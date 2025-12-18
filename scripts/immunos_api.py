#!/usr/bin/env python3
"""
IMMUNOS API Routes

RESTful API endpoints for all IMMUNOS cell types and components.
"""

import sys
import json
import hashlib
import requests
from pathlib import Path
from datetime import datetime, timedelta
from flask import request, jsonify, g, current_app
from threading import Thread

# Import IMMUNOS modules
sys.path.insert(0, str(Path(__file__).parent))

try:
    from immunos_nk_scan import NKCellScanner, Anomaly
    from immunos_memory import MemoryStore
    from immunos_todo import TodoStore
except ImportError as e:
    print(f"⚠ Warning: Could not import IMMUNOS modules: {e}")


# Utility Functions
# =================

def api_response(success=True, data=None, error=None):
    """Standardized API response format"""
    return jsonify({
        'success': success,
        'data': data,
        'error': error,
        'timestamp': datetime.now().isoformat()
    })


def generate_anomaly_hash(file_path: str, line_number: int, pattern: str) -> str:
    """Generate unique hash for an anomaly"""
    hash_input = f"{file_path}:{line_number}:{pattern}"
    return hashlib.sha256(hash_input.encode()).hexdigest()[:16]


def get_db():
    """Get database connection"""
    return g.get('db')


# System Health & Status Endpoints
# =================================

def register_routes(app, socketio):
    """Register all API routes with the Flask app"""

    @app.route('/api/health')
    def api_health():
        """Get overall system health"""
        try:
            db = get_db()

            # Get latest health score
            row = db.execute('''
                SELECT * FROM health_scores
                ORDER BY date DESC LIMIT 1
            ''').fetchone()

            if row:
                overall_score = row['overall_score']
                status = 'secure' if overall_score >= 80 else 'warning' if overall_score >= 50 else 'critical'
                threat_level = 'low' if overall_score >= 80 else 'medium' if overall_score >= 50 else 'high'

                health_data = {
                    'overall_health': overall_score,
                    'status': status,
                    'threat_level': threat_level,
                    'last_scan': row['created_at'],
                    'breakdown': {
                        'security': row['security_score'],
                        'quality': row['quality_score'],
                        'structure': row['structure_score'],
                        'docs': row['docs_score'],
                        'tasks': row['task_score']
                    },
                    'counts': {
                        'critical': row['critical_count'],
                        'warning': row['warning_count'],
                        'info': row['info_count']
                    }
                }
            else:
                # No health data yet
                health_data = {
                    'overall_health': 100,
                    'status': 'unknown',
                    'threat_level': 'unknown',
                    'last_scan': None,
                    'breakdown': {},
                    'counts': {'critical': 0, 'warning': 0, 'info': 0}
                }

            return api_response(data=health_data)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/status')
    def api_status():
        """Get status of all IMMUNOS components"""
        try:
            db = get_db()
            base_path = Path(current_app.config['BASE_PATH'])

            # NK Cell status
            nk_status = db.execute('''
                SELECT * FROM scan_history
                WHERE status = 'completed'
                ORDER BY completed_at DESC LIMIT 1
            ''').fetchone()

            # Memory status
            memory_path = base_path / '.immunos' / 'memory'
            memory_index_path = memory_path / 'index.json'
            memory_count = 0
            if memory_index_path.exists():
                with open(memory_index_path) as f:
                    index = json.load(f)
                    memory_count = len(index.get('memories', []))

            # Todo status
            todo_path = base_path / 'todo' / '.index' / 'todos.json'
            todo_stats = {'active': 0, 'completed': 0}
            if todo_path.exists():
                with open(todo_path) as f:
                    todos = json.load(f)
                    todo_stats['active'] = sum(1 for t in todos if t.get('status') != 'done')
                    todo_stats['completed'] = sum(1 for t in todos if t.get('status') == 'done')

            # Snapshot status
            snapshot_path = base_path / '.immunos' / 'snapshots'
            snapshot_count = len(list(snapshot_path.glob('snap_*.json'))) if snapshot_path.exists() else 0

            # Ollama status
            ollama_url = current_app.config.get('OLLAMA_URL', 'http://localhost:11434')
            ollama_status = 'offline'
            try:
                resp = requests.get(f"{ollama_url}/api/tags", timeout=2)
                ollama_status = 'online' if resp.status_code == 200 else 'offline'
            except:
                pass

            components = {
                'nk_cell': {
                    'status': 'active' if nk_status else 'idle',
                    'last_run': nk_status['completed_at'] if nk_status else None,
                    'anomalies': {
                        'high': nk_status['high_severity'] if nk_status else 0,
                        'medium': nk_status['medium_severity'] if nk_status else 0,
                        'low': nk_status['low_severity'] if nk_status else 0
                    }
                },
                't_cell': {
                    'status': 'active',
                    'memories': {'active': memory_count}
                },
                'todo': {
                    'status': 'active',
                    'stats': todo_stats
                },
                'snapshots': {
                    'status': 'active',
                    'count': snapshot_count
                },
                'ollama': {
                    'status': ollama_status,
                    'url': ollama_url
                }
            }

            return api_response(data={'components': components})

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/stats')
    def api_stats():
        """Get dashboard statistics"""
        try:
            db = get_db()

            # Recent anomalies
            anomalies = db.execute('''
                SELECT COUNT(*) as total,
                       SUM(CASE WHEN severity = 'HIGH' THEN 1 ELSE 0 END) as high,
                       SUM(CASE WHEN severity = 'MEDIUM' THEN 1 ELSE 0 END) as medium,
                       SUM(CASE WHEN severity = 'LOW' THEN 1 ELSE 0 END) as low
                FROM anomaly_actions
                WHERE action NOT IN ('fixed', 'false_positive')
            ''').fetchone()

            stats = {
                'critical': anomalies['high'] or 0,
                'warnings': anomalies['medium'] or 0,
                'info': anomalies['low'] or 0,
                'fixed': db.execute('SELECT COUNT(*) as c FROM anomaly_actions WHERE action = "fixed"').fetchone()['c']
            }

            return api_response(data=stats)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    # NK Cell (Scanner) Endpoints
    # ============================

    @app.route('/api/nk/scan', methods=['POST'])
    def api_nk_scan():
        """Trigger a new NK Cell scan"""
        try:
            data = request.json or {}
            scan_type = data.get('scan_type', 'full')
            use_ollama = data.get('use_ollama', False)
            trigger = data.get('trigger', 'manual')

            db = get_db()

            # Create scan record
            cursor = db.execute('''
                INSERT INTO scan_history (scan_type, started_at, status, trigger)
                VALUES (?, ?, 'running', ?)
            ''', (scan_type, datetime.now(), trigger))
            scan_id = cursor.lastrowid
            db.commit()

            # Run scan in background thread
            def run_scan():
                from flask import current_app as app_instance
                import sqlite3

                try:
                    base_path = Path(app_instance.config['BASE_PATH'])
                    scanner = NKCellScanner()

                    # Progress callback
                    def progress_callback(percent, status):
                        socketio.emit('scan_progress', {
                            'scan_id': scan_id,
                            'percent': percent,
                            'status': status
                        })

                    # Scan files (simplified - would need full implementation)
                    anomalies = []
                    for py_file in base_path.glob('**/*.py'):
                        if any(skip in str(py_file) for skip in ['.venv', 'node_modules', '.git']):
                            continue

                        try:
                            content = py_file.read_text()
                            scanner.check_security_patterns(str(py_file), content)
                        except:
                            pass

                    anomalies = scanner.anomalies

                    # Calculate severity counts
                    high_count = sum(1 for a in anomalies if a.severity == 'HIGH')
                    medium_count = sum(1 for a in anomalies if a.severity == 'MEDIUM')
                    low_count = sum(1 for a in anomalies if a.severity == 'LOW')

                    # Update scan record with direct connection
                    db_conn = sqlite3.connect(app_instance.config['DATABASE'])
                    db_conn.execute('''
                        UPDATE scan_history
                        SET completed_at = ?,
                            status = 'completed',
                            results_json = ?,
                            total_anomalies = ?,
                            high_severity = ?,
                            medium_severity = ?,
                            low_severity = ?
                        WHERE id = ?
                    ''', (datetime.now().isoformat(), json.dumps([vars(a) for a in anomalies]),
                          len(anomalies), high_count, medium_count, low_count, scan_id))
                    db_conn.commit()
                    db_conn.close()

                    # Emit completion event
                    socketio.emit('scan_completed', {
                        'scan_id': scan_id,
                        'total': len(anomalies),
                        'high': high_count,
                        'medium': medium_count,
                        'low': low_count
                    })

                except Exception as e:
                    # Update with error using direct connection
                    try:
                        db_conn = sqlite3.connect(app_instance.config['DATABASE'])
                        db_conn.execute('''
                            UPDATE scan_history
                            SET status = 'failed', error = ?
                            WHERE id = ?
                        ''', (str(e), scan_id))
                        db_conn.commit()
                        db_conn.close()
                    except:
                        pass

                    socketio.emit('scan_error', {'scan_id': scan_id, 'error': str(e)})

            # Start background scan
            thread = Thread(target=run_scan)
            thread.daemon = True
            thread.start()

            return api_response(data={'scan_id': scan_id, 'status': 'started'})

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/nk/results')
    def api_nk_results():
        """Get latest scan results"""
        try:
            db = get_db()
            scan = db.execute('''
                SELECT * FROM scan_history
                WHERE status = 'completed'
                ORDER BY completed_at DESC LIMIT 1
            ''').fetchone()

            if not scan:
                return api_response(data=None)

            results = {
                'scan_id': scan['id'],
                'scan_type': scan['scan_type'],
                'completed_at': scan['completed_at'],
                'total': scan['total_anomalies'],
                'high': scan['high_severity'],
                'medium': scan['medium_severity'],
                'low': scan['low_severity'],
                'anomalies': json.loads(scan['results_json']) if scan['results_json'] else []
            }

            return api_response(data=results)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/nk/anomalies')
    def api_nk_anomalies():
        """List anomalies with filtering"""
        try:
            severity = request.args.get('severity')
            status = request.args.get('status')  # 'unresolved', 'fixed', 'ignored'

            db = get_db()
            query = 'SELECT * FROM anomaly_actions WHERE 1=1'
            params = []

            if severity:
                query += ' AND severity = ?'
                params.append(severity.upper())

            if status == 'unresolved':
                query += ' AND action NOT IN ("fixed", "false_positive")'
            elif status:
                query += ' AND action = ?'
                params.append(status)

            query += ' ORDER BY action_date DESC'

            rows = db.execute(query, params).fetchall()
            anomalies = [dict(row) for row in rows]

            return api_response(data=anomalies)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/nk/anomalies/<anomaly_id>/action', methods=['POST'])
    def api_nk_anomaly_action(anomaly_id):
        """Mark anomaly action (fix/ignore/false_positive)"""
        try:
            data = request.json
            action = data.get('action')  # 'fixed', 'ignored', 'false_positive'
            note = data.get('note', '')

            if action not in ['fixed', 'ignored', 'false_positive']:
                return api_response(success=False, error='Invalid action'), 400

            db = get_db()
            db.execute('''
                UPDATE anomaly_actions
                SET action = ?, note = ?, action_date = ?
                WHERE anomaly_hash = ?
            ''', (action, note, datetime.now(), anomaly_id))
            db.commit()

            return api_response(data={'anomaly_id': anomaly_id, 'action': action})

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    # T Cell (Memory) Endpoints
    # ==========================

    @app.route('/api/memory/list')
    def api_memory_list():
        """List memories with filtering"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            memory_store = MemoryStore(base_path)

            # Get all memories from index
            memories = memory_store.list_memories()

            return api_response(data=memories)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/memory/<memory_id>')
    def api_memory_get(memory_id):
        """Get specific memory"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            memory_store = MemoryStore(base_path)

            mem = memory_store.get_memory(memory_id)
            if not mem:
                return api_response(success=False, error='Memory not found'), 404

            return api_response(data=mem)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/memory/decay', methods=['POST'])
    def api_memory_decay():
        """Trigger memory decay"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            memory_store = MemoryStore(base_path)

            # Run decay process
            result = memory_store.decay_memories()

            return api_response(data={'expired_count': result.get('cleaned', 0)})

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    # Ollama (LLM) Endpoints
    # ======================

    @app.route('/api/ollama/status')
    def api_ollama_status():
        """Get Ollama server status"""
        try:
            ollama_url = current_app.config.get('OLLAMA_URL', 'http://localhost:11434')

            try:
                resp = requests.get(f"{ollama_url}/api/tags", timeout=2)
                if resp.status_code == 200:
                    models = resp.json().get('models', [])
                    return api_response(data={
                        'status': 'online',
                        'url': ollama_url,
                        'models': models
                    })
            except:
                pass

            return api_response(data={'status': 'offline', 'url': ollama_url})

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    @app.route('/api/ollama/test', methods=['POST'])
    def api_ollama_test():
        """Test a model with a prompt"""
        try:
            data = request.json
            model = data.get('model')
            prompt = data.get('prompt')

            if not model or not prompt:
                return api_response(success=False, error='Model and prompt required'), 400

            ollama_url = current_app.config.get('OLLAMA_URL', 'http://localhost:11434')

            start_time = datetime.now()
            resp = requests.post(
                f"{ollama_url}/api/generate",
                json={'model': model, 'prompt': prompt, 'stream': False},
                timeout=30
            )
            end_time = datetime.now()

            if resp.status_code == 200:
                result = resp.json()
                response_time = (end_time - start_time).total_seconds() * 1000

                return api_response(data={
                    'response': result.get('response'),
                    'response_time_ms': response_time,
                    'model': model
                })

            return api_response(success=False, error=f'Ollama error: {resp.status_code}'), 500

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    # Activity Log Endpoint
    # =====================

    @app.route('/api/activity')
    def api_activity():
        """Get recent activity"""
        try:
            limit = int(request.args.get('limit', 20))

            db = get_db()
            rows = db.execute('''
                SELECT * FROM activity_log
                ORDER BY created_at DESC
                LIMIT ?
            ''', (limit,)).fetchall()

            activities = [dict(row) for row in rows]

            return api_response(data=activities)

        except Exception as e:
            return api_response(success=False, error=str(e)), 500


    print("✓ API routes registered successfully")

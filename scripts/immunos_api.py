#!/usr/bin/env python3
"""
IMMUNOS API Routes

RESTful API endpoints for all IMMUNOS cell types and components.
"""

import os
import sys
import json
import hashlib
import requests
import random
import subprocess
import time
from pathlib import Path
from datetime import datetime, timedelta
from flask import request, jsonify, g, current_app
from threading import Thread, Lock

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


def load_json_file(path: Path) -> dict:
    """Load JSON file if present."""
    if not path.exists():
        return {}
    with path.open(encoding='utf-8') as handle:
        return json.load(handle)


def load_orchestrator_config(base_path: Path) -> dict:
    """Load orchestrator config, creating defaults if missing."""
    config_dir = base_path / ".immunos" / "config"
    config_path = config_dir / "orchestrator.json"
    default_config = {
        "connectivity": "online",
        "online_backend": {
            "provider": "claude_code",
            "model": "",
            "base_url": "",
            "api_key_env": "ANTHROPIC_AUTH_TOKEN"
        },
        "offline_backend": {
            "provider": "ollama",
            "model": "qwen2.5-coder:7b",
            "base_url": "http://localhost:11434",
            "api_key_env": ""
        }
    }

    if not config_path.exists():
        config_dir.mkdir(parents=True, exist_ok=True)
        with config_path.open("w", encoding="utf-8") as handle:
            json.dump(default_config, handle, indent=2)
        return default_config

    try:
        with config_path.open(encoding="utf-8") as handle:
            data = json.load(handle)
    except json.JSONDecodeError:
        data = {}

    merged = default_config.copy()
    merged.update(data)
    merged["online_backend"].update(data.get("online_backend", {}))
    merged["offline_backend"].update(data.get("offline_backend", {}))
    return merged


def save_orchestrator_config(base_path: Path, config: dict) -> dict:
    """Persist orchestrator config and return merged config."""
    config_dir = base_path / ".immunos" / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    merged = load_orchestrator_config(base_path)

    merged.update({k: v for k, v in config.items() if k in {"connectivity", "online_backend", "offline_backend"}})
    if "online_backend" in config:
        merged["online_backend"].update(config["online_backend"])
    if "offline_backend" in config:
        merged["offline_backend"].update(config["offline_backend"])

    with (config_dir / "orchestrator.json").open("w", encoding="utf-8") as handle:
        json.dump(merged, handle, indent=2)
    return merged


def resolve_orchestrator_backend(config: dict) -> dict:
    """Select active backend based on connectivity."""
    connectivity = config.get("connectivity", "online")
    online_backend = config.get("online_backend") or {}
    offline_backend = config.get("offline_backend") or {}

    def backend_available(backend: dict) -> tuple[bool, str]:
        provider = backend.get("provider")
        if provider == "ollama":
            return True, ""

        if provider == "claude_code":
            api_key_env = backend.get("api_key_env") or "ANTHROPIC_AUTH_TOKEN"
            if not os.getenv(api_key_env):
                return False, "missing_api_key"
            model = backend.get("model") or os.getenv("ANTHROPIC_MODEL")
            if not model:
                return False, "missing_model"
            return True, ""

        if provider == "chatgpt":
            api_key_env = backend.get("api_key_env") or "OPENAI_API_KEY"
            if not os.getenv(api_key_env):
                return False, "missing_api_key"
            model = backend.get("model") or os.getenv("OPENAI_MODEL")
            if not model:
                return False, "missing_model"
            return True, ""

        if provider == "openrouter":
            api_key_env = backend.get("api_key_env") or "OPENROUTER_API_KEY"
            if not os.getenv(api_key_env):
                return False, "missing_api_key"
            model = backend.get("model") or os.getenv("OPENROUTER_MODEL")
            if not model:
                return False, "missing_model"
            return True, ""

        if provider == "local_server":
            if not backend.get("base_url"):
                return False, "missing_base_url"
            if not backend.get("model"):
                return False, "missing_model"
            return True, ""

        return False, "unknown_provider"

    selected = online_backend if connectivity == "online" else offline_backend
    fallback_reason = ""

    if connectivity == "online":
        ok, reason = backend_available(selected)
        if not ok:
            selected = offline_backend
            connectivity = "offline"
            fallback_reason = reason

    return {
        "connectivity": connectivity,
        "provider": selected.get("provider"),
        "model": selected.get("model"),
        "base_url": selected.get("base_url"),
        "api_key_env": selected.get("api_key_env"),
        "fallback_reason": fallback_reason
    }


def summarize_agent_state(agents_dir: Path) -> dict:
    """Summarize tracked agent artifacts for dashboard display."""
    domains = {
        "emotion": {},
        "hallucination": {
            "bcell": agents_dir / "hallucination_bcell.json",
            "nk": agents_dir / "hallucination_nk.json"
        },
        "network": {
            "nk_enhanced": agents_dir / "network_nk_enhanced.json"
        },
        "code": {},
        "research": {
            "bcell": agents_dir / "research_bcell.json",
            "nk": agents_dir / "research_nk.json"
        }
    }

    summary = {}
    for domain, paths in domains.items():
        stats = {
            "bcell_patterns": 0,
            "nk_detectors": 0,
            "nk_self_patterns": 0,
            "trained": False,
            "last_trained": None
        }
        latest_mtime = 0.0

        bcell_path = paths.get("bcell")
        if bcell_path and bcell_path.exists():
            data = load_json_file(bcell_path)
            stats["bcell_patterns"] = len(data.get("patterns", []))
            stats["trained"] = True
            latest_mtime = max(latest_mtime, bcell_path.stat().st_mtime)

        nk_path = paths.get("nk")
        if nk_path and nk_path.exists():
            data = load_json_file(nk_path)
            stats["nk_detectors"] = len(data.get("detectors", []))
            stats["nk_self_patterns"] = len(data.get("self_patterns", []))
            stats["trained"] = True
            latest_mtime = max(latest_mtime, nk_path.stat().st_mtime)

        nk_enhanced_path = paths.get("nk_enhanced")
        if nk_enhanced_path and nk_enhanced_path.exists():
            data = load_json_file(nk_enhanced_path)
            detector_sets = data.get("detector_sets", {})
            stats["nk_detectors"] = sum(len(item.get("detectors", [])) for item in detector_sets.values())
            self_patterns = data.get("self_patterns_by_class", {})
            stats["nk_self_patterns"] = sum(len(items) for items in self_patterns.values())
            stats["trained"] = True
            latest_mtime = max(latest_mtime, nk_enhanced_path.stat().st_mtime)

        if latest_mtime:
            stats["last_trained"] = datetime.fromtimestamp(latest_mtime).strftime("%Y-%m-%d %H:%M")

        summary[domain] = stats

    return summary


def load_todos_index(base_path: Path) -> dict:
    """Load the todo index file if present."""
    todo_index = base_path / "todo" / ".index" / "todos.json"
    if not todo_index.exists():
        return {}
    try:
        with todo_index.open(encoding="utf-8") as handle:
            return json.load(handle)
    except json.JSONDecodeError:
        return {}


def load_thymus_intake(path: Path) -> list:
    """Load thymus intake queue from disk."""
    if not path.exists():
        return []
    try:
        with path.open(encoding="utf-8") as handle:
            data = json.load(handle)
            if isinstance(data, list):
                return data
    except json.JSONDecodeError:
        return []
    return []


def save_thymus_intake(path: Path, items: list) -> None:
    """Save thymus intake queue to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(items, handle, indent=2)


def load_agent_templates(base_path: Path) -> dict:
    """Load agent templates for the Agent Foundry."""
    template_path = base_path / "immunos-mcp" / "src" / "immunos_mcp" / "config" / "agent_templates.json"
    if not template_path.exists():
        return {"version": "1.0", "templates": []}
    try:
        with template_path.open(encoding="utf-8") as handle:
            data = json.load(handle)
            if isinstance(data, dict):
                return data
    except json.JSONDecodeError:
        pass
    return {"version": "1.0", "templates": []}


def list_foundry_agents(foundry_dir: Path, limit: int = 25) -> list:
    """List recent Agent Foundry stubs."""
    if not foundry_dir.exists():
        return []
    files = sorted(foundry_dir.glob("agent_*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    agents = []
    for path in files[:limit]:
        try:
            with path.open(encoding="utf-8") as handle:
                data = json.load(handle)
                if isinstance(data, dict):
                    data["path"] = str(path)
                    agents.append(data)
        except (OSError, json.JSONDecodeError):
            continue
    return agents


def save_foundry_agent(foundry_dir: Path, payload: dict) -> dict:
    """Save a foundry agent stub to disk."""
    foundry_dir.mkdir(parents=True, exist_ok=True)
    agent_id = payload.get("id") or f"agent_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    payload["id"] = agent_id
    path = foundry_dir / f"{agent_id}.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    payload["path"] = str(path)
    return payload


def add_thymus_intake(path: Path, entry: dict) -> dict:
    """Append an entry to the thymus intake queue."""
    items = load_thymus_intake(path)
    items.insert(0, entry)
    save_thymus_intake(path, items)
    return entry


def update_thymus_entry(path: Path, entry_id: str, updates: dict) -> dict:
    """Update a thymus intake entry by id."""
    items = load_thymus_intake(path)
    updated = None
    for item in items:
        if item.get("id") == entry_id:
            item.update(updates)
            updated = item
            break
    save_thymus_intake(path, items)
    return updated or {}


def next_thymus_entry(path: Path) -> dict:
    """Return the next queued thymus entry."""
    items = load_thymus_intake(path)
    for item in items:
        if item.get("status") == "queued":
            return item
    return {}


def summarize_issues(index: dict, project: str = None) -> dict:
    """Summarize issues from the todo index."""
    todos = index.get("todos", {})
    stats = {
        "total": 0,
        "active": 0,
        "by_status": {},
        "by_priority": {},
        "overdue": 0,
        "due_today": 0,
        "due_this_week": 0
    }

    today = datetime.now().date()
    week_from_now = today + timedelta(days=7)

    for todo in todos.values():
        if project:
            tags = todo.get("tags", []) or []
            if todo.get("project") != project and project not in tags:
                continue
        stats["total"] += 1

        status = todo.get("status", "unknown")
        priority = todo.get("priority", "medium")
        stats["by_status"][status] = stats["by_status"].get(status, 0) + 1
        stats["by_priority"][priority] = stats["by_priority"].get(priority, 0) + 1

        if status != "done":
            stats["active"] += 1

        due_str = todo.get("due")
        if status != "done" and due_str:
            try:
                due_date = datetime.fromisoformat(due_str).date()
            except ValueError:
                continue

            if due_date < today:
                stats["overdue"] += 1
            elif due_date == today:
                stats["due_today"] += 1
            elif due_date <= week_from_now:
                stats["due_this_week"] += 1

    return stats


def list_issues(index: dict, status: str = None, limit: int = 25, project: str = None) -> list:
    """Return a sorted list of issues."""
    todos = list(index.get("todos", {}).values())
    if project:
        filtered = []
        for todo in todos:
            tags = todo.get("tags", []) or []
            if todo.get("project") == project or project in tags:
                filtered.append(todo)
        todos = filtered
    if status:
        todos = [todo for todo in todos if todo.get("status") == status]

    todos.sort(key=lambda item: item.get("created", ""), reverse=True)
    return todos[:limit]


# System Health & Status Endpoints
# =================================

def register_routes(app, socketio):
    """Register all API routes with the Flask app"""
    thymus_lock = Lock()
    thymus_state = {"running": False, "current_id": None, "paused": False}

    def build_training_command(entry: dict, base_path: Path, agents_dir: Path) -> dict:
        """Build training command based on thymus entry."""
        domain = entry.get("domain")
        dataset_path = entry.get("dataset_path") or None
        sample_count = entry.get("sample_count")
        max_samples = str(int(sample_count)) if sample_count else None

        env = os.environ.copy()
        env["PYTHONPATH"] = f"{base_path / 'immunos-mcp' / 'src'}{os.pathsep}{env.get('PYTHONPATH', '')}"

        if domain == "hallucination":
            cmd = [sys.executable, "-m", "immunos_mcp.training.hallucination_trainer"]
            if dataset_path:
                cmd += ["--truthfulqa", dataset_path]
            if max_samples:
                cmd += ["--max-samples", max_samples]
            cmd += ["--save-bcell", str(agents_dir / "hallucination_bcell.json")]
            cmd += ["--save-nk", str(agents_dir / "hallucination_nk.json")]
            return {"cmd": cmd, "env": env, "dataset": dataset_path or "TruthfulQA", "domain": domain}

        if domain == "network":
            cmd = [sys.executable, "-m", "immunos_mcp.training.network_trainer"]
            if dataset_path:
                cmd += ["--kdd", dataset_path]
            if max_samples:
                cmd += ["--max-samples", max_samples]
            cmd += ["--save-nk", str(agents_dir / "network_nk_enhanced.json")]
            return {"cmd": cmd, "env": env, "dataset": dataset_path or "KDDTrain+", "domain": domain}

        if domain == "research":
            cmd = [sys.executable, "-m", "immunos_mcp.training.research_trainer"]
            if dataset_path:
                cmd += ["--claims", dataset_path]
            if max_samples:
                cmd += ["--max-samples", max_samples]
            cmd += ["--save-bcell", str(agents_dir / "research_bcell.json")]
            cmd += ["--save-nk", str(agents_dir / "research_nk.json")]
            return {"cmd": cmd, "env": env, "dataset": dataset_path or "SciFact", "domain": domain}

        if domain == "code":
            cmd = [sys.executable, "-m", "immunos_mcp.training.code_trainer"]
            if dataset_path:
                cmd += ["--code-dir", dataset_path]
            if max_samples:
                cmd += ["--max-samples", max_samples]
            cmd += ["--save-bcell", str(agents_dir / "code_bcell.json")]
            cmd += ["--save-nk", str(agents_dir / "code_nk.json")]
            return {"cmd": cmd, "env": env, "dataset": dataset_path or "Code Samples", "domain": domain}

        if domain == "emotion":
            cmd = [sys.executable, "-m", "immunos_mcp.training.emotion_trainer"]
            if dataset_path:
                cmd += ["--features", dataset_path]
            if max_samples:
                cmd += ["--max-samples", max_samples]
            cmd += ["--save-bcell", str(agents_dir / "emotion_bcell.json")]
            cmd += ["--save-nk", str(agents_dir / "emotion_nk.json")]
            return {"cmd": cmd, "env": env, "dataset": dataset_path or "Emotion Features", "domain": domain}

        return {"cmd": None, "env": env, "dataset": dataset_path or "manual", "domain": domain}

    def run_thymus_job(entry: dict, base_path: Path):
        intake_path = base_path / ".immunos" / "thymus" / "intake.json"
        agents_dir = base_path / "immunos-mcp" / ".immunos" / "agents"
        agents_dir.mkdir(parents=True, exist_ok=True)

        command_info = build_training_command(entry, base_path, agents_dir)
        cmd = command_info["cmd"]

        update_thymus_entry(intake_path, entry["id"], {
            "status": "running" if cmd else "manual",
            "started_at": datetime.now().isoformat(),
            "dataset_label": command_info["dataset"]
        })

        socketio.emit('training_started', {
            'domain': entry.get("domain", "unknown"),
            'dataset': command_info["dataset"],
            'num_samples': entry.get("sample_count") or 0,
            'timestamp': datetime.now().isoformat()
        })

        if not cmd:
            update_thymus_entry(intake_path, entry["id"], {
                "status": "manual",
                "completed_at": datetime.now().isoformat(),
                "error": "No training pipeline available"
            })
            socketio.emit('training_complete', {
                'domain': entry.get("domain", "unknown"),
                'detectors_created': 0,
                'final_accuracy': 0.0,
                'duration_seconds': 0,
                'timestamp': datetime.now().isoformat()
            })
            with thymus_lock:
                thymus_state["running"] = False
                thymus_state["current_id"] = None
            return

        start_time = time.time()
        progress = 0
        try:
            process = subprocess.Popen(
                cmd,
                cwd=str(base_path),
                env=command_info["env"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            while process.poll() is None:
                time.sleep(2)
                progress = min(95, progress + 5)
                socketio.emit('training_progress', {
                    'domain': entry.get("domain", "unknown"),
                    'current': progress,
                    'total': 100,
                    'percent': progress,
                    'detectors_created': 0,
                    'current_accuracy': 0.0,
                    'timestamp': datetime.now().isoformat()
                })

            stdout, stderr = process.communicate()
            duration = int(time.time() - start_time)
            status = "completed" if process.returncode == 0 else "failed"

            domain_stats = summarize_agent_state(agents_dir)
            domain_entry = domain_stats.get(entry.get("domain", ""), {})

            update_thymus_entry(intake_path, entry["id"], {
                "status": status,
                "completed_at": datetime.now().isoformat(),
                "output": stdout[-2000:] if stdout else "",
                "error": stderr[-2000:] if stderr else "",
                "detectors": domain_entry.get("nk_detectors", 0),
                "bcell_patterns": domain_entry.get("bcell_patterns", 0),
                "nk_self_patterns": domain_entry.get("nk_self_patterns", 0)
            })

            socketio.emit('training_complete', {
                'domain': entry.get("domain", "unknown"),
                'detectors_created': 0,
                'final_accuracy': 0.0,
                'duration_seconds': duration,
                'timestamp': datetime.now().isoformat()
            })
        finally:
            with thymus_lock:
                thymus_state["running"] = False
                thymus_state["current_id"] = None

    def maybe_start_thymus_queue(base_path: Path):
        intake_path = base_path / ".immunos" / "thymus" / "intake.json"
        with thymus_lock:
            if thymus_state["running"] or thymus_state.get("paused"):
                return
            next_entry = next_thymus_entry(intake_path)
            if not next_entry:
                return
            thymus_state["running"] = True
            thymus_state["current_id"] = next_entry.get("id")

        thread = Thread(target=run_thymus_job, args=(next_entry, base_path), daemon=True)
        thread.start()

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


    # ========================================================================
    # MODEL MANAGEMENT API (NEW - for dashboard)
    # ========================================================================

    # Import new modules
    try:
        from immunos_model_manager import OllamaModelManager as ModelManager
        from immunos_token_tracker import TokenTracker
        from immunos_routing import ModelRouter
    except ImportError as e:
        print(f"⚠ Warning: Could not import new model management modules: {e}")
        ModelManager = None
        TokenTracker = None
        ModelRouter = None

    @app.route('/api/models/list')
    def api_models_list():
        """GET - List all available Ollama models (enhanced)"""
        try:
            if ModelManager:
                manager = ModelManager()
                models = manager.list_models()
                return jsonify(models)
            else:
                # Fallback to basic check
                return api_ollama_status()
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/models/<model_name>/load', methods=['POST'])
    def api_model_load(model_name):
        """POST - Load model into memory"""
        try:
            if not ModelManager:
                return jsonify({'error': 'Model manager not available'}), 503

            manager = ModelManager()
            success = manager.load_model(model_name)

            if success:
                # Emit WebSocket event
                socketio.emit('model_loaded', {
                    'model': model_name,
                    'timestamp': datetime.now().isoformat()
                })
                return jsonify({'success': True, 'model': model_name})

            return jsonify({'success': False, 'error': 'Failed to load model'}), 500
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/models/<model_name>/unload', methods=['POST'])
    def api_model_unload(model_name):
        """POST - Unload model from memory"""
        try:
            if not ModelManager:
                return jsonify({'error': 'Model manager not available'}), 503

            manager = ModelManager()
            success = manager.unload_model(model_name)

            if success:
                socketio.emit('model_unloaded', {
                    'model': model_name,
                    'timestamp': datetime.now().isoformat()
                })
                return jsonify({'success': True, 'model': model_name})

            return jsonify({'success': False, 'error': 'Failed to unload model'}), 500
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/models/status')
    def api_models_status():
        """GET - Get status of all models"""
        try:
            if not ModelManager:
                return jsonify({'error': 'Model manager not available'}), 503

            manager = ModelManager()
            status = manager.get_all_status()
            return jsonify(status)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/ollama/start', methods=['POST'])
    def api_ollama_start():
        """POST - Start Ollama server"""
        try:
            if not ModelManager:
                return jsonify({'error': 'Model manager not available'}), 503

            manager = ModelManager()
            success = manager.start_ollama_server()

            if success:
                socketio.emit('ollama_server_started', {
                    'timestamp': datetime.now().isoformat()
                })
                return jsonify({'success': True})

            return jsonify({'success': False, 'error': 'Failed to start server'}), 500
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/ollama/stop', methods=['POST'])
    def api_ollama_stop():
        """POST - Stop Ollama server"""
        try:
            if not ModelManager:
                return jsonify({'error': 'Model manager not available'}), 503

            manager = ModelManager()
            success = manager.stop_ollama_server()

            if success:
                socketio.emit('ollama_server_stopped', {
                    'timestamp': datetime.now().isoformat()
                })
                return jsonify({'success': True})

            return jsonify({'success': False, 'error': 'Failed to stop server'}), 500
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    # ========================================================================
    # TOKEN TRACKING API
    # ========================================================================

    @app.route('/api/tokens/session/<session_id>')
    def api_tokens_session(session_id):
        """GET - Token usage for specific Claude session"""
        try:
            if not TokenTracker:
                return jsonify({'error': 'Token tracker not available'}), 503

            tracker = TokenTracker()
            usage = tracker.get_session_usage(session_id)
            return jsonify(usage)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/tokens/session/current')
    def api_tokens_session_current():
        """GET - Token usage for current session"""
        try:
            if not TokenTracker:
                return jsonify({'error': 'Token tracker not available'}), 503

            tracker = TokenTracker()
            session_id = tracker.get_or_create_session()
            usage = tracker.get_session_usage(session_id)
            return jsonify(usage)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/tokens/comparison')
    def api_tokens_comparison():
        """GET - Claude vs Ollama comparison (Query param: ?hours=24)"""
        try:
            if not TokenTracker:
                return jsonify({'error': 'Token tracker not available'}), 503

            hours = request.args.get('hours', 24, type=int)
            tracker = TokenTracker()
            data = tracker.get_provider_comparison(hours=hours)
            return jsonify(data)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/tokens/models')
    def api_tokens_models():
        """GET - Per-model token breakdown"""
        try:
            if not TokenTracker:
                return jsonify({'error': 'Token tracker not available'}), 503

            tracker = TokenTracker()
            models = tracker.get_model_breakdown()
            return jsonify(models)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/tokens/routing')
    def api_tokens_routing():
        """GET - Model routing statistics"""
        try:
            if not TokenTracker:
                return jsonify({'error': 'Token tracker not available'}), 503

            tracker = TokenTracker()
            stats = tracker.get_routing_stats()
            return jsonify(stats)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    # ========================================================================
    # ROUTING API
    # ========================================================================

    @app.route('/api/routing/config')
    def api_routing_config():
        """GET - Current routing configuration"""
        try:
            if not ModelRouter:
                return jsonify({'error': 'Model router not available'}), 503

            router = ModelRouter()
            config = router.get_routing_config()
            return jsonify(config)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/routing/config/update', methods=['POST'])
    def api_routing_config_update():
        """POST - Update routing configuration"""
        try:
            if not ModelRouter:
                return jsonify({'error': 'Model router not available'}), 503

            router = ModelRouter()
            config = request.get_json()
            success = router.update_routing_config(config)

            if success:
                return jsonify({'success': True})

            return jsonify({'success': False, 'error': 'Failed to update config'}), 500
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/routing/decide', methods=['POST'])
    def api_routing_decide():
        """POST - Get routing decision for a task"""
        try:
            if not ModelRouter:
                return jsonify({'error': 'Model router not available'}), 503

            router = ModelRouter()
            data = request.get_json()

            model, reason = router.route_task(
                task_type=data.get('task_type', 'unknown'),
                task_classification=data.get('task_classification', 'routine'),
                estimated_tokens=data.get('estimated_tokens')
            )

            return jsonify({'model': model, 'reason': reason})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/routing/status')
    def api_routing_status():
        """GET - Current routing status with recommendations"""
        try:
            if not ModelRouter:
                return jsonify({'error': 'Model router not available'}), 503

            router = ModelRouter()
            status = router.get_session_status()
            return jsonify(status)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/domains/status')
    def api_domains_status():
        """GET - Summary of tracked agent artifacts per domain"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            agents_dir = base_path / "immunos-mcp" / ".immunos" / "agents"
            domains = summarize_agent_state(agents_dir)
            return jsonify({'domains': domains})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/orchestrator/status')
    def api_orchestrator_status():
        """GET - Current orchestrator backend configuration"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            config = load_orchestrator_config(base_path)
            active = resolve_orchestrator_backend(config)
            return jsonify({'config': config, 'active_backend': active})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/orchestrator/config', methods=['POST'])
    def api_orchestrator_config():
        """POST - Update orchestrator configuration"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            payload = request.get_json() or {}
            config = save_orchestrator_config(base_path, payload)
            active = resolve_orchestrator_backend(config)
            return jsonify({'success': True, 'config': config, 'active_backend': active})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/thymus/intake')
    def api_thymus_intake_list():
        """GET - Thymus intake queue"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            intake_path = base_path / ".immunos" / "thymus" / "intake.json"
            limit = request.args.get("limit", 20, type=int)
            items = load_thymus_intake(intake_path)
            return jsonify({'items': items[:limit]})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/thymus/intake', methods=['POST'])
    def api_thymus_intake_create():
        """POST - Add entry to thymus intake queue"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            intake_path = base_path / ".immunos" / "thymus" / "intake.json"
            payload = request.get_json() or {}
            domain = payload.get("domain", "unknown")
            dataset_path = payload.get("dataset_path", "")
            sample_count = payload.get("sample_count")
            if not isinstance(sample_count, (int, float)):
                sample_count = None
            notes = payload.get("notes", "")

            entry = {
                "id": hashlib.sha256(f"{domain}:{dataset_path}:{datetime.now().isoformat()}".encode()).hexdigest()[:12],
                "domain": domain,
                "dataset_path": dataset_path,
                "sample_count": sample_count,
                "notes": notes,
                "status": "queued",
                "created_at": datetime.now().isoformat()
            }
            add_thymus_intake(intake_path, entry)
            maybe_start_thymus_queue(base_path)
            return jsonify({'success': True, 'entry': entry})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/thymus/status')
    def api_thymus_status():
        """GET - Thymus queue status"""
        try:
            with thymus_lock:
                status = {
                    "running": thymus_state.get("running", False),
                    "paused": thymus_state.get("paused", False),
                    "current_id": thymus_state.get("current_id")
                }
            return jsonify(status)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/thymus/control', methods=['POST'])
    def api_thymus_control():
        """POST - Control thymus queue (pause/resume/run_next)"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            payload = request.get_json() or {}
            action = payload.get("action")

            with thymus_lock:
                if action == "pause":
                    thymus_state["paused"] = True
                elif action == "resume":
                    thymus_state["paused"] = False
                elif action == "run_next":
                    thymus_state["paused"] = False
                else:
                    return jsonify({'error': 'Unknown action'}), 400

            if action in {"resume", "run_next"}:
                maybe_start_thymus_queue(base_path)

            return jsonify({
                "success": True,
                "running": thymus_state.get("running", False),
                "paused": thymus_state.get("paused", False),
                "current_id": thymus_state.get("current_id")
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/issues/summary')
    def api_issues_summary():
        """GET - Summary of IMMUNOS issues from todo system"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            index = load_todos_index(base_path)
            project = request.args.get("project")
            summary = summarize_issues(index, project=project)
            return jsonify({'summary': summary})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/issues/list')
    def api_issues_list():
        """GET - Recent IMMUNOS issues from todo system"""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            index = load_todos_index(base_path)
            status = request.args.get("status")
            limit = request.args.get("limit", 25, type=int)
            project = request.args.get("project")
            items = list_issues(index, status=status, limit=limit, project=project)
            return jsonify({'items': items})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/spleen/summary')
    def api_spleen_summary():
        """GET - Global anomaly and issue summary."""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            project = request.args.get("project")

            index = load_todos_index(base_path)
            issues_summary = summarize_issues(index, project=project)

            db = get_db()
            unresolved = db.execute('''
                SELECT COUNT(*) as total,
                       SUM(CASE WHEN severity = 'HIGH' THEN 1 ELSE 0 END) as high,
                       SUM(CASE WHEN severity = 'MEDIUM' THEN 1 ELSE 0 END) as medium,
                       SUM(CASE WHEN severity = 'LOW' THEN 1 ELSE 0 END) as low
                FROM anomaly_actions
                WHERE action NOT IN ('fixed', 'false_positive')
            ''').fetchone()

            resolved = db.execute('''
                SELECT COUNT(*) as total
                FROM anomaly_actions
                WHERE action IN ('fixed', 'false_positive')
            ''').fetchone()

            last_scan = db.execute('''
                SELECT completed_at
                FROM scan_history
                WHERE status = 'completed'
                ORDER BY completed_at DESC LIMIT 1
            ''').fetchone()

            summary = {
                "anomalies": {
                    "unresolved_total": unresolved["total"] if unresolved else 0,
                    "high": unresolved["high"] if unresolved else 0,
                    "medium": unresolved["medium"] if unresolved else 0,
                    "low": unresolved["low"] if unresolved else 0,
                    "resolved_total": resolved["total"] if resolved else 0
                },
                "issues": issues_summary,
                "last_scan": last_scan["completed_at"] if last_scan else None
            }

            return jsonify(summary)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/agents/templates')
    def api_agents_templates():
        """GET - List Agent Foundry templates."""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            data = load_agent_templates(base_path)
            return jsonify({'version': data.get('version'), 'templates': data.get('templates', [])})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/agents/foundry', methods=['GET', 'POST'])
    def api_agents_foundry():
        """GET - List Agent Foundry stubs. POST - Create a new stub."""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            foundry_dir = base_path / ".immunos" / "foundry"

            if request.method == 'GET':
                limit = request.args.get("limit", 25, type=int)
                items = list_foundry_agents(foundry_dir, limit=limit)
                return jsonify({'items': items})

            payload = request.get_json() or {}
            template_id = payload.get("template_id")
            name = payload.get("name") or "Agent Stub"
            domain = payload.get("domain") or "generic"
            notes = payload.get("notes") or ""

            templates = load_agent_templates(base_path).get("templates", [])
            template = next((item for item in templates if item.get("id") == template_id), None)
            if not template:
                return jsonify({'error': 'Template not found'}), 400

            stub = {
                "id": f"agent_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
                "template_id": template_id,
                "name": name,
                "domain": domain,
                "role": template.get("role", "unknown"),
                "description": template.get("description", ""),
                "defaults": template.get("defaults", {}),
                "notes": notes,
                "created_at": datetime.now().isoformat()
            }
            saved = save_foundry_agent(foundry_dir, stub)
            return jsonify({'item': saved})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/kb/index')
    def api_kb_index():
        """GET - List available KB pages."""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            kb_path = base_path / "docs" / "kb"
            if not kb_path.exists():
                return jsonify({'pages': []})

            pages = []
            for path in sorted(kb_path.glob("*.md")):
                name = path.stem
                title = name.replace("-", " ").title()
                try:
                    with path.open(encoding="utf-8") as handle:
                        for line in handle:
                            line = line.strip()
                            if line.startswith("# "):
                                title = line[2:].strip()
                                break
                except OSError:
                    pass
                pages.append({"name": name, "title": title})

            return jsonify({'pages': pages})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/kb/page')
    def api_kb_page():
        """GET - Return KB page content."""
        try:
            base_path = Path(current_app.config['BASE_PATH'])
            kb_path = base_path / "docs" / "kb"
            name = request.args.get("name", "")
            if not name.isidentifier() and "-" in name:
                name = name.replace("/", "").replace("..", "")

            path = (kb_path / f"{name}.md").resolve()
            if kb_path not in path.parents or not path.exists():
                return jsonify({'error': 'KB page not found'}), 404

            content = path.read_text(encoding="utf-8")
            title = name.replace("-", " ").title()
            for line in content.splitlines():
                if line.startswith("# "):
                    title = line[2:].strip()
                    break

            return jsonify({'name': name, 'title': title, 'content': content})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/events/detection', methods=['POST'])
    def api_events_detection():
        """POST - Emit detection event to dashboard"""
        try:
            payload = request.get_json() or {}
            payload.setdefault('timestamp', datetime.now().isoformat())
            socketio.emit('detection_event', payload)
            return jsonify({'success': True, 'event': payload})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/events/simulate', methods=['POST'])
    def api_events_simulate():
        """POST - Simulate a detection event for the dashboard map"""
        try:
            payload = request.get_json() or {}
            domain = payload.get('domain') or random.choice(['hallucination', 'network', 'research'])
            result = payload.get('result') or random.choice(['self', 'non_self', 'danger', 'uncertain'])
            confidence = round(random.uniform(0.45, 0.95), 2)
            danger_signal = round(random.uniform(0.0, 1.0), 2)

            event = {
                'domain': domain,
                'result': result,
                'confidence': confidence,
                'matched_detectors': payload.get('matched_detectors', []),
                'danger_signal': danger_signal,
                'timestamp': datetime.now().isoformat()
            }
            socketio.emit('detection_event', event)
            return jsonify({'success': True, 'event': event})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    # ========================================================================
    # CHAT API
    # ========================================================================

    @app.route('/api/chat', methods=['POST'])
    def api_chat():
        """POST - Chat with IMMUNOS models"""
        try:
            from immunos_chat import ImmunosChat
            data = request.json
            user_input = data.get('message', '')
            task_type = data.get('task_type', 'analysis')
            task_classification = data.get('task_classification', 'routine')
            mode = data.get('mode', 'orchestrator')
            domain = data.get('domain')
            connectivity = data.get('connectivity')

            chat = ImmunosChat()
            result = chat.chat(
                user_input,
                task_type=task_type,
                task_classification=task_classification,
                mode=mode,
                domain=domain,
                connectivity=connectivity
            )
            return jsonify(result)
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    # ========================================================================
    # HANDOFF API
    # ========================================================================

    @app.route('/api/handoff/check')
    def api_handoff_check():
        """GET - Check if handoff needed"""
        try:
            from immunos_handoff import ContextHandoff
            handoff = ContextHandoff()
            should_handoff, reason = handoff.should_handoff()
            session = handoff.tracker.get_session_usage(handoff.current_session_id)

            return jsonify({
                'should_handoff': should_handoff,
                'reason': reason,
                'token_usage': session.get('total_tokens', 0),
                'threshold_warning': 150000,
                'threshold_critical': 180000
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    @app.route('/api/handoff/create', methods=['POST'])
    def api_handoff_create():
        """POST - Create context handoff"""
        try:
            from immunos_handoff import ContextHandoff
            data = request.json
            handoff = ContextHandoff()

            filepath = handoff.save_handoff(
                conversation_history=data.get('conversation', []),
                current_task=data.get('current_task', 'Unknown'),
                files_being_worked_on=data.get('files', []),
                next_steps=data.get('next_steps', []),
                context=data.get('context', {}),
                reason=data.get('reason', 'manual')
            )

            return jsonify({'success': True, 'filepath': filepath})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    # ========================================================================
    # DASHBOARD ROUTES
    # ========================================================================

    @app.route('/monitor')
    def monitor_dashboard():
        """Real-time monitoring dashboard"""
        try:
            return app.send_static_file('monitor.html')
        except:
            from flask import render_template
            try:
                return render_template('monitor.html')
            except:
                return "Monitor dashboard under construction.", 404

    @app.route('/docs')
    def docs_page():
        """IMMUNOS Documentation and Knowledge Base"""
        from flask import render_template
        try:
            return render_template('docs.html')
        except:
            return "Documentation under construction.", 404

    print("✓ API routes registered successfully")
    print("  - Model Management: /api/models/* (7 endpoints)")
    print("  - Token Tracking: /api/tokens/* (5 endpoints)")
    print("  - Routing: /api/routing/* (4 endpoints)")
    print("  - Chat: /api/chat (1 endpoint)")
    print("  - Handoff: /api/handoff/* (2 endpoints)")
    print("  - Monitor: /monitor (dashboard)")
    print("  - Docs: /docs (knowledge base)")

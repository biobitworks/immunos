#!/usr/bin/env python3
"""
IMMUNOS Model Manager
=====================

Manages Ollama model lifecycle and monitoring for the IMMUNOS dashboard.

Features:
- Model discovery and status
- Load/unload models
- Ollama server management
- Database logging of events
"""

import os
import sys
import json
import time
import requests
import subprocess
import sqlite3
import psutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, asdict


@dataclass
class ModelStatus:
    """Status of an Ollama model"""
    name: str
    loaded: bool
    ram_mb: int
    last_used: Optional[str]
    total_requests: int
    avg_latency_ms: float


@dataclass
class ServerStatus:
    """Status of the Ollama server"""
    running: bool
    pid: Optional[int]
    uptime_seconds: int
    total_requests: int


class OllamaModelManager:
    """
    Manages Ollama model lifecycle and monitoring

    Responsibilities:
    - List available models
    - Load/unload models
    - Track model status and metrics
    - Start/stop Ollama server
    - Log events to database
    """

    def __init__(self, ollama_url: str = "http://localhost:11434", db_path: str = None):
        self.ollama_url = ollama_url
        self.db_path = db_path or str(Path.home() / ".immunos" / "db" / "immunos.db")

        # In-memory tracking of loaded models
        self.models_loaded: Dict[str, Dict[str, Any]] = {}

        # Server process info
        self.server_pid: Optional[int] = None
        self.server_start_time: Optional[float] = None

        # Initialize database connection
        self._init_db()

    # ========================================================================
    # DATABASE OPERATIONS
    # ========================================================================

    def _init_db(self):
        """Initialize database connection and ensure tables exist"""
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)

        with sqlite3.connect(self.db_path) as conn:
            # Create ollama_server_log table if not exists
            conn.execute("""
                CREATE TABLE IF NOT EXISTS ollama_server_log (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    event_type TEXT NOT NULL,
                    status TEXT NOT NULL,
                    pid INTEGER,
                    cpu_percent REAL,
                    memory_mb REAL,
                    details TEXT
                )
            """)

            # Create llm_model_usage table if not exists
            conn.execute("""
                CREATE TABLE IF NOT EXISTS llm_model_usage (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    model_name TEXT NOT NULL,
                    component TEXT NOT NULL,
                    prompt_tokens INTEGER DEFAULT 0,
                    response_tokens INTEGER DEFAULT 0,
                    response_time_ms INTEGER,
                    provider TEXT NOT NULL,
                    routing_reason TEXT
                )
            """)

            conn.commit()

    def _log_server_event(self, event_type: str, status: str, pid: Optional[int] = None,
                         cpu_percent: Optional[float] = None, memory_mb: Optional[float] = None,
                         details: str = None):
        """Log server event to database"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO ollama_server_log
                (timestamp, event_type, status, pid, cpu_percent, memory_mb, details)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                datetime.now().isoformat(),
                event_type,
                status,
                pid,
                cpu_percent,
                memory_mb,
                details
            ))
            conn.commit()

    def _get_model_usage_stats(self, model_name: str) -> Dict[str, Any]:
        """Get usage statistics for a model from database"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT
                    COUNT(*) as total_requests,
                    AVG(response_time_ms) as avg_latency_ms,
                    MAX(timestamp) as last_used,
                    SUM(prompt_tokens + response_tokens) as total_tokens
                FROM llm_model_usage
                WHERE model_name = ? AND provider = 'ollama'
            """, (model_name,))

            row = cursor.fetchone()
            if row and row[0] > 0:
                return {
                    'total_requests': row[0],
                    'avg_latency_ms': round(row[1] or 0, 2),
                    'last_used': row[2],
                    'total_tokens': row[3] or 0
                }

            return {
                'total_requests': 0,
                'avg_latency_ms': 0,
                'last_used': None,
                'total_tokens': 0
            }

    # ========================================================================
    # MODEL DISCOVERY
    # ========================================================================

    def list_models(self) -> List[Dict[str, Any]]:
        """
        List all available Ollama models

        Returns:
            List of model info dicts with name, size, digest, modified
        """
        try:
            response = requests.get(f"{self.ollama_url}/api/tags", timeout=5)
            response.raise_for_status()

            data = response.json()
            models = data.get('models', [])

            # Enrich with usage stats
            result = []
            for model in models:
                model_name = model.get('name', 'unknown')
                stats = self._get_model_usage_stats(model_name)

                result.append({
                    'name': model_name,
                    'size': model.get('size', 0),
                    'digest': model.get('digest', ''),
                    'modified': model.get('modified_at', ''),
                    **stats
                })

            return result

        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Failed to list models: {e}")
            return []

    # ========================================================================
    # MODEL LIFECYCLE
    # ========================================================================

    def load_model(self, model_name: str) -> bool:
        """
        Load a model into memory

        This is done by sending an empty generation request, which forces
        Ollama to load the model.

        Args:
            model_name: Name of the model to load (e.g., 'qwen2.5-coder:7b')

        Returns:
            True if successful, False otherwise
        """
        try:
            print(f"[INFO] Loading model: {model_name}")

            # Send empty generation request to load model
            response = requests.post(
                f"{self.ollama_url}/api/generate",
                json={
                    'model': model_name,
                    'prompt': '',
                    'stream': False
                },
                timeout=30
            )
            response.raise_for_status()

            # Track in memory
            self.models_loaded[model_name] = {
                'loaded_at': time.time(),
                'ram_mb': self._estimate_model_ram(model_name),
                'last_used': datetime.now().isoformat()
            }

            # Log event
            self._log_server_event(
                event_type='model_load',
                status='success',
                details=json.dumps({'model': model_name})
            )

            print(f"[SUCCESS] Model loaded: {model_name}")
            return True

        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Failed to load model {model_name}: {e}")
            self._log_server_event(
                event_type='model_load',
                status='error',
                details=json.dumps({'model': model_name, 'error': str(e)})
            )
            return False

    def unload_model(self, model_name: str) -> bool:
        """
        Unload a model from memory

        Note: Ollama doesn't have a direct unload API. This is a placeholder
        for future functionality or manual process management.

        Args:
            model_name: Name of the model to unload

        Returns:
            True if successful, False otherwise
        """
        try:
            print(f"[INFO] Unloading model: {model_name}")

            # Remove from tracking
            if model_name in self.models_loaded:
                del self.models_loaded[model_name]

            # Log event
            self._log_server_event(
                event_type='model_unload',
                status='success',
                details=json.dumps({'model': model_name})
            )

            print(f"[SUCCESS] Model unloaded: {model_name}")
            return True

        except Exception as e:
            print(f"[ERROR] Failed to unload model {model_name}: {e}")
            self._log_server_event(
                event_type='model_unload',
                status='error',
                details=json.dumps({'model': model_name, 'error': str(e)})
            )
            return False

    def _estimate_model_ram(self, model_name: str) -> int:
        """
        Estimate RAM usage for a model

        This is a rough estimate based on model size.

        Args:
            model_name: Name of the model

        Returns:
            Estimated RAM in MB
        """
        # Common model sizes (approximate)
        size_map = {
            '7b': 4096,
            '13b': 8192,
            '34b': 20480,
            '70b': 40960
        }

        for size_key, ram_mb in size_map.items():
            if size_key in model_name.lower():
                return ram_mb

        # Default estimate
        return 4096

    # ========================================================================
    # MODEL STATUS
    # ========================================================================

    def get_model_status(self, model_name: str) -> Dict[str, Any]:
        """
        Get status of a specific model

        Args:
            model_name: Name of the model

        Returns:
            Status dict with loaded, ram_mb, last_used, total_requests, avg_latency_ms
        """
        loaded = model_name in self.models_loaded
        stats = self._get_model_usage_stats(model_name)

        return {
            'name': model_name,
            'loaded': loaded,
            'ram_mb': self.models_loaded.get(model_name, {}).get('ram_mb', 0),
            'last_used': stats['last_used'],
            'total_requests': stats['total_requests'],
            'avg_latency_ms': stats['avg_latency_ms']
        }

    def get_all_status(self) -> Dict[str, Dict[str, Any]]:
        """
        Get status of all models

        Returns:
            Dict mapping model names to their status
        """
        models = self.list_models()
        result = {}

        for model in models:
            model_name = model['name']
            result[model_name] = self.get_model_status(model_name)

        return result

    # ========================================================================
    # SERVER MANAGEMENT
    # ========================================================================

    def is_server_running(self) -> bool:
        """Check if Ollama server is running"""
        try:
            response = requests.get(f"{self.ollama_url}/api/tags", timeout=2)
            return response.status_code == 200
        except requests.exceptions.RequestException:
            return False

    def get_server_pid(self) -> Optional[int]:
        """Get PID of running Ollama server process"""
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                # Look for ollama serve process
                cmdline = proc.info.get('cmdline', [])
                if cmdline and 'ollama' in ' '.join(cmdline).lower() and 'serve' in ' '.join(cmdline).lower():
                    return proc.info['pid']
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue

        return None

    def start_ollama_server(self) -> bool:
        """
        Start the Ollama server

        Returns:
            True if started successfully, False otherwise
        """
        try:
            # Check if already running
            if self.is_server_running():
                print("[INFO] Ollama server is already running")
                return True

            print("[INFO] Starting Ollama server...")

            # Start server process
            process = subprocess.Popen(
                ['ollama', 'serve'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                start_new_session=True
            )

            self.server_pid = process.pid
            self.server_start_time = time.time()

            # Wait for server to be ready (max 10 seconds)
            for _ in range(10):
                time.sleep(1)
                if self.is_server_running():
                    print(f"[SUCCESS] Ollama server started (PID: {self.server_pid})")

                    # Log event
                    self._log_server_event(
                        event_type='start',
                        status='success',
                        pid=self.server_pid
                    )

                    return True

            print("[ERROR] Ollama server failed to start within timeout")
            return False

        except FileNotFoundError:
            print("[ERROR] 'ollama' command not found. Is Ollama installed?")
            return False

        except Exception as e:
            print(f"[ERROR] Failed to start Ollama server: {e}")
            self._log_server_event(
                event_type='start',
                status='error',
                details=json.dumps({'error': str(e)})
            )
            return False

    def stop_ollama_server(self) -> bool:
        """
        Stop the Ollama server gracefully

        Returns:
            True if stopped successfully, False otherwise
        """
        try:
            pid = self.get_server_pid()

            if not pid:
                print("[INFO] No Ollama server process found")
                return True

            print(f"[INFO] Stopping Ollama server (PID: {pid})...")

            # Send SIGTERM for graceful shutdown
            process = psutil.Process(pid)
            process.terminate()

            # Wait up to 5 seconds for graceful shutdown
            try:
                process.wait(timeout=5)
            except psutil.TimeoutExpired:
                print("[WARNING] Server did not stop gracefully, forcing shutdown")
                process.kill()

            # Log event
            self._log_server_event(
                event_type='stop',
                status='success',
                pid=pid
            )

            self.server_pid = None
            self.server_start_time = None

            print("[SUCCESS] Ollama server stopped")
            return True

        except psutil.NoSuchProcess:
            print("[INFO] Ollama server process already stopped")
            return True

        except Exception as e:
            print(f"[ERROR] Failed to stop Ollama server: {e}")
            self._log_server_event(
                event_type='stop',
                status='error',
                details=json.dumps({'error': str(e)})
            )
            return False

    def get_server_status(self) -> Dict[str, Any]:
        """
        Get status of the Ollama server

        Returns:
            Status dict with running, pid, uptime_seconds, total_requests
        """
        running = self.is_server_running()
        pid = self.get_server_pid() if running else None

        uptime_seconds = 0
        if running and self.server_start_time:
            uptime_seconds = int(time.time() - self.server_start_time)

        # Get total requests from database
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT COUNT(*) FROM llm_model_usage WHERE provider = 'ollama'
            """)
            total_requests = cursor.fetchone()[0]

        return {
            'running': running,
            'pid': pid,
            'uptime_seconds': uptime_seconds,
            'total_requests': total_requests
        }


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for model manager"""
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Model Manager')
    parser.add_argument('action', choices=['list', 'load', 'unload', 'status', 'start-server', 'stop-server'],
                       help='Action to perform')
    parser.add_argument('model', nargs='?', help='Model name (for load/unload/status actions)')
    parser.add_argument('--url', default='http://localhost:11434', help='Ollama URL')
    parser.add_argument('--db', help='Database path')

    args = parser.parse_args()

    manager = OllamaModelManager(ollama_url=args.url, db_path=args.db)

    if args.action == 'list':
        models = manager.list_models()
        print(f"\nAvailable models ({len(models)}):\n")
        for model in models:
            print(f"  {model['name']}")
            print(f"    Size: {model['size'] / (1024**3):.1f} GB")
            print(f"    Requests: {model['total_requests']}")
            print(f"    Avg Latency: {model['avg_latency_ms']:.0f} ms")
            print()

    elif args.action == 'load':
        if not args.model:
            print("[ERROR] Model name required for load action")
            sys.exit(1)

        success = manager.load_model(args.model)
        sys.exit(0 if success else 1)

    elif args.action == 'unload':
        if not args.model:
            print("[ERROR] Model name required for unload action")
            sys.exit(1)

        success = manager.unload_model(args.model)
        sys.exit(0 if success else 1)

    elif args.action == 'status':
        if args.model:
            status = manager.get_model_status(args.model)
            print(f"\nModel: {status['name']}")
            print(f"Loaded: {status['loaded']}")
            print(f"RAM: {status['ram_mb']} MB")
            print(f"Last Used: {status['last_used'] or 'Never'}")
            print(f"Requests: {status['total_requests']}")
            print(f"Avg Latency: {status['avg_latency_ms']:.0f} ms")
        else:
            status = manager.get_server_status()
            print(f"\nOllama Server Status:")
            print(f"Running: {status['running']}")
            print(f"PID: {status['pid'] or 'N/A'}")
            print(f"Uptime: {status['uptime_seconds']} seconds")
            print(f"Total Requests: {status['total_requests']}")

    elif args.action == 'start-server':
        success = manager.start_ollama_server()
        sys.exit(0 if success else 1)

    elif args.action == 'stop-server':
        success = manager.stop_ollama_server()
        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()

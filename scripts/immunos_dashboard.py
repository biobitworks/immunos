#!/usr/bin/env python3
"""
IMMUNOS Dashboard - Main Application

A professional antivirus-style web interface for the IMMUNOS system.
Provides real-time monitoring, threat detection, and interactive controls
for all IMMUNOS cell types.

Usage:
    python immunos_dashboard.py [--port 5000] [--host 127.0.0.1]

Then open: http://localhost:5000
"""

import os
import sys
import sqlite3
import argparse
from pathlib import Path
from datetime import datetime
from flask import Flask, render_template, g
from flask_socketio import SocketIO, emit
from flask_cors import CORS

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

# Initialize Flask app
app = Flask(__name__,
            template_folder=str(Path(__file__).parent / 'templates'),
            static_folder=str(Path(__file__).parent / 'static'))

# Configuration
app.config['SECRET_KEY'] = 'immunos-dashboard-secret-key-change-in-production'
app.config['DATABASE'] = str(Path(__file__).parent.parent / '.immunos' / 'db' / 'dashboard.db')
app.config['BASE_PATH'] = str(Path(__file__).parent.parent)

# Enable CORS for development
CORS(app)

# Initialize SocketIO
socketio = SocketIO(app, cors_allowed_origins="*")


# Database Functions
# ==================

def get_db():
    """Get database connection for current request"""
    if 'db' not in g:
        g.db = sqlite3.connect(
            app.config['DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row
    return g.db


@app.teardown_appcontext
def close_db(error):
    """Close database connection after request"""
    db = g.pop('db', None)
    if db is not None:
        db.close()


def init_db():
    """Initialize the database with schema"""
    db_path = Path(app.config['DATABASE'])
    db_path.parent.mkdir(parents=True, exist_ok=True)

    schema_path = db_path.parent / 'schema.sql'

    # Check if schema exists
    if not schema_path.exists():
        print(f"❌ Error: Schema file not found at {schema_path}")
        print("Please ensure .immunos/db/schema.sql exists")
        sys.exit(1)

    # Create/update database
    conn = sqlite3.connect(str(db_path))
    with open(schema_path, 'r') as f:
        conn.executescript(f.read())
    conn.commit()
    conn.close()

    print(f"✓ Database initialized at {db_path}")


def get_setting(key: str, default=None):
    """Get a setting value from database"""
    db = get_db()
    row = db.execute('SELECT value, value_type FROM settings WHERE key = ?', (key,)).fetchone()

    if row is None:
        return default

    value, value_type = row['value'], row['value_type']

    # Convert based on type
    if value_type == 'bool':
        return value.lower() == 'true'
    elif value_type == 'int':
        return int(value)
    elif value_type == 'json':
        import json
        return json.loads(value)
    else:
        return value


def set_setting(key: str, value, value_type: str = 'string'):
    """Set a setting value in database"""
    import json

    # Convert value to string
    if value_type == 'bool':
        str_value = 'true' if value else 'false'
    elif value_type == 'json':
        str_value = json.dumps(value)
    else:
        str_value = str(value)

    db = get_db()
    db.execute('''
        INSERT INTO settings (key, value, value_type, updated_at)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(key) DO UPDATE SET
            value = excluded.value,
            value_type = excluded.value_type,
            updated_at = excluded.updated_at
    ''', (key, str_value, value_type, datetime.now()))
    db.commit()


def log_activity(event_type: str, event_action: str, title: str,
                 description: str = None, severity: str = 'info', metadata: dict = None):
    """Log an activity event"""
    import json

    db = get_db()
    db.execute('''
        INSERT INTO activity_log (event_type, event_action, title, description, severity, metadata_json)
        VALUES (?, ?, ?, ?, ?, ?)
    ''', (event_type, event_action, title, description, severity,
          json.dumps(metadata) if metadata else None))
    db.commit()

    # Emit to WebSocket clients
    socketio.emit('activity', {
        'type': event_type,
        'action': event_action,
        'title': title,
        'description': description,
        'severity': severity,
        'timestamp': datetime.now().isoformat()
    })


# Main Routes
# ===========

@app.route('/')
def index():
    """Main dashboard page"""
    return render_template('immunos_dashboard.html')


@app.route('/scans')
def scans():
    """Scans page (NK Cell focus)"""
    return render_template('scans.html')


@app.route('/memory')
def memory():
    """Memory browser page (T Cell focus)"""
    return render_template('memory.html')


@app.route('/llm')
def llm():
    """LLM interface page"""
    return render_template('llm.html')


@app.route('/todos')
def todos():
    """Redirect to todo dashboard"""
    from flask import redirect
    return redirect('http://localhost:5001', code=302)


# WebSocket Events
# ================

@socketio.on('connect')
def handle_connect():
    """Handle client connection"""
    print(f'✓ Client connected: {datetime.now().strftime("%H:%M:%S")}')
    emit('connected', {'message': 'Connected to IMMUNOS Dashboard'})


@socketio.on('disconnect')
def handle_disconnect():
    """Handle client disconnection"""
    print(f'✓ Client disconnected: {datetime.now().strftime("%H:%M:%S")}')


@socketio.on('subscribe')
def handle_subscribe(data):
    """Handle subscription to event types"""
    event_types = data.get('events', [])
    print(f'✓ Client subscribed to: {", ".join(event_types)}')
    emit('subscribed', {'events': event_types})


# Import API routes
# =================

try:
    from immunos_api import register_routes
    register_routes(app, socketio)
    print("✓ API routes registered")
except ImportError as e:
    print(f"⚠ Warning: Could not import API routes: {e}")
    print("  The dashboard will run but API endpoints will not be available")


# CLI Entry Point
# ===============

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='IMMUNOS Dashboard')
    parser.add_argument('--port', type=int, default=5000,
                       help='Port to run dashboard on (default: 5000)')
    parser.add_argument('--host', type=str, default='127.0.0.1',
                       help='Host to bind to (default: 127.0.0.1)')
    parser.add_argument('--init-db', action='store_true',
                       help='Initialize database and exit')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode')

    args = parser.parse_args()

    # Initialize database
    print("\n🧬 IMMUNOS Dashboard")
    print("=" * 60)

    if args.init_db:
        print("Initializing database...")
        init_db()
        print("✓ Database initialization complete")
        sys.exit(0)

    # Check if database exists, if not initialize it
    if not Path(app.config['DATABASE']).exists():
        print("Database not found, initializing...")
        init_db()

    print(f"Starting server on http://{args.host}:{args.port}")
    print("=" * 60)
    print()
    print(f"📊 Main Dashboard:  http://localhost:{args.port}")
    print(f"🔍 Scans & Threats: http://localhost:{args.port}/scans")
    print(f"🧠 Memory Browser:  http://localhost:{args.port}/memory")
    print(f"🤖 LLM Interface:   http://localhost:{args.port}/llm")
    print(f"✅ Todo System:     http://localhost:5001 (separate server)")
    print()
    print("Press Ctrl+C to stop the server")
    print()

    # Run app with SocketIO
    socketio.run(
        app,
        host=args.host,
        port=args.port,
        debug=args.debug,
        allow_unsafe_werkzeug=True  # For development only
    )

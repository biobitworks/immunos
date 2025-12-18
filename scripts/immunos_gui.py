#!/usr/bin/env python3
"""
IMMUNOS Web GUI

A web-based interface for the IMMUNOS todo system.

Usage:
    python immunos_gui.py [--port 5000]

Then open: http://localhost:5000
"""

from flask import Flask, render_template, request, jsonify, redirect, url_for
from pathlib import Path
import sys
import json
from datetime import datetime

# Import todo system
sys.path.insert(0, str(Path(__file__).parent))
from immunos_todo import TodoStore, Status, Priority

app = Flask(__name__,
            template_folder=str(Path(__file__).parent / 'templates'),
            static_folder=str(Path(__file__).parent / 'static'))

# Initialize todo store
store = TodoStore(Path('/Users/byron/projects'))


@app.route('/')
def index():
    """Main dashboard"""
    stats = store._calculate_stats()

    # Get todos by status
    inbox_todos = store.list_todos(status=Status.INBOX.value)
    next_todos = store.list_todos(status=Status.NEXT.value)
    waiting_todos = store.list_todos(status=Status.WAITING.value)
    someday_todos = store.list_todos(status=Status.SOMEDAY.value)

    # Get overdue
    overdue_todos = store.list_todos(overdue=True)

    # Get due today
    due_today = store.list_todos(due_today=True)

    return render_template('dashboard.html',
                         stats=stats,
                         inbox=inbox_todos,
                         next=next_todos,
                         waiting=waiting_todos,
                         someday=someday_todos,
                         overdue=overdue_todos,
                         due_today=due_today)


@app.route('/api/todos', methods=['GET'])
def get_todos():
    """Get todos (API endpoint)"""
    status = request.args.get('status')
    project = request.args.get('project')
    priority = request.args.get('priority')
    overdue = request.args.get('overdue') == 'true'

    todos = store.list_todos(
        status=status,
        project=project,
        priority=priority,
        overdue=overdue
    )

    return jsonify([{
        'id': t.id,
        'title': t.title,
        'status': t.status,
        'priority': t.priority,
        'due': t.due,
        'project': t.project,
        'tags': t.tags,
        'description': t.description,
        'created': t.created,
        'completed': t.completed
    } for t in todos])


@app.route('/api/todos', methods=['POST'])
def create_todo():
    """Create a new todo"""
    data = request.json

    todo = store.add_todo(
        title=data['title'],
        status=data.get('status', Status.INBOX.value),
        due=data.get('due'),
        priority=data.get('priority', Priority.MEDIUM.value),
        project=data.get('project'),
        tags=data.get('tags', []),
        description=data.get('description')
    )

    return jsonify({'success': True, 'todo_id': todo.id})


@app.route('/api/todos/<todo_id>', methods=['PUT'])
def update_todo(todo_id):
    """Update a todo"""
    data = request.json

    # Update fields
    kwargs = {}
    if 'title' in data:
        kwargs['title'] = data['title']
    if 'status' in data:
        kwargs['status'] = data['status']
    if 'priority' in data:
        kwargs['priority'] = data['priority']
    if 'due' in data:
        kwargs['due'] = data['due']
    if 'project' in data:
        kwargs['project'] = data['project']
    if 'tags' in data:
        kwargs['tags'] = data['tags']
    if 'description' in data:
        kwargs['description'] = data['description']

    todo = store.update_todo(todo_id, **kwargs)

    if not todo:
        return jsonify({'success': False, 'error': 'Todo not found'}), 404

    return jsonify({'success': True})


@app.route('/api/todos/<todo_id>/move', methods=['POST'])
def move_todo(todo_id):
    """Move todo to new status"""
    data = request.json
    new_status = data.get('status')

    todo = store.move_todo(todo_id, new_status)

    if not todo:
        return jsonify({'success': False, 'error': 'Todo not found'}), 404

    return jsonify({'success': True})


@app.route('/api/todos/<todo_id>/complete', methods=['POST'])
def complete_todo(todo_id):
    """Mark todo as completed"""
    todo = store.complete_todo(todo_id)

    if not todo:
        return jsonify({'success': False, 'error': 'Todo not found'}), 404

    return jsonify({'success': True})


@app.route('/api/stats', methods=['GET'])
def get_stats():
    """Get statistics"""
    stats = store._calculate_stats()
    return jsonify(stats)


def create_templates():
    """Create HTML templates"""
    template_dir = Path(__file__).parent / 'templates'
    template_dir.mkdir(exist_ok=True)

    # Create dashboard template
    dashboard_html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IMMUNOS - GTD Todo System</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            background: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }

        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }

        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }

        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }

        header p {
            opacity: 0.9;
            font-size: 1.1em;
        }

        .stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }

        .stat-card {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            border-left: 4px solid #667eea;
        }

        .stat-card h3 {
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
            margin-bottom: 10px;
        }

        .stat-card .value {
            font-size: 2.5em;
            font-weight: bold;
            color: #667eea;
        }

        .stat-card.danger {
            border-left-color: #e74c3c;
        }

        .stat-card.danger .value {
            color: #e74c3c;
        }

        .stat-card.warning {
            border-left-color: #f39c12;
        }

        .stat-card.warning .value {
            color: #f39c12;
        }

        .stat-card.success {
            border-left-color: #27ae60;
        }

        .stat-card.success .value {
            color: #27ae60;
        }

        .sections {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
        }

        .section {
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow: hidden;
        }

        .section-header {
            background: #f8f9fa;
            padding: 15px 20px;
            border-bottom: 2px solid #e9ecef;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .section-header h2 {
            font-size: 1.3em;
            color: #2c3e50;
        }

        .section-header .badge {
            background: #667eea;
            color: white;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.9em;
            font-weight: bold;
        }

        .section-header.inbox .badge { background: #95a5a6; }
        .section-header.next .badge { background: #3498db; }
        .section-header.waiting .badge { background: #f39c12; }
        .section-header.someday .badge { background: #9b59b6; }
        .section-header.overdue .badge { background: #e74c3c; }

        .todo-list {
            padding: 10px;
            max-height: 500px;
            overflow-y: auto;
        }

        .todo-item {
            padding: 15px;
            border-left: 3px solid #e9ecef;
            margin-bottom: 10px;
            background: #f8f9fa;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.2s;
        }

        .todo-item:hover {
            background: #e9ecef;
            border-left-color: #667eea;
            transform: translateX(5px);
        }

        .todo-item.high {
            border-left-color: #e74c3c;
        }

        .todo-item.medium {
            border-left-color: #f39c12;
        }

        .todo-item.low {
            border-left-color: #27ae60;
        }

        .todo-title {
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 5px;
        }

        .todo-meta {
            font-size: 0.85em;
            color: #7f8c8d;
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
        }

        .todo-meta span {
            display: inline-flex;
            align-items: center;
            gap: 5px;
        }

        .priority-badge {
            padding: 2px 8px;
            border-radius: 10px;
            font-size: 0.75em;
            font-weight: bold;
            text-transform: uppercase;
        }

        .priority-high {
            background: #fee;
            color: #e74c3c;
        }

        .priority-medium {
            background: #fef5e7;
            color: #f39c12;
        }

        .priority-low {
            background: #eafaf1;
            color: #27ae60;
        }

        .empty-state {
            text-align: center;
            padding: 40px 20px;
            color: #95a5a6;
        }

        .empty-state i {
            font-size: 3em;
            margin-bottom: 10px;
            display: block;
        }

        .btn {
            display: inline-block;
            padding: 10px 20px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-size: 1em;
            transition: all 0.2s;
            text-decoration: none;
        }

        .btn:hover {
            background: #5568d3;
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }

        .add-todo-btn {
            position: fixed;
            bottom: 30px;
            right: 30px;
            width: 60px;
            height: 60px;
            border-radius: 50%;
            background: #667eea;
            color: white;
            font-size: 2em;
            border: none;
            cursor: pointer;
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
            transition: all 0.2s;
        }

        .add-todo-btn:hover {
            background: #5568d3;
            transform: scale(1.1);
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 IMMUNOS</h1>
            <p>Getting Things Done - Todo System</p>
        </header>

        <div class="stats">
            <div class="stat-card">
                <h3>Total Active</h3>
                <div class="value">{{ stats.total - stats.by_status.get('done', 0) }}</div>
            </div>
            <div class="stat-card danger">
                <h3>Overdue</h3>
                <div class="value">{{ stats.overdue }}</div>
            </div>
            <div class="stat-card warning">
                <h3>Due Today</h3>
                <div class="value">{{ stats.due_today }}</div>
            </div>
            <div class="stat-card success">
                <h3>Due This Week</h3>
                <div class="value">{{ stats.due_this_week }}</div>
            </div>
        </div>

        {% if overdue %}
        <div class="section">
            <div class="section-header overdue">
                <h2>🔴 Overdue</h2>
                <span class="badge">{{ overdue|length }}</span>
            </div>
            <div class="todo-list">
                {% for todo in overdue %}
                <div class="todo-item {{ todo.priority }}">
                    <div class="todo-title">{{ todo.title }}</div>
                    <div class="todo-meta">
                        <span class="priority-badge priority-{{ todo.priority }}">{{ todo.priority }}</span>
                        {% if todo.project %}<span>📁 {{ todo.project }}</span>{% endif %}
                        {% if todo.due %}<span>📅 {{ todo.due[:10] }}</span>{% endif %}
                    </div>
                </div>
                {% endfor %}
            </div>
        </div>
        {% endif %}

        <div class="sections">
            <div class="section">
                <div class="section-header inbox">
                    <h2>📥 Inbox</h2>
                    <span class="badge">{{ inbox|length }}</span>
                </div>
                <div class="todo-list">
                    {% if inbox %}
                        {% for todo in inbox %}
                        <div class="todo-item {{ todo.priority }}">
                            <div class="todo-title">{{ todo.title }}</div>
                            <div class="todo-meta">
                                <span class="priority-badge priority-{{ todo.priority }}">{{ todo.priority }}</span>
                                {% if todo.project %}<span>📁 {{ todo.project }}</span>{% endif %}
                                {% if todo.due %}<span>📅 {{ todo.due[:10] }}</span>{% endif %}
                            </div>
                        </div>
                        {% endfor %}
                    {% else %}
                        <div class="empty-state">
                            <i>✓</i>
                            <div>Inbox is empty!</div>
                        </div>
                    {% endif %}
                </div>
            </div>

            <div class="section">
                <div class="section-header next">
                    <h2>⚡ Next Actions</h2>
                    <span class="badge">{{ next|length }}</span>
                </div>
                <div class="todo-list">
                    {% if next %}
                        {% for todo in next %}
                        <div class="todo-item {{ todo.priority }}">
                            <div class="todo-title">{{ todo.title }}</div>
                            <div class="todo-meta">
                                <span class="priority-badge priority-{{ todo.priority }}">{{ todo.priority }}</span>
                                {% if todo.project %}<span>📁 {{ todo.project }}</span>{% endif %}
                                {% if todo.due %}<span>📅 {{ todo.due[:10] }}</span>{% endif %}
                            </div>
                        </div>
                        {% endfor %}
                    {% else %}
                        <div class="empty-state">
                            <i>✓</i>
                            <div>No next actions</div>
                        </div>
                    {% endif %}
                </div>
            </div>

            <div class="section">
                <div class="section-header waiting">
                    <h2>⏸️ Waiting</h2>
                    <span class="badge">{{ waiting|length }}</span>
                </div>
                <div class="todo-list">
                    {% if waiting %}
                        {% for todo in waiting %}
                        <div class="todo-item {{ todo.priority }}">
                            <div class="todo-title">{{ todo.title }}</div>
                            <div class="todo-meta">
                                <span class="priority-badge priority-{{ todo.priority }}">{{ todo.priority }}</span>
                                {% if todo.project %}<span>📁 {{ todo.project }}</span>{% endif %}
                                {% if todo.due %}<span>📅 {{ todo.due[:10] }}</span>{% endif %}
                            </div>
                        </div>
                        {% endfor %}
                    {% else %}
                        <div class="empty-state">
                            <i>✓</i>
                            <div>Nothing waiting</div>
                        </div>
                    {% endif %}
                </div>
            </div>

            <div class="section">
                <div class="section-header someday">
                    <h2>💭 Someday/Maybe</h2>
                    <span class="badge">{{ someday|length }}</span>
                </div>
                <div class="todo-list">
                    {% if someday %}
                        {% for todo in someday %}
                        <div class="todo-item {{ todo.priority }}">
                            <div class="todo-title">{{ todo.title }}</div>
                            <div class="todo-meta">
                                <span class="priority-badge priority-{{ todo.priority }}">{{ todo.priority }}</span>
                                {% if todo.project %}<span>📁 {{ todo.project }}</span>{% endif %}
                            </div>
                        </div>
                        {% endfor %}
                    {% else %}
                        <div class="empty-state">
                            <i>✓</i>
                            <div>No future plans</div>
                        </div>
                    {% endif %}
                </div>
            </div>
        </div>

        <button class="add-todo-btn" onclick="alert('Add todo functionality coming soon!')">+</button>
    </div>
</body>
</html>'''

    with open(template_dir / 'dashboard.html', 'w') as f:
        f.write(dashboard_html)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Web GUI')
    parser.add_argument('--port', type=int, default=5000, help='Port to run on')
    parser.add_argument('--host', type=str, default='127.0.0.1', help='Host to bind to')

    args = parser.parse_args()

    # Create templates
    print("Creating templates...")
    create_templates()

    print(f"\n🧬 IMMUNOS Web GUI")
    print(f"{'='*60}")
    print(f"Starting server on http://{args.host}:{args.port}")
    print(f"{'='*60}\n")
    print(f"Open your browser to: http://localhost:{args.port}")
    print(f"\nPress Ctrl+C to stop the server\n")

    app.run(host=args.host, port=args.port, debug=True)

---
source: /Users/byron/projects/scripts/templates/dashboard.html
relative: scripts/templates/dashboard.html
generated_at: 2025-12-23 10:28
---

```html
<!DOCTYPE html>
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
</html>
```

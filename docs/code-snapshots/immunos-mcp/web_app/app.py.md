---
source: /Users/byron/projects/immunos-mcp/web_app/app.py
relative: immunos-mcp/web_app/app.py
generated_at: 2025-12-23 10:28
---

```python
"""
IMMUNOS-MCP Portfolio Demo - Flask Web Application
Multi-domain AI pattern recognition and anomaly detection
"""

import json
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from database import Pattern, Analysis, UserSubmission, Metric, get_stats

# Initialize Flask app
app = Flask(__name__)
app.config['SECRET_KEY'] = 'immunos-mcp-demo-secret-key-2025'
CORS(app)

# Database setup
DB_PATH = 'data/immunos.db'
engine = create_engine(f'sqlite:///{DB_PATH}')
Session = sessionmaker(bind=engine)


# ============================================================================
# HOMEPAGE & NAVIGATION
# ============================================================================

@app.route('/')
def index():
    """Landing page with project overview"""
    session = Session()
    stats = get_stats(session)
    session.close()

    return render_template('index.html', stats=stats)


@app.route('/algorithms')
def algorithms():
    """Algorithm visualization page"""
    return render_template('algorithms.html')


@app.route('/training')
def training():
    """Live training interface"""
    session = Session()
    stats = get_stats(session)
    session.close()

    return render_template('training.html', stats=stats)


# ============================================================================
# DEMO PAGES
# ============================================================================

@app.route('/security')
def security():
    """Code security scanner demo"""
    return render_template('security.html')


@app.route('/emotion')
def emotion():
    """Emotion detection demo"""
    return render_template('emotion.html')


@app.route('/spam')
def spam():
    """Email/spam classification demo"""
    return render_template('spam.html')


@app.route('/network')
def network():
    """Network intrusion detection demo"""
    return render_template('network.html')


# ============================================================================
# API ENDPOINTS
# ============================================================================

@app.route('/api/stats')
def api_stats():
    """Get database statistics"""
    session = Session()
    stats = get_stats(session)
    session.close()

    return jsonify(stats)


@app.route('/api/scan', methods=['POST'])
def api_scan():
    """
    Code security scanning endpoint

    Request:
        {
            "code": "def execute_command(cmd):\n    os.system(cmd)",
            "pipeline": "standard"  # quick|standard|deep
        }

    Response:
        {
            "predicted_class": "vulnerable",
            "confidence": 0.95,
            "cwe_types": ["CWE-78"],
            "agents_used": ["bcell", "nkcell"],
            "execution_time": 0.124
        }
    """
    data = request.json
    code = data.get('code', '')
    pipeline = data.get('pipeline', 'standard')

    # TODO: Implement actual scanning with IMMUNOS agents
    # For now, return mock response
    result = {
        'predicted_class': 'safe',
        'confidence': 0.85,
        'cwe_types': [],
        'agents_used': ['bcell'],
        'execution_time': 0.05
    }

    # Store analysis
    session = Session()
    analysis = Analysis(
        domain='security',
        input_data=code,
        predicted_label=result['predicted_class'],
        confidence=result['confidence'],
        agent_used=','.join(result['agents_used']),
        execution_time=result['execution_time']
    )
    session.add(analysis)
    session.commit()
    session.close()

    return jsonify(result)


@app.route('/api/detect_emotion', methods=['POST'])
def api_detect_emotion():
    """
    Emotion detection endpoint

    Request:
        {
            "text": "I am so happy today!",
            "modality": "text"  # text|facial
        }

    Response:
        {
            "dominant_emotion": "joy",
            "probabilities": {"joy": 0.92, "sadness": 0.03, ...},
            "confidence": 0.92
        }
    """
    data = request.json
    text = data.get('text', '')
    modality = data.get('modality', 'text')

    # TODO: Implement emotion detection
    result = {
        'dominant_emotion': 'neutral',
        'probabilities': {
            'joy': 0.15,
            'sadness': 0.10,
            'anger': 0.05,
            'fear': 0.05,
            'neutral': 0.65
        },
        'confidence': 0.65
    }

    # Store analysis
    session = Session()
    analysis = Analysis(
        domain='emotion',
        input_data=text,
        predicted_label=result['dominant_emotion'],
        confidence=result['confidence'],
        agent_used='bcell',
        execution_time=0.03
    )
    session.add(analysis)
    session.commit()
    session.close()

    return jsonify(result)


@app.route('/api/classify_email', methods=['POST'])
def api_classify_email():
    """
    Email classification endpoint

    Request:
        {
            "subject": "Meeting tomorrow",
            "body": "Hi, let's meet at 2pm"
        }

    Response:
        {
            "classification": "ham",
            "confidence": 0.94,
            "indicators": ["professional language", "specific time"]
        }
    """
    data = request.json
    subject = data.get('subject', '')
    body = data.get('body', '')
    email_text = f"{subject}\n{body}"

    # TODO: Implement spam classification
    result = {
        'classification': 'ham',
        'confidence': 0.90,
        'indicators': ['professional language']
    }

    # Store analysis
    session = Session()
    analysis = Analysis(
        domain='spam',
        input_data=email_text,
        predicted_label=result['classification'],
        confidence=result['confidence'],
        agent_used='bcell',
        execution_time=0.02
    )
    session.add(analysis)
    session.commit()
    session.close()

    return jsonify(result)


@app.route('/api/analyze_traffic', methods=['POST'])
def api_analyze_traffic():
    """
    Network traffic analysis endpoint

    Request:
        {
            "traffic_data": {...},  # Network flow features
            "format": "json"  # json|csv
        }

    Response:
        {
            "classification": "normal",
            "attack_type": null,
            "confidence": 0.88,
            "anomaly_score": 0.12
        }
    """
    data = request.json
    traffic_data = data.get('traffic_data', {})

    # TODO: Implement network intrusion detection
    result = {
        'classification': 'normal',
        'attack_type': None,
        'confidence': 0.88,
        'anomaly_score': 0.12
    }

    # Store analysis
    session = Session()
    analysis = Analysis(
        domain='network',
        input_data=json.dumps(traffic_data),
        predicted_label=result['classification'],
        confidence=result['confidence'],
        agent_used='nkcell',
        execution_time=0.08
    )
    session.add(analysis)
    session.commit()
    session.close()

    return jsonify(result)


@app.route('/api/submit_training', methods=['POST'])
def api_submit_training():
    """
    Submit new training example

    Request:
        {
            "domain": "security",
            "label": "safe",
            "data": "def add(a, b): return a + b"
        }

    Response:
        {
            "success": true,
            "message": "Training example submitted",
            "pattern_id": 123
        }
    """
    data = request.json
    domain = data.get('domain')
    label = data.get('label')
    input_data = data.get('data')

    # Validate inputs
    if not all([domain, label, input_data]):
        return jsonify({
            'success': False,
            'message': 'Missing required fields'
        }), 400

    # Store submission
    session = Session()
    submission = UserSubmission(
        domain=domain,
        data=input_data,
        user_label=label,
        accepted=False  # Requires admin approval
    )
    session.add(submission)
    session.commit()

    submission_id = submission.id
    session.close()

    return jsonify({
        'success': True,
        'message': 'Training example submitted successfully',
        'submission_id': submission_id
    })


# ============================================================================
# ERROR HANDLERS
# ============================================================================

@app.errorhandler(404)
def not_found(error):
    """Handle 404 errors"""
    return render_template('404.html'), 404


@app.errorhandler(500)
def server_error(error):
    """Handle 500 errors"""
    return render_template('500.html'), 500


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("🧬 IMMUNOS-MCP Portfolio Demo")
    print("=" * 60)
    print("\n✅ Flask server starting...")
    print(f"📁 Database: {DB_PATH}")
    print(f"🌐 Access at: http://localhost:5000")
    print("\nPress Ctrl+C to stop the server\n")

    app.run(
        host='0.0.0.0',
        port=5000,
        debug=True
    )

```

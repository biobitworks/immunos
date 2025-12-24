---
source: /Users/byron/projects/immunos-mcp/web_app/templates/index.html
relative: immunos-mcp/web_app/templates/index.html
generated_at: 2025-12-23 10:28
---

```html
{% extends "base.html" %}

{% block title %}IMMUNOS-MCP - Artificial Immune System for Multi-Domain Pattern Recognition{% endblock %}

{% block content %}
<!-- Hero Section -->
<div class="jumbotron bg-gradient text-white p-5 rounded-3 mb-5" style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);">
    <div class="container">
        <h1 class="display-4 fw-bold">
            <i class="fas fa-shield-virus"></i> IMMUNOS-MCP
        </h1>
        <p class="lead">Artificial Immune System for Multi-Domain Pattern Recognition & Anomaly Detection</p>
        <hr class="my-4" style="border-color: rgba(255,255,255,0.3);">
        <p>Inspired by biological immune systems | Powered by QML-AiNet & Opt-AiNet algorithms | 4 domains demonstrated</p>
        <div class="mt-4">
            <a class="btn btn-light btn-lg me-2" href="/algorithms" role="button">
                <i class="fas fa-chart-line"></i> View Algorithms
            </a>
            <a class="btn btn-outline-light btn-lg" href="/training" role="button">
                <i class="fas fa-graduation-cap"></i> Try Live Training
            </a>
        </div>
    </div>
</div>

<!-- Quick Stats -->
<div class="row mb-5">
    <div class="col-md-3">
        <div class="card text-center">
            <div class="card-body">
                <i class="fas fa-database fa-3x text-primary mb-3"></i>
                <h3 class="card-title">{{ stats.total_patterns }}</h3>
                <p class="card-text text-muted">Training Patterns</p>
            </div>
        </div>
    </div>
    <div class="col-md-3">
        <div class="card text-center">
            <div class="card-body">
                <i class="fas fa-search fa-3x text-success mb-3"></i>
                <h3 class="card-title">{{ stats.total_analyses }}</h3>
                <p class="card-text text-muted">Analyses Performed</p>
            </div>
        </div>
    </div>
    <div class="col-md-3">
        <div class="card text-center">
            <div class="card-body">
                <i class="fas fa-paper-plane fa-3x text-info mb-3"></i>
                <h3 class="card-title">{{ stats.total_submissions }}</h3>
                <p class="card-text text-muted">User Submissions</p>
            </div>
        </div>
    </div>
    <div class="col-md-3">
        <div class="card text-center">
            <div class="card-body">
                <i class="fas fa-flask fa-3x text-warning mb-3"></i>
                <h3 class="card-title">4</h3>
                <p class="card-text text-muted">Use Cases</p>
            </div>
        </div>
    </div>
</div>

<!-- Use Cases -->
<h2 class="mb-4"><i class="fas fa-microscope"></i> Multi-Domain Demonstrations</h2>
<div class="row mb-5">
    <!-- Code Security -->
    <div class="col-md-6 mb-4">
        <div class="card h-100 shadow-sm hover-card">
            <div class="card-body">
                <div class="d-flex align-items-center mb-3">
                    <i class="fas fa-code fa-2x text-danger me-3"></i>
                    <h4 class="card-title mb-0">Code Security Scanner</h4>
                </div>
                <p class="card-text">
                    Detect vulnerabilities in source code using B Cell pattern matching.
                    Identifies CWE patterns including SQL injection, command injection, and XSS.
                </p>
                <ul class="list-unstyled">
                    <li><i class="fas fa-check text-success"></i> Real-time code analysis</li>
                    <li><i class="fas fa-check text-success"></i> CWE classification</li>
                    <li><i class="fas fa-check text-success"></i> {{ stats.patterns_by_domain.security }} training patterns</li>
                </ul>
                <a href="/security" class="btn btn-danger">
                    <i class="fas fa-play"></i> Launch Demo
                </a>
            </div>
        </div>
    </div>

    <!-- Emotion Detection -->
    <div class="col-md-6 mb-4">
        <div class="card h-100 shadow-sm hover-card">
            <div class="card-body">
                <div class="d-flex align-items-center mb-3">
                    <i class="fas fa-smile fa-2x text-primary me-3"></i>
                    <h4 class="card-title mb-0">Emotion Detection</h4>
                </div>
                <p class="card-text">
                    Analyze text and facial expressions for emotional content.
                    Multi-class classification with confidence scores.
                </p>
                <ul class="list-unstyled">
                    <li><i class="fas fa-check text-success"></i> Text-based analysis</li>
                    <li><i class="fas fa-check text-success"></i> 7 emotion categories</li>
                    <li><i class="fas fa-check text-success"></i> {{ stats.patterns_by_domain.emotion }} training patterns</li>
                </ul>
                <a href="/emotion" class="btn btn-primary">
                    <i class="fas fa-play"></i> Launch Demo
                </a>
            </div>
        </div>
    </div>

    <!-- Email Classification -->
    <div class="col-md-6 mb-4">
        <div class="card h-100 shadow-sm hover-card">
            <div class="card-body">
                <div class="d-flex align-items-center mb-3">
                    <i class="fas fa-envelope fa-2x text-success me-3"></i>
                    <h4 class="card-title mb-0">Email/Phishing Classifier</h4>
                </div>
                <p class="card-text">
                    Classify emails as spam or legitimate using pattern recognition.
                    Identifies phishing attempts and suspicious content.
                </p>
                <ul class="list-unstyled">
                    <li><i class="fas fa-check text-success"></i> Spam detection</li>
                    <li><i class="fas fa-check text-success"></i> Phishing identification</li>
                    <li><i class="fas fa-check text-success"></i> {{ stats.patterns_by_domain.spam }} training patterns</li>
                </ul>
                <a href="/spam" class="btn btn-success">
                    <i class="fas fa-play"></i> Launch Demo
                </a>
            </div>
        </div>
    </div>

    <!-- Network Intrusion -->
    <div class="col-md-6 mb-4">
        <div class="card h-100 shadow-sm hover-card">
            <div class="card-body">
                <div class="d-flex align-items-center mb-3">
                    <i class="fas fa-network-wired fa-2x text-warning me-3"></i>
                    <h4 class="card-title mb-0">Network Intrusion Detection</h4>
                </div>
                <p class="card-text">
                    Detect network attacks using NK Cell negative selection.
                    Identifies DoS, Probe, R2L, and U2R attacks.
                </p>
                <ul class="list-unstyled">
                    <li><i class="fas fa-check text-success"></i> Anomaly detection</li>
                    <li><i class="fas fa-check text-success"></i> Attack classification</li>
                    <li><i class="fas fa-check text-success"></i> {{ stats.patterns_by_domain.network }} training patterns</li>
                </ul>
                <a href="/network" class="btn btn-warning">
                    <i class="fas fa-play"></i> Launch Demo
                </a>
            </div>
        </div>
    </div>
</div>

<!-- System Architecture -->
<h2 class="mb-4"><i class="fas fa-sitemap"></i> System Architecture</h2>
<div class="row mb-5">
    <div class="col-md-8 offset-md-2">
        <div class="card">
            <div class="card-body text-center">
                <h5 class="card-title">Multi-Agent Immune System</h5>
                <div class="architecture-diagram mt-4 mb-4">
                    <div class="row">
                        <div class="col-md-4">
                            <div class="agent-box bg-primary text-white p-3 rounded">
                                <i class="fas fa-user-shield fa-2x mb-2"></i>
                                <h6>B Cell Agent</h6>
                                <small>Pattern Matching</small>
                            </div>
                        </div>
                        <div class="col-md-4">
                            <div class="agent-box bg-danger text-white p-3 rounded">
                                <i class="fas fa-exclamation-triangle fa-2x mb-2"></i>
                                <h6>NK Cell Agent</h6>
                                <small>Anomaly Detection</small>
                            </div>
                        </div>
                        <div class="col-md-4">
                            <div class="agent-box bg-success text-white p-3 rounded">
                                <i class="fas fa-network-wired fa-2x mb-2"></i>
                                <h6>Orchestrator</h6>
                                <small>Multi-Domain Coordination</small>
                            </div>
                        </div>
                    </div>
                </div>
                <p class="text-muted">
                    Inspired by the human immune system's ability to recognize patterns and detect anomalies
                </p>
            </div>
        </div>
    </div>
</div>

<!-- Key Features -->
<h2 class="mb-4"><i class="fas fa-star"></i> Key Features</h2>
<div class="row mb-5">
    <div class="col-md-4">
        <div class="feature-box">
            <i class="fas fa-brain fa-3x text-primary mb-3"></i>
            <h5>Adaptive Learning</h5>
            <p>Self-adjusting pattern recognition based on QML-AiNet and Opt-AiNet algorithms</p>
        </div>
    </div>
    <div class="col-md-4">
        <div class="feature-box">
            <i class="fas fa-layer-group fa-3x text-success mb-3"></i>
            <h5>Multi-Domain</h5>
            <p>Single system handles security, emotion, spam, and network analysis</p>
        </div>
    </div>
    <div class="col-md-4">
        <div class="feature-box">
            <i class="fas fa-tachometer-alt fa-3x text-danger mb-3"></i>
            <h5>Real-Time Performance</h5>
            <p>Sub-second response times for all classification tasks</p>
        </div>
    </div>
</div>

<!-- Call to Action -->
<div class="row">
    <div class="col-md-12">
        <div class="card bg-light">
            <div class="card-body text-center p-5">
                <h3><i class="fas fa-rocket"></i> Ready to Explore?</h3>
                <p class="lead">Try the interactive demos or train the system with your own data</p>
                <div class="mt-4">
                    <a href="/security" class="btn btn-danger btn-lg me-2">
                        <i class="fas fa-code"></i> Code Scanner
                    </a>
                    <a href="/emotion" class="btn btn-primary btn-lg me-2">
                        <i class="fas fa-smile"></i> Emotion Detector
                    </a>
                    <a href="/training" class="btn btn-success btn-lg">
                        <i class="fas fa-graduation-cap"></i> Live Training
                    </a>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}

{% block extra_scripts %}
<style>
    .hover-card {
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .hover-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 8px 16px rgba(0,0,0,0.1);
    }
    .feature-box {
        text-align: center;
        padding: 20px;
    }
    .agent-box {
        margin: 10px 0;
    }
    .bg-gradient {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    }
</style>
{% endblock %}

```

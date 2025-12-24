// ============================================================================
// IMMUNOS Real-Time Monitor - Main JavaScript
// ============================================================================

const escapeHtml = (value) => {
    const text = value === null || value === undefined ? '' : String(value);
    return text.replace(/[&<>"']/g, (char) => ({
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
    })[char]);
};

const Monitor = {
    // WebSocket connection
    socket: null,
    connected: false,

    // State
    currentTab: 'overview',
    currentDomain: 'emotion',
    models: {},
    trainingJobs: {},
    sessionId: null,

    // Charts
    activityChart: null,
    tokensChart: null,

    // Event history for activity chart (last 20 minutes)
    eventHistory: Array(20).fill(0),
    timeLabels: [],

    // Domain information
    domainInfo: {
        emotion: {
            name: 'Emotion Recognition',
            description: 'Detects anomalous emotional expressions in facial images using CNN features, HOG descriptors, and LBP patterns. Trained on CK+ and FER2013 datasets.'
        },
        hallucination: {
            name: 'Hallucination Detection',
            description: 'Identifies factual inconsistencies and hallucinations in LLM outputs using BERT embeddings, linguistic markers, and entailment scoring. Uses TruthfulQA and HaluEval datasets.'
        },
        network: {
            name: 'Network Intrusion',
            description: 'Detects network intrusions and anomalous traffic patterns using NSL-KDD features (41 dimensions). Identifies DoS, Probe, R2L, and U2R attacks.'
        },
        code: {
            name: 'Code Vulnerability',
            description: 'Scans source code for security vulnerabilities using AST analysis, dataflow patterns, and known vulnerability signatures. Trained on DiverseVul dataset.'
        },
        research: {
            name: 'Research Verification',
            description: 'Verifies scientific claims against evidence using SciBERT embeddings, citation analysis, and entailment models. Uses SciFact dataset.'
        }
    },

    // ========================================================================
    // INITIALIZATION
    // ========================================================================

    init() {
        console.log('[IMMUNOS] Initializing monitor...');
        this.setupTabs();
        this.setupDomainButtons();
        this.setupEventListeners();
        this.setupCharts();
        this.connectWebSocket();
        this.loadInitialData();
        this.startPeriodicUpdates();
        this.initTimeLabels();
    },

    initTimeLabels() {
        const now = new Date();
        for (let i = 19; i >= 0; i--) {
            const time = new Date(now.getTime() - i * 60000);
            this.timeLabels.push(time.toLocaleTimeString('en-US', { hour: '2-digit', minute: '2-digit' }));
        }
    },

    // ========================================================================
    // WEBSOCKET
    // ========================================================================

    connectWebSocket() {
        console.log('[IMMUNOS] Connecting to WebSocket...');
        this.socket = io();

        this.socket.on('connect', () => {
            this.connected = true;
            this.updateConnectionStatus();
            console.log('[IMMUNOS] WebSocket connected');
        });

        this.socket.on('disconnect', () => {
            this.connected = false;
            this.updateConnectionStatus();
            console.log('[IMMUNOS] WebSocket disconnected');
        });

        // Training events
        this.socket.on('training_started', (data) => this.handleTrainingStarted(data));
        this.socket.on('training_progress', (data) => this.handleTrainingProgress(data));
        this.socket.on('training_complete', (data) => this.handleTrainingComplete(data));

        // Model events
        this.socket.on('model_switched', (data) => this.handleModelSwitched(data));
        this.socket.on('model_loaded', (data) => this.handleModelLoaded(data));
        this.socket.on('model_unloaded', (data) => this.handleModelUnloaded(data));

        // Detection events
        this.socket.on('detection_event', (data) => this.handleDetectionEvent(data));

        // Resource events
        this.socket.on('resource_update', (data) => this.handleResourceUpdate(data));

        // Server events
        this.socket.on('ollama_server_started', (data) => this.handleOllamaServerStarted(data));
        this.socket.on('ollama_server_stopped', (data) => this.handleOllamaStopped(data));
    },

    updateConnectionStatus() {
        const statusDiv = document.getElementById('connection-status');
        if (this.connected) {
            statusDiv.innerHTML = `
                <div class="w-3 h-3 bg-green-500 rounded-full animate-pulse"></div>
                <span class="text-sm text-green-400">Connected</span>
            `;
        } else {
            statusDiv.innerHTML = `
                <div class="w-3 h-3 bg-red-500 rounded-full"></div>
                <span class="text-sm text-red-400">Disconnected</span>
            `;
        }
    },

    // ========================================================================
    // INITIAL DATA LOAD
    // ========================================================================

    async loadInitialData() {
        try {
            await Promise.all([
                this.loadModelsStatus(),
                this.loadTokensComparison(),
                this.loadSessionStatus(),
                this.loadRoutingStatus()
            ]);
            console.log('[IMMUNOS] Initial data loaded');
        } catch (error) {
            console.error('[IMMUNOS] Error loading initial data:', error);
        }
    },

    async loadModelsStatus() {
        try {
            const response = await fetch('/api/models/status');
            if (!response.ok) throw new Error('Failed to load models');
            this.models = await response.json();
            this.renderModels();
        } catch (error) {
            console.error('[IMMUNOS] Error loading models:', error);
        }
    },

    async loadTokensComparison() {
        try {
            const response = await fetch('/api/tokens/comparison?hours=24');
            if (!response.ok) throw new Error('Failed to load token comparison');
            const data = await response.json();
            this.updateTokensTab(data);
        } catch (error) {
            console.error('[IMMUNOS] Error loading tokens comparison:', error);
        }
    },

    async loadSessionStatus() {
        try {
            const response = await fetch('/api/tokens/session/current');
            if (!response.ok) throw new Error('Failed to load session');
            const data = await response.json();
            this.sessionId = data.session_id;
            this.updateSessionTokens(data.total_tokens);
        } catch (error) {
            console.error('[IMMUNOS] Error loading session:', error);
        }
    },

    async loadRoutingStatus() {
        try {
            const response = await fetch('/api/routing/status');
            if (!response.ok) throw new Error('Failed to load routing status');
            const data = await response.json();
            // Update system status
            this.updateSystemStatus(data);
        } catch (error) {
            console.error('[IMMUNOS] Error loading routing status:', error);
        }
    },

    // ========================================================================
    // EVENT HANDLERS
    // ========================================================================

    handleTrainingStarted(data) {
        console.log('[IMMUNOS] Training started:', data);
        this.trainingJobs[data.domain] = {
            domain: data.domain,
            dataset: data.dataset,
            num_samples: data.num_samples,
            current: 0,
            total: data.num_samples,
            detectors_created: 0,
            current_accuracy: 0,
            started_at: data.timestamp
        };
        this.renderTrainingJob(data.domain);
        this.addEventToFeed('training_started', data);
    },

    handleTrainingProgress(data) {
        if (this.trainingJobs[data.domain]) {
            Object.assign(this.trainingJobs[data.domain], data);
            this.updateTrainingJobProgress(data.domain);
        }
    },

    handleTrainingComplete(data) {
        console.log('[IMMUNOS] Training complete:', data);
        if (this.trainingJobs[data.domain]) {
            delete this.trainingJobs[data.domain];
            this.removeTrainingJob(data.domain);
        }
        this.addEventToFeed('training_complete', data);
    },

    handleModelSwitched(data) {
        console.log('[IMMUNOS] Model switched:', data);
        this.addEventToFeed('model_switched', data);
        if (data.session_tokens) {
            this.updateSessionTokens(data.session_tokens);
        }
    },

    handleModelLoaded(data) {
        console.log('[IMMUNOS] Model loaded:', data);
        this.addEventToFeed('model_loaded', data);
        this.loadModelsStatus();  // Refresh models
    },

    handleModelUnloaded(data) {
        console.log('[IMMUNOS] Model unloaded:', data);
        this.addEventToFeed('model_unloaded', data);
        this.loadModelsStatus();  // Refresh models
    },

    handleDetectionEvent(data) {
        console.log('[IMMUNOS] Detection event:', data);
        this.addEventToFeed('detection_event', data);
    },

    handleResourceUpdate(data) {
        console.log('[IMMUNOS] Resource update:', data);
        // Update resource metrics if needed
    },

    handleOllamaServerStarted(data) {
        console.log('[IMMUNOS] Ollama server started:', data);
        this.addEventToFeed('ollama_server_started', data);
        this.updateOllamaServerStatus(true);
    },

    handleOllamaStopped(data) {
        console.log('[IMMUNOS] Ollama server stopped:', data);
        this.addEventToFeed('ollama_server_stopped', data);
        this.updateOllamaServerStatus(false);
    },

    // ========================================================================
    // MODEL CONTROLS
    // ========================================================================

    async loadModel(modelName) {
        try {
            const response = await fetch(`/api/models/${modelName}/load`, {
                method: 'POST'
            });
            if (response.ok) {
                console.log(`[IMMUNOS] Loaded ${modelName}`);
            } else {
                console.error(`[IMMUNOS] Failed to load ${modelName}`);
            }
        } catch (error) {
            console.error('[IMMUNOS] Error loading model:', error);
        }
    },

    async unloadModel(modelName) {
        try {
            const response = await fetch(`/api/models/${modelName}/unload`, {
                method: 'POST'
            });
            if (response.ok) {
                console.log(`[IMMUNOS] Unloaded ${modelName}`);
            } else {
                console.error(`[IMMUNOS] Failed to unload ${modelName}`);
            }
        } catch (error) {
            console.error('[IMMUNOS] Error unloading model:', error);
        }
    },

    async startOllama() {
        try {
            const response = await fetch('/api/ollama/start', { method: 'POST' });
            if (response.ok) {
                console.log('[IMMUNOS] Ollama server starting...');
            } else {
                console.error('[IMMUNOS] Failed to start Ollama');
            }
        } catch (error) {
            console.error('[IMMUNOS] Error starting Ollama:', error);
        }
    },

    async stopOllama() {
        try {
            const response = await fetch('/api/ollama/stop', { method: 'POST' });
            if (response.ok) {
                console.log('[IMMUNOS] Ollama server stopping...');
            } else {
                console.error('[IMMUNOS] Failed to stop Ollama');
            }
        } catch (error) {
            console.error('[IMMUNOS] Error stopping Ollama:', error);
        }
    },

    // ========================================================================
    // RENDERING METHODS
    // ========================================================================

    renderModels() {
        const container = document.getElementById('models-container');
        if (!container) return;

        if (Object.keys(this.models).length === 0) {
            container.innerHTML = `
                <div class="col-span-3 text-center text-gray-500 text-sm py-8">
                    No models available. Start Ollama server to see models.
                </div>
            `;
            return;
        }

        container.innerHTML = '';
        Object.entries(this.models).forEach(([name, status]) => {
            const card = this.createModelCard(name, status);
            container.appendChild(card);
        });

        // Update active models count
        const activeCount = Object.values(this.models).filter(m => m.loaded).length;
        const activeModelsElem = document.getElementById('active-models');
        if (activeModelsElem) {
            activeModelsElem.textContent = activeCount;
        }
    },

    createModelCard(name, status) {
        const card = document.createElement('div');
        card.className = 'bg-gray-800 rounded-lg p-6 border border-gray-700';

        const statusClass = status.loaded ? 'text-green-400' : 'text-gray-500';
        const statusText = status.loaded ? '● Loaded' : '○ Unloaded';
        const lastUsed = status.last_used ? new Date(status.last_used).toLocaleString() : 'Never';
        const safeName = escapeHtml(name);
        const ramMb = escapeHtml(status.ram_mb ?? 0);
        const totalRequests = escapeHtml(status.total_requests ?? 0);
        const avgLatency = escapeHtml((status.avg_latency_ms || 0).toFixed(0));
        const safeLastUsed = escapeHtml(lastUsed);

        card.innerHTML = `
            <div class="flex justify-between items-start mb-4">
                <div>
                    <h3 class="font-bold text-lg">${safeName}</h3>
                </div>
                <span class="${statusClass} text-sm">${statusText}</span>
            </div>
            <div class="grid grid-cols-2 gap-3 text-sm mb-4">
                <div>
                    <span class="text-gray-400">RAM:</span>
                    <span class="ml-2 font-mono">${ramMb} MB</span>
                </div>
                <div>
                    <span class="text-gray-400">Requests:</span>
                    <span class="ml-2 font-mono">${totalRequests}</span>
                </div>
                <div class="col-span-2">
                    <span class="text-gray-400">Latency:</span>
                    <span class="ml-2 font-mono">${avgLatency} ms</span>
                </div>
                <div class="col-span-2">
                    <span class="text-gray-400">Last Used:</span>
                    <span class="ml-2 font-mono text-xs">${safeLastUsed}</span>
                </div>
            </div>
        `;

        const actions = document.createElement('div');
        actions.className = 'flex gap-2';
        const button = document.createElement('button');
        button.type = 'button';
        button.className = `flex-1 px-4 py-2 ${status.loaded ? 'bg-red-600 hover:bg-red-700' : 'bg-green-600 hover:bg-green-700'} rounded-lg text-sm font-medium transition-colors`;
        button.textContent = status.loaded ? 'Unload' : 'Load';
        button.addEventListener('click', () => {
            if (status.loaded) {
                this.unloadModel(name);
            } else {
                this.loadModel(name);
            }
        });
        actions.appendChild(button);
        card.appendChild(actions);

        return card;
    },

    addEventToFeed(eventType, data) {
        const feed = document.getElementById('events-feed');
        if (!feed) return;

        // Remove "no events" message if present
        const noEventsMsg = feed.querySelector('.text-center.text-gray-500');
        if (noEventsMsg) {
            feed.innerHTML = '';
        }

        const event = this.createEventElement(eventType, data);
        feed.insertBefore(event, feed.firstChild);

        // Keep only last 20 events
        while (feed.children.length > 20) {
            feed.removeChild(feed.lastChild);
        }

        // Track for activity chart
        this.trackEvent(eventType);
    },

    createEventElement(eventType, data) {
        const div = document.createElement('div');
        div.className = 'flex items-start space-x-3 p-3 bg-gray-700/50 rounded-lg';

        let icon = '●';
        let color = 'text-blue-400';
        let message = '';
        const safeDomain = escapeHtml(data?.domain || 'unknown');
        const safeFrom = escapeHtml(data?.from_model || 'unknown');
        const safeTo = escapeHtml(data?.to_model || 'unknown');
        const safeReason = escapeHtml(data?.reason || '');
        const safeModel = escapeHtml(data?.model || 'unknown');

        switch (eventType) {
            case 'training_started':
                icon = '🎯';
                color = 'text-blue-400';
                message = `Training started: ${safeDomain} (${data.num_samples ?? 0} samples)`;
                break;
            case 'training_complete':
                icon = '✅';
                color = 'text-green-400';
                message = `Training complete: ${safeDomain} (${data.detectors_created ?? 0} detectors, ${((data.final_accuracy || 0) * 100).toFixed(1)}% accuracy)`;
                break;
            case 'model_switched':
                icon = '🔄';
                color = 'text-yellow-400';
                message = `Model switched: ${safeFrom} → ${safeTo} (${safeReason || 'unknown'})`;
                break;
            case 'detection_event':
                icon = '🔍';
                color = data.result === 'non_self' ? 'text-red-400' : 'text-green-400';
                message = `Detection: ${safeDomain} - ${escapeHtml(data.result || 'unknown')} (${((data.confidence || 0) * 100).toFixed(1)}%)`;
                break;
            case 'model_loaded':
                icon = '⬆️';
                color = 'text-green-400';
                message = `Model loaded: ${safeModel}`;
                break;
            case 'model_unloaded':
                icon = '⬇️';
                color = 'text-gray-400';
                message = `Model unloaded: ${safeModel}`;
                break;
            case 'ollama_server_started':
                icon = '🚀';
                color = 'text-green-400';
                message = 'Ollama server started';
                break;
            case 'ollama_server_stopped':
                icon = '🛑';
                color = 'text-red-400';
                message = 'Ollama server stopped';
                break;
            default:
                message = escapeHtml(JSON.stringify(data));
        }

        const timestamp = new Date(data.timestamp || Date.now()).toLocaleTimeString();

        div.innerHTML = `
            <div class="${color} text-lg">${icon}</div>
            <div class="flex-1">
                <div class="text-sm">${message}</div>
                <div class="text-xs text-gray-500 mt-1">${escapeHtml(timestamp)}</div>
            </div>
        `;

        return div;
    },

    trackEvent(eventType) {
        // Increment current minute's event count
        this.eventHistory[this.eventHistory.length - 1]++;
        this.updateActivityChart();
    },

    updateSessionTokens(tokens) {
        const tokenElem = document.getElementById('session-tokens');
        if (tokenElem) {
            tokenElem.textContent = tokens.toLocaleString();
        }

        // Update progress bar
        const percent = (tokens / 200000) * 100;
        const progressBar = document.getElementById('token-progress-bar');
        if (progressBar) {
            progressBar.style.width = `${percent}%`;

            // Change color based on threshold
            if (tokens >= 180000) {
                progressBar.className = 'h-full bg-red-500 transition-all duration-300';
            } else if (tokens >= 150000) {
                progressBar.className = 'h-full bg-yellow-500 transition-all duration-300';
            } else {
                progressBar.className = 'h-full bg-blue-500 transition-all duration-300';
            }
        }

        // Update text
        const progressText = document.getElementById('token-progress-text');
        if (progressText) {
            progressText.textContent = `${tokens.toLocaleString()} / 200,000`;
        }

        // Show alerts
        const alertsDiv = document.getElementById('token-alerts');
        if (alertsDiv) {
            if (tokens >= 180000) {
                alertsDiv.innerHTML = '<div class="bg-red-900/50 border border-red-700 rounded-lg p-4 text-sm">🚨 <strong>Auto-Switch Active:</strong> Routing all tasks to Ollama to conserve tokens</div>';
            } else if (tokens >= 150000) {
                alertsDiv.innerHTML = '<div class="bg-yellow-900/50 border border-yellow-700 rounded-lg p-4 text-sm">⚠️ <strong>Warning:</strong> Approaching token limit. Routine tasks being routed to Ollama.</div>';
            } else {
                alertsDiv.innerHTML = '';
            }
        }
    },

    updateTokensTab(data) {
        // Claude stats
        const claudeTokens = document.getElementById('claude-tokens');
        const claudeRequests = document.getElementById('claude-requests');
        const claudeCost = document.getElementById('claude-cost');
        const claudeLatency = document.getElementById('claude-latency');

        if (claudeTokens) claudeTokens.textContent = (data.claude?.tokens || 0).toLocaleString();
        if (claudeRequests) claudeRequests.textContent = (data.claude?.requests || 0).toLocaleString();
        if (claudeCost) claudeCost.textContent = `$${(data.claude?.cost_usd || 0).toFixed(2)}`;
        if (claudeLatency) claudeLatency.textContent = `${(data.claude?.avg_latency_ms || 0).toFixed(0)} ms`;

        // Ollama stats
        const ollamaTokens = document.getElementById('ollama-tokens');
        const ollamaRequests = document.getElementById('ollama-requests');
        const ollamaCost = document.getElementById('ollama-cost');
        const ollamaLatency = document.getElementById('ollama-latency');

        if (ollamaTokens) ollamaTokens.textContent = (data.ollama?.tokens || 0).toLocaleString();
        if (ollamaRequests) ollamaRequests.textContent = (data.ollama?.requests || 0).toLocaleString();
        if (ollamaCost) ollamaCost.textContent = `$${(data.ollama?.cost_usd || 0).toFixed(2)}`;
        if (ollamaLatency) ollamaLatency.textContent = `${(data.ollama?.avg_latency_ms || 0).toFixed(0)} ms`;

        // Savings
        const totalSavings = document.getElementById('total-savings');
        const savingsPercent = document.getElementById('savings-percent');
        const savings24h = document.getElementById('savings-24h');

        if (totalSavings) totalSavings.textContent = `$${(data.savings_usd || 0).toFixed(2)}`;
        if (savingsPercent) savingsPercent.textContent = `${(data.savings_percent || 0).toFixed(1)}% reduction in costs`;
        if (savings24h) savings24h.textContent = `$${(data.savings_usd || 0).toFixed(2)}`;

        // Total requests
        const totalRequests = document.getElementById('total-requests');
        if (totalRequests) {
            const total = (data.claude?.requests || 0) + (data.ollama?.requests || 0);
            totalRequests.textContent = total.toLocaleString();
        }
    },

    updateOllamaServerStatus(running) {
        const statusElem = document.getElementById('ollama-server-status');
        const ollamaStatus = document.getElementById('ollama-status');

        if (statusElem) {
            statusElem.textContent = running ? 'Running' : 'Stopped';
            statusElem.className = running ? 'text-green-400' : 'text-red-400';
        }

        if (ollamaStatus) {
            ollamaStatus.textContent = running ? '● Running' : '● Stopped';
            ollamaStatus.className = running ? 'text-green-400' : 'text-red-400';
        }
    },

    updateSystemStatus(data) {
        // Update system status card in overview tab
        // This can be expanded with more routing status info
    },

    // ========================================================================
    // TRAINING JOB RENDERING
    // ========================================================================

    renderTrainingJob(domain) {
        const container = document.getElementById('training-jobs-container');
        if (!container) return;

        // Remove "no active jobs" message
        const noJobsMsg = container.querySelector('.text-center.text-gray-500');
        if (noJobsMsg) {
            container.innerHTML = '';
        }

        const job = this.trainingJobs[domain];
        const jobCard = document.createElement('div');
        jobCard.id = `training-job-${domain}`;
        jobCard.className = 'bg-gray-800 rounded-lg p-6 border border-gray-700';

        const percent = job.total > 0 ? (job.current / job.total) * 100 : 0;
        const safeDomain = escapeHtml(domain);
        const safeDataset = escapeHtml(job.dataset || '');

        jobCard.innerHTML = `
            <div class="flex justify-between items-start mb-4">
                <div>
                    <h4 class="font-bold text-lg">${safeDomain}</h4>
                    <p class="text-sm text-gray-400">${safeDataset}</p>
                </div>
                <div class="text-right">
                    <div class="text-sm text-gray-400">Progress</div>
                    <div class="text-xl font-bold text-blue-400">${percent.toFixed(1)}%</div>
                </div>
            </div>
            <div class="mb-4">
                <div class="h-4 bg-gray-700 rounded-full overflow-hidden">
                    <div class="h-full bg-blue-500 transition-all duration-300" style="width: ${percent}%"></div>
                </div>
            </div>
            <div class="grid grid-cols-3 gap-4 text-sm">
                <div>
                    <span class="text-gray-400">Samples:</span>
                    <span class="ml-2 font-mono">${job.current} / ${job.total}</span>
                </div>
                <div>
                    <span class="text-gray-400">Detectors:</span>
                    <span class="ml-2 font-mono">${job.detectors_created}</span>
                </div>
                <div>
                    <span class="text-gray-400">Accuracy:</span>
                    <span class="ml-2 font-mono">${(job.current_accuracy * 100).toFixed(1)}%</span>
                </div>
            </div>
        `;

        container.appendChild(jobCard);
    },

    updateTrainingJobProgress(domain) {
        const jobCard = document.getElementById(`training-job-${domain}`);
        if (!jobCard) return;

        const job = this.trainingJobs[domain];
        const percent = job.total > 0 ? (job.current / job.total) * 100 : 0;

        // Update progress bar
        const progressBar = jobCard.querySelector('.bg-blue-500');
        if (progressBar) {
            progressBar.style.width = `${percent}%`;
        }

        // Update percentage text
        const percentText = jobCard.querySelector('.text-blue-400');
        if (percentText) {
            percentText.textContent = `${percent.toFixed(1)}%`;
        }

        // Update stats
        const stats = jobCard.querySelectorAll('.font-mono');
        if (stats[0]) stats[0].textContent = `${job.current} / ${job.total}`;
        if (stats[1]) stats[1].textContent = job.detectors_created;
        if (stats[2]) stats[2].textContent = `${(job.current_accuracy * 100).toFixed(1)}%`;
    },

    removeTrainingJob(domain) {
        const jobCard = document.getElementById(`training-job-${domain}`);
        if (jobCard) {
            jobCard.remove();
        }

        // Show "no active jobs" if container is empty
        const container = document.getElementById('training-jobs-container');
        if (container && container.children.length === 0) {
            container.innerHTML = `
                <div class="text-center text-gray-500 text-sm py-8 bg-gray-800 rounded-lg border border-gray-700">
                    No active training jobs
                </div>
            `;
        }
    },

    // ========================================================================
    // CHARTS
    // ========================================================================

    setupCharts() {
        // Activity chart (events per minute)
        const activityCtx = document.getElementById('activity-chart');
        if (activityCtx) {
            this.activityChart = new Chart(activityCtx.getContext('2d'), {
                type: 'line',
                data: {
                    labels: this.timeLabels,
                    datasets: [{
                        label: 'Events/min',
                        data: this.eventHistory,
                        borderColor: '#3b82f6',
                        backgroundColor: 'rgba(59, 130, 246, 0.1)',
                        tension: 0.4,
                        fill: true
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: {
                            display: false
                        }
                    },
                    scales: {
                        y: {
                            beginAtZero: true,
                            grid: {
                                color: '#374151'
                            },
                            ticks: {
                                color: '#9ca3af'
                            }
                        },
                        x: {
                            grid: {
                                color: '#374151'
                            },
                            ticks: {
                                color: '#9ca3af'
                            }
                        }
                    }
                }
            });
        }

        // Tokens chart
        const tokensCtx = document.getElementById('tokens-chart');
        if (tokensCtx) {
            this.tokensChart = new Chart(tokensCtx.getContext('2d'), {
                type: 'line',
                data: {
                    labels: [],
                    datasets: [
                        {
                            label: 'Claude',
                            data: [],
                            borderColor: '#3b82f6',
                            backgroundColor: 'rgba(59, 130, 246, 0.1)',
                            tension: 0.4,
                            fill: true
                        },
                        {
                            label: 'Ollama',
                            data: [],
                            borderColor: '#10b981',
                            backgroundColor: 'rgba(16, 185, 129, 0.1)',
                            tension: 0.4,
                            fill: true
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: {
                            labels: {
                                color: '#9ca3af'
                            }
                        }
                    },
                    scales: {
                        y: {
                            beginAtZero: true,
                            grid: {
                                color: '#374151'
                            },
                            ticks: {
                                color: '#9ca3af'
                            }
                        },
                        x: {
                            grid: {
                                color: '#374151'
                            },
                            ticks: {
                                color: '#9ca3af'
                            }
                        }
                    }
                }
            });
        }
    },

    updateActivityChart() {
        if (this.activityChart) {
            this.activityChart.data.datasets[0].data = this.eventHistory;
            this.activityChart.update('none'); // Update without animation for performance
        }
    },

    // ========================================================================
    // TAB MANAGEMENT
    // ========================================================================

    setupTabs() {
        const tabButtons = document.querySelectorAll('.tab-btn');
        tabButtons.forEach(btn => {
            btn.addEventListener('click', () => {
                this.switchTab(btn.dataset.tab);
            });
        });
    },

    switchTab(tabName) {
        this.currentTab = tabName;

        // Update button states
        document.querySelectorAll('.tab-btn').forEach(btn => {
            if (btn.dataset.tab === tabName) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
        });

        // Update panel visibility
        document.querySelectorAll('.tab-panel').forEach(panel => {
            panel.classList.add('hidden');
        });

        const activePanel = document.getElementById(`${tabName}-tab`);
        if (activePanel) {
            activePanel.classList.remove('hidden');
        }
    },

    // ========================================================================
    // DOMAIN MANAGEMENT
    // ========================================================================

    setupDomainButtons() {
        const domainButtons = document.querySelectorAll('.domain-btn');
        domainButtons.forEach(btn => {
            btn.addEventListener('click', () => {
                this.switchDomain(btn.dataset.domain);
            });
        });
    },

    switchDomain(domainName) {
        this.currentDomain = domainName;

        // Update button states
        document.querySelectorAll('.domain-btn').forEach(btn => {
            if (btn.dataset.domain === domainName) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
        });

        // Update domain content
        this.renderDomainContent(domainName);
    },

    renderDomainContent(domainName) {
        const info = this.domainInfo[domainName];
        if (!info) return;

        const titleElem = document.getElementById('domain-title');
        const descElem = document.getElementById('domain-description');

        if (titleElem) titleElem.textContent = info.name;
        if (descElem) descElem.textContent = info.description;

        // TODO: Load domain-specific stats from API
        // For now, show placeholder data
    },

    // ========================================================================
    // EVENT LISTENERS
    // ========================================================================

    setupEventListeners() {
        // Ollama server controls
        const btnOllamaStart = document.getElementById('btn-ollama-start');
        const btnOllamaStop = document.getElementById('btn-ollama-stop');

        if (btnOllamaStart) {
            btnOllamaStart.addEventListener('click', () => {
                this.startOllama();
            });
        }

        if (btnOllamaStop) {
            btnOllamaStop.addEventListener('click', () => {
                this.stopOllama();
            });
        }
    },

    // ========================================================================
    // PERIODIC UPDATES
    // ========================================================================

    startPeriodicUpdates() {
        // Update every 5 seconds
        setInterval(() => {
            this.loadTokensComparison();
            this.loadModelsStatus();
        }, 5000);

        // Shift event history every minute
        setInterval(() => {
            this.eventHistory.shift();
            this.eventHistory.push(0);

            // Update time labels
            const now = new Date();
            this.timeLabels.shift();
            this.timeLabels.push(now.toLocaleTimeString('en-US', { hour: '2-digit', minute: '2-digit' }));

            if (this.activityChart) {
                this.activityChart.data.labels = this.timeLabels;
                this.activityChart.update();
            }
        }, 60000); // Every minute
    }
};

// Initialize when DOM ready
document.addEventListener('DOMContentLoaded', () => {
    console.log('[IMMUNOS] DOM loaded, initializing monitor...');
    Monitor.init();
});

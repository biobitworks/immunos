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
    currentTab: 'chat',
    currentDomain: 'emotion',
    models: {},
    trainingJobs: {},
    domainStats: {},
    issuesSummary: {},
    issuesList: [],
    issuesProject: 'immunos',
    orchestratorConfig: null,
    thymusQueue: [],
    kbPages: [],
    foundryTemplates: [],
    foundryAgents: [],
    spleenSummary: {},
    sessionId: null,

    // Charts
    activityChart: null,
    tokensChart: null,
    thymusChart: null,

    agentGraph: null,
    agentGraphCanvas: null,
    agentGraphCtx: null,
    agentGraphConfig: { showLabels: true },
    agentGraphAnimating: false,
    lastDetection: null,

    // Event history for activity chart (last 20 minutes)
    eventHistory: Array(20).fill(0),
    timeLabels: [],

    // Domain information
    domainInfo: {
        emotion: {
            name: 'Emotion Recognition',
            description: 'Detects anomalous emotional context using local feature extractors and AIS pattern matching. Targets CK+ / FER2013-style signals.'
        },
        hallucination: {
            name: 'Hallucination Detection',
            description: 'Tags hallucinations with B Cell + NK using SimpleText embeddings and TruthfulQA/HaluEval training data.'
        },
        network: {
            name: 'Network Intrusion',
            description: 'Detects intrusions using Enhanced NK on NSL-KDD (41-dim features) and negative selection.'
        },
        code: {
            name: 'Code Vulnerability',
            description: 'Scans code for vulnerabilities using AST features, negative selection, and curated safe/vuln patterns.'
        },
        research: {
            name: 'Research Verification',
            description: 'Verifies scientific claims using SciFact with B Cell + NK and SimpleText embeddings.'
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
        this.setupThymusChart();
        this.setupThymusIntake();
        this.setupFoundry();
        this.setupAgentGraph();
        this.setupChat();
        this.setupQuickChat();
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
                this.loadOllamaStatus(),
                this.loadRoutingStatus(),
                this.loadDomainStatus(),
                this.loadIssues(),
                this.loadSpleenSummary(),
                this.loadOrchestratorStatus(),
                this.loadThymusQueue(),
                this.loadThymusStatus(),
                this.loadKbIndex(),
                this.loadFoundryTemplates(),
                this.loadFoundryHistory()
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

    async loadOllamaStatus() {
        try {
            const response = await fetch('/api/ollama/status');
            if (!response.ok) throw new Error('Failed to load Ollama status');
            const payload = await response.json();
            const status = payload.data?.status;
            this.updateOllamaServerStatus(status === 'online');
        } catch (error) {
            console.error('[IMMUNOS] Error loading Ollama status:', error);
            this.updateOllamaServerStatus(false);
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

    async loadDomainStatus() {
        try {
            const response = await fetch('/api/domains/status');
            if (!response.ok) throw new Error('Failed to load domain status');
            const data = await response.json();
            this.domainStats = data;
            this.updateDomainStats(this.currentDomain);
            this.updateSystemStatus(data);
            this.updateThymusOverview();
        } catch (error) {
            console.error('[IMMUNOS] Error loading domain status:', error);
        }
    },

    async loadIssues() {
        try {
            const projectParam = this.issuesProject ? `?project=${encodeURIComponent(this.issuesProject)}` : '';
            const listParam = this.issuesProject ? `&project=${encodeURIComponent(this.issuesProject)}` : '';
            const [summaryResp, listResp] = await Promise.all([
                fetch(`/api/issues/summary${projectParam}`),
                fetch(`/api/issues/list?limit=10${listParam}`)
            ]);

            if (!summaryResp.ok || !listResp.ok) {
                throw new Error('Failed to load issues');
            }

            const summaryData = await summaryResp.json();
            const listData = await listResp.json();
            this.issuesSummary = summaryData.summary || {};
            this.issuesList = listData.items || [];
            this.renderIssues();
        } catch (error) {
            console.error('[IMMUNOS] Error loading issues:', error);
        }
    },

    async loadSpleenSummary() {
        try {
            const projectParam = this.issuesProject ? `?project=${encodeURIComponent(this.issuesProject)}` : '';
            const response = await fetch(`/api/spleen/summary${projectParam}`);
            if (!response.ok) throw new Error('Failed to load spleen summary');
            const data = await response.json();
            this.spleenSummary = data || {};
            this.updateSpleenSummary();
        } catch (error) {
            console.error('[IMMUNOS] Error loading spleen summary:', error);
        }
    },

    updateSpleenSummary() {
        const anomalies = this.spleenSummary?.anomalies || {};
        const issues = this.spleenSummary?.issues || {};
        const lastScan = this.spleenSummary?.last_scan;

        const unresolvedElem = document.getElementById('spleen-unresolved');
        const highElem = document.getElementById('spleen-high');
        const resolvedElem = document.getElementById('spleen-resolved');
        const openIssuesElem = document.getElementById('spleen-open-issues');
        const lastScanElem = document.getElementById('spleen-last-scan');

        if (unresolvedElem) unresolvedElem.textContent = anomalies.unresolved_total ?? 0;
        if (highElem) highElem.textContent = anomalies.high ?? 0;
        if (resolvedElem) resolvedElem.textContent = anomalies.resolved_total ?? 0;
        if (openIssuesElem) openIssuesElem.textContent = issues.active ?? 0;
        if (lastScanElem) {
            lastScanElem.textContent = lastScan ? new Date(lastScan).toLocaleString() : '—';
        }
    },

    async loadOrchestratorStatus() {
        try {
            const response = await fetch('/api/orchestrator/status');
            if (!response.ok) throw new Error('Failed to load orchestrator status');
            const data = await response.json();
            this.orchestratorConfig = data;
            this.updateOrchestratorUI();
        } catch (error) {
            console.error('[IMMUNOS] Error loading orchestrator status:', error);
        }
    },

    async loadThymusQueue() {
        try {
            const response = await fetch('/api/thymus/intake?limit=10');
            if (!response.ok) throw new Error('Failed to load thymus queue');
            const data = await response.json();
            this.thymusQueue = data.items || [];
            this.renderThymusQueue();
        } catch (error) {
            console.error('[IMMUNOS] Error loading thymus queue:', error);
        }
    },

    async loadKbIndex() {
        try {
            const response = await fetch('/api/kb/index');
            if (!response.ok) throw new Error('Failed to load KB index');
            const data = await response.json();
            this.kbPages = data.pages || [];
            this.renderKbList();
            if (this.kbPages.length) {
                this.loadKbPage(this.kbPages[0].name);
            }
        } catch (error) {
            console.error('[IMMUNOS] Error loading KB index:', error);
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
        this.updateTrainingJobCount();
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
        this.updateTrainingJobCount();
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
        this.updateAgentGraph(data);
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

    async simulateDetection() {
        try {
            const response = await fetch('/api/events/simulate', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ domain: this.currentDomain })
            });
            if (!response.ok) {
                console.error('[IMMUNOS] Failed to simulate detection');
            }
        } catch (error) {
            console.error('[IMMUNOS] Error simulating detection:', error);
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
            case 'thymus_intake':
                icon = '🧬';
                color = 'text-blue-400';
                message = 'Thymus intake queued';
                break;
            case 'orchestrator_update':
                icon = '🧭';
                color = 'text-blue-400';
                message = 'Orchestrator configuration updated';
                break;
            default:
                message = escapeHtml(JSON.stringify(data));
        }

        const timestamp = new Date(data.timestamp || Date.now()).toLocaleTimeString();

        div.innerHTML = `
            <div class="${color} text-lg">${icon}</div>
            <div class="flex-1">
                <div class="text-sm">${message}</div>
                <div class="text-xs text-gray-500 mt-1">${timestamp}</div>
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
        const chatStatus = document.getElementById('chat-ollama-status');
        const offlineBanner = document.getElementById('chat-ollama-offline');

        if (statusElem) {
            statusElem.textContent = running ? 'Running' : 'Stopped';
            statusElem.className = running ? 'text-green-400' : 'text-red-400';
        }

        if (ollamaStatus) {
            ollamaStatus.textContent = running ? '● Running' : '● Stopped';
            ollamaStatus.className = running ? 'text-green-400' : 'text-red-400';
        }

        if (chatStatus) {
            chatStatus.textContent = running ? 'Running' : 'Stopped';
            chatStatus.className = running ? 'text-xs text-green-400' : 'text-xs text-red-400';
        }

        if (offlineBanner) {
            if (running) {
                offlineBanner.classList.add('hidden');
            } else {
                offlineBanner.classList.remove('hidden');
            }
        }
    },

    updateSystemStatus(data) {
        // Update system status card in overview tab
        const trainedDomainsElem = document.getElementById('trained-domains');
        if (trainedDomainsElem && data?.domains) {
            const total = Object.keys(this.domainInfo || {}).length;
            const trained = Object.values(data.domains || {}).filter(domain => domain.trained).length;
            trainedDomainsElem.textContent = `${trained} / ${total}`;
        }
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

    updateTrainingJobCount() {
        const countElem = document.getElementById('training-jobs');
        if (!countElem) return;
        countElem.textContent = Object.keys(this.trainingJobs).length.toString();
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

    setupThymusChart() {
        const thymusCtx = document.getElementById('thymus-training-chart');
        if (!thymusCtx) return;

        this.thymusChart = new Chart(thymusCtx.getContext('2d'), {
            type: 'bar',
            data: {
                labels: [],
                datasets: [
                    {
                        label: 'B Cell Patterns',
                        data: [],
                        backgroundColor: 'rgba(59, 130, 246, 0.6)'
                    },
                    {
                        label: 'NK Detectors',
                        data: [],
                        backgroundColor: 'rgba(16, 185, 129, 0.6)'
                    },
                    {
                        label: 'NK Self Patterns',
                        data: [],
                        backgroundColor: 'rgba(168, 85, 247, 0.6)'
                    }
                ]
            },
            options: {
                indexAxis: 'y',
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
    },

    updateActivityChart() {
        if (this.activityChart) {
            this.activityChart.data.datasets[0].data = this.eventHistory;
            this.activityChart.update('none'); // Update without animation for performance
        }
    },

    updateThymusOverview() {
        const stats = this.domainStats?.domains || {};
        const domainKeys = Object.keys(this.domainInfo || {});

        const labels = [];
        const patterns = [];
        const detectors = [];
        const selfPatterns = [];

        let trainedCount = 0;
        let totalPatterns = 0;
        let totalDetectors = 0;

        domainKeys.forEach((key) => {
            const info = this.domainInfo[key];
            const domainStats = stats[key] || {};
            const label = info?.name ? info.name.split(' ')[0] : key;

            labels.push(label);
            const bcell = domainStats.bcell_patterns ?? 0;
            const nk = domainStats.nk_detectors ?? 0;
            const self = domainStats.nk_self_patterns ?? 0;

            patterns.push(bcell);
            detectors.push(nk);
            selfPatterns.push(self);

            totalPatterns += bcell;
            totalDetectors += nk;
            if (domainStats.trained) {
                trainedCount += 1;
            }
        });

        const trainedElem = document.getElementById('thymus-trained-count');
        if (trainedElem) trainedElem.textContent = trainedCount.toString();

        const patternElem = document.getElementById('thymus-total-patterns');
        if (patternElem) patternElem.textContent = totalPatterns.toString();

        const detectorElem = document.getElementById('thymus-total-detectors');
        if (detectorElem) detectorElem.textContent = totalDetectors.toString();

        if (this.thymusChart) {
            this.thymusChart.data.labels = labels;
            this.thymusChart.data.datasets[0].data = patterns;
            this.thymusChart.data.datasets[1].data = detectors;
            this.thymusChart.data.datasets[2].data = selfPatterns;
            this.thymusChart.update();
        }

        this.renderThymusCards(domainKeys, stats);
    },

    renderThymusCards(domainKeys, stats) {
        const container = document.getElementById('thymus-domain-cards');
        if (!container) return;

        const cards = domainKeys.map((key) => {
            const info = this.domainInfo[key];
            const domainStats = stats[key] || {};
            const name = escapeHtml(info?.name || key);
            const trained = domainStats.trained ? 'Trained' : 'Untrained';
            const lastTrained = escapeHtml(domainStats.last_trained || 'Never');

            return `
                <div class="bg-gray-800 rounded-lg p-4 border border-gray-700">
                    <div class="flex justify-between items-center mb-2">
                        <h4 class="font-semibold text-sm">${name}</h4>
                        <span class="text-xs ${domainStats.trained ? 'text-green-400' : 'text-gray-500'}">${trained}</span>
                    </div>
                    <div class="grid grid-cols-3 gap-2 text-xs text-gray-400 mb-2">
                        <div>Patterns: <span class="font-mono text-blue-400">${domainStats.bcell_patterns ?? 0}</span></div>
                        <div>Detectors: <span class="font-mono text-green-400">${domainStats.nk_detectors ?? 0}</span></div>
                        <div>Self: <span class="font-mono text-purple-400">${domainStats.nk_self_patterns ?? 0}</span></div>
                    </div>
                    <div class="text-xs text-gray-500">Last trained: <span class="font-mono">${lastTrained}</span></div>
                </div>
            `;
        }).join('');

        container.innerHTML = cards || `
            <div class="col-span-3 text-center text-gray-500 text-sm py-6">
                No training stats available yet.
            </div>
        `;
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

        if (tabName === 'orchestrator') {
            this.startAgentGraph();
            setTimeout(() => this.resizeAgentGraph(), 50);
        } else {
            this.agentGraphAnimating = false;
        }

        if (tabName === 'thymus' && this.thymusChart) {
            setTimeout(() => this.thymusChart.resize(), 50);
        }
    },

    // ========================================================================
    // AGENT NETWORK MAP
    // ========================================================================

    setupAgentGraph() {
        const canvas = document.getElementById('orchestrator-canvas');
        if (!canvas) return;

        this.agentGraphCanvas = canvas;
        this.agentGraphCtx = canvas.getContext('2d');
        this.agentGraph = this.buildAgentGraph();
        this.resizeAgentGraph();

        window.addEventListener('resize', () => this.resizeAgentGraph());

        const toggleLabels = document.getElementById('toggle-agent-labels');
        if (toggleLabels) {
            toggleLabels.addEventListener('click', () => {
                this.agentGraphConfig.showLabels = !this.agentGraphConfig.showLabels;
            });
        }

        const resetBtn = document.getElementById('reset-agent-map');
        if (resetBtn) {
            resetBtn.addEventListener('click', () => this.clearAgentHighlights());
        }

        if (this.currentTab === 'orchestrator') {
            this.startAgentGraph();
        }
    },

    buildAgentGraph() {
        const nodes = [
            { id: 'domain_emotion', label: 'Emotion', x: 0.12, y: 0.18, type: 'domain' },
            { id: 'domain_hallucination', label: 'Hallucination', x: 0.12, y: 0.35, type: 'domain' },
            { id: 'domain_network', label: 'Network', x: 0.12, y: 0.52, type: 'domain' },
            { id: 'domain_code', label: 'Code', x: 0.12, y: 0.69, type: 'domain' },
            { id: 'domain_research', label: 'Research', x: 0.12, y: 0.86, type: 'domain' },
            { id: 'input', label: 'Input', x: 0.28, y: 0.5, type: 'input' },
            { id: 'dendritic', label: 'Dendritic', x: 0.38, y: 0.32, type: 'agent' },
            { id: 'orchestrator', label: 'Orchestrator', x: 0.52, y: 0.5, type: 'orchestrator' },
            { id: 'bcell', label: 'B Cell', x: 0.66, y: 0.22, type: 'agent' },
            { id: 'nk', label: 'NK Cell', x: 0.66, y: 0.78, type: 'agent' },
            { id: 'memory', label: 'Memory', x: 0.78, y: 0.35, type: 'agent' },
            { id: 'tcell', label: 'T Cell', x: 0.78, y: 0.65, type: 'agent' },
            { id: 'decision_self', label: 'Self', x: 0.92, y: 0.2, type: 'decision_self' },
            { id: 'decision_non_self', label: 'Non-Self', x: 0.92, y: 0.45, type: 'decision_non_self' },
            { id: 'decision_danger', label: 'Danger', x: 0.92, y: 0.7, type: 'decision_danger' },
            { id: 'decision_uncertain', label: 'Uncertain', x: 0.92, y: 0.9, type: 'decision_uncertain' }
        ];

        const edges = [
            { from: 'domain_emotion', to: 'orchestrator' },
            { from: 'domain_hallucination', to: 'orchestrator' },
            { from: 'domain_network', to: 'orchestrator' },
            { from: 'domain_code', to: 'orchestrator' },
            { from: 'domain_research', to: 'orchestrator' },
            { from: 'input', to: 'dendritic' },
            { from: 'dendritic', to: 'orchestrator' },
            { from: 'orchestrator', to: 'bcell' },
            { from: 'orchestrator', to: 'nk' },
            { from: 'orchestrator', to: 'memory' },
            { from: 'orchestrator', to: 'tcell' },
            { from: 'bcell', to: 'orchestrator' },
            { from: 'nk', to: 'orchestrator' },
            { from: 'memory', to: 'orchestrator' },
            { from: 'tcell', to: 'orchestrator' },
            { from: 'orchestrator', to: 'decision_self' },
            { from: 'orchestrator', to: 'decision_non_self' },
            { from: 'orchestrator', to: 'decision_danger' },
            { from: 'orchestrator', to: 'decision_uncertain' }
        ];

        return { nodes, edges };
    },

    resizeAgentGraph() {
        if (!this.agentGraphCanvas || !this.agentGraphCtx) return;
        const rect = this.agentGraphCanvas.getBoundingClientRect();
        const ratio = window.devicePixelRatio || 1;
        this.agentGraphCanvas.width = rect.width * ratio;
        this.agentGraphCanvas.height = rect.height * ratio;
        this.agentGraphCtx.setTransform(ratio, 0, 0, ratio, 0, 0);
    },

    startAgentGraph() {
        if (this.agentGraphAnimating) return;
        this.agentGraphAnimating = true;
        this.renderAgentGraph();
    },

    renderAgentGraph() {
        if (!this.agentGraphAnimating) return;
        if (!this.agentGraphCanvas || !this.agentGraphCtx || !this.agentGraph) return;
        this.drawAgentGraph();
        requestAnimationFrame(() => this.renderAgentGraph());
    },

    drawAgentGraph() {
        const ctx = this.agentGraphCtx;
        const canvas = this.agentGraphCanvas;
        if (!ctx || !canvas || !this.agentGraph) return;

        const width = canvas.clientWidth;
        const height = canvas.clientHeight;
        ctx.clearRect(0, 0, width, height);

        const now = Date.now();
        const colors = {
            orchestrator: '#60a5fa',
            agent: '#22d3ee',
            domain: '#a855f7',
            input: '#94a3b8',
            decision_self: '#34d399',
            decision_non_self: '#facc15',
            decision_danger: '#f87171',
            decision_uncertain: '#fb923c'
        };

        const getNode = (id) => this.agentGraph.nodes.find(node => node.id === id);

        // Draw edges
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.35)';
        ctx.lineWidth = 1;
        this.agentGraph.edges.forEach(edge => {
            const from = getNode(edge.from);
            const to = getNode(edge.to);
            if (!from || !to) return;
            ctx.beginPath();
            ctx.moveTo(from.x * width, from.y * height);
            ctx.lineTo(to.x * width, to.y * height);
            ctx.stroke();
        });

        // Draw nodes
        this.agentGraph.nodes.forEach(node => {
            const baseColor = colors[node.type] || '#e5e7eb';
            const age = node.lastActive ? (now - node.lastActive) : 99999;
            const active = age < 3500 ? (1 - age / 3500) : 0;
            const radius = node.type.startsWith('decision') ? 14 : (node.type === 'orchestrator' ? 18 : 12);
            const x = node.x * width;
            const y = node.y * height;

            if (active > 0) {
                ctx.beginPath();
                ctx.strokeStyle = baseColor;
                ctx.globalAlpha = 0.3 + (0.5 * active);
                ctx.lineWidth = 6;
                ctx.arc(x, y, radius + 8, 0, Math.PI * 2);
                ctx.stroke();
            }

            ctx.beginPath();
            ctx.fillStyle = baseColor;
            ctx.globalAlpha = 0.5 + (0.5 * (node.type === 'orchestrator' ? 1 : 0.7)) + active * 0.4;
            ctx.arc(x, y, radius, 0, Math.PI * 2);
            ctx.fill();
            ctx.globalAlpha = 1;

            if (this.agentGraphConfig.showLabels) {
                ctx.fillStyle = '#e5e7eb';
                ctx.font = '12px sans-serif';
                ctx.textAlign = 'center';
                ctx.fillText(node.label, x, y + radius + 14);
            }
        });
    },

    updateAgentGraph(data) {
        if (!data || !this.agentGraph) return;

        const result = (data.result || '').toLowerCase();
        const resultMap = {
            self: 'decision_self',
            non_self: 'decision_non_self',
            danger: 'decision_danger',
            uncertain: 'decision_uncertain'
        };

        this.markAgentActive('orchestrator', 1);
        this.markAgentActive('dendritic', 0.6);
        this.markAgentActive('bcell', 0.7);
        this.markAgentActive('nk', 0.9);
        this.markAgentActive('memory', 0.5);
        this.markAgentActive('tcell', 0.4);

        if (data.domain) {
            this.markAgentActive(`domain_${data.domain}`, 1);
        }

        const decisionNode = resultMap[result];
        if (decisionNode) {
            this.markAgentActive(decisionNode, 1);
        }

        this.updateLastDetectionPanel(data);
    },

    markAgentActive(nodeId, strength) {
        const node = this.agentGraph.nodes.find(n => n.id === nodeId);
        if (!node) return;
        node.lastActive = Date.now();
        node.activity = strength || 1;
    },

    clearAgentHighlights() {
        if (!this.agentGraph) return;
        this.agentGraph.nodes.forEach(node => {
            node.lastActive = null;
            node.activity = 0;
        });
    },

    updateLastDetectionPanel(data) {
        if (!data) return;
        const domainElem = document.getElementById('last-detection-domain');
        const resultElem = document.getElementById('last-detection-result');
        const confElem = document.getElementById('last-detection-confidence');
        const dangerElem = document.getElementById('last-detection-danger');
        const detElem = document.getElementById('last-detection-detectors');
        const timeElem = document.getElementById('last-detection-time');

        if (domainElem) domainElem.textContent = data.domain || 'unknown';

        if (resultElem) {
            const result = data.result || 'unknown';
            const colorClass = {
                self: 'text-green-400',
                non_self: 'text-yellow-300',
                danger: 'text-red-400',
                uncertain: 'text-orange-400'
            }[result] || 'text-gray-300';
            resultElem.textContent = result;
            resultElem.className = `font-mono ${colorClass}`;
        }

        if (confElem) {
            const confidence = typeof data.confidence === 'number' ? `${(data.confidence * 100).toFixed(1)}%` : '—';
            confElem.textContent = confidence;
        }

        if (dangerElem) {
            const danger = typeof data.danger_signal === 'number' ? data.danger_signal.toFixed(2) : '—';
            dangerElem.textContent = danger;
        }

        if (detElem) {
            const detectors = Array.isArray(data.matched_detectors) ? data.matched_detectors.length : 0;
            detElem.textContent = detectors ? `${detectors} matched` : '—';
        }

        if (timeElem) {
            const time = data.timestamp ? new Date(data.timestamp).toLocaleTimeString() : '—';
            timeElem.textContent = time;
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

        this.updateDomainStats(domainName);
    },

    updateDomainStats(domainName) {
        const stats = this.domainStats?.domains?.[domainName];
        if (!stats) return;

        const patternsElem = document.getElementById('domain-detections');
        const detectorsElem = document.getElementById('domain-confidence');
        const selfElem = document.getElementById('domain-detectors');
        const trainedElem = document.getElementById('domain-last-detection');

        if (patternsElem) patternsElem.textContent = stats.bcell_patterns ?? 0;
        if (detectorsElem) detectorsElem.textContent = stats.nk_detectors ?? 0;
        if (selfElem) selfElem.textContent = stats.nk_self_patterns ?? 0;
        if (trainedElem) trainedElem.textContent = stats.last_trained || 'Never';
    },

    setupThymusIntake() {
        const submitBtn = document.getElementById('thymus-submit');
        const pauseBtn = document.getElementById('thymus-pause');
        const resumeBtn = document.getElementById('thymus-resume');
        const runNextBtn = document.getElementById('thymus-run-next');
        if (submitBtn) {
            submitBtn.addEventListener('click', () => this.submitThymusIntake());
        }

        if (pauseBtn) {
            pauseBtn.addEventListener('click', () => this.controlThymusQueue('pause'));
        }
        if (resumeBtn) {
            resumeBtn.addEventListener('click', () => this.controlThymusQueue('resume'));
        }
        if (runNextBtn) {
            runNextBtn.addEventListener('click', () => this.controlThymusQueue('run_next'));
        }
    },

    setupFoundry() {
        const createBtn = document.getElementById('foundry-create');
        if (createBtn) {
            createBtn.addEventListener('click', () => this.submitFoundry());
        }
    },

    async submitThymusIntake() {
        const domainElem = document.getElementById('thymus-domain');
        const datasetElem = document.getElementById('thymus-dataset');
        const samplesElem = document.getElementById('thymus-samples');
        const notesElem = document.getElementById('thymus-notes');

        const domain = domainElem ? domainElem.value : '';
        const datasetPath = datasetElem ? datasetElem.value.trim() : '';
        const samples = samplesElem ? samplesElem.value.trim() : '';
        const notes = notesElem ? notesElem.value.trim() : '';

        if (!domain || (!datasetPath && !notes)) {
            this.addEventToFeed('thymus_intake', { timestamp: new Date().toISOString() });
            return;
        }

        try {
            const response = await fetch('/api/thymus/intake', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    domain,
                    dataset_path: datasetPath,
                    sample_count: samples ? Number(samples) : null,
                    notes
                })
            });

            if (!response.ok) {
                throw new Error('Failed to queue intake');
            }

            if (datasetElem) datasetElem.value = '';
            if (samplesElem) samplesElem.value = '';
            if (notesElem) notesElem.value = '';

            await this.loadThymusQueue();
            await this.loadThymusStatus();
            this.addEventToFeed('thymus_intake', { timestamp: new Date().toISOString() });
        } catch (error) {
            console.error('[IMMUNOS] Error submitting thymus intake:', error);
        }
    },

    async loadFoundryTemplates() {
        try {
            const response = await fetch('/api/agents/templates');
            if (!response.ok) throw new Error('Failed to load templates');
            const data = await response.json();
            this.foundryTemplates = data.templates || [];
            this.renderFoundryTemplates();
        } catch (error) {
            console.error('[IMMUNOS] Error loading foundry templates:', error);
        }
    },

    renderFoundryTemplates() {
        const listElem = document.getElementById('foundry-templates');
        const selectElem = document.getElementById('foundry-template');

        if (selectElem) {
            selectElem.innerHTML = '';
        }

        if (!this.foundryTemplates.length) {
            if (listElem) {
                listElem.innerHTML = '<div class="text-center text-gray-500 text-sm py-6">No templates available.</div>';
            }
            if (selectElem) {
                selectElem.innerHTML = '<option value="">No templates</option>';
            }
            return;
        }

        if (listElem) {
            listElem.innerHTML = this.foundryTemplates.map((template) => {
                const name = escapeHtml(template.name || 'Untitled');
                const description = escapeHtml(template.description || 'No description');
                const role = escapeHtml(template.role || 'unknown');
                return `
                    <div class="bg-gray-700/40 rounded-lg p-3 border border-gray-700">
                        <div class="font-semibold text-gray-200">${name}</div>
                        <div class="text-xs text-gray-400 mt-1">${description}</div>
                        <div class="text-xs text-gray-500 mt-2">Role: ${role}</div>
                    </div>
                `;
            }).join('');
        }

        if (selectElem) {
            selectElem.innerHTML = this.foundryTemplates.map((template) => `
                <option value="${escapeHtml(template.id)}">${escapeHtml(template.name || 'Template')}</option>
            `).join('');
        }
    },

    async loadFoundryHistory() {
        try {
            const response = await fetch('/api/agents/foundry?limit=8');
            if (!response.ok) throw new Error('Failed to load foundry history');
            const data = await response.json();
            this.foundryAgents = data.items || [];
            this.renderFoundryHistory();
        } catch (error) {
            console.error('[IMMUNOS] Error loading foundry history:', error);
        }
    },

    renderFoundryHistory() {
        const historyElem = document.getElementById('foundry-history');
        if (!historyElem) return;

        if (!this.foundryAgents.length) {
            historyElem.innerHTML = '<div class="text-center text-gray-500 text-sm py-4">No stubs yet.</div>';
            return;
        }

        historyElem.innerHTML = this.foundryAgents.map((agent) => {
            const created = agent.created_at ? new Date(agent.created_at).toLocaleString() : '—';
            const name = escapeHtml(agent.name || agent.id || 'Agent');
            const description = escapeHtml(agent.description || 'No description');
            const role = escapeHtml(agent.role || 'unknown');
            const domain = escapeHtml(agent.domain || 'generic');
            const notes = escapeHtml(agent.notes || 'No notes');
            return `
                <div class="bg-gray-700/40 rounded-lg p-3 border border-gray-700">
                    <div class="flex justify-between items-center">
                        <div class="font-semibold text-gray-200">${name}</div>
                        <div class="text-xs text-gray-500">${escapeHtml(created)}</div>
                    </div>
                    <div class="text-xs text-gray-400 mt-1">${description}</div>
                    <div class="text-xs text-gray-500 mt-2">Role: ${role} | Domain: ${domain}</div>
                    <div class="text-xs text-gray-500 mt-1">${notes}</div>
                </div>
            `;
        }).join('');
    },

    async submitFoundry() {
        const templateElem = document.getElementById('foundry-template');
        const nameElem = document.getElementById('foundry-name');
        const domainElem = document.getElementById('foundry-domain');
        const notesElem = document.getElementById('foundry-notes');
        const resultElem = document.getElementById('foundry-result');

        const templateId = templateElem ? templateElem.value : '';
        const name = nameElem ? nameElem.value.trim() : '';
        const domain = domainElem ? domainElem.value.trim() : '';
        const notes = notesElem ? notesElem.value.trim() : '';

        if (!templateId) {
            if (resultElem) {
                resultElem.textContent = 'Select a template to create an agent stub.';
            }
            return;
        }

        try {
            const response = await fetch('/api/agents/foundry', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    template_id: templateId,
                    name: name || 'Agent Stub',
                    domain: domain || 'generic',
                    notes
                })
            });

            if (!response.ok) {
                throw new Error('Failed to create agent stub');
            }

            const data = await response.json();
            const saved = data.item || {};
            if (resultElem) {
                resultElem.textContent = saved.path ? `Saved: ${saved.path}` : 'Agent stub created.';
            }
            if (nameElem) nameElem.value = '';
            if (domainElem) domainElem.value = '';
            if (notesElem) notesElem.value = '';
            await this.loadFoundryHistory();
        } catch (error) {
            console.error('[IMMUNOS] Error creating agent stub:', error);
            if (resultElem) {
                resultElem.textContent = 'Failed to create agent stub.';
            }
        }
    },

    async loadThymusStatus() {
        try {
            const response = await fetch('/api/thymus/status');
            if (!response.ok) throw new Error('Failed to load thymus status');
            const data = await response.json();
            this.updateThymusStatus(data);
        } catch (error) {
            console.error('[IMMUNOS] Error loading thymus status:', error);
        }
    },

    updateThymusStatus(status) {
        const statusElem = document.getElementById('thymus-queue-status');
        if (!statusElem) return;
        const label = status.paused ? 'paused' : (status.running ? 'running' : 'idle');
        statusElem.textContent = label;
    },

    async controlThymusQueue(action) {
        try {
            const response = await fetch('/api/thymus/control', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action })
            });
            if (!response.ok) throw new Error('Failed to control queue');
            const data = await response.json();
            this.updateThymusStatus(data);
        } catch (error) {
            console.error('[IMMUNOS] Error controlling thymus queue:', error);
        }
    },

    renderThymusQueue() {
        const container = document.getElementById('thymus-queue');
        const countElem = document.getElementById('thymus-queue-count');
        if (!container) return;

        if (countElem) {
            countElem.textContent = `${this.thymusQueue.length} entries`;
        }

        if (!this.thymusQueue.length) {
            container.innerHTML = '<div class="text-center text-gray-500 text-sm py-6">No queued items yet.</div>';
            return;
        }

        container.innerHTML = this.thymusQueue.map(item => {
            const domain = escapeHtml(item.domain || 'unknown');
            const dataset = escapeHtml(item.dataset_path || 'N/A');
            const samplesValue = item.sample_count;
            const samples = samplesValue === 0 || samplesValue ? `Samples: ${samplesValue}` : 'Samples: —';
            const notes = escapeHtml(item.notes || 'No notes');
            const created = item.created_at ? new Date(item.created_at).toLocaleString() : '';
            const status = escapeHtml(item.status || 'queued');

            return `
                <div class="bg-gray-700/40 rounded-lg p-3 border border-gray-700">
                    <div class="flex justify-between items-center mb-1">
                        <span class="text-sm font-semibold">${domain}</span>
                        <span class="text-xs text-gray-400">${escapeHtml(created)}</span>
                    </div>
                    <div class="text-xs text-gray-300">Path: <span class="font-mono">${dataset}</span></div>
                    <div class="text-xs text-gray-400 mt-1">${escapeHtml(samples)}</div>
                    <div class="text-xs text-gray-500 mt-1">Status: ${status}</div>
                    <div class="text-xs text-gray-500 mt-1">${notes}</div>
                </div>
            `;
        }).join('');

        this.renderThymusHistory();
    },

    renderThymusHistory() {
        const tbody = document.getElementById('training-history-tbody');
        if (!tbody) return;

        const history = (this.thymusQueue || []).filter(item => item.status && item.status !== 'queued');
        if (!history.length) {
            tbody.innerHTML = '<tr><td colspan="6" class="text-center text-gray-500 py-4">No training history</td></tr>';
            return;
        }

        tbody.innerHTML = history.map(item => {
            const domain = escapeHtml(item.domain || 'unknown');
            const dataset = escapeHtml(item.dataset_label || item.dataset_path || 'N/A');
            const detectors = escapeHtml(item.detectors ?? '—');
            const accuracy = typeof item.accuracy === 'number' ? `${item.accuracy}%` : '—';
            const started = item.started_at ? new Date(item.started_at) : null;
            const completed = item.completed_at ? new Date(item.completed_at) : null;
            const duration = started && completed ? `${Math.max(1, Math.round((completed - started) / 1000))}s` : '—';
            const dateLabel = completed ? completed.toLocaleString() : (started ? started.toLocaleString() : '');

            return `
                <tr class="border-b border-gray-700/50">
                    <td class="py-2">${domain}</td>
                    <td class="py-2">${dataset}</td>
                    <td class="py-2">${detectors}</td>
                    <td class="py-2">${escapeHtml(accuracy)}</td>
                    <td class="py-2">${escapeHtml(duration)}</td>
                    <td class="py-2">${escapeHtml(dateLabel)}</td>
                </tr>
            `;
        }).join('');
    },

    renderKbList() {
        const list = document.getElementById('kb-list');
        if (!list) return;

        if (!this.kbPages.length) {
            list.innerHTML = '<div class="text-center text-gray-500 text-sm py-6">No KB pages found.</div>';
            return;
        }

        list.innerHTML = this.kbPages.map(page => {
            const name = escapeHtml(page.name || '');
            const title = escapeHtml(page.title || page.name || 'Untitled');
            return `
                <button class="w-full text-left px-3 py-2 rounded bg-gray-700/40 hover:bg-gray-700 text-sm" data-kb="${name}">
                    ${title}
                </button>
            `;
        }).join('');

        list.querySelectorAll('button[data-kb]').forEach(btn => {
            btn.addEventListener('click', () => {
                const name = btn.getAttribute('data-kb');
                if (name) {
                    this.loadKbPage(name);
                }
            });
        });
    },

    async loadKbPage(name) {
        const content = document.getElementById('kb-content');
        const title = document.getElementById('kb-title');
        if (content) {
            content.textContent = 'Loading...';
        }

        try {
            const response = await fetch(`/api/kb/page?name=${encodeURIComponent(name)}`);
            if (!response.ok) throw new Error('Failed to load KB page');
            const data = await response.json();
            if (title) title.textContent = data.title || name;
            if (content) content.textContent = data.content || '';
        } catch (error) {
            if (content) content.textContent = `Error: ${error.message}`;
        }
    },

    renderIssues() {
        const summary = this.issuesSummary || {};
        const activeElem = document.getElementById('issues-active');
        const overdueElem = document.getElementById('issues-overdue');
        const dueTodayElem = document.getElementById('issues-due-today');
        const dueWeekElem = document.getElementById('issues-due-week');

        if (activeElem) activeElem.textContent = summary.active ?? 0;
        if (overdueElem) overdueElem.textContent = summary.overdue ?? 0;
        if (dueTodayElem) dueTodayElem.textContent = summary.due_today ?? 0;
        if (dueWeekElem) dueWeekElem.textContent = summary.due_this_week ?? 0;

        const list = document.getElementById('issues-list');
        if (!list) return;

        if (!this.issuesList || this.issuesList.length === 0) {
            list.innerHTML = `
                <div class="text-center text-gray-500 text-sm py-6 col-span-2">
                    No issues tracked yet.
                </div>
            `;
            return;
        }

        list.innerHTML = this.issuesList.map(item => {
            const title = escapeHtml(item.title || 'Untitled');
            const status = escapeHtml(item.status || 'unknown');
            const priority = escapeHtml(item.priority || 'medium');
            const due = item.due ? new Date(item.due).toLocaleDateString() : '—';
            const tags = (item.tags || []).slice(0, 3).map(tag => `#${escapeHtml(tag)}`).join(' ');
            const itemId = escapeHtml(item.id || '');

            return `
                <div class="bg-gray-700/40 rounded-lg p-4 border border-gray-700">
                    <div class="flex justify-between items-center mb-2">
                        <div class="font-semibold text-sm">${title}</div>
                        <span class="text-xs text-gray-400">${itemId}</span>
                    </div>
                    <div class="flex justify-between text-xs text-gray-400 mb-2">
                        <span>Status: <span class="text-gray-200">${status}</span></span>
                        <span>Priority: <span class="text-gray-200">${priority}</span></span>
                    </div>
                    <div class="flex justify-between text-xs text-gray-500">
                        <span>Due: ${escapeHtml(due)}</span>
                        <span>${tags}</span>
                    </div>
                </div>
            `;
        }).join('');
    },

    setupChat() {
        const sendBtn = document.getElementById('chat-send');
        const input = document.getElementById('chat-input');
        const saveBtn = document.getElementById('orchestrator-save');
        const ollamaStart = document.getElementById('chat-ollama-start');
        const ollamaStop = document.getElementById('chat-ollama-stop');
        const offlineStart = document.getElementById('chat-ollama-offline-start');

        if (sendBtn && input) {
            sendBtn.addEventListener('click', () => this.sendChatMessage());
            input.addEventListener('keydown', (event) => {
                if (event.key === 'Enter' && (event.metaKey || event.ctrlKey)) {
                    event.preventDefault();
                    this.sendChatMessage();
                }
            });
        }

        if (saveBtn) {
            saveBtn.addEventListener('click', () => this.saveOrchestratorConfig());
        }

        if (ollamaStart) {
            ollamaStart.addEventListener('click', () => this.startOllama());
        }

        if (ollamaStop) {
            ollamaStop.addEventListener('click', () => this.stopOllama());
        }

        if (offlineStart) {
            offlineStart.addEventListener('click', () => this.startOllama());
        }
    },

    appendChatMessage(role, content, meta = {}) {
        const container = document.getElementById('chat-messages');
        if (!container) return;

        const placeholder = container.querySelector('.text-center.text-gray-500');
        if (placeholder) {
            container.innerHTML = '';
        }

        const bubble = document.createElement('div');
        const isUser = role === 'user';
        bubble.className = `rounded-lg p-3 ${isUser ? 'bg-blue-600/20 border border-blue-600/40' : 'bg-gray-700/50 border border-gray-700'}`;
        const metaText = meta.model ? `${meta.model}${meta.reason ? ` • ${meta.reason}` : ''}` : '';
        if (metaText) {
            const metaLine = document.createElement('div');
            metaLine.className = 'text-xs text-gray-400 mb-1';
            metaLine.textContent = metaText;
            bubble.appendChild(metaLine);
        }
        const contentEl = document.createElement('div');
        contentEl.className = 'text-sm text-gray-100 whitespace-pre-wrap';
        contentEl.textContent = content;
        bubble.appendChild(contentEl);

        container.appendChild(bubble);
        container.scrollTop = container.scrollHeight;
    },

    setupQuickChat() {
        const sendBtn = document.getElementById('quick-chat-send');
        const input = document.getElementById('quick-chat-input');

        if (sendBtn && input) {
            sendBtn.addEventListener('click', () => this.sendQuickChatMessage());
            input.addEventListener('keydown', (event) => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    this.sendQuickChatMessage();
                }
            });
        }
    },

    appendQuickChatMessage(role, content) {
        const feed = document.getElementById('quick-chat-feed');
        if (!feed) return;

        const placeholder = feed.querySelector('.text-center.text-gray-500');
        if (placeholder) {
            feed.innerHTML = '';
        }

        const bubble = document.createElement('div');
        bubble.className = role === 'user'
            ? 'rounded-lg bg-blue-600/20 border border-blue-600/40 p-3 text-sm'
            : 'rounded-lg bg-gray-700/50 border border-gray-700 p-3 text-sm';
        bubble.textContent = content;
        feed.appendChild(bubble);
        feed.scrollTop = feed.scrollHeight;
    },

    async sendQuickChatMessage() {
        const input = document.getElementById('quick-chat-input');
        if (!input) return;
        const message = input.value.trim();
        if (!message) return;

        input.value = '';
        this.appendQuickChatMessage('user', message);

        try {
            const response = await fetch('/api/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    message,
                    task_type: 'analysis',
                    task_classification: 'routine',
                    mode: 'orchestrator'
                })
            });

            if (!response.ok) {
                throw new Error('Quick chat failed');
            }

            const data = await response.json();
            this.appendQuickChatMessage('assistant', data.response || 'No response');
        } catch (error) {
            this.appendQuickChatMessage('assistant', `Error: ${error.message}`);
        }
    },

    async sendChatMessage() {
        const input = document.getElementById('chat-input');
        if (!input) return;

        const message = input.value.trim();
        if (!message) return;

        const modeElem = document.getElementById('chat-mode');
        const domainElem = document.getElementById('chat-domain');
        const connectivityElem = document.getElementById('orchestrator-connectivity');

        const mode = modeElem ? modeElem.value : 'orchestrator';
        const domain = domainElem ? domainElem.value : '';
        const connectivity = connectivityElem ? connectivityElem.value : null;

        input.value = '';
        this.appendChatMessage('user', message);

        try {
            const payload = {
                message,
                task_type: domain || 'analysis',
                task_classification: 'routine',
                mode,
                domain: domain || null
            };

            if (mode === 'orchestrator' && connectivity) {
                payload.connectivity = connectivity;
            }

            const response = await fetch('/api/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });

            if (!response.ok) {
                throw new Error('Chat request failed');
            }

            const data = await response.json();
            this.appendChatMessage('assistant', data.response || 'No response', {
                model: data.model || 'unknown',
                reason: data.routing_reason || ''
            });
        } catch (error) {
            this.appendChatMessage('assistant', `Error: ${error.message}`);
        }
    },

    updateOrchestratorUI() {
        const config = this.orchestratorConfig?.config;
        const active = this.orchestratorConfig?.active_backend;
        if (!config) return;

        const connectivityElem = document.getElementById('orchestrator-connectivity');
        if (connectivityElem) connectivityElem.value = config.connectivity || 'online';

        const onlineProvider = document.getElementById('orchestrator-online-provider');
        const onlineModel = document.getElementById('orchestrator-online-model');
        const onlineUrl = document.getElementById('orchestrator-online-url');

        if (onlineProvider) onlineProvider.value = config.online_backend?.provider || 'claude_code';
        if (onlineModel) onlineModel.value = config.online_backend?.model || '';
        if (onlineUrl) onlineUrl.value = config.online_backend?.base_url || '';

        const offlineProvider = document.getElementById('orchestrator-offline-provider');
        const offlineModel = document.getElementById('orchestrator-offline-model');
        const offlineUrl = document.getElementById('orchestrator-offline-url');

        if (offlineProvider) offlineProvider.value = config.offline_backend?.provider || 'ollama';
        if (offlineModel) offlineModel.value = config.offline_backend?.model || '';
        if (offlineUrl) offlineUrl.value = config.offline_backend?.base_url || '';

        const activeElem = document.getElementById('orchestrator-active');
        if (activeElem && active) {
            let label = `${active.provider || 'orchestrator'}${active.model ? `:${active.model}` : ''} (${active.connectivity})`;
            if (active.fallback_reason) {
                label += ` fallback:${active.fallback_reason}`;
            }
            activeElem.textContent = label;
        }

        const fallbackElem = document.getElementById('orchestrator-fallback');
        if (fallbackElem) {
            if (active?.fallback_reason) {
                fallbackElem.textContent = `Orchestrator fallback active (${active.fallback_reason}). Using offline backend until online credentials are available.`;
                fallbackElem.classList.remove('hidden');
            } else {
                fallbackElem.classList.add('hidden');
            }
        }
    },

    resolveApiKeyEnv(provider) {
        switch (provider) {
            case 'chatgpt':
                return 'OPENAI_API_KEY';
            case 'openrouter':
                return 'OPENROUTER_API_KEY';
            case 'claude_code':
                return 'ANTHROPIC_AUTH_TOKEN';
            default:
                return '';
        }
    },

    async saveOrchestratorConfig() {
        const connectivityElem = document.getElementById('orchestrator-connectivity');
        const onlineProvider = document.getElementById('orchestrator-online-provider');
        const onlineModel = document.getElementById('orchestrator-online-model');
        const onlineUrl = document.getElementById('orchestrator-online-url');
        const offlineProvider = document.getElementById('orchestrator-offline-provider');
        const offlineModel = document.getElementById('orchestrator-offline-model');
        const offlineUrl = document.getElementById('orchestrator-offline-url');

        const payload = {
            connectivity: connectivityElem ? connectivityElem.value : 'online',
            online_backend: {
                provider: onlineProvider ? onlineProvider.value : 'claude_code',
                model: onlineModel ? onlineModel.value : '',
                base_url: onlineUrl ? onlineUrl.value : '',
                api_key_env: this.resolveApiKeyEnv(onlineProvider ? onlineProvider.value : '')
            },
            offline_backend: {
                provider: offlineProvider ? offlineProvider.value : 'ollama',
                model: offlineModel ? offlineModel.value : '',
                base_url: offlineUrl ? offlineUrl.value : '',
                api_key_env: this.resolveApiKeyEnv(offlineProvider ? offlineProvider.value : '')
            }
        };

        try {
            const response = await fetch('/api/orchestrator/config', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });

            if (!response.ok) {
                throw new Error('Failed to save orchestrator config');
            }

            const data = await response.json();
            this.orchestratorConfig = data;
            this.updateOrchestratorUI();
        } catch (error) {
            this.addEventToFeed('orchestrator_update', { timestamp: new Date().toISOString() });
            console.error('[IMMUNOS] Error saving orchestrator config:', error);
        }
    },

    // ========================================================================
    // EVENT LISTENERS
    // ========================================================================

    setupEventListeners() {
        // Ollama server controls
        const btnOllamaStart = document.getElementById('btn-ollama-start');
        const btnOllamaStop = document.getElementById('btn-ollama-stop');
        const simulateBtn = document.getElementById('simulate-detection');

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

        if (simulateBtn) {
            simulateBtn.addEventListener('click', () => {
                this.simulateDetection();
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
            this.loadOllamaStatus();
        }, 5000);

        // Update domain stats every 30 seconds
        setInterval(() => {
            this.loadDomainStatus();
            this.loadIssues();
            this.loadThymusQueue();
            this.loadThymusStatus();
            this.loadSpleenSummary();
            this.loadFoundryHistory();
        }, 30000);

        // Refresh orchestrator status every 60 seconds
        setInterval(() => {
            this.loadOrchestratorStatus();
        }, 60000);

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

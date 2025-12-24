/**
 * IMMUNOS Dashboard - Main JavaScript Controller
 *
 * Handles WebSocket connections, API calls, and real-time updates
 */

class ImmunosApp {
    constructor() {
        this.socket = null;
        this.connected = false;
        this.refreshInterval = null;
    }

    /**
     * Initialize the application
     */
    init() {
        console.log('🧬 Initializing IMMUNOS Dashboard...');

        // Connect to WebSocket
        this.connectWebSocket();

        // Setup auto-refresh
        this.startAutoRefresh();

        // Setup event listeners
        this.setupEventListeners();

        console.log('✓ IMMUNOS Dashboard initialized');
    }

    /**
     * Connect to WebSocket server
     */
    connectWebSocket() {
        try {
            this.socket = io();

            this.socket.on('connect', () => {
                console.log('✓ WebSocket connected');
                this.connected = true;
                this.updateConnectionStatus(true);

                // Subscribe to events
                this.socket.emit('subscribe', {
                    events: ['scan_progress', 'scan_completed', 'scan_error', 'activity', 'health_updated']
                });
            });

            this.socket.on('disconnect', () => {
                console.log('✗ WebSocket disconnected');
                this.connected = false;
                this.updateConnectionStatus(false);
            });

            // Scan progress updates
            this.socket.on('scan_progress', (data) => {
                this.handleScanProgress(data);
            });

            // Scan completed
            this.socket.on('scan_completed', (data) => {
                this.handleScanCompleted(data);
            });

            // Scan error
            this.socket.on('scan_error', (data) => {
                this.handleScanError(data);
            });

            // Activity updates
            this.socket.on('activity', (data) => {
                this.handleActivity(data);
            });

            // Health updates
            this.socket.on('health_updated', (data) => {
                this.handleHealthUpdate(data);
            });

        } catch (error) {
            console.error('WebSocket connection error:', error);
            this.updateConnectionStatus(false);
        }
    }

    /**
     * Update connection status indicator
     */
    updateConnectionStatus(connected) {
        const statusDot = document.querySelector('.status-dot');
        const statusText = document.querySelector('.status-text');

        if (statusDot) {
            if (connected) {
                statusDot.classList.remove('offline');
                statusText.textContent = 'Connected';
            } else {
                statusDot.classList.add('offline');
                statusText.textContent = 'Disconnected';
            }
        }
    }

    /**
     * Setup event listeners
     */
    setupEventListeners() {
        // No global event listeners needed yet
    }

    /**
     * Start auto-refresh for dashboard data
     */
    startAutoRefresh() {
        // Refresh every 30 seconds
        this.refreshInterval = setInterval(() => {
            if (this.connected) {
                this.refreshDashboard();
            }
        }, 30000);
    }

    /**
     * Refresh dashboard data
     */
    async refreshDashboard() {
        try {
            // Only refresh if on dashboard page
            if (window.location.pathname === '/') {
                await this.loadDashboardData();
            }
        } catch (error) {
            console.error('Dashboard refresh error:', error);
        }
    }

    /**
     * Load dashboard data (called from dashboard page)
     */
    async loadDashboardData() {
        try {
            // Load health status
            await this.loadHealth();

            // Load stats
            await this.loadStats();

            // Load component status
            await this.loadComponentStatus();

            // Load activity
            await this.loadActivity();

        } catch (error) {
            console.error('Error loading dashboard data:', error);
        }
    }

    /**
     * Load health status
     */
    async loadHealth() {
        try {
            const response = await fetch('/api/health');
            const result = await response.json();

            if (result.success && result.data) {
                this.updateHealthDisplay(result.data);
            }
        } catch (error) {
            console.error('Error loading health:', error);
        }
    }

    /**
     * Update health display
     */
    updateHealthDisplay(health) {
        const banner = document.getElementById('threat-banner');
        const title = document.getElementById('threat-title');
        const subtitle = document.getElementById('threat-subtitle');
        const scoreValue = document.getElementById('score-value');
        const scoreFill = document.getElementById('score-fill');

        if (!banner) return;

        // Update banner class
        banner.className = 'threat-banner';
        if (health.status === 'warning') {
            banner.classList.add('warning');
        } else if (health.status === 'critical') {
            banner.classList.add('critical');
        }

        // Update title
        if (title) {
            if (health.status === 'secure') {
                title.textContent = 'System Secure';
            } else if (health.status === 'warning') {
                title.textContent = 'Warnings Detected';
            } else if (health.status === 'critical') {
                title.textContent = 'Critical Issues';
            } else {
                title.textContent = 'System Status Unknown';
            }
        }

        // Update subtitle
        if (subtitle) {
            if (health.counts) {
                const parts = [];
                if (health.counts.critical > 0) parts.push(`${health.counts.critical} critical`);
                if (health.counts.warning > 0) parts.push(`${health.counts.warning} warnings`);
                if (health.counts.info > 0) parts.push(`${health.counts.info} info`);

                subtitle.textContent = parts.length > 0
                    ? parts.join(', ')
                    : 'No issues detected';
            }
        }

        // Update score
        if (scoreValue) {
            scoreValue.textContent = health.overall_health || 100;
        }

        // Update score circle
        if (scoreFill) {
            const score = health.overall_health || 100;
            const circumference = 314;
            const offset = circumference - (circumference * score / 100);
            scoreFill.style.strokeDashoffset = offset;

            // Update color based on score
            if (score >= 80) {
                scoreFill.style.stroke = '#27ae60';
            } else if (score >= 50) {
                scoreFill.style.stroke = '#f39c12';
            } else {
                scoreFill.style.stroke = '#e74c3c';
            }
        }
    }

    /**
     * Load stats
     */
    async loadStats() {
        try {
            const response = await fetch('/api/stats');
            const result = await response.json();

            if (result.success && result.data) {
                this.updateStatsDisplay(result.data);
            }
        } catch (error) {
            console.error('Error loading stats:', error);
        }
    }

    /**
     * Update stats display
     */
    updateStatsDisplay(stats) {
        const criticalEl = document.getElementById('stat-critical');
        const warningsEl = document.getElementById('stat-warnings');
        const infoEl = document.getElementById('stat-info');
        const fixedEl = document.getElementById('stat-fixed');

        if (criticalEl) criticalEl.textContent = stats.critical || 0;
        if (warningsEl) warningsEl.textContent = stats.warnings || 0;
        if (infoEl) infoEl.textContent = stats.info || 0;
        if (fixedEl) fixedEl.textContent = stats.fixed || 0;
    }

    /**
     * Load component status
     */
    async loadComponentStatus() {
        try {
            const response = await fetch('/api/status');
            const result = await response.json();

            if (result.success && result.data) {
                this.updateComponentStatus(result.data.components);
            }
        } catch (error) {
            console.error('Error loading component status:', error);
        }
    }

    /**
     * Update component status
     */
    updateComponentStatus(components) {
        // NK Cell
        if (components.nk_cell) {
            const nk = components.nk_cell;
            const statusEl = document.getElementById('nk-status');
            const lastScanEl = document.getElementById('nk-last-scan');

            if (statusEl) {
                const badge = statusEl.querySelector('.status-badge');
                if (badge) {
                    badge.className = 'status-badge ' + nk.status;
                    badge.textContent = nk.status.charAt(0).toUpperCase() + nk.status.slice(1);
                }
            }

            if (lastScanEl && nk.last_run) {
                lastScanEl.textContent = this.formatDate(nk.last_run);
            }

            // Update anomaly counts
            if (nk.anomalies) {
                const anomaliesEl = document.getElementById('nk-anomalies');
                if (anomaliesEl) {
                    const highEl = anomaliesEl.querySelector('.anomaly-item.high .anomaly-count');
                    const mediumEl = anomaliesEl.querySelector('.anomaly-item.medium .anomaly-count');
                    const lowEl = anomaliesEl.querySelector('.anomaly-item.low .anomaly-count');

                    if (highEl) highEl.textContent = nk.anomalies.high || 0;
                    if (mediumEl) mediumEl.textContent = nk.anomalies.medium || 0;
                    if (lowEl) lowEl.textContent = nk.anomalies.low || 0;
                }
            }
        }

        // T Cell Memory
        if (components.t_cell) {
            const memoryCountEl = document.getElementById('memory-count');
            if (memoryCountEl) {
                memoryCountEl.textContent = components.t_cell.memories?.active || 0;
            }
        }

        // Todo
        if (components.todo) {
            const activeEl = document.getElementById('todo-active');
            const completedEl = document.getElementById('todo-completed');

            if (activeEl) activeEl.textContent = components.todo.stats?.active || 0;
            if (completedEl) completedEl.textContent = components.todo.stats?.completed || 0;
        }

        // Snapshots
        if (components.snapshots) {
            const countEl = document.getElementById('snapshot-count');
            if (countEl) countEl.textContent = components.snapshots.count || 0;
        }
    }

    /**
     * Load activity log
     */
    async loadActivity() {
        try {
            const response = await fetch('/api/activity?limit=10');
            const result = await response.json();

            if (result.success && result.data) {
                this.updateActivityDisplay(result.data);
            }
        } catch (error) {
            console.error('Error loading activity:', error);
        }
    }

    /**
     * Update activity display
     */
    updateActivityDisplay(activities) {
        const timelineEl = document.getElementById('activity-timeline');
        if (!timelineEl) return;

        if (!activities || activities.length === 0) {
            timelineEl.innerHTML = `
                <div class="activity-empty">
                    <div class="empty-icon">📭</div>
                    <p>No recent activity</p>
                </div>
            `;
            return;
        }

        // TODO: Render activity items
        // For now, keep the empty state
    }

    /**
     * Trigger NK Cell scan
     */
    async triggerScan(scanType = 'full') {
        try {
            this.showLoading(true);
            this.showToast('Starting scan...', 'info');

            const response = await fetch('/api/nk/scan', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    scan_type: scanType,
                    use_ollama: false,
                    trigger: 'manual'
                })
            });

            const result = await response.json();

            if (result.success) {
                this.showToast('Scan started successfully', 'success');
            } else {
                this.showToast(`Scan failed: ${result.error}`, 'error');
            }

            this.showLoading(false);

        } catch (error) {
            console.error('Error triggering scan:', error);
            this.showToast('Failed to start scan', 'error');
            this.showLoading(false);
        }
    }

    /**
     * Run memory decay
     */
    async runMemoryDecay() {
        try {
            this.showLoading(true);

            const response = await fetch('/api/memory/decay', {
                method: 'POST'
            });

            const result = await response.json();

            if (result.success) {
                const expired = result.data?.expired_count || 0;
                this.showToast(`Memory decay complete. ${expired} memories expired.`, 'success');
            } else {
                this.showToast(`Decay failed: ${result.error}`, 'error');
            }

            this.showLoading(false);

        } catch (error) {
            console.error('Error running decay:', error);
            this.showToast('Failed to run memory decay', 'error');
            this.showLoading(false);
        }
    }

    /**
     * Handle scan progress updates from WebSocket
     */
    handleScanProgress(data) {
        console.log('Scan progress:', data);
        // Update scan status in UI
        this.showToast(`Scan progress: ${data.percent}%`, 'info');
    }

    /**
     * Handle scan completion from WebSocket
     */
    handleScanCompleted(data) {
        console.log('Scan completed:', data);
        this.showToast(`Scan completed! Found ${data.total} anomalies`, 'success');

        // Refresh dashboard
        this.loadDashboardData();
    }

    /**
     * Handle scan error from WebSocket
     */
    handleScanError(data) {
        console.error('Scan error:', data);
        this.showToast(`Scan failed: ${data.error}`, 'error');
    }

    /**
     * Handle activity updates from WebSocket
     */
    handleActivity(data) {
        console.log('Activity:', data);
        // Could update activity timeline here
    }

    /**
     * Handle health updates from WebSocket
     */
    handleHealthUpdate(data) {
        console.log('Health update:', data);
        this.loadHealth();
    }

    /**
     * Show toast notification
     */
    showToast(message, type = 'info') {
        const container = document.getElementById('toast-container');
        if (!container) return;

        const toast = document.createElement('div');
        toast.className = `toast ${type}`;
        const icon = document.createElement('div');
        icon.className = 'toast-icon';
        icon.textContent = this.getToastIcon(type);

        const text = document.createElement('div');
        text.className = 'toast-message';
        text.textContent = message;

        toast.appendChild(icon);
        toast.appendChild(text);

        container.appendChild(toast);

        // Auto-remove after 5 seconds
        setTimeout(() => {
            toast.remove();
        }, 5000);
    }

    /**
     * Get toast icon based on type
     */
    getToastIcon(type) {
        const icons = {
            success: '✅',
            error: '❌',
            warning: '⚠️',
            info: 'ℹ️'
        };
        return icons[type] || icons.info;
    }

    /**
     * Show/hide loading overlay
     */
    showLoading(show) {
        const overlay = document.getElementById('loading-overlay');
        if (overlay) {
            overlay.style.display = show ? 'flex' : 'none';
        }
    }

    /**
     * Format date for display
     */
    formatDate(dateStr) {
        if (!dateStr) return 'Never';

        const date = new Date(dateStr);
        const now = new Date();
        const diff = now - date;

        // Less than 1 minute
        if (diff < 60000) {
            return 'Just now';
        }

        // Less than 1 hour
        if (diff < 3600000) {
            const minutes = Math.floor(diff / 60000);
            return `${minutes} minute${minutes > 1 ? 's' : ''} ago`;
        }

        // Less than 1 day
        if (diff < 86400000) {
            const hours = Math.floor(diff / 3600000);
            return `${hours} hour${hours > 1 ? 's' : ''} ago`;
        }

        // Format as date
        return date.toLocaleDateString();
    }
}

// Export for use in templates
window.ImmunosApp = ImmunosApp;

// IMMUNOS Minimal Dashboard - Chat & Model Control
// State management
let socket = null;
let currentModel = '';
let isConnected = false;

// DOM elements
const elements = {
    connectionStatus: document.getElementById('connectionStatus'),
    connectionText: document.getElementById('connectionText'),
    ollamaStatus: document.getElementById('ollamaStatus'),
    ollamaText: document.getElementById('ollamaText'),
    startOllama: document.getElementById('startOllama'),
    stopOllama: document.getElementById('stopOllama'),
    modelSelect: document.getElementById('modelSelect'),
    modelsList: document.getElementById('modelsList'),
    chatMessages: document.getElementById('chatMessages'),
    chatInput: document.getElementById('chatInput'),
    sendButton: document.getElementById('sendButton'),
    kbContent: document.getElementById('kbContent')
};

// Selected KB files for context
let selectedKBFiles = new Set();

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    initializeWebSocket();
    initializeEventListeners();
    checkOllamaStatus();
    loadModels();
    loadKnowledgeBase();
    loadProjectList();
    checkRecoveryStatus();
});

// WebSocket Connection
function initializeWebSocket() {
    socket = io();

    socket.on('connect', () => {
        isConnected = true;
        updateConnectionStatus(true);
    });

    socket.on('disconnect', () => {
        isConnected = false;
        updateConnectionStatus(false);
    });

    socket.on('chat_response', (data) => {
        addMessage('assistant', data.response, data.model);
    });

    socket.on('chat_error', (data) => {
        addMessage('error', data.error);
    });
}

// Event Listeners
function initializeEventListeners() {
    // Send button
    elements.sendButton.addEventListener('click', sendMessage);

    // Enter key to send
    elements.chatInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            sendMessage();
        }
    });

    // Model selection
    elements.modelSelect.addEventListener('change', (e) => {
        currentModel = e.target.value;
        console.log('Selected model:', currentModel);

        // Auto-highlight memory panel for T Cell (context-related) model
        const memoryPanel = document.getElementById('memoryPanel');
        if (memoryPanel) {
            if (currentModel.includes('qwen2.5:7b')) {
                memoryPanel.classList.add('ring-2', 'ring-blue-500');
                memoryPanel.querySelector('h2').textContent = 'T Cell Memory (Active)';
            } else {
                memoryPanel.classList.remove('ring-2', 'ring-blue-500');
                memoryPanel.querySelector('h2').textContent = 'T Cell Memory';
            }
        }
    });

    // Ollama controls
    elements.startOllama.addEventListener('click', startOllama);
    elements.stopOllama.addEventListener('click', stopOllama);
}

// Chat Functions
async function sendMessage() {
    const message = elements.chatInput.value.trim();
    if (!message) return;

    if (!currentModel) {
        addMessage('error', 'Please select a model first');
        return;
    }

    // Add user message to chat
    addMessage('user', message);

    // Clear input
    elements.chatInput.value = '';

    // Get KB context if any files selected
    const kbContext = await getSelectedKBContext();
    const fullMessage = kbContext ? `${kbContext}\n\n### User Question:\n${message}` : message;

    // Show context indicator
    const kbFilesCount = document.querySelectorAll('.kb-checkbox:checked').length;
    const uploadedCount = uploadedFiles.filter(f => f.selected).length;
    const totalContextFiles = kbFilesCount + uploadedCount;
    if (totalContextFiles > 0) {
        addMessage('assistant', `[Using ${totalContextFiles} context file(s): ${kbFilesCount} KB + ${uploadedCount} uploaded]`, 'System');
    }

    // Send to server
    try {
        const response = await fetch('/api/chat', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                message: fullMessage,
                model: currentModel
            })
        });
        const data = await response.json();
        if (data.error) {
            addMessage('error', data.error);
        }
    } catch (error) {
        addMessage('error', `Failed to send message: ${error.message}`);
    }
}

function addMessage(type, content, model = null) {
    const messagesContainer = elements.chatMessages;

    // Remove "no messages" placeholder
    if (messagesContainer.querySelector('.text-center')) {
        messagesContainer.innerHTML = '';
    }

    const messageDiv = document.createElement('div');
    messageDiv.className = 'message';

    if (type === 'user') {
        messageDiv.innerHTML = `
            <div class="flex justify-end">
                <div class="bg-blue-600 rounded-lg px-4 py-2 max-w-md">
                    <p class="text-sm">${escapeHtml(content)}</p>
                </div>
            </div>
        `;
    } else if (type === 'assistant') {
        messageDiv.innerHTML = `
            <div class="flex justify-start">
                <div class="bg-gray-700 rounded-lg px-4 py-2 max-w-md">
                    ${model ? `<p class="text-xs text-gray-400 mb-1">${escapeHtml(model)}</p>` : ''}
                    <p class="text-sm whitespace-pre-wrap">${escapeHtml(content)}</p>
                </div>
            </div>
        `;
    } else if (type === 'error') {
        messageDiv.innerHTML = `
            <div class="flex justify-center">
                <div class="bg-red-900 border border-red-700 rounded-lg px-4 py-2">
                    <p class="text-sm text-red-200">${escapeHtml(content)}</p>
                </div>
            </div>
        `;
    }

    messagesContainer.appendChild(messageDiv);
    messagesContainer.scrollTop = messagesContainer.scrollHeight;
}

// Ollama Functions
async function checkOllamaStatus() {
    try {
        // Check Ollama directly (faster, more reliable)
        const response = await fetch('http://localhost:11434/api/tags');
        if (response.ok) {
            updateOllamaStatus(true);
            return;
        }
        updateOllamaStatus(false);
    } catch (error) {
        updateOllamaStatus(false);
        console.error('Error checking Ollama status:', error);
    }
}

async function startOllama() {
    try {
        elements.startOllama.disabled = true;
        elements.startOllama.textContent = 'Starting...';

        // Try Flask API first
        const response = await fetch('/api/ollama/start', { method: 'POST' });
        const data = await response.json();

        if (data.success) {
            // Wait for Ollama to be ready
            await new Promise(resolve => setTimeout(resolve, 2000));
            await checkOllamaStatus();
            loadModels();
        } else {
            addMessage('error', `Failed to start Ollama: ${data.error || 'Unknown error'}`);
        }
    } catch (error) {
        addMessage('error', `Error starting Ollama: ${error.message}. Try running 'ollama serve' manually.`);
    } finally {
        elements.startOllama.disabled = false;
        elements.startOllama.textContent = 'Start Ollama';
    }
}

async function stopOllama() {
    try {
        elements.stopOllama.disabled = true;
        elements.stopOllama.textContent = 'Stopping...';

        const response = await fetch('/api/ollama/stop', { method: 'POST' });
        const data = await response.json();

        if (data.success) {
            updateOllamaStatus(false);
            elements.modelSelect.innerHTML = '<option value="">Ollama stopped</option>';
        } else {
            addMessage('error', `Failed to stop Ollama: ${data.error}`);
        }
    } catch (error) {
        addMessage('error', `Error stopping Ollama: ${error.message}`);
    } finally {
        elements.stopOllama.disabled = false;
        elements.stopOllama.textContent = 'Stop Ollama';
    }
}

// IMMUNOS model role mappings (Biological Immune System Metaphor)
// NK Cell = Quick Scanner | B Cell = Code Verifier | Dendritic = Research Analyst | T Cell = Memory Agent
const MODEL_ROLES = {
    'qwen2.5-coder': { role: 'B Cell', desc: 'Code Implementation' },
    'deepseek-r1': { role: 'Dendritic', desc: 'Research Analysis' },
    'qwen2.5:1.5b': { role: 'NK Cell', desc: 'Quick Tasks' },
    'qwen2.5:7b': { role: 'T Cell', desc: 'Memory & Context' },
    'llama': { role: 'General', desc: 'General Purpose' },
    'mistral': { role: 'General', desc: 'General Purpose' },
    'phi': { role: 'NK Cell', desc: 'Quick Tasks' },
    'gemma': { role: 'General', desc: 'General Purpose' }
};

function getModelRole(modelName) {
    // Check for exact match first
    if (MODEL_ROLES[modelName]) return MODEL_ROLES[modelName];

    // Check for partial matches
    for (const [key, value] of Object.entries(MODEL_ROLES)) {
        if (modelName.includes(key)) return value;
    }

    // Embedding models (not for chat)
    if (modelName.includes('embed') || modelName.includes('bge')) {
        return { role: 'Embed', desc: 'Embeddings Only' };
    }

    return { role: 'Agent', desc: 'General Agent' };
}

async function loadModels() {
    try {
        // Fetch directly from Ollama API
        const response = await fetch('http://localhost:11434/api/tags');
        const data = await response.json();

        // Clear both lists
        if (elements.modelsList) elements.modelsList.innerHTML = '';
        if (elements.modelSelect) elements.modelSelect.innerHTML = '<option value="">Select model for chat...</option>';

        if (data.models && data.models.length > 0) {
            // Filter chat models (exclude embedding models)
            const chatModels = data.models.filter(m =>
                !m.name.includes('embed') && !m.name.includes('bge')
            );

            chatModels.forEach(model => {
                const sizeGB = (model.size / 1e9).toFixed(1);
                const roleInfo = getModelRole(model.name);

                // Add to checkbox list
                if (elements.modelsList) {
                    const item = document.createElement('label');
                    item.className = 'flex items-center gap-2 bg-gray-700 rounded p-2 cursor-pointer hover:bg-gray-600';
                    item.innerHTML = `
                        <input type="checkbox" class="model-checkbox" data-model="${model.name}" checked>
                        <span class="text-xs">
                            <span class="text-blue-400 font-medium">[${roleInfo.role}]</span>
                            <span class="text-gray-300">${model.name.split(':')[0]}</span>
                            <span class="text-gray-500">(${sizeGB}GB)</span>
                        </span>
                    `;
                    elements.modelsList.appendChild(item);
                }

                // Add to dropdown
                if (elements.modelSelect) {
                    const option = document.createElement('option');
                    option.value = model.name;
                    option.textContent = `[${roleInfo.role}] ${model.name} (${sizeGB}GB)`;
                    elements.modelSelect.appendChild(option);
                }
            });

            if (chatModels.length === 0 && elements.modelsList) {
                elements.modelsList.innerHTML = '<p class="text-gray-500 text-sm col-span-4">No chat models available</p>';
            }
        } else {
            if (elements.modelsList) {
                elements.modelsList.innerHTML = '<p class="text-gray-500 text-sm col-span-4">No models available</p>';
            }
        }
    } catch (error) {
        console.error('Error loading models:', error);
        if (elements.modelsList) {
            elements.modelsList.innerHTML = '<p class="text-red-400 text-sm col-span-4">Error loading models</p>';
        }
    }
}

// Status Updates
function updateConnectionStatus(connected) {
    if (connected) {
        elements.connectionStatus.className = 'w-3 h-3 rounded-full bg-green-500';
        elements.connectionText.textContent = 'Connected';
        elements.connectionText.className = 'text-sm text-green-400';
    } else {
        elements.connectionStatus.className = 'w-3 h-3 rounded-full bg-red-500';
        elements.connectionText.textContent = 'Disconnected';
        elements.connectionText.className = 'text-sm text-red-400';
    }
}

function updateOllamaStatus(running) {
    if (running) {
        elements.ollamaStatus.className = 'w-3 h-3 rounded-full bg-green-500';
        elements.ollamaText.textContent = 'Ollama Running';
        elements.ollamaText.className = 'text-sm text-green-400';
    } else {
        elements.ollamaStatus.className = 'w-3 h-3 rounded-full bg-red-500';
        elements.ollamaText.textContent = 'Ollama Stopped';
        elements.ollamaText.className = 'text-sm text-red-400';
    }
}

// Utility Functions
function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// Knowledge Base & Project Functions
let kbEnabled = true;
let uploadedFiles = [];
let projects = JSON.parse(localStorage.getItem('immunos_projects') || '{}');

// File Upload
function handleFileUpload(event) {
    const files = event.target.files;
    processFiles(files);
}

function handleFileDrop(event) {
    event.preventDefault();
    const files = event.dataTransfer.files;
    processFiles(files);
}

function processFiles(files) {
    for (const file of files) {
        const reader = new FileReader();
        reader.onload = (e) => {
            uploadedFiles.push({
                name: file.name,
                content: e.target.result,
                selected: true
            });
            renderUploadedFiles();
            updateKBCount();
        };
        reader.readAsText(file);
    }
}

function renderUploadedFiles() {
    const container = document.getElementById('uploadedFiles');
    const list = document.getElementById('uploadedFilesList');
    if (!container || !list) return;

    if (uploadedFiles.length === 0) {
        container.classList.add('hidden');
        return;
    }

    container.classList.remove('hidden');
    list.innerHTML = uploadedFiles.map((f, i) => `
        <label class="flex items-center gap-2 bg-gray-700 rounded p-2 cursor-pointer hover:bg-gray-600">
            <input type="checkbox" class="uploaded-checkbox" data-index="${i}" ${f.selected ? 'checked' : ''} onchange="toggleUploadedFile(${i})">
            <span class="text-sm truncate">📄 ${f.name}</span>
            <button onclick="removeUploadedFile(${i})" class="text-red-400 hover:text-red-300 ml-auto">&times;</button>
        </label>
    `).join('');
}

function toggleUploadedFile(index) {
    uploadedFiles[index].selected = !uploadedFiles[index].selected;
    updateKBCount();
}

function removeUploadedFile(index) {
    uploadedFiles.splice(index, 1);
    renderUploadedFiles();
    updateKBCount();
}

// Project Management
function loadProjectList() {
    const select = document.getElementById('projectSelect');
    if (!select) return;

    select.innerHTML = '<option value="">Select Project...</option>';

    // Sort by last used (most recent first)
    const sortedProjects = Object.entries(projects).sort((a, b) => {
        const aTime = a[1].lastUsed ? new Date(a[1].lastUsed) : new Date(0);
        const bTime = b[1].lastUsed ? new Date(b[1].lastUsed) : new Date(0);
        return bTime - aTime;
    });

    sortedProjects.forEach(([name, project]) => {
        const opt = document.createElement('option');
        opt.value = name;
        // Show project name and abbreviated working directory
        const workDir = project.workingDirectory || '';
        const shortDir = workDir.length > 30 ? '...' + workDir.slice(-27) : workDir;
        opt.textContent = shortDir ? `${name} (${shortDir})` : name;
        opt.title = workDir; // Full path on hover
        select.appendChild(opt);
    });
}

function saveProject() {
    const name = prompt('Project name:');
    if (!name) return;

    const workingDir = prompt('Working directory (local path or server URL):', '/Users/byron/projects');
    if (workingDir === null) return; // Cancelled

    const selectedModels = Array.from(document.querySelectorAll('.model-checkbox:checked')).map(cb => cb.dataset.model);
    const selectedKB = Array.from(document.querySelectorAll('.kb-checkbox:checked')).map(cb => cb.dataset.name);
    const activeModel = document.getElementById('modelSelect')?.value || '';

    projects[name] = {
        models: selectedModels,
        kb: selectedKB,
        activeModel: activeModel,
        uploadedFiles: uploadedFiles.map(f => ({ name: f.name, content: f.content })),
        workingDirectory: workingDir,
        createdAt: new Date().toISOString(),
        lastUsed: new Date().toISOString()
    };

    localStorage.setItem('immunos_projects', JSON.stringify(projects));
    loadProjectList();
    document.getElementById('projectSelect').value = name;
    addMessage('assistant', `Project "${name}" saved\nWorking Directory: ${workingDir}`, 'System');
}

function loadProject() {
    const select = document.getElementById('projectSelect');
    const name = select?.value;
    if (!name || !projects[name]) {
        addMessage('error', 'Select a project first');
        return;
    }

    const project = projects[name];

    // Restore models
    document.querySelectorAll('.model-checkbox').forEach(cb => {
        cb.checked = project.models.includes(cb.dataset.model);
    });

    // Restore KB
    document.querySelectorAll('.kb-checkbox').forEach(cb => {
        cb.checked = project.kb.includes(cb.dataset.name);
    });

    // Restore active model
    if (project.activeModel && document.getElementById('modelSelect')) {
        document.getElementById('modelSelect').value = project.activeModel;
        currentModel = project.activeModel;
    }

    // Restore uploaded files
    uploadedFiles = project.uploadedFiles || [];
    renderUploadedFiles();
    updateKBCount();

    // Update last used
    projects[name].lastUsed = new Date().toISOString();
    localStorage.setItem('immunos_projects', JSON.stringify(projects));

    // Show project info
    const workDir = project.workingDirectory || 'Not set';
    const lastUsed = project.lastUsed ? new Date(project.lastUsed).toLocaleString() : 'Unknown';
    addMessage('assistant', `Project "${name}" loaded\nWorking Directory: ${workDir}\nLast Used: ${lastUsed}`, 'System');
}

function newProject() {
    // Clear selections
    document.querySelectorAll('.model-checkbox').forEach(cb => cb.checked = true);
    document.querySelectorAll('.kb-checkbox').forEach(cb => cb.checked = false);
    uploadedFiles = [];
    renderUploadedFiles();
    updateKBCount();
    document.getElementById('projectSelect').value = '';
    addMessage('assistant', 'New project started', 'System');
}

async function loadKnowledgeBase() {
    try {
        const response = await fetch('/api/kb/index');
        const data = await response.json();

        if (!elements.kbContent) return;
        elements.kbContent.innerHTML = '';

        if (data.pages && data.pages.length > 0) {
            data.pages.forEach(page => {
                const item = document.createElement('label');
                item.className = 'flex items-center gap-3 bg-gray-700 rounded p-3 cursor-pointer hover:bg-gray-600';
                item.innerHTML = `
                    <input type="checkbox" class="kb-checkbox w-4 h-4" data-name="${page.name}" onchange="updateKBCount()">
                    <span class="text-xl">${getKBIcon(page.name)}</span>
                    <div class="flex-1">
                        <div class="text-sm font-medium">${page.title || page.name}</div>
                        <div class="text-xs text-gray-400">${page.name}.md</div>
                    </div>
                `;
                elements.kbContent.appendChild(item);
            });
        } else {
            elements.kbContent.innerHTML = '<p class="text-gray-500 text-sm">No KB files found</p>';
        }
    } catch (error) {
        console.error('Error loading KB:', error);
        if (elements.kbContent) {
            elements.kbContent.innerHTML = '<p class="text-red-400 text-sm">Error loading KB</p>';
        }
    }
}

function updateKBCount() {
    const kbCount = document.querySelectorAll('.kb-checkbox:checked').length;
    const uploadCount = uploadedFiles.filter(f => f.selected).length;
    const totalCount = kbCount + uploadCount;

    const kbCountEl = document.getElementById('kbCount');
    if (kbCountEl) kbCountEl.textContent = `${totalCount} files`;

    // Auto-enable KB if files selected
    const kbEnabledEl = document.getElementById('kbEnabled');
    if (kbEnabledEl && totalCount > 0) {
        kbEnabledEl.checked = true;
        kbEnabled = true;
    }
}

function toggleKBPanel() {
    const kbEnabledEl = document.getElementById('kbEnabled');
    kbEnabled = kbEnabledEl ? kbEnabledEl.checked : false;
}

function openKBModal() {
    const modal = document.getElementById('kbModal');
    if (modal) modal.classList.remove('hidden');
}

function closeKBModal() {
    const modal = document.getElementById('kbModal');
    if (modal) modal.classList.add('hidden');
}

function selectAllKB() {
    document.querySelectorAll('.kb-checkbox').forEach(cb => cb.checked = true);
    updateKBCount();
}

function clearAllKB() {
    document.querySelectorAll('.kb-checkbox').forEach(cb => cb.checked = false);
    updateKBCount();
}

// Get selected KB content for context
async function getSelectedKBContext() {
    const kbEnabledEl = document.getElementById('kbEnabled');
    if (!kbEnabledEl || !kbEnabledEl.checked) return '';

    let context = '';

    // Include KB files
    const checkboxes = document.querySelectorAll('.kb-checkbox:checked');
    if (checkboxes.length > 0) {
        context += '### Knowledge Base Context:\n\n';
        for (const checkbox of checkboxes) {
            try {
                const response = await fetch(`/api/kb/page?name=${encodeURIComponent(checkbox.dataset.name)}`);
                const data = await response.json();
                if (data.content) {
                    context += `--- ${data.title || checkbox.dataset.name} ---\n${data.content}\n\n`;
                }
            } catch (e) {
                console.error('Error loading KB:', e);
            }
        }
    }

    // Include uploaded files
    const selectedUploads = uploadedFiles.filter(f => f.selected);
    if (selectedUploads.length > 0) {
        context += '### Uploaded Files:\n\n';
        for (const file of selectedUploads) {
            context += `--- ${file.name} ---\n${file.content}\n\n`;
        }
    }

    return context;
}

function getKBIcon(filename) {
    if (filename.includes('model')) return '🧠';
    if (filename.includes('architecture')) return '🏗️';
    if (filename.includes('immune') || filename.includes('immunos')) return '🛡️';
    if (filename.includes('reference')) return '📚';
    if (filename.includes('handoff')) return '🤝';
    return '📄';
}

async function viewKBFile(name) {
    try {
        const response = await fetch(`/api/kb/page?name=${encodeURIComponent(name)}`);
        const data = await response.json();

        if (data.content) {
            // Display in chat area as system message
            const preview = data.content.substring(0, 1500);
            addMessage('assistant', `**${data.title || name}**\n\n${preview}${data.content.length > 1500 ? '\n\n...(truncated)' : ''}`, 'KB');
        } else if (data.error) {
            addMessage('error', `KB Error: ${data.error}`);
        }
    } catch (error) {
        addMessage('error', `Failed to load ${name}`);
    }
}

// Recovery & Snapshot Functions
async function runRecovery() {
    const btn = document.getElementById('recoveryBtn');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.textContent = 'Recovering...';

    try {
        const response = await fetch('/api/recovery/run', { method: 'POST' });
        const data = await response.json();

        if (data.success) {
            addMessage('assistant', `Context recovered successfully!\n\n${data.output}`, 'Recovery');

            // Show recovery content summary in chat
            if (data.recovery_content) {
                const lines = data.recovery_content.split('\n').slice(0, 20).join('\n');
                addMessage('assistant', `**Recovery Summary:**\n${lines}\n\n...(truncated)`, 'Context');
            }
        } else {
            addMessage('error', `Recovery failed: ${data.error || 'Unknown error'}`);
        }
    } catch (error) {
        addMessage('error', `Recovery error: ${error.message}`);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function createSnapshot() {
    const btn = document.getElementById('snapshotBtn');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.textContent = 'Creating...';

    try {
        const summary = prompt('Snapshot summary:', 'Dashboard session checkpoint');
        if (!summary) {
            btn.disabled = false;
            btn.textContent = originalText;
            return;
        }

        const response = await fetch('/api/snapshot/create', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ summary: summary, trigger: 'manual' })
        });
        const data = await response.json();

        if (data.success) {
            addMessage('assistant', `Snapshot created!\n${data.output}`, 'Snapshot');
        } else {
            addMessage('error', `Snapshot failed: ${data.error || 'Unknown error'}`);
        }
    } catch (error) {
        addMessage('error', `Snapshot error: ${error.message}`);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

// Check recovery status on load
async function checkRecoveryStatus() {
    try {
        const response = await fetch('/api/recovery/status');
        const data = await response.json();

        if (data.latest_snapshot) {
            const modified = new Date(data.latest_snapshot.modified);
            const hoursAgo = Math.round((Date.now() - modified.getTime()) / (1000 * 60 * 60));

            if (hoursAgo > 4) {
                addMessage('assistant', `Last snapshot was ${hoursAgo} hours ago (${data.latest_snapshot.name}). Consider creating a new one.`, 'System');
            }
        }
    } catch (error) {
        console.error('Error checking recovery status:', error);
    }
}

// Periodic status checks
setInterval(checkOllamaStatus, 30000); // Check every 30 seconds

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
    chatMessages: document.getElementById('chatMessages'),
    chatInput: document.getElementById('chatInput'),
    sendButton: document.getElementById('sendButton')
};

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    initializeWebSocket();
    initializeEventListeners();
    checkOllamaStatus();
    loadModels();
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
    });

    // Ollama controls
    elements.startOllama.addEventListener('click', startOllama);
    elements.stopOllama.addEventListener('click', stopOllama);
}

// Chat Functions
function sendMessage() {
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

    // Send to server
    fetch('/api/chat', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            message: message,
            model: currentModel
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            addMessage('error', data.error);
        }
    })
    .catch(error => {
        addMessage('error', `Failed to send message: ${error.message}`);
    });
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
        const response = await fetch('/api/ollama/status');
        const data = await response.json();
        updateOllamaStatus(data.running);
    } catch (error) {
        updateOllamaStatus(false);
        console.error('Error checking Ollama status:', error);
    }
}

async function startOllama() {
    try {
        elements.startOllama.disabled = true;
        elements.startOllama.textContent = 'Starting...';

        const response = await fetch('/api/ollama/start', { method: 'POST' });
        const data = await response.json();

        if (data.success) {
            updateOllamaStatus(true);
            loadModels();
        } else {
            addMessage('error', `Failed to start Ollama: ${data.error}`);
        }
    } catch (error) {
        addMessage('error', `Error starting Ollama: ${error.message}`);
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

async function loadModels() {
    try {
        const response = await fetch('/api/models/status');
        const data = await response.json();

        elements.modelSelect.innerHTML = '';

        if (data.models && data.models.length > 0) {
            // Add placeholder option
            const placeholderOption = document.createElement('option');
            placeholderOption.value = '';
            placeholderOption.textContent = 'Select a model...';
            elements.modelSelect.appendChild(placeholderOption);

            // Add model options
            data.models.forEach(model => {
                const option = document.createElement('option');
                option.value = model.name;
                option.textContent = `${model.name} (${model.size})`;
                elements.modelSelect.appendChild(option);
            });
        } else {
            const option = document.createElement('option');
            option.value = '';
            option.textContent = 'No models available';
            elements.modelSelect.appendChild(option);
        }
    } catch (error) {
        console.error('Error loading models:', error);
        elements.modelSelect.innerHTML = '<option value="">Error loading models</option>';
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

// Periodic status checks
setInterval(checkOllamaStatus, 30000); // Check every 30 seconds

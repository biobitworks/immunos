---
project: immunos-mcp
source: Offline-Deployment-Architecture.md
type: code-mirror
language: md
size: 20757
modified: 2025-11-26T12:55:03.499227
hash: a8ef4d1a95ee92457deb0dbfe078c201
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `Offline-Deployment-Architecture.md`
> **Size**: 20757 bytes
> **Modified**: 2025-11-26
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```markdown
# Offline Deployment Architecture

> **IMMUNOS-MCP for Air-Gapped and Private Networks**

## Executive Summary

This document describes how to deploy IMMUNOS-MCP in **closed, offline, and air-gapped environments** where internet access is restricted or prohibited. Common scenarios include:

- **Corporate Security Teams** - No external internet for compliance
- **Government/Defense** - Air-gapped networks for classified work
- **Financial Institutions** - Isolated security analysis systems
- **Healthcare** - HIPAA-compliant offline analysis
- **Laptop Deployments** - Security analysts working offline

**Key Requirements Met:**
✅ No internet connectivity required
✅ Runs on standard laptops (8GB+ RAM)
✅ Works on internal wireless networks
✅ Uses local Ollama models (no API keys)
✅ All processing stays on-premises

---

## 🏗️ Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    AIR-GAPPED ENVIRONMENT                        │
│                    (No Internet Access)                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐         ┌──────────────┐                     │
│  │   Analyst    │         │   Analyst    │                     │
│  │   Laptop 1   │◄───────►│   Laptop 2   │                     │
│  │              │  WiFi   │              │                     │
│  │ - Ollama     │         │ - Ollama     │                     │
│  │ - IMMUNOS    │         │ - IMMUNOS    │                     │
│  │ - Scanner    │         │ - Scanner    │                     │
│  └──────┬───────┘         └──────┬───────┘                     │
│         │                        │                              │
│         │      Internal WiFi     │                              │
│         │      (No Internet)     │                              │
│         │                        │                              │
│  ┌──────┴────────────────────────┴───────┐                     │
│  │      Shared Network Storage           │                     │
│  │      (Optional)                        │                     │
│  │  - Model cache                         │                     │
│  │  - Vulnerability databases             │                     │
│  │  - Scan results repository             │                     │
│  └────────────────────────────────────────┘                     │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 💻 Laptop Configuration

### Minimum Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **CPU** | 4 cores | 8+ cores | Apple M1/M2/M3 ideal |
| **RAM** | 8 GB | 16-32 GB | More = larger models |
| **Storage** | 20 GB | 50+ GB | For models + data |
| **OS** | macOS 11+ | macOS 14+ | Linux also supported |
| **GPU** | Optional | Apple Silicon | Significantly faster |

### Recommended Models by Hardware

#### 8GB RAM Laptops (Entry Level)
```bash
# Smallest effective models
ollama pull qwen3:0.6b        # 500MB - Ultra fast
ollama pull phi4:14b           # 8GB - Best quality for size
ollama pull smollm2:1.7b       # 1GB - Very fast embeddings

# Embeddings
ollama pull nomic-embed-text   # 274MB - Code embeddings
```

**Expected Performance:**
- Code scanning: 1-2 seconds per file
- QML learning: 5-10 seconds per model
- Memory usage: 2-4 GB

#### 16GB RAM Laptops (Standard)
```bash
# Balanced models
ollama pull qwen2.5-coder:7b   # 4.7GB - Code analysis
ollama pull deepseek-r1:7b     # 4.7GB - Reasoning
ollama pull phi4:14b           # 8GB - General purpose

# Embeddings
ollama pull nomic-embed-text   # Better code understanding
```

**Expected Performance:**
- Code scanning: <1 second per file
- QML learning: 2-5 seconds per model
- Memory usage: 6-10 GB

#### 32GB+ RAM Laptops (Advanced)
```bash
# Powerful models
ollama pull deepseek-coder-v2:16b  # 9.4GB - GPT-4 level
ollama pull qwq:32b                # 20GB - Advanced reasoning
ollama pull qwen2.5-coder:32b      # 19GB - Best code analysis

# Embeddings
ollama pull bge-m3             # Multi-lingual embeddings
```

**Expected Performance:**
- Code scanning: <0.5 seconds per file
- QML learning: 1-2 seconds per model
- Memory usage: 12-24 GB

---

## 📦 Offline Installation Process

### Step 1: Download on Internet-Connected Machine

```bash
# On a machine WITH internet access

# 1. Download Ollama installer
curl -fsSL https://ollama.com/install.sh -o ollama-install.sh

# 2. Download models
ollama pull qwen3:1.8b
ollama pull nomic-embed-text
ollama pull qwen2.5-coder:7b

# 3. Export models
mkdir ollama-models-export
cp -r ~/.ollama/models/* ollama-models-export/

# 4. Clone IMMUNOS-MCP repository
git clone https://github.com/biobitworks/immunos-mcp.git

# 5. Download Python dependencies
cd immunos-mcp
python3 -m venv .venv
source .venv/bin/activate
pip download -r requirements.txt -d offline-packages/

# 6. Create transfer package
tar -czf immunos-offline-bundle.tar.gz \
    ollama-install.sh \
    ollama-models-export/ \
    immunos-mcp/ \
    offline-packages/
```

### Step 2: Transfer to Air-Gapped System

```bash
# Use approved transfer method:
# - USB drive (approved by security)
# - Internal file server
# - Secure network transfer zone
```

### Step 3: Install on Offline Machine

```bash
# Extract bundle
tar -xzf immunos-offline-bundle.tar.gz

# 1. Install Ollama (if not installed)
bash ollama-install.sh

# 2. Import models
cp -r ollama-models-export/* ~/.ollama/models/

# 3. Verify models
ollama list

# 4. Install IMMUNOS-MCP
cd immunos-mcp
python3 -m venv .venv
source .venv/bin/activate
pip install --no-index --find-links=offline-packages/ -e .

# 5. Test installation
python examples/code_security_scanner/scanner.py
```

---

## 🔧 Optimized Configuration

### Configuration File: `.immunos-offline.toml`

```toml
[deployment]
mode = "offline"
network_type = "internal_only"  # or "standalone"
allow_telemetry = false

[models]
# Primary models (must be pre-installed)
code_analyzer = "qwen2.5-coder:7b"
reasoning_engine = "deepseek-r1:7b"
embedder = "nomic-embed-text"

# Fallback to simple embeddings if Ollama unavailable
use_simple_embeddings_fallback = true

[performance]
# Laptop optimization
max_concurrent_scans = 2
batch_size = 10
enable_gpu = true  # Auto-detect Apple Silicon
cache_embeddings = true
cache_location = "~/.immunos/cache"

[scanner]
# Use lightweight agents for laptops
use_enhanced_nk = false  # Use standard NK Cell (less memory)
bcell_population_size = 30  # Smaller population
nk_detectors = 50  # Fewer detectors

[qml]
# Optimize QML-AiNet for laptops
population_size = 20
max_generations = 30
clone_multiplier = 5

[network]
# For multi-laptop deployments
enable_result_sharing = true
shared_storage_path = "/Volumes/SecurityTeam/immunos-results"
sync_interval_minutes = 15
```

### Load Configuration

```python
# examples/offline_scanner.py
import toml
from pathlib import Path

def load_offline_config():
    """Load offline-optimized configuration."""
    config_path = Path.home() / ".immunos-offline.toml"

    if config_path.exists():
        return toml.load(config_path)

    # Default offline config
    return {
        "deployment": {"mode": "offline"},
        "models": {
            "code_analyzer": "qwen2.5-coder:7b",
            "use_simple_embeddings_fallback": True
        },
        "performance": {
            "max_concurrent_scans": 2,
            "enable_gpu": True,
            "cache_embeddings": True
        }
    }

config = load_offline_config()
```

---

## 🌐 Internal Network Deployment

### Scenario: Security Team with 5 Laptops

```
Security Analyst Team (Air-Gapped Network)
├── Analyst 1 (Lead) - MacBook Pro 32GB
│   ├── Full model suite (deepseek-coder-v2:16b)
│   ├── Central result aggregator
│   └── Training data curator
│
├── Analyst 2-4 (Standard) - MacBook Air 16GB
│   ├── Standard models (qwen2.5-coder:7b)
│   ├── Active scanning
│   └── Report to shared storage
│
└── Analyst 5 (Mobile) - MacBook Air 8GB
    ├── Lightweight models (qwen3:1.8b)
    ├── Quick triage scanning
    └── Sync results when connected
```

### Shared Storage Structure

```
/Volumes/SecurityTeam/immunos-shared/
├── models/
│   ├── qwen2.5-coder-7b/      # Shared model cache
│   ├── nomic-embed-text/
│   └── deepseek-r1-7b/
│
├── datasets/
│   ├── vulnerability-patterns/ # Updated patterns
│   ├── safe-patterns/
│   └── custom-rules/
│
├── results/
│   ├── 2025-11-26/
│   │   ├── analyst1-scans/
│   │   ├── analyst2-scans/
│   │   └── aggregated-report.json
│   └── historical/
│
└── cache/
    └── embeddings/             # Shared embedding cache
```

### Network Sync Script

```python
# examples/network_sync.py
"""
Sync scan results across air-gapped team network.
"""

import shutil
import json
from pathlib import Path
from datetime import datetime

class OfflineTeamSync:
    def __init__(self, shared_path: str, analyst_id: str):
        self.shared_path = Path(shared_path)
        self.analyst_id = analyst_id
        self.local_results = Path.home() / ".immunos" / "results"

    def push_results(self):
        """Push local scan results to shared storage."""
        today = datetime.now().strftime("%Y-%m-%d")
        target = self.shared_path / "results" / today / f"{self.analyst_id}-scans"
        target.mkdir(parents=True, exist_ok=True)

        # Copy new results
        for result_file in self.local_results.glob("*.json"):
            if not (target / result_file.name).exists():
                shutil.copy(result_file, target)
                print(f"✓ Synced: {result_file.name}")

    def pull_patterns(self):
        """Pull updated vulnerability patterns from shared storage."""
        shared_patterns = self.shared_path / "datasets" / "vulnerability-patterns"
        local_patterns = Path("datasets")

        # Update if newer
        for pattern_file in shared_patterns.glob("*.py"):
            shared_mtime = pattern_file.stat().st_mtime
            local_file = local_patterns / pattern_file.name

            if not local_file.exists() or \
               local_file.stat().st_mtime < shared_mtime:
                shutil.copy(pattern_file, local_file)
                print(f"✓ Updated pattern: {pattern_file.name}")

    def get_team_stats(self):
        """Aggregate team scanning statistics."""
        today = datetime.now().strftime("%Y-%m-%d")
        team_results = self.shared_path / "results" / today

        stats = {
            "scans_today": 0,
            "vulnerabilities_found": 0,
            "analysts_active": set()
        }

        for analyst_dir in team_results.glob("*-scans"):
            analyst = analyst_dir.name.split("-")[0]
            stats["analysts_active"].add(analyst)

            for result_file in analyst_dir.glob("*.json"):
                with open(result_file) as f:
                    result = json.load(f)
                    stats["scans_today"] += 1
                    if result.get("is_vulnerable"):
                        stats["vulnerabilities_found"] += 1

        stats["analysts_active"] = len(stats["analysts_active"])
        return stats

# Usage
sync = OfflineTeamSync("/Volumes/SecurityTeam/immunos-shared", "analyst2")
sync.push_results()
sync.pull_patterns()
print(sync.get_team_stats())
```

---

## 🚀 Performance Optimization

### 1. Embedding Cache Strategy

```python
# src/embeddings/cached_embedder.py
"""
Persistent embedding cache for offline deployments.
Reduces computation in air-gapped environments.
"""

import pickle
import hashlib
from pathlib import Path
import numpy as np

class CachedEmbedder:
    def __init__(self, base_embedder, cache_dir="~/.immunos/cache"):
        self.embedder = base_embedder
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.hits = 0
        self.misses = 0

    def embed(self, code: str) -> np.ndarray:
        """Embed with caching."""
        # Hash code for cache key
        code_hash = hashlib.sha256(code.encode()).hexdigest()
        cache_file = self.cache_dir / f"{code_hash}.pkl"

        # Check cache
        if cache_file.exists():
            with open(cache_file, 'rb') as f:
                self.hits += 1
                return pickle.load(f)

        # Generate embedding
        embedding = self.embedder.embed(code)

        # Save to cache
        with open(cache_file, 'wb') as f:
            pickle.dump(embedding, f)

        self.misses += 1
        return embedding

    def get_stats(self):
        """Cache hit rate statistics."""
        total = self.hits + self.misses
        hit_rate = self.hits / total if total > 0 else 0
        return {
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate": hit_rate,
            "cache_size_mb": sum(f.stat().st_size for f in self.cache_dir.glob("*.pkl")) / 1024 / 1024
        }
```

### 2. Batch Processing

```python
# examples/offline_batch_scanner.py
"""
Efficient batch scanning for offline environments.
"""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

class OfflineBatchScanner:
    def __init__(self, scanner, max_workers=None):
        self.scanner = scanner
        # Use fewer workers on laptops
        self.max_workers = max_workers or min(4, multiprocessing.cpu_count())

    def scan_directory(self, directory: Path, extensions=(".py",)):
        """Scan entire directory in parallel."""
        files = []
        for ext in extensions:
            files.extend(directory.rglob(f"*{ext}"))

        print(f"Found {len(files)} files to scan")

        # Parallel scanning
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._scan_file, f): f
                for f in files
            }

            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    print(f"Error scanning {futures[future]}: {e}")

        return results

    def _scan_file(self, file_path):
        """Scan single file."""
        with open(file_path) as f:
            code = f.read()

        result = self.scanner.scan_code(code)
        result.file_path = str(file_path)
        return result
```

---

## 📊 Resource Monitoring

### Laptop Resource Monitor

```python
# examples/resource_monitor.py
"""
Monitor resource usage during offline scanning.
"""

import psutil
import time

class ResourceMonitor:
    def __init__(self):
        self.start_time = time.time()
        self.peak_memory = 0
        self.peak_cpu = 0

    def update(self):
        """Update resource metrics."""
        process = psutil.Process()

        # Memory
        mem = process.memory_info().rss / 1024 / 1024  # MB
        self.peak_memory = max(self.peak_memory, mem)

        # CPU
        cpu = process.cpu_percent()
        self.peak_cpu = max(self.peak_cpu, cpu)

        return {
            "memory_mb": mem,
            "peak_memory_mb": self.peak_memory,
            "cpu_percent": cpu,
            "peak_cpu_percent": self.peak_cpu,
            "runtime_seconds": time.time() - self.start_time
        }

    def get_system_info(self):
        """Get system capabilities."""
        return {
            "total_memory_gb": psutil.virtual_memory().total / 1024**3,
            "available_memory_gb": psutil.virtual_memory().available / 1024**3,
            "cpu_count": psutil.cpu_count(),
            "cpu_freq_mhz": psutil.cpu_freq().current if psutil.cpu_freq() else None
        }

    def recommend_config(self):
        """Recommend configuration based on hardware."""
        info = self.get_system_info()

        if info["total_memory_gb"] < 12:
            return {
                "model": "qwen3:1.8b",
                "population_size": 20,
                "max_workers": 2,
                "note": "Low memory - using lightweight config"
            }
        elif info["total_memory_gb"] < 20:
            return {
                "model": "qwen2.5-coder:7b",
                "population_size": 30,
                "max_workers": 4,
                "note": "Standard config"
            }
        else:
            return {
                "model": "deepseek-coder-v2:16b",
                "population_size": 50,
                "max_workers": 8,
                "note": "High-performance config"
            }
```

---

## 🔒 Security Considerations

### Air-Gap Verification

```python
def verify_air_gap():
    """Verify system is truly offline."""
    import socket

    # Test external connectivity
    try:
        socket.create_connection(("8.8.8.8", 53), timeout=3)
        print("⚠️  WARNING: Internet connectivity detected!")
        print("   This system should be air-gapped.")
        return False
    except (socket.timeout, socket.error):
        print("✓ Air-gap verified - no external connectivity")
        return True

# Run at startup
if __name__ == "__main__":
    verify_air_gap()
```

### Data Retention Policy

```python
# Auto-cleanup old results
from datetime import datetime, timedelta

def cleanup_old_results(days=30):
    """Remove results older than retention period."""
    cutoff = datetime.now() - timedelta(days=days)
    results_dir = Path.home() / ".immunos" / "results"

    removed = 0
    for result_file in results_dir.glob("*.json"):
        if datetime.fromtimestamp(result_file.stat().st_mtime) < cutoff:
            result_file.unlink()
            removed += 1

    print(f"Cleaned up {removed} old results (>{days} days)")
```

---

## 📝 Quick Start Guide

### Individual Laptop Setup

```bash
# 1. Install (one-time)
./offline-install.sh

# 2. Run scanner
python examples/code_security_scanner/demo.py

# 3. Scan a project
python examples/offline_batch_scanner.py /path/to/code
```

### Team Network Setup

```bash
# 1. Mount shared storage
# (varies by IT infrastructure)

# 2. Configure sync
cp .immunos-offline.toml.example ~/.immunos-offline.toml
# Edit shared_storage_path

# 3. Start scanning
python examples/offline_scanner.py

# 4. Sync results (runs automatically or manual)
python examples/network_sync.py
```

---

## 🎯 Best Practices

1. **Model Selection**: Start small (qwen3:1.8b), upgrade as needed
2. **Caching**: Enable embedding cache to speed up repeat scans
3. **Batch Processing**: Scan directories, not individual files
4. **Resource Monitoring**: Watch memory usage, especially on 8GB laptops
5. **Regular Sync**: Push results to shared storage daily
6. **Pattern Updates**: Pull new vulnerability patterns weekly
7. **Testing**: Validate models after offline import
8. **Backups**: Keep model exports backed up on approved media

---

## 📚 Additional Resources

- **Hardware Requirements**: See sections above
- **Model Comparison**: Check Offline-Deployment-Architecture.md
- **Troubleshooting**: See TROUBLESHOOTING.md
- **Network Setup**: Contact IT for shared storage configuration

---

## 🆘 Support for Offline Environments

Since you're offline, refer to:
- Local documentation: `docs/` directory
- Offline help: `python -m immunos.help`
- Team knowledge base: Check shared storage `/docs/`

---

**Document Version:** 1.0
**Last Updated:** 2025-11-26
**Target:** Air-Gapped and Private Network Deployments

```

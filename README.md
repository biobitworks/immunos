# Byron Projects

Master repository for all software engineering projects featuring immune system-inspired AI agents, thermal analysis applications, and aging biology research.

## 🎯 Overview

This workspace consolidates multiple Python and full-stack projects into a single, well-organized repository with unified dependency management and development workflow.

**Obsidian Vault**: This directory is configured as an Obsidian knowledge base with wiki-linking and network graphs.

**Quick Access**: Open [[HOME|HOME.md]] for project dashboards and navigation.

## 📁 Projects

### 🧬 Aging Biology Research (NEW!)
**Longevity Science & Cellular Senescence Research**

Comprehensive knowledge base for aging biology research, integrating databases, tracking researchers/organizations, and building analysis infrastructure.

- **Type**: Research Project + Knowledge Base
- **Purpose**: Systems-level aging biology analysis, data integration, collaboration discovery
- **Scope**: GenAge (307 genes), CellAge (950 genes), protein interactions, spatial transcriptomics
- **Tools**: HuggingFace data hosting, BioGRID networks, spatialLIBD queries

**Dashboard**: [[research/INDEX|Aging Research Dashboard]]

**Key Features**:
- Researcher tracking (Peter Lidsky - pathogen control theory)
- 4 organizations (CityUHK, Monarch Initiative, Future House, longevity-db)
- 3 databases documented (HuggingFace, spatialLIBD, BioGRID)
- Query systems for spatial transcriptomics
- Ready for HuggingFace dataset publication

**Documentation**: [research/README.md](research/README.md)

---

### 🧬 immunos81
**Classical Artificial Immune System for Pattern Recognition**

Implementation of the Immunos-81 algorithm (Hunt & Cooke, 2000) for transparent, interpretable classification.

- **Type**: Python Package
- **Purpose**: Medical diagnosis, pattern recognition, interpretable ML
- **Architecture**: T-cells, B-cells, Clones with affinity/avidity calculations
- **Strategies**: SHA (Simple Highest Avidity), RHA (Relative Highest Avidity)
- **Features**: Online learning, mixed data types, transparent decision-making

**Quick Start**:
```bash
cd immunos81
python examples/simple_demo.py
```

**Documentation**: [immunos81/README.md](immunos81/README.md)

---

### 🤖 IMMUNOS
**Artificial Immune System Orchestrator (Local + Airgapped)**

Multi-agent AIS system where agents play biological immune system roles for adaptive pattern recognition and anomaly detection.

- **Type**: Orchestrator + Dashboard
- **Purpose**: Code security, content moderation, anomaly detection, local AI routing
- **Architecture**: B Cells (pattern matching), NK Cells (anomaly detection), dendritic + memory + T-cell roles
- **Features**:
  - Online/offline routing across multiple model providers
  - NegSl-AIS methodology with adaptive thresholds
  - Thymus training queue + Spleen Summary
  - Agent Foundry (Bone Marrow) for new agent stubs

**Quick Start**:
```bash
python3 scripts/immunos_dashboard.py --port 5001 --host 127.0.0.1
python3 scripts/immunos_chat.py
```

**Documentation**: [IMMUNOS_QUICKSTART.md](IMMUNOS_QUICKSTART.md), [docs/kb/README.md](docs/kb/README.md)

---

### 🌡️ rockbeatspaper
**Full-Stack Thermal Data Analysis Application**

React + TypeScript frontend with Python/Jupyter backend for thermal data analysis using the THRML-HACK package.

- **Type**: Full-Stack Web Application
- **Purpose**: Thermal data analysis and visualization
- **Frontend**: React, TanStack Router, Tailwind CSS
- **Backend**: tRPC, Prisma ORM, Python/Jupyter
- **Infrastructure**: Docker, PostgreSQL, MinIO, Redis
- **Analysis**: JAX-based thermal computing

**Quick Start**:
```bash
cd rockbeatspaper
pnpm install
uv pip install -r requirements.txt
pnpm dev
```

**Documentation**: [rockbeatspaper/README.md](rockbeatspaper/README.md)

---

### 🧪 experiments
**Research & Experimental Code**

Ad-hoc scripts, research explorations, and prototypes.

**Contents**:
- `benford_pi_simple.py`: Benford's Law analysis on pi digits
- `benford_weights_simple.py`: JAX-based weight training

**Finding**: Pi digits do NOT follow Benford's Law (uniform ~11% per digit vs expected log distribution).

---

## 🚀 Quick Start

### Installation

```bash
# Clone or navigate to projects
cd /Users/byron/projects/

# Install workspace dependencies
uv sync

# Or install per-project
cd immunos81 && uv pip install -e .
cd immunos-mcp && uv pip install -e ".[dev,mcp,llm]"
cd rockbeatspaper && uv pip install -r requirements.txt && pnpm install
```

### Running Examples

```bash
# immunos81 - Pattern recognition demo
python immunos81/examples/simple_demo.py

# immunos-mcp - Multi-agent immune system
python immunos-mcp/examples/simple_demo.py
python immunos-mcp/examples/nk_comparison_demo.py

# rockbeatspaper - Start dev server
cd rockbeatspaper && pnpm dev
```

---

## 🛠️ Technology Stack

### Core Python
- **NumPy**: Numerical computing
- **Pandas**: Data manipulation
- **SciPy**: Scientific computing
- **JAX**: Accelerated numerical computing
- **Matplotlib, Plotly**: Visualization

### Machine Learning & AI
- **scikit-learn**: Benchmarking
- **Pydantic**: Data validation
- **Anthropic/OpenAI APIs**: LLM integration (optional)
- **ChromaDB**: Vector database (optional)

### Web & Infrastructure
- **React + TypeScript**: Frontend
- **tRPC**: Type-safe APIs
- **Prisma**: Database ORM
- **Docker**: Containerization
- **PostgreSQL, Redis, MinIO**: Data storage

### Development Tools
- **UV**: Package manager
- **pytest**: Testing
- **black, ruff**: Code formatting
- **Jupyter**: Interactive notebooks

---

## 📚 Research Background

### Immunos-81 (Hunt & Cooke, 2000)
Original artificial immune system for medical diagnosis featuring:
- T-cell and B-cell architecture
- Affinity and avidity calculations
- SHA/RHA recognition strategies
- Transparent, interpretable decisions

### Negative Selection Algorithm (Forrest et al., 1994)
- Self/non-self discrimination
- Anomaly detection without labeled examples
- Inspired by thymic T-cell selection

### Danger Theory (Matzinger, 2002)
- Context-aware immune responses
- Multi-signal integration (PAMP, danger, safe)
- Applied to dendritic cell algorithms

### NegSl-AIS (2025)
Modern negative selection with:
- Adaptive threshold calculation
- Per-class detector generation
- 96.48% arousal, 98.63% valence accuracy on MAHNOB-HCI dataset

---

## 🎓 Key Insights

### Enhancement Lessons (immunos81)

**Attempted**: Multiple clones per class, affinity maturation, lateral inhibition
**Result**: Accuracy decreased from 63.7% → 53.8%
**Finding**: Techniques from optimization (CSA) and anomaly detection (NSA) don't transfer well to supervised classification
**Lesson**: Context matters - apply algorithms to their intended problem domains

### NegSl-AIS Integration Success (immunos-mcp)

**Implemented**: Adaptive thresholds, per-class detectors, enhanced validation
**Method**: min_distance threshold, 15-25 detectors per class
**Result**: Improved accuracy and F1 score
**Validation**: d(detector, nearest_self) > τ (more rigorous than fixed similarity)

---

## 🗺️ Roadmap

### immunos81
- [ ] Cross-validation implementation
- [ ] Feature scaling and normalization
- [ ] Close gap to paper's 83.2% accuracy

### immunos-mcp
- [ ] Complete Dendritic Cell agent
- [ ] Memory agent with vector database
- [ ] T Killer Cell (validation)
- [ ] T Regulatory Cell (calibration)
- [ ] Orchestrator MCP server
- [ ] Real LLM embeddings integration

### rockbeatspaper
- [ ] Complete frontend UI
- [ ] Thermal analysis pipeline integration
- [ ] Production deployment

### General
- [ ] Comprehensive unit tests
- [ ] Performance benchmarks
- [ ] API documentation
- [ ] GitHub repository setup

---

## 📖 Documentation

- **immunos81**: [README.md](immunos81/README.md), [BENCHMARK_ANALYSIS.md](immunos81/BENCHMARK_ANALYSIS.md), [ENHANCEMENTS_SUMMARY.md](immunos81/ENHANCEMENTS_SUMMARY.md)
- **immunos-mcp**: [README.md](immunos-mcp/README.md), [GETTING_STARTED.md](immunos-mcp/GETTING_STARTED.md)
- **rockbeatspaper**: [README.md](rockbeatspaper/README.md)
- **Workspace Context**: [.claude](.claude) - Claude AI tracking and context

---

## 🤝 Development Workflow

### Virtual Environments
- **Primary**: `/Users/byron/projects/.venv/` (shared workspace)
- **Project-specific**: Each project maintains its own if needed
- **Python Version**: >=3.10 (some projects require >=3.14)

### Dependency Management
- **Primary**: UV with workspace configuration
- **Secondary**: pip with requirements.txt
- **Node**: pnpm for TypeScript projects

### Git & Version Control
Each project can be uploaded to GitHub individually:
- `byron/immunos81` - Standalone classifier
- `byron/immunos-mcp` - MCP server
- `byron/rockbeatspaper` - Full-stack app

---

## 📄 License

Individual projects may have their own licenses. Check each project's README for details.

---

## 👤 Author

**Byron** - Software Engineer & AI Researcher

---

## 🔗 Resources

- **Immunos-81 Paper**: Hunt, J. E., & Cooke, D. E. (2000). Learning using an artificial immune system. *Journal of Network and Computer Applications*, 19(2), 189-212.
- **Negative Selection**: Forrest, S., Perelson, A. S., Allen, L., & Cherukuri, R. (1994). Self-nonself discrimination in a computer. *IEEE Symposium on Security and Privacy*.
- **Model Context Protocol**: https://github.com/modelcontextprotocol
- **THRML-HACK**: https://github.com/PrimeIntellect-ai/THRML-HACK

---

**Last Updated**: 2025-12-04
**Workspace Version**: 1.0.0

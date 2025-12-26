# Stage 1 - Inventory (Read-only)

## Intent
Provide a high-level inventory of `/Users/byron/projects` with emphasis on immunOS replication + immunos-mcp productization. This is a top-down scan; deeper file-level review will be scoped in later stages.

## Sources Consulted (initial pass)
- `/Users/byron/projects/README.md`
- `/Users/byron/projects/IMMUNOS_STARTUP_GUIDE.md`
- `/Users/byron/projects/IMMUNOS_QUICKSTART.md`
- `/Users/byron/projects/MULTI_MODEL_SYSTEM_COMPLETE.md`
- `/Users/byron/projects/pyproject.toml`
- `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`
- `/Users/byron/projects/docs/reference/publications-log.md`
- `/Users/byron/projects/data/README.md`
- `/Users/byron/projects/scripts/README.md`
- `/Users/byron/projects/immunos-mcp/README.md`
- `/Users/byron/projects/immunos-preprint/README.md`

## Workspace Top-Level Structure
Directories (selected):
- `immunos-mcp/`: AIS multi-agent system, optional MCP wrapper.
- `immunos81/`: Immunos-81 classifier implementation.
- `immunos-preprint/`: main preprint workspace (includes replication draft).
- `immunos-replication-preprint/`: redirect stub only.
- `data/`: centralized datasets (GenAge, CellAge, etc.).
- `data/immunos_data/`: research datasets (SciFact, TruthfulQA, HaluEval, DiverseVul, etc.).
- `docs/`: reference notes, code snapshots, KB.
- `research/`: literature, datasets, collaborations, organizations.
- `scripts/`: automation + Obsidian tooling.
- `daily/`: day-by-day journals.
- `.immunos/`: agent memory, journals, recovery, and model context files (local-only).

## Key Projects (Relevant to current goals)
- immunOS replication experiment
  - Data: `data/immunos_data/research/scifact/`
  - Tracking: `docs/reference/scifact-baseline-replication.md`
  - Preprint: `immunos-preprint/immunos-preprint-v1.md`
- immunos-mcp
  - AIS multi-agent system with optional MCP packaging
  - Target: public release + VM/USB builder distribution

## Journals, Notes, and Context Stores
- Daily journals: `daily/YYYY-MM-DD.md`
- IMMUNOS journals: `.immunos/journal/`
- Reference notes: `docs/reference/`
- Literature notes: `research/literature/`
- Code snapshots: `docs/code-snapshots/`
- Obsidian vault metadata: `.obsidian/`

## Software Stack (from docs)
Core:
- Python (>=3.10 workspace, with 3.8 legacy env for SciFact)
- uv (workspace/package manager)
- pytest, black, ruff
- Jupyter

ML/AI:
- scikit-learn, numpy, scipy
- torch, transformers, tokenizers (SciFact)
- scispacy + spacy
- ChromaDB (optional, for memory agent)

Web/Infra:
- React + TypeScript
- tRPC, Prisma
- Docker
- PostgreSQL, Redis, MinIO

Ops/Tooling:
- Obsidian (knowledge base)
- conda (legacy SciFact env)

## AI Models Mentioned
LLMs / assistants:
- ChatGPT (Codex CLI)
- Claude Code
- Ollama models: Qwen 2.5 Coder 7B, Qwen 2.5 1.5B, DeepSeek R1 14B

SciFact models:
- `rationale_roberta_large_scifact`
- `label_roberta_large_fever_scifact`

Other model variants referenced in SciFact scripts:
- roberta_large, roberta_base, scibert, biomed_roberta_base

## Publications (Core)
Source of record: `docs/reference/publications-log.md`
- NegSl-AIS (Umair et al., 2025)
- SciFact dataset paper (Wadden et al., 2020)
- CliVER (Ismail et al., 2024)
- SEED (Pasumarthi et al., 2022)
- Negative Selection Algorithms surveys and related AIS literature

## Data Assets (Core)
- SciFact: `data/immunos_data/research/scifact/`
- Hallucination: `data/immunos_data/hallucination/`
- Code vulnerabilities: `data/immunos_data/code/`
- Aging biology datasets (GenAge, CellAge): `data/genage/`, `data/cellage/`

## Known Gaps / Next Questions
- Confirm authoritative list of software/models actually used in production vs optional.
- Define what ships with immunos-mcp vs optional addons.
- Establish publishable artifact list for preprint submission.
- Canonical preprint location is `immunos-preprint/` (replication draft included).

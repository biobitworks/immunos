# Requirements and Assets Register (Draft)

## Software & Tooling
- Python (workspace >=3.10), legacy SciFact env (Python 3.8)
- uv workspace (primary package manager)
- conda (legacy SciFact env)
- Node + pnpm (rockbeatspaper)
- Docker (infra)
- PostgreSQL, Redis, MinIO
- Prisma, tRPC
- Jupyter
- Obsidian (knowledge base)

## AI Models & Frameworks
- LLMs: ChatGPT (Codex CLI), Claude Code, Ollama models (Qwen 2.5 Coder 7B, Qwen 2.5 1.5B, DeepSeek R1 14B)
- SciFact models (transformer baselines):
  - rationale_roberta_large_scifact
  - label_roberta_large_fever_scifact
  - roberta_large, roberta_base, scibert, biomed_roberta_base

## Datasets (Core)
- SciFact: `data/immunos_data/research/scifact/`
- TruthfulQA: `data/immunos_data/hallucination/truthfulqa/`
- HaluEval: `data/immunos_data/hallucination/halueval/`
- DiverseVul: `data/immunos_data/code/diversevul/`
- GenAge + CellAge: `data/genage/`, `data/cellage/`

## Publications (Core)
Source of record: `docs/reference/publications-log.md`
- NegSl-AIS (Umair et al., 2025)
- SciFact dataset (Wadden et al., 2020)
- CliVER (Ismail et al., 2024)
- SEED (Pasumarthi et al., 2022)
- NSA/AIS surveys and supporting AIS literature

## Notes, Journals, and Logs
- Daily journals: `daily/`
- IMMUNOS journals (local): `.immunos/journal/`
- Reference logs: `docs/reference/`
- Literature notes: `research/literature/`
- Code snapshots: `docs/code-snapshots/`
- SciFact replication logs: `docs/reference/scifact-baseline-replication.md`


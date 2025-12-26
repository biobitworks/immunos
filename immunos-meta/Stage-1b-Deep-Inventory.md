# Stage 1b - Deep Inventory (Read-only)

## Summary
This stage provides a deeper directory-level inventory of `/Users/byron/projects` without modifying existing project files. It focuses on structure, size, counts, and key assets.

## Workspace Size Overview (Top-level)
- Large: `data/` (~3.0G), `.venv/` (~852M), `research/` (~565M), `rockbeatspaper/` (~397M), `papers/` (~323M), `docs/` (~264M), `immunos-mcp/` (~212M)
- Medium: `.immunos/` (~31M), `immunos81/` (~33M)
- Small: `immunos-preprint/` (~28K), `immunos-replication-preprint/` (redirect stub), `scripts/` (~1.9M), `daily/` (~124K)

## File/Directory Counts (selected)
- `immunos-mcp/`: 684 dirs, 3631 files
- `rockbeatspaper/`: 577 dirs, 3957 files
- `research/`: 350 dirs, 1899 files
- `docs/`: 263 dirs, 1299 files
- `data/`: 124 dirs, 302 files
- `immunos81/`: 86 dirs, 167 files
- `scripts/`: 10 dirs, 72 files
- `immunos-preprint/`: primary preprint workspace
- `immunos-replication-preprint/`: redirect stub only

## Project-Level Inventory (Core)

### immunOS replication + preprint
- `immunos-preprint/`: main manuscript + replication draft + journal
- `docs/reference/scifact-baseline-replication.md`: canonical metrics log
- `docs/reference/publications-log.md`: non-BibTeX citations
- `data/immunos_data/research/scifact/`: SciFact dataset + scripts

### immunos-mcp
- `immunos-mcp/`: AIS multi-agent system with optional MCP wrapper
- Source layout in `immunos-mcp/src/`
- Docs in `immunos-mcp/docs/`

### immunos81
- `immunos81/`: Immunos-81 AIS classifier

## Data Assets
- Central datasets: `data/` (GenAge, CellAge)
- AIS research datasets: `data/immunos_data/` (SciFact, TruthfulQA, HaluEval, DiverseVul)

## Tooling / Environments
- `pyproject.toml` (uv workspace)
- Requirements found:
  - `data/immunos_data/research/scifact/requirements.txt`
  - `data/immunos_data/hallucination/truthfulqa/requirements.txt`
  - `rockbeatspaper/requirements.txt`
  - proteomics experiments requirements in `research/experiments/`
- Dockerfile found: `rockbeatspaper/docker/Dockerfile`

## Licensing Files (detected)
- `data/immunos_data/hallucination/truthfulqa/LICENSE`
- `data/immunos_data/hallucination/halueval/LICENSE`
- `data/immunos_data/research/scifact/LICENSE.md`

## Bibliography Files
- `immunos-preprint/bibliography.bib`
- `research/citations/immunos-references.bib`
- `research/citations/genage-references.bib`
- `research/citations/cellage-references.bib`
- `papers/references/papers.bib`
- `prion-clock/citations/prion-clock-references.bib`

## Journals & Context Stores
- `daily/`: day-by-day logs
- `.immunos/journal/`: IMMUNOS daily journals (local)
- `docs/reference/`: experiment logs + reference notes
- `research/literature/`: literature archive and notes
- `docs/code-snapshots/`: Obsidian-readable code mirrors

## Observations (PM / Eng)
- Preprint content is now centralized in `immunos-preprint/` (replication draft included).
- Large volume of documentation + code snapshots increases maintenance burden; automation is present but needs governance.
- Data and environment files are spread across project-specific directories; standardization should reduce friction.

## Stage 1b Prompts (for agents)
- "Inventory deep scan completed. Validate canonical preprint location and dataset provenance. Report any missing licenses or dataset terms."
- "Identify minimal reproducibility bundle for SciFact replication (scripts, env, models, data pointers)."

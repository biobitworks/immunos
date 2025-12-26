# Stage 2 - Architecture & Standards (Draft)

## Principles (Industry Best Practice)
- Separation of concerns: core platform vs user projects vs data/compute.
- Reproducibility: lockfile + pinned datasets + versioned models.
- Security: no secrets in repo; explicit boundary for PII and user data.
- Portability: predictable paths and environment bootstrapping.
- Observability: logs, metrics, run manifests, and provenance.

## Proposed Workspace Architecture (User > Project > Environment)

### 1) User Layer
- Global config: `~/.config/immunos/` (or `~/.immunos/`)
- User secrets: `~/.config/immunos/secrets/` (never in repo)
- Global cache: `~/.cache/immunos/` (model weights, dataset cache)

### 2) Project Layer (Logical Project)
- Project root (repo or workspace): `/Users/<user>/projects/<project>`
- Project config: `<project>/immunos.project.yaml`
- Project state: `<project>/.immunos/` (run manifests, logs, indexes)
- Project data: `<project>/data/` (dataset pointers + small artifacts)

### 3) Environment Layer (Runtime)
- One environment per project (default)
  - Python: `uv` + `pyproject.toml` + `uv.lock`
  - Node: `package.json` + lockfile
- Optional per-run containers for reproducibility

## Data Governance & PII Boundaries
- **Repo-safe**: configs, schemas, small datasets, code, docs.
- **Local-only**: PII, secrets, large datasets, model weights.
- **Paths**:
  - `data/` for shared datasets
  - `models/` for model checkpoints (local-only by default)
  - `.immunos/` for memory/logs (local-only by default)

## Directory Documentation Standard (Required)
- Every top-level directory has `README.md` with:
  - `## Summary` (1–3 sentences)
  - `## Directory Map` (tree or abbreviated list)
- Major subfolders (docs/data/src/tests/scripts) also include README with a map.
- Use `INDEX.md` for navigational dashboards (Obsidian-friendly).

## Workspace Organization Plan (Current → Target)
- Current: mixed monorepo + standalone repos with shared datasets in `data/`.
- Target: consistent separation by layer:
  - `core/` (platform libraries)
  - `apps/` (UIs and services)
  - `experiments/` (research, benchmarks)
  - `data/` (datasets, small artifacts)
  - `models/` (weights, local-only by default)
  - `runs/` (metrics, manifests, logs; local-only)
- Keep `.immunos/` as the authoritative run/log store for multi-agent coordination.

## Run Manifests & Provenance
- Each experiment run should emit:
  - `run.json` (command, env hash, dataset hash, model hash)
  - `metrics.json` (summary metrics)
  - `artifacts/` (plots, tables)
- Store in `runs/<project>/<date>/` or `.immunos/runs/`.

## Papers Organization (DOI → Human Search)
- DOI folders are canonical storage.
- Add index files for searchability:
  - `papers/index.csv` and `papers/index.json` (generated from metadata)
  - `papers/INDEX.md` (Obsidian dashboard)
- Next step (design only; not executed):
  - Create `papers/by-year/<year>/` index files or alias folders.
  - Optionally add `papers/by-title/<slug>.md` as pointers to DOI folders.

## Core vs Optional (What Ships)

### Core (ships with immunos-mcp)
- CLI + orchestrator
- Local agents (B/NK/Dendritic/Memory/T-cell stubs)
- Default configs + minimal demo datasets
- Documentation + offline quickstart
- Basic UI (local dashboard)

### Optional Add-ons
- LLM provider integrations (OpenAI, Anthropic, etc.)
- Vector DB integrations (Chroma, Postgres+pgvector)
- Large datasets (SciFact, TruthfulQA, etc.)
- Specialized analyzers (code scanning, biomedical NLP)

## Standardized Repo Layout (Recommended)
```
immunos/
  core/                 # core platform libraries
  cli/                  # CLI entrypoints
  agents/               # AIS agents
  adapters/             # LLM providers + integrations
  ui/                   # optional UI
  configs/              # default config templates
  docs/                 # product docs
  examples/             # minimal demo project
  scripts/              # maintenance
  data/                 # small demo datasets
  .immunos/             # local run state (gitignored)
```

## Environment Bootstrapping
- `immunos init` to create project environment
- `immunos doctor` to validate environment
- `immunos run` for reproducible runs (creates run manifest)

## Standards & Policies
- **Versioning**: semantic versioning for core packages
- **Releases**: tag + changelog + migration notes
- **Docs**: every feature must link to a how-to + a reference
- **Testing**: unit tests for algorithms + integration tests for pipelines
- **Telemetry**: local-only by default; opt-in for analytics

## Open Questions
- Confirm top-level canonical preprint location.
- Decide if immunos-mcp is a standalone repo or monorepo workspace.
- Confirm which datasets can be bundled vs download-on-demand.

# Stage 1c - README Audit (Read-only)

## Goal
Identify which directories are missing README files and which existing READMEs lack a directory map + summary. No changes were made to existing project files.

## Scope
- Scan depth: 2 levels under `/Users/byron/projects` (excluding hidden directories).
- Result: many missing README files; several existing READMEs lack a directory map and/or summary.

## Missing Top-Level READMEs
- `/Users/byron/projects/daily/`
- `/Users/byron/projects/docs/`
- `/Users/byron/projects/projects/`
- `/Users/byron/projects/research/`
- `/Users/byron/projects/static/`
- `/Users/byron/projects/templates/`
- `/Users/byron/projects/todo/`

## Missing READMEs (Selected Subdirectories)
- `/Users/byron/projects/data/immunos_data/`
- `/Users/byron/projects/docs/code-snapshots/`
- `/Users/byron/projects/docs/papers/`
- `/Users/byron/projects/docs/presentations/`
- `/Users/byron/projects/immunos-mcp/{config,docs,examples,logs,scripts,src,tests,web_app}`
- `/Users/byron/projects/immunos81/{core,dashboard,data,engine,examples,jax_accel,strategies,tests,utils}`
- `/Users/byron/projects/immunos-preprint/{data,figures,scripts,supplementary,journal}`
- `/Users/byron/projects/rockbeatspaper/{docker,environments,prisma,public,scripts,src}`
- `/Users/byron/projects/prion-clock/{chat-history,citations,documents,figures,osf-submission,zenodo-submission}`

## Large Collections (Not Practical to Add README per Leaf)
- `/Users/byron/projects/papers/*` includes many DOI-named folders.
  - Recommendation: keep a single `papers/README.md` plus optional index files.
  - Add README only for subfolders that are active projects (`papers/scripts`, `papers/references`, etc.).

## Existing READMEs Missing Directory Maps (Sample)
- `/Users/byron/projects/data/genage/README.md`
- `/Users/byron/projects/data/cellage/README.md`
- `/Users/byron/projects/docs/reference/README.md`
- `/Users/byron/projects/docs/kb/README.md`
- `/Users/byron/projects/immunos-replication-preprint/README.md` (redirect stub)
- `/Users/byron/projects/portfolio/README.md`
- `/Users/byron/projects/research/{citations,collaborations,databases,datasets,literature,organizations,researchers}/README.md`

## Existing READMEs Missing Summary (Sample)
- `/Users/byron/projects/immunos-preprint/README.md`
- `/Users/byron/projects/portfolio/README.md`
- `/Users/byron/projects/data/proteomics-claims-2025/README.md`

## Recommendation (Best Practice)
- Create README files for all **top-level** directories.
- Create README files for **project roots** and **major subfolders** (code, docs, data, models, experiments, scripts).
- Avoid mandatory README for every leaf directory (e.g., each paper DOI folder). Use index files instead.
- Standardize README template: summary + directory map + key links + ownership/notes.

## Decision Needed
Confirm the scope of README creation:
1) Top-level directories only
2) Top-level + project roots (immunos-mcp, immunos81, preprints, rockbeatspaper, prion-clock)
3) Top-level + project roots + major subfolders (docs/data/src/tests/etc.)

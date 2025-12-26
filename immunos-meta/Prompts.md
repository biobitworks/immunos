# Stage Prompts

## Stage 1 Prompt (Inventory)
"Review /Users/byron/projects top-level structure and core docs. Summarize key projects, data assets, journals, AI models, software stack, and publications. Do not edit existing files; write inventory notes in immunos-meta."

## Stage 1b Prompt (Deep Inventory)
"Perform a deeper directory-level inventory of /Users/byron/projects. Capture sizes, file/dir counts, key licenses, requirements, Docker files, and bibliography locations. Do not modify existing project files; write a Stage-1b report in immunos-meta."

## Stage 1c Prompt (README Audit)
"Audit directories for missing README files and missing directory maps/summaries. Produce a Stage-1c report in immunos-meta, and request scope confirmation before creating new READMEs."

## Stage 1d Prompt (Unexpected Changes)
"Inspect unexpected file changes (diffs) and document them in immunos-meta. Do not revert or edit unless explicitly approved."

## Stage 2 Prompt (Architecture & Standards)
"Propose a unified workspace architecture: user > project > environment. Define core vs optional components, data boundaries, and versioning standards. Deliver folder structure and policy guidelines (PII, datasets, logs, models)."

## Stage 3 Prompt (Delivery & Licensing)
"Design a delivery plan for immunOS replication + preprint submission and immunos-mcp productization. Include VM/USB builder website strategy, release artifacts, and licensing options with ROI tradeoffs."

## Stage 4 Prompt (Final Review)
"Produce a final executive summary: gaps, risks, recommended next steps, and staged roadmap."

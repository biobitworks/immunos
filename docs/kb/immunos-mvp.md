---
project: immunos
type: documentation
tags: [immunos, mvp, presentation]
last_updated: 2025-12-22
---

# IMMUNOS MVP

## Problem
Local and airgapped teams need trustworthy, multi-model analysis without cloud dependencies, with clear traceability of self vs non-self detections.

## MVP Scope
- Orchestrator routing across online and offline backends.
- Core AIS agents: dendritic, B cell, NK cell, memory, T cell validator.
- Thymus queue for negative selection training.
- Spleen Summary for global anomaly triage.
- Dashboard chat as the primary operator interface.

## MVP Demo Flow
1) Start dashboard and verify orchestrator routing.
2) Queue a Thymus intake and run training.
3) Run a scan that triggers anomaly events.
4) Show Spleen Summary and open issues list.
5) Switch to offline model and re-run a query.

## Out of Scope (MVP)
- Multi-tenant auth and permissions.
- Full ticketing integration.
- Cloud-only model routing.
- Large-scale distributed training.

## Success Metrics
- Orchestrator routes 90 percent of tasks correctly by domain.
- Thymus jobs complete end-to-end with saved detectors.
- Spleen Summary reflects anomalies within 60 seconds of scan completion.
- Offline fallback works without manual edits.

## Next Milestones
- Agent Foundry creates and stores agent stubs.
- MALT intake connectors for files and API sources.
- Issue sentinel agent auto-tags priorities.

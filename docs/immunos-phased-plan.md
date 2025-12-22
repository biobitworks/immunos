# IMMUNOS Phased Plan (2025-12-22)

## Where we are now
- Chat is the primary dashboard surface (default landing tab).
- Thymus tab exists with intake UI + queue execution for hallucination/network/research.
- Issues (todo) surface in the Overview tab, filtered to `project=immunos`.
- Orchestrator supports online/offline fallback and is configurable from the Chat tab.
- Models + Tokens merged into a single System tab.
- Quick Chat exists on the Overview tab.
- Spleen Summary panel provides global anomaly + issue posture.
- Agent Foundry (Bone Marrow) creates agent stubs from templates.
- KB tab exposes internal docs (feature map, MVP, references).

## Current phase
Phase 2 (Training Expansion) has started. Phase 1 (Dashboard Stabilization) is mostly complete.

## Phased plan

### Phase 0 - Baseline setup (done)
- [x] Monitor dashboard is running on `http://localhost:5001/monitor`
- [x] Orchestrator map + event stream wired
- [x] Issues list from todo system
- [x] Thymus overview chart and agent status cards

### Phase 1 - Dashboard stabilization (mostly done)
- [x] Chat front and center with fallback messaging
- [x] System tab consolidates Models + Tokens
- [x] Quick Chat on Overview
- [x] Add UI hints for missing models / missing Ollama
- [x] Add “health” badge for Thymus queue state (idle/running/failed)
- [x] Add Spleen Summary (global anomaly triage)
- [x] Add Agent Foundry (Bone Marrow) stub UI

### Phase 2 - Training expansion (active)
- [x] Thymus intake queue + auto-run for hallucination/network/research
- [x] Add training pipelines for `code` domain
- [x] Add training pipelines for `emotion` domain
- [ ] Parse training outputs → populate detectors/accuracy in history table
- [x] Add queue controls (pause/resume/run next)

### Phase 3 - Data collection & telemetry
- [ ] Standardize logging schema across orchestrator + domain agents
- [ ] Add daily/weekly summary export to `daily/`
- [ ] Add issue auto-tagging based on detection events

### Phase 4 - Orchestrator routing maturity
- [ ] Add backend profiles (ChatGPT, OpenRouter, local server) in config
- [ ] Add automatic online/offline detection
- [ ] Add cost/latency based routing strategy

### Phase 5 - Packaging & airgapped readiness
- [ ] Vendor all static assets (no CDN dependencies)
- [ ] Offline install script for Ubuntu USB images
- [ ] “Lite” vs “Full” runtime profiles
- [ ] Document airgapped update + data import flow

### Phase 6 - Validation & QA
- [ ] Add smoke tests for `/api/*`
- [ ] Add training validation tests for each domain
- [ ] Regression checks on dashboard rendering

## Immediate next steps (recommended order)
1) Add training result parsing and update Thymus history table.
2) Add MALT intake sources and filters for boundary inputs.
3) Add issue auto-tagging based on detection events.

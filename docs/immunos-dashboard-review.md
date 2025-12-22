# IMMUNOS Dashboard Review and Plan (2025-12-22)

## Scope (files reviewed)
- /Users/byron/projects/scripts/immunos_dashboard.py
- /Users/byron/projects/scripts/immunos_api.py
- /Users/byron/projects/scripts/static/js/monitor.js
- /Users/byron/projects/scripts/templates/monitor.html
- /Users/byron/projects/scripts/immunos_model_manager.py
- /Users/byron/projects/scripts/immunos_token_tracker.py
- /Users/byron/projects/scripts/immunos_routing.py
- /Users/byron/projects/scripts/immunos_todo.py
- /Users/byron/projects/scripts/immunos_data_logger.py

## Quick health check (local)
- CLI curl to `http://localhost:5001` failed with connection refused.
- `lsof` showed a Python process listening on `127.0.0.1:5001`.
- Action: verify from browser or restart the server; confirm the process is the same one running `scripts/immunos_dashboard.py` and that it is bound to `127.0.0.1`.

## Dashboard components that appear wired correctly
- `/monitor` route is defined in `/Users/byron/projects/scripts/immunos_api.py` and renders `monitor.html`.
- `monitor.html` loads `static/js/monitor.js` and includes the Orchestrator Map tab.
- API endpoints referenced by the UI exist:
  - `/api/models/status`
  - `/api/tokens/session/current`
  - `/api/tokens/comparison`
  - `/api/routing/status`
  - `/api/domains/status`
  - `/api/events/simulate`
- WebSocket events expected by the UI exist in `/Users/byron/projects/scripts/immunos_universal.py` when initialized with a SocketIO instance.

## Gaps / risks found
- CDN dependencies (Tailwind, Chart.js, Socket.IO) will not work airgapped. These should be vendored locally for offline installs.
- Training events only emit when `immunos_universal` is executed inside the dashboard server context; there is no training API route yet, so the Training tab can remain empty in normal usage.
- `SECRET_KEY` is hardcoded and CORS is wide open (fine for dev, not for shared machines).
- `/api/domains/status` reads artifacts from `immunos-mcp/.immunos/agents`, which couples the dashboard to that repo layout.
- Token stats are empty unless `immunos_token_tracker.py` is actively logging.

## Installation needs (dashboard and monitor)

### Necessary (minimum)
- Python deps: `flask`, `flask_socketio`, `flask_cors`, `requests`, `psutil`, `numpy`
- App entry point: `/Users/byron/projects/scripts/immunos_dashboard.py`
- Routes: `/Users/byron/projects/scripts/immunos_api.py`
- Templates/assets: `/Users/byron/projects/scripts/templates/monitor.html`, `/Users/byron/projects/scripts/static/js/monitor.js`
- Database: `/Users/byron/projects/.immunos/db/schema.sql` and `/Users/byron/projects/.immunos/db/dashboard.db`

### Optional (can be skipped for quick start)
- Model manager and Ollama status endpoints (`/Users/byron/projects/scripts/immunos_model_manager.py`)
- Token tracking and routing analytics (`/Users/byron/projects/scripts/immunos_token_tracker.py`)
- Docs UI (`/docs`) and Chat API (`/api/chat`)
- Orchestrator map simulation (nice to have, not required to run)

### Good practice (recommended)
- Vendor JS/CSS locally to support airgapped runs.
- Configure secrets and URLs via env vars (`IMMUNOS_SECRET_KEY`, `IMMUNOS_OLLAMA_URL`, `IMMUNOS_BASE_PATH`).
- Add auth or shared-token gating for `/api/*` in multi-user setups.
- Add health checks for `ollama`, DB writeability, and agent artifact availability.
- Provide a pinned `requirements.txt` for `scripts/` services.

## Data collection + issue tracking direction
- Use `/Users/byron/projects/scripts/immunos_data_logger.py` to log interactions into `/Users/byron/projects/.immunos/db/immunos.db`.
- Use `/Users/byron/projects/scripts/immunos_todo.py` as a lightweight ticketing system by tagging issues and treating `todo` entries as tickets.
- Optional: add an `issues` table + `/api/issues/*` endpoints and a dashboard tab for IT triage.

## Phased plan (checklist)

### Phase 0 - Verify and stabilize (today)
- [ ] Confirm the dashboard responds at `http://localhost:5001/monitor`.
- [ ] Verify `/api/health`, `/api/domains/status`, `/api/tokens/session/current`.
- [ ] Confirm `.immunos/db/schema.sql` and `.immunos/db/dashboard.db` exist and are writable.
- [ ] Open devtools console in the browser for JS errors on `monitor.html`.

### Phase 1 - Data collection baseline (fast)
- [ ] Wire `immunos_data_logger.py` into model calls (routing + chat + detect).
- [ ] Define a shared schema for `task_type`, `task_classification`, `routing_reason`, `confidence_score`.
- [ ] Add a daily export summary to `/Users/byron/projects/daily/`.

### Phase 2 - Issue and ticketing MVP
- [ ] Decide: extend `immunos_todo.py` or add a new `issues` table.
- [ ] Add CLI wrapper for issues (create, update, resolve).
- [ ] Add `/api/issues` and a dashboard tab with filters (priority, owner, status).
- [ ] Log issue events into the activity feed for audit.

### Phase 3 - Orchestrator and training wiring
- [ ] Ensure orchestrator emits to `/api/events/detection`.
- [ ] Add API endpoints to start training jobs and emit `training_*` events.
- [ ] Add per-domain health score computed from logged metrics.

### Phase 4 - Packaging and airgapped readiness
- [ ] Vendor all static assets locally (no CDN calls).
- [ ] Provide an offline install script for Ubuntu USB images.
- [ ] Add a minimal "lite" profile vs "full" profile configuration.
- [ ] Document startup and validation steps in `/Users/byron/projects/docs/`.

## Update log
### 2025-12-22
- Renamed Training tab to Thymus and added a training overview chart plus per-domain cards.
- Added IMMUNOS issues summary/list to the overview page (from todo system).
- Added Spleen Summary panel for global anomaly + issue posture.
- Added Agent Foundry (Bone Marrow) UI for template-based agent stubs.
- Added KB tab with in-dashboard documentation viewer.
- Added orchestrator chat tab with online/offline backend toggles and saved config.
- Promoted Chat to the first tab and default landing view; removed empty Rigor tab.
- Added Thymus intake queue UI + API to capture dataset paths/notes.
- Merged Models + Tokens into a single System tab.
- Added Quick Chat on Overview for fast orchestrator prompts.
- Thymus intake now triggers training queue execution with progress events.

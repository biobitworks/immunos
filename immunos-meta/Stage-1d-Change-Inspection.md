# Stage 1d - Unexpected Changes Review (Inspect, Document, Ignore)

## Intent
Inspect unexpected modifications and document them. No changes were reverted.

## Files with unexpected diffs (root repo)

### `scripts/immunos_api.py`
- Adds ticketing + sentinel API endpoints:
  - GET /api/tickets, GET /api/tickets/<id>, PATCH /api/tickets/<id>
  - GET /api/tickets/stats
  - POST /api/sentinel/start
  - GET /api/sentinel/patterns
  - POST /api/tokens/estimate
- Bumps API route summary for token tracking from 5 to 6 endpoints.

### `scripts/immunos_dashboard.py`
- Adds `/tickets` route rendering `tickets.html`.

### `scripts/templates/base.html`
- Adds navigation item for Tickets page.
- Updates cache-busting query params (`dashboard.css`, `app.js`).
- Version label changed to `v20251226-sentinel`.

## Submodule/Repo Status (README additions only)
- `immunos-mcp`: README updated + new subfolder READMEs
- `immunos81`: README updated + new subfolder READMEs
- `prion-clock`: README updated + new subfolder READMEs
- `rockbeatspaper`: README updated + new subfolder READMEs
- `experiments`: README updated

## Data/Repo Status
- `data/immunos_data/research/scifact` appears as an untracked submodule in root status; contains its own `.git` directory.

## Action
- Documented only. Left changes untouched per instruction.

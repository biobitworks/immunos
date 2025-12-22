---
project: immunos
type: documentation
tags: [immunos, handoff, prompts, claude-code]
last_updated: 2025-12-22
---

# Handoff and Model Switching

This page captures safe handoff workflows between Claude Code, ChatGPT, and local models.
Do not include API keys or proprietary data in handoff notes.

## When to Handoff

- Token usage warning in the dashboard (150k+).
- Offline/airgapped workflows.
- Switching between providers (Claude Code -> ChatGPT -> Ollama).

## Manual Handoff Command

```bash
python3 scripts/immunos_handoff.py save --task "Current work" --steps "Next step 1" "Next step 2"
```

## Prompt Template (Claude Code / ChatGPT)

Paste this into the new model when handing off:

```
PROJECT: IMMUNOS
GOAL: Continue work on IMMUNOS dashboard + AIS orchestration.

CONTEXT:
- Dashboard runs at http://localhost:5001/monitor
- Orchestrator config stored in .immunos/config/orchestrator.json
- Training queue under .immunos/thymus/intake.json

FILES TO REVIEW:
- scripts/immunos_api.py
- scripts/static/js/monitor.js
- scripts/templates/monitor.html
- docs/immunos-phased-plan.md
- docs/immunos-dashboard-review.md
- docs/kb/README.md

WHAT CHANGED RECENTLY:
- Spleen Summary panel + /api/spleen/summary
- Agent Foundry (Bone Marrow) templates + /api/agents/foundry
- KB pages: feature map, MVP, references

NEXT STEPS:
1) Verify monitor updates render after restart.
2) Extend training history parsing in Thymus tab.
3) Add MALT intake sources and filters.
```

## Offline Fallback Checklist

- Start Ollama and confirm models list in System tab.
- Ensure offline backend is selected in the chat tab.
- Verify chat responses and update token counters.

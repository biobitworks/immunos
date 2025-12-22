---
project: immunos
type: presentation
tags: [immunos, mvp, slides]
last_updated: 2025-12-22
---

# IMMUNOS MVP (Slides Draft)

## Slide 1: MVP Goal
Run a multi-agent AIS locally with offline fallback and visible self vs non-self signals.

## Slide 2: MVP Features
- Orchestrator routing (online/offline)
- Core AIS agents (B cell, NK cell, dendritic, memory, T cell)
- Thymus training queue
- Spleen Summary (global anomaly view)
- Dashboard-first operations

## Slide 3: Demo Flow
1) Start dashboard
2) Queue Thymus intake
3) Trigger scan
4) Review Spleen Summary
5) Switch offline model

## Slide 4: MVP Metrics
- Correct routing by domain
- Training jobs complete with saved detectors
- Anomaly summary updates within 60 seconds

## Slide 5: Next After MVP
- Agent Foundry expansion
- MALT intake connectors
- Issue sentinel auto-tagging

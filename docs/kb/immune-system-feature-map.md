---
project: immunos
type: documentation
tags: [immunos, feature-map, immune-system]
last_updated: 2025-12-22
---

# Immune System Feature Map

This maps immune system structures to IMMUNOS components, UI surfaces, and code targets.

Legend:
- existing: already in repo or dashboard
- planned: recommended next step
- optional: later or nice-to-have

## Primary Lymphoid Organs

- Bone Marrow (agent genesis) -> Agent Foundry for new agent templates and detector seeding.
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/config/`, `immunos-mcp/src/immunos_mcp/training/`
  - UI: Foundry panel for templates and stub creation
- Thymus (negative selection) -> Training and selection pipeline and queue.
  - Status: existing
  - Code targets: `immunos-mcp/src/immunos_mcp/training/`
  - UI: Thymus tab (intake + queue + charts)

## Secondary Lymphoid Organs

- Spleen (global filtering) -> System-wide anomaly triage and resolution state.
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/orchestrator/`
  - UI: Spleen Summary panel for global risk and anomalies
- Lymph Nodes (local hubs) -> Domain stations (code, network, research, hallucination).
  - Status: partial
  - Code targets: `immunos-mcp/src/immunos_mcp/training/`
  - UI: per-domain performance, detector coverage, alerts
- MALT (entry protections) -> Intake connectors and boundary sensors.
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/core/antigen.py`
  - UI: MALT intake sources and filters

## Endocrine + Metabolic Coupling

- Adipose Tissue (context reservoir) -> long-lived embeddings and sensitivity tuning.
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/embeddings/`, `immunos-mcp/src/immunos_mcp/agents/memory_agent.py`
  - UI: storage health and retention policies
- HPA/HPT axes (resource control) -> budget and stress response (tokens, CPU/GPU).
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/orchestrator/`
  - UI: resource governor and auto-routing

## Barrier + Support Structures

- Skin (perimeter barrier) -> input validation, sanitization, PII screening.
  - Status: planned
  - Code targets: `immunos-mcp/src/immunos_mcp/core/antigen.py`
  - UI: input policy and redaction rules
- Liver (acute response) -> emergency playbooks and artifact cleanup.
  - Status: optional
  - Code targets: `immunos-mcp/src/immunos_mcp/orchestrator/`
  - UI: incident log and response toggles
- Lymphatic vessels (transport) -> event bus and queue health.
  - Status: partial
  - Code targets: `immunos-mcp/src/immunos_mcp/training/`
  - UI: event timeline and queue status

## Low-Friction Next Features

1) Agent Foundry stub UI + config templates.
2) Spleen Summary panel with global anomaly counts.
3) MALT intake sources list and filters.

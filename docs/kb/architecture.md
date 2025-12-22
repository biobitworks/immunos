# Architecture

## Overview
- IMMUNOS is a multi-agent artificial immune system with an orchestrator coordinating domain agents.
- Primary domains: hallucination, network, research, code, emotion.

## Core components
- Orchestrator (routing + fallback)
- Agents: B Cell, NK Cell, Memory, T Cell, Dendritic
- Thymus (training intake + execution)
- Monitor dashboard

## Key paths
- Orchestrator config: `/Users/byron/projects/.immunos/config/orchestrator.json`
- Agent artifacts: `/Users/byron/projects/immunos-mcp/.immunos/agents`
- Thymus queue: `/Users/byron/projects/.immunos/thymus/intake.json`

## Update notes
- 2025-12-22: Initial KB structure created.

---
date: 2025-12-01 (planned)
tags: [planning, next-session, mcp, datasets]
---

# Next Session Plan

## Goals
1. **Collect benchmark/validation datasets** (1000+ code samples)
2. **Implement MCP server** (VS Code immune agents integration)
3. **Expand training data** (improve agent accuracy)

## Session Structure (2-3 hours)

### Phase 1: Data Collection (45 min)
**Datasets to gather:**
- SARD (Software Assurance Reference Dataset)
- Juliet Test Suite (118+ CWE categories)
- OWASP Benchmark
- Safe code: Python stdlib, Django, Flask, NumPy
- Vulnerable code: CVE PoCs, OWASP apps

**Actions:**
1. Create `datasets/` directory structure
2. Download and extract samples
3. Write conversion scripts: External → Antigen format
4. Label: safe (500+), vulnerable (500+)

### Phase 2: MCP Server Implementation (60 min)
**Files to create:**
- `immunos-mcp/src/servers/immunos_mcp_server.py` (~450 lines)
- `~/.continue/mcpServers/immunos.yaml` (~15 lines)
- `immunos-mcp/src/servers/pattern_manager.py` (~200 lines)

**Tools to expose:**
- `/immune-bcell` → B Cell pattern matching
- `/immune-nk` → NK Cell anomaly detection
- `/immune-full` → Full 6-agent analysis

**Integration:**
```
VS Code → Continue.dev → MCP Server → 6 Immune Agents → Ollama
```

### Phase 3: Training Expansion (30 min)
1. Train B Cell on 1000+ labeled samples
2. Train NK Cell on 500+ safe samples (self)
3. Save patterns with `pattern_manager.py`
4. Update MCP server to load saved patterns
5. Validate on held-out test set

### Phase 4: Documentation (15 min)
1. Update `daily/2025-12-01.md`
2. Create experiment log for validation
3. Update `NEXT_SESSION.md` → archive
4. Git commit

## Expected Deliverables

- ✅ 1000+ labeled code samples
- ✅ MCP server working in VS Code
- ✅ 3 tools callable via `/immune-*` commands
- ✅ Improved agent accuracy (target: >85%)
- ✅ Benchmark validation results
- ✅ Documentation in Obsidian
- ✅ Git commit

## Success Criteria

1. Can run `/immune-full` in VS Code Continue on selected code
2. Get comprehensive analysis with risk levels + explanations
3. Agents trained on 1000+ samples show >85% accuracy
4. Patterns persist across server restarts
5. Complete documentation for reproducibility

## Resources

**Current State:**
- Continue.dev installed and configured
- 3 Ollama models running locally
- Obsidian vault set up (70+ files)
- 6 LLM-enhanced agents implemented

**Documentation:**
- [[docs/VS-Code-LLM-Integration]] - Continue setup
- [[research/experiments/qml-ainet-validation-2025-11-30]] - Current accuracy
- [[projects/immunos-mcp/diagrams/system-overview]] - Architecture

**Research from today:**
- MCP integration architecture designed
- 3 tools planned with detailed specs
- Pattern persistence strategy defined

## Preparation Before Session

Optional prep work:
1. Download SARD dataset: https://samate.nist.gov/SARD/
2. Clone Juliet: https://samate.nist.gov/SARD/test-suites/112
3. Review MCP docs: https://modelcontextprotocol.io/
4. Install MCP SDK: `uv pip install mcp>=1.3.1`

## Notes

- Session built on today's Obsidian vault setup
- Continue.dev already tested and working
- MCP architecture researched (see research agent output)
- Estimated time: 2-3 hours for full implementation

---

**Created**: 2025-11-30
**Status**: Ready to execute
**Priority**: High (enables real immune analysis in VS Code)

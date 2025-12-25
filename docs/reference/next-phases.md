# Next Phases Plan

**Created**: 2025-12-24
**Status**: Active

---

## Phase 1: Environment Finalization (Immediate)

| Task | Owner | Status |
|------|-------|--------|
| Run `conda init zsh` to enable conda commands | User | Pending |
| Verify TorchEEG can load MAHNOB-HCI | Agent | Pending |
| Test proteomics env with existing notebooks | Agent | Pending |

**Command**:
```bash
/opt/anaconda3/bin/conda init zsh
# Then restart terminal
conda activate immunos-ml
python -c "import torcheeg; print('TorchEEG:', torcheeg.__version__)"
```

---

## Phase 2: MAHNOB-HCI Dataset Access

| Task | Owner | Status |
|------|-------|--------|
| Draft email to iBUG team for dataset access | Agent | Pending |
| Check Internet Archive for MAHNOB portal form | Agent | Pending |
| Prepare fallback to DEAP/DREAMER if blocked | Agent | Pending |

**Contact**: iBUG group at Imperial College London
**Portal**: https://mahnob-db.eu/ (HCI-Tagging subpage 404)

---

## Phase 3: Preprint Development

### immunOS Replication Preprint (existing)
| Task | Status |
|------|--------|
| SciFact baseline metrics complete | ✓ |
| Citation Relevance table added | ✓ |
| Open dev F1 tables in place | ✓ |
| Add abstract + author info | Pending |
| Zenodo DOI reservation | Pending |

### immunOS Main Preprint (new scaffold)
| Task | Status |
|------|--------|
| Project scaffold created | ✓ |
| Define scope vs replication preprint | Pending |
| Outline main contributions | Pending |

---

## Phase 4: Prion Clock Submission

| Task | Status |
|------|--------|
| All Opus 4.5 critical fixes complete | ✓ |
| 30/30 citations verified | ✓ |
| Add author attribution | Pending |
| Final proofread | Pending |
| Submit to bioRxiv | Pending |

**Location**: `/Users/byron/projects/prion-clock/`

---

## Phase 5: IMMUNOS Blog Security Fixes

5 critical issues before publication:

| Issue | Status |
|-------|--------|
| Hardcoded secrets | Pending |
| CORS configuration | Pending |
| Missing functions | Pending |
| Input validation | Pending |
| Error handling | Pending |

**Location**: `/Users/byron/projects/scripts/immunos_dashboard.py`

---

## Phase 6: Orchestrator Enhancement

| Task | Status |
|------|--------|
| Terminal smoke tests logged | ✓ |
| Model roles mapped | ✓ |
| Add LLM-enhanced agents to main orchestrator | Pending |
| Create interactive REPL mode | Pending |

---

## Agent Assignments

| Phase | Primary Agent | Tools |
|-------|---------------|-------|
| 1 | Claude Code | Bash, verify envs |
| 2 | General-purpose | WebFetch, email draft |
| 3 | Claude Code | Edit preprints |
| 4 | Claude Code | Final edits, bioRxiv prep |
| 5 | General-purpose | Security audit, code fixes |
| 6 | General-purpose | Orchestrator development |

---

## Recommended Priority Order

1. **Phase 1** - Quick win, unblocks ML work
2. **Phase 4** - Prion Clock ready, low effort to submit
3. **Phase 3** - Preprints need abstracts for Zenodo
4. **Phase 2** - Dataset access may take time (external dependency)
5. **Phase 5** - Security fixes before blog goes public
6. **Phase 6** - Enhancement, not blocking

---

## Quick Commands

```bash
# Morning startup
python3 scripts/immunos_recover.py
cat .immunos/recovery/CONTEXT_RECOVERY.md

# Check what's next
cat docs/reference/next-phases.md

# Process any new inbox files
python3 scripts/immunos_inbox.py --verbose

# Activate ML environment
conda activate immunos-ml
```

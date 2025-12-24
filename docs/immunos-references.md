---
project: immunos
type: references
last_updated: 2025-12-22
---

# IMMUNOS References and Reproducibility Sources

This file lists **only sources available in `/Users/byron/projects` or explicitly provided by the user**.
No external claims are added without a direct link.

## Core AIS / IMMUNOS Papers (DOI links)

- NegSl-AIS baseline paper (Results in Engineering, 2025): https://doi.org/10.1016/j.rineng.2025.106601
- IMMUNOS-81 original paper: https://doi.org/10.1016/B978-0-12-354840-7.50037-6
- Negative Selection Algorithms review: https://doi.org/10.1162/106365601750190975
- V-Detector disease diagnosis (user-provided abstract): https://doi.org/10.1016/B978-0-12-803468-2.00011-4
- Applied Soft Computing AIS paper (user-provided): https://doi.org/10.1016/j.asoc.2021.108129
- Immunological computation (NCBI Bookshelf, Chapter 23): https://www.ncbi.nlm.nih.gov/books/NBK459484/

## SciFact Baselines (data + code)

- Paper (arXiv): https://arxiv.org/abs/2004.14974
- Dataset + code (local): `/Users/byron/projects/data/immunos_data/research/scifact`
  - Dataset tarball: https://scifact.s3-us-west-2.amazonaws.com/release/latest/data.tar.gz
  - Baseline scripts:
    - `script/rationale-selection.sh`
    - `script/label-prediction.sh`
    - `script/pipeline.sh`
- Leaderboard: https://leaderboard.allenai.org/scifact/
- SOTA reference mentioned in dataset README: https://github.com/dwadden/multivers

## Reported Baseline Metrics (quoted from provided abstract)

The following metrics are **quoted** from the user-provided abstract and are not independently verified here:

- V-Detector disease diagnosis results:
  - Breast Cancer Wisconsin: 98.95% detection
  - BUPA Liver Disorder: 74.44% detection
  - Biomedical data (unspecified): 71.64% detection
- NegSl-AIS multimodal emotion classification (Results in Engineering, 2025):
  - Dataset: MAHNOB-HCI (data availability: on request per paper excerpt)
  - Arousal accuracy: 96.48%
  - Valence accuracy: 98.63%
  - Overall accuracy: 94%
  - Cohen's Kappa: 0.919
  - MCC: 0.920

## NegSl-AIS Baseline (local copy)

- Local paper: `/Users/byron/projects/papers/umair2025negsl/1-s2.0-S2590123025026702-main.pdf`
- Local summary: `/Users/byron/projects/papers/umair2025negsl/README.md`
- Analysis notes: `/Users/byron/projects/.immunos/journal/IMMUNOS_FOUNDATION_PAPER.md`
- Baseline summary: `/Users/byron/projects/docs/reference/negsl-ais-2025-baseline.md`

## IMMUNOS Evaluation Harness

- Script: `/Users/byron/projects/scripts/immunos_training_eval.py`
- Latest publication-grade k-fold report (gold evidence + mxbai, balanced):
  - `/Users/byron/projects/docs/training-eval-2025-12-22-research-ollama-mxbai-kfold-final.md`

## Baseline Replication TODO

- TODO-20251223-001: Replicate baselines from AIS/SciFact papers (download data, run scripts, record metrics, compare).
- SciFact replication log: `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`

## Gaps / To Verify

- DOI links and code/data repos for Claude-provided references in `/Users/byron/projects/` (not yet linked here).
- Any DOI for the SciFact arXiv paper (if exists, add direct DOI link).

## Zotero / BibTeX

- NegSl-AIS reference list (DOI-only, for Zotero import):
  - `/Users/byron/projects/research/citations/immunos-references.bib`
- Publications log (non-BibTeX):
  - `/Users/byron/projects/docs/reference/publications-log.md`

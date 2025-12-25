# SciFact Baseline Replication Log

## Sources
- Repo: `/Users/byron/projects/data/immunos_data/research/scifact`
- Dataset root: `/Users/byron/projects/data/immunos_data/research/scifact/data`
- Metrics (oracle dev): `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`
- Metrics (open dev): `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev_open.json`
- Eval log (oracle dev): `/Users/byron/projects/data/immunos_data/research/scifact/logs/eval_dev.log`

## Environment (x86_64)
- Miniconda: `/Users/byron/miniconda3-x86`
- Env: `scifact-py38-legacy` (Python 3.8, osx-64)
- Pinned requirements installed (legacy resolver): numpy 1.18.2, scipy 1.4.1, scikit-learn 0.22.2.post1, torch 1.5.0, transformers 2.7.0, tokenizers 0.5.2, requests 2.23.0, urllib3 1.26.5.
- SciSpaCy stack: scispacy 0.2.5, spacy 2.3.9, en_core_sci_sm 0.2.5.
- sentencepiece 0.1.85 built locally; dylibs copied into env `lib/`.

## Oracle Dev Results

Metrics from `logs/metrics_dev.json`:

| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.794 | 0.590 | 0.677 |
| Sentence label | 0.713 | 0.530 | 0.608 |
| Abstract label (only) | 0.910 | 0.675 | 0.775 |
| Abstract label (rationalized) | 0.852 | 0.632 | 0.725 |

Notes:
- Evaluator warns when rationales exceed 3 sentences; it truncates to 3 for abstract-level scoring.

## Open Dev Results

Metrics from `logs/metrics_dev_open.json`:

| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.525 | 0.437 | 0.477 |
| Sentence label | 0.469 | 0.391 | 0.426 |
| Abstract label (only) | 0.553 | 0.474 | 0.510 |
| Abstract label (rationalized) | 0.525 | 0.450 | 0.485 |

## Open Retrieval (dev) Complete

Step A (TF-IDF retrieval) completed:

- Output: `/Users/byron/projects/data/immunos_data/research/scifact/prediction/abstract_retrieval_open.jsonl`
- Log: `/Users/byron/projects/data/immunos_data/research/scifact/logs/open_retrieval.log`

Step B (rationale selection) completed:

- Output: `/Users/byron/projects/data/immunos_data/research/scifact/prediction/rationale_selection_open.jsonl`
- Log: `/Users/byron/projects/data/immunos_data/research/scifact/logs/open_rationale.log`

Step C (label prediction + merge + metrics) completed:

- Output: `/Users/byron/projects/data/immunos_data/research/scifact/prediction/label_prediction_open.jsonl`
- Output: `/Users/byron/projects/data/immunos_data/research/scifact/prediction/merged_predictions_open.jsonl`
- Logs: `/Users/byron/projects/data/immunos_data/research/scifact/logs/open_label.log`, `/Users/byron/projects/data/immunos_data/research/scifact/logs/open_autorun.log`

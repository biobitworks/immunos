# IMMUNOS Training Evaluation (2025-12-22 23:17)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=False)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train labels: {'unknown': 304, 'contradict': 173, 'support': 332}
  - Test labels: {'unknown': 112, 'support': 124, 'contradict': 64}
  - Train/Test: 809 / 300
  - B Cell accuracy: 41.3% (macro-F1 19.5%)
  - NK accuracy: 53.7% (precision 57.7%, recall 78.4%)

- Research (SciFact K-Fold)
  - Split: kfold-5 (balanced=False)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=3)
  - NK mode: classic
  - B Cell accuracy: 41.2% ± 0.5%
  - B Cell macro-F1: 20.1% ± 0.7%
  - NK accuracy: 53.5% ± 1.8%
  - NK F1: 66.3% ± 1.2%

## Detailed Metrics

### Research_Dev
- bcell_metrics: `{
  "accuracy": 0.41333333333333333,
  "macro_precision": 0.13777777777777778,
  "macro_recall": 0.3333333333333333,
  "macro_f1": 0.19496855345911948,
  "total": 300,
  "correct": 124,
  "per_class_accuracy": {
    "contradict": 0.0,
    "support": 1.0,
    "unknown": 0.0
  },
  "per_class_precision": {
    "contradict": 0.0,
    "support": 0.41333333333333333,
    "unknown": 0.0
  },
  "per_class_recall": {
    "contradict": 0.0,
    "support": 1.0,
    "unknown": 0.0
  },
  "per_class_f1": {
    "contradict": 0.0,
    "support": 0.5849056603773585,
    "unknown": 0.0
  }
}`
- nk_metrics: `{
  "accuracy": 0.5366666666666666,
  "precision": 0.5774058577405857,
  "recall": 0.7840909090909091,
  "f1": 0.6650602409638554,
  "tp": 138,
  "fp": 101,
  "tn": 23,
  "fn": 38
}`

### Research_Kfold

# IMMUNOS Training Evaluation (2025-12-22 23:04)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train labels: {'unknown': 173, 'support': 173, 'contradict': 173}
  - Test labels: {'unknown': 112, 'support': 124, 'contradict': 64}
  - Train/Test: 519 / 300
  - B Cell accuracy: 70.7% (macro-F1 64.8%)
  - NK accuracy: 53.7% (precision 57.1%, recall 84.1%)

- Research (SciFact Test)
  - Split: train/test (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train labels: {'unknown': 173, 'support': 173, 'contradict': 173}
  - Test labels: {'unknown': 300}
  - Warning: test split has a single label; test metrics are not publication-grade.
  - Train/Test: 519 / 300
  - B Cell accuracy: 94.3% (macro-F1 32.4%)
  - NK accuracy: 100.0% (precision 100.0%, recall 100.0%)

## Detailed Metrics

### Research_Dev
- bcell_metrics: `{
  "accuracy": 0.7066666666666667,
  "macro_precision": 0.6543290043290043,
  "macro_recall": 0.645497311827957,
  "macro_f1": 0.6481684560432189,
  "total": 300,
  "correct": 212,
  "per_class_accuracy": {
    "contradict": 0.28125,
    "support": 0.717741935483871,
    "unknown": 0.9375
  },
  "per_class_precision": {
    "contradict": 0.32727272727272727,
    "support": 0.6357142857142857,
    "unknown": 1.0
  },
  "per_class_recall": {
    "contradict": 0.28125,
    "support": 0.717741935483871,
    "unknown": 0.9375
  },
  "per_class_f1": {
    "contradict": 0.3025210084033613,
    "support": 0.6742424242424243,
    "unknown": 0.967741935483871
  }
}`
- nk_metrics: `{
  "accuracy": 0.5366666666666666,
  "precision": 0.5714285714285714,
  "recall": 0.8409090909090909,
  "f1": 0.6804597701149425,
  "tp": 148,
  "fp": 111,
  "tn": 13,
  "fn": 28
}`

### Research_Test
- bcell_metrics: `{
  "accuracy": 0.9433333333333334,
  "macro_precision": 0.3333333333333333,
  "macro_recall": 0.31444444444444447,
  "macro_f1": 0.3236134934248142,
  "total": 300,
  "correct": 283,
  "per_class_accuracy": {
    "contradict": 0.0,
    "support": 0.0,
    "unknown": 0.9433333333333334
  },
  "per_class_precision": {
    "contradict": 0.0,
    "support": 0.0,
    "unknown": 1.0
  },
  "per_class_recall": {
    "contradict": 0.0,
    "support": 0.0,
    "unknown": 0.9433333333333334
  },
  "per_class_f1": {
    "contradict": 0.0,
    "support": 0.0,
    "unknown": 0.9708404802744426
  }
}`
- nk_metrics: `{
  "accuracy": 1.0,
  "precision": 1.0,
  "recall": 1.0,
  "f1": 1.0,
  "tp": 300,
  "fp": 0,
  "tn": 0,
  "fn": 0
}`

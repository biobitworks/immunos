# IMMUNOS Training Evaluation (2025-12-23 00:05)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (bge-m3)
  - Evidence: gold (max=3)
  - NK mode: classic
  - Train labels: {'unknown': 173, 'support': 173, 'contradict': 173}
  - Test labels: {'unknown': 112, 'support': 124, 'contradict': 64}
  - Train/Test: 519 / 300
  - B Cell accuracy: 70.0% (macro-F1 65.9%)
  - NK accuracy: 54.0% (precision 57.3%, recall 84.7%)

- Research (SciFact K-Fold)
  - Split: kfold-3 (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (bge-m3)
  - Evidence: gold (max=3)
  - NK mode: classic
  - B Cell accuracy: 72.5% ± 2.7%
  - B Cell macro-F1: 69.3% ± 2.8%
  - NK accuracy: 56.7% ± 0.4%
  - NK F1: 71.5% ± 0.4%

## Detailed Metrics

### Research_Dev
- bcell_metrics: `{
  "accuracy": 0.7,
  "macro_precision": 0.665556831228473,
  "macro_recall": 0.6546658986175116,
  "macro_f1": 0.6593650185702643,
  "total": 300,
  "correct": 210,
  "per_class_accuracy": {
    "contradict": 0.375,
    "support": 0.6693548387096774,
    "unknown": 0.9196428571428571
  },
  "per_class_precision": {
    "contradict": 0.3582089552238806,
    "support": 0.6384615384615384,
    "unknown": 1.0
  },
  "per_class_recall": {
    "contradict": 0.375,
    "support": 0.6693548387096774,
    "unknown": 0.9196428571428571
  },
  "per_class_f1": {
    "contradict": 0.36641221374045807,
    "support": 0.6535433070866141,
    "unknown": 0.9581395348837208
  }
}`
- nk_metrics: `{
  "accuracy": 0.54,
  "precision": 0.573076923076923,
  "recall": 0.8465909090909091,
  "f1": 0.6834862385321101,
  "tp": 149,
  "fp": 111,
  "tn": 13,
  "fn": 27
}`

### Research_Kfold

# IMMUNOS Training Evaluation (2025-12-22 23:58)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: retrieval (max=3)
  - NK mode: classic
  - Train labels: {'unknown': 173, 'support': 173, 'contradict': 173}
  - Test labels: {'unknown': 112, 'support': 124, 'contradict': 64}
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.3% (macro-F1 38.3%)
  - NK accuracy: 51.7% (precision 56.4%, recall 77.3%)

- Research (SciFact K-Fold)
  - Split: kfold-3 (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
- Evidence: retrieval (max=3)
  - NK mode: classic
  - B Cell accuracy: 42.3% ± 2.0%
  - B Cell macro-F1: 38.7% ± 4.0%
  - NK accuracy: 56.3% ± 0.6%
  - NK F1: 70.3% ± 0.6%

## Detailed Metrics

### Research_Dev
- bcell_metrics: `{
  "accuracy": 0.3933333333333333,
  "macro_precision": 0.3859521737951122,
  "macro_recall": 0.3842405913978495,
  "macro_f1": 0.382932865461306,
  "total": 300,
  "correct": 118,
  "per_class_accuracy": {
    "contradict": 0.328125,
    "support": 0.3870967741935484,
    "unknown": 0.4375
  },
  "per_class_precision": {
    "contradict": 0.25925925925925924,
    "support": 0.39344262295081966,
    "unknown": 0.5051546391752577
  },
  "per_class_recall": {
    "contradict": 0.328125,
    "support": 0.3870967741935484,
    "unknown": 0.4375
  },
  "per_class_f1": {
    "contradict": 0.28965517241379307,
    "support": 0.3902439024390244,
    "unknown": 0.4688995215311005
  }
}`
- nk_metrics: `{
  "accuracy": 0.5166666666666667,
  "precision": 0.5643153526970954,
  "recall": 0.7727272727272727,
  "f1": 0.6522781774580335,
  "tp": 136,
  "fp": 105,
  "tn": 19,
  "fn": 40
}`

### Research_Kfold

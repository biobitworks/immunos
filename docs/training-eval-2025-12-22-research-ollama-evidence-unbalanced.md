# IMMUNOS Training Evaluation (2025-12-22 22:49)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=False)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (nomic-embed-text)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train/Test: 809 / 300
  - B Cell accuracy: 41.3%
  - NK accuracy: 53.7% (precision 57.8%, recall 77.8%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.41333333333333333,
  "total": 300,
  "correct": 124,
  "per_class_accuracy": {
    "unknown": 0.0,
    "support": 1.0,
    "contradict": 0.0
  }
}`
- nk_metrics: `{
  "accuracy": 0.5366666666666666,
  "precision": 0.5780590717299579,
  "recall": 0.7784090909090909,
  "f1": 0.6634382566585956,
  "tp": 137,
  "fp": 100,
  "tn": 24,
  "fn": 39
}`

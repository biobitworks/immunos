# IMMUNOS Training Evaluation (2025-12-22 22:31)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 41.3% (precision 50.0%, recall 1.1%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.39,
  "total": 300,
  "correct": 117,
  "per_class_accuracy": {
    "unknown": 0.17857142857142858,
    "support": 0.5806451612903226,
    "contradict": 0.390625
  }
}`
- nk_metrics: `{
  "accuracy": 0.41333333333333333,
  "precision": 0.5,
  "recall": 0.011363636363636364,
  "f1": 0.022222222222222223,
  "tp": 2,
  "fp": 2,
  "tn": 122,
  "fn": 174
}`

# IMMUNOS Training Evaluation (2025-12-22 22:28)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 51.0% (precision 56.5%, recall 71.6%)

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
  "accuracy": 0.51,
  "precision": 0.5650224215246636,
  "recall": 0.7159090909090909,
  "f1": 0.631578947368421,
  "tp": 126,
  "fp": 97,
  "tn": 27,
  "fn": 50
}`

# IMMUNOS Training Evaluation (2025-12-22 22:39)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Evidence: off (max=0)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 52.0% (precision 56.2%, recall 82.4%)

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
  "accuracy": 0.52,
  "precision": 0.562015503875969,
  "recall": 0.8238636363636364,
  "f1": 0.6682027649769585,
  "tp": 145,
  "fp": 113,
  "tn": 11,
  "fn": 31
}`

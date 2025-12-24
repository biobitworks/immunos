# IMMUNOS Training Evaluation (2025-12-22 22:40)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Evidence: off (max=0)
  - NK mode: enhanced (threshold=min_distance, detectors/class=25)
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 41.3% (precision 0.0%, recall 0.0%)

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
  "precision": 0.0,
  "recall": 0.0,
  "f1": 0.0,
  "tp": 0,
  "fp": 0,
  "tn": 124,
  "fn": 176
}`

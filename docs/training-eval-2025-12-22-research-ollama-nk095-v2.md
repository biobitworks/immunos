# IMMUNOS Training Evaluation (2025-12-22 22:39)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Evidence: off (max=0)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 54.0% (precision 57.0%, recall 88.1%)

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
  "accuracy": 0.54,
  "precision": 0.5698529411764706,
  "recall": 0.8806818181818182,
  "f1": 0.6919642857142856,
  "tp": 155,
  "fp": 117,
  "tn": 7,
  "fn": 21
}`

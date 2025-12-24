# IMMUNOS Training Evaluation (2025-12-22 22:40)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Embedder: ollama (nomic-embed-text)
  - Evidence: off (max=0)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 39.0%
  - NK accuracy: 53.7% (precision 56.5%, recall 91.5%)

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
  "accuracy": 0.5366666666666666,
  "precision": 0.5649122807017544,
  "recall": 0.9147727272727273,
  "f1": 0.6984815618221257,
  "tp": 161,
  "fp": 124,
  "tn": 0,
  "fn": 15
}`

# IMMUNOS Training Evaluation (2025-12-22 22:17)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - Train/Test: 519 / 300
  - B Cell accuracy: 33.7%
  - NK accuracy: 54.3% (precision 57.3%, recall 86.9%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.33666666666666667,
  "total": 300,
  "correct": 101,
  "per_class_accuracy": {
    "unknown": 0.33035714285714285,
    "support": 0.28225806451612906,
    "contradict": 0.453125
  }
}`
- nk_metrics: `{
  "accuracy": 0.5433333333333333,
  "precision": 0.5730337078651685,
  "recall": 0.8693181818181818,
  "f1": 0.6907449209932279,
  "tp": 153,
  "fp": 114,
  "tn": 10,
  "fn": 23
}`

# IMMUNOS Training Evaluation (2025-12-22 22:09)

## Summary

- Hallucination (TruthfulQA)
  - Train/Test: 320 / 80
  - B Cell accuracy: 47.5%
  - NK accuracy: 30.0% (precision 23.3%, recall 17.5%)

- Research (SciFact)
  - Train/Test: 239 / 61
  - B Cell accuracy: 41.0%
  - NK accuracy: 52.5% (precision 56.1%, recall 88.9%)

- Network (NSL-KDD)
  - Train/Test: 1000 / 500
  - NK accuracy: 45.8% (precision 0.0%, recall 0.0%)

## Detailed Metrics

### Hallucination
- bcell_metrics: `{
  "accuracy": 0.475,
  "total": 80,
  "correct": 38,
  "per_class_accuracy": {
    "truthful": 0.85,
    "hallucinated": 0.1
  }
}`
- nk_metrics: `{
  "accuracy": 0.3,
  "precision": 0.23333333333333334,
  "recall": 0.175,
  "f1": 0.2,
  "tp": 7,
  "fp": 23,
  "tn": 17,
  "fn": 33
}`

### Research
- bcell_metrics: `{
  "accuracy": 0.4098360655737705,
  "total": 61,
  "correct": 25,
  "per_class_accuracy": {
    "unknown": 0.0,
    "support": 1.0,
    "contradict": 0.0
  }
}`
- nk_metrics: `{
  "accuracy": 0.5245901639344263,
  "precision": 0.5614035087719298,
  "recall": 0.8888888888888888,
  "f1": 0.6881720430107527,
  "tp": 32,
  "fp": 25,
  "tn": 0,
  "fn": 4
}`

### Network
- nk_metrics: `{
  "accuracy": 0.458,
  "precision": 0.0,
  "recall": 0.0,
  "f1": 0.0,
  "tp": 0,
  "fp": 0,
  "tn": 229,
  "fn": 271
}`

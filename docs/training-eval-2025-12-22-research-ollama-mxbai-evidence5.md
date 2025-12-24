# IMMUNOS Training Evaluation (2025-12-22 22:52)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=5)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 70.7%
  - NK accuracy: 53.3% (precision 57.0%, recall 83.5%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.7066666666666667,
  "total": 300,
  "correct": 212,
  "per_class_accuracy": {
    "unknown": 0.9375,
    "support": 0.7096774193548387,
    "contradict": 0.296875
  }
}`
- nk_metrics: `{
  "accuracy": 0.5333333333333333,
  "precision": 0.5697674418604651,
  "recall": 0.8352272727272727,
  "f1": 0.6774193548387097,
  "tp": 147,
  "fp": 111,
  "tn": 13,
  "fn": 29
}`

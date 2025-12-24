# IMMUNOS Training Evaluation (2025-12-22 22:51)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=hybrid, strategy=sha, embed_weight=0.7
  - Embedder: ollama (mxbai-embed-large)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 70.7%
  - NK accuracy: 53.7% (precision 57.1%, recall 84.1%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.7066666666666667,
  "total": 300,
  "correct": 212,
  "per_class_accuracy": {
    "unknown": 0.9375,
    "support": 0.717741935483871,
    "contradict": 0.28125
  }
}`
- nk_metrics: `{
  "accuracy": 0.5366666666666666,
  "precision": 0.5714285714285714,
  "recall": 0.8409090909090909,
  "f1": 0.6804597701149425,
  "tp": 148,
  "fp": 111,
  "tn": 13,
  "fn": 28
}`

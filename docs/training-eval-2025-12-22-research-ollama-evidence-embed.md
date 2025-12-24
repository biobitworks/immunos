# IMMUNOS Training Evaluation (2025-12-22 22:47)

## Summary

- Research (SciFact)
  - Split: train/dev (balanced=True)
  - B Cell: affinity=embedding, strategy=rha, embed_weight=1.0
  - Embedder: ollama (nomic-embed-text)
  - Evidence: on (max=3)
  - NK mode: classic
  - Train/Test: 519 / 300
  - B Cell accuracy: 36.0%
  - NK accuracy: 53.7% (precision 57.2%, recall 83.5%)

## Detailed Metrics

### Research
- bcell_metrics: `{
  "accuracy": 0.36,
  "total": 300,
  "correct": 108,
  "per_class_accuracy": {
    "unknown": 0.008928571428571428,
    "support": 0.6854838709677419,
    "contradict": 0.34375
  }
}`
- nk_metrics: `{
  "accuracy": 0.5366666666666666,
  "precision": 0.5719844357976653,
  "recall": 0.8352272727272727,
  "f1": 0.6789838337182448,
  "tp": 147,
  "fp": 110,
  "tn": 14,
  "fn": 29
}`

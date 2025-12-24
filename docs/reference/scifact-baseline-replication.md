# SciFact Baseline Replication

## Sources
- Metrics: `/Users/byron/projects/data/immunos_data/research/scifact/logs/metrics_dev.json`
- Dataset root: `/Users/byron/projects/data/immunos_data/research/scifact`

## Dev Metrics (from metrics_dev.json)
| Task | Precision | Recall | F1 |
| --- | --- | --- | --- |
| Sentence selection | 0.794 | 0.590 | 0.677 |
| Sentence label | 0.713 | 0.530 | 0.608 |
| Abstract label (only) | 0.910 | 0.675 | 0.775 |
| Abstract label (rationalized) | 0.852 | 0.632 | 0.725 |

## Notes
- Baseline reports sentence selection, sentence labeling, and abstract-level label metrics.
- Metrics are logged in the SciFact replication run under the dataset log directory.

---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/updated_analysis/reports/analysis_results.json
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/updated_analysis/reports/analysis_results.json
generated_at: 2025-12-23 10:28
---

```json
{
  "timestamp": "2025-09-28T01:57:04.945227",
  "dataset_info": {
    "samples": 44,
    "proteins": 5853,
    "ups_proteins_validated": 132
  },
  "claims": {
    "claim1_ups_proteins": {
      "evaluation": "SUPPORTED",
      "total_ups": 132,
      "significantly_changed": 38,
      "percent_significant": 28.8,
      "evidence": "Only 38/132 (28.8%) UPS proteins significantly changed"
    },
    "claim2_sqstm1": {
      "evaluation": "STRONGLY_SUPPORTED",
      "log2_fc": 3.412941152590909,
      "fold_change": 10.65117843635957,
      "pvalue": 9.293505048467219e-08,
      "evidence": "SQSTM1 shows 10.7-fold upregulation (p=9.29e-08)"
    },
    "claim5_mitophagy": {
      "evaluation": "STRONGLY_SUPPORTED",
      "proteins_analyzed": 5,
      "significantly_changed": 3,
      "key_proteins": {
        "SQSTM1": 3.412941152590909,
        "NBR1": 1.4874093532727244
      },
      "evidence": "3/5 mitophagy proteins significantly changed"
    }
  }
}
```

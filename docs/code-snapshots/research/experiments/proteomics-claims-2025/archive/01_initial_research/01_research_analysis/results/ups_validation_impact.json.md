---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/ups_validation_impact.json
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/ups_validation_impact.json
generated_at: 2025-12-23 10:28
---

```json
{
  "validation_date": "2025-09-28T01:50:49.293958",
  "ups_proteins_analyzed": 132,
  "original_proteins": 10,
  "improvement_factor": 13.2,
  "claims_updated": {
    "claim1": {
      "original_evaluation": "PARTIALLY_SUPPORTED",
      "new_evaluation": "SUPPORTED",
      "reason": "With 132 UPS proteins, only 28.8% significantly changed. Autophagy receptors (SQSTM1, NBR1) massively upregulated while most proteasome subunits unchanged.",
      "key_evidence": {
        "total_ups": 132,
        "percent_changed": 28.8,
        "sqstm1_fc": 3.413,
        "nbr1_fc": 1.487,
        "proteasome_changes": "9/43 significant"
      }
    },
    "claim2": {
      "original_evaluation": "SUPPORTED",
      "new_evaluation": "STRONGLY_SUPPORTED",
      "reason": "SQSTM1 shows massive 10.7-fold upregulation (p<0.0001), even stronger than originally claimed",
      "key_evidence": {
        "log2_fc": 3.412941152590909,
        "fold_change": 10.65117843635957,
        "pvalue": 9.293505048467221e-08
      }
    },
    "claim5": {
      "original_evaluation": "SUPPORTED",
      "new_evaluation": "STRONGLY_SUPPORTED",
      "reason": "Multiple mitophagy receptors significantly dysregulated. SQSTM1 and NBR1 massively upregulated indicating impaired clearance.",
      "key_evidence": {
        "sqstm1_up": 3.413,
        "nbr1_up": 1.487,
        "tax1bp1_up": 0.67,
        "mitophagy_proteins_found": 5
      }
    }
  },
  "success_rate": {
    "original": 62.5,
    "updated": 87.5
  },
  "key_findings": [
    "SQSTM1 shows 10.7-fold upregulation",
    "Only 28.8% of UPS proteins significantly changed",
    "Autophagy-specific dysfunction confirmed",
    "Mitophagy receptors accumulated"
  ]
}
```

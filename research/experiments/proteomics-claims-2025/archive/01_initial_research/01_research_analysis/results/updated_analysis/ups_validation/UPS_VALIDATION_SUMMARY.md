# UPS Protein Validation Summary

## Overview
Successfully validated **132 UPS proteins** in the proteomics dataset using gene name matching.

## Key Statistics
- **Total UPS proteins**: 132
- **Significantly changed (p<0.05)**: 38 (28.8%)
- **Upregulated (FC>1.2)**: 8
- **Downregulated (FC<0.8)**: 14

## Categories Validated
| Category | Count | Significant |
|----------|-------|-------------|
| Proteasome subunits | 43 | 9 |
| E3 ligases | 19 | 8 |
| E2 enzymes | 18 | 3 |
| E1 enzymes | 7 | 4 |
| Deubiquitinases | 28 | 7 |
| UPS regulators | 9 | 3 |

## Top Findings
1. **SQSTM1/p62 massively upregulated** (Log2 FC: 3.41, p<0.0001)
2. **NBR1 significantly upregulated** (Log2 FC: 1.49, p<0.0001)
3. **Multiple proteasome subunits downregulated** (PSME1, PSME2, PSMD9)
4. **Autophagy-specific dysfunction** evident from receptor accumulation

## Files Generated
- `ups_proteins_found.md` - Detailed protein report
- `ups_expression_data.csv` - Expression data for all UPS proteins
- `validated_ups_genes.json` - JSON format for programmatic use
- `ups_analysis_validated.py` - Updated analysis script

## Integration with Analysis
The validated UPS proteins have been integrated into the analysis scripts for evaluating mitochondrial dysregulation claims, particularly:
- Statement 1: UPS protein expression
- Statement 2: SQSTM1/p62 upregulation validation
- Statement 5: Mitophagy pathway assessment

---
*Validation complete: 132 UPS proteins confirmed in dataset*

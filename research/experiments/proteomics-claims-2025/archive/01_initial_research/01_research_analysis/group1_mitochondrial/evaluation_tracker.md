# Finding Group Evaluation Tracker

## Project Information
- **Dataset**: pool_processed_v2.h5ad (Proteomic dataset - Alzheimer's disease)
- **Context PDFs**:
  - Late-Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure.pdf
  - Sequential_Failure_of_Proteostasis_Mechanisms.pdf
- **Colab Notebook**: https://colab.research.google.com/drive/1qeVYENLzDEnU-gb0jOV5Sx2qmu-SRHbq?usp=sharing

## Dataset Characteristics
- Mini pools of 10 neurons from AD cases
- MC1 quantification (misfolded tau)
- Pseudotime ordering
- Metadata: age at death, MC1 score, pseudotime, tau status (positive/negative)
- Data is log2 transformed

---

## ANALYSIS STATEMENTS (8 statements)

### Statement 1: UPS Proteins
**Claim**: Targeted analyses show no significant UPS protein alterations across tau-positive versus tau-negative neurons.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 2: SQSTM1 Upregulation
**Claim**: SQSTM1 (p62) is massively upregulated (log2FC = 3.413, FDR = 1.76 × 10^-8) and increases with pseudotime (β = 4.951, FDR < 0.001).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 3: SQSTM1 with Autophagy Proteins
**Claim**: The upregulation SQSTM1 (p62) is accompanied by upregulation of BECN1 and CTSD but downregulation of ATG12, ULK1, and CTSL, while none of 11 UPS proteins are significantly changed.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 4: Protein Correlations
**Claim**: Of TAX1BP1, CAT, VDAC1, CYCS, ATP5F1A, UQCRC2, COX4I1, PRDX1, KEAP1, and TFRC, all correlate with SQSTM1 and pseudotime (p < 0.05), with 8/10 and 5/10 remaining significant after Bonferroni correction (α = 0.0025), respectively.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 5: SQSTM1-VDAC1 Global Correlation
**Claim**: The global correlation between SQSTM1 and VDAC1 across the full dataset is negligible (r = 0.0536, p = 0.730).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 6: Running Correlation Analysis
**Claim**: A running correlation between SQSTM1 and VDAC along pseudotime (sliding window n = 20) shifts from significantly negative early (mean r = -0.417 for pseudotime < 0.33) to positive late (mean r = 0.478 for pseudotime > 0.67), with a strong trend (r = 0.851, p = 6.98 × 10^-8).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 7: CYCS Biphasic Pattern
**Claim**: CYCS (cytochrome c) shows a biphasic pattern: stable/elevated at lower MC1 (MC1 < 2.5: 15.259 ± 0.167) but sharply reduced at high MC1 (≥ 3.0: 14.798 ± 0.191; t = 5.474, p = 5.097 × 10^-5, Cohen's d = -2.58).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 8: CYCS-VATPase Synchrony
**Claim**: CYCS and V-ATPase decline in synchrony (CYCS–VATPase r = 0.696, p < 0.001).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

## INTERPRETATION STATEMENTS (3 statements)

### Interpretation 1: Late-stage Dysregulation
**Context**: Stage-resolved metrics and dynamic coupling analyses together with running-correlation analyses show the SQSTM1–VDAC1 relationship shifting from negative early to positive late.

**Claim**: These results collectively support a late-stage mitochondrial dysregulation with mitophagy failure that is synchronized with lysosomal decompensation.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

### Interpretation 2: [Title]
**Claim**: [To be provided]

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

### Interpretation 3: [Title]
**Claim**: [To be provided]

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

## Summary Statistics
- Total Statements Evaluated: ___/11
- Analysis Statements: ___/8
  - Supported: ___
  - Refuted: ___
  - Unsure: ___
- Interpretation Statements: ___/3
  - Supported: ___
  - Refuted: ___
  - Unsure: ___

## Notes
_[Any additional observations or challenges encountered during evaluation]_
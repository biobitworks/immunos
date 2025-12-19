# Finding Group 2 Evaluation Tracker

## Project Information
- **Dataset**: pool_processed_v2.h5ad (Proteomic dataset - Alzheimer's disease)
- **Context PDF**: Sequential_Failure_of_Proteostasis_Mechanisms.pdf
- **Colab Notebook**: https://colab.research.google.com/drive/1xKk6BlJmjIA0iw3mIrJ-WJvOJwt8Ci1F?usp=sharing

## Dataset Characteristics
- Mini pools of 10 neurons from AD cases
- MC1 quantification (misfolded tau)
- Pseudotime ordering
- Metadata: age at death, PMI, PatientID, MC1 score, pseudotime, tau status
- Data is log2 transformed
- 5,853 proteins analyzed

---

## ANALYSIS STATEMENTS (7 statements)

### Statement 1: Differential Expression Overview
**Claim**: A covariate-controlled differential expression analysis (age, PMI, and PatientID) across 5,853 proteins (BH-FDR) identified 2,115 proteins (36.14%) significantly altered between tau-positive and tau-negative neurons.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your covariate-controlled DE analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 2: SQSTM1 Top Upregulated
**Claim**: The autophagy receptor SQSTM1 was the most upregulated protein (+3.41 log2) in the covariate-controlled differential expression analysis.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 3: Collagen Decreases
**Claim**: Marked decreases were observed in the covariate-controlled differential expression analysis with collagens (COL1A1, COL1A2, COL6A2) showing the largest declines (> −4.0 log2).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 4: V-ATPase-MC1 Correlations
**Claim**: Within tau-positive cells, individual V-ATPase subunits (ATP6V1A/B2/C1/H) correlated strongly and negatively with MC1 (mean r = −0.723; range −0.581 to −0.791; all p < 0.005).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 5: V-ATPase Score Segmented Analysis
**Claim**: The "V-ATPase Score" was calculated as the arithmetic mean of the log2-scaled expression values of the detected V-ATPase subunit genes. A segmented analysis of this V-ATPase Score against MC1 located a critical threshold at MC1 = 2.831 where a near-stable phase (slope −0.0029) transitions into a steep decline (slope −0.4255).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your segmented regression analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 6: V-ATPase Biphasic Behavior
**Claim**: The V-ATPase system showed significant biphasic behavior, but with a later breakpoint (0.654) than the proteasome (difference 0.282 pseudotime units), and an early decline followed by an apparent late "recovery" (slopes −0.202 then +0.304).

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your biphasic analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 7: V-ATPase vs Proteasome Timing
**Claim**: Along the disease trajectory, the V-ATPase breakpoint occurs later (pseudotime ~0.65) than the proteasome's, placing lysosomal failure downstream of proteasomal dysfunction.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Analysis Code**:
```python
# Your breakpoint comparison analysis here
```

**Explanation**:
_[Your detailed explanation with evidence]_

---

### Statement 8: [Empty]
**Claim**: [No claim provided]

**Evaluation**: N/A

**Analysis Code**:
```python
# N/A
```

**Explanation**:
_Statement 8 was not provided_

---

## INTERPRETATION STATEMENTS (3 statements)

### Interpretation 1: Failed-Compensation Model
**Context**: Individual V-ATPase subunits (ATP6V1A/B2/C1/H) are correlated strongly and negatively with MC1, with expression highest at low MC1 and progressively reduced as misfolded-tau burden rose.

**Claim**: These results support a failed-compensation model under increasing pathology.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

### Interpretation 2: Selective Survival
**Context**: The V-ATPase system also showed significant biphasic behavior, but with a later breakpoint (0.654) than the proteasome (difference 0.282 pseudotime units), and an early decline followed by an apparent late "recovery".

**Claim**: These results reflect selective survival rather than restoration of function, linking acidification competence to neuronal persistence.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

### Interpretation 3: Progressive Dysfunction
**Context**: Along the disease trajectory, an early proteasome inflection followed by a later collapse of lysosomal acidification coincides with a broad proteomic downshift, including reduced SLC2A1 (glucose transport), LMNA (nuclear architecture), and collagens.

**Claim**: These results associate tau accumulation to progressive dysfunction and differential survival outcomes.

**Evaluation**: [ SUPPORTED / REFUTED / UNSURE ]

**Explanation**:
_[Your interpretation based on report context only - no data analysis]_

---

## Summary Statistics
- Total Statements Evaluated: ___/10
- Analysis Statements: ___/7
  - Supported: ___
  - Refuted: ___
  - Unsure: ___
- Interpretation Statements: ___/3
  - Supported: ___
  - Refuted: ___
  - Unsure: ___

## Notes
_[Any additional observations or challenges encountered during evaluation]_
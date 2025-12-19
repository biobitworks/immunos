# 🏋️ Practice Exercise 5: Complete Integration Challenge

## 🎯 Learning Goal
Integrate all your skills to perform a comprehensive proteomics analysis from raw data to biological insights.

---

## 📋 Your Mission

You're a biologist studying Alzheimer's disease. Your lab just received proteomics data from 150 brain samples. Your PI asks: **"What's happening to the protein quality control system in tau pathology, and when does it fail?"**

### Your Resources
- Proteomics data (H5AD format)
- Sample metadata (age, sex, tau status, pseudotime)
- Your newly acquired analysis skills!

---

## 🔬 Part 1: Initial Data Analysis (30 min)

### Task 1.1: Data Loading and QC

```python
# Your code here:
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# TODO: Load and explore data
adata = sc.read_h5ad('alzheimer_proteomics.h5ad')

# TODO: Basic QC checks
# - Check dimensions
# - Look for missing values
# - Examine metadata
# - Check for outliers

print(f"Data shape: {adata.shape}")
print(f"Metadata columns: {list(adata.obs.columns)}")
```

**Initial observations:**
- [ ] Total proteins: _____
- [ ] Total samples: _____
- [ ] Missing data issues: _____
- [ ] Potential outliers: _____

### Task 1.2: Focus on UPS/Autophagy Proteins

```python
# Your code here:
# Protein quality control systems
ups_proteins = ['PSMA1', 'PSMA2', 'PSMA3', 'PSMB1', 'PSMB2', 'PSMB5', 'PSMD1', 'PSMD2']
autophagy_proteins = ['SQSTM1', 'MAP1LC3B', 'BECN1', 'ATG5', 'ATG7', 'GABARAPL2']
chaperones = ['HSPA1A', 'HSPA5', 'HSP90AA1', 'HSPB1']

qc_proteins = ups_proteins + autophagy_proteins + chaperones

# TODO: Extract these proteins
# TODO: Check which are present in data
# TODO: Create focused dataset
```

**Quality control proteins found:**
- [ ] UPS proteins present: _____ / 8
- [ ] Autophagy proteins present: _____ / 6
- [ ] Chaperones present: _____ / 4

---

## 🔬 Part 2: Differential Expression Analysis (45 min)

### Task 2.1: Simple Comparison

```python
# Your code here:
# TODO: Compare tau+ vs tau- for QC proteins
# Perform t-tests
# Calculate fold changes
# Apply FDR correction

de_results = []
for protein in qc_proteins:
    if protein in adata.var_names:
        # TODO: Implement DE analysis
        pass

de_df = pd.DataFrame(de_results)
```

### Task 2.2: Adjusted Analysis

```python
# Your code here:
# TODO: Repeat with covariate adjustment
# Adjust for age, sex, PMI
# Compare results to simple analysis

from statsmodels.formula.api import ols

adjusted_results = []
for protein in qc_proteins:
    if protein in adata.var_names:
        # TODO: Implement adjusted analysis
        pass
```

**DE Results Summary:**

| System | Upregulated | Downregulated | Unchanged |
|--------|-------------|---------------|-----------|
| UPS | | | |
| Autophagy | | | |
| Chaperones | | | |

---

## 🔬 Part 3: Temporal Dynamics (45 min)

### Task 3.1: Pseudotime Analysis

```python
# Your code here:
# TODO: Order cells by disease progression
# Use provided pseudotime or calculate

# TODO: Track QC protein changes over pseudotime
for protein in ['SQSTM1', 'PSMA1', 'HSPA5']:
    if protein in adata.var_names:
        # TODO: Plot expression vs pseudotime
        # TODO: Fit smoothed trajectory
        pass
```

### Task 3.2: System Coordination

```python
# Your code here:
# TODO: Calculate average expression for each system
ups_avg = []  # Average UPS expression over pseudotime
autophagy_avg = []  # Average autophagy expression
chaperone_avg = []  # Average chaperone expression

# TODO: Sliding window correlation between systems
# When do they become coordinated?
```

**Temporal Findings:**
- [ ] UPS peaks at pseudotime: _____
- [ ] Autophagy peaks at pseudotime: _____
- [ ] Systems become coupled at: _____
- [ ] Critical transition at: _____

---

## 🔬 Part 4: Pathway Integration (30 min)

### Task 4.1: Pathway Enrichment

```python
# Your code here:
# TODO: Take all significant proteins from DE
# TODO: Perform GO enrichment
# TODO: Check for protein degradation pathways

import gseapy as gp

# For upregulated proteins
up_proteins = de_df[de_df['log2fc'] > 0.5]['protein'].tolist()
go_up = gp.enrichr(gene_list=up_proteins,
                   gene_sets='GO_Biological_Process_2021',
                   outdir=None)
```

### Task 4.2: Network Analysis

```python
# Your code here:
# TODO: Build protein interaction network
# TODO: Identify modules
# TODO: Find hub proteins

import networkx as nx

G = nx.Graph()
# TODO: Add QC proteins as nodes
# TODO: Add interactions (from STRING or correlation)
# TODO: Identify central proteins
```

**Network insights:**
- [ ] Most connected protein: _____
- [ ] Number of modules: _____
- [ ] Key interactions: _____

---

## 🔬 Part 5: Biological Interpretation (30 min)

### Task 5.1: Timeline of Failure

Based on your analyses, create a timeline:

```python
# Your interpretation:
timeline = {
    'Early (0-0.3)': 'Your interpretation here',
    'Middle (0.3-0.6)': 'Your interpretation here',
    'Late (0.6-1.0)': 'Your interpretation here'
}

# TODO: Support with specific data
# Which proteins change when?
# What systems fail first?
```

### Task 5.2: Mechanistic Model

```python
# Draw your model:
"""
Your mechanistic model here:

Healthy → [What happens?] → Stressed → [What happens?] → Failed

Key proteins involved at each stage:
1. _______________
2. _______________
3. _______________
"""
```

---

## 🔬 Part 6: Clinical Translation (20 min)

### Task 6.1: Biomarker Candidates

```python
# Your code here:
# TODO: Identify best biomarker candidates
# Criteria: Large effect size, robust to adjustment, early change

biomarker_scores = []
for protein in qc_proteins:
    if protein in adata.var_names:
        score = 0
        # TODO: Score based on:
        # - Effect size (Cohen's d)
        # - Significance after adjustment
        # - Early pseudotime change
        # - Low variance within groups
        pass
```

**Top 3 Biomarker Candidates:**
1. Protein: _____ | Score: _____ | Rationale: _____
2. Protein: _____ | Score: _____ | Rationale: _____
3. Protein: _____ | Score: _____ | Rationale: _____

### Task 6.2: Therapeutic Windows

```python
# Based on your temporal analysis:
therapeutic_windows = {
    'Prevention': 'Pseudotime 0-___: Target _____',
    'Early intervention': 'Pseudotime ___-___: Target _____',
    'Late stage': 'Pseudotime ___-1.0: Target _____'
}

print("Therapeutic recommendations:")
for window, recommendation in therapeutic_windows.items():
    print(f"{window}: {recommendation}")
```

---

## 📊 Part 7: Scientific Communication (20 min)

### Task 7.1: Create Publication Figure

```python
# Your code here:
# TODO: Create a 4-panel figure summarizing findings

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Volcano plot of QC proteins
ax = axes[0, 0]
# TODO: Volcano plot highlighting UPS, autophagy, chaperones

# Panel B: Temporal dynamics
ax = axes[0, 1]
# TODO: Line plot of system averages over pseudotime

# Panel C: Network at early vs late
ax = axes[1, 0]
# TODO: Network visualization

# Panel D: Model diagram
ax = axes[1, 1]
# TODO: Schematic of your mechanistic model

plt.suptitle('Protein Quality Control Failure in Tau Pathology')
plt.tight_layout()
```

### Task 7.2: Write Abstract

Write a 150-word abstract summarizing your findings:

```
Title: _________________________________

Background: ____________________________
___________________________________________

Methods: _________________________________
___________________________________________

Results: _________________________________
___________________________________________
___________________________________________

Conclusions: _____________________________
___________________________________________
___________________________________________
```

---

## 🎯 Integration Checklist

Confirm you've successfully integrated all skills:

**Data Analysis:**
- [ ] Loaded and QC'd proteomics data
- [ ] Handled missing values appropriately
- [ ] Identified outliers

**Statistics:**
- [ ] Performed appropriate statistical tests
- [ ] Applied multiple testing correction
- [ ] Calculated effect sizes
- [ ] Adjusted for covariates

**Temporal Analysis:**
- [ ] Analyzed changes over pseudotime
- [ ] Identified critical transitions
- [ ] Tracked system coordination

**Biological Interpretation:**
- [ ] Connected statistics to biology
- [ ] Proposed mechanistic model
- [ ] Identified clinical implications

**Communication:**
- [ ] Created clear visualizations
- [ ] Wrote scientific summary
- [ ] Made actionable recommendations

---

## 🤔 Reflection Questions

1. **What was the most surprising finding?**
   Your answer: _______________

2. **What additional data would help confirm your model?**
   Your answer: _______________

3. **How would you validate these findings experimentally?**
   Your answer: _______________

4. **What are the clinical implications?**
   Your answer: _______________

---

## 🎓 Self-Assessment Rubric

Rate yourself (1-5) on each skill:

| Skill | Score | Evidence |
|-------|-------|----------|
| Data manipulation | _/5 | |
| Statistical analysis | _/5 | |
| Temporal analysis | _/5 | |
| Network analysis | _/5 | |
| Biological interpretation | _/5 | |
| Scientific communication | _/5 | |

**Total: ___/30**

- 25-30: Expert level! Ready for independent research
- 20-24: Proficient - Minor gaps to address
- 15-19: Competent - Review specific topics
- <15: Keep practicing - Review earlier exercises

---

## 🚀 What's Next?

### If you scored 25-30:
- Apply these skills to your own research
- Explore advanced methods (machine learning, multi-omics)
- Consider contributing to open-source tools

### If you scored 20-24:
- Review sections where you struggled
- Try the analysis with different parameters
- Practice with additional datasets

### If you scored <20:
- Revisit earlier exercises
- Focus on fundamentals
- Don't hesitate to review tutorials

---

## 📝 Solution Overview

<details>
<summary>Click to see solution approach (try yourself first!)</summary>

### Expected Key Findings:

1. **Differential Expression:**
   - SQSTM1: 2.5-fold increase (p < 0.001)
   - Proteasome subunits: Mixed (some up, some down)
   - HSP70/90: Increased (stress response)

2. **Temporal Pattern:**
   - Early: Proteasome upregulation (compensation)
   - Middle: Autophagy activation (SQSTM1 surge)
   - Late: System collapse (everything dysregulated)

3. **Critical Transition:**
   - Occurs around pseudotime 0.45
   - Marked by system coupling
   - Point of no return for intervention

4. **Mechanistic Model:**
   ```
   Healthy neurons
        ↓ (tau accumulation)
   Proteasome stress (early compensation)
        ↓ (overwhelmed)
   Autophagy activation (SQSTM1 up)
        ↓ (capacity exceeded)
   System coupling (loss of independence)
        ↓
   Complete failure (cell death pathway)
   ```

5. **Clinical Translation:**
   - Biomarkers: SQSTM1 (early), PSMA1 (progression)
   - Therapeutic window: Before pseudotime 0.45
   - Target: Boost proteasome early, support autophagy middle

### Technical Excellence Indicators:
- Proper multiple testing correction
- Covariate adjustment performed
- Bootstrap confidence intervals
- Network metrics calculated
- Clear visualizations
- Biological plausibility
</details>

---

**🎉 Congratulations! You've completed the comprehensive proteomics analysis challenge!**

You've successfully integrated data analysis, statistics, temporal dynamics, and biological interpretation to answer a real research question. These are the exact skills you'll use in actual proteomics research.

**Remember:** Every expert was once a beginner. You've come incredibly far!

---

[Return to Exercise Overview](README.md) | [Return to Main Guide](../README.md)
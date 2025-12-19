# 🏋️ Practice Exercise 1: Basic Differential Expression Analysis

## 🎯 Learning Goal
Practice performing a basic differential expression analysis comparing tau-positive and tau-negative neurons.

---

## 📋 Your Task

You've been given proteomics data from Alzheimer's disease brain samples. Your goal is to identify proteins that differ between tau-positive and tau-negative neurons.

### Dataset Information
- **Samples**: 150 neurons (75 tau-positive, 75 tau-negative)
- **Proteins**: Focus on top 100 most variable proteins
- **Format**: H5AD file with expression matrix and metadata

---

## 🔬 Exercise Steps

### Step 1: Data Loading and Exploration
Load the data and answer these questions:
1. How many total proteins are in the dataset?
2. What is the distribution of tau-positive vs tau-negative samples?
3. What metadata columns are available?

```python
# Your code here:
import scanpy as sc
import pandas as pd
import numpy as np

# Load data
adata = sc.read_h5ad('path_to_data.h5ad')

# TODO: Complete the exploration
```

**Questions to answer:**
- [ ] Total number of proteins: _____
- [ ] Number of tau+ samples: _____
- [ ] Number of tau- samples: _____
- [ ] List metadata columns: _____

### Step 2: Quality Control
Check data quality:
1. Are there any missing values?
2. What is the distribution of total protein expression per sample?
3. Are there obvious outliers?

```python
# Your code here:
# TODO: Check for missing values
# TODO: Plot distribution of total expression
# TODO: Identify potential outliers
```

### Step 3: Simple Differential Expression
Perform t-tests for the top 10 most variable proteins:

```python
# Your code here:
from scipy import stats

# Get top 10 variable proteins
# Hint: Calculate variance for each protein across samples

top_proteins = []  # TODO: Fill this

results = []
for protein in top_proteins:
    # TODO: Extract expression for this protein
    # TODO: Split by tau status
    # TODO: Perform t-test
    # TODO: Store results
    pass

# Create results dataframe
```

**Fill in your results:**

| Protein | Mean Tau+ | Mean Tau- | T-statistic | P-value | Significant? |
|---------|-----------|-----------|-------------|---------|--------------|
| | | | | | |
| | | | | | |

### Step 4: Multiple Testing Correction
Apply FDR correction to your p-values:

```python
# Your code here:
from statsmodels.stats.multitest import fdrcorrection

# TODO: Apply FDR correction
# TODO: How many proteins remain significant?
```

**Questions:**
- [ ] How many proteins were significant before correction?
- [ ] How many remain significant after FDR correction?
- [ ] What does this tell you about multiple testing?

### Step 5: Effect Size Calculation
Calculate Cohen's d for significant proteins:

```python
# Your code here:
def calculate_cohens_d(group1, group2):
    # TODO: Implement Cohen's d
    pass

# TODO: Calculate for significant proteins
```

### Step 6: Visualization
Create a volcano plot:

```python
# Your code here:
import matplotlib.pyplot as plt

# TODO: Create volcano plot
# X-axis: log2 fold change
# Y-axis: -log10(p-value)
# Color significant points differently
```

---

## 🤔 Reflection Questions

1. **Why do we need multiple testing correction?**
   Your answer: _______________

2. **What's the difference between statistical significance and biological relevance?**
   Your answer: _______________

3. **How would batch effects impact your results?**
   Your answer: _______________

4. **What additional analyses would strengthen your findings?**
   Your answer: _______________

---

## 🎯 Success Criteria

Your analysis is complete when you can:
- [ ] Load and explore proteomics data
- [ ] Perform basic quality control
- [ ] Conduct differential expression analysis
- [ ] Apply multiple testing correction
- [ ] Calculate effect sizes
- [ ] Create informative visualizations
- [ ] Interpret results biologically

---

## 💡 Hints

<details>
<summary>Hint 1: Finding Variable Proteins</summary>

```python
# Calculate variance for each protein
protein_variances = np.var(adata.X, axis=0)
# Get indices of top 10
top_indices = np.argsort(protein_variances)[-10:]
top_proteins = adata.var_names[top_indices]
```
</details>

<details>
<summary>Hint 2: T-test Implementation</summary>

```python
tau_pos = adata[adata.obs['tau_status'] == 'tau_positive', protein].X.flatten()
tau_neg = adata[adata.obs['tau_status'] == 'tau_negative', protein].X.flatten()
t_stat, p_val = stats.ttest_ind(tau_pos, tau_neg)
```
</details>

<details>
<summary>Hint 3: Cohen's d Formula</summary>

```python
def calculate_cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    return (np.mean(group1) - np.mean(group2)) / pooled_std
```
</details>

---

## 📊 Expected Outcomes

If done correctly, you should find:
- Several proteins with significant differences (p < 0.05)
- Fewer proteins significant after FDR correction
- A mix of upregulated and downregulated proteins
- Effect sizes ranging from small to large

---

## 🚀 Extension Challenges

### Challenge 1: Bootstrap Confidence Intervals
Instead of t-tests, use bootstrap to calculate confidence intervals for the difference in means.

### Challenge 2: Permutation Testing
Implement a permutation test to validate your t-test results.

### Challenge 3: Power Analysis
Calculate the statistical power of your tests given the sample size and effect sizes.

---

## 📝 Solution Overview

<details>
<summary>Click to see solution approach (try yourself first!)</summary>

### Key Steps:
1. **Data Loading**: Use scanpy to load H5AD file
2. **QC**: Check for missing values, outliers using boxplots
3. **DE Analysis**: T-tests with proper handling of multiple comparisons
4. **Visualization**: Volcano plot highlighting significant proteins
5. **Interpretation**: Consider both statistical and biological significance

### Common Pitfalls:
- Forgetting multiple testing correction
- Ignoring effect sizes
- Not checking assumptions (normality, equal variance)
- Over-interpreting p-values

### Best Practices:
- Always visualize your data first
- Report effect sizes alongside p-values
- Use appropriate statistical tests
- Document your analysis steps
</details>

---

## 📚 Related Resources

- [Statistical Methods Guide](../06_resources_and_support/statistical_methods_reference.md)
- [Visualization Best Practices](../06_resources_and_support/visualization_guide.md)
- [Troubleshooting Common Issues](../06_resources_and_support/troubleshooting_guide.md)

---

**Ready for the next challenge? Continue to [Exercise 2: Pathway Analysis](exercise_2_pathway_analysis.md)**
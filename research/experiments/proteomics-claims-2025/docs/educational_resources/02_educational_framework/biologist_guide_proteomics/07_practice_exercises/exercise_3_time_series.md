# 🏋️ Practice Exercise 3: Temporal Dynamics Analysis

## 🎯 Learning Goal
Learn to analyze how protein relationships change over disease progression using sliding window correlation analysis.

---

## 📋 Your Task

Analyze the temporal dynamics of protein-protein correlations across disease pseudotime, focusing on the SQSTM1-VDAC1 relationship and expanding to network-level analysis.

### Dataset Context
- **Samples**: 150 neurons ordered by disease pseudotime (0 = healthy, 1 = severe)
- **Focus proteins**: SQSTM1 (autophagy), VDAC1 (mitochondria)
- **Goal**: Track correlation changes over disease progression

---

## 🔬 Exercise Steps

### Step 1: Data Preparation and Pseudotime Ordering

```python
# Your code here:
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
import matplotlib.pyplot as plt

# Load data
adata = sc.read_h5ad('path_to_data.h5ad')

# TODO: Order samples by pseudotime
# Hint: Use disease markers to infer progression
# Options: diffusion pseudotime, monocle, or provided pseudotime

pseudotime = []  # TODO: Calculate or extract pseudotime

# Sort data by pseudotime
sort_idx = np.argsort(pseudotime)
adata_sorted = adata[sort_idx, :]
```

**Questions:**
- [ ] How did you determine pseudotime ordering?
- [ ] What is the distribution of pseudotime values?
- [ ] Do tau+ and tau- samples separate along pseudotime?

### Step 2: Sliding Window Correlation

Implement sliding window correlation analysis:

```python
# Your code here:
def sliding_window_correlation(expr1, expr2, pseudotime, window_size=0.3, step=0.05):
    """
    Calculate correlation in sliding windows across pseudotime
    """
    correlations = []
    window_centers = []
    n_samples = []

    # TODO: Implement sliding window
    # Start from window_size/2, end at 1 - window_size/2
    # For each window position:
    #   - Select samples within window
    #   - Calculate correlation
    #   - Store results

    return window_centers, correlations, n_samples

# Extract SQSTM1 and VDAC1 expression
sqstm1_expr = adata_sorted[:, 'SQSTM1'].X.flatten()
vdac1_expr = adata_sorted[:, 'VDAC1'].X.flatten()

# Calculate sliding window correlations
centers, corrs, samples = sliding_window_correlation(
    sqstm1_expr, vdac1_expr, pseudotime
)
```

**Record your findings:**
- [ ] Correlation at early pseudotime (0-0.2): _____
- [ ] Correlation at middle pseudotime (0.4-0.6): _____
- [ ] Correlation at late pseudotime (0.8-1.0): _____

### Step 3: Statistical Testing

Test for significant correlation changes:

```python
# Your code here:
# TODO: Test if correlation changes are significant

# Option 1: Mann-Kendall trend test
from scipy.stats import kendalltau
tau, p_value = kendalltau(centers, corrs)
print(f"Trend test: tau={tau:.3f}, p={p_value:.4f}")

# Option 2: Compare early vs late correlations
early_samples = adata_sorted[pseudotime < 0.3, :]
late_samples = adata_sorted[pseudotime > 0.7, :]

# TODO: Bootstrap confidence intervals for difference
def bootstrap_correlation_diff(early_data, late_data, n_bootstrap=1000):
    # TODO: Implement bootstrap
    pass
```

**Statistical Results:**
- [ ] Is there a significant trend? (p-value): _____
- [ ] Early vs Late difference: _____ (95% CI: [____, ____])

### Step 4: Identify Critical Transitions

Find points where correlation changes rapidly:

```python
# Your code here:
# TODO: Calculate rate of change in correlation
correlation_changes = np.diff(corrs) / np.diff(centers)

# TODO: Identify change points
# Option 1: Find maximum rate of change
# Option 2: Use change point detection algorithm

critical_points = []  # TODO: Identify critical transitions

# Visualize
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Correlation trajectory
ax = axes[0]
ax.plot(centers, corrs, 'b-', linewidth=2)
ax.fill_between(centers, corrs - 0.1, corrs + 0.1, alpha=0.3)  # Add CI
ax.set_ylabel('Correlation')
ax.set_title('SQSTM1-VDAC1 Correlation Dynamics')

# Rate of change
ax = axes[1]
ax.plot(centers[1:], correlation_changes, 'r-')
ax.set_xlabel('Pseudotime')
ax.set_ylabel('Rate of Change')

plt.tight_layout()
plt.show()
```

**Critical Transitions:**
- [ ] Number of transitions identified: _____
- [ ] Pseudotime of major transition: _____
- [ ] Biological interpretation: _____

### Step 5: Expand to Multiple Protein Pairs

Analyze temporal dynamics for multiple protein pairs:

```python
# Your code here:
# Key protein pairs to analyze
protein_pairs = [
    ('SQSTM1', 'VDAC1'),  # Autophagy-Mitochondria
    ('SQSTM1', 'MAP1LC3B'),  # Autophagy markers
    ('VDAC1', 'ATP5A1'),  # Mitochondrial proteins
    ('MAPT', 'SQSTM1'),  # Tau-Autophagy
]

# TODO: Calculate correlations for all pairs
pair_dynamics = {}
for protein1, protein2 in protein_pairs:
    # TODO: Apply sliding window analysis
    # TODO: Store results
    pass

# TODO: Create heatmap showing correlation evolution
# Rows = protein pairs, Columns = pseudotime windows
```

**Multi-protein Findings:**
| Protein Pair | Early Corr | Late Corr | Change | Pattern |
|--------------|------------|-----------|--------|---------|
| SQSTM1-VDAC1 | | | | |
| SQSTM1-MAP1LC3B | | | | |
| VDAC1-ATP5A1 | | | | |
| MAPT-SQSTM1 | | | | |

### Step 6: Network-Level Analysis

Analyze how the entire correlation network evolves:

```python
# Your code here:
# TODO: Build correlation networks at different timepoints
timepoints = [0.2, 0.5, 0.8]
networks = []

for tp in timepoints:
    # Select samples near this timepoint
    window_samples = np.abs(pseudotime - tp) < 0.15
    subset = adata_sorted[window_samples, :]

    # TODO: Calculate correlation matrix
    # TODO: Threshold to create network
    # TODO: Calculate network metrics

    network_metrics = {
        'density': 0,  # TODO: Calculate
        'clustering': 0,  # TODO: Calculate
        'modularity': 0,  # TODO: Calculate
    }
    networks.append(network_metrics)

# TODO: Plot network evolution
```

**Network Evolution:**
- [ ] Early network density: _____
- [ ] Late network density: _____
- [ ] Main change in network structure: _____

---

## 🤔 Analysis Questions

1. **What biological process does the SQSTM1-VDAC1 correlation represent?**
   Your answer: _______________

2. **Why might correlations increase during disease progression?**
   Your answer: _______________

3. **What does a critical transition point indicate biologically?**
   Your answer: _______________

4. **How would you validate these temporal patterns?**
   Your answer: _______________

---

## 💡 Hints

<details>
<summary>Hint 1: Pseudotime Calculation</summary>

```python
# Using diffusion pseudotime
sc.pp.neighbors(adata)
sc.tl.diffmap(adata)
sc.tl.dpt(adata, root=healthy_cell_index)
pseudotime = adata.obs['dpt_pseudotime']
```
</details>

<details>
<summary>Hint 2: Sliding Window Implementation</summary>

```python
def sliding_window_correlation(expr1, expr2, pseudotime, window_size=0.3, step=0.05):
    window_centers = np.arange(window_size/2, 1 - window_size/2, step)
    correlations = []

    for center in window_centers:
        # Get samples in window
        in_window = np.abs(pseudotime - center) <= window_size/2
        if np.sum(in_window) >= 10:  # Minimum samples
            corr, _ = stats.spearmanr(expr1[in_window], expr2[in_window])
            correlations.append(corr)
        else:
            correlations.append(np.nan)

    return window_centers, correlations
```
</details>

<details>
<summary>Hint 3: Change Point Detection</summary>

```python
# Using ruptures package
import ruptures as rpt

# Detect change points
algo = rpt.Pelt(model="rbf").fit(corrs)
change_points = algo.predict(pen=10)

# Or simple derivative method
grad = np.gradient(corrs)
peaks = np.where(np.abs(grad) > np.std(grad) * 2)[0]
```
</details>

---

## 🎯 Success Criteria

You've mastered this exercise when you can:
- [ ] Order samples by disease progression
- [ ] Implement sliding window correlation
- [ ] Identify significant temporal trends
- [ ] Detect critical transition points
- [ ] Analyze multiple protein relationships
- [ ] Interpret biological meaning

---

## 📊 Expected Results

You should observe:
- **Progressive correlation increase** between SQSTM1-VDAC1
- **Critical transition** around pseudotime 0.4-0.5
- **Network densification** over disease progression
- **Module coupling** in late disease

---

## 🚀 Extension Challenges

### Challenge 1: Bootstrapped Confidence Bands
Add confidence intervals to your correlation trajectory using bootstrap.

### Challenge 2: Multivariate Dynamics
Use partial correlations to control for confounding proteins.

### Challenge 3: Trajectory Clustering
Identify proteins with similar temporal patterns using clustering.

---

## 📝 Solution Overview

<details>
<summary>Click to see solution approach (try yourself first!)</summary>

### Key Findings:

1. **SQSTM1-VDAC1 Dynamics:**
   - Early: r ≈ 0.15 (independent function)
   - Middle: r ≈ 0.45 (coordinating)
   - Late: r ≈ 0.75 (locked together)
   - Critical transition at pseudotime ≈ 0.45

2. **Biological Interpretation:**
   - Early: Normal autophagy and mitochondrial function
   - Middle: Stress response coordination
   - Late: System failure and coupling

3. **Network Changes:**
   - Density increases 3-fold
   - Modular structure lost
   - Hub proteins emerge

### Technical Notes:
- Use Spearman correlation (robust to outliers)
- Window size affects smoothness vs resolution
- Bootstrap for confidence intervals
- Multiple testing correction for many pairs
</details>

---

**Great progress! Continue to [Exercise 4: Covariate Adjustment](exercise_4_covariate_adjustment.md)**
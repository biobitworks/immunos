# Statistical Methods Reference Guide

## Overview
This guide provides detailed explanations of all statistical methods used in the proteomic analysis evaluations, with connections to ISLP (Introduction to Statistical Learning with Python) concepts.

---

## 1. Differential Expression Analysis

### Basic T-Test
**Purpose**: Compare mean expression between two groups (tau+ vs tau-)

```python
from scipy.stats import ttest_ind
t_statistic, p_value = ttest_ind(group1, group2)
```

**Assumptions**:
- Normal distribution (check with Q-Q plots)
- Equal variance (use Welch's t-test if violated)
- Independent observations

**ISLP Connection**: Chapter 13 - Hypothesis Testing

### Log2 Fold Change
**Formula**: `log2FC = mean(log2(tau+)) - mean(log2(tau-))`

**Interpretation**:
- log2FC = 1 → 2-fold increase
- log2FC = 3.413 → 2^3.413 = 10.7-fold increase
- log2FC = -1 → 2-fold decrease

**Note**: If data is already log2-transformed, simply subtract means

---

## 2. Multiple Testing Correction

### False Discovery Rate (FDR) - Benjamini-Hochberg

**Problem**: Testing 5,853 proteins → expect 293 false positives at α=0.05

**Solution**: Control FDR instead of family-wise error rate

```python
from statsmodels.stats.multitest import multipletests
rejected, pvals_corrected, _, _ = multipletests(pvalues, method='fdr_bh', alpha=0.05)
```

**Steps**:
1. Sort p-values: p(1) ≤ p(2) ≤ ... ≤ p(m)
2. Find largest i where: p(i) ≤ (i/m) × α
3. Reject hypotheses 1, ..., i

**Interpretation**: FDR = 0.05 means 5% of discoveries are expected to be false

**ISLP Connection**: Chapter 13.3 - Multiple Testing

---

## 3. Covariate-Controlled Analysis

### Linear Model with Covariates

**Model**: `Expression ~ TauStatus + Age + PMI + PatientID`

```python
from statsmodels.formula.api import ols
model = ols('expression ~ tau_status + age + PMI + C(PatientID)', data=df)
results = model.fit()
tau_effect = results.params['tau_status[T.positive]']
```

**Why Use Covariates**:
- Remove confounding effects
- Increase statistical power
- Get unbiased estimate of tau effect

**Interpretation of Coefficients**:
- β₁ (tau_status): Effect of tau after controlling for other variables
- β₂ (age): Change in expression per year of age
- β₃ (PMI): Change in expression per hour post-mortem

**ISLP Connection**: Chapter 3 - Linear Regression

---

## 4. Correlation Analysis

### Pearson Correlation
**Formula**: `r = Σ[(xi - x̄)(yi - ȳ)] / √[Σ(xi - x̄)² × Σ(yi - ȳ)²]`

```python
from scipy.stats import pearsonr
correlation, p_value = pearsonr(x, y)
```

**Assumptions**:
- Linear relationship
- Normal distribution
- No outliers

**Interpretation**:
- r = 0.7 to 1.0: Strong positive
- r = 0.3 to 0.7: Moderate positive
- r = -0.3 to 0.3: Weak/no relationship
- r = -0.7 to -0.3: Moderate negative
- r = -1.0 to -0.7: Strong negative

### Spearman Correlation
**Use when**: Data is not normally distributed or relationship is monotonic but not linear

```python
from scipy.stats import spearmanr
rho, p_value = spearmanr(x, y)
```

---

## 5. Sliding Window Analysis

### Running Correlation
**Purpose**: Detect dynamic relationship changes over time/pseudotime

```python
def sliding_window_correlation(data, window_size=20, step=1):
    correlations = []
    for i in range(0, len(data) - window_size + 1, step):
        window = data[i:i+window_size]
        corr = pearsonr(window.x, window.y)[0]
        correlations.append(corr)
    return correlations
```

**Key Parameters**:
- Window size: Balance between stability and resolution
- Step size: Overlap between windows

**Interpretation**:
- Changing correlation sign: Relationship reversal
- Increasing correlation: Strengthening relationship
- Correlation trend: Overall trajectory

**ISLP Connection**: Chapter 7 - Moving Beyond Linearity (Local Methods)

---

## 6. Segmented Regression (Piecewise Linear)

### Breakpoint Detection
**Purpose**: Find threshold where relationship changes

```python
from scipy.optimize import curve_fit

def piecewise_linear(x, x0, y0, b1, b2):
    return np.piecewise(x, [x < x0],
                        [lambda x: y0 + b1*(x-x0),
                         lambda x: y0 + b2*(x-x0)])

params, _ = curve_fit(piecewise_linear, x_data, y_data)
breakpoint = params[0]
```

**Interpretation**:
- x0: Breakpoint location
- b1: Slope before breakpoint
- b2: Slope after breakpoint

**Applications**:
- Threshold effects
- Phase transitions
- Non-linear relationships

---

## 7. Effect Size Measures

### Cohen's d
**Formula**: `d = (mean1 - mean2) / pooled_std`

```python
pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
cohens_d = (mean1 - mean2) / pooled_std
```

**Interpretation**:
- |d| < 0.2: Negligible
- |d| = 0.2-0.5: Small
- |d| = 0.5-0.8: Medium
- |d| > 0.8: Large

**Advantage**: Independent of sample size (unlike p-values)

---

## 8. Bootstrap Confidence Intervals

### Non-parametric Confidence Intervals

```python
def bootstrap_ci(data, statistic, n_boot=10000, alpha=0.05):
    bootstrap_stats = []
    for _ in range(n_boot):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrap_stats.append(statistic(sample))

    ci_lower = np.percentile(bootstrap_stats, 100*alpha/2)
    ci_upper = np.percentile(bootstrap_stats, 100*(1-alpha/2))
    return ci_lower, ci_upper
```

**Advantages**:
- No distribution assumptions
- Works for any statistic
- Robust to outliers

**ISLP Connection**: Chapter 5 - Resampling Methods

---

## 9. Model Evaluation Metrics

### R-squared (Coefficient of Determination)
**Formula**: `R² = 1 - (SS_residual / SS_total)`

**Interpretation**:
- R² = 0.7: Model explains 70% of variance
- Higher is better (but beware overfitting)

### AIC (Akaike Information Criterion)
**Formula**: `AIC = 2k - 2ln(L)`

**Use**: Model selection (lower is better)

---

## 10. Bonferroni Correction

### Conservative Multiple Testing Adjustment

```python
bonferroni_alpha = 0.05 / n_tests
significant = p_value < bonferroni_alpha
```

**When to Use**:
- Small number of tests
- Need strict control of Type I error

**Example**: Testing 10 proteins against 2 variables = 20 tests
- Bonferroni α = 0.05/20 = 0.0025

---

## 11. Temporal/Pseudotime Analysis

### Linear Regression with Pseudotime

```python
import statsmodels.api as sm
X = sm.add_constant(pseudotime)
model = sm.OLS(expression, X)
results = model.fit()
beta = results.params[1]  # Slope (rate of change)
```

**Interpretation of β**:
- β > 0: Expression increases over disease progression
- β < 0: Expression decreases over disease progression
- |β|: Rate of change

---

## 12. Statistical Power Considerations

### Sample Size Effects

**Factors Affecting Power**:
1. Effect size (larger → more power)
2. Sample size (more → more power)
3. Significance level (higher α → more power)
4. Variance (lower → more power)

**Rule of Thumb**: Need ~20 samples per group for moderate effect (d=0.5)

---

## Common Pitfalls and Solutions

### 1. P-hacking
**Problem**: Testing until significant
**Solution**: Pre-specify analyses, use FDR correction

### 2. Outliers
**Problem**: Single points driving results
**Solution**: Use robust methods (median, MAD), check with/without outliers

### 3. Non-independence
**Problem**: Multiple measurements from same patient
**Solution**: Mixed effects models, cluster-robust standard errors

### 4. Missing Data
**Problem**: Biased results if not random
**Solution**: Multiple imputation, sensitivity analysis

### 5. Batch Effects
**Problem**: Technical variation masking biological signal
**Solution**: Include batch as covariate, ComBat correction

---

## Evaluation Guidelines

### When to Call "SUPPORTED"
- Statistical significance achieved (FDR < 0.05)
- Effect size within expected range (±20%)
- Multiple lines of evidence converge

### When to Call "REFUTED"
- No statistical significance
- Effect size opposite direction or >50% different
- Contradictory evidence

### When to Call "UNSURE"
- Missing data/proteins
- Borderline significance
- Technical issues prevent analysis

---

## Python Libraries Summary

```python
# Core statistical testing
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu, pearsonr, spearmanr

# Multiple testing correction
from statsmodels.stats.multitest import multipletests

# Linear models
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Machine learning
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor

# Data manipulation
import pandas as pd
import numpy as np

# Single-cell analysis
import scanpy as sc

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
```

---

## References

1. **ISLP Book**: James, Witten, Hastie, Tibshirani, Taylor (2023) "An Introduction to Statistical Learning with Applications in Python"

2. **Multiple Testing**: Benjamini & Hochberg (1995) "Controlling the false discovery rate"

3. **Effect Sizes**: Cohen (1988) "Statistical Power Analysis for the Behavioral Sciences"

4. **Linear Models**: Faraway (2014) "Linear Models with R"

5. **Proteomics Statistics**: Kammers et al. (2015) "Detecting significant changes in protein abundance"

6. **Single-Cell Analysis**: Wolf et al. (2018) "SCANPY: large-scale single-cell gene expression data analysis"
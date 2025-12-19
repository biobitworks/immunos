# 🏋️ Practice Exercise 4: Covariate-Adjusted Analysis

## 🎯 Learning Goal
Master the critical skill of adjusting for confounding variables to reveal true biological signals in your proteomics data.

---

## 📋 Your Task

You've discovered that several proteins appear differentially expressed between tau+ and tau- neurons. However, your colleague points out that tau+ samples are from older patients. Are your findings real or just age effects? Let's find out!

### The Challenge
- Initial finding: 50 proteins significantly different (p < 0.05)
- Confounding factor: Age differs between groups (tau+: mean age 78, tau-: mean age 65)
- Question: Which proteins are truly associated with tau pathology?

---

## 🔬 Exercise Steps

### Step 1: Explore the Confounding Problem

```python
# Your code here:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Load data
data = pd.read_csv('proteomics_with_covariates.csv')  # Adjust path

# TODO: Check group balance
print("Group Statistics:")
print(data.groupby('tau_status')[['age', 'sex', 'pmi']].describe())

# TODO: Visualize the confounding
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Age distribution by group
ax = axes[0]
# TODO: Create violin plot of age by tau status

# Example protein vs age
ax = axes[1]
# TODO: Scatter plot of a protein vs age, colored by tau status
```

**Observations:**
- [ ] Age difference between groups: _____ years
- [ ] Other imbalanced covariates: _____
- [ ] Visual evidence of confounding: _____

### Step 2: Simple vs Adjusted Comparison

Compare results with and without adjustment:

```python
# Your code here:
# Select a protein known to change with age (e.g., GFAP)
protein_name = 'GFAP'  # Astrocyte marker

# TODO: Simple t-test (unadjusted)
tau_pos = data[data['tau_status'] == 1][protein_name]
tau_neg = data[data['tau_status'] == 0][protein_name]
t_stat_simple, p_val_simple = stats.ttest_ind(tau_pos, tau_neg)

print(f"Simple t-test: t={t_stat_simple:.3f}, p={p_val_simple:.4f}")

# TODO: Linear model with age adjustment
model_adjusted = ols(f'{protein_name} ~ tau_status + age', data=data).fit()
print("\nAdjusted model:")
print(model_adjusted.summary())

# Extract adjusted effect
tau_effect_adjusted = model_adjusted.params['tau_status']
tau_pval_adjusted = model_adjusted.pvalues['tau_status']

print(f"\nComparison:")
print(f"Unadjusted effect: {tau_pos.mean() - tau_neg.mean():.3f}")
print(f"Adjusted effect: {tau_effect_adjusted:.3f}")
```

**Results for GFAP:**
- [ ] Unadjusted p-value: _____
- [ ] Adjusted p-value: _____
- [ ] Change in effect size: _____
- [ ] Interpretation: _____

### Step 3: Systematic Covariate-Adjusted DE

Perform adjusted analysis for all proteins:

```python
# Your code here:
proteins = [col for col in data.columns if col not in ['tau_status', 'age', 'sex', 'pmi', 'batch']]

results = []
for protein in proteins[:100]:  # Start with 100 proteins
    try:
        # TODO: Fit adjusted model
        formula = f'{protein} ~ tau_status + age + sex + pmi'
        model = ols(formula, data=data).fit()

        # TODO: Also fit unadjusted model for comparison
        model_simple = ols(f'{protein} ~ tau_status', data=data).fit()

        # TODO: Store results
        results.append({
            'protein': protein,
            'unadj_coef': model_simple.params['tau_status'],
            'unadj_pval': model_simple.pvalues['tau_status'],
            'adj_coef': model.params['tau_status'],
            'adj_pval': model.pvalues['tau_status'],
            'age_coef': model.params['age'],
            'age_pval': model.pvalues['age'],
            'r_squared': model.rsquared
        })
    except:
        continue

results_df = pd.DataFrame(results)

# TODO: Apply FDR correction
from statsmodels.stats.multitest import fdrcorrection
results_df['unadj_padj'] = fdrcorrection(results_df['unadj_pval'])[1]
results_df['adj_padj'] = fdrcorrection(results_df['adj_pval'])[1]
```

**Summary Statistics:**
- [ ] Significant before adjustment (FDR < 0.05): _____
- [ ] Significant after adjustment (FDR < 0.05): _____
- [ ] Lost significance after adjustment: _____
- [ ] Gained significance after adjustment: _____

### Step 4: Categorize Proteins by Adjustment Impact

```python
# Your code here:
# TODO: Categorize proteins
def categorize_protein(row):
    if row['unadj_padj'] < 0.05 and row['adj_padj'] < 0.05:
        return 'Robust to adjustment'
    elif row['unadj_padj'] < 0.05 and row['adj_padj'] >= 0.05:
        return 'Confounded (lost significance)'
    elif row['unadj_padj'] >= 0.05 and row['adj_padj'] < 0.05:
        return 'Masked (gained significance)'
    else:
        return 'Not significant'

results_df['category'] = results_df.apply(categorize_protein, axis=1)

# TODO: Visualize categories
category_counts = results_df['category'].value_counts()
plt.figure(figsize=(8, 6))
# Create pie chart or bar plot
```

**Protein Categories:**
| Category | Count | Examples |
|----------|-------|----------|
| Robust to adjustment | | |
| Confounded | | |
| Masked by confounding | | |
| Not significant | | |

### Step 5: Diagnostic Checks

Verify model assumptions:

```python
# Your code here:
# Pick a significant protein after adjustment
example_protein = results_df[results_df['adj_padj'] < 0.05].iloc[0]['protein']

# TODO: Refit model for diagnostics
model = ols(f'{example_protein} ~ tau_status + age + sex + pmi', data=data).fit()

# TODO: Check assumptions
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Residuals vs Fitted
ax = axes[0, 0]
ax.scatter(model.fittedvalues, model.resid, alpha=0.5)
ax.axhline(y=0, color='red', linestyle='--')
ax.set_xlabel('Fitted Values')
ax.set_ylabel('Residuals')
ax.set_title('Residuals vs Fitted')

# Q-Q plot
ax = axes[0, 1]
stats.probplot(model.resid, dist="norm", plot=ax)
ax.set_title('Q-Q Plot')

# TODO: Add other diagnostic plots

# Check for multicollinearity
from statsmodels.stats.outliers_influence import variance_inflation_factor
X = data[['age', 'sex', 'pmi']]
vif = pd.DataFrame()
vif['Variable'] = X.columns
vif['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
print("\nVariance Inflation Factors:")
print(vif)
```

**Diagnostic Results:**
- [ ] Residuals normally distributed? _____
- [ ] Homoscedasticity satisfied? _____
- [ ] Any high VIF values (>10)? _____
- [ ] Model assumptions met? _____

### Step 6: Sensitivity Analysis

Test robustness to covariate selection:

```python
# Your code here:
# TODO: Leave-one-covariate-out analysis
sensitivity_results = []
covariates = ['age', 'sex', 'pmi']

for protein in results_df.nsmallest(10, 'adj_pval')['protein']:
    protein_sensitivity = {'protein': protein}

    # Full model
    full_model = ols(f'{protein} ~ tau_status + age + sex + pmi', data=data).fit()
    protein_sensitivity['full_coef'] = full_model.params['tau_status']

    # Leave each covariate out
    for leave_out in covariates:
        remaining = [c for c in covariates if c != leave_out]
        formula = f'{protein} ~ tau_status + {" + ".join(remaining)}'
        reduced_model = ols(formula, data=data).fit()
        protein_sensitivity[f'without_{leave_out}'] = reduced_model.params['tau_status']

    sensitivity_results.append(protein_sensitivity)

sensitivity_df = pd.DataFrame(sensitivity_results)

# TODO: Visualize sensitivity
# Create heatmap or bar plots showing coefficient changes
```

**Sensitivity Findings:**
- [ ] Most influential covariate: _____
- [ ] Proteins stable across models: _____
- [ ] Proteins sensitive to covariate choice: _____

---

## 🤔 Interpretation Questions

1. **Why did some proteins lose significance after adjustment?**
   Your answer: _______________

2. **What does it mean when a protein gains significance after adjustment?**
   Your answer: _______________

3. **Should we always adjust for every available covariate?**
   Your answer: _______________

4. **How do you decide which covariates to include?**
   Your answer: _______________

5. **What's the difference between a confounder and a mediator?**
   Your answer: _______________

---

## 💡 Hints

<details>
<summary>Hint 1: Checking Balance</summary>

```python
# Check continuous variables
for var in ['age', 'pmi', 'rin']:
    tau_pos = data[data['tau_status'] == 1][var]
    tau_neg = data[data['tau_status'] == 0][var]
    t_stat, p_val = stats.ttest_ind(tau_pos, tau_neg)
    print(f"{var}: Tau+ = {tau_pos.mean():.1f}, Tau- = {tau_neg.mean():.1f}, p = {p_val:.4f}")

# Check categorical variables
print(pd.crosstab(data['sex'], data['tau_status'], normalize='columns'))
```
</details>

<details>
<summary>Hint 2: Interpreting Adjusted Coefficients</summary>

```python
# The tau_status coefficient in adjusted model represents:
# The difference in protein expression between tau+ and tau-
# AFTER accounting for differences in age, sex, and PMI

# If coefficient changes dramatically after adjustment:
# Original effect was partly/wholly due to confounders

# If coefficient stays similar:
# Effect is independent of confounders (robust finding)
```
</details>

<details>
<summary>Hint 3: Model Selection</summary>

```python
# Compare models using AIC/BIC
models = {
    'Simple': ols(f'{protein} ~ tau_status', data=data).fit(),
    'Age only': ols(f'{protein} ~ tau_status + age', data=data).fit(),
    'Full': ols(f'{protein} ~ tau_status + age + sex + pmi', data=data).fit()
}

for name, model in models.items():
    print(f"{name}: AIC = {model.aic:.1f}, BIC = {model.bic:.1f}")
# Lower AIC/BIC = better model
```
</details>

---

## 🎯 Success Criteria

You've mastered covariate adjustment when you can:
- [ ] Identify confounding variables
- [ ] Build and interpret adjusted models
- [ ] Compare adjusted vs unadjusted results
- [ ] Perform diagnostic checks
- [ ] Conduct sensitivity analyses
- [ ] Decide which covariates to include

---

## 📊 Expected Insights

After proper adjustment, you should find:
- **~30-40% of proteins** lose significance (were age-confounded)
- **~5-10% of proteins** gain significance (suppressor effects)
- **~50% of proteins** remain robust
- **Age** is typically the strongest confounder
- **Batch effects** can be substantial if present

---

## 🚀 Extension Challenges

### Challenge 1: Interaction Effects
Test if the tau effect differs by age or sex (tau × age interaction).

### Challenge 2: Non-linear Relationships
Use splines or polynomial terms for age if relationship is non-linear.

### Challenge 3: Propensity Score Matching
Instead of regression adjustment, try matching tau+ and tau- samples on covariates.

---

## 📝 Solution Overview

<details>
<summary>Click to see solution approach (try yourself first!)</summary>

### Key Findings:

1. **Major Confounders:**
   - Age accounts for 30-50% of apparent differences
   - Batch effects present but smaller (~10%)
   - Sex has minimal confounding effect

2. **Protein Categories:**
   - Robust: True tau-related proteins (SQSTM1, MAP2)
   - Confounded: Age-related proteins (GFAP, inflammatory markers)
   - Masked: Proteins where age opposed tau effect

3. **Best Practices:**
   - Always check covariate balance
   - Include strong confounders (age, batch)
   - Don't over-adjust (too many covariates)
   - Validate with sensitivity analysis

### Common Mistakes to Avoid:
- Ignoring obvious confounders
- Adjusting for mediators
- Including highly correlated covariates
- Not checking model assumptions
</details>

---

**Excellent work on mastering covariate adjustment! Continue to [Exercise 5: Integration Challenge](exercise_5_integration.md)**
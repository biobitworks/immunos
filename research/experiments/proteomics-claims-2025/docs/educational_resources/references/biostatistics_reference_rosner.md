# Biostatistics Reference: Rosner's Fundamentals Applied to Proteomic Analysis

## Overview
This guide integrates concepts from Bernard Rosner's "Fundamentals of Biostatistics" (8th Edition) with the proteomic analysis evaluations. Rosner's text provides the theoretical foundation for many statistical methods used in biological research.

---

## Chapter Connections to Our Analysis

### Chapter 5: Estimation (Rosner)
**Application**: Confidence intervals for fold changes and correlations

```python
# Example: 95% CI for log2 fold change
import numpy as np
from scipy import stats

def calculate_ci_mean_difference(group1, group2, alpha=0.05):
    """
    Calculate confidence interval for difference of means
    Following Rosner Chapter 5
    """
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)

    # Pooled variance (assuming equal variances)
    pooled_var = ((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)
    se_diff = np.sqrt(pooled_var * (1/n1 + 1/n2))

    # t-distribution with df = n1 + n2 - 2
    df = n1 + n2 - 2
    t_critical = stats.t.ppf(1 - alpha/2, df)

    diff = mean1 - mean2
    ci_lower = diff - t_critical * se_diff
    ci_upper = diff + t_critical * se_diff

    return diff, (ci_lower, ci_upper)
```

**Rosner Concept**: Two-sample t-test confidence intervals
**Our Application**: SQSTM1 fold change confidence intervals

---

### Chapter 6: Hypothesis Testing (Rosner)
**Application**: All differential expression analyses

**Key Concepts from Rosner**:
1. **Type I Error (α)**: False positive rate
2. **Type II Error (β)**: False negative rate
3. **Power (1-β)**: Probability of detecting true effect
4. **P-value interpretation**: Probability of observing data given H₀ is true

```python
def power_analysis_two_sample(effect_size, n1, n2, alpha=0.05):
    """
    Power analysis for two-sample t-test
    Based on Rosner Chapter 6
    """
    from scipy import stats
    import numpy as np

    # Non-centrality parameter
    ncp = effect_size * np.sqrt(n1 * n2 / (n1 + n2))

    # Critical value
    df = n1 + n2 - 2
    t_critical = stats.t.ppf(1 - alpha/2, df)

    # Power calculation
    power = 1 - stats.nct.cdf(t_critical, df, ncp) + stats.nct.cdf(-t_critical, df, ncp)

    return power
```

**Our Application**: Determining if we have sufficient power to detect claimed fold changes

---

### Chapter 7: One-Sample Inference (Rosner)
**Application**: Testing if correlations differ from expected values

```python
def test_correlation_hypothesis(observed_r, expected_r, n, alpha=0.05):
    """
    Test if observed correlation differs from expected
    Using Fisher's z-transformation (Rosner Chapter 7)
    """
    import numpy as np
    from scipy import stats

    # Fisher's z-transformation
    z_obs = 0.5 * np.log((1 + observed_r) / (1 - observed_r))
    z_exp = 0.5 * np.log((1 + expected_r) / (1 - expected_r))

    # Standard error
    se_z = 1 / np.sqrt(n - 3)

    # Test statistic
    z_stat = (z_obs - z_exp) / se_z

    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

    return z_stat, p_value
```

**Our Application**: Testing if observed correlations match expected values (e.g., r = 0.851 in sliding window)

---

### Chapter 8: Two-Sample Inference (Rosner)
**Application**: Comparing tau-positive vs tau-negative neurons

**Rosner's Key Points**:
1. **Paired vs Independent samples**: Our data is independent
2. **Equal vs Unequal variances**: Use Welch's t-test if variances differ
3. **Normal vs Non-normal data**: Use Mann-Whitney U if non-normal

```python
def comprehensive_two_sample_test(group1, group2):
    """
    Comprehensive two-sample testing following Rosner principles
    """
    from scipy import stats
    import numpy as np

    results = {}

    # 1. Test for equal variances (Levene's test)
    levene_stat, levene_p = stats.levene(group1, group2)
    equal_var = levene_p > 0.05

    # 2. Test for normality (Shapiro-Wilk)
    _, norm_p1 = stats.shapiro(group1)
    _, norm_p2 = stats.shapiro(group2)
    normal_data = (norm_p1 > 0.05) and (norm_p2 > 0.05)

    # 3. Choose appropriate test
    if normal_data:
        if equal_var:
            # Student's t-test
            t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=True)
            test_used = "Student's t-test"
        else:
            # Welch's t-test
            t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=False)
            test_used = "Welch's t-test"
    else:
        # Mann-Whitney U test
        u_stat, p_val = stats.mannwhitneyu(group1, group2, alternative='two-sided')
        test_used = "Mann-Whitney U test"

    results['test_used'] = test_used
    results['p_value'] = p_val
    results['equal_variance'] = equal_var
    results['normal_data'] = normal_data

    return results
```

---

### Chapter 10: Regression and Correlation (Rosner)
**Application**: Pseudotime correlation analysis and covariate adjustment

**Key Rosner Concepts**:
1. **Simple Linear Regression**: Y = α + βX + ε
2. **Multiple Regression**: Y = β₀ + β₁X₁ + β₂X₂ + ... + ε
3. **Correlation vs Causation**: Correlation doesn't imply causation

```python
def regression_diagnostics(X, y):
    """
    Regression diagnostics following Rosner Chapter 10
    """
    import matplotlib.pyplot as plt
    from scipy import stats
    import numpy as np

    # Fit regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)

    # Predicted values and residuals
    y_pred = slope * X + intercept
    residuals = y - y_pred

    # Diagnostic plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Residuals vs fitted
    axes[0,0].scatter(y_pred, residuals)
    axes[0,0].axhline(y=0, color='red', linestyle='--')
    axes[0,0].set_xlabel('Fitted Values')
    axes[0,0].set_ylabel('Residuals')
    axes[0,0].set_title('Residuals vs Fitted')

    # 2. Q-Q plot of residuals
    stats.probplot(residuals, dist="norm", plot=axes[0,1])
    axes[0,1].set_title('Q-Q Plot of Residuals')

    # 3. Scale-location plot
    axes[1,0].scatter(y_pred, np.sqrt(np.abs(residuals)))
    axes[1,0].set_xlabel('Fitted Values')
    axes[1,0].set_ylabel('√|Residuals|')
    axes[1,0].set_title('Scale-Location Plot')

    # 4. Original data with regression line
    axes[1,1].scatter(X, y, alpha=0.5)
    axes[1,1].plot(X, y_pred, 'red', linewidth=2)
    axes[1,1].set_xlabel('X')
    axes[1,1].set_ylabel('Y')
    axes[1,1].set_title(f'Regression Line (R² = {r_value**2:.3f})')

    plt.tight_layout()
    plt.show()

    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'p_value': p_value,
        'std_error': std_err
    }
```

**Our Application**: SQSTM1 vs pseudotime regression (β = 4.951)

---

### Chapter 11: Multiple Regression (Rosner)
**Application**: Covariate-controlled differential expression

**Rosner's Key Teaching Points**:
1. **Confounding variables**: Variables that affect both predictor and outcome
2. **Partial correlation**: Correlation after removing effects of other variables
3. **Model selection**: Forward, backward, stepwise selection

```python
def covariate_adjusted_analysis(expression, tau_status, age, pmi, patient_id):
    """
    Covariate-adjusted analysis following Rosner Chapter 11
    """
    import pandas as pd
    import statsmodels.api as sm
    from statsmodels.formula.api import ols

    # Create dataframe
    df = pd.DataFrame({
        'expression': expression,
        'tau_status': tau_status,
        'age': age,
        'pmi': pmi,
        'patient_id': patient_id
    })

    # Model without covariates (crude analysis)
    model_crude = ols('expression ~ tau_status', data=df)
    results_crude = model_crude.fit()

    # Model with covariates (adjusted analysis)
    model_adjusted = ols('expression ~ tau_status + age + pmi + C(patient_id)', data=df)
    results_adjusted = model_adjusted.fit()

    # Extract coefficients
    crude_effect = results_crude.params['tau_status[T.positive]']
    adjusted_effect = results_adjusted.params['tau_status[T.positive]']

    # Confounding assessment
    percent_change = 100 * (adjusted_effect - crude_effect) / crude_effect

    print(f"Crude tau effect: {crude_effect:.3f}")
    print(f"Adjusted tau effect: {adjusted_effect:.3f}")
    print(f"Percent change: {percent_change:.1f}%")

    if abs(percent_change) > 10:
        print("Significant confounding detected!")

    return results_crude, results_adjusted
```

---

### Chapter 12: Analysis of Variance (ANOVA) (Rosner)
**Application**: Comparing multiple groups (e.g., MC1 tertiles)

```python
def mc1_tertile_analysis(expression, mc1_scores):
    """
    Compare expression across MC1 tertiles using ANOVA
    Following Rosner Chapter 12
    """
    import numpy as np
    from scipy import stats

    # Create tertiles
    tertile_cutoffs = np.percentile(mc1_scores, [33.33, 66.67])

    group1 = expression[mc1_scores <= tertile_cutoffs[0]]  # Low MC1
    group2 = expression[(mc1_scores > tertile_cutoffs[0]) &
                       (mc1_scores <= tertile_cutoffs[1])]  # Medium MC1
    group3 = expression[mc1_scores > tertile_cutoffs[1]]   # High MC1

    # One-way ANOVA
    f_stat, p_value = stats.f_oneway(group1, group2, group3)

    # Post-hoc pairwise comparisons (if ANOVA significant)
    if p_value < 0.05:
        # Bonferroni correction for 3 comparisons
        alpha_corrected = 0.05 / 3

        t1, p1 = stats.ttest_ind(group1, group2)
        t2, p2 = stats.ttest_ind(group1, group3)
        t3, p3 = stats.ttest_ind(group2, group3)

        print(f"ANOVA: F = {f_stat:.3f}, p = {p_value:.3e}")
        print(f"Low vs Medium: p = {p1:.3e} {'*' if p1 < alpha_corrected else ''}")
        print(f"Low vs High: p = {p2:.3e} {'*' if p2 < alpha_corrected else ''}")
        print(f"Medium vs High: p = {p3:.3e} {'*' if p3 < alpha_corrected else ''}")

    return f_stat, p_value
```

**Our Application**: CYCS expression across MC1 groups

---

### Chapter 13: Nonparametric Methods (Rosner)
**Application**: When normality assumptions are violated

**Rosner's Nonparametric Arsenal**:
1. **Wilcoxon Signed-Rank**: Paired data, non-normal
2. **Mann-Whitney U**: Independent samples, non-normal
3. **Kruskal-Wallis**: Multiple groups, non-normal
4. **Spearman Correlation**: Non-linear relationships

```python
def robust_differential_expression(group1, group2):
    """
    Robust DE analysis using nonparametric methods
    Following Rosner Chapter 13
    """
    from scipy import stats
    import numpy as np

    # Test normality
    _, norm_p1 = stats.shapiro(group1[:50])  # Shapiro works best with n < 50
    _, norm_p2 = stats.shapiro(group2[:50])

    if norm_p1 > 0.05 and norm_p2 > 0.05:
        # Use parametric test
        t_stat, p_val = stats.ttest_ind(group1, group2)
        test_type = "t-test"
        effect_size = (np.mean(group1) - np.mean(group2)) / np.sqrt((np.var(group1) + np.var(group2))/2)
    else:
        # Use nonparametric test
        u_stat, p_val = stats.mannwhitneyu(group1, group2, alternative='two-sided')
        test_type = "Mann-Whitney U"
        # Rank-biserial correlation as effect size
        n1, n2 = len(group1), len(group2)
        effect_size = 1 - (2 * u_stat) / (n1 * n2)

    return {
        'test_type': test_type,
        'p_value': p_val,
        'effect_size': effect_size,
        'normal_data': norm_p1 > 0.05 and norm_p2 > 0.05
    }
```

---

### Chapter 16: Logistic Regression (Rosner)
**Application**: Predicting tau status from protein expression

```python
def tau_status_prediction(protein_expression, tau_status):
    """
    Predict tau status using logistic regression
    Following Rosner Chapter 16
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import classification_report, roc_auc_score
    import numpy as np

    # Convert tau status to binary
    y = (tau_status == 'positive').astype(int)
    X = protein_expression.reshape(-1, 1)

    # Fit logistic regression
    model = LogisticRegression()
    model.fit(X, y)

    # Predictions
    y_pred_proba = model.predict_proba(X)[:, 1]
    y_pred = model.predict(X)

    # Performance metrics
    auc = roc_auc_score(y, y_pred_proba)

    print(f"AUC: {auc:.3f}")
    print(classification_report(y, y_pred))

    return model, auc
```

**Our Application**: How well can SQSTM1 expression predict tau status?

---

## Study Design Considerations (Rosner Chapter 4)

### Sample Size Calculation
```python
def sample_size_two_sample(effect_size, power=0.8, alpha=0.05):
    """
    Sample size calculation for two-sample t-test
    Following Rosner Chapter 4
    """
    from scipy import stats
    import numpy as np

    # Critical values
    z_alpha = stats.norm.ppf(1 - alpha/2)
    z_beta = stats.norm.ppf(power)

    # Sample size per group
    n = 2 * ((z_alpha + z_beta) / effect_size) ** 2

    return int(np.ceil(n))

# Example: Detect fold change of 1.5 (log2FC = 0.585)
n_needed = sample_size_two_sample(0.585)
print(f"Need {n_needed} samples per group to detect 1.5-fold change")
```

---

## Rosner's Approach to P-value Interpretation

### The Rosner P-value Framework:
- **p < 0.001**: Very strong evidence against H₀
- **0.001 ≤ p < 0.01**: Strong evidence against H₀
- **0.01 ≤ p < 0.05**: Moderate evidence against H₀
- **0.05 ≤ p < 0.10**: Weak evidence against H₀
- **p ≥ 0.10**: Little or no evidence against H₀

### Applied to Our Statements:
```python
def interpret_p_value_rosner(p_value):
    """
    Interpret p-values using Rosner's framework
    """
    if p_value < 0.001:
        return "Very strong evidence against H₀"
    elif p_value < 0.01:
        return "Strong evidence against H₀"
    elif p_value < 0.05:
        return "Moderate evidence against H₀"
    elif p_value < 0.10:
        return "Weak evidence against H₀"
    else:
        return "Little or no evidence against H₀"

# Example: SQSTM1 p-value interpretation
sqstm1_p = 1.76e-8
print(f"SQSTM1 p-value ({sqstm1_p:.2e}): {interpret_p_value_rosner(sqstm1_p)}")
```

---

## Key Biostatistical Principles from Rosner

### 1. Always Check Assumptions
```python
def check_test_assumptions(data1, data2=None):
    """
    Comprehensive assumption checking
    """
    from scipy import stats

    assumptions = {}

    # Normality (Shapiro-Wilk)
    _, p_norm1 = stats.shapiro(data1[:50])
    assumptions['group1_normal'] = p_norm1 > 0.05

    if data2 is not None:
        _, p_norm2 = stats.shapiro(data2[:50])
        assumptions['group2_normal'] = p_norm2 > 0.05

        # Equal variances (Levene)
        _, p_levene = stats.levene(data1, data2)
        assumptions['equal_variances'] = p_levene > 0.05

    return assumptions
```

### 2. Report Effect Sizes
```python
def comprehensive_comparison(group1, group2, protein_name="Protein"):
    """
    Complete comparison following Rosner principles
    """
    from scipy import stats
    import numpy as np

    # Descriptive statistics
    n1, n2 = len(group1), len(group2)
    mean1, mean2 = np.mean(group1), np.mean(group2)
    std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)

    # Statistical test
    t_stat, p_val = stats.ttest_ind(group1, group2)

    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
    cohens_d = (mean1 - mean2) / pooled_std

    # Confidence interval for difference
    se_diff = pooled_std * np.sqrt(1/n1 + 1/n2)
    df = n1 + n2 - 2
    t_crit = stats.t.ppf(0.975, df)
    diff = mean1 - mean2
    ci_lower = diff - t_crit * se_diff
    ci_upper = diff + t_crit * se_diff

    print(f"{protein_name} Analysis")
    print(f"Group 1: {mean1:.3f} ± {std1:.3f} (n={n1})")
    print(f"Group 2: {mean2:.3f} ± {std2:.3f} (n={n2})")
    print(f"Difference: {diff:.3f} (95% CI: {ci_lower:.3f}, {ci_upper:.3f})")
    print(f"Cohen's d: {cohens_d:.3f}")
    print(f"t-statistic: {t_stat:.3f}")
    print(f"p-value: {p_val:.3e}")
    print(f"Interpretation: {interpret_p_value_rosner(p_val)}")
```

---

## Integration with Our Analysis

### For Each Statement Evaluation:

1. **Check Assumptions** using Rosner's diagnostic approaches
2. **Choose Appropriate Test** based on data characteristics
3. **Report Effect Sizes** along with p-values
4. **Interpret Results** using Rosner's p-value framework
5. **Consider Biological Significance** beyond statistical significance

This biostatistical foundation from Rosner ensures our proteomic analyses are statistically sound and properly interpreted in the biological context.
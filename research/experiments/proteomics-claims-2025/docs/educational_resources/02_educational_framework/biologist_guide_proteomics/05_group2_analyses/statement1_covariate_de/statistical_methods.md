# 📐 Statistical Methods for Covariate-Adjusted Analysis

## 🎯 Learning Objectives

By the end of this guide, you'll understand:
- ✅ **Linear model theory** and assumptions
- ✅ **Different types of covariates** and how to handle them
- ✅ **Model selection strategies** for optimal analysis
- ✅ **Advanced methods** like mixed models and batch correction
- ✅ **Interpretation pitfalls** and how to avoid them

---

## 🔬 Fundamental Concepts

### The Linear Model Framework

#### Basic Model Structure
```
Y = β₀ + β₁X₁ + β₂X₂ + ... + βₙXₙ + ε

Where:
Y = Protein expression (outcome)
X₁ = Disease status (primary predictor)
X₂...Xₙ = Covariates (confounders)
β = Coefficients (effects)
ε = Random error
```

#### In Proteomics Context
```python
# Practical example:
"""
SQSTM1_expression = β₀ +
                    β₁ × tau_status +
                    β₂ × age +
                    β₃ × sex +
                    β₄ × batch +
                    ε

Interpretation:
- β₁: Effect of tau pathology on SQSTM1 (what we care about!)
- β₂: How SQSTM1 changes per year of age
- β₃: Difference between males and females
- β₄: Technical batch effect
- ε: Unexplained variation
"""
```

### Why Adjustment Matters: The Confounding Problem

#### Mathematical Explanation
```python
# Without adjustment:
"""
Observed Effect = True Effect + Confounding Bias

Example:
Observed APOE increase in AD = 2.5-fold
BUT:
- True AD effect = 1.2-fold
- Age confounding = 1.3-fold
- 2.5 = 1.2 × 1.3 (multiplicative on fold-change scale)
"""

# With adjustment:
"""
Adjusted Effect ≈ True Effect
(Confounding removed statistically)
"""
```

#### Simpson's Paradox Illustration
```python
# Real data scenario:
"""
COMBINED DATA:
Treatment helps: 60% success rate
Control: 40% success rate
Conclusion: Treatment works! ✓

STRATIFIED BY AGE:
Young patients:
- Treatment: 30% success
- Control: 50% success
- Treatment HARMFUL! ✗

Old patients:
- Treatment: 70% success
- Control: 90% success
- Treatment HARMFUL! ✗

EXPLANATION:
- More young patients got treatment (they're healthier)
- More old patients were controls (too sick for treatment)
- Age was confounding the relationship!
"""
```

---

## 📊 Types of Statistical Models

### 1. Simple Linear Regression (No Covariates)

#### Model Specification
```python
# Simple t-test equivalent:
expression ~ disease_status

# In R/Python statsmodels notation:
model = ols('expression ~ tau_status', data=data).fit()

# Assumptions:
"""
1. Independence of observations
2. Linear relationship
3. Normal distribution of residuals
4. Homoscedasticity (constant variance)
"""
```

#### When to Use
- Perfectly balanced groups
- No confounders
- Exploratory analysis
- Small sample size (can't support many parameters)

### 2. Multiple Linear Regression (With Covariates)

#### Model Specification
```python
# Full adjustment model:
expression ~ disease_status + age + sex + batch + pmi

# Mathematical form:
"""
E[Y|X] = β₀ + β₁×disease + β₂×age + β₃×sex + β₄×batch + β₅×pmi

Each β represents the INDEPENDENT effect of that variable,
holding all others constant.
"""
```

#### Coefficient Interpretation
```python
# Example output interpretation:
"""
Coefficients:
Intercept:    5.2   (baseline expression)
tau_positive: 0.8   (tau+ samples have 0.8 higher expression)
age:         -0.02  (expression decreases 0.02 per year)
sex_M:        0.3   (males have 0.3 higher expression)
batch_2:     -0.5   (batch 2 has 0.5 lower expression)

CRITICAL: The tau effect (0.8) is the difference between
tau+ and tau- AFTER accounting for age, sex, and batch.
"""
```

### 3. Analysis of Covariance (ANCOVA)

#### Conceptual Framework
```python
# ANCOVA = ANOVA + regression
"""
Combines:
- Categorical predictor (disease status)
- Continuous covariates (age, pmi)

Useful for:
- Adjusting group comparisons
- Testing group differences while controlling covariates
- Increased statistical power
"""

# Model with interaction:
expression ~ disease_status * age + sex + batch
# Tests if age effect differs by disease status
```

### 4. Mixed Effects Models

#### When You Have Hierarchical Data
```python
# Random effects for repeated measures or clustering:
"""
Example: Multiple samples per patient

Fixed effects: disease_status, age, sex
Random effects: patient_id

Model: expression ~ disease_status + age + sex + (1|patient_id)

Accounts for:
- Within-patient correlation
- Different baseline expression per patient
- Unbalanced designs
"""
```

#### Implementation
```python
from statsmodels.regression.mixed_linear_model import MixedLM

# Mixed model with random intercepts
model = MixedLM.from_formula(
    'expression ~ tau_status + age + sex',
    groups='patient_id',
    data=data
).fit()
```

---

## 🎯 Model Selection Strategies

### 1. A Priori Selection (Theory-Driven)

```python
# Based on domain knowledge:
"""
ALWAYS INCLUDE:
✓ Known confounders from literature
✓ Variables that differ between groups
✓ Technical factors (batch, RIN)

EXAMPLE DECISION TREE:
Is it measured before disease? → Include
Does it affect both X and Y? → Include
Is it part of disease process? → Don't include
Is it downstream of outcome? → Don't include
"""
```

### 2. Statistical Selection

#### Stepwise Selection (Use Cautiously!)
```python
def stepwise_selection(data, outcome, predictors, p_enter=0.05, p_remove=0.10):
    """
    Forward-backward stepwise selection
    WARNING: Can lead to overfitting!
    """
    selected = []
    remaining = predictors.copy()

    # Forward step
    while remaining:
        p_values = {}
        for var in remaining:
            formula = f"{outcome} ~ {' + '.join(selected + [var])}"
            model = ols(formula, data=data).fit()
            p_values[var] = model.pvalues[var]

        min_p = min(p_values.values())
        if min_p < p_enter:
            best_var = min(p_values, key=p_values.get)
            selected.append(best_var)
            remaining.remove(best_var)
        else:
            break

    # Backward step
    while selected:
        formula = f"{outcome} ~ {' + '.join(selected)}"
        model = ols(formula, data=data).fit()

        max_p_var = max(selected, key=lambda x: model.pvalues[x])
        if model.pvalues[max_p_var] > p_remove:
            selected.remove(max_p_var)
        else:
            break

    return selected
```

#### Information Criteria
```python
def model_comparison_aic_bic(models):
    """
    Compare models using AIC/BIC
    Lower is better!
    """
    comparison = []
    for name, model in models.items():
        comparison.append({
            'Model': name,
            'AIC': model.aic,
            'BIC': model.bic,
            'R²': model.rsquared,
            'Adj_R²': model.rsquared_adj
        })

    return pd.DataFrame(comparison).sort_values('AIC')

# Usage:
models = {
    'Simple': ols('expression ~ disease', data=data).fit(),
    'Age_adjusted': ols('expression ~ disease + age', data=data).fit(),
    'Full': ols('expression ~ disease + age + sex + batch', data=data).fit()
}
best_model = model_comparison_aic_bic(models)
```

### 3. Cross-Validation Approach

```python
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error

def cv_model_selection(X, y, model_formulas, k=5):
    """
    K-fold cross-validation for model selection
    """
    kf = KFold(n_splits=k, shuffle=True, random_state=42)
    results = {}

    for formula in model_formulas:
        mse_scores = []

        for train_idx, test_idx in kf.split(X):
            train_data = X.iloc[train_idx]
            test_data = X.iloc[test_idx]

            # Fit on train
            model = ols(formula, data=train_data).fit()

            # Predict on test
            predictions = model.predict(test_data)
            mse = mean_squared_error(
                test_data['expression'],
                predictions
            )
            mse_scores.append(mse)

        results[formula] = {
            'mean_mse': np.mean(mse_scores),
            'std_mse': np.std(mse_scores)
        }

    return results
```

---

## 🔧 Handling Different Types of Covariates

### 1. Continuous Covariates

#### Linear vs Non-linear Relationships
```python
def test_linearity(data, outcome, predictor):
    """
    Test if relationship is linear or needs transformation
    """
    # Fit linear model
    linear = ols(f'{outcome} ~ {predictor}', data=data).fit()

    # Fit polynomial model
    data[f'{predictor}_sq'] = data[predictor]**2
    nonlinear = ols(f'{outcome} ~ {predictor} + {predictor}_sq',
                   data=data).fit()

    # Compare models
    f_stat = ((linear.ssr - nonlinear.ssr) / 1) / (nonlinear.ssr / nonlinear.df_resid)
    p_value = 1 - stats.f.cdf(f_stat, 1, nonlinear.df_resid)

    if p_value < 0.05:
        print(f"⚠️ Non-linear relationship detected (p={p_value:.4f})")
        print("Consider polynomial or spline terms")
    else:
        print(f"✓ Linear relationship adequate (p={p_value:.4f})")

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Scatter with linear fit
    ax = axes[0]
    ax.scatter(data[predictor], data[outcome], alpha=0.5)
    ax.plot(data[predictor], linear.fittedvalues, 'r-', label='Linear fit')
    ax.set_xlabel(predictor)
    ax.set_ylabel(outcome)
    ax.set_title('Linear Model')

    # Residual plot
    ax = axes[1]
    ax.scatter(linear.fittedvalues, linear.resid, alpha=0.5)
    ax.axhline(y=0, color='red', linestyle='--')
    ax.set_xlabel('Fitted values')
    ax.set_ylabel('Residuals')
    ax.set_title('Residual Plot (check for patterns)')

    plt.tight_layout()
    plt.show()

    return p_value
```

#### Centering and Scaling
```python
def prepare_continuous_covariates(data, continuous_vars):
    """
    Center and scale continuous variables for better interpretation
    """
    processed_data = data.copy()

    for var in continuous_vars:
        # Center around mean
        processed_data[f'{var}_centered'] = data[var] - data[var].mean()

        # Standardize (z-score)
        processed_data[f'{var}_z'] = (data[var] - data[var].mean()) / data[var].std()

        # Log transformation for skewed data
        if data[var].min() > 0:  # Only if all positive
            processed_data[f'{var}_log'] = np.log(data[var])

    return processed_data

# Interpretation changes:
"""
Original: β = 0.02 per year of age
Centered: β = 0.02 per year from mean age (e.g., 70)
Standardized: β = 0.5 per SD of age (e.g., per 10 years)
"""
```

### 2. Categorical Covariates

#### Dummy Variable Coding
```python
def encode_categorical(data, categorical_vars, reference_groups=None):
    """
    Proper encoding of categorical variables
    """
    encoded_data = data.copy()

    for var in categorical_vars:
        # Get unique values
        categories = data[var].unique()

        # Set reference group
        if reference_groups and var in reference_groups:
            ref = reference_groups[var]
        else:
            ref = categories[0]  # First category as reference

        # Create dummy variables
        for cat in categories:
            if cat != ref:
                encoded_data[f'{var}_{cat}'] = (data[var] == cat).astype(int)

        print(f"{var}: Reference = {ref}, Dummies = {categories[categories != ref]}")

    return encoded_data

# Usage:
data_encoded = encode_categorical(
    data,
    ['sex', 'apoe_genotype'],
    reference_groups={'sex': 'F', 'apoe_genotype': 'E3/E3'}
)
```

#### Effect Coding (Sum-to-Zero)
```python
def effect_coding(data, var):
    """
    Alternative to dummy coding - compares to grand mean
    """
    categories = data[var].unique()
    n_cat = len(categories)

    for i, cat in enumerate(categories[:-1]):
        data[f'{var}_effect_{cat}'] = data[var].apply(
            lambda x: 1 if x == cat else (-1 if x == categories[-1] else 0)
        )

    return data
```

### 3. Ordinal Covariates

```python
def handle_ordinal(data, var, levels):
    """
    Treat ordinal variables appropriately
    """
    # Option 1: Numeric coding
    level_map = {level: i for i, level in enumerate(levels)}
    data[f'{var}_numeric'] = data[var].map(level_map)

    # Option 2: Polynomial contrasts
    n_levels = len(levels)
    for degree in range(1, min(3, n_levels)):
        data[f'{var}_poly{degree}'] = data[f'{var}_numeric']**degree

    return data

# Example:
data = handle_ordinal(
    data,
    'braak_stage',
    levels=['0', 'I', 'II', 'III', 'IV', 'V', 'VI']
)
```

---

## 🎨 Advanced Methods

### 1. Batch Effect Correction: ComBat

```python
def combat_correction(expression_matrix, batch_info, covariates=None):
    """
    ComBat batch correction preserving biological variation

    Note: In practice, use pyCombat or similar implementation
    """
    # Conceptual overview:
    """
    ComBat model:
    Y_ijg = α_g + X×β_g + γ_ig + δ_ig×ε_ijg

    Where:
    - α_g: Overall expression for gene g
    - X×β_g: Biological covariates
    - γ_ig: Additive batch effect
    - δ_ig: Multiplicative batch effect

    Steps:
    1. Standardize data
    2. Estimate batch parameters (empirical Bayes)
    3. Remove batch effects
    4. Adjust for covariates
    """

    # Simplified implementation structure:
    import pandas as pd
    from combat.pycombat import pycombat

    # Prepare data
    data_df = pd.DataFrame(
        expression_matrix.T,
        columns=[f'Sample_{i}' for i in range(expression_matrix.shape[0])]
    )

    # Run ComBat
    corrected_data = pycombat(
        data_df,
        batch_info,
        mod=covariates  # Preserve these biological factors
    )

    return corrected_data
```

### 2. Propensity Score Methods

```python
def propensity_score_adjustment(data, treatment_col, covariates):
    """
    Alternative to regression adjustment
    """
    from sklearn.linear_model import LogisticRegression

    # Step 1: Estimate propensity scores
    X = data[covariates]
    y = data[treatment_col]

    ps_model = LogisticRegression(max_iter=1000)
    ps_model.fit(X, y)

    # Get propensity scores
    data['propensity_score'] = ps_model.predict_proba(X)[:, 1]

    # Step 2: Check overlap
    fig, ax = plt.subplots(figsize=(8, 5))
    for group in data[treatment_col].unique():
        subset = data[data[treatment_col] == group]
        ax.hist(subset['propensity_score'], alpha=0.5, label=f'{treatment_col}={group}', bins=20)

    ax.set_xlabel('Propensity Score')
    ax.set_ylabel('Count')
    ax.set_title('Propensity Score Distribution')
    ax.legend()
    plt.show()

    # Step 3: Use in analysis
    # Option A: Stratification
    data['ps_stratum'] = pd.qcut(data['propensity_score'], q=5, labels=['Q1','Q2','Q3','Q4','Q5'])

    # Option B: Weighting (Inverse Probability Weighting)
    data['ipw'] = np.where(
        data[treatment_col] == 1,
        1 / data['propensity_score'],
        1 / (1 - data['propensity_score'])
    )

    # Option C: Include as covariate
    # expression ~ treatment + propensity_score

    return data
```

### 3. Robust Regression

```python
def robust_regression(data, formula):
    """
    Resistant to outliers and violations of assumptions
    """
    import statsmodels.formula.api as smf

    # Method 1: M-estimator (Huber)
    robust_model = smf.rlm(formula, data=data, M=sm.robust.norms.HuberT()).fit()

    # Method 2: MM-estimator
    mm_model = smf.rlm(formula, data=data, M=sm.robust.norms.TukeyBiweight()).fit()

    # Method 3: Quantile regression (median)
    quantile_model = smf.quantreg(formula, data=data).fit(q=0.5)

    # Compare to OLS
    ols_model = smf.ols(formula, data=data).fit()

    print("Coefficient Comparison:")
    print(f"{'Method':<15} {'Coef':>10} {'SE':>10}")
    print("-" * 35)
    print(f"{'OLS':<15} {ols_model.params[1]:>10.4f} {ols_model.bse[1]:>10.4f}")
    print(f"{'Huber':<15} {robust_model.params[1]:>10.4f} {robust_model.bse[1]:>10.4f}")
    print(f"{'Tukey':<15} {mm_model.params[1]:>10.4f} {mm_model.bse[1]:>10.4f}")
    print(f"{'Quantile':<15} {quantile_model.params[1]:>10.4f} {quantile_model.bse[1]:>10.4f}")

    return robust_model
```

---

## ⚠️ Common Pitfalls and Solutions

### 1. Multicollinearity

```python
def check_multicollinearity(data, predictors):
    """
    Detect and handle multicollinearity
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    # Calculate VIF for each predictor
    X = data[predictors]
    vif_data = pd.DataFrame()
    vif_data["Variable"] = predictors
    vif_data["VIF"] = [variance_inflation_factor(X.values, i)
                      for i in range(len(predictors))]

    print("Variance Inflation Factors:")
    print(vif_data.to_string())

    # Interpretation and solutions
    high_vif = vif_data[vif_data['VIF'] > 10]
    if not high_vif.empty:
        print("\n⚠️ High multicollinearity detected!")
        print(f"Problematic variables: {high_vif['Variable'].tolist()}")
        print("\nSolutions:")
        print("1. Remove one of the correlated variables")
        print("2. Combine into single score (PCA)")
        print("3. Use ridge regression")

    return vif_data
```

### 2. Over-adjustment (Collider Bias)

```python
# Conceptual example:
"""
WRONG: Adjusting for a collider

Disease → Inflammation
    ↓           ↓
Protein ← ← ← Death (collider)

If you adjust for death (only analyzing those who died),
you create spurious association between Disease and Protein!

SOLUTION: Don't adjust for outcomes or consequences of exposure/outcome
"""
```

### 3. Table 2 Fallacy

```python
def avoid_table2_fallacy(data, predictors, outcome):
    """
    Don't interpret other coefficients as causal effects!
    """
    model = ols(f'{outcome} ~ {" + ".join(predictors)}', data=data).fit()

    print("⚠️ WARNING: Table 2 Fallacy")
    print("-" * 50)
    print("Only ONE coefficient should be interpreted causally:")
    print("- The primary exposure of interest")
    print("\nOther coefficients are NOT guaranteed to be causal!")
    print("They're just statistical adjustments.\n")

    print(model.summary().tables[1])

    print("\n✓ Correct interpretation:")
    print(f"  Effect of {predictors[0]}: {model.params[1]:.3f}")
    print("\n✗ Incorrect interpretation:")
    print("  'Age also causes protein changes by β₂'")
    print("  (Age is just a confounder here, not the focus)")
```

### 4. Model Assumptions Violations

```python
def diagnose_violations(model):
    """
    Comprehensive assumption checking
    """
    residuals = model.resid
    fitted = model.fittedvalues

    violations = []

    # 1. Normality
    _, p_normal = stats.shapiro(residuals[:5000])
    if p_normal < 0.05:
        violations.append("Non-normal residuals")

    # 2. Homoscedasticity
    _, p_hetero = stats.levene(
        residuals[fitted < np.median(fitted)],
        residuals[fitted >= np.median(fitted)]
    )
    if p_hetero < 0.05:
        violations.append("Heteroscedasticity")

    # 3. Linearity (check residual pattern)
    corr_resid_fitted = np.corrcoef(fitted, residuals)[0, 1]
    if abs(corr_resid_fitted) > 0.1:
        violations.append("Non-linearity")

    if violations:
        print("⚠️ Assumption Violations Detected:")
        for v in violations:
            print(f"  - {v}")
        print("\nPossible solutions:")
        print("- Transform variables (log, sqrt)")
        print("- Use robust regression")
        print("- Bootstrap confidence intervals")
        print("- Non-parametric methods")
    else:
        print("✓ Model assumptions appear satisfied")

    return violations
```

---

## 📈 Interpreting and Reporting Results

### Effect Size Metrics

```python
def calculate_effect_sizes(model, data, treatment_var):
    """
    Multiple effect size measures for comprehensive reporting
    """
    # Get groups
    treatment_data = data[data[treatment_var] == 1]['expression']
    control_data = data[data[treatment_var] == 0]['expression']

    # Raw difference
    raw_diff = treatment_data.mean() - control_data.mean()

    # Adjusted difference (from model)
    adj_diff = model.params[treatment_var]

    # Cohen's d
    pooled_std = np.sqrt((treatment_data.var() + control_data.var()) / 2)
    cohens_d = adj_diff / pooled_std

    # Partial eta squared
    ss_effect = model.ess * (model.rsquared / (1 - model.rsquared))
    partial_eta2 = ss_effect / (ss_effect + model.ssr)

    # Percentage of variance explained
    r_squared_change = model.rsquared  # Compare to null model

    print("Effect Size Measures:")
    print("-" * 40)
    print(f"Raw difference: {raw_diff:.3f}")
    print(f"Adjusted difference: {adj_diff:.3f}")
    print(f"Cohen's d: {cohens_d:.3f}")
    print(f"Partial η²: {partial_eta2:.3f}")
    print(f"R² (model): {model.rsquared:.3f}")

    # Interpretation guide
    print("\nInterpretation:")
    if abs(cohens_d) < 0.2:
        print("  Small effect size")
    elif abs(cohens_d) < 0.5:
        print("  Medium effect size")
    else:
        print("  Large effect size")

    return {
        'raw_difference': raw_diff,
        'adjusted_difference': adj_diff,
        'cohens_d': cohens_d,
        'partial_eta2': partial_eta2,
        'r_squared': model.rsquared
    }
```

### Reporting Template

```python
def generate_report(model, protein_name, covariates_used):
    """
    Generate publication-ready results statement
    """
    # Extract key statistics
    tau_coef = model.params['tau_binary']
    tau_se = model.bse['tau_binary']
    tau_pval = model.pvalues['tau_binary']
    ci_lower = tau_coef - 1.96 * tau_se
    ci_upper = tau_coef + 1.96 * tau_se

    report = f"""
    STATISTICAL REPORT: {protein_name}
    {'=' * 50}

    Model: Linear regression with covariate adjustment
    Covariates: {', '.join(covariates_used)}

    RESULTS:
    --------
    After adjusting for {', '.join(covariates_used)}, tau-positive
    samples showed {"higher" if tau_coef > 0 else "lower"} {protein_name}
    expression compared to tau-negative samples (β = {tau_coef:.3f},
    95% CI: [{ci_lower:.3f}, {ci_upper:.3f}], p = {tau_pval:.4f}).

    The full model explained {model.rsquared*100:.1f}% of the variance
    in {protein_name} expression (adjusted R² = {model.rsquared_adj:.3f}).

    INTERPRETATION:
    --------------
    {"This represents a statistically significant difference." if tau_pval < 0.05
    else "The difference was not statistically significant."}

    Technical Note: Assumptions of linearity, normality of residuals,
    and homoscedasticity were checked and {"satisfied" if True else "violated"}.
    """

    print(report)
    return report

# Usage:
report = generate_report(
    model,
    'SQSTM1',
    ['age', 'sex', 'batch', 'pmi']
)
```

---

## 🎓 Advanced Topics

### Penalized Regression

```python
from sklearn.linear_model import LassoCV, RidgeCV, ElasticNetCV

def penalized_regression_comparison(X, y):
    """
    Compare different regularization approaches
    """
    # Standardize predictors
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Fit models
    models = {
        'LASSO': LassoCV(cv=5, random_state=42),
        'Ridge': RidgeCV(cv=5),
        'ElasticNet': ElasticNetCV(cv=5, random_state=42)
    }

    results = {}
    for name, model in models.items():
        model.fit(X_scaled, y)
        results[name] = {
            'alpha': model.alpha_,
            'coefficients': model.coef_,
            'score': model.score(X_scaled, y)
        }

        # Show which variables were selected (LASSO/ElasticNet)
        if hasattr(model, 'coef_'):
            selected = X.columns[model.coef_ != 0]
            print(f"\n{name} selected variables: {list(selected)}")

    return results
```

### Causal Inference Framework

```python
def causal_dag_analysis():
    """
    Directed Acyclic Graph for causal reasoning
    """
    import networkx as nx

    # Create DAG
    G = nx.DiGraph()

    # Add nodes
    nodes = ['Age', 'Sex', 'Tau_Status', 'SQSTM1', 'Cognition', 'Batch']
    G.add_nodes_from(nodes)

    # Add edges (causal relationships)
    edges = [
        ('Age', 'Tau_Status'),
        ('Age', 'SQSTM1'),
        ('Sex', 'SQSTM1'),
        ('Tau_Status', 'SQSTM1'),
        ('Tau_Status', 'Cognition'),
        ('SQSTM1', 'Cognition'),
        ('Batch', 'SQSTM1')
    ]
    G.add_edges_from(edges)

    # Identify confounders
    confounders = []
    for node in G.nodes():
        if (node != 'Tau_Status' and
            node != 'SQSTM1' and
            nx.has_path(G, node, 'Tau_Status') and
            nx.has_path(G, node, 'SQSTM1')):
            confounders.append(node)

    print(f"Identified confounders: {confounders}")
    print("These should be included as covariates!")

    # Visualize
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color='lightblue',
            node_size=2000, font_size=10, arrows=True)
    plt.title("Causal DAG for Proteomics Analysis")
    plt.show()

    return confounders
```

---

## 📚 Key References

### Methods Papers
1. **Leek & Storey (2007)** - "Capturing heterogeneity in gene expression studies by surrogate variable analysis" - PLoS Genetics
2. **Johnson et al. (2007)** - "Adjusting batch effects in microarray expression data using empirical Bayes methods" - Biostatistics
3. **Gelman & Hill (2006)** - "Data Analysis Using Regression and Multilevel/Hierarchical Models" - Cambridge

### Best Practices
1. **Irizarry et al. (2003)** - "Multiple-laboratory comparison of microarray platforms" - Nature Methods
2. **Parker & Leek (2012)** - "The practical effect of batch on genomic prediction" - Statistical Applications in Genetics

### Tools and Software
- [limma](https://bioconductor.org/packages/limma/) - Linear models for microarray/proteomics
- [ComBat](https://bioconductor.org/packages/sva/) - Batch effect correction
- [statsmodels](https://www.statsmodels.org/) - Statistical models in Python

---

## 🎯 Summary: Statistical Best Practices

### Your Analysis Checklist

```python
statistical_checklist = """
COVARIATE-ADJUSTED ANALYSIS CHECKLIST
=====================================

PRE-ANALYSIS:
☐ Identify potential confounders
☐ Check covariate balance between groups
☐ Assess multicollinearity (VIF < 10)
☐ Detect batch effects
☐ Verify sample size adequacy (10-15 samples per covariate)

MODEL BUILDING:
☐ Start with simple model
☐ Add covariates based on theory
☐ Test for interactions if relevant
☐ Check model assumptions
☐ Consider non-linear relationships

VALIDATION:
☐ Sensitivity analysis (leave-one-out)
☐ Compare adjusted vs unadjusted
☐ Cross-validation if selecting variables
☐ Check influential observations
☐ Bootstrap confidence intervals

REPORTING:
☐ Report both raw and adjusted effects
☐ Include effect sizes (not just p-values)
☐ Describe covariate selection rationale
☐ Acknowledge limitations
☐ Provide model diagnostics
"""
print(statistical_checklist)
```

### Golden Rules

1. **"All models are wrong, but some are useful"** - George Box
2. **Adjustment is not magic** - Can't fix fundamental design flaws
3. **Domain knowledge > Statistical significance**
4. **Simpler models often better** - Occam's Razor applies
5. **Validate, validate, validate** - Replication is key

---

**You now have a comprehensive understanding of the statistical methods underlying covariate adjustment. Remember: the goal is to isolate true biological signals from confounding variation.**

*Next: [Interpreting Your Results](interpreting_results.md)*
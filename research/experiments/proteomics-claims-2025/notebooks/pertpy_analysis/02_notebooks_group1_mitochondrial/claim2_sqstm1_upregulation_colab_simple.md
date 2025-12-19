# SQSTM1 Analysis - Simplified PertPy Version
## Testing: "SQSTM1 is massively upregulated (log2FC = 3.413)"

**🚀 Ultra-simple**: Find SQSTM1, test expression, get answer!

---

## Setup & Data Loading

```python
# 🔧 Setup
import os
IN_COLAB = 'google.colab' in str(get_ipython())

if IN_COLAB:
    !pip install -q pertpy scanpy pandas numpy matplotlib seaborn scipy scikit-learn
    from google.colab import files
    print("📁 Upload your pool_processed_v2.h5ad file:")
    uploaded = files.upload()
    data_file = list(uploaded.keys())[0]
else:
    data_file = "pool_processed_v2.h5ad"

import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')

print("✅ Setup complete!")
```

## Load Data & Find SQSTM1

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

# 🔍 Find SQSTM1 (p62)
protein_names = list(adata.var_names)
sqstm1_patterns = ['SQSTM1', 'P62', 'SEQUESTOSOME']

sqstm1_protein = None
for pattern in sqstm1_patterns:
    matches = [p for p in protein_names if pattern.upper() in p.upper()]
    if matches:
        sqstm1_protein = matches[0]
        break

if sqstm1_protein:
    print(f"✅ Found SQSTM1: {sqstm1_protein}")
    sqstm1_idx = protein_names.index(sqstm1_protein)
else:
    print("❌ SQSTM1 not found!")
    print("Available autophagy proteins:", [p for p in protein_names if any(x in p.upper() for x in ['NBR1', 'OPTN', 'LC3', 'CALCOCO'])][:5])
```

## Analyze SQSTM1 Expression

```python
if sqstm1_protein:
    # 📊 Extract SQSTM1 expression
    sqstm1_expr = adata.X[:, sqstm1_idx]

    # Split by tau status
    tau_pos_expr = sqstm1_expr[adata.obs['tau_positive'] == 1]
    tau_neg_expr = sqstm1_expr[adata.obs['tau_positive'] == 0]

    # 🧮 Calculate statistics
    mean_pos = np.mean(tau_pos_expr)
    mean_neg = np.mean(tau_neg_expr)
    log2fc_observed = mean_pos - mean_neg

    # Statistical tests
    t_stat, p_value = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

    # Results
    print("📈 SQSTM1 ANALYSIS RESULTS")
    print("="*40)
    print(f"Tau+ mean: {mean_pos:.3f}")
    print(f"Tau- mean: {mean_neg:.3f}")
    print(f"Observed log2FC: {log2fc_observed:.3f}")
    print(f"Linear fold change: {2**log2fc_observed:.2f}x")
    print(f"P-value: {p_value:.2e}")
    print(f"Samples: Tau+ n={len(tau_pos_expr)}, Tau- n={len(tau_neg_expr)}")
else:
    print("⚠️ Cannot analyze - SQSTM1 not found")
    log2fc_observed = 0
    p_value = 1
```

## Visualize Results

```python
if sqstm1_protein:
    # 📊 Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Box plot
    data_plot = [tau_neg_expr, tau_pos_expr]
    box_plot = ax1.boxplot(data_plot, labels=['Tau-', 'Tau+'], patch_artist=True)
    box_plot['boxes'][0].set_facecolor('lightblue')
    box_plot['boxes'][1].set_facecolor('lightcoral')
    ax1.set_ylabel('SQSTM1 Expression (log2)')
    ax1.set_title(f'SQSTM1 Expression by Tau Status\\nlog2FC = {log2fc_observed:.2f}, p = {p_value:.2e}')
    ax1.grid(True, alpha=0.3)

    # Comparison bar plot
    claimed_fc = 3.413
    bars = ax2.bar(['Claimed', 'Observed'], [claimed_fc, log2fc_observed],
                   color=['gray', 'green'], alpha=0.7)
    ax2.set_ylabel('Log2 Fold Change')
    ax2.set_title('Fold Change Comparison')
    ax2.grid(True, alpha=0.3, axis='y')

    # Add value labels
    for bar, val in zip(bars, [claimed_fc, log2fc_observed]):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.1,
                f'{val:.2f}', ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.show()

    # 📊 Distribution plot
    plt.figure(figsize=(8, 5))
    plt.hist(tau_neg_expr, alpha=0.6, label='Tau-', bins=20, color='blue')
    plt.hist(tau_pos_expr, alpha=0.6, label='Tau+', bins=20, color='red')
    plt.xlabel('SQSTM1 Expression (log2)')
    plt.ylabel('Number of Samples')
    plt.title('SQSTM1 Expression Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
```

## Pseudotime Analysis (if available)

```python
# 📈 Check for pseudotime data
has_pseudotime = 'pseudotime' in adata.obs.columns or 'Pseudotime' in adata.obs.columns

if has_pseudotime and sqstm1_protein:
    # Get pseudotime column
    if 'pseudotime' in adata.obs.columns:
        pseudotime = adata.obs['pseudotime']
    else:
        pseudotime = adata.obs['Pseudotime']

    # Remove NaN values
    valid_mask = ~(np.isnan(sqstm1_expr) | np.isnan(pseudotime))

    if valid_mask.sum() > 10:
        # Linear regression
        X = pseudotime[valid_mask].values.reshape(-1, 1)
        y = sqstm1_expr[valid_mask]

        reg = LinearRegression().fit(X, y)
        beta = reg.coef_[0]
        r_squared = reg.score(X, y)

        # Correlation
        corr, corr_p = stats.pearsonr(pseudotime[valid_mask], y)

        print("\\n⏰ PSEUDOTIME ANALYSIS")
        print("="*40)
        print(f"Beta coefficient: {beta:.3f}")
        print(f"R-squared: {r_squared:.3f}")
        print(f"Correlation: {corr:.3f} (p = {corr_p:.2e})")
        print(f"Claimed beta: 4.951")
        print(f"Difference: {abs(beta - 4.951):.3f}")

        # Quick plot
        plt.figure(figsize=(8, 5))
        colors = ['red' if x == 1 else 'blue' for x in adata.obs['tau_positive'][valid_mask]]
        plt.scatter(pseudotime[valid_mask], y, c=colors, alpha=0.6)

        # Regression line
        x_line = np.linspace(pseudotime[valid_mask].min(), pseudotime[valid_mask].max(), 100)
        y_line = reg.predict(x_line.reshape(-1, 1))
        plt.plot(x_line, y_line, 'k--', linewidth=2, label=f'β = {beta:.3f}')

        plt.xlabel('Pseudotime (Disease Progression)')
        plt.ylabel('SQSTM1 Expression')
        plt.title('SQSTM1 vs Disease Progression')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
    else:
        print("\\n⚠️ Insufficient pseudotime data")
        beta = None
else:
    print("\\n⚠️ No pseudotime data available")
    beta = None
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*50)
print("🔍 CLAIM EVALUATION")
print("="*50)
print("Claim: SQSTM1 massively upregulated (log2FC = 3.413)")
print()

if sqstm1_protein:
    claimed_fc = 3.413
    fc_diff_pct = abs(log2fc_observed - claimed_fc) / claimed_fc * 100

    # Evaluate fold change
    print(f"📊 FOLD CHANGE ANALYSIS:")
    print(f"Claimed: {claimed_fc:.3f} ({2**claimed_fc:.1f}x)")
    print(f"Observed: {log2fc_observed:.3f} ({2**log2fc_observed:.1f}x)")
    print(f"Difference: {fc_diff_pct:.1f}%")

    # Evaluate significance
    print(f"\\n📈 SIGNIFICANCE:")
    print(f"P-value: {p_value:.2e}")
    print(f"Significant: {'Yes' if p_value < 0.05 else 'No'}")

    # Overall verdict
    is_upregulated = log2fc_observed > 0.5
    is_significant = p_value < 0.05

    if is_upregulated and is_significant:
        if fc_diff_pct < 30:
            verdict = "✅ SUPPORTED"
            explanation = "SQSTM1 is significantly upregulated as claimed"
        else:
            verdict = "⚠️ PARTIALLY SUPPORTED"
            explanation = f"SQSTM1 upregulated but magnitude differs ({fc_diff_pct:.1f}% difference)"
    elif is_upregulated:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = "SQSTM1 upregulated but not statistically significant"
    else:
        verdict = "❌ REFUTED"
        explanation = "SQSTM1 not significantly upregulated"

    print(f"\\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Pseudotime verdict
    if beta is not None:
        if beta > 2:
            time_verdict = "✅ Increases with disease progression"
        elif beta > 0:
            time_verdict = "⚠️ Modest increase with progression"
        else:
            time_verdict = "❌ No clear progression relationship"
        print(f"⏰ Temporal pattern: {time_verdict}")

    # Biological interpretation
    print(f"\\n🧬 BIOLOGICAL MEANING:")
    if is_upregulated:
        print("• SQSTM1/p62 accumulation suggests autophagy dysfunction")
        print("• Impaired protein clearance in tau-positive neurons")
        print("• Supports proteostasis failure hypothesis")
    else:
        print("• No clear autophagy receptor accumulation")
        print("• Protein clearance mechanisms may be functional")

else:
    print("❌ Cannot evaluate - SQSTM1 not found in dataset")
```

## Save Results

```python
# 💾 Save results
if sqstm1_protein:
    results = {
        'protein': sqstm1_protein,
        'claimed_log2FC': 3.413,
        'observed_log2FC': log2fc_observed,
        'p_value': p_value,
        'verdict': verdict if 'verdict' in locals() else 'Unknown',
        'tau_pos_mean': mean_pos,
        'tau_neg_mean': mean_neg,
        'beta_pseudotime': beta if beta is not None else 'N/A'
    }

    results_df = pd.DataFrame([results])

    if IN_COLAB:
        results_df.to_csv('sqstm1_results.csv', index=False)
        files.download('sqstm1_results.csv')
        print("📁 Results downloaded!")
    else:
        results_df.to_csv('sqstm1_analysis.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ SQSTM1 analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified SQSTM1 analysis:
- ✅ **Finds SQSTM1** automatically in your dataset
- ✅ **Tests the claim** with clear statistics
- ✅ **Compares fold changes** (claimed vs observed)
- ✅ **Analyzes temporal patterns** if pseudotime available
- ✅ **Gives clear verdict** with biological interpretation

**Perfect for testing autophagy dysfunction claims!** 🔬
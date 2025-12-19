# UPS Protein Analysis - Simplified PertPy Version
## Testing: "No significant UPS protein alterations across tau-positive versus tau-negative neurons"

**🚀 Ultra-simple Colab version**: Just run the cells and get results!

---

## Setup & Data Loading

```python
# 🔧 Install packages and setup
import os
IN_COLAB = 'google.colab' in str(get_ipython())

if IN_COLAB:
    !pip install -q pertpy scanpy pandas numpy matplotlib seaborn scipy
    from google.colab import files
    print("📁 Upload your pool_processed_v2.h5ad file:")
    uploaded = files.upload()
    data_file = list(uploaded.keys())[0]
else:
    data_file = "pool_processed_v2.h5ad"

# 📚 Import libraries
import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("✅ Setup complete!")
```

## Load Data & Define Proteins

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Standardize tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🧬 UPS Proteins (simplified list of key proteins)
ups_proteins = [
    # Proteasome core
    'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
    'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7',

    # Proteasome regulatory
    'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
    'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD11', 'PSMD12',

    # E1/E2/E3 enzymes
    'UBA1', 'UBB', 'UBC', 'UBE2D1', 'UBE2D3', 'UBE2L3',
    'UBE3A', 'UBE3B', 'PARK2', 'HUWE1',

    # Deubiquitinases
    'UCHL1', 'UCHL5', 'USP7', 'USP14', 'USP5', 'USP9X',

    # UPS regulators
    'VCP', 'SQSTM1', 'UBQLN1', 'UBQLN2', 'NBR1', 'OPTN'
]

print(f"🎯 Analyzing {len(ups_proteins)} key UPS proteins")
```

## Find & Subset UPS Proteins

```python
# 🔍 Find available UPS proteins
protein_names = adata.var_names if 'gene_name' not in adata.var else adata.var['gene_name']
available_ups = [p for p in ups_proteins if p in protein_names]
missing = set(ups_proteins) - set(available_ups)

print(f"Found: {len(available_ups)}/{len(ups_proteins)} UPS proteins ({len(available_ups)/len(ups_proteins)*100:.1f}%)")
if missing and len(missing) <= 10:
    print(f"Missing: {list(missing)}")

# ✂️ Create UPS subset
ups_mask = [p in available_ups for p in protein_names]
adata_ups = adata[:, ups_mask].copy()
print(f"UPS subset: {adata_ups.shape}")
```

## Run Simple Differential Expression

```python
# 🧮 Simple t-test approach (most reliable)
results = []

for i, protein in enumerate(adata_ups.var_names):
    # Get expression data
    expr = adata_ups.X[:, i]
    tau_pos_expr = expr[adata_ups.obs['tau_positive'] == 1]
    tau_neg_expr = expr[adata_ups.obs['tau_positive'] == 0]

    # Calculate statistics
    mean_pos = np.mean(tau_pos_expr)
    mean_neg = np.mean(tau_neg_expr)
    log2fc = mean_pos - mean_neg

    # T-test
    t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

    results.append({
        'protein': protein,
        'log2FC': log2fc,
        'p_value': p_val,
        'mean_tau_pos': mean_pos,
        'mean_tau_neg': mean_neg
    })

# Convert to DataFrame and add FDR correction
df = pd.DataFrame(results)
from statsmodels.stats.multitest import fdrcorrection
df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
df['significant'] = df['p_adjusted'] < 0.05

print(f"🔬 Analysis complete!")
print(f"Significant proteins (FDR < 0.05): {df['significant'].sum()}/{len(df)}")
```

## Results & Visualization

```python
# 📊 Quick summary
n_sig = df['significant'].sum()
pct_sig = n_sig / len(df) * 100

print("📈 RESULTS SUMMARY")
print("="*40)
print(f"UPS proteins tested: {len(df)}")
print(f"Significantly altered: {n_sig} ({pct_sig:.1f}%)")

if n_sig > 0:
    sig_df = df[df['significant']].sort_values('p_adjusted')
    print(f"\\nTop changed proteins:")
    for _, row in sig_df.head(5).iterrows():
        direction = "↑" if row['log2FC'] > 0 else "↓"
        print(f"  {row['protein']:10} {direction} {abs(row['log2FC']):.2f} (p={row['p_adjusted']:.2e})")

# 🌋 Volcano plot
plt.figure(figsize=(8, 6))
colors = ['red' if sig else 'gray' for sig in df['significant']]
plt.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6)
plt.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
plt.axvline(0, color='black', linestyle='-', alpha=0.3)
plt.xlabel('Log2 Fold Change (Tau+ vs Tau-)')
plt.ylabel('-Log10(p-value)')
plt.title(f'UPS Proteins Volcano Plot\\n{n_sig}/{len(df)} significant')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*50)
print("🔍 CLAIM EVALUATION")
print("="*50)
print("Claim: 'No significant UPS protein alterations'")
print()

# Decision logic
if pct_sig < 5:
    verdict = "✅ SUPPORTED"
    explanation = f"Only {pct_sig:.1f}% of UPS proteins altered"
elif pct_sig < 15:
    verdict = "⚠️ PARTIALLY SUPPORTED"
    explanation = f"{pct_sig:.1f}% altered - borderline significant"
elif pct_sig < 30:
    verdict = "❌ PARTIALLY REFUTED"
    explanation = f"{pct_sig:.1f}% of UPS proteins significantly altered"
else:
    verdict = "❌ REFUTED"
    explanation = f"{pct_sig:.1f}% of UPS proteins show major alterations"

print(f"VERDICT: {verdict}")
print(f"Reason: {explanation}")

# 🧬 Biological interpretation
print(f"\\n🧬 BIOLOGICAL IMPACT:")
if n_sig > len(df) * 0.2:  # >20% altered
    print("• Major UPS dysfunction detected")
    print("• Protein degradation likely impaired")
    print("• May contribute to protein aggregation")
else:
    print("• UPS function appears largely preserved")
    print("• Protein degradation capacity maintained")

print(f"\\n📊 Quick Stats:")
print(f"• Proteins analyzed: {len(df)}/{len(ups_proteins)}")
print(f"• Effect sizes: {df['log2FC'].abs().mean():.3f} avg |log2FC|")
print(f"• Largest change: {df.loc[df['log2FC'].abs().idxmax(), 'protein']} ({df['log2FC'].abs().max():.2f})")
```

## Save Results

```python
# 💾 Save results
if IN_COLAB:
    df.to_csv('ups_results.csv', index=False)
    files.download('ups_results.csv')
    print("📁 Results downloaded!")
else:
    df.to_csv('ups_analysis_results.csv', index=False)
    print("📁 Results saved to ups_analysis_results.csv")

print("\\n✅ Analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified analysis:
- ✅ **Tests key UPS proteins** with robust statistics
- ✅ **Clear visualization** with volcano plot
- ✅ **Objective evaluation** of the claim
- ✅ **Biological interpretation** of results
- ✅ **Easy to run** in Google Colab

**Just upload your data file and run all cells!** 🚀
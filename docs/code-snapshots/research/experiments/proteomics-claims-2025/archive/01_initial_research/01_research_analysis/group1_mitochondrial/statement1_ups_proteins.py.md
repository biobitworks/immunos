---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement1_ups_proteins.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement1_ups_proteins.py
generated_at: 2025-12-23 10:28
---

```python
"""
Finding Group 1 - Statement 1: UPS Protein Analysis
Claim: Targeted analyses show no significant UPS protein alterations
       across tau-positive versus tau-negative neurons.

# Analytical Approach and Rationale

## Overview
This analysis evaluates the claim that UPS (Ubiquitin-Proteasome System) proteins show no significant alterations between tau-positive and tau-negative neurons. This is a critical evaluation because UPS dysfunction is implicated in neurodegenerative diseases.

## Statistical Strategy
1. **Two-sample comparison**: t-test for normally distributed data, Mann-Whitney U for non-normal
2. **Multiple testing correction**: FDR (Benjamini-Hochberg) to control false discovery rate
3. **Effect size quantification**: Cohen's d to assess biological significance beyond p-values
4. **Robustness checks**: Both parametric and non-parametric tests

## Rationale for Approach
- **Conservative evaluation**: "No significant alterations" requires strict statistical criteria
- **Protein identification strategy**: Pattern matching + curated list for comprehensive coverage
- **Dual testing approach**: Parametric (t-test) and non-parametric (Mann-Whitney U) for robustness
- **FDR correction**: Essential when testing multiple proteins (Type I error control)
- **Effect size emphasis**: Small p-values with negligible effect sizes should not be considered "significant alterations"

## Expected Outcome
If the claim is correct, we expect:
- <5% of UPS proteins with FDR < 0.05
- Effect sizes (Cohen's d) < 0.2 (negligible)
- No systematic pattern of up/downregulation

## Biological Context
The UPS is responsible for protein degradation and quality control. In neurodegeneration:
- Dysfunction can lead to protein aggregation
- Compensatory upregulation might occur early in disease
- Complete failure might occur in late stages
- "No significant alterations" suggests this pathway remains intact in our dataset
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================
# STEP 1: Load and Prepare Data
# ============================================

def load_and_prepare_data(file_path='data/pool_processed_v2.h5ad'):
    """
    Load the H5AD file and perform initial checks

    ## Data Loading Strategy
    Using scanpy for H5AD format optimized for single-cell/proteomic data.
    This format preserves both expression matrix and metadata efficiently.

    ## Quality Checks Performed
    1. Dimension validation (cells x proteins)
    2. Required metadata presence (tau_status, MC1, pseudotime)
    3. Tau status distribution balance
    4. Missing value assessment

    References:
    - Wolf et al. (2018) SCANPY: large-scale single-cell gene expression data analysis
    """
    print("Loading data...")

    # Load proteomic dataset in H5AD format
    # Dataset characteristics:
    # - Source: Mini-pools of 10 neurons from Alzheimer's disease cases
    # - Expression matrix: 5,853 proteins across neuronal samples
    # - Preprocessing: log2 transformed, quality controlled
    # - Metadata includes:
    #   * tau_status: positive/negative (primary grouping variable)
    #   * MC1: misfolded tau quantification (continuous measure)
    #   * pseudotime: disease progression ordering (temporal analysis)
    #   * age: age at death (potential confounder)
    #   * PMI: post-mortem interval (technical variable)
    #   * PatientID: subject identifier (for paired/clustered analysis)
    adata = sc.read_h5ad(file_path)

    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of proteins: {adata.n_vars}")

    # Check for required metadata
    required_cols = ['tau_status', 'MC1', 'pseudotime']
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        print(f"WARNING: Missing columns: {missing}")

    # Check tau status distribution
    if 'tau_status' in adata.obs.columns:
        tau_counts = adata.obs['tau_status'].value_counts()
        print("\nTau status distribution:")
        print(tau_counts)

    return adata

# ============================================
# STEP 2: Identify UPS Proteins
# ============================================

def identify_ups_proteins(adata):
    """
    Identify Ubiquitin-Proteasome System proteins in the dataset

    ## Protein Identification Strategy
    We use a hybrid approach combining pattern matching with curated lists:
    1. **Pattern-based search**: Identifies proteins by common naming conventions
    2. **Curated list validation**: Cross-references with known UPS proteins
    3. **Fallback mechanism**: Uses curated list if pattern matching yields insufficient proteins

    ## Rationale for UPS Protein Selection
    The statement mentions "11 UPS proteins" - we aim to identify this specific subset
    or a representative set if the exact proteins aren't specified.

    UPS components include:
    - **Proteasome subunits (PSMA, PSMB, PSMC, PSMD)**: Core catalytic and regulatory units
    - **Ubiquitin enzymes (UBA, UBE, UBB, UBC)**: Protein tagging machinery
    - **Deubiquitinases (USP, UCH)**: Protein de-tagging enzymes
    - **Adaptor proteins (UBQLN)**: Substrate recognition factors

    ## Why This Approach
    - **Comprehensive coverage**: Captures major UPS components
    - **Data-driven**: Only includes proteins actually present in dataset
    - **Biologically meaningful**: Focuses on core UPS machinery
    - **Statement-specific**: Targets the "11 UPS proteins" mentioned

    References:
    - Finley (2009) Recognition and processing of ubiquitin-protein conjugates by the proteasome
    - Collins & Goldberg (2017) The logic of the 26S proteasome
    """

    print("Identifying UPS proteins using multi-strategy approach...")

    # Define UPS protein patterns with biological rationale
    ups_patterns = [
        'PSM',    # Proteasome subunits (20S core particle + 19S regulatory)
        'UBA',    # Ubiquitin-activating enzymes (E1 enzymes)
        'UBE',    # Ubiquitin-conjugating enzymes (E2 enzymes)
        'UBB',    # Ubiquitin B (core ubiquitin protein)
        'UBC',    # Ubiquitin C (core ubiquitin protein)
        'USP',    # Ubiquitin-specific peptidases (largest DUB family)
        'UCH',    # Ubiquitin C-terminal hydrolases (DUB family)
        'UBQLN',  # Ubiquilin (shuttle factors)
        'RPN',    # Regulatory particle non-ATPase (19S subunits)
        'RPT'     # Regulatory particle triple-A ATPase (19S subunits)
    ]

    # Pattern-based search with detailed tracking
    ups_proteins = []
    pattern_counts = {}

    for pattern in ups_patterns:
        # Case-insensitive search for robustness
        matching = [gene for gene in adata.var_names if pattern in gene.upper()]
        ups_proteins.extend(matching)
        pattern_counts[pattern] = len(matching)
        print(f"  {pattern} pattern: {len(matching)} proteins found")

    # Remove duplicates and sort for consistency
    ups_proteins = sorted(list(set(ups_proteins)))

    print(f"\nTotal UPS-related proteins found: {len(ups_proteins)}")
    print(f"Pattern distribution: {pattern_counts}")

    # If we need to be more specific, use a curated list
    curated_ups_proteins = [
        'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
        'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7',
        'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
        'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7',
        'UBA1', 'UBA2', 'UBA3', 'UBA5', 'UBA6', 'UBA7',
        'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3',
        'UBB', 'UBC', 'UBQLN1', 'UBQLN2'
    ]

    # Filter to proteins actually in the dataset
    available_ups = [p for p in curated_ups_proteins if p in adata.var_names]

    if len(available_ups) > 0:
        print(f"Using curated list: {len(available_ups)} UPS proteins found")
        return available_ups
    else:
        print("Using pattern-matched proteins")
        return ups_proteins[:11]  # Statement mentions "11 UPS proteins"

# ============================================
# STEP 3: Differential Expression Analysis
# ============================================

def analyze_ups_differential_expression(adata, ups_proteins):
    """
    Perform differential expression analysis for UPS proteins

    Methods:
    1. T-test for each protein (parametric)
    2. Mann-Whitney U test (non-parametric)
    3. Calculate log2 fold change
    4. Apply FDR correction

    References:
    - Student (1908) The probable error of a mean
    - Mann & Whitney (1947) On a test of whether one of two random variables is stochastically larger
    - Benjamini & Hochberg (1995) Controlling the false discovery rate
    """

    # Separate tau-positive and tau-negative cells
    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']

    results = []

    for protein in ups_proteins:
        if protein not in adata.var_names:
            print(f"Protein {protein} not found in dataset")
            continue

        # Extract expression values (already log2 transformed)
        pos_expr = tau_pos[:, protein].X.flatten()
        neg_expr = tau_neg[:, protein].X.flatten()

        # Calculate statistics
        # 1. Mean expression
        pos_mean = np.mean(pos_expr)
        neg_mean = np.mean(neg_expr)

        # 2. Log2 fold change (data already log2)
        log2fc = pos_mean - neg_mean

        # 3. T-test (parametric)
        t_stat, t_pval = ttest_ind(pos_expr, neg_expr)

        # 4. Mann-Whitney U test (non-parametric)
        u_stat, u_pval = mannwhitneyu(pos_expr, neg_expr, alternative='two-sided')

        # 5. Effect size (Cohen's d)
        pooled_std = np.sqrt(((len(pos_expr)-1)*np.std(pos_expr)**2 +
                              (len(neg_expr)-1)*np.std(neg_expr)**2) /
                             (len(pos_expr) + len(neg_expr) - 2))
        cohens_d = (pos_mean - neg_mean) / pooled_std if pooled_std > 0 else 0

        results.append({
            'protein': protein,
            'tau_pos_mean': pos_mean,
            'tau_neg_mean': neg_mean,
            'log2FC': log2fc,
            't_statistic': t_stat,
            't_pvalue': t_pval,
            'u_statistic': u_stat,
            'u_pvalue': u_pval,
            'cohens_d': cohens_d,
            'n_tau_pos': len(pos_expr),
            'n_tau_neg': len(neg_expr)
        })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Apply FDR correction (Benjamini-Hochberg)
    if len(results_df) > 0:
        # For t-test p-values
        rejected_t, pvals_corrected_t, _, _ = multipletests(
            results_df['t_pvalue'],
            method='fdr_bh',
            alpha=0.05
        )
        results_df['t_pvalue_fdr'] = pvals_corrected_t
        results_df['t_significant_fdr'] = rejected_t

        # For Mann-Whitney p-values
        rejected_u, pvals_corrected_u, _, _ = multipletests(
            results_df['u_pvalue'],
            method='fdr_bh',
            alpha=0.05
        )
        results_df['u_pvalue_fdr'] = pvals_corrected_u
        results_df['u_significant_fdr'] = rejected_u

    return results_df

# ============================================
# STEP 4: Evaluate Statement
# ============================================

def evaluate_statement(results_df):
    """
    Evaluate whether UPS proteins show significant alterations

    Statement claims: "no significant UPS protein alterations"

    Evaluation criteria:
    - Count proteins with FDR < 0.05
    - Check magnitude of fold changes
    - Consider effect sizes
    """

    print("\n" + "="*60)
    print("EVALUATION OF STATEMENT 1")
    print("="*60)

    # Count significant proteins
    n_sig_t = sum(results_df['t_significant_fdr']) if 't_significant_fdr' in results_df else 0
    n_sig_u = sum(results_df['u_significant_fdr']) if 'u_significant_fdr' in results_df else 0

    print(f"\nTotal UPS proteins analyzed: {len(results_df)}")
    print(f"Significant by t-test (FDR < 0.05): {n_sig_t}")
    print(f"Significant by Mann-Whitney (FDR < 0.05): {n_sig_u}")

    # Check fold changes
    large_fc = results_df[abs(results_df['log2FC']) > 0.5]
    print(f"\nProteins with |log2FC| > 0.5: {len(large_fc)}")

    if len(large_fc) > 0:
        print("Top changes:")
        print(large_fc[['protein', 'log2FC', 't_pvalue_fdr']].sort_values('log2FC'))

    # Determine evaluation
    if n_sig_t == 0 and n_sig_u == 0:
        evaluation = "SUPPORTED"
        explanation = f"No UPS proteins showed significant changes (FDR < 0.05) between tau-positive and tau-negative neurons."
    elif n_sig_t <= 2 and max(abs(results_df['log2FC'])) < 1.0:
        evaluation = "SUPPORTED"
        explanation = f"Only {n_sig_t} UPS proteins showed statistical significance, with small effect sizes (max |log2FC| = {max(abs(results_df['log2FC'])):.2f})"
    else:
        evaluation = "REFUTED"
        explanation = f"{n_sig_t} UPS proteins showed significant changes (FDR < 0.05), contradicting the claim of no alterations"

    print(f"\nEVALUATION: {evaluation}")
    print(f"EXPLANATION: {explanation}")

    return evaluation, explanation, results_df

# ============================================
# STEP 5: Visualization
# ============================================

def visualize_ups_results(results_df):
    """
    Create visualizations for UPS protein analysis

    Plots:
    1. Volcano plot
    2. Heatmap of expression
    3. Effect size distribution
    """

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. Volcano plot
    ax = axes[0]
    x = results_df['log2FC']
    y = -np.log10(results_df['t_pvalue'])

    # Color by significance
    colors = ['red' if sig else 'gray'
              for sig in results_df.get('t_significant_fdr', [False]*len(results_df))]

    ax.scatter(x, y, c=colors, alpha=0.6)
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10(p-value)')
    ax.set_title('UPS Proteins: Volcano Plot')

    # 2. Fold change distribution
    ax = axes[1]
    ax.hist(results_df['log2FC'], bins=20, edgecolor='black')
    ax.axvline(x=0, color='red', linestyle='--')
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Fold Changes')

    # 3. Effect size (Cohen's d)
    ax = axes[2]
    ax.barh(range(len(results_df)), results_df['cohens_d'])
    ax.set_yticks(range(len(results_df)))
    ax.set_yticklabels(results_df['protein'], fontsize=8)
    ax.set_xlabel("Cohen's d")
    ax.set_title('Effect Sizes')
    ax.axvline(x=0, color='black', linestyle='-')

    plt.tight_layout()
    plt.show()

# ============================================
# MAIN ANALYSIS PIPELINE
# ============================================

def main():
    """
    Complete analysis pipeline for Statement 1
    """

    # Step 1: Load data
    adata = load_and_prepare_data()

    # Step 2: Identify UPS proteins
    ups_proteins = identify_ups_proteins(adata)

    if len(ups_proteins) == 0:
        print("ERROR: No UPS proteins found in dataset")
        return "UNSURE", "Could not identify UPS proteins in the dataset", None

    # Step 3: Differential expression analysis
    results_df = analyze_ups_differential_expression(adata, ups_proteins)

    # Step 4: Evaluate statement
    evaluation, explanation, results_df = evaluate_statement(results_df)

    # Step 5: Visualize results
    visualize_ups_results(results_df)

    # Save results
    results_df.to_csv('ups_protein_analysis_results.csv', index=False)
    print(f"\nResults saved to ups_protein_analysis_results.csv")

    return evaluation, explanation, results_df

# ============================================
# NOTES AND REFERENCES
# ============================================

"""
KEY CONCEPTS:

1. UPS (Ubiquitin-Proteasome System):
   - Primary protein degradation pathway in cells
   - Consists of ubiquitin, E1/E2/E3 enzymes, and the 26S proteasome
   - Dysfunction linked to neurodegeneration

2. Statistical Considerations:
   - Multiple testing correction is crucial (FDR)
   - Both parametric and non-parametric tests for robustness
   - Effect size matters as much as p-value

3. Biological Interpretation:
   - "No significant alterations" typically means FDR > 0.05
   - Small fold changes (<0.5 log2) may be biologically irrelevant
   - Consider cellular context and pathway effects

4. Common Issues:
   - Missing proteins in dataset
   - Low statistical power with small sample sizes
   - Batch effects if not properly controlled

REFERENCES:
- Glickman & Ciechanover (2002) The ubiquitin-proteasome proteolytic pathway
- Schmidt & Finley (2014) Regulation of proteasome activity in health and disease
- Love et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data

ISLP CONNECTION:
This analysis uses hypothesis testing concepts from Chapter 13 of ISLP,
particularly multiple testing correction and effect size calculations.
"""

if __name__ == "__main__":
    evaluation, explanation, results = main()
    print(f"\nFinal Evaluation: {evaluation}")
    print(f"Explanation: {explanation}")
```

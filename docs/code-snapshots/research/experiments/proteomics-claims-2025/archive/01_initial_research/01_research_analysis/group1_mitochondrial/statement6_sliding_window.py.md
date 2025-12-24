---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement6_sliding_window.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement6_sliding_window.py
generated_at: 2025-12-23 10:28
---

```python
"""
Finding Group 1 - Statement 6: Sliding Window Correlation Analysis
Claim: A running correlation between SQSTM1 and VDAC1 along pseudotime (sliding window n = 20)
       shifts from significantly negative early (mean r = -0.417 for pseudotime < 0.33) to
       positive late (mean r = 0.478 for pseudotime > 0.67), with a strong trend
       (r = 0.851, p = 6.98 × 10^-8)

This is an advanced time-series analysis showing dynamic relationship changes.
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# ============================================
# BIOLOGICAL CONTEXT
# ============================================

"""
SQSTM1-VDAC1 Relationship:
- SQSTM1/p62: Autophagy receptor, accumulates when autophagy fails
- VDAC1: Voltage-dependent anion channel, mitochondrial outer membrane
- Their relationship indicates mitophagy (mitochondrial autophagy) status

Dynamic correlation changes suggest:
- Early: Compensatory mechanism (negative correlation)
- Late: System failure (positive correlation)

References:
- Geisler et al. (2010) PINK1/Parkin-mediated mitophagy
- Narendra et al. (2008) Parkin is recruited selectively to impaired mitochondria
"""

# ============================================
# STEP 1: Sliding Window Implementation
# ============================================

def sliding_window_correlation(adata, protein1='SQSTM1', protein2='VDAC1',
                              window_size=20, step_size=1):
    """
    Calculate correlation in sliding windows along pseudotime

    Parameters:
    -----------
    window_size : int
        Number of cells in each window (n=20 as specified)
    step_size : int
        Step between windows (1 for maximum resolution)

    Theory:
    --------
    Sliding window analysis reveals temporal dynamics that
    static correlations miss. This is similar to:
    - Moving averages in time series (ISLP Ch 8)
    - Local regression (LOESS)
    - Change point detection
    """

    print("="*60)
    print("SLIDING WINDOW CORRELATION ANALYSIS")
    print("="*60)
    print(f"Proteins: {protein1} vs {protein2}")
    print(f"Window size: {window_size}")
    print(f"Step size: {step_size}")

    # Check proteins exist
    if protein1 not in adata.var_names or protein2 not in adata.var_names:
        print(f"ERROR: Proteins not found in dataset")
        return None

    # Sort cells by pseudotime
    print("\nSorting cells by pseudotime...")
    sorted_indices = np.argsort(adata.obs['pseudotime'].values)
    adata_sorted = adata[sorted_indices]

    # Extract expression data
    expr1 = adata_sorted[:, protein1].X.flatten()
    expr2 = adata_sorted[:, protein2].X.flatten()
    pseudotime_sorted = adata_sorted.obs['pseudotime'].values

    # Initialize results
    results = {
        'window_start_idx': [],
        'window_end_idx': [],
        'pseudotime_mean': [],
        'pseudotime_min': [],
        'pseudotime_max': [],
        'correlation': [],
        'p_value': [],
        'n_cells': [],
        'expr1_mean': [],
        'expr2_mean': [],
        'expr1_std': [],
        'expr2_std': []
    }

    # Slide window through data
    n_windows = 0
    for i in range(0, len(adata_sorted) - window_size + 1, step_size):
        window_indices = range(i, i + window_size)

        # Extract window data
        window_expr1 = expr1[window_indices]
        window_expr2 = expr2[window_indices]
        window_pseudotime = pseudotime_sorted[window_indices]

        # Calculate correlation
        if np.std(window_expr1) > 0 and np.std(window_expr2) > 0:
            corr, p_val = pearsonr(window_expr1, window_expr2)

            # Store results
            results['window_start_idx'].append(i)
            results['window_end_idx'].append(i + window_size - 1)
            results['pseudotime_mean'].append(np.mean(window_pseudotime))
            results['pseudotime_min'].append(np.min(window_pseudotime))
            results['pseudotime_max'].append(np.max(window_pseudotime))
            results['correlation'].append(corr)
            results['p_value'].append(p_val)
            results['n_cells'].append(window_size)
            results['expr1_mean'].append(np.mean(window_expr1))
            results['expr2_mean'].append(np.mean(window_expr2))
            results['expr1_std'].append(np.std(window_expr1))
            results['expr2_std'].append(np.std(window_expr2))

            n_windows += 1

    print(f"\nProcessed {n_windows} windows")

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    return results_df

# ============================================
# STEP 2: Phase Analysis
# ============================================

def analyze_phases(results_df):
    """
    Analyze early, middle, and late phases of correlation

    Statement specifies:
    - Early: pseudotime < 0.33, mean r = -0.417
    - Late: pseudotime > 0.67, mean r = 0.478
    """

    print("\n" + "="*60)
    print("PHASE ANALYSIS")
    print("="*60)

    # Define phases
    early_phase = results_df[results_df['pseudotime_mean'] < 0.33]
    middle_phase = results_df[(results_df['pseudotime_mean'] >= 0.33) &
                              (results_df['pseudotime_mean'] <= 0.67)]
    late_phase = results_df[results_df['pseudotime_mean'] > 0.67]

    phase_stats = {}

    # Early phase statistics
    if len(early_phase) > 0:
        early_corrs = early_phase['correlation'].values
        phase_stats['early'] = {
            'mean_correlation': np.mean(early_corrs),
            'std_correlation': np.std(early_corrs),
            'median_correlation': np.median(early_corrs),
            'min_correlation': np.min(early_corrs),
            'max_correlation': np.max(early_corrs),
            'n_windows': len(early_phase),
            'ci_lower': np.percentile(early_corrs, 2.5),
            'ci_upper': np.percentile(early_corrs, 97.5)
        }

        print(f"Early Phase (pseudotime < 0.33):")
        print(f"  Mean correlation: {phase_stats['early']['mean_correlation']:.3f}")
        print(f"  95% CI: [{phase_stats['early']['ci_lower']:.3f}, {phase_stats['early']['ci_upper']:.3f}]")
        print(f"  N windows: {phase_stats['early']['n_windows']}")
        print(f"  Expected: -0.417")

    # Middle phase statistics
    if len(middle_phase) > 0:
        middle_corrs = middle_phase['correlation'].values
        phase_stats['middle'] = {
            'mean_correlation': np.mean(middle_corrs),
            'std_correlation': np.std(middle_corrs),
            'n_windows': len(middle_phase)
        }

        print(f"\nMiddle Phase (0.33 ≤ pseudotime ≤ 0.67):")
        print(f"  Mean correlation: {phase_stats['middle']['mean_correlation']:.3f}")
        print(f"  N windows: {phase_stats['middle']['n_windows']}")

    # Late phase statistics
    if len(late_phase) > 0:
        late_corrs = late_phase['correlation'].values
        phase_stats['late'] = {
            'mean_correlation': np.mean(late_corrs),
            'std_correlation': np.std(late_corrs),
            'median_correlation': np.median(late_corrs),
            'min_correlation': np.min(late_corrs),
            'max_correlation': np.max(late_corrs),
            'n_windows': len(late_phase),
            'ci_lower': np.percentile(late_corrs, 2.5),
            'ci_upper': np.percentile(late_corrs, 97.5)
        }

        print(f"\nLate Phase (pseudotime > 0.67):")
        print(f"  Mean correlation: {phase_stats['late']['mean_correlation']:.3f}")
        print(f"  95% CI: [{phase_stats['late']['ci_lower']:.3f}, {phase_stats['late']['ci_upper']:.3f}]")
        print(f"  N windows: {phase_stats['late']['n_windows']}")
        print(f"  Expected: 0.478")

    # Test for significant difference between phases
    if len(early_phase) > 0 and len(late_phase) > 0:
        t_stat, p_val = stats.ttest_ind(early_corrs, late_corrs)
        phase_stats['phase_comparison'] = {
            't_statistic': t_stat,
            'p_value': p_val,
            'effect_size': (phase_stats['late']['mean_correlation'] -
                           phase_stats['early']['mean_correlation'])
        }

        print(f"\nPhase Comparison (Late vs Early):")
        print(f"  Difference: {phase_stats['phase_comparison']['effect_size']:.3f}")
        print(f"  T-statistic: {phase_stats['phase_comparison']['t_statistic']:.3f}")
        print(f"  P-value: {phase_stats['phase_comparison']['p_value']:.2e}")

    return phase_stats

# ============================================
# STEP 3: Trend Analysis
# ============================================

def analyze_trend(results_df):
    """
    Analyze the trend of correlation change over pseudotime

    Statement claims: r = 0.851, p = 6.98 × 10^-8
    This represents correlation of correlations with pseudotime
    """

    print("\n" + "="*60)
    print("TREND ANALYSIS")
    print("="*60)

    # Calculate trend (correlation of correlations with pseudotime)
    pseudotime_means = results_df['pseudotime_mean'].values
    correlations = results_df['correlation'].values

    # Remove any NaN values
    mask = ~(np.isnan(pseudotime_means) | np.isnan(correlations))
    pseudotime_means = pseudotime_means[mask]
    correlations = correlations[mask]

    if len(correlations) < 3:
        print("ERROR: Not enough windows for trend analysis")
        return None

    # Pearson correlation of correlation values with pseudotime
    trend_r, trend_p = pearsonr(pseudotime_means, correlations)

    # Spearman correlation (non-parametric)
    trend_rho, trend_p_spearman = spearmanr(pseudotime_means, correlations)

    # Linear regression for more details
    from sklearn.linear_model import LinearRegression
    lr = LinearRegression()
    X = pseudotime_means.reshape(-1, 1)
    lr.fit(X, correlations)
    slope = lr.coef_[0]
    intercept = lr.intercept_
    r_squared = lr.score(X, correlations)

    # Bootstrap confidence interval for trend
    n_bootstrap = 1000
    bootstrap_trends = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(len(correlations), len(correlations), replace=True)
        boot_pseudo = pseudotime_means[idx]
        boot_corr = correlations[idx]
        boot_r, _ = pearsonr(boot_pseudo, boot_corr)
        bootstrap_trends.append(boot_r)

    trend_ci_lower = np.percentile(bootstrap_trends, 2.5)
    trend_ci_upper = np.percentile(bootstrap_trends, 97.5)

    results = {
        'trend_correlation': trend_r,
        'trend_p_value': trend_p,
        'trend_spearman': trend_rho,
        'trend_p_spearman': trend_p_spearman,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'trend_ci_lower': trend_ci_lower,
        'trend_ci_upper': trend_ci_upper
    }

    print(f"Trend Statistics:")
    print(f"  Pearson r: {trend_r:.3f} (p = {trend_p:.2e})")
    print(f"  95% CI: [{trend_ci_lower:.3f}, {trend_ci_upper:.3f}]")
    print(f"  Spearman ρ: {trend_rho:.3f} (p = {trend_p_spearman:.2e})")
    print(f"  Linear slope: {slope:.3f}")
    print(f"  R²: {r_squared:.3f}")
    print(f"\nExpected: r = 0.851, p = 6.98e-08")
    print(f"Observed: r = {trend_r:.3f}, p = {trend_p:.2e}")

    return results

# ============================================
# STEP 4: Statistical Validation
# ============================================

def validate_results(results_df, phase_stats, trend_results):
    """
    Validate the sliding window results using permutation tests
    """

    print("\n" + "="*60)
    print("STATISTICAL VALIDATION")
    print("="*60)

    # 1. Check if phase differences are significant
    validation = {}

    if 'early' in phase_stats and 'late' in phase_stats:
        early_mean = phase_stats['early']['mean_correlation']
        late_mean = phase_stats['late']['mean_correlation']
        difference = late_mean - early_mean

        print(f"Phase Difference Validation:")
        print(f"  Early mean: {early_mean:.3f} (expected: -0.417)")
        print(f"  Late mean: {late_mean:.3f} (expected: 0.478)")
        print(f"  Observed difference: {difference:.3f}")
        print(f"  Expected difference: {0.478 - (-0.417):.3f}")

        validation['phase_match'] = (
            abs(early_mean - (-0.417)) < 0.15 and
            abs(late_mean - 0.478) < 0.15
        )

    # 2. Check trend significance
    if trend_results:
        validation['trend_match'] = (
            abs(trend_results['trend_correlation'] - 0.851) < 0.1 and
            trend_results['trend_p_value'] < 1e-6
        )

        print(f"\nTrend Validation:")
        print(f"  Trend correlation match: {validation['trend_match']}")

    # 3. Check window size effect
    print(f"\nWindow Size Analysis:")
    print(f"  Total windows: {len(results_df)}")
    print(f"  Windows with p < 0.05: {sum(results_df['p_value'] < 0.05)}")
    print(f"  Percentage significant: {100*sum(results_df['p_value'] < 0.05)/len(results_df):.1f}%")

    return validation

# ============================================
# STEP 5: Comprehensive Visualization
# ============================================

def visualize_sliding_window(results_df, phase_stats, trend_results):
    """
    Create comprehensive visualization of sliding window analysis
    """

    fig = plt.figure(figsize=(16, 12))

    # 1. Main sliding window plot
    ax1 = plt.subplot(3, 2, 1)
    pseudotime = results_df['pseudotime_mean'].values
    correlations = results_df['correlation'].values

    # Color by correlation strength
    colors = plt.cm.coolwarm((correlations + 1) / 2)  # Normalize to [0, 1]
    scatter = ax1.scatter(pseudotime, correlations, c=correlations,
                         cmap='coolwarm', vmin=-1, vmax=1, s=20, alpha=0.6)

    # Add phase boundaries
    ax1.axvline(x=0.33, color='gray', linestyle='--', alpha=0.5, label='Phase boundaries')
    ax1.axvline(x=0.67, color='gray', linestyle='--', alpha=0.5)

    # Add expected values
    ax1.axhline(y=-0.417, color='blue', linestyle=':', alpha=0.5, label='Expected early')
    ax1.axhline(y=0.478, color='red', linestyle=':', alpha=0.5, label='Expected late')

    # Add trend line
    if trend_results:
        x_trend = np.linspace(pseudotime.min(), pseudotime.max(), 100)
        y_trend = trend_results['slope'] * x_trend + trend_results['intercept']
        ax1.plot(x_trend, y_trend, 'black', linewidth=2, label='Trend line')

    ax1.set_xlabel('Pseudotime')
    ax1.set_ylabel('Correlation (SQSTM1 vs VDAC1)')
    ax1.set_title('Sliding Window Correlation Analysis')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax1, label='Correlation')

    # 2. Phase comparison boxplot
    ax2 = plt.subplot(3, 2, 2)
    phase_data = []
    phase_labels = []

    early = results_df[results_df['pseudotime_mean'] < 0.33]['correlation']
    middle = results_df[(results_df['pseudotime_mean'] >= 0.33) &
                       (results_df['pseudotime_mean'] <= 0.67)]['correlation']
    late = results_df[results_df['pseudotime_mean'] > 0.67]['correlation']

    if len(early) > 0:
        phase_data.append(early)
        phase_labels.append(f'Early\n(n={len(early)})')
    if len(middle) > 0:
        phase_data.append(middle)
        phase_labels.append(f'Middle\n(n={len(middle)})')
    if len(late) > 0:
        phase_data.append(late)
        phase_labels.append(f'Late\n(n={len(late)})')

    bp = ax2.boxplot(phase_data, labels=phase_labels)
    ax2.set_ylabel('Correlation')
    ax2.set_title('Correlation by Phase')
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    ax2.grid(True, alpha=0.3)

    # 3. Histogram of correlations
    ax3 = plt.subplot(3, 2, 3)
    ax3.hist(correlations, bins=30, edgecolor='black', alpha=0.7)
    ax3.axvline(x=0, color='black', linestyle='-', linewidth=2)
    ax3.axvline(x=correlations.mean(), color='red', linestyle='--',
               label=f'Mean = {correlations.mean():.3f}')
    ax3.set_xlabel('Correlation')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Correlations')
    ax3.legend()

    # 4. P-value distribution
    ax4 = plt.subplot(3, 2, 4)
    p_values = results_df['p_value'].values
    ax4.hist(p_values, bins=30, edgecolor='black', alpha=0.7)
    ax4.axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
    ax4.set_xlabel('P-value')
    ax4.set_ylabel('Frequency')
    ax4.set_title('P-value Distribution')
    ax4.legend()

    # 5. Moving average smooth
    ax5 = plt.subplot(3, 2, 5)
    window_smooth = 10
    if len(correlations) > window_smooth:
        smoothed = pd.Series(correlations).rolling(window=window_smooth, center=True).mean()
        ax5.plot(pseudotime, correlations, 'lightgray', alpha=0.5, label='Raw')
        ax5.plot(pseudotime, smoothed, 'blue', linewidth=2, label=f'Smoothed (window={window_smooth})')
    else:
        ax5.plot(pseudotime, correlations, 'blue', linewidth=2, label='Correlations')

    ax5.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax5.set_xlabel('Pseudotime')
    ax5.set_ylabel('Correlation')
    ax5.set_title('Smoothed Correlation Trajectory')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    # 6. Phase means with error bars
    ax6 = plt.subplot(3, 2, 6)
    if phase_stats:
        phases = []
        means = []
        errors = []
        colors_bar = []

        if 'early' in phase_stats:
            phases.append('Early\n(<0.33)')
            means.append(phase_stats['early']['mean_correlation'])
            errors.append(phase_stats['early']['std_correlation'])
            colors_bar.append('blue')

        if 'middle' in phase_stats:
            phases.append('Middle\n(0.33-0.67)')
            means.append(phase_stats['middle']['mean_correlation'])
            errors.append(phase_stats['middle']['std_correlation'])
            colors_bar.append('gray')

        if 'late' in phase_stats:
            phases.append('Late\n(>0.67)')
            means.append(phase_stats['late']['mean_correlation'])
            errors.append(phase_stats['late']['std_correlation'])
            colors_bar.append('red')

        x_pos = np.arange(len(phases))
        ax6.bar(x_pos, means, yerr=errors, capsize=5, color=colors_bar, alpha=0.7)
        ax6.set_xticks(x_pos)
        ax6.set_xticklabels(phases)
        ax6.set_ylabel('Mean Correlation ± SD')
        ax6.set_title('Phase Summary Statistics')
        ax6.axhline(y=0, color='black', linestyle='-', alpha=0.3)

        # Add expected values
        ax6.axhline(y=-0.417, color='blue', linestyle=':', alpha=0.5, label='Expected early')
        ax6.axhline(y=0.478, color='red', linestyle=':', alpha=0.5, label='Expected late')
        ax6.legend()

    plt.tight_layout()
    plt.show()

# ============================================
# STEP 6: Evaluation
# ============================================

def evaluate_statement(phase_stats, trend_results, validation):
    """
    Evaluate the statement based on all analyses
    """

    print("\n" + "="*60)
    print("STATEMENT EVALUATION")
    print("="*60)

    criteria_met = []
    criteria_failed = []

    # Check early phase
    if 'early' in phase_stats:
        early_match = abs(phase_stats['early']['mean_correlation'] - (-0.417)) < 0.15
        if early_match:
            criteria_met.append("Early phase correlation matches (-0.417)")
        else:
            criteria_failed.append(f"Early phase: observed {phase_stats['early']['mean_correlation']:.3f}, expected -0.417")

    # Check late phase
    if 'late' in phase_stats:
        late_match = abs(phase_stats['late']['mean_correlation'] - 0.478) < 0.15
        if late_match:
            criteria_met.append("Late phase correlation matches (0.478)")
        else:
            criteria_failed.append(f"Late phase: observed {phase_stats['late']['mean_correlation']:.3f}, expected 0.478")

    # Check trend
    if trend_results:
        trend_match = (abs(trend_results['trend_correlation'] - 0.851) < 0.15 and
                      trend_results['trend_p_value'] < 1e-6)
        if trend_match:
            criteria_met.append("Trend correlation matches (0.851)")
        else:
            criteria_failed.append(f"Trend: observed r={trend_results['trend_correlation']:.3f}, expected 0.851")

    # Determine evaluation
    total_criteria = len(criteria_met) + len(criteria_failed)
    if total_criteria == 0:
        evaluation = "UNSURE"
        explanation = "Unable to perform complete analysis"
    elif len(criteria_met) >= 2:
        evaluation = "SUPPORTED"
        explanation = f"Key criteria met: {'; '.join(criteria_met)}"
    elif len(criteria_met) == 1:
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"Some criteria met: {'; '.join(criteria_met)}. Failed: {'; '.join(criteria_failed)}"
    else:
        evaluation = "REFUTED"
        explanation = f"Criteria not met: {'; '.join(criteria_failed)}"

    print(f"Criteria Met ({len(criteria_met)}):")
    for c in criteria_met:
        print(f"  ✓ {c}")

    print(f"\nCriteria Failed ({len(criteria_failed)}):")
    for c in criteria_failed:
        print(f"  ✗ {c}")

    print(f"\nFINAL EVALUATION: {evaluation}")
    print(f"EXPLANATION: {explanation}")

    return evaluation, explanation

# ============================================
# MAIN ANALYSIS PIPELINE
# ============================================

def main():
    """
    Complete pipeline for sliding window correlation analysis
    """

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad('data/pool_processed_v2.h5ad')

    # Step 1: Sliding window analysis
    results_df = sliding_window_correlation(adata)

    if results_df is None or len(results_df) == 0:
        return "UNSURE", "Unable to perform sliding window analysis", None

    # Step 2: Phase analysis
    phase_stats = analyze_phases(results_df)

    # Step 3: Trend analysis
    trend_results = analyze_trend(results_df)

    # Step 4: Validation
    validation = validate_results(results_df, phase_stats, trend_results)

    # Step 5: Visualization
    visualize_sliding_window(results_df, phase_stats, trend_results)

    # Step 6: Evaluation
    evaluation, explanation = evaluate_statement(phase_stats, trend_results, validation)

    # Save results
    results_df.to_csv('sliding_window_results.csv', index=False)
    print(f"\nResults saved to sliding_window_results.csv")

    return evaluation, explanation

# ============================================
# KEY CONCEPTS AND NOTES
# ============================================

"""
ADVANCED CONCEPTS IN THIS ANALYSIS:

1. Sliding Window Analysis:
   - Reveals temporal dynamics invisible in static correlations
   - Window size affects resolution vs stability trade-off
   - Similar to LOESS regression (locally weighted smoothing)

2. Dynamic Correlation Interpretation:
   - Negative correlation: Compensatory mechanism
   - Positive correlation: Co-regulation or system failure
   - Transition point: Critical threshold in disease

3. Statistical Considerations:
   - Multiple testing: Many windows = many correlations
   - Temporal autocorrelation: Adjacent windows are not independent
   - Bootstrap for robust confidence intervals

4. ISLP Connections:
   - Chapter 7: Moving beyond linearity (local methods)
   - Chapter 8: Tree-based methods (for change point detection)
   - Chapter 5: Resampling methods (bootstrap)

5. Biological Significance:
   - Early negative correlation: SQSTM1 up when VDAC1 down (compensation)
   - Late positive correlation: Both up or both down together (system failure)
   - Transition: Marks disease progression milestone

TROUBLESHOOTING:
- Small window size: More noise, less stable correlations
- Large window size: Less temporal resolution
- Outliers: Can dominate correlation in small windows
- Edge effects: Fewer windows at extremes of pseudotime

ADVANCED EXTENSIONS:
1. Variable window size based on data density
2. Weighted correlations based on expression levels
3. Multivariate sliding window (include more proteins)
4. Change point detection algorithms
5. Time-varying network analysis
"""

if __name__ == "__main__":
    evaluation, explanation = main()
```

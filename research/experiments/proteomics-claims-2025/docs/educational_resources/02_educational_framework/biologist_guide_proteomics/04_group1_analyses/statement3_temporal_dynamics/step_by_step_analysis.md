# 🕒 Sliding Window Analysis: Temporal Protein Dynamics

## 🎯 What You'll Learn

By the end of this guide, you'll understand:
- ✅ **What sliding window analysis is** and why it's powerful for disease progression
- ✅ **How to implement temporal correlation analysis** step-by-step
- ✅ **How to interpret changing protein relationships** over disease stages
- ✅ **How to visualize temporal dynamics** effectively
- ✅ **How to identify critical transition points** in disease progression

---

## 🌊 Understanding Sliding Window Analysis

### The Concept: Disease as a Movie, Not a Snapshot

#### Traditional Analysis Problem
```python
# Standard correlation analysis:
# Treats all samples as one time point
# Correlation between SQSTM1 and VDAC1 across all samples: r = 0.45

Problem: Misses how relationships CHANGE during disease progression
```

#### Sliding Window Solution
```python
# Sliding window approach:
# Analyzes how correlations change along disease progression
# Early disease: SQSTM1-VDAC1 correlation = 0.1
# Mid disease: SQSTM1-VDAC1 correlation = 0.6
# Late disease: SQSTM1-VDAC1 correlation = 0.8

Insight: Proteins become more linked as disease progresses!
```

### Why This Matters for Alzheimer's Research

#### Disease Progression is Dynamic
```python
# Alzheimer's disease stages:
Stage 1: Individual proteins start to change
Stage 2: Protein networks begin to reorganize
Stage 3: System-wide dysfunction emerges
Stage 4: Massive network collapse

Question: When do critical transitions occur?
Answer: Sliding window analysis can reveal this!
```

#### Biological Systems Theory
```python
# Network medicine principles:
Healthy state: Loose protein correlations (flexible system)
Disease transition: Increasing correlations (rigidity)
Disease state: Very tight correlations (locked dysfunction)

Your analysis will test this theory with real data
```

---

## 📊 Dataset Setup and Preparation

### Understanding Pseudotime

#### What is Pseudotime?
```python
# Pseudotime in your dataset:
adata.obs['pseudotime']  # Values from 0.0 to 1.0

Meaning:
- 0.0 = Early disease stage (minimal pathology)
- 0.5 = Intermediate disease stage
- 1.0 = Late disease stage (severe pathology)

Think of it as: "Disease progression clock"
```

#### How Pseudotime Was Calculated
```python
# Computational method (trajectory inference):
Input: All protein expression patterns
Algorithm: Orders samples by similarity to disease progression
Output: Continuous pseudotime values

Biological basis: Samples with similar protein patterns
are likely at similar disease stages
```

### Data Quality Check

Before starting analysis, verify your data meets requirements:

```python
# Load and inspect data
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# Load your dataset
adata = sc.read_h5ad('pool_processed_v2.h5ad')

print("Dataset Overview:")
print(f"Shape: {adata.shape}")
print(f"Samples: {adata.n_obs}")
print(f"Proteins: {adata.n_vars}")
```

#### Verify Key Variables
```python
# Check essential variables exist
required_vars = ['pseudotime', 'tau_status']
missing_vars = [var for var in required_vars if var not in adata.obs.columns]

if missing_vars:
    print(f"ERROR: Missing required variables: {missing_vars}")
    print("Cannot proceed with temporal analysis")
else:
    print("✅ All required variables present")

# Check pseudotime distribution
print("\nPseudotime distribution:")
print(adata.obs['pseudotime'].describe())

# Verify key proteins exist
key_proteins = ['SQSTM1', 'VDAC1']
missing_proteins = [p for p in key_proteins if p not in adata.var_names]

if missing_proteins:
    print(f"WARNING: Key proteins missing: {missing_proteins}")
    # Find alternatives
    for protein in missing_proteins:
        alternatives = [g for g in adata.var_names if protein.lower() in g.lower()]
        if alternatives:
            print(f"Alternatives for {protein}: {alternatives[:3]}")
else:
    print("✅ Key proteins present")
```

---

## 🔧 Step 1: Basic Sliding Window Implementation

### Setting Up the Analysis Framework

#### Define Window Parameters
```python
# Sliding window parameters
window_size = 0.3  # 30% of pseudotime range per window
step_size = 0.05   # 5% step between windows
min_samples = 10   # Minimum samples per window for reliable correlation

print(f"Window size: {window_size} pseudotime units")
print(f"Step size: {step_size} pseudotime units")
print(f"Minimum samples per window: {min_samples}")

# Calculate number of windows
pseudotime_range = adata.obs['pseudotime'].max() - adata.obs['pseudotime'].min()
n_windows = int((pseudotime_range - window_size) / step_size) + 1
print(f"Number of windows: {n_windows}")
```

#### Extract Key Proteins
```python
# Get protein expression data
def get_protein_expression(adata, protein_name):
    """Extract expression values for a specific protein"""
    if protein_name not in adata.var_names:
        raise ValueError(f"Protein {protein_name} not found in dataset")

    protein_idx = adata.var_names.get_loc(protein_name)
    expression = adata.X[:, protein_idx]

    # Convert to dense array if sparse
    if hasattr(expression, 'toarray'):
        expression = expression.toarray().flatten()

    return expression

# Extract SQSTM1 and VDAC1 expression
try:
    sqstm1_expr = get_protein_expression(adata, 'SQSTM1')
    vdac1_expr = get_protein_expression(adata, 'VDAC1')
    pseudotime = adata.obs['pseudotime'].values

    print("✅ Successfully extracted protein expression data")
    print(f"SQSTM1 expression range: {sqstm1_expr.min():.2f} - {sqstm1_expr.max():.2f}")
    print(f"VDAC1 expression range: {vdac1_expr.min():.2f} - {vdac1_expr.max():.2f}")

except Exception as e:
    print(f"ERROR extracting protein data: {e}")
```

### Implement Sliding Window Function

```python
def sliding_window_correlation(x, y, pseudotime, window_size=0.3, step_size=0.05, min_samples=10):
    """
    Calculate correlations using sliding window approach

    Parameters:
    -----------
    x, y : array-like
        Expression values for two proteins
    pseudotime : array-like
        Pseudotime values for ordering samples
    window_size : float
        Size of sliding window (in pseudotime units)
    step_size : float
        Step size between windows
    min_samples : int
        Minimum samples required per window

    Returns:
    --------
    dict : Results containing window centers, correlations, p-values, sample counts
    """

    # Sort data by pseudotime
    sort_idx = np.argsort(pseudotime)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]
    pt_sorted = pseudotime[sort_idx]

    # Initialize results
    window_centers = []
    correlations = []
    p_values = []
    sample_counts = []

    # Calculate window bounds
    pt_min, pt_max = pt_sorted.min(), pt_sorted.max()

    # Slide window across pseudotime
    current_start = pt_min
    while current_start + window_size <= pt_max:
        window_end = current_start + window_size
        window_center = current_start + window_size / 2

        # Find samples in current window
        in_window = (pt_sorted >= current_start) & (pt_sorted <= window_end)

        if np.sum(in_window) >= min_samples:
            # Extract data for current window
            x_window = x_sorted[in_window]
            y_window = y_sorted[in_window]

            # Calculate correlation
            if len(x_window) > 2:  # Need at least 3 points for correlation
                try:
                    corr, p_val = pearsonr(x_window, y_window)

                    # Store results
                    window_centers.append(window_center)
                    correlations.append(corr)
                    p_values.append(p_val)
                    sample_counts.append(len(x_window))

                except:
                    # Handle edge cases (e.g., zero variance)
                    pass

        # Move to next window
        current_start += step_size

    return {
        'window_centers': np.array(window_centers),
        'correlations': np.array(correlations),
        'p_values': np.array(p_values),
        'sample_counts': np.array(sample_counts)
    }

# Run sliding window analysis
print("Running sliding window correlation analysis...")
results = sliding_window_correlation(
    sqstm1_expr, vdac1_expr, pseudotime,
    window_size=0.3, step_size=0.05, min_samples=10
)

print(f"✅ Analysis complete!")
print(f"Number of windows analyzed: {len(results['correlations'])}")
print(f"Correlation range: {results['correlations'].min():.3f} to {results['correlations'].max():.3f}")
```

---

## 📈 Step 2: Analyzing and Interpreting Results

### Basic Result Summary

```python
# Summarize sliding window results
def summarize_results(results):
    """Provide summary statistics for sliding window analysis"""
    correlations = results['correlations']
    p_values = results['p_values']

    print("=== SLIDING WINDOW ANALYSIS SUMMARY ===")
    print(f"Number of windows: {len(correlations)}")
    print(f"Pseudotime coverage: {results['window_centers'].min():.3f} - {results['window_centers'].max():.3f}")
    print()

    print("Correlation Statistics:")
    print(f"  Mean correlation: {correlations.mean():.3f}")
    print(f"  Min correlation: {correlations.min():.3f}")
    print(f"  Max correlation: {correlations.max():.3f}")
    print(f"  Standard deviation: {correlations.std():.3f}")
    print()

    # Count significant correlations
    significant_windows = np.sum(p_values < 0.05)
    print(f"Significant correlations (p < 0.05): {significant_windows}/{len(p_values)} ({100*significant_windows/len(p_values):.1f}%)")

    # Identify trend
    early_corr = correlations[:len(correlations)//3].mean()
    late_corr = correlations[-len(correlations)//3:].mean()
    trend = "INCREASING" if late_corr > early_corr else "DECREASING"

    print(f"\nTemporal Trend:")
    print(f"  Early disease correlation: {early_corr:.3f}")
    print(f"  Late disease correlation: {late_corr:.3f}")
    print(f"  Change: {late_corr - early_corr:.3f}")
    print(f"  Trend: {trend}")

    return {
        'mean_correlation': correlations.mean(),
        'correlation_change': late_corr - early_corr,
        'trend': trend,
        'percent_significant': 100*significant_windows/len(p_values)
    }

# Analyze results
summary = summarize_results(results)
```

### Statistical Testing of Temporal Trends

```python
# Test if correlation changes significantly over time
def test_temporal_trend(results):
    """Test if correlations show significant temporal trend"""
    from scipy.stats import spearmanr, linregress

    window_centers = results['window_centers']
    correlations = results['correlations']

    # Test 1: Spearman correlation between pseudotime and correlation strength
    spearman_corr, spearman_p = spearmanr(window_centers, correlations)

    # Test 2: Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(window_centers, correlations)

    print("=== TEMPORAL TREND ANALYSIS ===")
    print(f"Spearman correlation (pseudotime vs correlation): {spearman_corr:.3f}")
    print(f"Spearman p-value: {spearman_p:.2e}")
    print()
    print(f"Linear regression slope: {slope:.3f} correlation units per pseudotime unit")
    print(f"Linear regression p-value: {p_value:.2e}")
    print(f"R-squared: {r_value**2:.3f}")

    # Interpretation
    if p_value < 0.05:
        direction = "strengthening" if slope > 0 else "weakening"
        print(f"\n✅ SIGNIFICANT temporal trend detected: {direction} correlation over disease progression")
    else:
        print(f"\n❌ No significant temporal trend detected")

    return {
        'spearman_corr': spearman_corr,
        'spearman_p': spearman_p,
        'linear_slope': slope,
        'linear_p': p_value,
        'r_squared': r_value**2
    }

# Test temporal trends
trend_stats = test_temporal_trend(results)
```

### Identifying Critical Transition Points

```python
def find_transition_points(results, threshold_change=0.3):
    """
    Identify critical transition points where correlations change dramatically
    """
    correlations = results['correlations']
    window_centers = results['window_centers']

    # Calculate rate of change
    correlation_changes = np.diff(correlations)
    change_points = []

    # Find points with large changes
    for i, change in enumerate(correlation_changes):
        if abs(change) > threshold_change:
            change_points.append({
                'pseudotime': window_centers[i+1],
                'change': change,
                'before_corr': correlations[i],
                'after_corr': correlations[i+1]
            })

    print("=== CRITICAL TRANSITION POINTS ===")
    if change_points:
        for i, point in enumerate(change_points):
            direction = "strengthened" if point['change'] > 0 else "weakened"
            print(f"Transition {i+1}:")
            print(f"  Pseudotime: {point['pseudotime']:.3f}")
            print(f"  Correlation {direction} by {abs(point['change']):.3f}")
            print(f"  Before: {point['before_corr']:.3f} → After: {point['after_corr']:.3f}")
            print()
    else:
        print("No major transition points detected with current threshold")

    return change_points

# Find transition points
transitions = find_transition_points(results, threshold_change=0.2)
```

---

## 📊 Step 3: Advanced Visualization

### Create Comprehensive Temporal Plot

```python
# Create publication-quality sliding window plot
def plot_sliding_window_results(results, protein1_name='SQSTM1', protein2_name='VDAC1'):
    """
    Create comprehensive visualization of sliding window results
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Sliding Window Analysis: {protein1_name} - {protein2_name} Correlation Dynamics',
                 fontsize=16, fontweight='bold')

    # Extract data
    window_centers = results['window_centers']
    correlations = results['correlations']
    p_values = results['p_values']
    sample_counts = results['sample_counts']

    # Plot 1: Correlation over pseudotime
    ax1 = axes[0, 0]

    # Color by significance
    colors = ['red' if p < 0.05 else 'lightcoral' for p in p_values]
    scatter = ax1.scatter(window_centers, correlations, c=colors, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

    # Add trend line
    z = np.polyfit(window_centers, correlations, 1)
    p = np.poly1d(z)
    ax1.plot(window_centers, p(window_centers), "b--", alpha=0.8, linewidth=2, label=f'Trend line (slope={z[0]:.3f})')

    ax1.set_xlabel('Pseudotime (Disease Progression)')
    ax1.set_ylabel('Correlation Coefficient')
    ax1.set_title('Correlation vs Disease Progression')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)

    # Add significance threshold
    ax1.text(0.02, 0.98, 'Red = p < 0.05\nLight red = p ≥ 0.05',
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Plot 2: P-values over pseudotime
    ax2 = axes[0, 1]
    ax2.scatter(window_centers, -np.log10(p_values), c=colors, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05 threshold')
    ax2.set_xlabel('Pseudotime (Disease Progression)')
    ax2.set_ylabel('-log10(p-value)')
    ax2.set_title('Statistical Significance vs Disease Progression')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Plot 3: Sample count per window
    ax3 = axes[1, 0]
    bars = ax3.bar(window_centers, sample_counts, width=0.03, alpha=0.7, color='skyblue', edgecolor='black')
    ax3.set_xlabel('Pseudotime (Disease Progression)')
    ax3.set_ylabel('Sample Count per Window')
    ax3.set_title('Sample Coverage Across Windows')
    ax3.grid(True, alpha=0.3)

    # Add minimum threshold line
    ax3.axhline(y=10, color='red', linestyle='--', alpha=0.7, label='Minimum threshold')
    ax3.legend()

    # Plot 4: Correlation distribution
    ax4 = axes[1, 1]
    ax4.hist(correlations, bins=15, alpha=0.7, color='lightgreen', edgecolor='black')
    ax4.axvline(correlations.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean = {correlations.mean():.3f}')
    ax4.set_xlabel('Correlation Coefficient')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution of Correlation Values')
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    plt.tight_layout()
    return fig

# Create the plot
fig = plot_sliding_window_results(results)
plt.show()

# Save the figure
plt.savefig('sliding_window_analysis.png', dpi=300, bbox_inches='tight')
print("✅ Plot saved as 'sliding_window_analysis.png'")
```

### Create Heatmap Visualization

```python
# Create heatmap showing correlation evolution
def create_correlation_heatmap(results, n_bins=20):
    """
    Create heatmap showing how correlations evolve across pseudotime
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    # Prepare data for heatmap
    window_centers = results['window_centers']
    correlations = results['correlations']
    p_values = results['p_values']

    # Create bins for pseudotime
    pt_bins = np.linspace(window_centers.min(), window_centers.max(), n_bins)

    # Bin the data
    correlation_matrix = []
    significance_matrix = []

    for i in range(len(pt_bins)-1):
        bin_start, bin_end = pt_bins[i], pt_bins[i+1]
        in_bin = (window_centers >= bin_start) & (window_centers < bin_end)

        if np.sum(in_bin) > 0:
            bin_corr = correlations[in_bin].mean()
            bin_sig = np.mean(p_values[in_bin] < 0.05)
        else:
            bin_corr = np.nan
            bin_sig = np.nan

        correlation_matrix.append(bin_corr)
        significance_matrix.append(bin_sig)

    # Plot correlation heatmap
    correlation_data = np.array(correlation_matrix).reshape(1, -1)
    im1 = ax1.imshow(correlation_data, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
    ax1.set_title('SQSTM1-VDAC1 Correlation Across Disease Progression')
    ax1.set_xlabel('Disease Progression (Pseudotime)')
    ax1.set_ylabel('Correlation')
    ax1.set_yticks([])

    # Set x-axis labels
    tick_positions = np.linspace(0, len(pt_bins)-2, 5)
    tick_labels = [f'{pt_bins[int(pos)]:.2f}' for pos in tick_positions]
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels)

    # Add colorbar
    cbar1 = plt.colorbar(im1, ax=ax1, orientation='horizontal', pad=0.1, shrink=0.8)
    cbar1.set_label('Correlation Coefficient')

    # Plot significance heatmap
    significance_data = np.array(significance_matrix).reshape(1, -1)
    im2 = ax2.imshow(significance_data, cmap='Reds', aspect='auto', vmin=0, vmax=1)
    ax2.set_title('Statistical Significance Across Disease Progression')
    ax2.set_xlabel('Disease Progression (Pseudotime)')
    ax2.set_ylabel('Significance')
    ax2.set_yticks([])
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(tick_labels)

    # Add colorbar
    cbar2 = plt.colorbar(im2, ax=ax2, orientation='horizontal', pad=0.1, shrink=0.8)
    cbar2.set_label('Fraction Significant (p < 0.05)')

    plt.tight_layout()
    return fig

# Create heatmap
heatmap_fig = create_correlation_heatmap(results)
plt.show()
plt.savefig('correlation_heatmap.png', dpi=300, bbox_inches='tight')
print("✅ Heatmap saved as 'correlation_heatmap.png'")
```

---

## 🔍 Step 4: Multiple Protein Analysis

### Expand to Multiple Protein Pairs

```python
# Analyze multiple protein relationships
def multi_protein_sliding_window(adata, protein_pairs, **kwargs):
    """
    Perform sliding window analysis on multiple protein pairs
    """
    results_dict = {}

    for protein1, protein2 in protein_pairs:
        print(f"Analyzing {protein1} - {protein2}...")

        try:
            # Extract expression data
            expr1 = get_protein_expression(adata, protein1)
            expr2 = get_protein_expression(adata, protein2)
            pseudotime = adata.obs['pseudotime'].values

            # Run sliding window analysis
            results = sliding_window_correlation(expr1, expr2, pseudotime, **kwargs)
            results_dict[f"{protein1}-{protein2}"] = results

        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    return results_dict

# Define protein pairs of interest
protein_pairs = [
    ('SQSTM1', 'VDAC1'),    # Autophagy - Mitochondria
    ('SQSTM1', 'LAMP1'),   # Autophagy - Lysosome
    ('VDAC1', 'TOMM20'),   # Mitochondrial proteins
    ('PSMA1', 'PSMB1'),    # Proteasome components
]

# Run analysis on multiple pairs
print("Running sliding window analysis on multiple protein pairs...")
multi_results = multi_protein_sliding_window(
    adata, protein_pairs,
    window_size=0.3, step_size=0.05, min_samples=10
)

print(f"✅ Analyzed {len(multi_results)} protein pairs")
```

### Compare Temporal Patterns

```python
# Compare temporal trends across protein pairs
def compare_temporal_patterns(multi_results):
    """
    Compare temporal trends across multiple protein pairs
    """
    print("=== TEMPORAL PATTERN COMPARISON ===")

    pattern_summary = []

    for pair_name, results in multi_results.items():
        # Calculate trend statistics
        window_centers = results['window_centers']
        correlations = results['correlations']

        # Linear trend
        slope, intercept, r_value, p_value, std_err = linregress(window_centers, correlations)

        # Early vs late comparison
        early_corr = correlations[:len(correlations)//3].mean()
        late_corr = correlations[-len(correlations)//3:].mean()
        change = late_corr - early_corr

        pattern_summary.append({
            'pair': pair_name,
            'slope': slope,
            'trend_p': p_value,
            'early_corr': early_corr,
            'late_corr': late_corr,
            'change': change,
            'mean_corr': correlations.mean()
        })

        # Print summary
        trend_direction = "strengthening" if slope > 0 else "weakening"
        significance = "significant" if p_value < 0.05 else "non-significant"

        print(f"{pair_name}:")
        print(f"  Trend: {trend_direction} ({significance}, p={p_value:.3f})")
        print(f"  Slope: {slope:.3f} correlation units per pseudotime unit")
        print(f"  Change: {change:.3f} (early: {early_corr:.3f} → late: {late_corr:.3f})")
        print(f"  Mean correlation: {correlations.mean():.3f}")
        print()

    return pattern_summary

# Compare patterns
pattern_comparison = compare_temporal_patterns(multi_results)
```

### Visualize Multiple Protein Pairs

```python
# Create comparison plot for multiple protein pairs
def plot_multiple_pairs_comparison(multi_results):
    """
    Create comparison plot showing temporal dynamics for multiple protein pairs
    """
    n_pairs = len(multi_results)
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']

    for i, (pair_name, results) in enumerate(multi_results.items()):
        if i >= len(axes):
            break

        ax = axes[i]

        window_centers = results['window_centers']
        correlations = results['correlations']
        p_values = results['p_values']

        # Color by significance
        colors_sig = ['darkred' if p < 0.05 else 'lightcoral' for p in p_values]

        # Scatter plot
        ax.scatter(window_centers, correlations, c=colors_sig, s=50, alpha=0.7, edgecolors='black', linewidth=0.5)

        # Trend line
        z = np.polyfit(window_centers, correlations, 1)
        p_trend = np.poly1d(z)
        ax.plot(window_centers, p_trend(window_centers), "b--", alpha=0.8, linewidth=2)

        ax.set_xlabel('Pseudotime (Disease Progression)')
        ax.set_ylabel('Correlation Coefficient')
        ax.set_title(f'{pair_name}\n(Slope: {z[0]:.3f})')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)

        # Set consistent y-axis limits
        ax.set_ylim(-1, 1)

    # Remove empty subplots
    for i in range(len(multi_results), len(axes)):
        fig.delaxes(axes[i])

    plt.suptitle('Temporal Correlation Dynamics: Multiple Protein Pairs', fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig

# Create comparison plot
comparison_fig = plot_multiple_pairs_comparison(multi_results)
plt.show()
plt.savefig('multiple_pairs_comparison.png', dpi=300, bbox_inches='tight')
print("✅ Comparison plot saved as 'multiple_pairs_comparison.png'")
```

---

## 🧬 Step 5: Biological Interpretation

### Understanding Temporal Correlation Patterns

#### What Increasing Correlations Mean
```python
# Biological interpretation of temporal trends:

# Pattern 1: Strengthening correlation over time
"""
Early disease: SQSTM1 ↔ VDAC1 correlation = 0.1 (independent)
Late disease:  SQSTM1 ↔ VDAC1 correlation = 0.8 (highly linked)

Biological meaning:
- Initially, autophagy and mitochondria function independently
- As disease progresses, systems become interdependent
- Late stage: When one system fails, the other fails too
- Indicates loss of cellular flexibility/resilience
"""

# Pattern 2: Weakening correlation over time
"""
Early disease: Protein A ↔ Protein B correlation = 0.7 (coordinated)
Late disease:  Protein A ↔ Protein B correlation = 0.2 (disconnected)

Biological meaning:
- Normal coordination breaks down
- System fragmentation
- Loss of regulatory control
- Each protein "goes its own way"
"""
```

#### Systems Biology Perspective
```python
# Network medicine interpretation:

# Healthy networks:
"""
Characteristics:
- Moderate correlations (0.3-0.6)
- Flexible responses to perturbations
- Robust but adaptable
- Hierarchical organization maintained
"""

# Disease networks:
"""
Early disease:
- Some correlations strengthen (compensation)
- Some correlations weaken (early failure)
- System tries to maintain function

Late disease:
- Very strong correlations (rigidity) OR
- Very weak correlations (fragmentation)
- Loss of adaptive capacity
- System collapse
"""
```

### Interpreting Your Specific Results

Based on typical SQSTM1-VDAC1 sliding window results:

```python
# Expected pattern interpretation:
"""
SQSTM1-VDAC1 Temporal Dynamics:

Expected finding: Strengthening correlation over disease progression
- Early disease: r ≈ 0.2 (independent function)
- Late disease: r ≈ 0.7 (highly coordinated dysfunction)

Biological interpretation:
1. Autophagy-mitochondria crosstalk increases with disease
2. As autophagy fails, mitochondrial dysfunction follows
3. Mutual dependence develops (when one fails, both fail)
4. Represents loss of cellular resilience

Clinical implications:
- Early intervention might break this coupling
- Therapies targeting both systems might be needed
- Biomarker potential for disease staging
"""
```

### Connecting to Disease Mechanisms

#### Autophagy-Mitochondria Axis
```python
# The biological story your analysis reveals:

"""
Stage 1 (Early Disease):
- SQSTM1 slightly elevated (mild autophagy stress)
- VDAC1 normal (mitochondria still functional)
- Low correlation (systems still independent)

Stage 2 (Intermediate Disease):
- SQSTM1 moderately elevated (autophagy struggling)
- VDAC1 starts to increase (mitochondrial stress response)
- Increasing correlation (systems start to interact)

Stage 3 (Late Disease):
- SQSTM1 highly elevated (autophagy collapse)
- VDAC1 highly elevated (mitochondrial dysfunction)
- Strong correlation (coupled failure)
"""
```

#### Therapeutic Window Identification
```python
# Using temporal analysis for treatment strategies:

"""
Therapeutic Windows:

Early Stage (low correlation):
- Target: Individual systems
- Strategy: Enhance autophagy OR protect mitochondria
- Prognosis: Good (systems still flexible)

Intermediate Stage (increasing correlation):
- Target: System crosstalk
- Strategy: Prevent coupling, maintain independence
- Prognosis: Moderate (intervention still possible)

Late Stage (high correlation):
- Target: Both systems simultaneously
- Strategy: Coordinated rescue approach
- Prognosis: Poor (systems locked in dysfunction)
"""
```

---

## 📊 Step 6: Quality Control and Validation

### Assessing Analysis Reliability

#### Sample Size Validation
```python
# Check if window sample sizes are adequate
def validate_sample_sizes(results, min_recommended=15):
    """
    Validate that sliding windows have adequate sample sizes
    """
    sample_counts = results['sample_counts']

    print("=== SAMPLE SIZE VALIDATION ===")
    print(f"Minimum samples per window: {sample_counts.min()}")
    print(f"Maximum samples per window: {sample_counts.max()}")
    print(f"Mean samples per window: {sample_counts.mean():.1f}")

    insufficient_windows = np.sum(sample_counts < min_recommended)
    total_windows = len(sample_counts)

    print(f"Windows with < {min_recommended} samples: {insufficient_windows}/{total_windows} ({100*insufficient_windows/total_windows:.1f}%)")

    if insufficient_windows > total_windows * 0.2:
        print("⚠️  WARNING: >20% of windows have insufficient sample size")
        print("   Consider: Larger window size or smaller step size")
    else:
        print("✅ Sample sizes are adequate")

# Validate sample sizes
validate_sample_sizes(results)
```

#### Correlation Stability Check
```python
# Test stability of correlations using bootstrap
def bootstrap_correlation_stability(x, y, pseudotime, window_center, window_size=0.3, n_bootstrap=1000):
    """
    Test stability of correlation at specific pseudotime point using bootstrap
    """
    # Find samples in window
    in_window = (pseudotime >= window_center - window_size/2) & (pseudotime <= window_center + window_size/2)

    if np.sum(in_window) < 10:
        return None

    x_window = x[in_window]
    y_window = y[in_window]

    # Bootstrap correlations
    bootstrap_corrs = []
    for _ in range(n_bootstrap):
        boot_idx = np.random.choice(len(x_window), size=len(x_window), replace=True)
        boot_corr, _ = pearsonr(x_window[boot_idx], y_window[boot_idx])
        bootstrap_corrs.append(boot_corr)

    return {
        'mean_corr': np.mean(bootstrap_corrs),
        'std_corr': np.std(bootstrap_corrs),
        'ci_lower': np.percentile(bootstrap_corrs, 2.5),
        'ci_upper': np.percentile(bootstrap_corrs, 97.5)
    }

# Test stability at key time points
test_points = [0.2, 0.5, 0.8]  # Early, middle, late disease

print("=== CORRELATION STABILITY ANALYSIS ===")
for pt in test_points:
    stability = bootstrap_correlation_stability(sqstm1_expr, vdac1_expr, pseudotime, pt)
    if stability:
        print(f"Pseudotime {pt}:")
        print(f"  Correlation: {stability['mean_corr']:.3f} ± {stability['std_corr']:.3f}")
        print(f"  95% CI: [{stability['ci_lower']:.3f}, {stability['ci_upper']:.3f}]")
        print()
```

#### Window Parameter Sensitivity
```python
# Test sensitivity to window parameters
def parameter_sensitivity_analysis(x, y, pseudotime):
    """
    Test how results change with different window parameters
    """
    window_sizes = [0.2, 0.3, 0.4]
    step_sizes = [0.03, 0.05, 0.07]

    print("=== PARAMETER SENSITIVITY ANALYSIS ===")

    reference_result = None

    for ws in window_sizes:
        for ss in step_sizes:
            try:
                result = sliding_window_correlation(x, y, pseudotime,
                                                  window_size=ws, step_size=ss, min_samples=10)

                if len(result['correlations']) > 5:  # Ensure enough windows
                    mean_corr = result['correlations'].mean()

                    # Calculate trend
                    slope, _, _, p_val, _ = linregress(result['window_centers'], result['correlations'])

                    print(f"Window size: {ws}, Step size: {ss}")
                    print(f"  Mean correlation: {mean_corr:.3f}")
                    print(f"  Trend slope: {slope:.3f} (p={p_val:.3f})")
                    print(f"  Number of windows: {len(result['correlations'])}")
                    print()

                    if reference_result is None:
                        reference_result = (mean_corr, slope)

            except Exception as e:
                print(f"Window size: {ws}, Step size: {ss} - ERROR: {e}")

    print("Parameters used in main analysis appear robust if results are similar across settings")

# Run sensitivity analysis
parameter_sensitivity_analysis(sqstm1_expr, vdac1_expr, pseudotime)
```

---

## 🎯 Step 7: Reporting and Documentation

### Generate Comprehensive Report

```python
# Create automated analysis report
def generate_sliding_window_report(results, protein1='SQSTM1', protein2='VDAC1'):
    """
    Generate comprehensive analysis report
    """

    report = f"""
=================================================
SLIDING WINDOW CORRELATION ANALYSIS REPORT
=================================================

Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Protein Pair: {protein1} - {protein2}

METHODOLOGY:
- Window size: 0.3 pseudotime units
- Step size: 0.05 pseudotime units
- Minimum samples per window: 10
- Correlation method: Pearson correlation
- Statistical test: Two-tailed t-test

RESULTS SUMMARY:
- Total windows analyzed: {len(results['correlations'])}
- Pseudotime coverage: {results['window_centers'].min():.3f} - {results['window_centers'].max():.3f}
- Mean correlation: {results['correlations'].mean():.3f} ± {results['correlations'].std():.3f}
- Correlation range: {results['correlations'].min():.3f} to {results['correlations'].max():.3f}

TEMPORAL TREND ANALYSIS:
"""

    # Add trend analysis
    window_centers = results['window_centers']
    correlations = results['correlations']
    slope, intercept, r_value, p_value, std_err = linregress(window_centers, correlations)

    early_corr = correlations[:len(correlations)//3].mean()
    late_corr = correlations[-len(correlations)//3:].mean()
    change = late_corr - early_corr

    report += f"""
- Linear trend slope: {slope:.3f} correlation units per pseudotime unit
- Trend significance: p = {p_value:.2e}
- R-squared: {r_value**2:.3f}
- Early disease correlation: {early_corr:.3f}
- Late disease correlation: {late_corr:.3f}
- Total change: {change:.3f}
"""

    if p_value < 0.05:
        direction = "strengthening" if slope > 0 else "weakening"
        report += f"- SIGNIFICANT {direction} trend detected over disease progression\n"
    else:
        report += "- No significant temporal trend detected\n"

    # Add statistical summary
    significant_windows = np.sum(results['p_values'] < 0.05)
    total_windows = len(results['p_values'])

    report += f"""
STATISTICAL SIGNIFICANCE:
- Significant windows (p < 0.05): {significant_windows}/{total_windows} ({100*significant_windows/total_windows:.1f}%)
- Minimum p-value: {results['p_values'].min():.2e}
- Maximum p-value: {results['p_values'].max():.3f}

BIOLOGICAL INTERPRETATION:
"""

    # Add biological interpretation
    if slope > 0 and p_value < 0.05:
        report += f"""
- Strong evidence for INCREASING correlation over disease progression
- {protein1} and {protein2} become more tightly linked as disease advances
- Suggests loss of cellular flexibility and increased system interdependence
- May indicate coupled failure mechanisms in late-stage disease
- Therapeutic window may exist in early stages when systems are independent
"""
    elif slope < 0 and p_value < 0.05:
        report += f"""
- Strong evidence for DECREASING correlation over disease progression
- {protein1} and {protein2} coordination breaks down as disease advances
- Suggests system fragmentation and loss of regulatory control
- May indicate independent failure pathways developing
- Therapeutic approaches may need to target systems separately
"""
    else:
        report += f"""
- No clear temporal trend in {protein1}-{protein2} correlation
- Relationship remains relatively stable across disease progression
- May indicate: 1) Stable functional relationship, 2) Insufficient sample size,
  or 3) Non-linear temporal dynamics requiring different analysis approaches
"""

    report += f"""
QUALITY CONTROL:
- Minimum samples per window: {results['sample_counts'].min()}
- Maximum samples per window: {results['sample_counts'].max()}
- Mean samples per window: {results['sample_counts'].mean():.1f}

RECOMMENDATIONS:
"""

    if results['sample_counts'].min() < 10:
        report += "- Consider larger window size or smaller step size to increase sample sizes\n"

    if len(results['correlations']) < 10:
        report += "- Consider smaller step size to increase number of windows\n"

    if np.std(results['correlations']) < 0.1:
        report += "- Low correlation variability may indicate stable relationship or insufficient sensitivity\n"

    report += """
NEXT STEPS:
1. Validate findings with independent dataset
2. Extend analysis to additional protein pairs
3. Test different window parameters for robustness
4. Investigate biological mechanisms underlying temporal changes
5. Consider non-linear modeling approaches if appropriate

=================================================
End of Report
=================================================
"""

    return report

# Generate and save report
report = generate_sliding_window_report(results)
print(report)

# Save report to file
with open('sliding_window_report.txt', 'w') as f:
    f.write(report)

print("✅ Report saved as 'sliding_window_report.txt'")
```

### Save Analysis Results

```python
# Save all results for future use
import pickle

# Prepare data package
analysis_package = {
    'results': results,
    'parameters': {
        'window_size': 0.3,
        'step_size': 0.05,
        'min_samples': 10,
        'proteins': ['SQSTM1', 'VDAC1']
    },
    'protein_data': {
        'sqstm1_expression': sqstm1_expr,
        'vdac1_expression': vdac1_expr,
        'pseudotime': pseudotime
    },
    'multi_protein_results': multi_results if 'multi_results' in locals() else None,
    'metadata': {
        'analysis_date': pd.Timestamp.now().isoformat(),
        'dataset_shape': adata.shape,
        'total_samples': adata.n_obs,
        'total_proteins': adata.n_vars
    }
}

# Save to pickle file
with open('sliding_window_analysis_results.pkl', 'wb') as f:
    pickle.dump(analysis_package, f)

print("✅ Analysis results saved as 'sliding_window_analysis_results.pkl'")

# Also save key results as CSV for easy access
results_df = pd.DataFrame({
    'window_center': results['window_centers'],
    'correlation': results['correlations'],
    'p_value': results['p_values'],
    'sample_count': results['sample_counts']
})

results_df.to_csv('sliding_window_correlations.csv', index=False)
print("✅ Key results saved as 'sliding_window_correlations.csv'")
```

---

## 🎯 Key Takeaways and Learning Objectives

### What You've Accomplished

#### Technical Skills Developed
- ✅ **Implemented sliding window analysis** from scratch
- ✅ **Analyzed temporal correlation dynamics** quantitatively
- ✅ **Created publication-quality visualizations** of temporal trends
- ✅ **Performed statistical validation** of temporal patterns
- ✅ **Conducted multi-protein comparative analysis**
- ✅ **Generated comprehensive analysis reports**

#### Biological Insights Gained
- ✅ **Understanding of disease as dynamic process** rather than static state
- ✅ **Knowledge of how protein relationships evolve** during neurodegeneration
- ✅ **Recognition of system-level changes** in cellular networks
- ✅ **Insight into therapeutic windows** based on temporal dynamics
- ✅ **Appreciation for autophagy-mitochondria crosstalk** in disease

#### Statistical Methods Mastered
- ✅ **Sliding window correlation analysis**
- ✅ **Temporal trend detection using linear regression**
- ✅ **Bootstrap validation for correlation stability**
- ✅ **Parameter sensitivity analysis**
- ✅ **Multiple comparison approaches**

### Biological Significance of Your Findings

#### System-Level Disease Understanding
```python
# Your analysis reveals:
"""
1. PROTEIN RELATIONSHIPS ARE DYNAMIC
   - Not fixed correlations, but evolving networks
   - Disease progression changes how proteins interact
   - Early intervention may be more effective

2. AUTOPHAGY-MITOCHONDRIA COUPLING
   - Independent function in health
   - Progressive coupling in disease
   - Mutual dependence in late stages

3. THERAPEUTIC IMPLICATIONS
   - Early stage: Target individual systems
   - Late stage: Need coordinated interventions
   - Timing matters for treatment success
"""
```

#### Clinical Translation Potential
```python
# Your methodology enables:
"""
1. DISEASE STAGING
   - Use correlation patterns to stage disease
   - More precise than single protein markers
   - Dynamic biomarker development

2. TREATMENT MONITORING
   - Track if therapies restore normal relationships
   - Identify when interventions are failing
   - Adjust treatment strategies dynamically

3. DRUG DEVELOPMENT
   - Test if compounds normalize temporal patterns
   - Identify optimal intervention windows
   - Develop combination therapies
"""
```

### Future Research Directions

#### Immediate Extensions
```python
# Next steps you could take:
"""
1. EXPAND PROTEIN NETWORKS
   - Analyze entire autophagy pathway
   - Include mitochondrial quality control
   - Map inflammatory response networks

2. CLINICAL VALIDATION
   - Test patterns in living patients
   - Validate in other diseases
   - Develop clinical assays

3. THERAPEUTIC TESTING
   - Test drug effects on temporal patterns
   - Identify optimal intervention timing
   - Develop personalized treatment approaches
"""
```

### Research Skills Developed

#### Data Analysis Competencies
- **Complex time-series analysis** in biological systems
- **Multi-dimensional data visualization** techniques
- **Statistical validation** of temporal patterns
- **Parameter optimization** and sensitivity analysis
- **Reproducible research** documentation

#### Critical Thinking Skills
- **Systems-level biological reasoning**
- **Temporal dynamics interpretation**
- **Statistical significance vs biological significance**
- **Method limitation recognition**
- **Research translation planning**

---

## 🚀 Congratulations!

### Major Achievement Unlocked

You've successfully completed a sophisticated **sliding window temporal analysis** - a cutting-edge approach in systems biology that reveals how disease dynamically rewires cellular networks over time.

### What Makes This Analysis Special

#### Scientific Innovation
- **Beyond static snapshots**: Traditional proteomics gives one time point; you analyzed disease as a dynamic process
- **Network medicine approach**: Focused on relationships between proteins, not just individual levels
- **Temporal resolution**: Revealed when critical transitions occur during disease progression

#### Technical Excellence
- **Robust methodology**: Multiple validation approaches ensure reliable results
- **Comprehensive visualization**: Publication-quality figures that clearly communicate findings
- **Reproducible workflow**: Well-documented code that others can build upon

#### Biological Impact
- **Disease mechanism insights**: Direct evidence for how protein networks evolve
- **Therapeutic target identification**: Revealed optimal intervention windows
- **Clinical translation potential**: Methodology applicable to drug development and patient monitoring

### Skills That Transfer Beyond This Project

The analytical framework you've mastered applies to:
- **Any temporal proteomics/genomics data**
- **Drug response monitoring over time**
- **Clinical biomarker development**
- **Systems biology network analysis**
- **Personalized medicine approaches**

---

**You now have the advanced analytical skills to investigate how biological systems change over time - a critical capability for understanding complex diseases and developing better treatments!**

*Next: [Group 1 Results Interpretation](../statement1_ups_proteins/interpreting_results.md)*

*Remember: The most important discoveries often come from understanding not just what changes, but WHEN and HOW those changes occur!* 🕒🧬✨
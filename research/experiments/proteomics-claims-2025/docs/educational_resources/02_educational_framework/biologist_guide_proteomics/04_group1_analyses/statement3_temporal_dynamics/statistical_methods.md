# 📊 Statistical Methods for Temporal Dynamics Analysis

## 🎯 Learning Objectives

By the end of this guide, you'll understand:
- ✅ **Sliding window correlation** principles and implementation
- ✅ **Pseudotime analysis** for disease trajectory inference
- ✅ **Statistical testing** for temporal trends
- ✅ **Bootstrap methods** for confidence intervals
- ✅ **Change point detection** algorithms
- ✅ **Multiple testing considerations** in temporal analysis

---

## 🪟 Sliding Window Analysis Fundamentals

### Core Concept

Sliding window analysis divides continuous data into overlapping segments to track how statistical properties change over time or disease progression.

#### Mathematical Definition
```python
# Sliding Window Parameters:
"""
W = Window size (fraction of total range)
S = Step size (how much window moves)
N = Minimum samples per window

For each window position i:
  - Select samples in window
  - Calculate statistic (e.g., correlation)
  - Move window by step size
  - Repeat until end
"""
```

### Window Size Selection

#### Trade-offs in Window Size
```python
# Small Windows (W = 0.1-0.2):
"""
Advantages:
✓ High temporal resolution
✓ Detect rapid changes
✓ Local patterns visible

Disadvantages:
✗ Fewer samples per window
✗ Less statistical power
✗ Noisier estimates
"""

# Large Windows (W = 0.4-0.5):
"""
Advantages:
✓ More samples per window
✓ Greater statistical power
✓ Smoother estimates

Disadvantages:
✗ Low temporal resolution
✗ May miss rapid changes
✗ Averaging over transitions
"""
```

#### Optimal Window Size Formula
```python
def calculate_optimal_window_size(n_samples, n_timepoints, min_samples=10):
    """
    Empirical formula for window size selection

    Based on:
    - Total sample size
    - Desired resolution
    - Statistical requirements
    """
    # Ensure adequate samples per window
    min_window = min_samples / n_samples

    # Balance resolution and power
    optimal_window = 0.3  # Default: 30% of range

    # Adjust based on sample density
    if n_samples < 50:
        optimal_window = max(0.4, min_window)  # Larger windows
    elif n_samples > 200:
        optimal_window = 0.2  # Smaller windows OK

    return optimal_window
```

### Step Size Optimization

#### Overlap Considerations
```python
# High Overlap (Step = 0.01-0.05):
"""
✓ Smooth trajectories
✓ Don't miss transitions
✗ Highly correlated windows
✗ Multiple testing burden
"""

# Low Overlap (Step = 0.2-0.3):
"""
✓ Independent windows
✓ Fewer tests
✗ May miss transitions
✗ Choppy trajectories
"""

# Recommended: 80-90% overlap
step_size = window_size * 0.1  # 90% overlap
```

---

## 📈 Correlation Analysis Methods

### Pearson Correlation in Windows

#### Standard Approach
```python
def windowed_pearson_correlation(x, y, time, window_size, step_size):
    """
    Calculate Pearson correlation in sliding windows
    """
    correlations = []
    window_centers = []
    p_values = []

    time_min, time_max = time.min(), time.max()
    current_start = time_min

    while current_start + window_size <= time_max:
        # Select samples in window
        in_window = (time >= current_start) & (time < current_start + window_size)

        if sum(in_window) >= 10:  # Minimum samples
            # Calculate correlation
            r, p = pearsonr(x[in_window], y[in_window])

            correlations.append(r)
            p_values.append(p)
            window_centers.append(current_start + window_size/2)

        current_start += step_size

    return correlations, p_values, window_centers
```

### Spearman Correlation for Robustness

#### When to Use Spearman
```python
# Use Spearman when:
"""
1. Non-linear monotonic relationships
2. Outliers present
3. Non-normal distributions
4. Ordinal data
"""

def windowed_spearman_correlation(x, y, time, window_size, step_size):
    """
    Rank-based correlation for robust analysis
    """
    # Similar to Pearson but uses:
    rho, p = spearmanr(x[in_window], y[in_window])
```

### Partial Correlation Control

#### Controlling for Confounders
```python
def windowed_partial_correlation(x, y, z, time, window_size):
    """
    Correlation between x and y, controlling for z

    Useful for:
    - Removing age effects
    - Controlling for disease severity
    - Adjusting for technical factors
    """
    # For each window:
    # 1. Regress x on z, get residuals
    # 2. Regress y on z, get residuals
    # 3. Correlate residuals

    x_residuals = x - linear_regression(x, z)
    y_residuals = y - linear_regression(y, z)
    r_partial = pearsonr(x_residuals, y_residuals)
```

---

## 📊 Statistical Testing for Temporal Trends

### Testing for Monotonic Trends

#### Mann-Kendall Test
```python
def mann_kendall_trend_test(correlations, window_centers):
    """
    Test if correlations show significant trend over time

    H0: No monotonic trend
    H1: Monotonic trend exists
    """
    from scipy.stats import kendalltau

    # Kendall's tau between time and correlation
    tau, p_value = kendalltau(window_centers, correlations)

    # Interpretation
    if p_value < 0.05:
        if tau > 0:
            trend = "Increasing"
        else:
            trend = "Decreasing"
    else:
        trend = "No significant trend"

    return tau, p_value, trend
```

### Linear Trend Analysis

#### Regression Approach
```python
def test_linear_trend(correlations, window_centers):
    """
    Linear regression to quantify trend
    """
    from scipy import stats

    # Fit linear model
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        window_centers, correlations
    )

    # Calculate confidence interval for slope
    n = len(window_centers)
    t_crit = stats.t.ppf(0.975, n-2)  # 95% CI
    slope_ci = [slope - t_crit*std_err, slope + t_crit*std_err]

    return {
        'slope': slope,
        'slope_ci': slope_ci,
        'p_value': p_value,
        'r_squared': r_value**2
    }
```

### Permutation Testing

#### Non-parametric Significance
```python
def permutation_test_trend(x, y, time, n_permutations=1000):
    """
    Test if observed trend could occur by chance
    """
    # Calculate observed trend
    observed_corrs, _, centers = windowed_correlation(x, y, time)
    observed_slope, _, _, _, _ = linregress(centers, observed_corrs)

    # Permutation null distribution
    null_slopes = []
    for _ in range(n_permutations):
        # Shuffle time labels
        time_shuffled = np.random.permutation(time)

        # Calculate trend with shuffled time
        perm_corrs, _, centers = windowed_correlation(x, y, time_shuffled)
        perm_slope, _, _, _, _ = linregress(centers, perm_corrs)
        null_slopes.append(perm_slope)

    # P-value: fraction of null slopes more extreme than observed
    p_value = np.mean(np.abs(null_slopes) >= np.abs(observed_slope))

    return p_value, null_slopes
```

---

## 🎭 Bootstrap Confidence Intervals

### Bootstrap for Correlation Trajectories

#### Implementation
```python
def bootstrap_correlation_trajectory(x, y, time, n_bootstrap=1000):
    """
    Bootstrap confidence bands for temporal correlations
    """
    n_samples = len(x)
    bootstrap_trajectories = []

    for _ in range(n_bootstrap):
        # Resample with replacement
        boot_idx = np.random.choice(n_samples, n_samples, replace=True)
        x_boot = x[boot_idx]
        y_boot = y[boot_idx]
        time_boot = time[boot_idx]

        # Calculate correlation trajectory
        corrs, _, centers = windowed_correlation(x_boot, y_boot, time_boot)
        bootstrap_trajectories.append(corrs)

    # Calculate percentile confidence bands
    lower_band = np.percentile(bootstrap_trajectories, 2.5, axis=0)
    upper_band = np.percentile(bootstrap_trajectories, 97.5, axis=0)

    return lower_band, upper_band
```

### Block Bootstrap for Time Series

#### Preserving Temporal Structure
```python
def block_bootstrap_correlation(x, y, time, block_size=10):
    """
    Block bootstrap preserves temporal dependencies
    """
    n_samples = len(x)
    n_blocks = n_samples // block_size

    # Create blocks
    blocks = []
    for i in range(n_blocks):
        start = i * block_size
        end = start + block_size
        blocks.append((x[start:end], y[start:end], time[start:end]))

    # Resample blocks
    boot_blocks = np.random.choice(blocks, n_blocks, replace=True)

    # Reconstruct time series
    x_boot = np.concatenate([b[0] for b in boot_blocks])
    y_boot = np.concatenate([b[1] for b in boot_blocks])
    time_boot = np.concatenate([b[2] for b in boot_blocks])

    return x_boot, y_boot, time_boot
```

---

## 🎯 Change Point Detection

### Finding Critical Transitions

#### PELT Algorithm (Pruned Exact Linear Time)
```python
def detect_change_points(correlations, penalty=3):
    """
    Detect points where correlation changes significantly

    Uses changepoint detection to find:
    - Sudden transitions
    - Disease stage boundaries
    - Critical points
    """
    import ruptures as rpt

    # PELT with normal model
    algo = rpt.Pelt(model="normal", min_size=5).fit(correlations)
    change_points = algo.predict(pen=penalty)

    return change_points
```

### Binary Segmentation Method

#### Recursive Partitioning
```python
def binary_segmentation(data, min_segment=10, threshold=0.1):
    """
    Recursively find change points
    """
    def find_single_change_point(segment):
        best_split = None
        max_difference = 0

        for split in range(min_segment, len(segment)-min_segment):
            left_mean = np.mean(segment[:split])
            right_mean = np.mean(segment[split:])
            difference = abs(left_mean - right_mean)

            if difference > max_difference:
                max_difference = difference
                best_split = split

        return best_split, max_difference

    change_points = []
    segments = [(0, len(data))]

    while segments:
        start, end = segments.pop()
        segment = data[start:end]

        split, difference = find_single_change_point(segment)

        if split and difference > threshold:
            change_points.append(start + split)
            segments.append((start, start + split))
            segments.append((start + split, end))

    return sorted(change_points)
```

---

## 📉 Pseudotime Analysis

### Understanding Pseudotime

#### Concept
```python
# Pseudotime Definition:
"""
Computational ordering of samples along a trajectory
representing biological progression (e.g., disease)

Properties:
- Continuous scale (0 to 1)
- Captures biological progression
- Not necessarily chronological time
- Based on molecular signatures
"""
```

### Trajectory Inference Methods

#### Diffusion Pseudotime
```python
def calculate_diffusion_pseudotime(expression_matrix):
    """
    Order samples along disease trajectory
    """
    import scanpy as sc

    # Create AnnData object
    adata = sc.AnnData(expression_matrix)

    # Preprocessing
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)

    # Dimensionality reduction
    sc.pp.pca(adata)

    # Diffusion map
    sc.pp.neighbors(adata)
    sc.tl.diffmap(adata)

    # Calculate pseudotime
    # Root = earliest disease stage sample
    root_idx = find_healthiest_sample(adata)
    adata.uns['iroot'] = root_idx

    # Diffusion pseudotime
    sc.tl.dpt(adata)

    return adata.obs['dpt_pseudotime']
```

### Validation of Pseudotime

#### Biological Markers
```python
def validate_pseudotime(pseudotime, known_markers):
    """
    Validate pseudotime using known disease markers
    """
    validations = {}

    for marker in known_markers:
        # Correlation with known progression markers
        r, p = pearsonr(pseudotime, marker['expression'])

        expected_direction = marker['expected_trend']
        observed_direction = 'increasing' if r > 0 else 'decreasing'

        validations[marker['name']] = {
            'correlation': r,
            'p_value': p,
            'matches_expectation': expected_direction == observed_direction
        }

    return validations
```

---

## 📊 Multiple Testing Correction

### Challenges in Temporal Analysis

#### Multiple Comparisons Problem
```python
# Sources of Multiple Testing:
"""
1. Multiple windows (e.g., 50 windows)
2. Multiple protein pairs (e.g., 100 pairs)
3. Multiple statistics per window (correlation, p-value)

Total tests = 50 × 100 = 5,000 tests
Expected false positives at α=0.05: 250
"""
```

### Correction Strategies

#### FDR Control for Sliding Windows
```python
def correct_sliding_window_pvalues(p_values_matrix):
    """
    FDR correction accounting for dependence between windows

    p_values_matrix: shape (n_windows, n_protein_pairs)
    """
    from statsmodels.stats.multitest import multipletests

    # Flatten all p-values
    all_pvals = p_values_matrix.flatten()

    # Benjamini-Yekutieli (BY) for dependent tests
    # More conservative than BH but handles dependence
    rejected, corrected_pvals, _, _ = multipletests(
        all_pvals,
        method='fdr_by',
        alpha=0.05
    )

    # Reshape back
    corrected_matrix = corrected_pvals.reshape(p_values_matrix.shape)

    return corrected_matrix, rejected.reshape(p_values_matrix.shape)
```

### Permutation-Based FDR

#### Empirical Null Distribution
```python
def permutation_fdr(observed_pvals, null_pvals_list):
    """
    FDR estimation using permutation null

    More accurate for correlated tests
    """
    n_total = len(observed_pvals)

    # For each p-value threshold
    thresholds = np.linspace(0, 1, 100)
    fdrs = []

    for threshold in thresholds:
        # Observed discoveries
        n_observed = np.sum(observed_pvals <= threshold)

        # Expected false discoveries from null
        null_discoveries = [np.sum(null <= threshold) for null in null_pvals_list]
        expected_false = np.mean(null_discoveries)

        # FDR estimate
        fdr = expected_false / max(n_observed, 1)
        fdrs.append(min(fdr, 1.0))

    return thresholds, fdrs
```

---

## 🎨 Visualization of Temporal Statistics

### Correlation Heatmap Over Time

```python
def plot_temporal_correlation_heatmap(correlation_matrix, window_centers):
    """
    Visualize correlation evolution as heatmap
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(12, 3))

    # Heatmap with time on x-axis
    sns.heatmap(correlation_matrix.T,
                xticklabels=np.round(window_centers, 2),
                yticklabels=protein_pairs,
                cmap='RdBu_r',
                center=0,
                vmin=-1, vmax=1,
                cbar_kws={'label': 'Correlation'})

    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Protein Pairs')
    ax.set_title('Temporal Correlation Dynamics')

    return fig
```

### Statistical Significance Track

```python
def plot_significance_track(p_values, window_centers, alpha=0.05):
    """
    Show statistical significance over time
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # P-values track
    ax1.plot(window_centers, -np.log10(p_values), 'b-', linewidth=2)
    ax1.axhline(-np.log10(alpha), color='r', linestyle='--', label=f'p={alpha}')
    ax1.set_ylabel('-log10(p-value)')
    ax1.set_title('Statistical Significance Over Disease Progression')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Binary significance track
    significant = p_values < alpha
    ax2.fill_between(window_centers, 0, significant,
                     alpha=0.5, color='green', label='Significant')
    ax2.set_ylabel('Significant')
    ax2.set_xlabel('Pseudotime')
    ax2.set_ylim(-0.1, 1.1)
    ax2.legend()

    plt.tight_layout()
    return fig
```

---

## 🔧 Practical Implementation Tips

### Computational Optimization

```python
# Speed up sliding window analysis:
"""
1. Vectorize operations
2. Use numba JIT compilation
3. Parallel processing for protein pairs
4. Cache correlation matrices
5. Use sparse matrices when applicable
"""

from numba import jit
import multiprocessing as mp

@jit(nopython=True)
def fast_correlation(x, y):
    """JIT-compiled correlation calculation"""
    n = len(x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    num = np.sum((x - x_mean) * (y - y_mean))
    den = np.sqrt(np.sum((x - x_mean)**2) * np.sum((y - y_mean)**2))

    return num / den if den != 0 else 0
```

### Quality Control Checks

```python
def quality_control_temporal_analysis(results):
    """
    Verify temporal analysis quality
    """
    checks = {
        'sufficient_windows': len(results['windows']) >= 10,
        'adequate_samples': all(n >= 10 for n in results['sample_counts']),
        'correlation_range': -1 <= min(results['correlations']) and max(results['correlations']) <= 1,
        'no_constant_windows': np.std(results['correlations']) > 0.01,
        'temporal_coverage': (max(results['window_centers']) - min(results['window_centers'])) > 0.5
    }

    return checks
```

---

## 📚 Statistical Best Practices

### Key Recommendations

1. **Window Size**: Start with 30% of range, adjust based on results
2. **Step Size**: Use 10% of window size for smooth trajectories
3. **Minimum Samples**: At least 10-15 per window
4. **Multiple Testing**: Use FDR correction for all p-values
5. **Validation**: Bootstrap confidence intervals for key findings
6. **Robustness**: Compare Pearson and Spearman correlations

### Common Pitfalls to Avoid

- **Too small windows**: Noisy, unreliable correlations
- **Too large windows**: Miss important transitions
- **Ignoring dependencies**: Windows are not independent
- **Over-interpretation**: Correlation ≠ causation
- **Technical artifacts**: Batch effects can create false patterns

---

**Next: [Interpreting Temporal Results](interpreting_results.md)**

*Remember: The statistics are tools to reveal biology - always interpret results in biological context!*
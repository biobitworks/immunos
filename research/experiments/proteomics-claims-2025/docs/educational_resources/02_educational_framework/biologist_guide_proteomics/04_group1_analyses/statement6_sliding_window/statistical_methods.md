# 📊 Statistical Methods for Network Temporal Analysis

## 🎯 Learning Objectives

By the end of this guide, you'll understand:
- ✅ **Network correlation metrics** and their statistical properties
- ✅ **Multiple testing correction** for network-wide analyses
- ✅ **Graph theory statistics** applied to protein networks
- ✅ **Temporal network comparison** methods
- ✅ **Bootstrap approaches** for network confidence
- ✅ **Module detection** and statistical validation

---

## 🌐 Network Correlation Fundamentals

### Correlation Matrix Properties

#### Mathematical Definition
```python
# Correlation Matrix C for n proteins:
"""
C[i,j] = correlation between protein i and protein j
Properties:
- Symmetric: C[i,j] = C[j,i]
- Diagonal = 1: C[i,i] = 1
- Bounded: -1 ≤ C[i,j] ≤ 1
- Positive semi-definite
"""
```

### Statistical Challenges in Network Analysis

#### Dependency Structure
```python
# The multiple testing problem in networks:
"""
For n proteins:
- Number of unique correlations = n(n-1)/2
- Example: 50 proteins = 1,225 correlations
- Tests are NOT independent (transitive correlations)
- Standard FDR may be too conservative
"""
```

---

## 📈 Network-Specific Statistical Tests

### Testing Network Density Changes

```python
def test_network_density_change(density_early, density_late,
                                n_proteins, n_samples_early, n_samples_late):
    """
    Test if network density significantly changes between timepoints

    H0: density_early = density_late
    H1: density_early ≠ density_late
    """
    # Fisher's z-transformation for correlation-based density
    z_early = np.arctanh(density_early)
    z_late = np.arctanh(density_late)

    # Standard error
    se_early = 1 / np.sqrt(n_samples_early - 3)
    se_late = 1 / np.sqrt(n_samples_late - 3)

    # Z-statistic
    z_stat = (z_late - z_early) / np.sqrt(se_early**2 + se_late**2)

    # P-value (two-tailed)
    p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

    return z_stat, p_value
```

### Module Preservation Statistics

```python
def calculate_module_preservation(module_corr_reference, module_corr_test):
    """
    Calculate module preservation statistics between conditions

    Based on Langfelder et al. (2011) PMID: 21118823
    """
    # Z-summary statistic combines multiple preservation measures

    # 1. Density preservation
    density_ref = np.mean(np.abs(module_corr_reference))
    density_test = np.mean(np.abs(module_corr_test))
    z_density = (density_test - density_ref) / np.std(module_corr_reference.flatten())

    # 2. Connectivity preservation
    connectivity_ref = np.sum(np.abs(module_corr_reference), axis=0)
    connectivity_test = np.sum(np.abs(module_corr_test), axis=0)
    r_connectivity = pearsonr(connectivity_ref, connectivity_test)[0]
    z_connectivity = np.arctanh(r_connectivity) * np.sqrt(len(connectivity_ref) - 3)

    # 3. Combine statistics
    z_summary = (z_density + z_connectivity) / 2

    # Interpretation
    preservation = "Strong" if z_summary > 10 else "Moderate" if z_summary > 2 else "Weak"

    return {
        'z_summary': z_summary,
        'z_density': z_density,
        'z_connectivity': z_connectivity,
        'preservation': preservation
    }
```

---

## 🔄 Temporal Network Comparison Methods

### Graphlet-Based Network Comparison

```python
def graphlet_correlation_distance(network1, network2, graphlet_size=3):
    """
    Compare networks using graphlet correlation distance (GCD)

    Graphlets are small subgraphs that capture local topology
    """
    import networkx as nx

    # Count graphlet frequencies
    def count_graphlets(G, size=3):
        graphlet_counts = {}
        for nodes in nx.algorithms.clique.enumerate_all_cliques(G):
            if len(nodes) == size:
                subgraph = G.subgraph(nodes)
                # Convert to canonical form
                canonical = nx.graph_atlas_g()[nx.graph_hashing.weisfeiler_lehman_graph_hash(subgraph)]
                graphlet_counts[canonical] = graphlet_counts.get(canonical, 0) + 1
        return graphlet_counts

    counts1 = count_graphlets(network1, graphlet_size)
    counts2 = count_graphlets(network2, graphlet_size)

    # Normalize and compute correlation
    all_graphlets = set(counts1.keys()) | set(counts2.keys())
    vector1 = [counts1.get(g, 0) for g in all_graphlets]
    vector2 = [counts2.get(g, 0) for g in all_graphlets]

    if sum(vector1) > 0 and sum(vector2) > 0:
        correlation = pearsonr(vector1, vector2)[0]
        gcd = 1 - correlation  # Distance measure
    else:
        gcd = 1.0

    return gcd
```

### Spectral Distance Between Networks

```python
def spectral_distance(corr_matrix1, corr_matrix2):
    """
    Compare networks using eigenvalue distributions

    Captures global network properties
    """
    # Compute eigenvalues
    eigenvalues1 = np.linalg.eigvalsh(corr_matrix1)
    eigenvalues2 = np.linalg.eigvalsh(corr_matrix2)

    # Sort in descending order
    eigenvalues1 = np.sort(eigenvalues1)[::-1]
    eigenvalues2 = np.sort(eigenvalues2)[::-1]

    # Ensure same length
    min_len = min(len(eigenvalues1), len(eigenvalues2))
    eigenvalues1 = eigenvalues1[:min_len]
    eigenvalues2 = eigenvalues2[:min_len]

    # Calculate spectral distance
    spectral_dist = np.sqrt(np.sum((eigenvalues1 - eigenvalues2)**2))

    return spectral_dist
```

---

## 🎲 Bootstrap Methods for Networks

### Network Bootstrap with Correlation Stability

```python
def bootstrap_network_confidence(expression_data, n_bootstrap=1000,
                                 correlation_threshold=0.3):
    """
    Bootstrap confidence intervals for network edges
    """
    n_samples, n_proteins = expression_data.shape
    edge_stability = np.zeros((n_proteins, n_proteins))

    for b in range(n_bootstrap):
        # Resample with replacement
        boot_idx = np.random.choice(n_samples, n_samples, replace=True)
        boot_data = expression_data[boot_idx, :]

        # Calculate correlation matrix
        boot_corr = np.corrcoef(boot_data.T)

        # Count stable edges
        edge_stability += (np.abs(boot_corr) > correlation_threshold).astype(int)

    # Calculate edge confidence
    edge_confidence = edge_stability / n_bootstrap

    return edge_confidence

def identify_stable_edges(edge_confidence, stability_threshold=0.95):
    """
    Identify edges that appear in >95% of bootstraps
    """
    stable_edges = edge_confidence > stability_threshold

    # Set diagonal to False
    np.fill_diagonal(stable_edges, False)

    n_stable = np.sum(stable_edges) / 2  # Divide by 2 for symmetry

    return stable_edges, n_stable
```

### Temporal Stability Assessment

```python
def assess_temporal_stability(sliding_window_results, n_bootstrap=100):
    """
    Bootstrap assessment of temporal pattern stability
    """
    n_windows = len(sliding_window_results['window_centers'])
    n_proteins = sliding_window_results['correlation_matrices'][0].shape[0]

    # Store bootstrap trajectories
    bootstrap_trajectories = []

    for b in range(n_bootstrap):
        # Resample windows with replacement
        boot_indices = np.random.choice(n_windows, n_windows, replace=True)

        # Sort to maintain temporal order
        boot_indices = np.sort(boot_indices)

        # Extract resampled trajectory
        boot_trajectory = sliding_window_results['correlation_matrices'][boot_indices]
        bootstrap_trajectories.append(boot_trajectory)

    # Calculate confidence intervals for each window
    confidence_intervals = []
    for w in range(n_windows):
        window_bootstraps = [traj[w] for traj in bootstrap_trajectories]
        lower = np.percentile(window_bootstraps, 2.5, axis=0)
        upper = np.percentile(window_bootstraps, 97.5, axis=0)
        confidence_intervals.append((lower, upper))

    return confidence_intervals
```

---

## 🎯 Module Detection and Validation

### Statistical Validation of Modules

```python
def validate_network_modules(correlation_matrix, module_assignments,
                            n_permutations=1000):
    """
    Test if modules have significantly higher internal correlation
    than expected by chance
    """
    modules = np.unique(module_assignments)
    module_statistics = {}

    for module in modules:
        # Get module proteins
        module_idx = np.where(module_assignments == module)[0]

        if len(module_idx) < 2:
            continue

        # Calculate observed internal correlation
        module_corr = correlation_matrix[np.ix_(module_idx, module_idx)]
        observed_mean = np.mean(np.abs(module_corr[np.triu_indices_from(module_corr, k=1)]))

        # Permutation test
        null_means = []
        for _ in range(n_permutations):
            # Random proteins of same size
            random_idx = np.random.choice(len(module_assignments),
                                        len(module_idx), replace=False)
            random_corr = correlation_matrix[np.ix_(random_idx, random_idx)]
            null_mean = np.mean(np.abs(random_corr[np.triu_indices_from(random_corr, k=1)]))
            null_means.append(null_mean)

        # Calculate p-value
        p_value = np.mean(null_means >= observed_mean)

        # Effect size (Cohen's d)
        cohens_d = (observed_mean - np.mean(null_means)) / np.std(null_means)

        module_statistics[module] = {
            'observed_correlation': observed_mean,
            'expected_correlation': np.mean(null_means),
            'p_value': p_value,
            'cohens_d': cohens_d,
            'significant': p_value < 0.05
        }

    return module_statistics
```

### Temporal Module Stability

```python
def assess_module_temporal_stability(sliding_window_results, module_assignments):
    """
    Measure how stable modules are across disease progression
    """
    correlation_matrices = sliding_window_results['correlation_matrices']
    window_centers = sliding_window_results['window_centers']

    modules = np.unique(module_assignments)
    stability_scores = {module: [] for module in modules}

    for corr_mat in correlation_matrices:
        for module in modules:
            module_idx = np.where(module_assignments == module)[0]

            if len(module_idx) < 2:
                stability_scores[module].append(np.nan)
                continue

            # Within-module correlation
            within_corr = corr_mat[np.ix_(module_idx, module_idx)]
            within_mean = np.mean(np.abs(within_corr[np.triu_indices_from(within_corr, k=1)]))

            # Between-module correlation
            other_idx = np.where(module_assignments != module)[0]
            if len(other_idx) > 0:
                between_corr = corr_mat[np.ix_(module_idx, other_idx)]
                between_mean = np.mean(np.abs(between_corr))

                # Modularity score
                modularity = within_mean - between_mean
            else:
                modularity = within_mean

            stability_scores[module].append(modularity)

    # Calculate stability metrics
    stability_metrics = {}
    for module, scores in stability_scores.items():
        valid_scores = [s for s in scores if not np.isnan(s)]
        if valid_scores:
            stability_metrics[module] = {
                'mean_modularity': np.mean(valid_scores),
                'std_modularity': np.std(valid_scores),
                'trend': np.corrcoef(range(len(valid_scores)), valid_scores)[0, 1],
                'cv': np.std(valid_scores) / np.mean(valid_scores) if np.mean(valid_scores) > 0 else np.inf
            }

    return stability_metrics
```

---

## 📊 Multiple Testing Correction for Networks

### Network-Specific FDR Control

```python
def network_fdr_correction(p_value_matrix, method='BY'):
    """
    FDR correction accounting for correlation structure in networks

    Methods:
    - 'BY': Benjamini-Yekutieli (handles dependence)
    - 'BH': Benjamini-Hochberg (assumes independence)
    - 'spectral': Spectral method using eigenvalues
    """
    n_proteins = p_value_matrix.shape[0]

    if method in ['BY', 'BH']:
        # Extract upper triangle
        upper_triangle_idx = np.triu_indices(n_proteins, k=1)
        p_values = p_value_matrix[upper_triangle_idx]

        # Apply correction
        from statsmodels.stats.multitest import multipletests
        _, corrected_p, _, _ = multipletests(p_values, method=method.lower())

        # Reconstruct matrix
        corrected_matrix = np.ones_like(p_value_matrix)
        corrected_matrix[upper_triangle_idx] = corrected_p
        corrected_matrix = np.minimum(corrected_matrix, corrected_matrix.T)

    elif method == 'spectral':
        # Estimate effective number of tests using eigenvalues
        eigenvalues = np.linalg.eigvalsh(np.corrcoef(p_value_matrix))
        effective_n = np.sum(eigenvalues > 1)  # Kaiser criterion

        # Bonferroni with effective N
        corrected_matrix = p_value_matrix * effective_n
        corrected_matrix = np.minimum(corrected_matrix, 1.0)

    return corrected_matrix
```

### Permutation-Based FDR for Networks

```python
def permutation_network_fdr(observed_correlations, null_correlations_list,
                           alpha=0.05):
    """
    Estimate FDR using permutation null distribution

    More accurate for network data with complex dependencies
    """
    n_proteins = observed_correlations.shape[0]

    # Calculate p-values for each edge
    p_value_matrix = np.ones((n_proteins, n_proteins))

    for i in range(n_proteins):
        for j in range(i+1, n_proteins):
            obs_corr = observed_correlations[i, j]

            # Null distribution for this edge
            null_dist = [null_corr[i, j] for null_corr in null_correlations_list]

            # Two-tailed p-value
            p_val = np.mean(np.abs(null_dist) >= np.abs(obs_corr))
            p_value_matrix[i, j] = p_val
            p_value_matrix[j, i] = p_val

    # Estimate FDR at different thresholds
    thresholds = np.logspace(-4, 0, 100)
    fdr_estimates = []

    for threshold in thresholds:
        # Observed discoveries
        n_discoveries = np.sum(p_value_matrix < threshold) / 2

        # Expected false discoveries
        null_discoveries = []
        for null_corr in null_correlations_list:
            null_p = calculate_null_pvalues(null_corr, null_correlations_list)
            null_discoveries.append(np.sum(null_p < threshold) / 2)

        expected_false = np.mean(null_discoveries)

        # FDR estimate
        fdr = expected_false / max(n_discoveries, 1)
        fdr_estimates.append(min(fdr, 1.0))

    # Find threshold for target FDR
    fdr_threshold_idx = np.where(np.array(fdr_estimates) <= alpha)[0]
    if len(fdr_threshold_idx) > 0:
        optimal_threshold = thresholds[fdr_threshold_idx[-1]]
    else:
        optimal_threshold = 0

    # Apply threshold
    significant_edges = p_value_matrix < optimal_threshold

    return significant_edges, p_value_matrix, optimal_threshold
```

---

## 🔧 Implementation Best Practices

### Computational Efficiency

```python
# Optimize network calculations:
"""
1. Use sparse matrices for large networks
2. Parallelize correlation calculations
3. Cache intermediate results
4. Use approximate methods for very large networks
5. Consider GPU acceleration for matrix operations
"""

from scipy.sparse import csr_matrix
from joblib import Parallel, delayed

def parallel_correlation_calculation(data, n_jobs=-1):
    """
    Parallelize correlation matrix calculation
    """
    n_proteins = data.shape[1]

    def calculate_row(i, data):
        row = np.zeros(n_proteins)
        for j in range(n_proteins):
            if i <= j:
                row[j] = np.corrcoef(data[:, i], data[:, j])[0, 1]
            else:
                row[j] = row[i]  # Symmetric
        return row

    # Parallel computation
    rows = Parallel(n_jobs=n_jobs)(
        delayed(calculate_row)(i, data) for i in range(n_proteins)
    )

    return np.array(rows)
```

### Quality Control for Network Analysis

```python
def network_quality_control(correlation_matrix, expression_data):
    """
    Quality checks for network analysis
    """
    qc_results = {
        'matrix_symmetric': np.allclose(correlation_matrix, correlation_matrix.T),
        'diagonal_ones': np.allclose(np.diag(correlation_matrix), 1),
        'bounded_values': np.all((-1 <= correlation_matrix) & (correlation_matrix <= 1)),
        'no_nan': not np.any(np.isnan(correlation_matrix)),
        'positive_definite': np.all(np.linalg.eigvalsh(correlation_matrix) >= -1e-10),
        'sufficient_variance': np.all(np.std(expression_data, axis=0) > 0.01)
    }

    return qc_results
```

---

## 📚 Statistical Recommendations

### Key Guidelines

1. **Correlation Threshold**: Use 0.3-0.4 for biological relevance
2. **Multiple Testing**: Apply BY correction for network data
3. **Bootstrap Iterations**: Minimum 1000 for stable estimates
4. **Module Size**: Require ≥5 proteins per module
5. **Temporal Windows**: Ensure ≥15 samples per window
6. **Network Comparison**: Use multiple metrics (density, modularity, spectral)

### Common Pitfalls

- **Ignoring network dependencies** in multiple testing
- **Over-interpreting unstable edges** without bootstrap validation
- **Using inappropriate null models** for significance testing
- **Neglecting batch effects** in network construction
- **Assuming stationarity** in temporal networks

---

**Next: [Interpreting Network Results](interpreting_results.md)**

*Remember: Network statistics must account for the complex dependencies inherent in biological systems!*
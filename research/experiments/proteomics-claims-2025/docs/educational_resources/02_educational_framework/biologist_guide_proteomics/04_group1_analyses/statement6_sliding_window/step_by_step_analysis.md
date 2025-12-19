# 🪟 Advanced Sliding Window Analysis: Multi-Protein Networks

## 🎯 Learning Objectives

By the end of this tutorial, you'll be able to:
- ✅ **Implement sliding window analysis** for multiple protein pairs
- ✅ **Create network correlation matrices** over time
- ✅ **Identify co-regulated protein modules** during disease progression
- ✅ **Perform pathway-level temporal analysis**
- ✅ **Visualize complex temporal networks** effectively

---

## 🌐 Beyond Pairwise: Network-Level Temporal Analysis

### Why Analyze Multiple Proteins Together?

While the SQSTM1-VDAC1 analysis revealed autophagy-mitochondria coupling, real biological systems involve complex networks of interacting proteins.

#### Network Perspective
```python
# Single pair analysis:
"""
SQSTM1 ↔ VDAC1
Limited to one relationship
"""

# Network analysis:
"""
Autophagy Module: SQSTM1, LC3B, ATG5, BECN1
     ↕               ↕           ↕         ↕
Mitochondria: VDAC1, TOMM20, COX4I1, ATP5A1
     ↕               ↕           ↕         ↕
Proteasome: PSMA1, PSMB1, PSMD1, PSMD2

Reveals system-wide coordination
"""
```

---

## 📊 Step 1: Data Preparation for Network Analysis

### Load and Prepare Multi-Protein Dataset

```python
# Import required libraries
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from scipy.cluster import hierarchy
from sklearn.preprocessing import StandardScaler
import networkx as nx
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

print("🪟 Advanced Sliding Window Network Analysis")
print("=" * 50)

# Load the dataset
adata = sc.read_h5ad('pool_processed_v2.h5ad')
print(f"Dataset shape: {adata.shape}")
```

### Define Protein Networks of Interest

```python
# Define protein modules for network analysis
protein_modules = {
    'Autophagy': ['SQSTM1', 'MAP1LC3B', 'ATG5', 'ATG7', 'BECN1', 'ULK1'],
    'Mitochondria': ['VDAC1', 'VDAC2', 'TOMM20', 'COX4I1', 'ATP5F1A', 'NDUFA1'],
    'Proteasome': ['PSMA1', 'PSMB1', 'PSMB2', 'PSMD1', 'PSMD2', 'PSMC1'],
    'Stress_Response': ['HSPA8', 'HSP90AA1', 'HSPB1', 'DNAJB1', 'HSPA5', 'HSPH1'],
    'Synaptic': ['SYN1', 'SYP', 'SNAP25', 'STX1A', 'VAMP2', 'SYT1']
}

# Flatten to get all proteins
all_proteins = []
for module, proteins in protein_modules.items():
    all_proteins.extend(proteins)

print(f"Total proteins to analyze: {len(all_proteins)}")
print(f"Modules: {list(protein_modules.keys())}")
```

### Extract and Validate Protein Data

```python
def extract_protein_network_data(adata, protein_list):
    """
    Extract expression data for multiple proteins
    """
    available_proteins = []
    missing_proteins = []
    expression_data = {}

    for protein in protein_list:
        if protein in adata.var_names:
            idx = adata.var_names.get_loc(protein)
            expr = adata.X[:, idx]
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            expression_data[protein] = expr
            available_proteins.append(protein)
        else:
            # Try case-insensitive search
            alternatives = [p for p in adata.var_names if protein.lower() in p.lower()]
            if alternatives:
                print(f"  {protein} not found, using {alternatives[0]}")
                idx = adata.var_names.get_loc(alternatives[0])
                expr = adata.X[:, idx]
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                expression_data[alternatives[0]] = expr
                available_proteins.append(alternatives[0])
            else:
                missing_proteins.append(protein)

    print(f"\n✅ Found {len(available_proteins)} proteins")
    if missing_proteins:
        print(f"⚠️  Missing {len(missing_proteins)} proteins: {missing_proteins[:5]}...")

    return expression_data, available_proteins

# Extract protein data
protein_data, valid_proteins = extract_protein_network_data(adata, all_proteins)

# Get pseudotime
pseudotime = adata.obs['pseudotime'].values if 'pseudotime' in adata.obs else None
if pseudotime is None:
    print("⚠️  No pseudotime found, using disease score proxy")
    # Create pseudotime from disease markers if not available
    if 'MC1_score' in adata.obs:
        pseudotime = adata.obs['MC1_score'].values
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
    else:
        # Use tau status as simple proxy
        pseudotime = (adata.obs['tau_status'] == 'positive').astype(float)
        # Add noise for continuous values
        pseudotime += np.random.normal(0, 0.1, len(pseudotime))
        pseudotime = np.clip(pseudotime, 0, 1)

print(f"Pseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")
```

---

## 🔄 Step 2: Multi-Protein Sliding Window Implementation

### Network Correlation Matrix Calculation

```python
def calculate_network_correlations_in_window(protein_data, indices):
    """
    Calculate correlation matrix for all protein pairs in a window
    """
    proteins = list(protein_data.keys())
    n_proteins = len(proteins)
    corr_matrix = np.zeros((n_proteins, n_proteins))
    p_matrix = np.ones((n_proteins, n_proteins))

    for i, protein1 in enumerate(proteins):
        for j, protein2 in enumerate(proteins):
            if i <= j:  # Calculate upper triangle only
                expr1 = protein_data[protein1][indices]
                expr2 = protein_data[protein2][indices]

                if len(expr1) > 3:  # Need minimum samples
                    corr, pval = pearsonr(expr1, expr2)
                    corr_matrix[i, j] = corr
                    corr_matrix[j, i] = corr  # Symmetric
                    p_matrix[i, j] = pval
                    p_matrix[j, i] = pval
                else:
                    corr_matrix[i, j] = np.nan
                    corr_matrix[j, i] = np.nan

    return corr_matrix, p_matrix

def sliding_window_network_analysis(protein_data, pseudotime,
                                   window_size=0.3, step_size=0.05,
                                   min_samples=15):
    """
    Perform sliding window analysis on protein network
    """
    print(f"\n🔄 Running sliding window network analysis...")
    print(f"  Window size: {window_size}")
    print(f"  Step size: {step_size}")
    print(f"  Minimum samples: {min_samples}")

    # Initialize storage
    correlation_matrices = []
    p_value_matrices = []
    window_centers = []
    sample_counts = []

    # Sort by pseudotime
    sort_idx = np.argsort(pseudotime)
    pseudotime_sorted = pseudotime[sort_idx]

    # Apply sorting to protein data
    protein_data_sorted = {}
    for protein, expr in protein_data.items():
        protein_data_sorted[protein] = expr[sort_idx]

    # Sliding window
    pt_min, pt_max = pseudotime_sorted.min(), pseudotime_sorted.max()
    current_start = pt_min

    with tqdm(desc="Processing windows") as pbar:
        while current_start + window_size <= pt_max:
            window_end = current_start + window_size
            window_center = current_start + window_size / 2

            # Find samples in window
            in_window = np.where((pseudotime_sorted >= current_start) &
                                (pseudotime_sorted < window_end))[0]

            if len(in_window) >= min_samples:
                # Calculate correlation matrix for this window
                corr_mat, p_mat = calculate_network_correlations_in_window(
                    protein_data_sorted, in_window
                )

                correlation_matrices.append(corr_mat)
                p_value_matrices.append(p_mat)
                window_centers.append(window_center)
                sample_counts.append(len(in_window))

            current_start += step_size
            pbar.update(1)

    print(f"✅ Analyzed {len(window_centers)} windows")

    return {
        'correlation_matrices': np.array(correlation_matrices),
        'p_value_matrices': np.array(p_value_matrices),
        'window_centers': np.array(window_centers),
        'sample_counts': np.array(sample_counts),
        'proteins': list(protein_data.keys())
    }

# Run network analysis
network_results = sliding_window_network_analysis(
    protein_data, pseudotime,
    window_size=0.3, step_size=0.05, min_samples=15
)
```

---

## 📈 Step 3: Extract Key Network Metrics

### Calculate Network Properties Over Time

```python
def calculate_network_metrics(correlation_matrices, threshold=0.3):
    """
    Calculate network properties for each time window
    """
    metrics = {
        'mean_correlation': [],
        'network_density': [],
        'clustering_coefficient': [],
        'modularity': [],
        'n_edges': []
    }

    for corr_mat in correlation_matrices:
        # Mean absolute correlation
        upper_triangle = np.triu(corr_mat, k=1)
        non_diag = upper_triangle[upper_triangle != 0]
        metrics['mean_correlation'].append(np.nanmean(np.abs(non_diag)))

        # Create network from correlation matrix
        G = nx.Graph()
        n_proteins = corr_mat.shape[0]

        # Add edges for correlations above threshold
        edges_added = 0
        for i in range(n_proteins):
            for j in range(i+1, n_proteins):
                if abs(corr_mat[i, j]) > threshold:
                    G.add_edge(i, j, weight=abs(corr_mat[i, j]))
                    edges_added += 1

        metrics['n_edges'].append(edges_added)

        # Network density
        if G.number_of_nodes() > 1:
            metrics['network_density'].append(nx.density(G))

            # Clustering coefficient
            try:
                metrics['clustering_coefficient'].append(
                    nx.average_clustering(G, weight='weight')
                )
            except:
                metrics['clustering_coefficient'].append(0)

            # Modularity (requires community detection)
            try:
                from networkx.algorithms import community
                communities = community.greedy_modularity_communities(G)
                modularity = community.modularity(G, communities)
                metrics['modularity'].append(modularity)
            except:
                metrics['modularity'].append(0)
        else:
            metrics['network_density'].append(0)
            metrics['clustering_coefficient'].append(0)
            metrics['modularity'].append(0)

    return pd.DataFrame(metrics)

# Calculate metrics
network_metrics = calculate_network_metrics(network_results['correlation_matrices'])
network_metrics['window_center'] = network_results['window_centers']

print("\n📊 Network Metrics Summary:")
print(network_metrics.describe())
```

### Identify Module-Specific Patterns

```python
def analyze_module_correlations(network_results, protein_modules):
    """
    Track correlations within and between modules
    """
    proteins = network_results['proteins']
    correlation_matrices = network_results['correlation_matrices']
    window_centers = network_results['window_centers']

    # Map proteins to modules
    protein_to_module = {}
    for module, module_proteins in protein_modules.items():
        for protein in module_proteins:
            if protein in proteins:
                protein_to_module[protein] = module

    # Calculate within and between module correlations
    module_correlations = {
        'window_center': window_centers,
    }

    # Within-module correlations
    for module in protein_modules.keys():
        module_correlations[f'{module}_internal'] = []

    # Between-module correlations
    module_pairs = []
    for i, mod1 in enumerate(protein_modules.keys()):
        for mod2 in list(protein_modules.keys())[i+1:]:
            module_pairs.append(f'{mod1}-{mod2}')
            module_correlations[f'{mod1}-{mod2}'] = []

    # Process each window
    for window_idx, corr_mat in enumerate(correlation_matrices):
        # Within-module
        for module in protein_modules.keys():
            module_indices = [i for i, p in enumerate(proteins)
                            if p in protein_to_module and
                            protein_to_module[p] == module]

            if len(module_indices) > 1:
                module_corrs = []
                for i in module_indices:
                    for j in module_indices:
                        if i < j:
                            module_corrs.append(corr_mat[i, j])

                mean_corr = np.nanmean(module_corrs) if module_corrs else np.nan
                module_correlations[f'{module}_internal'].append(mean_corr)
            else:
                module_correlations[f'{module}_internal'].append(np.nan)

        # Between-module
        for pair in module_pairs:
            mod1, mod2 = pair.split('-')

            mod1_indices = [i for i, p in enumerate(proteins)
                          if p in protein_to_module and
                          protein_to_module[p] == mod1]
            mod2_indices = [i for i, p in enumerate(proteins)
                          if p in protein_to_module and
                          protein_to_module[p] == mod2]

            if mod1_indices and mod2_indices:
                between_corrs = []
                for i in mod1_indices:
                    for j in mod2_indices:
                        between_corrs.append(corr_mat[i, j])

                mean_corr = np.nanmean(between_corrs) if between_corrs else np.nan
                module_correlations[pair].append(mean_corr)
            else:
                module_correlations[pair].append(np.nan)

    return pd.DataFrame(module_correlations)

# Analyze modules
module_results = analyze_module_correlations(network_results, protein_modules)
print("\n🔬 Module Correlation Analysis:")
print(module_results.describe())
```

---

## 🎨 Step 4: Advanced Visualizations

### Create Temporal Network Evolution Plot

```python
def plot_network_evolution(network_results, timepoints=[0.2, 0.5, 0.8],
                          threshold=0.4):
    """
    Visualize network structure at different disease stages
    """
    fig, axes = plt.subplots(1, len(timepoints), figsize=(18, 6))

    proteins = network_results['proteins']
    window_centers = network_results['window_centers']
    correlation_matrices = network_results['correlation_matrices']

    # Color by module
    protein_to_module = {}
    module_colors = {
        'Autophagy': 'red',
        'Mitochondria': 'blue',
        'Proteasome': 'green',
        'Stress_Response': 'orange',
        'Synaptic': 'purple'
    }

    for module, module_proteins in protein_modules.items():
        for protein in module_proteins:
            if protein in proteins:
                protein_to_module[protein] = module

    node_colors = [module_colors.get(protein_to_module.get(p, 'Other'), 'gray')
                   for p in proteins]

    for idx, timepoint in enumerate(timepoints):
        # Find closest window
        closest_idx = np.argmin(np.abs(window_centers - timepoint))
        corr_mat = correlation_matrices[closest_idx]

        # Create network
        G = nx.Graph()
        for i, protein in enumerate(proteins):
            G.add_node(i, label=protein)

        # Add edges for significant correlations
        for i in range(len(proteins)):
            for j in range(i+1, len(proteins)):
                if abs(corr_mat[i, j]) > threshold:
                    G.add_edge(i, j, weight=abs(corr_mat[i, j]))

        # Draw network
        ax = axes[idx]
        pos = nx.spring_layout(G, k=2, seed=42)

        # Draw edges with varying width
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]

        nx.draw_networkx_edges(G, pos, ax=ax,
                              width=[w*3 for w in weights],
                              alpha=0.5)

        # Draw nodes
        nx.draw_networkx_nodes(G, pos, ax=ax,
                              node_color=node_colors,
                              node_size=300,
                              alpha=0.8)

        # Add labels for key proteins
        labels = {i: proteins[i][:5] for i in range(len(proteins))
                 if proteins[i] in ['SQSTM1', 'VDAC1', 'PSMA1']}
        nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=8)

        ax.set_title(f'Pseudotime = {timepoint:.2f}\n' +
                    f'Edges = {G.number_of_edges()}, ' +
                    f'Density = {nx.density(G):.3f}')
        ax.axis('off')

    plt.suptitle('Protein Network Evolution During Disease Progression',
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('network_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("✅ Network evolution plot saved")

# Create network evolution plot
plot_network_evolution(network_results)
```

### Create Module Correlation Heatmap

```python
def plot_module_correlation_heatmap(module_results):
    """
    Visualize module correlations over time as heatmap
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    # Within-module correlations
    within_columns = [col for col in module_results.columns
                     if '_internal' in col]
    within_data = module_results[within_columns].T

    sns.heatmap(within_data,
                xticklabels=np.round(module_results['window_center'], 2),
                yticklabels=[col.replace('_internal', '') for col in within_columns],
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                ax=ax1, cbar_kws={'label': 'Mean Correlation'})

    ax1.set_xlabel('Disease Progression (Pseudotime)')
    ax1.set_ylabel('Module')
    ax1.set_title('Within-Module Correlations Over Time')

    # Between-module correlations
    between_columns = [col for col in module_results.columns
                      if '-' in col and col != 'window_center']
    between_data = module_results[between_columns].T

    sns.heatmap(between_data,
                xticklabels=np.round(module_results['window_center'], 2),
                yticklabels=between_columns,
                cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                ax=ax2, cbar_kws={'label': 'Mean Correlation'})

    ax2.set_xlabel('Disease Progression (Pseudotime)')
    ax2.set_ylabel('Module Pairs')
    ax2.set_title('Between-Module Correlations Over Time')

    plt.tight_layout()
    plt.savefig('module_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("✅ Module correlation heatmap saved")

# Create module heatmap
plot_module_correlation_heatmap(module_results)
```

### Network Metrics Trajectory Plot

```python
def plot_network_metrics_trajectory(network_metrics):
    """
    Plot how network properties change over time
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    metrics_to_plot = [
        ('mean_correlation', 'Mean Absolute Correlation', axes[0, 0]),
        ('network_density', 'Network Density', axes[0, 1]),
        ('clustering_coefficient', 'Clustering Coefficient', axes[1, 0]),
        ('n_edges', 'Number of Edges', axes[1, 1])
    ]

    for metric, title, ax in metrics_to_plot:
        if metric in network_metrics.columns:
            ax.plot(network_metrics['window_center'],
                   network_metrics[metric],
                   'b-', linewidth=2)

            # Add trend line
            z = np.polyfit(network_metrics['window_center'],
                          network_metrics[metric], 1)
            p = np.poly1d(z)
            ax.plot(network_metrics['window_center'],
                   p(network_metrics['window_center']),
                   'r--', alpha=0.5, label=f'Trend (slope={z[0]:.3f})')

            ax.set_xlabel('Disease Progression (Pseudotime)')
            ax.set_ylabel(title)
            ax.set_title(f'{title} Over Disease Progression')
            ax.grid(True, alpha=0.3)
            ax.legend()

    plt.suptitle('Network Properties During Disease Progression',
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('network_metrics_trajectory.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("✅ Network metrics plot saved")

# Plot network metrics
plot_network_metrics_trajectory(network_metrics)
```

---

## 📊 Step 5: Statistical Analysis of Network Changes

### Test for Significant Temporal Trends

```python
def test_network_trends(network_metrics):
    """
    Statistical testing of network metric trends
    """
    from scipy.stats import kendalltau, linregress

    print("\n📈 Statistical Analysis of Network Trends:")
    print("-" * 50)

    results = {}

    for metric in network_metrics.columns:
        if metric != 'window_center':
            # Remove NaN values
            valid_idx = ~network_metrics[metric].isna()
            x = network_metrics.loc[valid_idx, 'window_center']
            y = network_metrics.loc[valid_idx, metric]

            if len(x) > 3:
                # Mann-Kendall trend test
                tau, p_kendall = kendalltau(x, y)

                # Linear regression
                slope, intercept, r_value, p_linear, std_err = linregress(x, y)

                results[metric] = {
                    'kendall_tau': tau,
                    'kendall_p': p_kendall,
                    'linear_slope': slope,
                    'linear_p': p_linear,
                    'r_squared': r_value**2
                }

                # Interpretation
                if p_kendall < 0.05:
                    trend = "INCREASING" if tau > 0 else "DECREASING"
                    print(f"\n{metric}:")
                    print(f"  Trend: {trend} (τ={tau:.3f}, p={p_kendall:.3e})")
                    print(f"  Linear slope: {slope:.3f} per pseudotime unit")
                    print(f"  R²: {r_value**2:.3f}")

    return pd.DataFrame(results).T

# Test trends
trend_results = test_network_trends(network_metrics)
```

### Identify Critical Network Transitions

```python
def find_network_transitions(network_metrics, metric='mean_correlation',
                            threshold_change=0.1):
    """
    Identify points where network properties change dramatically
    """
    values = network_metrics[metric].values
    window_centers = network_metrics['window_center'].values

    # Calculate rate of change
    changes = np.diff(values)
    change_points = []

    for i, change in enumerate(changes):
        if abs(change) > threshold_change:
            change_points.append({
                'pseudotime': window_centers[i+1],
                'change': change,
                'before': values[i],
                'after': values[i+1]
            })

    if change_points:
        print(f"\n🎯 Critical Transitions in {metric}:")
        for i, point in enumerate(change_points):
            direction = "increased" if point['change'] > 0 else "decreased"
            print(f"  Transition {i+1} at pseudotime {point['pseudotime']:.3f}:")
            print(f"    {metric} {direction} by {abs(point['change']):.3f}")
            print(f"    ({point['before']:.3f} → {point['after']:.3f})")
    else:
        print(f"No major transitions detected in {metric}")

    return change_points

# Find transitions
transitions = find_network_transitions(network_metrics)
```

---

## 💾 Step 6: Export Results

### Save Analysis Results

```python
def save_network_analysis_results(network_results, network_metrics,
                                  module_results, filename_prefix='network_sliding_window'):
    """
    Save all analysis results for future use
    """
    import pickle

    # Create results package
    results_package = {
        'network_results': network_results,
        'network_metrics': network_metrics,
        'module_results': module_results,
        'protein_modules': protein_modules,
        'analysis_parameters': {
            'window_size': 0.3,
            'step_size': 0.05,
            'min_samples': 15
        },
        'metadata': {
            'date': pd.Timestamp.now().isoformat(),
            'n_proteins': len(network_results['proteins']),
            'n_windows': len(network_results['window_centers']),
            'n_samples': len(pseudotime)
        }
    }

    # Save as pickle
    with open(f'{filename_prefix}_results.pkl', 'wb') as f:
        pickle.dump(results_package, f)
    print(f"✅ Complete results saved to {filename_prefix}_results.pkl")

    # Save key dataframes as CSV
    network_metrics.to_csv(f'{filename_prefix}_metrics.csv', index=False)
    module_results.to_csv(f'{filename_prefix}_modules.csv', index=False)
    print(f"✅ CSV files saved")

    # Save protein list
    with open(f'{filename_prefix}_proteins.txt', 'w') as f:
        for protein in network_results['proteins']:
            f.write(f"{protein}\n")
    print(f"✅ Protein list saved")

# Save all results
save_network_analysis_results(network_results, network_metrics, module_results)
```

---

## 🎯 Key Findings and Interpretation

### Network-Level Insights

Based on your analysis, you should observe:

1. **Progressive Network Rigidification**
   - Early disease: Sparse, flexible network
   - Late disease: Dense, rigid network
   - Critical transition around pseudotime 0.5

2. **Module-Specific Patterns**
   - Autophagy module: Early correlation increase
   - Mitochondrial module: Mid-stage coupling
   - Synaptic module: Late-stage disruption

3. **Cross-Module Communication**
   - Autophagy-Mitochondria: Progressive coupling
   - Proteasome-Stress Response: Coordinated activation
   - Synaptic-Everything: Late disconnection

### Biological Interpretation

Your network analysis reveals:
- **System-wide coordination** emerges during disease
- **Loss of modularity** as pathways become interdependent
- **Critical points** where multiple systems couple
- **Therapeutic opportunities** before network rigidification

---

## 🚀 Next Steps

### Immediate Extensions
1. Test different correlation thresholds
2. Compare network metrics between tau+ and tau- subgroups
3. Integrate with pathway databases
4. Perform community detection analysis

### Advanced Analyses
1. Dynamic network biomarkers
2. Machine learning on network features
3. Causal network inference
4. Perturbation simulations

---

**Congratulations! You've completed advanced network-level temporal analysis, revealing how entire protein systems reorganize during neurodegeneration.**

*Next: [Statistical Methods for Network Analysis](statistical_methods.md)*

*Remember: Networks capture the complexity of biology - your analysis shows how diseases transform from local to system-wide failures.*
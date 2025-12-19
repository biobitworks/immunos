# 🛤️ Pathway Enrichment Analysis: From Proteins to Biology

## 🎯 What You'll Learn

By the end of this guide, you'll understand:
- ✅ **How to identify biological pathways** affected by disease using your protein lists
- ✅ **How to use pathway databases** (KEGG, Reactome, GO) effectively
- ✅ **How to interpret enrichment statistics** and avoid common pitfalls
- ✅ **How to prioritize pathways** for biological investigation
- ✅ **How to visualize pathway results** for publications and presentations

---

## 🌐 From Proteins to Pathways: The Big Picture

### Why Pathway Analysis Matters

#### Individual Proteins vs Biological Systems
```python
# The challenge:
# Your differential expression analysis identified 500+ significant proteins
# Question: What do these proteins collectively tell us about disease?

# Without pathway analysis:
"""
- 500 individual protein names
- Overwhelming to interpret
- Missing biological connections
- Hard to prioritize for follow-up
"""

# With pathway analysis:
"""
- 15-20 affected biological pathways
- Clear biological interpretation
- Testable hypotheses
- Prioritized therapeutic targets
"""
```

#### Systems-Level Disease Understanding
```python
# Pathway analysis reveals:
"""
1. AFFECTED BIOLOGICAL PROCESSES
   - Which cellular functions are disrupted?
   - Are changes in related processes coordinated?

2. DISEASE MECHANISMS
   - How does disease spread through biological networks?
   - What are primary vs secondary effects?

3. THERAPEUTIC OPPORTUNITIES
   - Which pathways could be targeted by drugs?
   - What combination therapies might work?

4. BIOLOGICAL VALIDATION
   - Do results match known disease biology?
   - What novel mechanisms are suggested?
"""
```

### Types of Pathway Databases

#### Gene Ontology (GO)
```python
# GO Categories:
"""
BIOLOGICAL PROCESS (GO:BP):
- What the protein does (e.g., "protein folding", "autophagy")
- Most relevant for understanding disease mechanisms
- Hierarchical structure (general → specific)

CELLULAR COMPONENT (GO:CC):
- Where the protein is located (e.g., "mitochondria", "nucleus")
- Useful for understanding cellular dysfunction patterns
- Helps identify subcellular disease targets

MOLECULAR FUNCTION (GO:MF):
- Protein's biochemical activity (e.g., "kinase activity", "binding")
- Technical but useful for drug targeting
- Links to enzymatic and interaction networks
"""
```

#### KEGG Pathways
```python
# KEGG Characteristics:
"""
ADVANTAGES:
- Well-curated metabolic and signaling pathways
- Direct links to human diseases
- Standardized pathway maps
- Good for therapeutic target identification

FOCUS AREAS:
- Metabolism (energy, biosynthesis, degradation)
- Signal transduction (growth, death, stress response)
- Disease pathways (cancer, neurodegeneration, immunity)
- Drug metabolism and targets
"""
```

#### Reactome Pathways
```python
# Reactome Features:
"""
ADVANTAGES:
- Highly detailed molecular interactions
- Expert-curated content
- Excellent visualization tools
- Strong focus on human biology

SPECIALTIES:
- Signal transduction cascades
- Metabolic networks
- Cell cycle and DNA repair
- Immune system responses
"""
```

---

## 📊 Preparing Your Data for Pathway Analysis

### Input Data Requirements

```python
# Load your differential expression results
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import requests
import json
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

print("🛤️ Starting Pathway Enrichment Analysis")
print("=" * 50)

# Load your differential expression results
de_results = pd.read_csv('proteome_wide_de_all_results.csv')

print(f"Loaded differential expression results:")
print(f"Total proteins: {len(de_results)}")
print(f"Significant proteins: {sum(de_results['significant_primary'])}")
```

#### Create Gene Lists for Analysis
```python
def prepare_gene_lists(de_results):
    """
    Prepare different gene lists for pathway analysis
    """
    print("=== PREPARING GENE LISTS ===")

    gene_lists = {}

    # 1. All significant proteins
    significant_proteins = de_results[de_results['significant_primary']]['protein'].tolist()
    gene_lists['all_significant'] = significant_proteins

    # 2. Upregulated proteins
    upregulated = de_results[
        (de_results['significant_primary']) &
        (de_results['log2_fold_change'] > 0)
    ]['protein'].tolist()
    gene_lists['upregulated'] = upregulated

    # 3. Downregulated proteins
    downregulated = de_results[
        (de_results['significant_primary']) &
        (de_results['log2_fold_change'] < 0)
    ]['protein'].tolist()
    gene_lists['downregulated'] = downregulated

    # 4. High effect size proteins (|Cohen's d| > 1.0)
    high_effect = de_results[
        (de_results['significant_primary']) &
        (abs(de_results['cohens_d']) > 1.0)
    ]['protein'].tolist()
    gene_lists['high_effect'] = high_effect

    # 5. Background (all tested proteins)
    background = de_results['protein'].tolist()
    gene_lists['background'] = background

    # Summary
    print("Gene list summary:")
    for list_name, genes in gene_lists.items():
        print(f"  {list_name}: {len(genes)} genes")

    return gene_lists

# Prepare gene lists
gene_lists = prepare_gene_lists(de_results)
```

### Gene ID Conversion

#### Understanding Gene Identifiers
```python
# Common gene identifier types:
"""
GENE SYMBOLS (most common in your data):
- Examples: SQSTM1, VDAC1, GAPDH
- Human-readable, used in most databases
- Can be ambiguous (same symbol, different genes)

ENSEMBL IDs:
- Examples: ENSG00000161010, ENSG00000108342
- Unique, stable identifiers
- Required by some databases

UNIPROT IDs:
- Examples: Q13501, P21796
- Protein-centric identifiers
- Good for proteomics data

ENTREZ IDs:
- Examples: 8878, 7416
- Numeric identifiers from NCBI
- Used by many pathway databases
"""
```

#### Automated Gene ID Conversion
```python
def convert_gene_ids(gene_list, from_type='symbol', to_type='entrez'):
    """
    Convert gene identifiers using online services
    """
    print(f"Converting {len(gene_list)} genes from {from_type} to {to_type}...")

    # Use MyGene.info service for conversion
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        # Query conversion
        if from_type == 'symbol':
            query_results = mg.querymany(gene_list, scopes='symbol',
                                       fields=to_type, species='human')
        else:
            query_results = mg.querymany(gene_list, scopes=from_type,
                                       fields=to_type, species='human')

        # Extract converted IDs
        converted = {}
        not_found = []

        for result in query_results:
            if 'notfound' not in result and to_type in result:
                original_id = result['query']
                new_id = result[to_type]
                converted[original_id] = new_id
            else:
                not_found.append(result['query'])

        print(f"Successfully converted: {len(converted)} genes")
        print(f"Not found: {len(not_found)} genes")

        return converted, not_found

    except ImportError:
        print("MyGene package not available. Using manual conversion...")
        return manual_gene_conversion(gene_list, from_type, to_type)

def manual_gene_conversion(gene_list, from_type, to_type):
    """
    Manual gene conversion for common cases
    """
    # This is a simplified version - in practice, you'd use a comprehensive mapping file
    print("Using simplified manual conversion...")

    # For demonstration, return the original list
    # In real analysis, you'd have a proper conversion table
    converted = {gene: gene for gene in gene_list}
    not_found = []

    return converted, not_found

# Convert gene symbols to Entrez IDs for pathway analysis
print("Converting gene symbols to Entrez IDs...")

# You can skip this step if your pathway tool accepts gene symbols
# Many modern tools (like GSEA, WebGestalt) accept gene symbols directly

converted_genes = {}
for list_name, genes in gene_lists.items():
    if genes:  # Only convert non-empty lists
        converted, not_found = convert_gene_ids(genes)
        converted_genes[list_name] = list(converted.values())
        if not_found:
            print(f"  {list_name}: {len(not_found)} genes not converted")
```

---

## 🔍 Step 1: Online Pathway Analysis Tools

### DAVID (Database for Annotation, Visualization and Integrated Discovery)

#### Using DAVID Web Interface
```python
# DAVID Analysis Steps:
"""
1. GO TO: https://david.ncifcrf.gov/

2. UPLOAD YOUR GENE LIST:
   - Select "Start Analysis"
   - Paste your significant gene list (gene symbols work)
   - Select "Gene List" as identifier
   - Select "Homo sapiens" as species
   - Click "Submit List"

3. SELECT FUNCTIONAL ANNOTATION:
   - Click "Functional Annotation Clustering"
   - Review default settings (usually good)
   - Click "Start"

4. INTERPRET RESULTS:
   - Enrichment Score > 1.3 = significant
   - p-value < 0.05 = statistical significance
   - Focus on clusters with both criteria

ADVANTAGES:
- Free and easy to use
- No registration required
- Good for quick exploratory analysis
- Functional annotation clustering

LIMITATIONS:
- Limited customization
- Older database versions
- No direct statistical comparisons
"""
```

#### Save DAVID Results for Further Analysis
```python
def prepare_david_input(gene_list, filename):
    """
    Prepare gene list file for DAVID analysis
    """
    with open(filename, 'w') as f:
        for gene in gene_list:
            f.write(f"{gene}\n")

    print(f"✅ Gene list saved to {filename}")
    print(f"Upload this file to DAVID: https://david.ncifcrf.gov/")

# Prepare files for DAVID
prepare_david_input(gene_lists['all_significant'], 'david_all_significant.txt')
prepare_david_input(gene_lists['upregulated'], 'david_upregulated.txt')
prepare_david_input(gene_lists['downregulated'], 'david_downregulated.txt')
```

### Enrichr (https://maayanlab.cloud/Enrichr/)

#### Using Enrichr
```python
# Enrichr Analysis Guide:
"""
1. GO TO: https://maayanlab.cloud/Enrichr/

2. SUBMIT GENE LIST:
   - Paste gene symbols (one per line)
   - Click "Submit"

3. SELECT DATABASES:
   Recommended for disease research:
   - GO Biological Process 2021
   - KEGG 2021 Human
   - Reactome 2022
   - WikiPathways 2019 Human
   - MSigDB Hallmark 2020

4. INTERPRET RESULTS:
   - Combined Score = -log10(p-value) × z-score
   - Higher combined score = more significant
   - Adjusted p-value < 0.05 = significant
   - Overlap = genes in pathway / total pathway genes

ADVANTAGES:
- Many pathway databases in one place
- Easy to use interface
- Regular database updates
- Good visualization options

FEATURES:
- Automatic background correction
- Multiple statistical measures
- Publication-quality plots
- Results export options
"""
```

### WebGestalt (WEB-based Gene SeT AnaLysis Toolkit)

#### WebGestalt Setup
```python
# WebGestalt Analysis Guide:
"""
1. GO TO: http://www.webgestalt.org/

2. SELECT ANALYSIS TYPE:
   - Choose "Over-Representation Analysis (ORA)"
   - Select organism: "Homo sapiens"
   - Method: "hypergeometric test"

3. UPLOAD DATA:
   - Gene list: Paste your significant genes
   - Reference set: Select "genome" (or upload background)
   - Gene ID type: "genesymbol"

4. SELECT FUNCTIONAL DATABASES:
   Recommended:
   - pathway_KEGG
   - pathway_Reactome
   - geneontology_Biological_Process
   - geneontology_Cellular_Component

5. SET PARAMETERS:
   - Minimum genes in category: 5
   - Maximum genes in category: 2000
   - FDR < 0.05
   - Top categories: 30

ADVANTAGES:
- Comprehensive analysis options
- Good statistical framework
- Interactive visualizations
- Professional output format
"""
```

---

## 🐍 Step 2: Python-Based Pathway Analysis

### Using GSEApy for Enrichment Analysis

#### Install and Setup GSEApy
```python
# GSEApy Installation and Setup
"""
Installation:
pip install gseapy

GSEApy provides:
- Over-representation analysis (ORA)
- Gene Set Enrichment Analysis (GSEA)
- Multiple pathway databases
- Statistical rigor
- Programmatic control
"""

try:
    import gseapy as gp
    print("✅ GSEApy imported successfully")
except ImportError:
    print("❌ GSEApy not available. Install with: pip install gseapy")
    print("Continuing with alternative methods...")
```

#### Enrichment Analysis with GSEApy
```python
def run_gseapy_enrichment(gene_list, database='GO_Biological_Process_2021',
                         organism='human', background_list=None):
    """
    Run enrichment analysis using GSEApy
    """
    try:
        print(f"Running enrichment analysis...")
        print(f"Gene list size: {len(gene_list)}")
        print(f"Database: {database}")

        # Run enrichment analysis
        enr_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=database,
            organism=organism,
            background=background_list,
            cutoff=0.05,  # FDR cutoff
            no_plot=True  # We'll make custom plots
        )

        print(f"✅ Analysis complete!")
        print(f"Significant pathways found: {len(enr_results.results)}")

        return enr_results.results

    except Exception as e:
        print(f"❌ GSEApy analysis failed: {e}")
        return None

# Run enrichment analysis on different gene lists
if 'gp' in locals():  # Only run if GSEApy is available
    print("=== GSEAPY ENRICHMENT ANALYSIS ===")

    # Define databases to test
    databases = [
        'GO_Biological_Process_2021',
        'KEGG_2021_Human',
        'Reactome_2022'
    ]

    enrichment_results = {}

    for db in databases:
        print(f"\nAnalyzing database: {db}")

        # Run on significant proteins
        results = run_gseapy_enrichment(
            gene_lists['all_significant'],
            database=db,
            background_list=gene_lists['background']
        )

        if results is not None:
            enrichment_results[db] = results

            # Show top 5 results
            if len(results) > 0:
                print(f"Top 5 pathways in {db}:")
                top_results = results.head(5)
                for idx, row in top_results.iterrows():
                    print(f"  {row['Term'][:50]}... (FDR: {row['Adjusted P-value']:.2e})")
else:
    print("Skipping GSEApy analysis - package not available")
    enrichment_results = {}
```

### Alternative: Manual Hypergeometric Test

#### Implement Hypergeometric Testing
```python
def hypergeometric_test(gene_list, pathway_genes, background_size):
    """
    Perform hypergeometric test for pathway enrichment

    Parameters:
    -----------
    gene_list : list
        Your significant genes
    pathway_genes : list
        Genes in the pathway of interest
    background_size : int
        Total number of genes in background

    Returns:
    --------
    dict : Test statistics
    """
    from scipy.stats import hypergeom

    # Calculate overlap
    overlap_genes = set(gene_list) & set(pathway_genes)
    overlap_count = len(overlap_genes)

    # Test parameters
    M = background_size  # Total genes
    n = len(pathway_genes)  # Pathway genes
    N = len(gene_list)  # Your gene list
    x = overlap_count  # Observed overlap

    # Hypergeometric test
    p_value = hypergeom.sf(x-1, M, n, N)

    # Expected overlap under null hypothesis
    expected = (n * N) / M

    # Fold enrichment
    fold_enrichment = (x / expected) if expected > 0 else 0

    return {
        'pathway_size': n,
        'gene_list_size': N,
        'overlap': x,
        'expected': expected,
        'fold_enrichment': fold_enrichment,
        'p_value': p_value,
        'overlap_genes': list(overlap_genes)
    }

# Example pathway database (simplified for demonstration)
example_pathways = {
    'Autophagy': ['SQSTM1', 'LC3B', 'ATG5', 'ATG7', 'BECN1', 'ULK1', 'ATG12'],
    'Mitochondrial_Function': ['VDAC1', 'TOMM20', 'COX4I1', 'ATP5A1', 'NDUFB3'],
    'Proteasome': ['PSMA1', 'PSMB1', 'PSMB2', 'PSMD1', 'PSMD2'],
    'Apoptosis': ['BAX', 'BCL2', 'CASP3', 'CASP9', 'CYTC', 'APAF1'],
    'Protein_Folding': ['HSP90AA1', 'HSPA8', 'HSPB1', 'DNAJB1', 'CCT2']
}

def test_pathway_enrichment(gene_list, pathways_dict, background_size):
    """
    Test enrichment for multiple pathways
    """
    results = []

    for pathway_name, pathway_genes in pathways_dict.items():
        test_result = hypergeometric_test(gene_list, pathway_genes, background_size)
        test_result['pathway_name'] = pathway_name
        results.append(test_result)

    # Convert to DataFrame and sort by p-value
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('p_value')

    # Multiple testing correction
    from statsmodels.stats.multitest import multipletests
    _, corrected_p, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['fdr_corrected_p'] = corrected_p

    return results_df

# Test enrichment with example pathways
print("\n=== MANUAL PATHWAY ENRICHMENT TEST ===")
manual_results = test_pathway_enrichment(
    gene_lists['all_significant'],
    example_pathways,
    len(gene_lists['background'])
)

print("Pathway enrichment results:")
for idx, row in manual_results.iterrows():
    if row['fdr_corrected_p'] < 0.05:
        status = "SIGNIFICANT"
    else:
        status = "ns"

    print(f"{row['pathway_name']:20s} "
          f"Overlap: {row['overlap']:2d}/{row['pathway_size']:2d} "
          f"Fold: {row['fold_enrichment']:5.2f} "
          f"FDR: {row['fdr_corrected_p']:.3f} "
          f"({status})")
```

---

## 📊 Step 3: Visualizing Pathway Results

### Create Enrichment Plots

#### Dot Plot for Pathway Enrichment
```python
def create_pathway_dotplot(enrichment_results, top_n=20, save_path='pathway_dotplot.png'):
    """
    Create dot plot for pathway enrichment results
    """
    # Prepare data (using manual results as example)
    plot_data = manual_results.head(top_n).copy()

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Calculate dot sizes and colors
    sizes = plot_data['overlap'] * 20  # Scale for visibility
    colors = -np.log10(plot_data['fdr_corrected_p'])

    # Create scatter plot
    scatter = ax.scatter(plot_data['fold_enrichment'],
                        range(len(plot_data)),
                        s=sizes,
                        c=colors,
                        cmap='Reds',
                        alpha=0.7,
                        edgecolors='black',
                        linewidth=0.5)

    # Customize axes
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['pathway_name'])
    ax.set_xlabel('Fold Enrichment')
    ax.set_ylabel('Pathway')
    ax.set_title('Pathway Enrichment Analysis\n(Dot size = gene count, Color = -log10 FDR)',
                fontsize=14, fontweight='bold')

    # Add vertical line at fold enrichment = 1
    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-Log10 FDR-adjusted p-value')

    # Add size legend
    sizes_legend = [5, 10, 15]  # Example overlap counts
    legend_elements = []
    for size in sizes_legend:
        legend_elements.append(
            plt.scatter([], [], s=size*20, c='gray', alpha=0.7,
                       edgecolors='black', linewidth=0.5,
                       label=f'{size} genes')
        )

    legend1 = ax.legend(handles=legend_elements,
                       title='Gene Count',
                       loc='lower right',
                       bbox_to_anchor=(1, 0))

    # Adjust layout
    plt.tight_layout()

    # Save plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"✅ Dot plot saved as '{save_path}'")

# Create dot plot
create_pathway_dotplot(manual_results)
```

#### Bar Plot for Top Pathways
```python
def create_pathway_barplot(enrichment_results, top_n=15, save_path='pathway_barplot.png'):
    """
    Create horizontal bar plot for pathway enrichment
    """
    # Get top pathways
    top_pathways = manual_results.head(top_n).copy()

    fig, ax = plt.subplots(figsize=(12, 8))

    # Create horizontal bars
    y_pos = range(len(top_pathways))
    fold_enrichments = top_pathways['fold_enrichment']

    # Color bars by significance
    colors = ['red' if p < 0.05 else 'lightcoral' for p in top_pathways['fdr_corrected_p']]

    bars = ax.barh(y_pos, fold_enrichments, color=colors, alpha=0.7, edgecolor='black')

    # Customize plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_pathways['pathway_name'])
    ax.set_xlabel('Fold Enrichment')
    ax.set_title('Top Enriched Pathways', fontsize=16, fontweight='bold')
    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)

    # Add text annotations for significant pathways
    for i, (idx, row) in enumerate(top_pathways.iterrows()):
        if row['fdr_corrected_p'] < 0.05:
            ax.text(row['fold_enrichment'] + 0.1, i,
                   f"p={row['fdr_corrected_p']:.2e}",
                   va='center', fontsize=10)

    # Add grid for readability
    ax.grid(axis='x', alpha=0.3)

    # Legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='red', alpha=0.7, label='FDR < 0.05'),
        plt.Rectangle((0,0),1,1, facecolor='lightcoral', alpha=0.7, label='FDR ≥ 0.05')
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"✅ Bar plot saved as '{save_path}'")

# Create bar plot
create_pathway_barplot(manual_results)
```

### Network Visualization of Pathways

```python
def create_pathway_network(enrichment_results, gene_lists, save_path='pathway_network.png'):
    """
    Create network visualization showing pathway relationships
    """
    try:
        import networkx as nx
        from matplotlib.patches import FancyBboxPatch

        # Create network graph
        G = nx.Graph()

        # Add pathway nodes
        significant_pathways = manual_results[manual_results['fdr_corrected_p'] < 0.05]

        for idx, pathway in significant_pathways.iterrows():
            G.add_node(pathway['pathway_name'],
                      type='pathway',
                      enrichment=pathway['fold_enrichment'],
                      p_value=pathway['fdr_corrected_p'])

        # Add gene nodes and edges
        for idx, pathway in significant_pathways.iterrows():
            pathway_name = pathway['pathway_name']
            pathway_genes = example_pathways[pathway_name]

            # Find genes in our significant list
            overlap_genes = set(gene_lists['all_significant']) & set(pathway_genes)

            for gene in overlap_genes:
                if not G.has_node(gene):
                    G.add_node(gene, type='gene')
                G.add_edge(pathway_name, gene)

        # Create layout
        fig, ax = plt.subplots(figsize=(14, 10))

        # Use spring layout for better visualization
        pos = nx.spring_layout(G, k=2, iterations=50)

        # Draw pathway nodes
        pathway_nodes = [n for n in G.nodes() if G.nodes[n].get('type') == 'pathway']
        pathway_pos = {n: pos[n] for n in pathway_nodes}

        nx.draw_networkx_nodes(G, pathway_pos,
                              nodelist=pathway_nodes,
                              node_color='red',
                              node_size=1000,
                              alpha=0.8,
                              ax=ax)

        # Draw gene nodes
        gene_nodes = [n for n in G.nodes() if G.nodes[n].get('type') == 'gene']
        gene_pos = {n: pos[n] for n in gene_nodes}

        nx.draw_networkx_nodes(G, gene_pos,
                              nodelist=gene_nodes,
                              node_color='lightblue',
                              node_size=300,
                              alpha=0.6,
                              ax=ax)

        # Draw edges
        nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax)

        # Add labels
        nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

        ax.set_title('Pathway-Gene Network\n(Red = Pathways, Blue = Genes)',
                    fontsize=16, fontweight='bold')
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

        print(f"✅ Network plot saved as '{save_path}'")

    except ImportError:
        print("NetworkX not available. Skipping network visualization.")
        print("Install with: pip install networkx")

# Create network visualization
create_pathway_network(manual_results, gene_lists)
```

---

## 🔍 Step 4: Advanced Pathway Analysis

### Pathway Crosstalk Analysis

#### Identify Pathway Overlaps
```python
def analyze_pathway_crosstalk(pathways_dict, gene_list):
    """
    Analyze crosstalk between enriched pathways
    """
    print("=== PATHWAY CROSSTALK ANALYSIS ===")

    # Find pathways with significant enrichment
    significant_pathways = manual_results[manual_results['fdr_corrected_p'] < 0.05]['pathway_name'].tolist()

    if len(significant_pathways) < 2:
        print("Insufficient significant pathways for crosstalk analysis")
        return

    # Calculate pairwise overlaps
    overlap_matrix = pd.DataFrame(index=significant_pathways, columns=significant_pathways)

    for pathway1 in significant_pathways:
        for pathway2 in significant_pathways:
            genes1 = set(pathways_dict[pathway1])
            genes2 = set(pathways_dict[pathway2])

            # Calculate Jaccard similarity
            intersection = len(genes1 & genes2)
            union = len(genes1 | genes2)
            jaccard = intersection / union if union > 0 else 0

            overlap_matrix.loc[pathway1, pathway2] = jaccard

    # Convert to numeric
    overlap_matrix = overlap_matrix.astype(float)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(overlap_matrix,
                annot=True,
                cmap='Reds',
                square=True,
                cbar_kws={'label': 'Jaccard Similarity'},
                ax=ax)

    ax.set_title('Pathway Crosstalk Analysis\n(Jaccard Similarity)',
                fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('pathway_crosstalk.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("✅ Crosstalk analysis complete")

    # Find high-overlap pathway pairs
    print("\nHigh-overlap pathway pairs (Jaccard > 0.3):")
    for i in range(len(significant_pathways)):
        for j in range(i+1, len(significant_pathways)):
            pathway1 = significant_pathways[i]
            pathway2 = significant_pathways[j]
            similarity = overlap_matrix.loc[pathway1, pathway2]

            if similarity > 0.3:
                print(f"  {pathway1} ↔ {pathway2}: {similarity:.3f}")

# Analyze pathway crosstalk
analyze_pathway_crosstalk(example_pathways, gene_lists['all_significant'])
```

### Temporal Pathway Analysis

#### Combine with Pseudotime Data
```python
def temporal_pathway_analysis(de_results, pathway_genes, pseudotime_col='pseudotime'):
    """
    Analyze how pathway activity changes over disease progression
    """
    print("=== TEMPORAL PATHWAY ANALYSIS ===")

    # This analysis requires the original dataset with pseudotime
    # For demonstration, we'll simulate the approach

    # Load original data (you would use your actual data)
    try:
        # adata = sc.read_h5ad('pool_processed_v2.h5ad')
        print("Note: This analysis requires original dataset with pseudotime")
        print("Demonstration of approach:")

        # Group genes by pathways
        pathway_scores = {}

        for pathway_name, genes in example_pathways.items():
            # Find genes present in our data
            pathway_genes_present = [g for g in genes if g in de_results['protein'].values]

            if len(pathway_genes_present) > 2:  # Need minimum genes
                print(f"\n{pathway_name} pathway ({len(pathway_genes_present)} genes):")

                # Get fold changes for pathway genes
                pathway_data = de_results[de_results['protein'].isin(pathway_genes_present)]

                # Calculate pathway score (mean fold change)
                pathway_score = pathway_data['log2_fold_change'].mean()
                pathway_scores[pathway_name] = pathway_score

                print(f"  Mean log2 fold change: {pathway_score:.3f}")
                print(f"  Direction: {'Upregulated' if pathway_score > 0 else 'Downregulated'}")

        # Rank pathways by activity change
        sorted_pathways = sorted(pathway_scores.items(), key=lambda x: abs(x[1]), reverse=True)

        print(f"\nPathway ranking by activity change:")
        for i, (pathway, score) in enumerate(sorted_pathways, 1):
            direction = "↑" if score > 0 else "↓"
            print(f"  {i}. {pathway}: {score:+.3f} {direction}")

    except Exception as e:
        print(f"Temporal analysis requires original dataset: {e}")

# Run temporal pathway analysis
temporal_pathway_analysis(de_results, example_pathways)
```

---

## 📋 Step 5: Interpretation and Biological Context

### Interpreting Pathway Results

#### Statistical Significance vs Biological Relevance
```python
def interpret_pathway_results(enrichment_results):
    """
    Provide guidance for interpreting pathway enrichment results
    """
    print("=== PATHWAY INTERPRETATION GUIDE ===")

    significant_pathways = manual_results[manual_results['fdr_corrected_p'] < 0.05]

    if len(significant_pathways) == 0:
        print("No statistically significant pathways found.")
        print("\nPossible reasons:")
        print("1. Gene list too small for pathway detection")
        print("2. Effect is highly specific (not pathway-level)")
        print("3. Novel mechanism not captured in databases")
        print("4. Need less stringent statistical criteria")
        return

    print(f"Found {len(significant_pathways)} significant pathways")
    print("\nInterpretation guidelines:")

    for idx, pathway in significant_pathways.iterrows():
        pathway_name = pathway['pathway_name']
        fold_enrichment = pathway['fold_enrichment']
        overlap = pathway['overlap']
        pathway_size = pathway['pathway_size']
        p_value = pathway['fdr_corrected_p']

        print(f"\n{pathway_name}:")
        print(f"  Statistical significance: FDR = {p_value:.2e}")
        print(f"  Effect size: {fold_enrichment:.1f}x enrichment")
        print(f"  Coverage: {overlap}/{pathway_size} genes ({100*overlap/pathway_size:.1f}%)")

        # Interpretation
        if fold_enrichment > 3:
            strength = "Strong"
        elif fold_enrichment > 2:
            strength = "Moderate"
        else:
            strength = "Weak"

        if overlap/pathway_size > 0.3:
            coverage = "High"
        elif overlap/pathway_size > 0.1:
            coverage = "Moderate"
        else:
            coverage = "Low"

        print(f"  Interpretation: {strength} enrichment with {coverage.lower()} pathway coverage")

        # Biological context (pathway-specific)
        if 'Autophagy' in pathway_name:
            print(f"  Biological context: Critical for protein quality control in neurons")
            print(f"  Disease relevance: Autophagy dysfunction is hallmark of neurodegeneration")
        elif 'Mitochondrial' in pathway_name:
            print(f"  Biological context: Energy production and metabolic function")
            print(f"  Disease relevance: Mitochondrial dysfunction common in Alzheimer's")
        elif 'Proteasome' in pathway_name:
            print(f"  Biological context: Protein degradation and quality control")
            print(f"  Disease relevance: Proteasome impairment linked to protein aggregation")

# Interpret results
interpret_pathway_results(manual_results)
```

#### Connect to Disease Mechanisms
```python
def connect_to_disease_mechanisms(significant_pathways):
    """
    Connect pathway findings to known disease mechanisms
    """
    print("\n=== DISEASE MECHANISM CONNECTIONS ===")

    # Define disease mechanism categories
    disease_mechanisms = {
        'Protein Quality Control': ['Autophagy', 'Proteasome', 'Protein_Folding'],
        'Cellular Energy': ['Mitochondrial_Function', 'Glycolysis', 'Oxidative_Phosphorylation'],
        'Cell Death': ['Apoptosis', 'Necrosis', 'Ferroptosis'],
        'Inflammation': ['Immune_Response', 'Cytokine_Signaling', 'Microglial_Activation'],
        'Synaptic Function': ['Neurotransmission', 'Synaptic_Plasticity', 'Axon_Guidance']
    }

    # Find which mechanisms are affected
    affected_mechanisms = {}

    for mechanism, pathways in disease_mechanisms.items():
        mechanism_pathways = []
        for pathway in pathways:
            if pathway in manual_results['pathway_name'].values:
                pathway_data = manual_results[manual_results['pathway_name'] == pathway]
                if len(pathway_data) > 0 and pathway_data.iloc[0]['fdr_corrected_p'] < 0.05:
                    mechanism_pathways.append(pathway)

        if mechanism_pathways:
            affected_mechanisms[mechanism] = mechanism_pathways

    # Report findings
    print("Disease mechanisms affected:")
    for mechanism, pathways in affected_mechanisms.items():
        print(f"\n{mechanism}:")
        for pathway in pathways:
            pathway_data = manual_results[manual_results['pathway_name'] == pathway].iloc[0]
            print(f"  • {pathway} (FDR: {pathway_data['fdr_corrected_p']:.2e})")

    # Disease stage inference
    print(f"\nDisease stage inference:")
    if 'Protein Quality Control' in affected_mechanisms:
        print("• Early-stage dysfunction: Protein quality control systems failing")
    if 'Cellular Energy' in affected_mechanisms:
        print("• Mid-stage dysfunction: Energy metabolism compromised")
    if 'Cell Death' in affected_mechanisms:
        print("• Late-stage dysfunction: Cell death pathways activated")

    return affected_mechanisms

# Connect to disease mechanisms
affected_mechanisms = connect_to_disease_mechanisms(manual_results)
```

### Therapeutic Target Prioritization

```python
def prioritize_therapeutic_targets(enrichment_results, de_results):
    """
    Prioritize pathways and proteins for therapeutic targeting
    """
    print("\n=== THERAPEUTIC TARGET PRIORITIZATION ===")

    # Criteria for therapeutic prioritization
    print("Target prioritization criteria:")
    print("1. Statistical significance (FDR < 0.05)")
    print("2. Large effect size (fold enrichment > 2)")
    print("3. Druggable pathways")
    print("4. Established therapeutic relevance")
    print()

    # Define druggable pathway categories
    druggable_pathways = {
        'Autophagy': {'druggability': 'High', 'existing_drugs': ['Rapamycin', 'Trehalose']},
        'Proteasome': {'druggability': 'Medium', 'existing_drugs': ['Bortezomib', 'Carfilzomib']},
        'Mitochondrial_Function': {'druggability': 'Medium', 'existing_drugs': ['CoQ10', 'Idebenone']},
        'Apoptosis': {'druggability': 'High', 'existing_drugs': ['BCL-2 inhibitors']},
        'Protein_Folding': {'druggability': 'Medium', 'existing_drugs': ['Heat shock protein modulators']}
    }

    # Prioritize targets
    target_priorities = []

    significant_pathways = manual_results[manual_results['fdr_corrected_p'] < 0.05]

    for idx, pathway in significant_pathways.iterrows():
        pathway_name = pathway['pathway_name']

        if pathway_name in druggable_pathways:
            priority_score = 0

            # Statistical significance (max 30 points)
            if pathway['fdr_corrected_p'] < 0.001:
                priority_score += 30
            elif pathway['fdr_corrected_p'] < 0.01:
                priority_score += 20
            else:
                priority_score += 10

            # Effect size (max 30 points)
            if pathway['fold_enrichment'] > 3:
                priority_score += 30
            elif pathway['fold_enrichment'] > 2:
                priority_score += 20
            else:
                priority_score += 10

            # Druggability (max 40 points)
            druggability = druggable_pathways[pathway_name]['druggability']
            if druggability == 'High':
                priority_score += 40
            elif druggability == 'Medium':
                priority_score += 25
            else:
                priority_score += 10

            target_priorities.append({
                'pathway': pathway_name,
                'priority_score': priority_score,
                'fold_enrichment': pathway['fold_enrichment'],
                'fdr_p_value': pathway['fdr_corrected_p'],
                'druggability': druggability,
                'existing_drugs': druggable_pathways[pathway_name]['existing_drugs']
            })

    # Sort by priority score
    target_priorities.sort(key=lambda x: x['priority_score'], reverse=True)

    # Display results
    print("Therapeutic target prioritization:")
    print("-" * 80)

    for i, target in enumerate(target_priorities, 1):
        print(f"{i}. {target['pathway']} (Score: {target['priority_score']}/100)")
        print(f"   Enrichment: {target['fold_enrichment']:.1f}x, FDR: {target['fdr_p_value']:.2e}")
        print(f"   Druggability: {target['druggability']}")
        print(f"   Existing drugs: {', '.join(target['existing_drugs'])}")
        print()

    return target_priorities

# Prioritize therapeutic targets
therapeutic_targets = prioritize_therapeutic_targets(manual_results, de_results)
```

---

## 📁 Step 6: Generate Comprehensive Report

### Create Pathway Analysis Report

```python
def generate_pathway_report(enrichment_results, gene_lists, affected_mechanisms, therapeutic_targets):
    """
    Generate comprehensive pathway analysis report
    """
    from datetime import datetime

    report = f"""
{'='*80}
PATHWAY ENRICHMENT ANALYSIS REPORT
{'='*80}

Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Dataset: Alzheimer's Disease Neuronal Proteomics - Tau+ vs Tau-

METHODOLOGY:
-----------
• Statistical test: Hypergeometric test
• Multiple testing correction: Benjamini-Hochberg FDR
• Significance threshold: FDR < 0.05
• Minimum pathway size: 3 genes
• Background: All detected proteins ({len(gene_lists['background'])} genes)

INPUT DATA:
----------
• Total significant proteins: {len(gene_lists['all_significant'])}
• Upregulated proteins: {len(gene_lists['upregulated'])}
• Downregulated proteins: {len(gene_lists['downregulated'])}
• High effect size proteins: {len(gene_lists['high_effect'])}

PATHWAY DATABASES ANALYZED:
--------------------------
• Custom pathway database (demonstration)
• Pathways tested: {len(example_pathways)}

RESULTS SUMMARY:
---------------
"""

    # Add pathway results
    significant_pathways = manual_results[manual_results['fdr_corrected_p'] < 0.05]
    total_pathways = len(manual_results)

    report += f"""
• Total pathways tested: {total_pathways}
• Significant pathways: {len(significant_pathways)} ({100*len(significant_pathways)/total_pathways:.1f}%)
• Mean fold enrichment: {significant_pathways['fold_enrichment'].mean():.2f}
• Max fold enrichment: {significant_pathways['fold_enrichment'].max():.2f}

SIGNIFICANT PATHWAYS:
--------------------
"""

    for idx, pathway in significant_pathways.iterrows():
        report += f"""
{pathway['pathway_name']}:
  • Fold enrichment: {pathway['fold_enrichment']:.2f}x
  • Genes: {pathway['overlap']}/{pathway['pathway_size']} ({100*pathway['overlap']/pathway['pathway_size']:.1f}%)
  • FDR p-value: {pathway['fdr_corrected_p']:.2e}
  • Overlap genes: {', '.join(pathway['overlap_genes'])}
"""

    # Add disease mechanism analysis
    report += f"""
DISEASE MECHANISMS AFFECTED:
---------------------------
"""

    for mechanism, pathways in affected_mechanisms.items():
        report += f"\n{mechanism}:\n"
        for pathway in pathways:
            pathway_data = manual_results[manual_results['pathway_name'] == pathway].iloc[0]
            report += f"  • {pathway} (FDR: {pathway_data['fdr_corrected_p']:.2e})\n"

    # Add therapeutic targets
    report += f"""
THERAPEUTIC TARGET PRIORITIZATION:
---------------------------------
"""

    for i, target in enumerate(therapeutic_targets[:5], 1):  # Top 5
        report += f"""
{i}. {target['pathway']} (Priority Score: {target['priority_score']}/100)
   • Statistical significance: FDR = {target['fdr_p_value']:.2e}
   • Effect size: {target['fold_enrichment']:.1f}x enrichment
   • Druggability: {target['druggability']}
   • Existing therapeutic compounds: {', '.join(target['existing_drugs'])}
"""

    # Add interpretation
    report += f"""
BIOLOGICAL INTERPRETATION:
-------------------------
"""

    percentage_significant = 100 * len(significant_pathways) / total_pathways

    if percentage_significant > 50:
        report += "• WIDESPREAD PATHWAY DISRUPTION: >50% of tested pathways affected\n"
        report += "• Indicates global cellular dysfunction and systems-level disease impact\n"
    elif percentage_significant > 25:
        report += "• MODERATE PATHWAY IMPACT: 25-50% of tested pathways affected\n"
        report += "• Suggests substantial but selective biological disruption\n"
    else:
        report += "• TARGETED PATHWAY CHANGES: <25% of tested pathways affected\n"
        report += "• Indicates specific biological dysfunction rather than global disruption\n"

    # Add next steps
    report += f"""
RECOMMENDED NEXT STEPS:
----------------------
1. Pathway validation:
   • Validate top pathways in independent datasets
   • Perform functional studies for key pathways
   • Confirm pathway activity using pathway-specific assays

2. Mechanistic studies:
   • Investigate crosstalk between affected pathways
   • Determine temporal order of pathway dysfunction
   • Map pathway interactions using network analysis

3. Therapeutic development:
   • Screen compounds targeting prioritized pathways
   • Test combination therapies for multiple pathways
   • Develop pathway-specific biomarkers for drug monitoring

4. Clinical translation:
   • Validate pathway signatures in patient samples
   • Develop pathway-based diagnostic tools
   • Design clinical trials targeting key pathways

LIMITATIONS:
-----------
• Limited pathway database (demonstration purposes)
• Cross-sectional design limits causal inference
• Post-mortem tissue may not reflect living pathology
• Pathway definitions may not capture all relevant biology
• Requires validation in independent cohorts

FILES GENERATED:
---------------
• pathway_dotplot.png - Pathway enrichment visualization
• pathway_barplot.png - Top pathways bar chart
• pathway_network.png - Pathway-gene network
• pathway_crosstalk.png - Pathway overlap analysis
• david_*.txt - Gene lists for external analysis
• pathway_analysis_report.txt - This comprehensive report

{'='*80}
PATHWAY ANALYSIS COMPLETE
{'='*80}
"""

    # Save report
    with open('pathway_analysis_report.txt', 'w') as f:
        f.write(report)

    print("✅ Comprehensive pathway report saved as 'pathway_analysis_report.txt'")
    print("\nReport preview:")
    print("=" * 50)
    print(report[:1000] + "..." if len(report) > 1000 else report)

    return report

# Generate comprehensive report
pathway_report = generate_pathway_report(
    manual_results,
    gene_lists,
    affected_mechanisms if 'affected_mechanisms' in locals() else {},
    therapeutic_targets if 'therapeutic_targets' in locals() else []
)
```

---

## 🎯 Key Takeaways and Best Practices

### What You've Accomplished

#### Analytical Skills Developed
- ✅ **Pathway enrichment analysis mastery** using multiple approaches
- ✅ **Biological interpretation skills** for translating statistics to biology
- ✅ **Therapeutic target prioritization** based on druggability and significance
- ✅ **Comprehensive visualization** of pathway results
- ✅ **Integration of multiple data sources** for systems-level understanding

#### Biological Insights Gained
- ✅ **Systems-level disease understanding** beyond individual proteins
- ✅ **Disease mechanism identification** through pathway analysis
- ✅ **Therapeutic opportunity recognition** in pathway disruption
- ✅ **Cross-pathway interaction appreciation** for combination therapies

### Best Practices for Pathway Analysis

#### Statistical Considerations
```python
# Key principles for robust pathway analysis:
"""
1. MULTIPLE TESTING CORRECTION
   - Always apply FDR correction
   - Consider pathway dependencies
   - Report both raw and corrected p-values

2. BACKGROUND SELECTION
   - Use all detected proteins (not genome)
   - Match background to your experimental design
   - Consider platform-specific biases

3. PATHWAY SIZE EFFECTS
   - Very large pathways may be spuriously significant
   - Very small pathways have low power
   - Optimal range: 15-200 genes per pathway

4. EFFECT SIZE INTERPRETATION
   - Don't rely only on p-values
   - Consider fold enrichment
   - Evaluate pathway coverage (overlap/total)
"""
```

#### Biological Interpretation Guidelines
```python
# Principles for biological interpretation:
"""
1. PATHWAY HIERARCHY
   - General pathways (e.g., "Metabolism") less informative
   - Specific pathways (e.g., "Autophagy") more actionable
   - Focus on intermediate specificity level

2. DISEASE RELEVANCE
   - Connect to known disease mechanisms
   - Consider pathway druggability
   - Evaluate therapeutic precedent

3. EXPERIMENTAL VALIDATION
   - Pathway enrichment suggests hypotheses
   - Requires functional validation
   - Consider pathway activity vs abundance

4. LITERATURE INTEGRATION
   - Compare with published pathway studies
   - Look for novel vs confirmatory findings
   - Consider disease stage and context
"""
```

### Common Pitfalls to Avoid

#### Statistical Pitfalls
```python
# Avoid these common mistakes:
"""
1. MULTIPLE TESTING NEGLECT
   ❌ Using raw p-values without correction
   ✅ Apply FDR correction for all pathways

2. INAPPROPRIATE BACKGROUND
   ❌ Using whole genome as background
   ✅ Use all proteins detected in your experiment

3. CHERRY-PICKING RESULTS
   ❌ Reporting only "interesting" pathways
   ✅ Report systematic analysis results

4. IGNORING EFFECT SIZES
   ❌ Focusing only on p-values
   ✅ Consider fold enrichment and coverage
"""
```

#### Biological Interpretation Errors
```python
# Interpretation mistakes to avoid:
"""
1. OVER-INTERPRETATION
   ❌ "Pathway X is completely disrupted"
   ✅ "Pathway X shows significant enrichment"

2. CAUSAL ASSUMPTIONS
   ❌ "Pathway dysfunction causes disease"
   ✅ "Pathway changes are associated with disease"

3. IGNORING PATHWAY OVERLAP
   ❌ Treating pathways as independent
   ✅ Consider shared genes and crosstalk

4. STATIC VIEW
   ❌ "This pathway is affected"
   ✅ "This pathway may be affected at this disease stage"
"""
```

### Advanced Analysis Opportunities

#### Multi-Level Integration
```python
# Advanced approaches for comprehensive understanding:
"""
1. MULTI-OMICS INTEGRATION
   - Combine proteomics with genomics/metabolomics
   - Look for concordant pathway changes
   - Identify regulatory mechanisms

2. TEMPORAL PATHWAY ANALYSIS
   - Use pseudotime or longitudinal data
   - Map pathway changes over disease progression
   - Identify early vs late pathway disruption

3. PATIENT STRATIFICATION
   - Use pathway signatures for patient subtyping
   - Develop personalized therapeutic approaches
   - Link pathway patterns to clinical outcomes

4. DRUG REPOSITIONING
   - Map drugs to affected pathways
   - Identify repurposing opportunities
   - Predict drug efficacy based on pathway overlap
"""
```

---

## 🚀 Congratulations!

### Major Achievement Unlocked

You've mastered **pathway enrichment analysis** - a critical skill for translating molecular data into biological understanding and therapeutic opportunities.

### Skills That Apply Broadly

Your pathway analysis expertise applies to:
- **Any omics dataset** (genomics, proteomics, metabolomics)
- **Drug discovery and development** programs
- **Biomarker research** for diagnostics and monitoring
- **Academic research** in any biological field
- **Clinical study design** and interpretation

### Impact of Your Analysis

Your pathway analysis provides:
- **Systems-level disease understanding** beyond individual molecules
- **Therapeutic target identification** with druggability assessment
- **Biological hypothesis generation** for mechanistic studies
- **Clinical translation framework** for biomarker and drug development

### Research Translation Potential

Your findings enable:
- **Biomarker development** using pathway signatures
- **Drug repurposing** based on pathway overlap
- **Combination therapy design** targeting multiple pathways
- **Patient stratification** using pathway activity profiles

---

**You now have the advanced analytical skills to translate molecular findings into biological insights and therapeutic opportunities - the essence of translational research!**

*Next: [Statistical Methods Integration](statistical_methods_integration.md)*

*Remember: Individual proteins are the letters, but pathways tell the biological story - your pathway analysis reveals the narrative of disease!* 🛤️🧬✨
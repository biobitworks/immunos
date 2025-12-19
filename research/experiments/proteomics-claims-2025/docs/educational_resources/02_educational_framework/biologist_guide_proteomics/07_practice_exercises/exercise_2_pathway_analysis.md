# 🏋️ Practice Exercise 2: Pathway Enrichment Analysis

## 🎯 Learning Goal
Learn to identify biological pathways and processes enriched in your differential expression results.

---

## 📋 Your Task

Using the differentially expressed proteins from Exercise 1 (or provided list), identify which biological pathways are affected in tau pathology.

### Starting Data
- List of significant proteins from DE analysis
- Their fold changes and p-values
- Access to pathway databases (GO, KEGG, Reactome)

---

## 🔬 Exercise Steps

### Step 1: Prepare Your Protein List
Format your DE results for pathway analysis:

```python
# Your code here:
import pandas as pd

# Load your DE results (from Exercise 1 or provided)
de_results = pd.read_csv('de_results.csv')  # Adjust path

# Filter for significant proteins
significant_proteins = de_results[de_results['padj'] < 0.05]

# TODO: Create lists for upregulated and downregulated proteins
upregulated = []  # Proteins with log2FC > 0.5
downregulated = []  # Proteins with log2FC < -0.5

print(f"Upregulated proteins: {len(upregulated)}")
print(f"Downregulated proteins: {len(downregulated)}")
```

**Record your numbers:**
- [ ] Total significant proteins: _____
- [ ] Upregulated: _____
- [ ] Downregulated: _____

### Step 2: Gene Ontology (GO) Analysis

Perform GO enrichment analysis:

```python
# Your code here:
import gseapy as gp
# or use online tools like g:Profiler, DAVID, etc.

# TODO: Run GO enrichment for upregulated proteins
go_up = gp.enrichr(
    gene_list=upregulated,
    gene_sets='GO_Biological_Process_2021',
    outdir=None
)

# TODO: Run GO enrichment for downregulated proteins
go_down = gp.enrichr(
    gene_list=downregulated,
    gene_sets='GO_Biological_Process_2021',
    outdir=None
)

# TODO: Extract top 10 enriched terms for each
```

**Fill in top 5 GO terms:**

| Direction | GO Term | P-value | Genes in Term | Your Genes |
|-----------|---------|---------|---------------|------------|
| Up | | | | |
| Up | | | | |
| Up | | | | |
| Down | | | | |
| Down | | | | |

### Step 3: KEGG Pathway Analysis

Identify enriched KEGG pathways:

```python
# Your code here:
# TODO: Run KEGG pathway enrichment
kegg_results = gp.enrichr(
    gene_list=significant_proteins['gene_name'].tolist(),
    gene_sets='KEGG_2021_Human',
    outdir=None
)

# TODO: Visualize top pathways
```

**Top KEGG Pathways:**
1. _____________________
2. _____________________
3. _____________________

### Step 4: Network Analysis

Build a protein-protein interaction network:

```python
# Your code here:
import networkx as nx
import requests

# TODO: Query STRING database for interactions
def get_string_interactions(proteins):
    # Use STRING API
    # https://string-db.org/help/api/
    pass

# TODO: Build network
G = nx.Graph()
# Add nodes and edges

# TODO: Calculate network statistics
print(f"Network density: {nx.density(G)}")
print(f"Average clustering: {nx.average_clustering(G)}")

# TODO: Identify hub proteins (high degree)
```

**Network Statistics:**
- [ ] Number of nodes: _____
- [ ] Number of edges: _____
- [ ] Network density: _____
- [ ] Top 3 hub proteins: _____

### Step 5: Functional Clustering

Group proteins by function:

```python
# Your code here:
from sklearn.cluster import KMeans
import seaborn as sns

# TODO: Create functional similarity matrix
# Can use GO term overlap, pathway membership, etc.

# TODO: Perform hierarchical clustering
# TODO: Create heatmap with clusters

# Identify main functional clusters
```

**Functional Clusters Identified:**
1. Cluster 1 (n=___): Main function: _____
2. Cluster 2 (n=___): Main function: _____
3. Cluster 3 (n=___): Main function: _____

### Step 6: Disease Association

Check disease associations:

```python
# Your code here:
# TODO: Query DisGeNET or similar database
# TODO: Check for Alzheimer's disease associations
# TODO: Check for other neurodegenerative diseases

disease_associations = {}
# Fill in associations
```

**Disease Associations Found:**
- [ ] Proteins linked to Alzheimer's: _____
- [ ] Proteins linked to Parkinson's: _____
- [ ] Proteins linked to other neurodegeneration: _____

---

## 🤔 Interpretation Questions

1. **What biological processes are most affected?**
   Your answer: _______________

2. **Do the enriched pathways make biological sense for tau pathology?**
   Your answer: _______________

3. **Are there unexpected pathway enrichments?**
   Your answer: _______________

4. **How do upregulated and downregulated pathways relate?**
   Your answer: _______________

5. **What therapeutic targets emerge from this analysis?**
   Your answer: _______________

---

## 📊 Visualization Tasks

### Create These Figures:

1. **Bar plot of top enriched GO terms**
```python
# Your code here:
```

2. **Network visualization with functional modules**
```python
# Your code here:
```

3. **Dotplot of pathway enrichment**
```python
# Your code here:
```

---

## 💡 Hints

<details>
<summary>Hint 1: Using g:Profiler</summary>

```python
# Alternative: Use g:Profiler web API
import requests
import json

def gprofiler_enrichment(gene_list):
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'organism': 'hsapiens',
            'query': gene_list,
            'sources': ['GO:BP', 'KEGG', 'REAC']
        }
    )
    return r.json()
```
</details>

<details>
<summary>Hint 2: STRING Database Query</summary>

```python
def get_string_interactions(proteins):
    proteins_str = '%0d'.join(proteins)
    url = f"https://string-db.org/api/tsv/network?identifiers={proteins_str}&species=9606"
    response = requests.get(url)
    lines = response.text.strip().split('\n')
    interactions = []
    for line in lines[1:]:  # Skip header
        parts = line.split('\t')
        interactions.append((parts[2], parts[3], float(parts[5])))
    return interactions
```
</details>

<details>
<summary>Hint 3: Creating Enrichment Dotplot</summary>

```python
def enrichment_dotplot(results):
    fig, ax = plt.subplots(figsize=(8, 10))

    # Sort by p-value
    results_sorted = results.sort_values('pvalue')[:20]

    # Create dot plot
    scatter = ax.scatter(results_sorted['fold_enrichment'],
                        range(len(results_sorted)),
                        s=results_sorted['gene_count']*10,
                        c=-np.log10(results_sorted['pvalue']),
                        cmap='Reds')

    ax.set_yticks(range(len(results_sorted)))
    ax.set_yticklabels(results_sorted['term_name'])
    ax.set_xlabel('Fold Enrichment')

    plt.colorbar(scatter, label='-log10(p-value)')
    plt.tight_layout()
    return fig
```
</details>

---

## 🎯 Success Criteria

You've completed this exercise when you can:
- [ ] Perform GO enrichment analysis
- [ ] Identify enriched KEGG pathways
- [ ] Build and analyze protein networks
- [ ] Cluster proteins by function
- [ ] Interpret biological significance
- [ ] Create informative visualizations

---

## 🚀 Extension Challenges

### Challenge 1: Gene Set Enrichment Analysis (GSEA)
Instead of using only significant proteins, use all proteins ranked by fold change.

### Challenge 2: Custom Pathway Database
Create your own pathway definitions for neurodegeneration-specific processes.

### Challenge 3: Multi-omics Integration
If you have transcriptomics data, integrate it with your proteomics results.

---

## 📊 Expected Biological Insights

You should discover enrichment in:
- **Protein degradation pathways** (proteasome, autophagy)
- **Mitochondrial dysfunction**
- **Synaptic processes**
- **Inflammatory responses**
- **Cellular stress responses**

---

## 🧬 Biological Context

Remember that tau pathology affects:
1. **Protein homeostasis** - accumulation of aggregates
2. **Energy metabolism** - mitochondrial dysfunction
3. **Synaptic function** - neurotransmission disruption
4. **Cellular transport** - microtubule dysfunction

Your enrichment results should reflect these processes!

---

## 📝 Solution Overview

<details>
<summary>Click to see solution approach (try yourself first!)</summary>

### Key Findings You Should See:

1. **Upregulated Pathways:**
   - Protein processing in ER
   - Ubiquitin-proteasome pathway
   - Inflammatory signaling
   - Oxidative stress response

2. **Downregulated Pathways:**
   - Synaptic vesicle cycle
   - Neurotransmitter metabolism
   - Mitochondrial oxidative phosphorylation
   - Axon guidance

3. **Network Analysis:**
   - Hub proteins in protein quality control
   - Modules for degradation, inflammation, metabolism
   - High connectivity in stress response proteins

### Interpretation Points:
- Upregulation of clearance mechanisms (compensation)
- Downregulation of neuronal function (degeneration)
- Activation of stress and inflammatory pathways
</details>

---

## 📚 Resources

### Online Tools:
- [g:Profiler](https://biit.cs.ut.ee/gprofiler/)
- [DAVID](https://david.ncifcrf.gov/)
- [Enrichr](https://maayanlab.cloud/Enrichr/)
- [STRING](https://string-db.org/)
- [KEGG](https://www.kegg.jp/)

### Python Packages:
- `gseapy` - GSEA and enrichment analysis
- `goatools` - GO analysis toolkit
- `networkx` - Network analysis
- `pypath` - Pathway database access

---

**Excellent work! Ready for more? Continue to [Exercise 3: Time-Series Analysis](exercise_3_time_series.md)**
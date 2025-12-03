---
name: "BioGRID"
type: "database"
category: "protein-interactions"
maintained_by: "The BioGRID Consortium"
website: "https://thebiogrid.org/"
repository: "https://github.com/BioGRID"
license: "MIT"
programming_language: "Multiple"
installation: "web-api"
status: "active"
last_release: "2024"
documentation_quality: "excellent"
data_size: "~500 MB (tab-delimited download)"
added_date: "2025-12-03"
last_updated: "2025-12-03"
tags:
  - database
  - protein-interactions
  - genetic-interactions
  - systems-biology
  - aging-biology
---

# BioGRID - Biological General Repository for Interaction Datasets

## Quick Summary

Comprehensive database of protein-protein, genetic, and chemical interactions. Contains 2.9+ million interactions from 87,534+ publications. Free access via web interface, REST API, and bulk downloads. Critical resource for aging gene network analysis.

## Overview

### Purpose and Scope
BioGRID curates physical and genetic interactions, protein post-translational modifications, and chemical associations from major model organisms and humans. Essential for understanding biological networks and pathways relevant to aging.

### Key Features
- **2.9+ million interactions**: Protein, genetic, chemical
- **1.1+ million PTMs**: Post-translational modifications
- **87,534+ publications**: Literature-curated data
- **Multiple organisms**: Human, mouse, yeast, fly, worm, etc.
- **Monthly updates**: Continuously expanding
- **Free access**: Open data, multiple formats
- **REST API**: Programmatic access
- **Web interface**: Interactive exploration

### Version Information
- **Current Version**: 5.0.252 (as of 2024)
- **Update Frequency**: Monthly curation updates
- **Status**: Actively maintained and expanding
- **Coverage**: Major model organisms + humans

## Technical Details

### Data Statistics (Build 5.0.252)

**⚠️ DATA SIZE TRACKING**

| Component | Count | Notes |
|-----------|-------|-------|
| **Protein/genetic interactions** | 2,905,263 | Total interactions |
| **Non-redundant interactions** | 2,203,605 | Unique interactions |
| **Chemical interactions** | 31,540 | Drug-protein associations |
| **PTMs** | 1,128,339 | Post-translational modifications |
| **Publications curated** | 87,534 | Literature sources |
| **Download size** | **~500 MB** | Tab-delimited format |

**Storage recommendation**: Download specific organism/interaction type rather than full database if storage is concern.

### Data Schema

#### Interaction Types
1. **Physical interactions**: Protein-protein binding
2. **Genetic interactions**: Functional relationships
3. **Chemical associations**: Drug-target interactions
4. **PTMs**: Phosphorylation, ubiquitination, etc.

#### Key Fields
- **Interactor A/B**: Gene/protein identifiers
- **Interaction type**: Physical, genetic, chemical
- **Experimental system**: Method used to detect
- **Publication**: PubMed ID
- **Organism**: Species (Human, Mouse, etc.)
- **Throughput**: Low/high-throughput
- **Score**: Confidence metrics (if available)

### Supported Organisms
- **Humans** (Homo sapiens)
- **Mouse** (Mus musculus)
- **Yeast** (S. cerevisiae)
- **Fly** (D. melanogaster)
- **Worm** (C. elegans)
- **Arabidopsis** (A. thaliana)
- **Zebrafish** (D. rerio)
- **Others** (50+ organisms total)

## Access Methods

### Web Interface
**URL**: https://thebiogrid.org/

**Search features**:
- Gene/protein name search
- Interaction browsing
- Network visualization
- Publication-based queries
- Organism filtering

**No download required** ✅

### REST API

**Base URL**: https://webservice.thebiogrid.org/

**Authentication**: Requires free API key from https://webservice.thebiogrid.org/

**Example query**:
```bash
# Get interactions for TP53
curl "https://webservice.thebiogrid.org/interactions?geneList=TP53&searchNames=true&format=json&accessKey=YOUR_KEY"
```

**Python example**:
```python
import requests

API_KEY = "your_api_key_here"
url = "https://webservice.thebiogrid.org/interactions"

params = {
    "accessKey": API_KEY,
    "format": "json",
    "geneList": "TP53|FOXO3|IGF1R",  # Pipe-separated
    "searchNames": "true",
    "includeInteractors": "true",
    "taxId": 9606  # Human
}

response = requests.get(url, params=params)
interactions = response.json()

print(f"Found {len(interactions)} interactions")
```

### Bulk Download

**Download page**: https://downloads.thebiogrid.org/BioGRID/

**Files available**:
- **BIOGRID-ALL**: All interactions (~500 MB)
- **BIOGRID-ORGANISM**: Species-specific subsets
- **BIOGRID-CHEMICALS**: Chemical interactions only
- **BIOGRID-PTMS**: Post-translational modifications

**Formats**:
- Tab-delimited (tab3)
- PSI-MI XML
- JSON

**Download example**:
```bash
# Download latest human interactions
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip

# Extract
unzip BIOGRID-ORGANISM-Homo_sapiens-*.tab3.zip
```

## Data Access Questions

### Basic Access Question

**Question**: How many protein-protein interaction partners does the growth hormone receptor gene (GHR) have in BioGRID?

**Expected Answer**:
```
12-15 partners (varies by build)
```

**How to answer**:
1. Visit https://thebiogrid.org/
2. Search for "GHR"
3. Select Homo sapiens entry
4. Count unique interaction partners shown
5. Or use API:

```python
import requests

params = {
    "accessKey": "YOUR_KEY",
    "format": "json",
    "geneList": "GHR",
    "searchNames": "true",
    "taxId": 9606
}

response = requests.get("https://webservice.thebiogrid.org/interactions", params=params)
data = response.json()

# Get unique partners
partners = set()
for interaction in data.values():
    if interaction['OFFICIAL_SYMBOL_A'] == 'GHR':
        partners.add(interaction['OFFICIAL_SYMBOL_B'])
    else:
        partners.add(interaction['OFFICIAL_SYMBOL_A'])

print(f"GHR has {len(partners)} interaction partners")
```

### Real-World Question

**Question**: For the top 10 GenAge human longevity genes, identify which ones form a highly-connected interaction network (>5 interactions among themselves) and list the shared interactors that connect them.

**Expected Answer**:
```
Highly connected genes: TP53, IGF1R, FOXO3, SIRT1, AKT1
Shared interactors: MTOR, PIK3CA, IRS1, EGFR
Network density: 0.65
```

**How to answer**:
1. Load top 10 GenAge genes: TP53, IGF1R, FOXO3, SIRT1, APOE, TERT, GHR, WRN, LMNA, AKT1
2. Query BioGRID for interactions between these genes
3. Build interaction network
4. Calculate network density
5. Identify genes with >5 intra-network connections
6. Find shared interaction partners

```python
import requests
import pandas as pd
from collections import Counter

API_KEY = "your_key"
genage_top10 = ["TP53", "IGF1R", "FOXO3", "SIRT1", "APOE",
                "TERT", "GHR", "WRN", "LMNA", "AKT1"]

# Get all interactions for these genes
all_interactions = []
for gene in genage_top10:
    params = {
        "accessKey": API_KEY,
        "format": "json",
        "geneList": gene,
        "searchNames": "true",
        "taxId": 9606,
        "includeInteractors": "true"
    }
    response = requests.get("https://webservice.thebiogrid.org/interactions",
                           params=params)
    all_interactions.extend(response.json().values())

# Build network: interactions within GenAge top 10
network_edges = []
shared_interactors = []

for interaction in all_interactions:
    geneA = interaction['OFFICIAL_SYMBOL_A']
    geneB = interaction['OFFICIAL_SYMBOL_B']

    if geneA in genage_top10 and geneB in genage_top10:
        network_edges.append((geneA, geneB))
    elif geneA in genage_top10:
        shared_interactors.append(geneB)
    elif geneB in genage_top10:
        shared_interactors.append(geneA)

# Count connections per gene
gene_connections = Counter([g for edge in network_edges for g in edge])
highly_connected = [gene for gene, count in gene_connections.items() if count > 5]

# Top shared interactors
common_partners = Counter(shared_interactors).most_common(10)

print(f"Highly connected genes: {', '.join(highly_connected)}")
print(f"Top shared interactors: {', '.join([p[0] for p in common_partners[:4]])}")
print(f"Network density: {len(network_edges) / (len(genage_top10) * (len(genage_top10)-1) / 2):.2f}")
```

## Integration with Our Research

### Relevance to Aging Biology ⭐⭐⭐⭐⭐

**Critical resource** for analyzing aging gene networks:
- Map GenAge gene interaction networks
- Identify CellAge senescence pathway connections
- Find shared interactors among aging genes
- Discover novel aging-related genes via network proximity

### GenAge Integration

**Use case**: Build interaction networks for aging genes

```python
import pandas as pd
import requests

# Load GenAge
genage = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv')

# Query BioGRID for all GenAge genes
aging_network = {}
for symbol in genage['symbol'][:50]:  # Start with top 50
    params = {
        "accessKey": API_KEY,
        "geneList": symbol,
        "searchNames": "true",
        "taxId": 9606,
        "format": "json"
    }
    response = requests.get("https://webservice.thebiogrid.org/interactions",
                           params=params)
    aging_network[symbol] = response.json()

# Analyze network properties
# - Hub genes (high degree)
# - Network modules
# - Pathway enrichment
```

### CellAge Integration

**Use case**: Map senescence gene interaction networks

```python
# Load CellAge
cellage = pd.read_csv('/Users/byron/projects/data/cellage/cellage3.tsv', sep='\t')

# Get senescence-inducing vs inhibiting genes
inducers = cellage[cellage['Senescence Effect'] == 'Induces']['Gene symbol']
inhibitors = cellage[cellage['Senescence Effect'] == 'Inhibits']['Gene symbol']

# Query interactions
# Are senescence inducers connected?
# Do they share common pathways?
# Are there master regulators?
```

### Cross-Reference with Other Databases

**STRING-DB comparison**:
- BioGRID: Experimentally validated interactions
- STRING: Predicted + experimental interactions
- **Strategy**: Use BioGRID for high-confidence, STRING for broader network

**Monarch Initiative integration**:
- BioGRID: Protein interactions
- Monarch: Phenotype associations
- **Combined**: Interaction → Phenotype connections

### Analysis Opportunities

1. **Aging gene network topology**
   - Are aging genes more connected than expected?
   - Network hubs and bottlenecks
   - Pathway enrichment analysis

2. **Senescence interaction patterns**
   - Do senescence genes form modules?
   - Shared regulators of senescence
   - Drug targets in senescence pathways

3. **Cross-species comparison**
   - Conservation of aging gene interactions
   - Model organism validation
   - Evolutionary analysis

4. **Drug target identification**
   - Chemical interactions with aging genes
   - Network-based drug discovery
   - Polypharmacology opportunities

## Usage Examples

### Example 1: Query Single Gene via Web

**Goal**: Find interaction partners of TP53

**Steps**:
1. Visit https://thebiogrid.org/
2. Search "TP53"
3. Select "Homo sapiens"
4. View interaction table
5. Export if needed

**Time**: 2 minutes

### Example 2: API Query for Multiple Genes

**Goal**: Get interactions for top aging genes

```python
import requests
import json

API_KEY = "your_key_here"
aging_genes = ["TP53", "FOXO3", "IGF1R", "SIRT1"]

url = "https://webservice.thebiogrid.org/interactions"

for gene in aging_genes:
    params = {
        "accessKey": API_KEY,
        "format": "json",
        "geneList": gene,
        "searchNames": "true",
        "taxId": 9606
    }

    response = requests.get(url, params=params)
    data = response.json()

    print(f"{gene}: {len(data)} interactions")
```

### Example 3: Build Aging Gene Network

**Goal**: Create network graph of GenAge gene interactions

```python
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load GenAge genes
genage = pd.read_csv('data/genage/human/genage_human.csv')
aging_genes = set(genage['symbol'].tolist())

# Query BioGRID (pseudo-code, requires API key)
# Get all interactions where both genes are in GenAge

G = nx.Graph()

# Add edges from BioGRID results
for interaction in biogrid_results:
    geneA = interaction['OFFICIAL_SYMBOL_A']
    geneB = interaction['OFFICIAL_SYMBOL_B']

    if geneA in aging_genes and geneB in aging_genes:
        G.add_edge(geneA, geneB)

# Analyze network
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Density: {nx.density(G):.3f}")

# Find hubs
degree_centrality = nx.degree_centrality(G)
hubs = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
print(f"Network hubs: {[gene for gene, cent in hubs]}")

# Visualize
nx.draw(G, with_labels=True, node_color='lightblue',
        node_size=500, font_size=8)
plt.savefig('aging_gene_network.png', dpi=300, bbox_inches='tight')
```

## Citation and Attribution

### Primary Citation

**BibTeX**:
```bibtex
@article{oughtred2021biogrid,
  author = {Oughtred, Rose and others},
  title = {The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions},
  journal = {Protein Science},
  volume = {30},
  number = {1},
  pages = {187--200},
  year = {2021},
  doi = {10.1002/pro.3978},
  pmid = {33070389}
}
```

### How to Cite in Publications

```
Data from BioGRID database (version 5.0.252, https://thebiogrid.org/)
(Oughtred et al., 2021)
```

### License
**MIT License** - Free for academic and commercial use

## Performance and Limitations

### Advantages
- **Comprehensive**: 2.9M+ interactions
- **Curated**: Literature-based, not predicted
- **Free**: Open access, no paywall
- **API**: Programmatic access available
- **Updated**: Monthly curation updates
- **Multiple organisms**: Cross-species analysis

### Limitations
- **Literature bias**: Well-studied genes have more interactions
- **Coverage gaps**: Some genes/pathways under-represented
- **Experimental bias**: Certain methods over-represented
- **False positives**: Some interactions may not occur in vivo
- **Dynamic interactions**: Doesn't capture temporal/spatial context
- **API rate limits**: Throttling on bulk queries

### BioGRID vs Alternatives

| Feature | BioGRID | STRING | IntAct | HPRD |
|---------|---------|--------|--------|------|
| **Data source** | Curated | Predicted+Curated | Curated | Curated |
| **Human interactions** | 500K+ | 11M+ | 150K+ | 40K+ |
| **Confidence scores** | Limited | Yes | Yes | Yes |
| **Organisms** | 50+ | 5,000+ | 100+ | Human only |
| **Update frequency** | Monthly | Annually | Quarterly | Discontinued |
| **API** | Yes | Yes | Yes | No |
| **Best for** | High-confidence | Predictions | Quality | Legacy |

**Recommendation**: Use BioGRID for experimentally validated interactions, STRING for broader network context.

## Documentation and Support

### Official Resources
- **Website**: https://thebiogrid.org/
- **Downloads**: https://downloads.thebiogrid.org/
- **API docs**: https://wiki.thebiogrid.org/doku.php/biogridrest
- **Wiki**: https://wiki.thebiogrid.org/
- **Publications**: https://thebiogrid.org/publications

### Getting API Key
1. Visit https://webservice.thebiogrid.org/
2. Register for free account
3. Generate API key
4. Use in API requests

### Support Channels
- **Email**: info@thebiogrid.org
- **Twitter**: @BioGRID
- **GitHub**: https://github.com/BioGRID
- **Publications**: Cite for questions

## Related Resources

### Complementary Databases
- **[[string-db|STRING-DB]]** - Predicted protein interactions
- **[[monarch-kg|Monarch Initiative]]** - Phenotype-genotype networks
- **IntAct** - Molecular interaction database
- **HPRD** - Human protein reference database (legacy)

### Our Databases
- **[[../datasets/genage|GenAge]]** - Query aging gene interactions
- **[[../datasets/cellage|CellAge]]** - Senescence pathway networks

## Next Steps

### Setup
- [ ] Register for BioGRID API key at https://webservice.thebiogrid.org/
- [ ] Test API access with sample query
- [ ] Download human interaction subset (~50 MB)

### Analysis Projects
- [ ] Map GenAge gene interaction network
- [ ] Identify senescence pathway modules in CellAge
- [ ] Find shared interactors among aging genes
- [ ] Network-based drug target discovery

### Integration
- [ ] Create script: `scripts/query_biogrid.py`
- [ ] Integrate with GenAge/CellAge analysis
- [ ] Compare with STRING-DB results
- [ ] Visualize aging gene networks

## Notes

### Strategic Importance
**Very High** for aging biology network analysis:
- Map gene interaction networks
- Identify pathway modules
- Find drug targets
- Validate aging gene connections

### Storage Considerations
- **Full database**: ~500 MB (manageable)
- **Human subset**: ~50 MB (small)
- **API queries**: No local storage needed

**Recommendation**: Start with API queries, download human subset if needed for extensive analysis.

### Best Practices
1. **Use API** for exploratory queries
2. **Download subset** for bulk analysis
3. **Cite properly** in publications
4. **Validate interactions** with literature
5. **Cross-reference** with other databases

### Last Updated
2025-12-03 - Initial documentation with access questions

---

**Created**: 2025-12-03
**Last Modified**: 2025-12-03
**Status**: Active - High priority for network analysis
**API Key**: Required for programmatic access
**Priority**: Very High

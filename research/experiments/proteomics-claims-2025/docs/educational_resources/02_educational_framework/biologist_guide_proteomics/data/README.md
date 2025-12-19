# 📁 Sample Data and Resources

## 🎯 Overview

This directory contains sample datasets and resources for the proteomics analysis exercises. These datasets are designed to be educational while reflecting real-world proteomics analysis challenges.

---

## 📊 Available Datasets

### 1. Main Dataset: `alzheimer_proteomics.h5ad`
**Primary dataset for all exercises**

#### Description
- **Source**: Post-mortem human brain tissue (frontal cortex)
- **Disease**: Alzheimer's disease
- **Technology**: Tandem mass tag (TMT) proteomics
- **Processing**: Log2-transformed, normalized intensities

#### Data Structure
```
Samples (obs): 150 neurons
├── tau_status: 'tau_positive' or 'tau_negative'
├── age: Age at death (years)
├── sex: 'M' or 'F'
├── pmi: Post-mortem interval (hours)
├── batch: Processing batch (1-3)
├── braak_stage: Braak tau pathology stage (0-VI)
├── pseudotime: Disease progression score (0-1)
└── cognitive_score: Pre-mortem cognitive assessment

Proteins (var): 5,853 proteins
├── gene_name: Official gene symbol
├── protein_name: Protein description
├── uniprot_id: UniProt accession
├── molecular_weight: Protein MW (kDa)
├── cellular_component: Main subcellular location
├── biological_process: Primary GO biological process
└── variance: Variance across samples
```

#### Key Proteins of Interest
| Symbol | Name | Function | Exercise Focus |
|--------|------|----------|----------------|
| SQSTM1 | Sequestosome-1 | Autophagy receptor | All exercises |
| VDAC1 | Voltage-dependent anion channel 1 | Mitochondrial transport | Temporal analysis |
| MAP1LC3B | LC3B | Autophagy marker | Exercise 2, 5 |
| PSMA1-7 | Proteasome subunits | Protein degradation | Exercise 1, 5 |
| HSPA5 | Heat shock 70kDa protein 5 | ER chaperone | Exercise 4, 5 |
| MAPT | Microtubule-associated protein tau | Cytoskeleton | Exercise 5 |

---

### 2. Subset Dataset: `proteomics_subset.csv`
**Simplified dataset for beginners**

#### Description
- **Samples**: Same 150 neurons
- **Proteins**: Top 100 most variable proteins
- **Format**: CSV for easy loading
- **Use**: Exercise 1 and initial exploration

#### Columns
```
Sample_ID, tau_status, age, sex, pmi, batch,
SQSTM1, VDAC1, MAP1LC3B, PSMA1, PSMA2, ... (100 proteins)
```

---

### 3. Pathway Reference: `pathway_definitions.csv`
**Curated pathway annotations**

#### Description
- **Pathways**: Protein quality control pathways
- **Source**: GO, KEGG, Reactome
- **Use**: Exercise 2 pathway analysis

#### Structure
```
pathway_id, pathway_name, protein_list, database_source
GO:0006511, "Ubiquitin-dependent protein catabolic process", "PSMA1;PSMA2;...", "GO"
hsa04141, "Protein processing in ER", "HSPA5;CALR;...", "KEGG"
```

---

### 4. Network Data: `protein_interactions.csv`
**Protein-protein interactions**

#### Description
- **Source**: STRING database (confidence > 0.7)
- **Proteins**: Focus on quality control proteins
- **Use**: Exercise 2, 3, 5 network analysis

#### Structure
```
protein_a, protein_b, confidence_score, interaction_type
SQSTM1, MAP1LC3B, 0.85, "binding"
PSMA1, PSMA2, 0.95, "complex"
```

---

## 🔗 Data Access Methods

### Option 1: Local Download (Recommended)

```python
# Download data files to your local machine
import urllib.request
import os

base_url = "https://github.com/your-repo/proteomics-data/"
data_files = [
    "alzheimer_proteomics.h5ad",
    "proteomics_subset.csv",
    "pathway_definitions.csv",
    "protein_interactions.csv"
]

# Create data directory
os.makedirs("data", exist_ok=True)

# Download files
for file in data_files:
    url = base_url + file
    local_path = f"data/{file}"
    urllib.request.urlretrieve(url, local_path)
    print(f"Downloaded: {file}")
```

### Option 2: Google Colab
```python
# Mount Google Drive and access shared folder
from google.colab import drive
drive.mount('/content/drive')

# Copy data to Colab environment
!cp "/content/drive/MyDrive/proteomics_course_data/*" ./data/

# Or download directly
!wget "https://github.com/your-repo/proteomics-data/alzheimer_proteomics.h5ad"
```

### Option 3: Jupyter Hub / Cloud Environment
```python
# If using institutional server with shared storage
import shutil

# Copy from shared location
shared_path = "/shared/proteomics_course/data/"
local_path = "./data/"
shutil.copytree(shared_path, local_path)
```

---

## 📋 Data Loading Examples

### Loading the Main Dataset
```python
import scanpy as sc
import pandas as pd

# Load H5AD file
adata = sc.read_h5ad('data/alzheimer_proteomics.h5ad')

# Basic exploration
print(f"Data shape: {adata.shape}")
print(f"Sample metadata: {list(adata.obs.columns)}")
print(f"Protein metadata: {list(adata.var.columns)}")

# Access expression matrix
expression_matrix = adata.X  # Sparse or dense array
protein_names = adata.var_names  # Gene symbols
sample_metadata = adata.obs  # Sample information
```

### Loading Subset Dataset
```python
# Load CSV file
df = pd.read_csv('data/proteomics_subset.csv', index_col=0)

# Separate metadata and expression
metadata_cols = ['tau_status', 'age', 'sex', 'pmi', 'batch']
metadata = df[metadata_cols]
expression = df.drop(columns=metadata_cols)

print(f"Samples: {len(df)}")
print(f"Proteins: {len(expression.columns)}")
```

### Loading Pathway Data
```python
# Load pathway definitions
pathways = pd.read_csv('data/pathway_definitions.csv')

# Create pathway dictionary
pathway_dict = {}
for _, row in pathways.iterrows():
    pathway_dict[row['pathway_id']] = {
        'name': row['pathway_name'],
        'proteins': row['protein_list'].split(';'),
        'source': row['database_source']
    }
```

---

## 🔍 Data Quality Information

### Completeness
- **Missing values**: <1% (primarily for PMI data)
- **Batch coverage**: Balanced across tau status
- **Age distribution**: Normal distribution (mean 72 ± 12 years)

### Technical Details
- **Normalization**: TMT reporter ion intensities, log2-transformed
- **Filtering**: Proteins detected in >80% of samples
- **Batch correction**: None applied (for educational purposes)
- **Outliers**: Identified but retained for learning

### Known Issues (Educational Features)
1. **Age confounding**: Tau+ samples slightly older (realistic scenario)
2. **Batch effects**: Present but moderate (practice detection)
3. **Missing data**: Some PMI values missing (practice imputation)
4. **Outliers**: 2-3 samples with extreme values (practice handling)

---

## 🛡️ Data Privacy and Ethics

### Synthetic Enhancement
While based on real proteomics principles, this dataset includes:
- ✅ Synthetic noise to protect privacy
- ✅ Simulated samples to balance groups
- ✅ Realistic but non-identifiable metadata
- ✅ Educational modifications for learning

### Usage Guidelines
- ✅ Educational and research use only
- ✅ Not for clinical applications
- ✅ Cite this educational resource if used
- ❌ Do not redistribute without permission

---

## 🔧 Troubleshooting Data Issues

### File Not Found
```python
import os
print("Current directory:", os.getcwd())
print("Files in data/:", os.listdir('data/'))

# Check if file exists
if os.path.exists('data/alzheimer_proteomics.h5ad'):
    print("Main dataset found!")
else:
    print("Please download the dataset first")
```

### Memory Issues
```python
# For large datasets, load subset first
adata = sc.read_h5ad('data/alzheimer_proteomics.h5ad')

# Work with top variable proteins
top_proteins = adata.var.nlargest(1000, 'variance').index
adata_subset = adata[:, top_proteins]

# Or specific proteins
proteins_of_interest = ['SQSTM1', 'VDAC1', 'MAP1LC3B']
adata_focused = adata[:, proteins_of_interest]
```

### Loading Errors
```python
# If H5AD loading fails, try CSV backup
try:
    adata = sc.read_h5ad('data/alzheimer_proteomics.h5ad')
except:
    print("H5AD failed, loading CSV backup...")
    df = pd.read_csv('data/proteomics_subset.csv', index_col=0)
```

---

## 📊 Data Generation Scripts

For transparency, here's how the educational dataset was created:

<details>
<summary>Click to see data generation approach</summary>

```python
def generate_educational_dataset():
    """
    Create realistic proteomics dataset for education
    """
    np.random.seed(42)  # Reproducibility

    # Sample metadata
    n_samples = 150
    metadata = pd.DataFrame({
        'tau_status': ['tau_positive'] * 75 + ['tau_negative'] * 75,
        'age': np.random.normal(72, 12, n_samples),
        'sex': np.random.choice(['M', 'F'], n_samples),
        'pmi': np.random.exponential(8, n_samples),
        'batch': np.random.choice([1, 2, 3], n_samples)
    })

    # Add realistic confounding (tau+ slightly older)
    tau_pos_idx = metadata['tau_status'] == 'tau_positive'
    metadata.loc[tau_pos_idx, 'age'] += 3

    # Protein expression matrix
    n_proteins = 5853
    base_expression = np.random.normal(5, 2, (n_samples, n_proteins))

    # Add biological effects
    # SQSTM1: upregulated in tau+
    sqstm1_idx = protein_names.index('SQSTM1')
    base_expression[tau_pos_idx, sqstm1_idx] += 1.5

    # Age effects
    age_sensitive = np.random.choice(n_proteins, 500)
    for idx in age_sensitive:
        base_expression[:, idx] += 0.02 * (metadata['age'] - 70)

    # Batch effects
    batch_effects = np.random.normal(0, 0.3, (3, n_proteins))
    for i, batch in enumerate([1, 2, 3]):
        batch_mask = metadata['batch'] == batch
        base_expression[batch_mask, :] += batch_effects[i, :]

    return base_expression, metadata

# Note: This is conceptual - actual dataset uses more sophisticated modeling
```
</details>

---

## 📚 Related Resources

### Documentation Links
- [H5AD Format Guide](https://anndata.readthedocs.io/en/latest/)
- [Scanpy Tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html)
- [Pandas Data Loading](https://pandas.pydata.org/docs/user_guide/io.html)

### External Databases
- [UniProt](https://www.uniprot.org/) - Protein information
- [STRING](https://string-db.org/) - Protein interactions
- [Gene Ontology](http://geneontology.org/) - Functional annotations
- [KEGG](https://www.kegg.jp/) - Pathway database

---

**Ready to start analyzing? Return to [Practice Exercises](../07_practice_exercises/) or [Main Guide](../README.md)**
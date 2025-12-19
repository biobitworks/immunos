# 📊 Understanding Your Proteomics Dataset

## 🎯 What You'll Learn

By the end of this guide, you'll understand:
- ✅ **What's in the proteomics dataset** we'll be analyzing
- ✅ **How the data was generated** using mass spectrometry
- ✅ **What each number means** in the dataset
- ✅ **How samples are organized** and what metadata we have
- ✅ **What biological questions** we can answer with this data

---

## 🧬 The Big Picture: What Are We Analyzing?

### Our Research Context

#### The Biological Question
**Central Question**: How do protein levels change in neurons affected by Alzheimer's disease?

**Specific Focus**:
- **Tau-positive neurons**: Neurons with tau pathology (more diseased)
- **Tau-negative neurons**: Neurons with less tau pathology (relatively healthier)
- **5,853 proteins**: Nearly the entire expressed neuronal proteome

#### Why This Matters
- **Disease mechanisms**: Understanding how proteins change reveals how disease progresses
- **Therapeutic targets**: Proteins that change dramatically could be drug targets
- **Biomarkers**: Protein changes could predict or track disease progression
- **Systems biology**: See how entire cellular networks are affected

### The Dataset We're Using

#### Source and Generation
- **Human brain tissue**: Post-mortem samples from Alzheimer's disease patients
- **Single-neuron resolution**: Individual neurons isolated and analyzed
- **Mass spectrometry**: Proteins identified and quantified using advanced MS
- **Pooled samples**: Small groups of neurons (mini-pools) to get enough material

#### What Makes This Dataset Special
- **Disease-relevant**: Real human neurons affected by Alzheimer's
- **Comprehensive**: Nearly complete proteome coverage
- **High resolution**: Individual neuron variation captured
- **Well-controlled**: Multiple covariates measured and accounted for

---

## 📁 Dataset Structure and Format

### File Format: H5AD (AnnData)

#### What is H5AD?
- **Hierarchical Data Format**: Efficiently stores large multidimensional arrays
- **AnnData structure**: Standard format for single-cell/proteomics data
- **Compressed**: Much smaller than equivalent CSV files
- **Fast access**: Quick loading and subsetting

#### Why H5AD for Proteomics?
```python
# Traditional CSV approach (problematic):
# File size: 150 samples × 5,853 proteins = 878,950 values
# As CSV: ~50-100 MB, slow to load, memory intensive

# H5AD approach (better):
# Compressed: ~10-20 MB
# Fast loading: Seconds instead of minutes
# Metadata integrated: All info in one file
# Subsetting: Can load only parts of data
```

### Dataset Dimensions

#### Overall Structure
```python
# Dataset shape: (n_observations, n_variables)
# Observations = neuronal samples (rows)
# Variables = proteins (columns)

Expected dimensions: (~150 samples, 5,853 proteins)
```

#### What Each Dimension Represents
- **Rows (Observations)**: Individual neuronal samples (mini-pools)
- **Columns (Variables)**: Proteins detected and quantified
- **Values**: Log2-transformed protein expression levels

### Data Organization

#### Main Data Matrix (X)
```python
# adata.X contains the expression matrix
# Shape: (150 samples, 5,853 proteins)
# Values: Log2-transformed protein intensities
# Range: Typically 10-25 (log2 scale)
```

#### Sample Metadata (obs)
```python
# adata.obs contains information about each sample
Key columns:
- tau_status: 'positive' or 'negative' (our main comparison)
- MC1_score: Quantitative measure of tau pathology
- pseudotime: Disease progression estimate
- age: Patient age at death
- sex: Patient sex
- PMI: Post-mortem interval
- patient_id: Which patient the sample came from
- batch: Technical batch for processing
```

#### Protein Information (var)
```python
# adata.var contains information about each protein
Key columns:
- gene_names: Standard gene symbols (e.g., 'SQSTM1', 'VDAC1')
- protein_ids: UniProt identifiers
- highly_variable: Proteins with high variance across samples
- mean_expression: Average expression level
- detection_rate: Fraction of samples where protein detected
```

---

## 🔬 How the Data Was Generated

### Sample Collection and Preparation

#### Brain Tissue Collection
1. **Post-mortem brain tissue** from Alzheimer's disease patients
2. **Specific brain regions** most affected by tau pathology
3. **Fresh-frozen preservation** to maintain protein integrity
4. **Quality control** for RNA/protein preservation

#### Single-Neuron Isolation
```
Brain Tissue → Nuclei Isolation → FACS Sorting → Mini-Pool Creation
     |               |              |              |
  Homogenize    Separate by     Sort neurons    Pool 10-20
   tissue       cell type        by markers     neurons
```

#### Why Mini-Pools?
- **Sensitivity**: Single neurons don't have enough protein for detection
- **Reliability**: Pooling reduces technical noise
- **Biology**: Still captures cell-to-cell variation
- **Practical**: Balances resolution with detectability

### Mass Spectrometry Analysis

#### Sample Processing
1. **Protein extraction**: Release proteins from neuronal pools
2. **Digestion**: Break proteins into peptides (smaller pieces)
3. **Cleanup**: Remove contaminants that interfere with MS
4. **Labeling**: Add tags for quantification (TMT/iTRAQ)

#### Mass Spectrometry Detection
```
Peptides → Ionization → Mass Analysis → Detection → Quantification
    |           |            |           |           |
  LC-MS/MS   Electric    m/z ratio   Intensity   Protein levels
           acceleration  measurement   counting
```

#### Data Processing Pipeline
1. **Raw data**: MS instrument output
2. **Peptide identification**: Match spectra to sequences
3. **Protein inference**: Combine peptides into proteins
4. **Quantification**: Calculate protein abundances
5. **Normalization**: Correct for technical variation
6. **Quality control**: Filter unreliable measurements

---

## 📈 Understanding the Data Values

### Expression Scale: Log2 Transformation

#### Why Log2 Transform?
```python
# Original protein intensities:
Raw values: [1000, 2000, 4000, 8000, 16000]
Problem: Huge range, right-skewed distribution

# After log2 transformation:
Log2 values: [10, 11, 12, 13, 14]
Benefit: Normal distribution, manageable range
```

#### Interpreting Log2 Values
```python
# Log2 expression examples:
log2_value = 10  → original_intensity = 2^10 = 1,024
log2_value = 15  → original_intensity = 2^15 = 32,768
log2_value = 20  → original_intensity = 2^20 = 1,048,576

# Fold change calculations:
# If protein A = 12 and protein B = 10:
# Fold change = 2^(12-10) = 2^2 = 4-fold higher
```

#### Typical Value Ranges
- **Low expression**: 8-12 (barely detectable)
- **Medium expression**: 12-18 (clearly detected)
- **High expression**: 18-22 (abundant proteins)
- **Very high expression**: 22+ (extremely abundant, e.g., housekeeping proteins)

### Expression Patterns

#### What Determines Protein Levels?
1. **Biological factors**: Cell type, disease state, age
2. **Technical factors**: Batch effects, sample quality
3. **Detection limits**: Some proteins too low to measure
4. **Dynamic range**: MS can detect 4-5 orders of magnitude

#### Expected Patterns in Our Data
```python
# Housekeeping proteins (should be stable):
# Ribosomal proteins, metabolic enzymes
# Expected: Similar levels across all samples

# Disease-related proteins (our focus):
# SQSTM1, tau, inflammatory markers
# Expected: Different between tau+ and tau- neurons

# Cell-type markers:
# Neuronal-specific proteins
# Expected: High in all samples (since all are neurons)
```

---

## 🏷️ Sample Metadata Explained

### Primary Variables

#### Tau Status (Our Main Comparison)
```python
# adata.obs['tau_status']
Values: 'positive' or 'negative'
Meaning:
- 'positive': Neurons with significant tau pathology
- 'negative': Neurons with minimal tau pathology
Biological significance: This is our primary disease comparison
```

#### MC1 Score (Disease Severity)
```python
# adata.obs['MC1_score']
Values: Continuous scale (0-10)
Meaning: Quantitative measure of misfolded tau
Use:
- Correlation with protein changes
- Disease severity assessment
- Continuous analysis (not just positive/negative)
```

#### Pseudotime (Disease Progression)
```python
# adata.obs['pseudotime']
Values: 0.0 to 1.0
Meaning: Computational estimate of disease stage
Use:
- 0.0 = early disease stage
- 1.0 = late disease stage
- Sliding window analysis
```

### Technical Covariates

#### Patient-Level Variables
```python
# adata.obs['patient_id']
Values: Patient identifiers
Use: Control for individual differences

# adata.obs['age']
Values: Age at death (years)
Use: Age affects protein expression

# adata.obs['sex']
Values: 'M' or 'F'
Use: Sex differences in protein expression
```

#### Sample Quality Variables
```python
# adata.obs['PMI']
Values: Post-mortem interval (hours)
Meaning: Time between death and tissue preservation
Use: Longer PMI can affect protein detection

# adata.obs['batch']
Values: Processing batch identifiers
Use: Technical batch effects correction
```

### Why Covariates Matter

#### Confounding Variables
```python
# Example problem:
# If tau-positive patients are older on average:
# Protein differences could be due to:
# 1. Tau pathology (what we want to study)
# 2. Age effects (confounding variable)

# Solution: Statistical control
# protein_level ~ tau_status + age + sex + PMI + batch
# This separates tau effects from age effects
```

---

## 🧪 Protein Information

### Protein Naming and Identification

#### Gene Names (Primary Identifiers)
```python
# adata.var_names (protein identifiers)
Examples:
- 'SQSTM1': Sequestosome 1 (our key protein)
- 'VDAC1': Voltage-dependent anion channel 1
- 'ACTB': Beta-actin (housekeeping protein)
- 'GAPDH': Glyceraldehyde-3-phosphate dehydrogenase
```

#### Protein Classification
```python
# Major protein categories in our dataset:
1. Structural proteins: Cytoskeleton, membranes
2. Metabolic enzymes: Energy production, biosynthesis
3. Regulatory proteins: Kinases, phosphatases, transcription factors
4. Quality control: Proteases, chaperones, autophagy
5. Signaling proteins: Receptors, ion channels
6. Disease-related: Amyloid, tau, inflammation
```

### Protein Detection Quality

#### Detection Rate
```python
# adata.var['detection_rate']
Values: 0.0 to 1.0
Meaning: Fraction of samples where protein detected
Quality thresholds:
- >0.8: Reliably detected (good for analysis)
- 0.5-0.8: Moderately detected (use with caution)
- <0.5: Poorly detected (exclude from analysis)
```

#### Expression Variability
```python
# adata.var['highly_variable']
Values: True/False
Meaning: Proteins with high variance across samples
Use: These proteins show interesting biological variation
```

---

## 🎯 Research Questions We Can Answer

### Group 1 Analyses (Focused)

#### 1. UPS System Analysis
```python
Question: "Are protein quality control systems intact?"
Data needed: UPS protein expression levels
Comparison: Tau-positive vs tau-negative neurons
Expected proteins: ~50 UPS components
```

#### 2. SQSTM1 Upregulation
```python
Question: "Is SQSTM1 really increased 10.7-fold?"
Data needed: SQSTM1 expression levels
Analysis: Bootstrap confidence intervals
Interpretation: Autophagy dysfunction severity
```

#### 3. Temporal Dynamics
```python
Question: "How do protein relationships change over disease progression?"
Data needed: SQSTM1, VDAC1, pseudotime
Analysis: Sliding window correlations
Insight: System failure patterns
```

### Group 2 Analyses (Broad)

#### 1. Proteome-Wide Impact
```python
Question: "What percentage of proteins are affected by disease?"
Data needed: All 5,853 proteins
Analysis: Differential expression with FDR correction
Expected: ~36% of proteins significantly changed
```

#### 2. Covariate Effects
```python
Question: "How do technical factors affect our results?"
Data needed: Expression + age + sex + PMI + batch
Analysis: Linear mixed models
Goal: Separate biological from technical effects
```

#### 3. Pathway Analysis
```python
Question: "Which biological pathways are most affected?"
Data needed: Protein fold changes + pathway databases
Analysis: Gene set enrichment analysis
Insight: Systems-level disease mechanisms
```

---

## 🔍 Data Quality Assessment

### What to Check Before Analysis

#### Overall Data Quality
```python
# Questions to ask:
1. Do we have the expected number of samples and proteins?
2. Are expression values in reasonable ranges?
3. Do technical replicates cluster together?
4. Are there obvious batch effects?
5. Do positive controls behave as expected?
```

#### Sample Quality Indicators
```python
# Good quality samples:
- Detected proteins: >3,000 proteins
- Expression range: 8-25 log2 units
- Housekeeping proteins: Stable across samples
- Batch effects: Minimal

# Poor quality samples:
- Detected proteins: <2,000 proteins
- Extreme values: Outside normal ranges
- Missing data: Many proteins undetected
- Technical artifacts: Unusual patterns
```

#### Protein Quality Indicators
```python
# Reliable proteins:
- Detection rate: >80% of samples
- Biological variation: Higher than technical variation
- Known biology: Consistent with literature

# Unreliable proteins:
- Detection rate: <50% of samples
- High technical noise: Batch-dependent
- Impossible values: Outside physical limits
```

---

## 📊 Initial Data Exploration Strategy

### Step 1: Basic Data Inspection
```python
# Load and examine dataset structure
import scanpy as sc
adata = sc.read_h5ad('pool_processed_v2.h5ad')

print(f"Dataset shape: {adata.shape}")
print(f"Sample metadata: {adata.obs.columns.tolist()}")
print(f"Protein metadata: {adata.var.columns.tolist()}")
```

### Step 2: Sample Overview
```python
# Examine sample distribution
print("Tau status distribution:")
print(adata.obs['tau_status'].value_counts())

print("\nAge distribution:")
print(adata.obs['age'].describe())

print("\nSex distribution:")
print(adata.obs['sex'].value_counts())
```

### Step 3: Expression Overview
```python
# Basic expression statistics
import numpy as np

print("Expression value distribution:")
print(f"Min: {np.min(adata.X):.2f}")
print(f"Max: {np.max(adata.X):.2f}")
print(f"Mean: {np.mean(adata.X):.2f}")
print(f"Median: {np.median(adata.X):.2f}")
```

### Step 4: Key Protein Check
```python
# Verify our key proteins are present
key_proteins = ['SQSTM1', 'VDAC1', 'PSMA1', 'PSMB1']
present_proteins = [p for p in key_proteins if p in adata.var_names]
missing_proteins = [p for p in key_proteins if p not in adata.var_names]

print(f"Key proteins present: {present_proteins}")
print(f"Key proteins missing: {missing_proteins}")
```

---

## 🎯 Learning Objectives Check

After reading this overview, you should understand:

### Dataset Fundamentals
- [ ] What H5AD format contains and why it's used
- [ ] Dataset dimensions and what they represent
- [ ] How samples and proteins are organized
- [ ] What log2 transformation means for interpretation

### Biological Context
- [ ] Why we're comparing tau-positive vs tau-negative neurons
- [ ] What kinds of proteins are in the dataset
- [ ] How mass spectrometry generates protein data
- [ ] What biological questions we can answer

### Technical Details
- [ ] What metadata variables mean and why they matter
- [ ] How to assess data quality
- [ ] What value ranges indicate good vs poor data
- [ ] Why covariate control is important

### Analysis Preparation
- [ ] How to load and inspect the dataset
- [ ] What to check before starting analysis
- [ ] How to identify key proteins for analysis
- [ ] Why different analyses require different approaches

---

## 🚀 Next Steps

### Immediate Actions
1. **Download the dataset** if you haven't already
2. **Set up your analysis environment** (local Jupyter or Google Colab)
3. **Load the data** and run basic inspection commands
4. **Verify data quality** using the checks described above

### Prepare for Analysis
- **Understand the biological questions** we're investigating
- **Review statistical concepts** needed for proteomics analysis
- **Plan your analysis strategy** based on research questions
- **Set up organized file structure** for reproducible analysis

### Continue Learning
- **Next guide**: Hands-on data exploration
- **Then**: Group 1 focused analyses (UPS, SQSTM1, temporal)
- **Finally**: Group 2 proteome-wide analyses

---

**You now have the foundation to understand what's in your proteomics dataset and why it's structured the way it is!** The next step is hands-on exploration to see these concepts in action.

*Next: [Hands-On Data Exploration](exploring_your_data.md)*

*Remember: Understanding your data deeply is the foundation of all good analysis - take time to explore before rushing into statistics!* 📊🧬
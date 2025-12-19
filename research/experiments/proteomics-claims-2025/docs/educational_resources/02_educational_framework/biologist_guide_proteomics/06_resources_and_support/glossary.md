# 📚 Proteomics and Statistics Glossary

## 🎯 How to Use This Glossary

This comprehensive glossary defines key terms you'll encounter in proteomics analysis:
- ✅ **Biological concepts** - proteins, pathways, disease mechanisms
- ✅ **Technical terms** - mass spectrometry, data formats, software
- ✅ **Statistical concepts** - p-values, effect sizes, multiple testing
- ✅ **Computational terms** - algorithms, data structures, programming

**Format**: Each term includes definition, context, and often an example to help you understand practical usage.

---

## 🧬 Biological Terms

### A

**Abundance**
The amount or concentration of a protein in a sample. In proteomics, abundance is typically measured as signal intensity from mass spectrometry and expressed in relative terms (e.g., protein A is 2-fold more abundant than protein B).

**Alzheimer's Disease**
A progressive neurodegenerative disorder characterized by accumulation of amyloid plaques and tau tangles in the brain, leading to cognitive decline and neuronal death.

**Amyloid**
Misfolded proteins that aggregate into insoluble fibrils. In Alzheimer's disease, amyloid-beta peptides form extracellular plaques that disrupt neuronal function.

**AnnData**
A data structure used in single-cell analysis (adapted for proteomics) that stores expression data (X), sample metadata (obs), and feature metadata (var) in a single object.

**Autophagy**
A cellular process that degrades and recycles damaged proteins and organelles. Dysfunction in autophagy is implicated in neurodegenerative diseases. Key proteins include SQSTM1 (p62) and LC3B.

### B

**Batch Effect**
Systematic technical variation in data that arises from processing samples in different batches or at different times, rather than from biological differences.

**Biomarker**
A measurable biological indicator of a disease state, treatment response, or biological process. Protein biomarkers can be diagnostic, prognostic, or predictive.

**Bootstrap**
A statistical resampling method that estimates the sampling distribution of a statistic by repeatedly sampling with replacement from the original data.

### C

**Cohen's d**
A standardized measure of effect size that quantifies the difference between two groups in terms of standard deviation units. Small (0.2), medium (0.5), and large (0.8) effects.

**Confidence Interval (CI)**
A range of values that likely contains the true population parameter with a specified probability (usually 95%). Provides information about the precision of estimates.

**Covariate**
A variable that is measured alongside the main variables of interest and may influence the outcome. Examples include age, sex, and batch in proteomics studies.

**Cross-validation**
A method for assessing how well a statistical model will generalize to new data by partitioning data into training and testing sets.

### D

**Differential Expression**
The comparison of protein abundance between different conditions (e.g., disease vs. healthy) to identify proteins that are significantly changed.

**Dynamic Range**
The range of protein concentrations that can be detected and quantified by a measurement technique. Mass spectrometry typically has a dynamic range of 4-6 orders of magnitude.

### E

**Effect Size**
A quantitative measure of the magnitude of a phenomenon. Unlike p-values, effect sizes are not dependent on sample size and indicate practical significance.

**Enrichment Analysis**
A statistical method to determine whether a set of genes/proteins is significantly associated with biological pathways or functional categories.

**Experimental Design**
The structure of an experiment, including how samples are collected, grouped, and analyzed to test specific hypotheses while controlling for confounding factors.

### F

**False Discovery Rate (FDR)**
The expected proportion of false positives among all discoveries. FDR control is crucial in proteomics to account for testing thousands of proteins simultaneously.

**Fold Change**
The ratio of expression levels between two conditions. A 2-fold change means one condition has twice the abundance of the other.

**Functional Annotation**
The process of assigning biological functions, pathways, and other properties to proteins based on experimental evidence or computational prediction.

### G

**Gene Ontology (GO)**
A standardized vocabulary for describing gene and protein functions, including biological processes, cellular components, and molecular functions.

**GSEA (Gene Set Enrichment Analysis)**
A computational method that determines whether a defined set of genes/proteins shows statistically significant, concordant differences between two biological states.

### H

**H5AD Format**
A file format based on HDF5 that efficiently stores large-scale biological data including expression matrices and metadata. Commonly used with scanpy.

**Housekeeping Proteins**
Proteins that are consistently expressed across different conditions and are essential for basic cellular functions. Examples include ACTB (beta-actin) and GAPDH.

**Hypergeometric Test**
A statistical test used in enrichment analysis to determine if the overlap between two sets (e.g., significant proteins and pathway proteins) is greater than expected by chance.

### I

**Imputation**
The process of replacing missing values in a dataset with estimated values based on other available information.

**Interaction Effect**
When the effect of one variable depends on the level of another variable. In proteomics, this might occur when the effect of disease depends on patient age.

### J

**Jupyter Notebook**
An interactive computing environment that allows you to combine code, visualizations, and narrative text in a single document.

### K

**KEGG Pathway**
A collection of manually curated pathway maps representing molecular interaction and reaction networks for metabolism, genetic information processing, and cellular processes.

### L

**Label-free Quantification**
A mass spectrometry approach that determines protein abundance by measuring signal intensity without using isotopic labels.

**Linear Model**
A statistical model that assumes a linear relationship between predictor variables and the outcome. Commonly used in proteomics to control for covariates.

**Log Transformation**
Converting data to logarithmic scale, typically log2 in proteomics. This normalizes skewed distributions and makes fold changes symmetric.

### M

**Machine Learning**
Computational methods that can learn patterns from data and make predictions. Used in proteomics for biomarker discovery and classification.

**Mann-Whitney U Test**
A non-parametric statistical test that compares two independent groups without assuming normal distributions. Alternative to t-test for skewed data.

**Mass Spectrometry (MS)**
An analytical technique that measures the mass-to-charge ratio of ions to identify and quantify proteins and other molecules.

**Meta-analysis**
A statistical method that combines results from multiple independent studies to increase statistical power and generalizability.

**Missing Values**
Data points that are absent from the dataset, common in proteomics when proteins are below the detection limit.

**Multiple Testing Correction**
Statistical adjustment for conducting many statistical tests simultaneously to control the overall error rate.

### N

**Normalization**
The process of adjusting data to remove systematic technical variation while preserving biological differences.

**Null Hypothesis**
The assumption that there is no difference between groups or no effect of an intervention. Statistical tests aim to reject or fail to reject this hypothesis.

### O

**Outlier**
A data point that differs significantly from other observations, potentially due to measurement error or genuine biological variation.

**Over-representation Analysis (ORA)**
A method to identify biological pathways that contain more significant genes/proteins than expected by chance.

### P

**Pathway**
A series of interactions between molecules in a cell that leads to a certain biological outcome. Examples include metabolic pathways and signaling cascades.

**p-value**
The probability of observing data at least as extreme as what was observed, assuming the null hypothesis is true. Lower p-values suggest stronger evidence against the null hypothesis.

**Post-translational Modification (PTM)**
Chemical changes to proteins after they are synthesized, such as phosphorylation, ubiquitination, or acetylation.

**Power Analysis**
Statistical calculation to determine the sample size needed to detect an effect of a given size with specified confidence.

**Proteasome**
A protein complex responsible for degrading unneeded or damaged proteins in cells. Components include PSMA and PSMB subunits.

**Proteome**
The entire set of proteins expressed by an organism, tissue, or cell type at a given time under specific conditions.

**Proteomics**
The large-scale study of proteins, including their structure, function, modification, and interactions.

**Pseudotime**
A computational measure that orders cells or samples along a trajectory representing biological progression, such as disease development.

### Q

**Quality Control (QC)**
Procedures to assess and ensure the reliability and accuracy of experimental data before analysis.

**Quantile Normalization**
A normalization method that transforms the distribution of each sample to have the same quantiles.

### R

**Random Forest**
A machine learning algorithm that combines multiple decision trees to make predictions or classifications.

**Reactome**
A database of biological pathways that provides detailed molecular-level descriptions of cellular processes.

**Regression Analysis**
Statistical method for modeling the relationship between a dependent variable and one or more independent variables.

**Replication**
Repeating experiments to verify results. Technical replicates use the same sample multiple times; biological replicates use different samples.

**Robust Statistics**
Statistical methods that are less sensitive to outliers and violations of assumptions compared to standard methods.

### S

**Sample Size**
The number of observations or samples in a study. Larger sample sizes generally provide more statistical power to detect effects.

**Scanpy**
A Python package for analyzing single-cell data, often adapted for proteomics analysis due to its efficient data structures.

**Significance Level (α)**
The threshold for determining statistical significance, commonly set at 0.05 (5% chance of false positive).

**Single-cell Proteomics**
Techniques to measure protein expression in individual cells, providing insights into cellular heterogeneity.

**SQSTM1 (p62)**
A protein that serves as an autophagy receptor, targeting damaged proteins and organelles for degradation. Accumulates when autophagy is impaired.

**Standard Deviation**
A measure of the amount of variation or dispersion in a dataset.

**Statistical Power**
The probability of correctly rejecting a false null hypothesis (detecting a true effect when one exists).

### T

**Tau Protein**
A microtubule-associated protein that becomes hyperphosphorylated and forms tangles in Alzheimer's disease.

**t-test**
A statistical test that compares the means of two groups, assuming normal distributions and equal variances.

**Temporal Analysis**
Analysis of how biological systems change over time or across disease progression stages.

**Transcriptomics**
The study of RNA transcripts, which provides information about gene expression but not protein abundance.

### U

**Ubiquitin-Proteasome System (UPS)**
A cellular mechanism for protein degradation involving tagging proteins with ubiquitin for destruction by the proteasome.

**Univariate Analysis**
Statistical analysis that examines one variable at a time, as opposed to multivariate analysis.

**Upregulation**
An increase in protein expression or activity compared to a control condition.

### V

**Variable**
A characteristic or property that can take different values. In proteomics, variables include protein abundances and sample metadata.

**Variance**
A measure of how spread out data points are from the mean. High variance indicates more dispersion.

**VDAC1**
Voltage-dependent anion channel 1, a mitochondrial outer membrane protein involved in metabolite transport and apoptosis.

**Violin Plot**
A visualization that combines aspects of box plots and kernel density plots to show the distribution of data.

### W

**Wilcoxon Test**
A non-parametric statistical test for comparing paired samples (Wilcoxon signed-rank) or independent samples (Wilcoxon rank-sum, equivalent to Mann-Whitney U).

### Z

**Z-score**
A standardized score that indicates how many standard deviations a value is from the mean. Used to identify outliers and normalize data.

---

## 📊 Statistical Concepts Explained

### Understanding P-values

**What they mean:**
- p < 0.001: Very strong evidence against null hypothesis
- p < 0.01: Strong evidence against null hypothesis
- p < 0.05: Moderate evidence against null hypothesis
- p > 0.05: Insufficient evidence to reject null hypothesis

**What they DON'T mean:**
- NOT the probability that the null hypothesis is true
- NOT the probability of replicating the result
- NOT a measure of effect size or biological importance

### Effect Size Guidelines

**Cohen's d (standardized difference):**
- 0.2: Small effect
- 0.5: Medium effect
- 0.8: Large effect

**Fold Change (ratio of means):**
- 1.0: No change
- 2.0: 2-fold increase (100% increase)
- 0.5: 2-fold decrease (50% of original)

**Correlation coefficients:**
- 0.1: Small correlation
- 0.3: Medium correlation
- 0.5: Large correlation

### Multiple Testing Correction

**Why needed:**
Testing thousands of proteins increases chance of false positives

**Common methods:**
- **Bonferroni**: Very conservative, controls family-wise error rate
- **FDR (Benjamini-Hochberg)**: Less conservative, controls false discovery rate
- **Permutation**: Empirical correction based on data resampling

**When to use:**
- Always when testing multiple hypotheses
- FDR preferred for exploratory proteomics
- Bonferroni for confirmatory studies

---

## 🔬 Technical Terms

### Data Formats

**H5AD**
- Hierarchical Data Format for AnnData
- Efficient storage of large matrices
- Includes expression data and metadata
- Compatible with scanpy and other tools

**CSV (Comma-Separated Values)**
- Simple text format
- Easy to read but inefficient for large datasets
- Compatible with Excel and most software

**HDF5**
- Hierarchical Data Format version 5
- Efficient for large scientific datasets
- Cross-platform compatibility

### Software and Tools

**Python**
- Programming language popular in data science
- Extensive libraries for scientific computing
- Free and open source

**Jupyter**
- Interactive notebook environment
- Combines code, visualizations, and text
- Excellent for exploratory data analysis

**Scanpy**
- Python package for single-cell analysis
- Adapted for proteomics data
- Efficient handling of large datasets

**Pandas**
- Python library for data manipulation
- Provides DataFrame structure
- Essential for data cleaning and analysis

**NumPy**
- Numerical computing library for Python
- Provides array operations and mathematical functions
- Foundation for most scientific Python packages

**SciPy**
- Scientific computing library
- Statistical functions and tests
- Optimization and signal processing

**Matplotlib/Seaborn**
- Python plotting libraries
- Create publication-quality figures
- Seaborn provides statistical visualizations

### Algorithms and Methods

**Principal Component Analysis (PCA)**
- Dimensionality reduction technique
- Identifies main sources of variation
- Useful for data visualization and quality control

**Hierarchical Clustering**
- Groups similar samples or proteins
- Creates tree-like structures (dendrograms)
- Helpful for identifying patterns

**K-means Clustering**
- Partitions data into k clusters
- Minimizes within-cluster variation
- Used for sample classification

**Random Forest**
- Ensemble machine learning method
- Combines multiple decision trees
- Good for feature selection and classification

---

## 🧪 Experimental Design Terms

### Study Types

**Cross-sectional Study**
- Observes subjects at a single time point
- Compares different groups simultaneously
- Cannot establish causation

**Longitudinal Study**
- Follows subjects over time
- Can track changes and disease progression
- More expensive but informative

**Case-Control Study**
- Compares diseased (cases) vs healthy (controls)
- Retrospective design
- Good for rare diseases

**Cohort Study**
- Follows a group of subjects over time
- Prospective design
- Can establish temporal relationships

### Sample Characteristics

**Biological Replicates**
- Different samples from the same condition
- Account for biological variation
- Essential for statistical inference

**Technical Replicates**
- Multiple measurements of the same sample
- Account for measurement error
- Help assess reproducibility

**Batch**
- Samples processed together at the same time
- Potential source of technical variation
- Must account for in analysis

**Confounding Variable**
- Variable associated with both exposure and outcome
- Can bias results if not controlled
- Examples: age, sex, batch effects

---

## 📈 Data Analysis Workflow Terms

### Data Preprocessing

**Quality Control**
- Assessment of data reliability
- Identification of outliers and problems
- Filtering of low-quality samples/features

**Normalization**
- Adjustment for systematic differences
- Makes samples comparable
- Various methods available (quantile, median, etc.)

**Imputation**
- Filling in missing values
- Various strategies (mean, median, model-based)
- Important for downstream analysis

**Transformation**
- Mathematical conversion of data
- Log transformation common in proteomics
- Stabilizes variance and normalizes distributions

### Statistical Analysis

**Univariate Analysis**
- Analysis of one variable at a time
- Simple but doesn't capture interactions
- Good starting point for exploration

**Multivariate Analysis**
- Analysis of multiple variables simultaneously
- Can model complex relationships
- More realistic but complex interpretation

**Supervised Learning**
- Machine learning with known outcomes
- Classification or regression problems
- Examples: disease prediction, biomarker discovery

**Unsupervised Learning**
- Machine learning without known outcomes
- Pattern discovery and clustering
- Examples: sample grouping, pathway discovery

### Results Interpretation

**Statistical Significance**
- Evidence against null hypothesis
- Determined by p-value threshold
- Doesn't necessarily mean biological importance

**Biological Significance**
- Practical importance of findings
- Considered alongside statistical significance
- Includes effect size and biological context

**Clinical Relevance**
- Potential impact on patient care
- Translation from research to practice
- Requires validation in clinical settings

---

## 🎯 Quick Reference Tables

### Statistical Test Selection

| Data Type | Normal Distribution | Non-Normal Distribution |
|-----------|-------------------|------------------------|
| Two groups | t-test | Mann-Whitney U |
| Paired samples | Paired t-test | Wilcoxon signed-rank |
| Multiple groups | ANOVA | Kruskal-Wallis |
| Correlation | Pearson | Spearman |

### Effect Size Interpretation

| Measure | Small | Medium | Large |
|---------|-------|---------|-------|
| Cohen's d | 0.2 | 0.5 | 0.8 |
| Correlation (r) | 0.1 | 0.3 | 0.5 |
| Fold Change | 1.2-1.5 | 1.5-2.0 | >2.0 |

### Sample Size Guidelines

| Effect Size | Power 80% | Power 90% |
|-------------|-----------|-----------|
| Small (d=0.2) | 400 per group | 525 per group |
| Medium (d=0.5) | 65 per group | 85 per group |
| Large (d=0.8) | 26 per group | 34 per group |

### P-value Interpretation

| P-value Range | Strength of Evidence | Common Interpretation |
|---------------|---------------------|----------------------|
| p < 0.001 | Very strong | Highly significant |
| 0.001 ≤ p < 0.01 | Strong | Very significant |
| 0.01 ≤ p < 0.05 | Moderate | Significant |
| 0.05 ≤ p < 0.1 | Weak | Marginally significant |
| p ≥ 0.1 | None | Not significant |

---

## 📖 Abbreviations and Acronyms

### Common Proteomics Abbreviations

- **AD**: Alzheimer's Disease
- **CI**: Confidence Interval
- **CV**: Coefficient of Variation
- **DE**: Differential Expression
- **ESI**: Electrospray Ionization
- **FC**: Fold Change
- **FDR**: False Discovery Rate
- **FWER**: Family-Wise Error Rate
- **GO**: Gene Ontology
- **GSEA**: Gene Set Enrichment Analysis
- **H5AD**: Hierarchical Data Format for AnnData
- **LC**: Liquid Chromatography
- **LC-MS/MS**: Liquid Chromatography-tandem Mass Spectrometry
- **LFQ**: Label-Free Quantification
- **MS**: Mass Spectrometry
- **ORA**: Over-Representation Analysis
- **PCA**: Principal Component Analysis
- **PMI**: Post-Mortem Interval
- **PTM**: Post-Translational Modification
- **QC**: Quality Control
- **RNA-seq**: RNA sequencing
- **RT**: Retention Time
- **SQSTM1**: Sequestosome 1 (also known as p62)
- **TMT**: Tandem Mass Tag
- **UPS**: Ubiquitin-Proteasome System
- **VDAC1**: Voltage-Dependent Anion Channel 1

### Statistical Abbreviations

- **ANOVA**: Analysis of Variance
- **AUC**: Area Under the Curve
- **CI**: Confidence Interval
- **CV**: Cross-Validation
- **FDR**: False Discovery Rate
- **IQR**: Interquartile Range
- **MAD**: Median Absolute Deviation
- **PCA**: Principal Component Analysis
- **ROC**: Receiver Operating Characteristic
- **SD**: Standard Deviation
- **SE**: Standard Error
- **SEM**: Standard Error of the Mean

### Software and File Format Abbreviations

- **API**: Application Programming Interface
- **CSV**: Comma-Separated Values
- **GUI**: Graphical User Interface
- **HDF**: Hierarchical Data Format
- **HTML**: HyperText Markup Language
- **IDE**: Integrated Development Environment
- **JSON**: JavaScript Object Notation
- **PDF**: Portable Document Format
- **SQL**: Structured Query Language
- **TSV**: Tab-Separated Values
- **URL**: Uniform Resource Locator
- **XML**: Extensible Markup Language

---

## 🔍 How to Use This Glossary Effectively

### For Learning
1. **Start with basics** - understand fundamental concepts before advanced terms
2. **Use cross-references** - many terms build on others
3. **Practice with examples** - apply definitions to real data
4. **Review regularly** - terminology becomes natural with repetition

### For Analysis
1. **Bookmark this page** - quick reference during analysis
2. **Check unfamiliar terms** - don't assume you know the meaning
3. **Verify statistical concepts** - ensure proper interpretation
4. **Share with collaborators** - establish common vocabulary

### For Communication
1. **Define terms** - when writing or presenting to mixed audiences
2. **Use standard terminology** - for consistency with literature
3. **Explain abbreviations** - spell out on first use
4. **Provide context** - adapt complexity to your audience

---

## 📚 Additional Resources for Deeper Understanding

### Online Databases and Tools
- **Gene Ontology**: http://geneontology.org/
- **KEGG Pathways**: https://www.genome.jp/kegg/pathway.html
- **Reactome**: https://reactome.org/
- **UniProt**: https://www.uniprot.org/
- **STRING**: https://string-db.org/

### Statistical Resources
- **Khan Academy Statistics**: Free online statistics course
- **Coursera Statistics**: University-level statistics courses
- **Cross Validated**: Stack Exchange for statistics questions
- **Nature Methods Statistics Guide**: Practical statistics for biologists

### Proteomics Literature
- **Journal of Proteome Research**: Premier proteomics journal
- **Molecular & Cellular Proteomics**: Technical and methodological advances
- **Nature Biotechnology**: High-impact biological technology
- **Bioinformatics**: Computational methods in biology

---

**This glossary is your reference companion throughout your proteomics journey. Don't hesitate to return here whenever you encounter unfamiliar terms - understanding the language is key to mastering the science!** 📚🔬✨

*Remember: Every expert was once a beginner who took the time to learn the terminology. Your investment in understanding these concepts will pay dividends throughout your research career!*
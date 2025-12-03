---
title: "Human Ageing Genomic Resources: updates on key databases in ageing research"
authors: [João Pedro de Magalhães, Zoya Abidi, Gabriel Arantes Dos Santos, Roberto A. Avelar, Diogo Barardo, Kasit Chatsirisupachai, Peter Clark, Evandro A. De-Souza, Emily J. Johnson, Inês Lopes, Guy Novoa, Ludovic Senez, Angelo Talay, Daniel Thornton, Paul Ka Po To]
year: 2024
venue: Nucleic Acids Research
type: journal-article
tags: [genage, aging-biology, longevity, genomics, database, bioinformatics, hagr, senescence]
relevance: data-source
status: referenced
project: aging-research
doi: 10.1093/nar/gkad927
url: https://doi.org/10.1093/nar/gkad927
pmid: 37933854
---

# Human Ageing Genomic Resources: GenAge Database (2024 Update)

## Full Citation

de Magalhães, J. P., Abidi, Z., Dos Santos, G. A., Avelar, R. A., Barardo, D., Chatsirisupachai, K., Clark, P., De-Souza, E. A., Johnson, E. J., Lopes, I., Novoa, G., Senez, L., Talay, A., Thornton, D., & To, P. K. P. (2024). Human Ageing Genomic Resources: updates on key databases in ageing research. *Nucleic Acids Research*, 52(D1), D900–D908. https://doi.org/10.1093/nar/gkad927

**BibTeX**: See [[../citations/genage-references.bib|genage-references.bib]]

---

## Summary

Comprehensive update on the Human Ageing Genomic Resources (HAGR), a collection of databases and tools for studying the biology of human aging. The 2024 version (Build 21) includes significant updates to six key databases, with GenAge being the primary database of aging-related genes.

**Key Numbers**:
- **GenAge Human**: 307 manually curated aging-related genes
- **GenAge Model Organisms**: 2,205 longevity/aging genes across 9 species
- **Gene Expression**: Meta-analysis from 127 datasets across mammals

---

## Key Databases in HAGR

### 1. **GenAge** (Primary Focus)

**Description**: Curated database of genes associated with aging and longevity

**Human Genes** (307 genes):
- Manually curated from literature
- Association with human aging process
- Extensive functional annotations
- Literature references for each gene
- **Downloaded**: [[../../data/genage/human/genage_human.csv|genage_human.csv]]

**Model Organisms** (2,205 genes):
- Yeast: *Saccharomyces cerevisiae* (911 genes)
- Worms: *Caenorhabditis elegans* (889 genes)
- Flies: *Drosophila melanogaster* (202 genes)
- Mice: *Mus musculus* (136 genes)
- Other organisms: 67 genes total
- **Downloaded**: [[../../data/genage/models/genage_models.csv|genage_models.csv]]

**Gene Expression Signatures**:
- Meta-analysis from 127 microarray/RNA-Seq datasets
- Species: Human, mouse, rat
- Genes consistently over/under-expressed with age
- Cross-species and cross-tissue analysis
- **Downloaded**: [[../../data/genage/expression/ageing_signatures.xlsx|ageing_signatures.xlsx]]

### 2. **DrugAge**

Compounds that extend or reduce lifespan in model organisms (database mentioned but not downloaded)

### 3. **GenDR**

Dietary restriction gene expression signatures

### 4. **LongevityMap**

Genetic variants associated with human longevity

### 5. **CellAge**

Genes associated with cell senescence

### 6. **AnAge**

Animal aging and longevity database (comparative biology)

---

## Methodology

### Data Curation Process
1. Manual literature review by domain experts
2. Criteria-based gene inclusion
3. Functional annotation from multiple sources
4. Cross-species comparative analysis
5. Regular updates (Build 21: August 2023)

### Quality Control
- Expert curation (not automated)
- Requirement for experimental evidence
- Literature citation for each entry
- Cross-validation with other databases

---

## Key Findings & Insights

### 1. Conservation of Aging Mechanisms

**Cross-species patterns**:
- Many aging genes conserved across species
- Yeast and worms most extensively studied (911 and 889 genes)
- Mouse models provide mammalian context (136 genes)

**Implications**:
- Fundamental aging mechanisms may be evolutionarily conserved
- Model organism research translates to human aging
- Comparative analysis reveals universal vs species-specific factors

### 2. Gene Expression Changes with Age

**Consistent patterns**:
- Meta-analysis identifies genes reliably changing with age
- Both up-regulated and down-regulated genes
- Tissue-specific and pan-tissue changes

**Biological themes**:
- Mitochondrial function decline
- Increased inflammation
- DNA damage response
- Proteostasis failure

### 3. Molecular Hallmarks of Aging

Genes in GenAge map to known hallmarks:
- Genomic instability
- Telomere attrition
- Epigenetic alterations
- Loss of proteostasis
- Deregulated nutrient sensing
- Mitochondrial dysfunction
- Cellular senescence
- Stem cell exhaustion
- Altered intercellular communication

---

## Data Access & Usage

### Download Information

**Website**: https://genomics.senescence.info/
**Download Page**: https://genomics.senescence.info/download.html
**Build**: 21 (released August 28, 2023)
**Downloaded**: 2025-12-02

**Local Data**:
- [[../../data/genage/human/|Human genes]]
- [[../../data/genage/models/|Model organisms]]
- [[../../data/genage/expression/|Expression signatures]]

### License

**Creative Commons Attribution 3.0 Unported License**
- Free for commercial, educational, and research use
- Attribution required in publications
- Full data export and reuse permitted

**Citation Requirement**: Must cite de Magalhães et al. (2024) NAR in publications

---

## Research Applications

### 1. Immunosenescence Research

**Opportunities**:
- Cross-reference aging genes with immune system genes
- Identify overlap between aging and immune decline
- Study inflammatory aging (inflammaging)
- Understand age-related immune dysfunction

**Potential Analysis**:
```python
# Identify genes in both GenAge and immune pathways
aging_genes = load_genage_human()
immune_genes = load_immune_pathways()
overlap = aging_genes.intersection(immune_genes)
```

### 2. IMMUNOS-MCP Integration

**Training Data**:
- Use aging patterns for pattern recognition training
- Biological validation of immune system models
- Cross-domain learning (aging + immunity)

**Anomaly Detection**:
- Identify premature aging signatures
- Detect abnormal aging patterns
- Compare to healthy aging profiles

**Potential Application**:
- Train B Cell agent on safe vs accelerated aging patterns
- NK Cell detection of age-related anomalies
- QML-AiNet learning of aging trajectories

### 3. Comparative Genomics

**Cross-species Analysis**:
- Study conservation of longevity mechanisms
- Identify species-specific aging factors
- Understand evolutionary pressures on aging

**Research Questions**:
- Why do different species age at different rates?
- What genes explain exceptional longevity?
- Can we identify anti-aging interventions?

### 4. Drug Discovery & Repurposing

**Target Identification**:
- Genes in aging pathways as drug targets
- Cross-reference with DrugAge for known compounds
- Repurposing existing drugs for aging

**Intervention Strategies**:
- Targeting multiple hallmarks simultaneously
- Personalized aging interventions
- Preventive vs therapeutic approaches

---

## Related Work & Citations

### Previous HAGR Publications

1. **Tacutu et al. (2018)** - Previous major update
   - [[tacutu-2018-hagr|Tacutu et al. (2018) NAR]]
   - Citation: Tacutu, R., et al. (2018). Nucleic Acids Research, 46(D1), D1083–D1090

2. **de Magalhães et al. (2009)** - Gene expression meta-analysis
   - Citation: de Magalhães, J. P., Curado, J., & Church, G. M. (2009). Bioinformatics, 25(7), 875–881
   - First large-scale meta-analysis of aging gene expression

### Complementary Databases

- **GenAge**: Aging genes (this paper)
- **DrugAge**: Longevity compounds
- **LongevityMap**: Human longevity variants
- **CellAge**: Cellular senescence
- **AnAge**: Comparative animal aging

---

## Technical Details

### Data Schema

**Human Genes Table**:
```
GenAge ID | symbol | name | entrez gene id | uniprot | why
```

**Model Organisms Table**:
```
GenAge ID | symbol | name | organism | entrez gene id |
avg lifespan change (max obsv) | lifespan effect | longevity influence
```

**Fields**:
- **symbol**: Standard gene symbol
- **name**: Full gene name
- **entrez gene id**: NCBI identifier for cross-referencing
- **organism**: Species (for model organisms)
- **why**: Rationale for inclusion (human genes)
- **lifespan effect**: Increase/Decrease/Mixed
- **longevity influence**: Pro-longevity, Anti-longevity, etc.

### Data Format

- **Primary Format**: Tab-delimited ASCII (CSV)
- **Expression Data**: Excel (.xlsx)
- **Compression**: ZIP archives
- **Size**: ~10-25 MB total (all datasets)

---

## Integration with Our Research

### Completed

- ✅ Downloaded all three GenAge datasets (2025-12-02)
- ✅ Created centralized data directory structure
- ✅ Generated analysis summary report
- ✅ Created BibTeX citations
- ✅ Documentation in Obsidian vault

### Analysis Performed

See: [[../experiments/genage-analysis-2025-12-02|GenAge Analysis Report]]

**Key Stats**:
- Human genes: 307 entries
- Model organisms: 2,205 entries (9 species)
- Expression signatures: 127 datasets analyzed

### Planned

- [ ] Map GenAge genes to immune system pathways
- [ ] Identify aging-immunity gene overlap
- [ ] Train IMMUNOS-MCP agents on aging patterns
- [ ] Comparative analysis across species
- [ ] Visualization of aging gene networks

---

## Critical Assessment

### Strengths

1. **Manual Curation**: Expert review ensures high quality
2. **Comprehensive Coverage**: 307 human + 2,205 model organism genes
3. **Well-Documented**: Extensive literature references
4. **Regularly Updated**: Active maintenance (Build 21, 2023)
5. **Open Access**: Free to use with attribution
6. **Multi-Species**: Cross-species comparative analysis possible

### Limitations

1. **Human Gene Count**: 307 genes may not capture full complexity of aging
2. **Selection Bias**: Focus on well-studied genes/organisms
3. **Causation vs Correlation**: Association doesn't prove causation
4. **Complexity**: Aging is multifactorial, genes act in networks
5. **Species Bias**: Yeast and *C. elegans* dominate model organism data

### Opportunities for Improvement

- Integration with proteomics and metabolomics data
- Machine learning to predict new aging genes
- Network analysis of gene interactions
- Phenotypic annotations for each gene
- Tissue-specific aging patterns

---

## Connection to IMMUNOS-MCP

### Biological Validation

**Immune System Parallels**:
- Both systems have pattern recognition (antigens vs aging markers)
- Both involve surveillance and response
- Both show network effects and coordination
- Both decline with age (immunosenescence)

**Computational Insights**:
- Aging genes → Training data for pattern recognition
- Biological networks → Inspire algorithm design
- Immune-aging overlap → Research direction

### Potential Research Directions

1. **Immunosenescence Analysis**
   - Identify genes in both GenAge and immune pathways
   - Study age-related immune decline
   - Predict immunosenescence risk

2. **Pattern Recognition**
   - Train B Cell agent on aging vs healthy patterns
   - NK Cell detection of abnormal aging
   - Multi-agent analysis of aging trajectories

3. **Comparative Biology**
   - Cross-species aging mechanisms
   - Evolutionary insights for algorithm design
   - Species-specific vs universal patterns

4. **Drug Discovery**
   - Cross-reference GenAge with DrugAge
   - Identify immunomodulatory anti-aging compounds
   - Repurposing for aging-related immune dysfunction

---

## Tags for Search

#genage #aging-biology #longevity #genomics #database #bioinformatics #hagr #senescence #lifespan #model-organisms #gene-expression #immunosenescence #data-source

---

## Metadata

**Relevance**: Primary data source for aging genomics research
**Status**: Referenced, data downloaded and analyzed
**DOI**: 10.1093/nar/gkad927
**PMID**: 37933854
**License**: CC BY 3.0
**Downloaded**: 2025-12-02
**Build**: 21 (August 2023)

---

**Related Documentation**:
- **Data**: [[../../data/genage/README|GenAge Data README]]
- **Analysis**: [[../experiments/genage-analysis-2025-12-02|Analysis Report]]
- **Citations**: [[../citations/genage-references.bib|BibTeX References]]
- **Datasets**:
  - [[../datasets/genage-human|Human Genes Documentation]]
  - [[../datasets/genage-models|Model Organisms Documentation]]
  - [[../datasets/genage-expression|Expression Signatures Documentation]]

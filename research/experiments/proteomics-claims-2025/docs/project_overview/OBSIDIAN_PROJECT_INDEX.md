# 🧬 Proteomics Analysis Project Index

## 📊 Main Analysis Directories

### [[pertpy_dge_analysis/README|PertPy DGE Analysis]]
Modern differential gene expression analysis using PertPy/PyDESeq2
- [[pertpy_dge_analysis/requirements|Requirements]]
- [[pertpy_dge_analysis/01_data_preparation/prepare_for_pertpy|Data Preparation]]

### [[ups_bias_analysis/README|UPS Bias Analysis]]
Comprehensive documentation of literature selection biases
- [[ups_bias_analysis/01_literature_comparison/CHERRY_PICKED_UPS_COMPARISON|Cherry-Picked Proteins]]
- [[ups_bias_analysis/01_literature_comparison/UPS_PROTEIN_COVERAGE_COMPARISON|Coverage Comparison]]
- [[ups_bias_analysis/01_literature_comparison/LITERATURE_COMPARISON_SQSTM1|SQSTM1 Literature]]

---

## 🔬 Group 1: Mitochondrial Dysregulation

### Contractor Notebooks
- [[contractor_notebook_1_sequential_failure|Sequential Failure Analysis (Notebook)]]
- [[contractor_notebook_1_sequential_failure.md|Sequential Failure Analysis (Markdown)]]
- [[contractor_notebook_2_mitochondrial_sqstm1|Mitochondrial & SQSTM1 (Notebook)]]
- [[contractor_notebook_2_mitochondrial_sqstm1.md|Mitochondrial & SQSTM1 (Markdown)]]

### PertPy Analyses
- [[pertpy_dge_analysis/02_group1_mitochondrial/claim1_ups_proteins|Claim 1: UPS Proteins]]
- [[pertpy_dge_analysis/02_group1_mitochondrial/claim2_sqstm1_upregulation|Claim 2: SQSTM1 Upregulation]]
- [[pertpy_dge_analysis/02_group1_mitochondrial/claim3_temporal_dynamics|Claim 3: Temporal Dynamics]]

### Original Analyses
- [[01_research_analysis/group1_mitochondrial/statement1_ups_proteins|Statement 1: UPS Proteins]]
- [[01_research_analysis/group1_mitochondrial/statement2_sqstm1_upregulation|Statement 2: SQSTM1]]
- [[01_research_analysis/group1_mitochondrial/evaluation_tracker|Evaluation Tracker]]

---

## 🧪 Group 2: Proteostasis Failure

### PertPy Analyses
- [[pertpy_dge_analysis/03_group2_proteostasis/claim1_vatpase_subunits|Claim 1: V-ATPase Subunits]]

### Original Analyses
- [[01_research_analysis/group2_proteostasis/evaluation_tracker_group2|Evaluation Tracker]]

---

## 📈 Comprehensive Analyses

### Summary Documents
- [[pertpy_dge_analysis/04_comprehensive_analysis/all_claims_summary|All Claims Summary]]
- [[ENHANCED_DOCUMENTATION_SUMMARY|Enhanced Documentation]]
- [[EXECUTION_GUIDE|Execution Guide]]
- [[ANALYSIS_DOCUMENTATION_TEMPLATE|Analysis Template]]

### Results
- [[01_research_analysis/results/PAPER_FIGURE_REPLICATION_SUMMARY|Figure Replication]]
- [[01_research_analysis/results/master_analysis/summary_report|Master Summary]]

---

## 🎓 Educational Framework

### Biologist Guide
- [[02_educational_framework/biologist_guide_proteomics/README|Main Guide]]
- [[02_educational_framework/biologist_guide_proteomics/QUICKSTART|Quick Start]]
- [[02_educational_framework/biologist_guide_proteomics/INDEX|Full Index]]
- [[02_educational_framework/biologist_guide_proteomics/SETUP_GUIDE|Setup Guide]]

### Tutorials
- [[02_educational_framework/biologist_guide_proteomics/00_getting_started/README_START_HERE|Getting Started]]
- [[02_educational_framework/biologist_guide_proteomics/01_background_knowledge/proteomics_basics|Proteomics Basics]]
- [[02_educational_framework/biologist_guide_proteomics/tutorials/uniprot_ups_integration|UniProt Integration]]

---

## 🗂️ UPS Bias Analysis Details

### Literature Comparisons
- [[ups_bias_analysis/01_literature_comparison/LITERATURE_COMPARISON_SQSTM1|SQSTM1 Literature Review]]
- [[ups_bias_analysis/01_literature_comparison/UPS_PROTEIN_COVERAGE_COMPARISON|Coverage Analysis]]
- [[ups_bias_analysis/01_literature_comparison/CHERRY_PICKED_UPS_COMPARISON|Cherry-Picking Impact]]

### Validated Proteins
- [[ups_bias_analysis/02_validated_proteins/ups_analysis_validated|Validated UPS List (132 proteins)]]

### Bias Analysis Notebooks
- [[ups_bias_analysis/04_bias_analysis_notebooks/01_go_term_bias_analysis|GO Term Bias]]
- [[ups_bias_analysis/04_bias_analysis_notebooks/02_cherry_picking_impact|Cherry-Picking Impact]]

### Data Files
- [[ups_bias_analysis/05_data_files/ups_132_proteins_list|Complete Protein List]]

---

## 🛠️ Utilities & Configuration

### Helper Functions
- [[pertpy_dge_analysis/06_utilities/pertpy_helpers|PertPy Helper Functions]]

### Configuration
- [[CLAUDE.md|Claude Configuration]]
- [[FILE_INVENTORY|File Inventory]]

### Documentation
- [[README|Main README]]
- [[04_documentation/DELIVERABLES|Deliverables]]

---

## 📊 Data Files

### Source Data
- `data/pool_processed_v2.h5ad` - Main proteomics dataset

### Prepared Data
- [[pertpy_dge_analysis/01_data_preparation/prepared_for_pertpy.h5ad|PertPy-Ready Data]]

---

## 🔍 Key Findings

### Major Discoveries
1. **UPS Disruption**: 28.8% of 132 proteins significantly changed (not "preserved")
2. **SQSTM1**: 1.32-fold upregulation (not 3.4-fold as claimed)
3. **Autophagy vs UPS**: 57% vs 29% disruption rate
4. **V-ATPase**: Differential expression confirmed
5. **Temporal Dynamics**: Progressive mitochondrial dysfunction validated

### Literature Bias Findings
- Cherry-picked studies miss **87%** of significant changes
- Single GO terms miss **73-93%** of UPS components
- Our 132 proteins is **2.6-13x** more comprehensive than typical studies

---

## 📝 Notebooks by Category

### Sequential Failure Analysis
- [[contractor_notebook_1_sequential_failure]]
- [[01_research_analysis/notebooks/Sequential_Failure_of_Proteostasis_Mechanisms_Notebook]]

### Mitochondrial Dysfunction
- [[contractor_notebook_2_mitochondrial_sqstm1]]
- [[01_research_analysis/notebooks/Late_Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure_Notebook]]

### Master Analysis
- [[01_research_analysis/notebooks/Master_Analysis_Notebook]]

---

## 🎯 Quick Navigation

### For Claim Evaluation
1. Start with [[pertpy_dge_analysis/README|PertPy Analysis Overview]]
2. Check [[pertpy_dge_analysis/04_comprehensive_analysis/all_claims_summary|All Claims Summary]]
3. Review specific claims in Group 1 or Group 2 folders

### For UPS Bias Analysis
1. Read [[ups_bias_analysis/README|UPS Bias Overview]]
2. Review [[ups_bias_analysis/01_literature_comparison/CHERRY_PICKED_UPS_COMPARISON|Cherry-Picking Impact]]
3. Check [[ups_bias_analysis/05_data_files/ups_132_proteins_list|Complete Protein List]]

### For Educational Content
1. Start with [[02_educational_framework/biologist_guide_proteomics/QUICKSTART|Quick Start Guide]]
2. Follow [[02_educational_framework/biologist_guide_proteomics/learning_path_guide|Learning Path]]
3. Use [[02_educational_framework/biologist_guide_proteomics/INDEX|Complete Index]]

---

## 📅 Project Timeline

- **Data**: pool_processed_v2.h5ad (44 neurons)
- **Analysis Date**: December 2024
- **Methods**: PertPy/PyDESeq2, comprehensive protein coverage
- **Claims Evaluated**: 16 total (8 Group 1, 8 Group 2)

---

## 🏷️ Tags for Obsidian

#proteomics #differential-expression #PertPy #PyDESeq2 #UPS #autophagy #mitochondrial #neurodegeneration #alzheimers #tau-pathology #bias-analysis #cherry-picking #comprehensive-analysis

---

## 🔗 External Links

- [PertPy Documentation](https://pertpy.readthedocs.io/)
- [PyDESeq2 GitHub](https://github.com/owkin/PyDESeq2)
- [UniProt Database](https://www.uniprot.org/)
- [BioGRID Database](https://thebiogrid.org/)

---

*Last Updated: December 2024*
*Total Files: 100+ notebooks and documents*
*Analysis Method: PertPy/PyDESeq2 with comprehensive protein coverage*
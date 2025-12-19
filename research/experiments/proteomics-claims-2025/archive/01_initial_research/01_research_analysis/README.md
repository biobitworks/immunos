# Research Analysis Section

## Overview
This directory contains the primary research deliverables for evaluating biological claims about proteostasis failure and mitochondrial dysfunction in neurodegeneration.

## Structure

### 📓 notebooks/
Core research notebooks implementing comprehensive statistical analyses:
- **Sequential_Failure_of_Proteostasis_Mechanisms_Notebook.ipynb**
  - Evaluates 8 claims about proteostasis system failure
  - Implements covariate-controlled differential expression
  - Analyzes V-ATPase, lysosomal, and transport proteins

- **Late_Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure_Notebook.ipynb**
  - Evaluates 8 claims about mitochondrial dysfunction
  - Analyzes UPS proteins and SQSTM1/p62 upregulation
  - Implements sliding window and biphasic pattern detection

- **Master_Analysis_Notebook.ipynb**
  - Comprehensive pipeline for all 16 biological claims
  - Unified analysis framework

### 📊 group1_mitochondrial/
Analysis scripts for Finding Group 1 (Mitochondrial Dysregulation):
- statement1_ups_proteins.py - UPS protein differential expression
- statement2_sqstm1_upregulation.py - SQSTM1/p62 validation
- statement3_becn1_sqstm1_correlation.py - Autophagy marker correlations
- statement4_becn1_reduction.py - BECN1 downregulation analysis
- statement5_impaired_pathway.py - Mitophagy pathway analysis
- statement6_sliding_window.py - Temporal pattern analysis
- statement7_biphasic_patterns.py - Biphasic expression detection
- statement8_tau_correlation.py - Tau pathology correlations

### 🧬 group2_proteostasis/
Analysis scripts for Finding Group 2 (Proteostasis Failure):
- statement1_v_atpase.py - V-ATPase subunit analysis
- statement2_atp6v0a1.py - ATP6V0A1 upregulation
- statement3_covariate_de.py - Covariate-controlled DE
- statement4_organellar_perturbation.py - Organellar dysfunction
- statement5_sos_response.py - Stress response analysis
- statement6_segmented_regression.py - Breakpoint detection
- statement7_temporal_ordering.py - Failure sequence analysis
- statement8_homeostasis_collapse.py - System-wide collapse patterns

### 🤖 ai_automation/
Automated analysis tools:
- ai_agent.py - Core AI agent for automated evaluation
- run_analysis.py - Main execution script
- utils.py - Helper functions
- validators.py - Data validation tools

## Key Features

### Statistical Methods
- **Differential Expression**: Covariate-controlled linear models
- **Multiple Testing Correction**: Benjamini-Hochberg FDR, Bonferroni
- **Temporal Analysis**: Sliding window correlations
- **Pattern Detection**: Segmented regression, biphasic detection
- **Effect Size Calculation**: Cohen's d, correlation coefficients

### Data Processing
- **Input**: pool_processed_v2.h5ad (44 samples × 5,853 proteins)
- **Format**: AnnData with log2-transformed expression values
- **Metadata**: Patient ID, age, PMI, tau status, MC1, pseudotime

## Usage

### Running Research Notebooks
```bash
# Navigate to notebooks directory
cd 01_research_analysis/notebooks/

# Launch Jupyter
jupyter notebook

# Open desired notebook and run all cells
```

### Running Individual Analyses
```bash
# Example: Run SQSTM1 analysis
python group1_mitochondrial/statement2_sqstm1_upregulation.py

# Run all Group 1 analyses
for script in group1_mitochondrial/statement*.py; do
    python "$script"
done
```

### Automated Analysis
```bash
# Run complete analysis with AI agent
python ai_automation/run_analysis.py \
    --data ../03_data/pool_processed_v2.h5ad \
    --group 1 --group 2 \
    --output results/
```

## Output Structure

Each analysis produces:
1. **Statistical Results**: P-values, effect sizes, confidence intervals
2. **Visualizations**: Publication-quality plots
3. **Evaluation**: SUPPORTED/REFUTED/UNSURE with rationale
4. **Detailed Report**: Comprehensive findings and biological interpretation

## Quality Assurance

All analyses include:
- ✅ Multiple testing correction
- ✅ Effect size reporting
- ✅ Confidence intervals
- ✅ Assumption checking
- ✅ Error handling for missing data
- ✅ Biological context interpretation

## Dependencies

```python
# Core requirements
scanpy >= 1.9.0
pandas >= 1.3.0
numpy >= 1.21.0
scipy >= 1.7.0
statsmodels >= 0.12.0
scikit-learn >= 0.24.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
```

## Navigation

- [[../INDEX|Back to Index]]
- [[../03_data/README|Data Documentation]]
- [[../04_documentation/analysis_documentation/README|Analysis Methods]]
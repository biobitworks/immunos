# Claude Code Configuration

## Project Overview
Bioinformatics Finding Group Evaluation Framework for proteomic data analysis in Alzheimer's disease research.

## Recommended Commands

### Data Analysis & Statistics
```bash
# Run full statistical analysis
python ai_automation/run_analysis.py --data data/pool_processed_v2.h5ad --group 1

# Run specific statement analysis
python group1_mitochondrial/statement1_ups_proteins.py

# Automated evaluation
python ai_automation/run_analysis.py --automation --output results/
```

### Code Quality & Testing
```bash
# Type checking
python -m mypy ai_automation/ai_agent.py

# Code formatting
python -m black ai_automation/ group1_mitochondrial/ group2_proteostasis/

# Linting
python -m pylint ai_automation/ai_agent.py

# Run tests (if test files exist)
python -m pytest tests/ -v
```

### Development Workflow
```bash
# Install dependencies
pip install pandas numpy scipy scanpy statsmodels scikit-learn matplotlib seaborn

# Check data integrity
python -c "import scanpy as sc; adata = sc.read_h5ad('data/pool_processed_v2.h5ad'); print(f'Data loaded: {adata.shape}')"

# Generate comprehensive report
python ai_automation/run_analysis.py --group 1 --group 2 --output analysis_results/
```

## Project Structure Context

### Core Analysis Files
- `group1_mitochondrial/`: Mitochondrial dysregulation analyses
- `group2_proteostasis/`: Proteostasis failure analyses
- `ai_automation/`: Automated analysis tools
- `data/`: Datasets and reference PDFs

### Key Dependencies
- **scanpy**: Single-cell analysis (proteomic data)
- **scipy.stats**: Statistical testing
- **statsmodels**: Advanced statistical models
- **pandas/numpy**: Data manipulation
- **matplotlib/seaborn**: Visualization

## AI Assistant Guidelines

### When Helping with Code:
1. **Statistical Rigor**: Always apply proper multiple testing correction
2. **Biological Context**: Consider the neurodegeneration research context
3. **Reproducibility**: Ensure all analyses are well-documented
4. **Error Handling**: Check for missing proteins/data
5. **Visualization**: Create publication-quality plots

### Common Tasks:
- Adding new statement analyses following existing patterns
- Debugging statistical tests and p-value calculations
- Optimizing data processing for large proteomic datasets
- Creating visualizations for biological data
- Implementing new statistical methods

### Code Patterns to Follow:
```python
# Standard analysis structure
def analyze_statement(adata, protein_list):
    """
    Analyze biological statement with proper statistics

    Parameters:
    -----------
    adata : AnnData
        Proteomic dataset
    protein_list : list
        Proteins to analyze

    Returns:
    --------
    dict : Analysis results with p-values, effect sizes, evaluation
    """
    # 1. Data validation
    # 2. Statistical analysis
    # 3. Multiple testing correction
    # 4. Effect size calculation
    # 5. Biological interpretation
    # 6. Return structured results
```

### Avoid These Patterns:
- Running statistical tests without multiple testing correction
- Ignoring effect sizes (don't rely only on p-values)
- Hard-coding protein names without checking if they exist
- Creating analyses without proper error handling
- Missing biological context in interpretations

## File-Specific Notes

### `ai_automation/ai_agent.py`
- Core AI agent for automated analysis
- Handles multiple analysis types
- Returns structured evaluation results

### `group1_mitochondrial/statement2_sqstm1_upregulation.py`
- Comprehensive SQSTM1 analysis example
- Shows proper differential expression workflow
- Good template for similar analyses

### `group1_mitochondrial/statement6_sliding_window.py`
- Advanced temporal analysis
- Complex visualization examples
- Statistical trend analysis

### `references/statistical_methods_reference.md`
- Complete guide to all statistical methods
- Reference for proper implementation
- Examples and best practices

## Common Issues & Solutions

### Missing Proteins
```python
if protein not in adata.var_names:
    alternatives = [g for g in adata.var_names if protein.lower() in g.lower()]
    if alternatives:
        print(f"Protein {protein} not found. Alternatives: {alternatives}")
    return None
```

### Memory Issues with Large Datasets
```python
# Process in chunks for large datasets
chunk_size = 1000
for i in range(0, adata.n_vars, chunk_size):
    chunk = adata[:, i:i+chunk_size]
    # Process chunk
```

### Statistical Power Warnings
```python
# Always check sample sizes
n_tau_pos = sum(adata.obs['tau_status'] == 'positive')
n_tau_neg = sum(adata.obs['tau_status'] == 'negative')
if min(n_tau_pos, n_tau_neg) < 20:
    print("WARNING: Small sample size may affect statistical power")
```

## Evaluation Criteria

### Statement Evaluation Guidelines
- **SUPPORTED**: p < 0.05 (FDR corrected) AND effect size within ±20% of claimed
- **REFUTED**: p > 0.05 OR effect size opposite direction OR >50% different
- **UNSURE**: Missing data, technical issues, or borderline results

### Quality Checks
- [ ] Multiple testing correction applied
- [ ] Effect sizes reported
- [ ] Confidence intervals calculated
- [ ] Assumptions checked (normality, equal variance)
- [ ] Biological interpretation provided
- [ ] Code properly documented

## Project Goals

Transform biological claims into rigorous statistical evaluations while maintaining:
- **Scientific accuracy**
- **Computational reproducibility**
- **Clear documentation**
- **Biological relevance**

This framework serves as a model for evidence-based evaluation of research findings in computational biology.
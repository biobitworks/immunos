# 🚀 Execution Guide for Bioinformatics Analysis Framework

## Quick Start

### 1. Install Dependencies
```bash
pip install pandas numpy scipy scanpy statsmodels scikit-learn matplotlib seaborn tqdm
```

### 2. Run the Master Jupyter Notebook
```bash
jupyter notebook Master_Analysis_Notebook.ipynb
```

## Alternative Execution Methods

### Option A: Individual Analysis Scripts
```bash
# UPS protein analysis
python group1_mitochondrial/statement1_ups_proteins.py

# SQSTM1 upregulation analysis
python group1_mitochondrial/statement2_sqstm1_upregulation.py

# Sliding window correlation analysis
python group1_mitochondrial/statement6_sliding_window.py

# Covariate-controlled differential expression
python group2_proteostasis/statement1_covariate_de.py
```

### Option B: AI Automation Pipeline
```bash
# Run automated analysis for both groups
python ai_automation/run_analysis.py --automation --output results/

# Run specific group
python ai_automation/run_analysis.py --group 1 --output results/

# Full automated pipeline
python ai_automation/run_analysis.py --group 1 --group 2 --automation --output comprehensive_results/
```

## 🔧 Fixed Issues

### ✅ Data Loading Paths
- **Fixed**: All Python files now use correct relative path: `data/pool_processed_v2.h5ad`
- **Was**: `../data/pool_processed_v2.h5ad` (incorrect)

### ✅ Import Dependencies
- **Verified**: All required imports are present
- **Added**: Missing tqdm imports where needed
- **Confirmed**: sklearn imports are properly placed

### ✅ File Structure
- **Enhanced**: All analysis files have comprehensive documentation
- **Standardized**: Consistent data loading patterns
- **Improved**: Error handling and validation

## 📊 Analysis Capabilities

### Group 1: Late-Stage Mitochondrial Dysregulation
1. **UPS Protein Analysis**: Rigorous statistical evaluation with dual testing
2. **SQSTM1 Upregulation**: Bootstrap validation of extreme fold changes
3. **Sliding Window Analysis**: Advanced temporal correlation dynamics

### Group 2: Sequential Failure of Proteostasis
1. **Covariate-Controlled DE**: Linear models with confounding variable control
2. **AI Automation**: Large-scale automated analysis pipelines

## 🎯 Performance Expectations

### Full Analysis (5,853 proteins)
- **Time**: 10-30 minutes depending on system
- **Memory**: ~2-4 GB RAM recommended
- **Output**: Comprehensive statistical results + visualizations

### Subset Analysis (500 proteins) - For Testing
- **Time**: 1-3 minutes
- **Memory**: ~500 MB RAM
- **Output**: Representative results for validation

## 📈 Expected Results

### Statement Evaluations
- **UPS Proteins**: Expected "SUPPORTED" (no significant alterations)
- **SQSTM1**: Expected "SUPPORTED" (~3.4 log2FC)
- **Sliding Window**: Expected "SUPPORTED" (early negative, late positive correlation)
- **Covariate DE**: Expected "SUPPORTED" (~36% significant proteins)

### Visualizations Generated
- **Volcano plots**: Differential expression results
- **Correlation plots**: Temporal relationship changes
- **Bootstrap distributions**: Confidence intervals
- **Phase analysis**: Disease progression patterns

## 🛠️ Troubleshooting

### Common Issues

#### 1. Data File Not Found
```
Error: FileNotFoundError: data/pool_processed_v2.h5ad
```
**Solution**: Ensure data file is in the `data/` directory

#### 2. Missing Dependencies
```
Error: ModuleNotFoundError: No module named 'scanpy'
```
**Solution**: Install all dependencies using the pip command above

#### 3. Memory Issues
```
Error: MemoryError
```
**Solution**:
- Reduce `max_proteins` parameter in analysis functions
- Use subset analysis: `max_proteins=500`
- Close other applications to free memory

#### 4. Slow Performance
**Solutions**:
- Use subset analysis for testing
- Run individual analyses instead of full pipeline
- Ensure data is on SSD if possible

### Getting Help

1. **Check CLAUDE.md**: Contains AI assistant guidelines and common patterns
2. **Review Analysis Templates**: Templates show expected patterns
3. **Run with Subsets**: Test with smaller protein sets first
4. **Check Dependencies**: Ensure all required packages are installed

## 📝 Output Files

### Generated Results
- **CSV files**: Statistical results for each analysis
- **PNG files**: Publication-quality visualizations
- **Log files**: Detailed execution logs
- **JSON files**: Structured evaluation summaries

### Key Output Locations
- `Master_Analysis_Notebook.ipynb`: Complete interactive analysis
- `results/`: Automated analysis outputs
- Individual script outputs in respective directories

## 🔬 Scientific Rigor

### Statistical Standards Implemented
- **Multiple testing correction**: FDR (Benjamini-Hochberg)
- **Effect size calculation**: Cohen's d for biological significance
- **Bootstrap confidence intervals**: Robust uncertainty quantification
- **Covariate control**: Advanced linear modeling
- **Temporal analysis**: Sliding window dynamics

### Documentation Standards
- **Analytical rationale**: Detailed justification for each method
- **Biological context**: Disease relevance and pathway significance
- **Statistical assumptions**: Validation and violation handling
- **Interpretation guidelines**: Clear evaluation criteria

This framework provides a complete, scientifically rigorous solution for evaluating biological statements using proteomic data.
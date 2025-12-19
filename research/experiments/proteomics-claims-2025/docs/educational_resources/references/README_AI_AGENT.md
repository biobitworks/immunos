# AI Agent for Bioinformatics Analysis

## Overview

This AI agent system automates the analysis of proteomic data for finding group evaluations. It's designed to handle complex statistical analyses, correlations, and differential expression studies on Alzheimer's disease datasets.

## Components

### 1. **ai_agent.py** - Core Agent Architecture
- Main `BioinformaticsAgent` class
- Handles different analysis types (DE, correlation, segmented regression, etc.)
- Generates evaluation results (SUPPORTED/REFUTED/UNSURE)
- Produces comprehensive reports

### 2. **analysis_automation.py** - Analysis Automation Module
- `AnalysisAutomation` class with specialized analysis methods
- Covariate-controlled differential expression
- Sliding window correlations
- Segmented regression and breakpoint detection
- Biphasic behavior analysis
- Composite score calculations (V-ATPase, proteasome)

### 3. **run_analysis.py** - Main Orchestrator
- Command-line interface for running analyses
- Configurations for both finding groups
- Automated report generation
- JSON summary output

## Installation

### Prerequisites

```bash
pip install pandas numpy scipy scanpy statsmodels scikit-learn matplotlib seaborn
```

## Usage

### Basic Usage

```bash
# Analyze both finding groups
python run_analysis.py --data pool_processed_v2.h5ad

# Analyze specific group
python run_analysis.py --data pool_processed_v2.h5ad --group 1

# Run with automated analysis module
python run_analysis.py --data pool_processed_v2.h5ad --automation

# Specify output directory
python run_analysis.py --data pool_processed_v2.h5ad --output results/
```

### Command Line Arguments

- `--data`: Path to H5AD data file (default: pool_processed_v2.h5ad)
- `--group`: Finding group to analyze (1 or 2, analyzes both if not specified)
- `--output`: Output directory for results (default: analysis_output)
- `--automation`: Run additional automated analysis

## Output Structure

```
analysis_output/
├── Late-Stage_Mitochondrial_Dysregulation/
│   ├── analysis_report.md          # Detailed analysis report
│   └── evaluation_summary.json     # Summary of evaluations
├── Sequential_Failure_of_Proteostasis/
│   ├── analysis_report.md
│   └── evaluation_summary.json
└── automated_analysis/              # If --automation flag used
    ├── differential_expression.csv
    ├── sqstm1_correlations.csv
    ├── pseudotime_correlations.csv
    ├── sliding_window_correlation.csv
    └── biphasic_analysis.json
```

## Analysis Types Supported

### 1. Differential Expression
- Compares tau-positive vs tau-negative neurons
- Supports covariate control (age, PMI, PatientID)
- FDR correction using Benjamini-Hochberg

### 2. Correlation Analysis
- Pearson and Spearman correlations
- Bonferroni correction for multiple testing
- Supports protein-protein and protein-metadata correlations

### 3. Segmented Regression
- Identifies breakpoints in continuous relationships
- Used for MC1 vs protein expression analysis
- Returns slopes before and after breakpoint

### 4. Sliding Window Correlation
- Analyzes changing correlations along pseudotime
- Configurable window size and step size
- Identifies early vs late phase behaviors

### 5. Biphasic Behavior Analysis
- Detects two-phase patterns in protein systems
- Compares V-ATPase and proteasome systems
- Determines temporal ordering of failures

## Evaluation Criteria

The agent evaluates statements as:
- **SUPPORTED**: Analysis results match expected values within tolerance
- **REFUTED**: Analysis results contradict the statement
- **UNSURE**: Unable to perform analysis or missing data

## Extending the Agent

### Adding New Analysis Types

1. Add method to `BioinformaticsAgent` class:
```python
def _analyze_custom(self, params: Dict) -> AnalysisResult:
    # Your analysis code
    return AnalysisResult(...)
```

2. Register in `analyze_statement` method:
```python
analysis_methods = {
    'custom_analysis': self._analyze_custom,
    # ...
}
```

3. Add configuration to finding group:
```python
{
    'type': 'custom_analysis',
    'params': {
        'statement_id': 'G1_S1',
        # your parameters
    }
}
```

## Important Notes

1. **Data Format**: Expects H5AD format with specific metadata columns (tau_status, MC1, pseudotime)
2. **Log2 Transform**: Data is assumed to be already log2 transformed
3. **Memory Usage**: Large datasets may require significant memory
4. **Computation Time**: Full analysis can take 10-30 minutes depending on dataset size

## Troubleshooting

### Common Issues

1. **Missing proteins**: Agent will mark as UNSURE if required proteins not found
2. **Missing metadata**: Check that tau_status, MC1, pseudotime columns exist
3. **Memory errors**: Consider analyzing subsets of proteins
4. **Convergence failures**: Segmented regression may fail with poor initial parameters

### Logging

Detailed logs are saved to `analysis_YYYYMMDD_HHMMSS.log` for debugging.

## Example Python Usage

```python
from ai_agent import BioinformaticsAgent
from analysis_automation import AnalysisAutomation

# Initialize agent
agent = BioinformaticsAgent('pool_processed_v2.h5ad')
agent.load_data()

# Run single analysis
result = agent.analyze_statement(
    'correlation',
    {
        'statement_id': 'test',
        'protein1': 'SQSTM1',
        'protein2': 'VDAC1',
        'expected_correlation': 0.05
    }
)

print(f"Evaluation: {result.evaluation.value}")
print(f"Evidence: {result.evidence}")

# Use automation module
automation = AnalysisAutomation(agent.adata)
de_results = automation.covariate_controlled_de()
print(f"Significant proteins: {sum(de_results['FDR'] < 0.05)}")
```

## Citation

If you use this agent in your research, please cite the original proteomic dataset and analysis methodology.

## Support

For issues or questions, please refer to the project documentation or contact the development team.
# Analysis Documentation Template

## Required Documentation Standards

Every analysis file must include the following components with detailed explanations:

### 1. **Header Docstring with Complete Analytical Rationale**

```python
"""
Finding Group X - Statement Y: [Title]
Claim: [Exact statement being evaluated]

# Analytical Approach and Rationale

## Overview
[Detailed biological and statistical context]

## Statistical Strategy
1. **Primary method**: [Method chosen and why]
2. **Validation approach**: [Secondary methods for robustness]
3. **Multiple testing**: [Correction method and justification]
4. **Effect size**: [Measures used and interpretation]

## Rationale for Approach
- **Method selection**: [Why these specific statistical tests]
- **Data preprocessing**: [Any transformations and justification]
- **Assumption handling**: [How violations are addressed]
- **Confounding control**: [Variables considered and methods]

## Expected Outcome
[Clear criteria for SUPPORTED/REFUTED/UNSURE]

## Biological Context
[Disease relevance and pathway significance]
"""
```

### 2. **Data Loading with Comprehensive Comments**

```python
def load_data(file_path='../data/pool_processed_v2.h5ad'):
    """
    ## Data Loading Strategy
    [Explanation of file format choice and optimization]

    ## Quality Checks Performed
    [List of validation steps]
    """

    # Load [dataset type] in [format] format
    # Dataset characteristics:
    # - Source: [detailed source description]
    # - Preprocessing: [transformations applied]
    # - Scope: [dimensions and coverage]
    # - Metadata includes:
    #   * [variable]: [purpose and usage]
    #   * [variable]: [purpose and usage]
    # - Quality control: [filters and validations applied]
    # - Biological context: [disease stage, tissue type, etc.]
    data = load_function(file_path)
```

### 3. **Method Implementation with Decision Rationale**

```python
def analysis_function(data, parameters):
    """
    ## Method Selection Rationale
    [Why this specific approach over alternatives]

    ## Statistical Assumptions
    [What we assume and how we validate]

    ## Parameter Choices
    [Justification for thresholds, corrections, etc.]
    """

    # Step 1: [Action]
    # Rationale: [Why this step is necessary]
    # Alternative considered: [Other options and why rejected]

    # Step 2: [Action]
    # Statistical consideration: [Power, bias, assumptions]
    # Implementation note: [Technical details]
```

### 4. **Evaluation Logic with Clear Criteria**

```python
def evaluate_statement(results):
    """
    ## Evaluation Criteria
    [Detailed thresholds and reasoning]

    ## Decision Tree
    [Logic flow for SUPPORTED/REFUTED/UNSURE]
    """

    # Primary criterion: [statistical significance]
    # Secondary criterion: [effect size magnitude]
    # Tertiary criterion: [biological plausibility]

    if [condition]:
        evaluation = "SUPPORTED"
        explanation = "[detailed reasoning]"
    elif [condition]:
        evaluation = "REFUTED"
        explanation = "[detailed reasoning]"
    else:
        evaluation = "UNSURE"
        explanation = "[detailed reasoning]"
```

### 5. **Visualization with Analytical Purpose**

```python
def create_visualizations(data, results):
    """
    ## Visualization Strategy
    [What each plot reveals and why it's needed]
    """

    # Plot 1: [Type] - Shows [specific aspect] to assess [criterion]
    # Rationale: [Why this visualization is optimal for the data]

    # Plot 2: [Type] - Validates [assumption] through [visual test]
    # Interpretation guide: [How to read the plot]
```

---

## File-Specific Requirements

### Group 1 - Mitochondrial Dysregulation Files

#### `statement1_ups_proteins.py` ✅ (Enhanced)
- Focus: Conservative evaluation of "no significant alterations"
- Key decisions: Dual testing (parametric/non-parametric), FDR correction
- Biological emphasis: UPS integrity in neurodegeneration

#### `statement2_sqstm1_upregulation.py` (Needs Enhancement)
Required additions:
```python
# At top:
"""
## Why SQSTM1 is Critical
- Autophagy receptor accumulating when system fails
- Massive upregulation (3.413 log2FC) suggests major dysfunction
- Pseudotime correlation indicates progressive accumulation

## Statistical Challenges
- Extreme fold change requires careful validation
- FDR in context of 5,853 proteins
- Pseudotime regression assumptions

## Method Selection
- Bootstrap CI for robust fold change estimation
- Linear regression for pseudotime (assumes linear relationship)
- Multiple validation approaches for such extreme claims
"""

# At data loading:
    # SQSTM1/p62 is central to autophagy-lysosome pathway
    # This protein accumulates when autophagy flux is impaired
    # Expected to be highly upregulated in disease progression
    # Log2 scale: 3.413 = 10.7-fold increase (massive biological change)
```

#### `statement6_sliding_window.py` (Needs Enhancement)
```python
"""
## Advanced Temporal Analysis Rationale
- Static correlations miss dynamic relationships
- Sliding window reveals system transitions
- Window size (n=20) balances resolution vs stability

## Biological Interpretation
- Negative early correlation: Compensatory mechanism
- Positive late correlation: System co-failure
- Transition timing: Critical disease threshold
"""
```

### Group 2 - Proteostasis Files

#### `statement1_covariate_de.py` (Needs Enhancement)
```python
"""
## Covariate Control Strategy
- Age: Known effect on protein expression
- PMI: Technical confounding factor
- PatientID: Individual biological variation

## 36.14% Significance Interpretation
- Massive proteomic remodeling in disease
- Much higher than expected by chance (5%)
- Indicates system-wide dysfunction
"""
```

---

## Implementation Priority

### Immediate Updates Needed:
1. **statement2_sqstm1_upregulation.py** - Add analytical rationale for extreme fold change
2. **statement6_sliding_window.py** - Explain temporal dynamics methodology
3. **statement1_covariate_de.py** - Justify covariate selection and interpretation

### Template Application:
Each file should include:
- [ ] Detailed header with analytical approach
- [ ] Data loading comments explaining dataset characteristics
- [ ] Method selection rationale at each step
- [ ] Clear evaluation criteria with biological context
- [ ] Visualization purpose and interpretation guides

### Code Comment Standards:
- **Every statistical test**: Why chosen over alternatives
- **Every threshold**: Biological or statistical justification
- **Every transformation**: Purpose and validation
- **Every assumption**: How checked and handled
- **Every result**: Biological interpretation

This documentation transforms the analysis scripts from technical tools into educational resources that demonstrate rigorous analytical thinking in computational biology.
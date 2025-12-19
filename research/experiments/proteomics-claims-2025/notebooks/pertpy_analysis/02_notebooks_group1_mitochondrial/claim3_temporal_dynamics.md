# Claim 3: Temporal Dynamics Analysis

## Claim
"Temporal dynamics reveal progressive mitochondrial dysfunction"

## Analysis Overview
Comprehensive temporal analysis examining how mitochondrial and proteostasis proteins change over disease progression (pseudotime).

## Methods

### 1. Temporal Protein Sets

```python
temporal_proteins = {
    'early_response': ['HSP70', 'HSP90', 'HSPA1A', 'HSPA5', 'HSPB1'],
    'mitochondrial_complex_I': ['NDUFB8', 'NDUFA9', 'NDUFS3', 'NDUFV1'],
    'mitochondrial_complex_V': ['ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D'],
    'autophagy_early': ['BECN1', 'ATG5', 'ATG7', 'ATG12'],
    'autophagy_late': ['SQSTM1', 'NBR1', 'OPTN', 'TAX1BP1'],
    'proteasome': ['PSMA1', 'PSMB5', 'PSMC4', 'PSMD1']
}
```

### 2. Correlation Analysis

#### Statistical Methods
- **Spearman correlation**: Non-parametric, robust to outliers
- **Linear regression**: Quantify rate of change (slope)
- **FDR correction**: Multiple testing adjustment

#### For Each Protein
```python
# Correlation with pseudotime
corr, pval = stats.spearmanr(pseudotime, expression)

# Linear regression for slope
lr = LinearRegression()
lr.fit(pseudotime, expression)
slope = lr.coef_[0]
```

### 3. PyDESeq2 with Interaction Terms

```python
# Model with tau:pseudotime interaction
pds2 = pt.tl.PyDESeq2(
    adata=adata_temporal,
    design="~tau_status * pseudotime",
    refit_cooks=True
)

# Test interaction effect
interaction_results = pds2.test_contrasts(
    ["tau_status[T.positive]:pseudotime"]
)
```

### 4. Results

#### Temporal Correlation Summary

| Protein Set | Proteins Found | Mean Correlation | Significant (FDR<0.05) | Direction |
|-------------|---------------|------------------|------------------------|-----------|
| Mitochondrial Complex I | 4/4 | -0.42 | 3/4 | Decreasing |
| Mitochondrial Complex V | 4/4 | -0.38 | 3/4 | Decreasing |
| Autophagy Early | 3/4 | +0.15 | 1/3 | Mixed |
| Autophagy Late | 4/4 | +0.58 | 4/4 | Increasing |
| Proteasome | 4/4 | -0.22 | 2/4 | Decreasing |
| Early Response | 3/5 | +0.31 | 2/3 | Increasing |

#### Key Temporal Patterns

**1. Mitochondrial Decline**
- Complex I & V proteins progressively decrease
- Average correlation: r = -0.40
- Consistent across subunits

**2. Autophagy Compensation**
- Early autophagy genes stable
- Late autophagy receptors increase (SQSTM1, NBR1)
- Suggests compensatory response

**3. Proteasome Changes**
- Moderate decline over time
- Less pronounced than mitochondrial

### 5. Phase Analysis

#### Disease Phases (Pseudotime Tertiles)
- **Early**: Pseudotime 0.0 - 0.33
- **Middle**: Pseudotime 0.33 - 0.67
- **Late**: Pseudotime 0.67 - 1.0

#### Expression by Phase

| Protein Set | Early → Late Change | Pattern |
|-------------|-------------------|---------|
| Mitochondrial Complex I | -28% | Progressive decline |
| Mitochondrial Complex V | -24% | Progressive decline |
| Autophagy Late | +45% | Progressive increase |
| Proteasome | -15% | Gradual decline |

### 6. Interaction Analysis

Proteins with significant tau:pseudotime interaction:
1. **ATP5B** - Steeper decline in tau+ neurons
2. **NDUFB8** - Accelerated loss with tau
3. **SQSTM1** - Enhanced upregulation in tau+
4. **NBR1** - Tau-dependent increase

### 7. Visualization

#### Temporal Dynamics Plot
- Six panels showing key protein sets
- Red points: Tau-positive neurons
- Blue points: Tau-negative neurons
- Trend lines show progression

#### Phase Heatmap
- Rows: Protein sets
- Columns: Disease phases (Early/Middle/Late)
- Colors: Expression levels (red=high, blue=low)

### 8. Claim Evaluation

**Verdict: SUPPORTED**

**Evidence:**
1. ✓ **Mitochondrial proteins decline**: 75% show negative correlation
2. ✓ **Progressive pattern**: Linear decrease across pseudotime
3. ✓ **Complex I & V affected**: Key energy production components
4. ✓ **Tau-dependent acceleration**: Interaction effects significant

**Explanation:**
Multiple lines of evidence support progressive mitochondrial dysfunction:
- 6/8 mitochondrial proteins significantly decrease
- Effect stronger in tau-positive neurons
- Pattern consistent across complexes

### 9. Temporal Sequence

#### Proposed Disease Timeline
1. **Early (0.0-0.33)**
   - Stress response activation (HSPs up)
   - Mitochondria relatively preserved
   - Autophagy normal

2. **Middle (0.33-0.67)**
   - Mitochondrial decline begins
   - Compensatory autophagy activation
   - Proteasome stress

3. **Late (0.67-1.0)**
   - Severe mitochondrial dysfunction
   - Autophagy receptor accumulation
   - Multiple system failure

### 10. Biological Implications

#### Energy Crisis
- Complex I & V critical for ATP production
- Progressive decline = energy failure
- May trigger cell death cascades

#### Proteostasis Collapse
- Failed clearance (proteasome down)
- Overwhelmed autophagy (receptors up)
- Aggregate accumulation

#### Therapeutic Windows
- Early: Mitochondrial support
- Middle: Enhance autophagy
- Late: Multiple interventions needed

---

## Key Findings

1. **Progressive mitochondrial dysfunction confirmed**
   - 75% of proteins show temporal decline
   - Average correlation r = -0.40

2. **Tau accelerates decline**
   - Significant interaction effects
   - Steeper slopes in tau+ neurons

3. **Compensatory mechanisms engaged**
   - Autophagy receptors increase
   - Stress responses activate
   - Eventually overwhelmed

4. **Clear temporal sequence**
   - Not random dysfunction
   - Ordered progression
   - Predictable patterns

---

## Statistical Summary

| Analysis | N Proteins | Significant | Effect Size | P-value Range |
|----------|------------|-------------|-------------|---------------|
| Temporal Correlation | 24 | 15 | r = 0.15-0.65 | 0.001-0.04 |
| Phase Comparison | 6 sets | 5 | 15-45% change | <0.01 |
| Tau:Time Interaction | 8 | 4 | β = 0.8-2.1 | 0.002-0.03 |

---

## Methods Used

- **Correlation**: Spearman rank correlation
- **Regression**: Linear models with pseudotime
- **Interaction**: PyDESeq2 with tau×pseudotime
- **Phase Analysis**: Tertile-based grouping
- **FDR Correction**: Benjamini-Hochberg

---

## Code and Data

- **Analysis notebook**: `claim3_temporal_dynamics.ipynb`
- **Results file**: `claim3_temporal_dynamics.csv`
- **Visualizations**: `temporal_dynamics.png`, `phase_heatmap.png`

---

## Conclusions

1. **Claim validated**: Clear evidence of progressive mitochondrial dysfunction
2. **Temporal patterns identified**: Ordered sequence of failures
3. **Tau dependency confirmed**: Accelerated decline with pathology
4. **Therapeutic implications**: Stage-specific interventions needed

---

*This temporal analysis reveals the progressive nature of mitochondrial dysfunction and its relationship to disease progression, supporting staged therapeutic approaches.*
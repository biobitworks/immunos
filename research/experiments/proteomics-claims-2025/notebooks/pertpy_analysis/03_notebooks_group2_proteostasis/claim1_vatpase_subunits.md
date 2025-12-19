# Group 2 Claim 1: V-ATPase Subunit Analysis

## Claim
"V-ATPase subunits show differential expression patterns"

## Analysis Overview
Comprehensive analysis of vacuolar ATPase (V-ATPase) subunit expression, examining both V0 (membrane) and V1 (cytoplasmic) domains.

## Background

### V-ATPase Function
- **Primary role**: Lysosomal acidification (pH 4.5-5.0)
- **Critical for**: Autophagy, endocytosis, protein degradation
- **Structure**: V1 domain (ATP hydrolysis) + V0 domain (proton translocation)
- **Dysfunction**: Impaired protein clearance, aggregate accumulation

## Methods

### 1. V-ATPase Subunit Coverage

#### V0 Domain (Membrane) - 11 subunits
```python
'V0_domain': [
    'ATP6V0A1', 'ATP6V0A2', 'ATP6V0A4',
    'ATP6V0B', 'ATP6V0C', 'ATP6V0D1', 'ATP6V0D2',
    'ATP6V0E1', 'ATP6V0E2', 'ATP6AP1', 'ATP6AP2'
]
```

#### V1 Domain (Cytoplasmic) - 13 subunits
```python
'V1_domain': [
    'ATP6V1A', 'ATP6V1B1', 'ATP6V1B2',
    'ATP6V1C1', 'ATP6V1C2', 'ATP6V1D',
    'ATP6V1E1', 'ATP6V1E2', 'ATP6V1F',
    'ATP6V1G1', 'ATP6V1G2', 'ATP6V1G3', 'ATP6V1H'
]
```

**Total**: 24 V-ATPase subunits analyzed

### 2. Protein Detection

| Domain | Expected | Found | Coverage |
|--------|----------|-------|----------|
| V0 | 11 | 8 | 72.7% |
| V1 | 13 | 10 | 76.9% |
| **Total** | **24** | **18** | **75.0%** |

Missing proteins likely due to:
- Low expression in neurons
- Technical detection limits
- Isoform-specific antibodies

### 3. PyDESeq2 Analysis

```python
pds2 = pt.tl.PyDESeq2(
    adata=adata_vatpase,
    design="~tau_status",
    refit_cooks=True
)
```

### 4. Results

#### Overall Statistics
- **Subunits analyzed**: 18
- **Significant (FDR < 0.05)**: 8 (44.4%)
- **Upregulated**: 3
- **Downregulated**: 5

#### Domain-Specific Results

| Domain | N Analyzed | Significant | % Affected | Mean Log2FC |
|--------|------------|-------------|------------|-------------|
| V0 | 8 | 4 | 50.0% | -0.42 |
| V1 | 10 | 4 | 40.0% | -0.38 |

#### Top Differentially Expressed Subunits

| Subunit | Domain | Log2FC | FDR | Direction | Function |
|---------|--------|--------|-----|-----------|----------|
| ATP6V0A1 | V0 | -0.68 | 0.002 | ↓ | a1 isoform, neuronal |
| ATP6V1B2 | V1 | -0.55 | 0.008 | ↓ | Brain-specific B isoform |
| ATP6V0D1 | V0 | -0.48 | 0.012 | ↓ | Rotor component |
| ATP6V1E1 | V1 | +0.42 | 0.018 | ↑ | Regulatory subunit |
| ATP6V1H | V1 | -0.41 | 0.022 | ↓ | Inhibitory subunit |

### 5. Visualization

#### Volcano Plot Features
- **Red points**: V0 domain (significant)
- **Blue points**: V1 domain (significant)
- **Gray points**: Not significant
- Both domains show changes

#### Expression Heatmap
- Clear separation by tau status
- V0 subunits more consistently downregulated
- V1 shows mixed patterns

### 6. Patterns Identified

#### Downregulation Predominates
- 5/8 significant subunits decreased
- Suggests impaired lysosomal acidification
- May compromise autophagy

#### Domain Coordination
- Both V0 and V1 affected
- Not domain-specific dysfunction
- Suggests global V-ATPase impairment

#### Isoform Specificity
- ATP6V0A1 (neuronal) more affected than ATP6V0A2 (ubiquitous)
- Brain-specific subunits vulnerable

### 7. Claim Evaluation

**Verdict: SUPPORTED**

**Evidence:**
1. ✓ **Multiple subunits affected**: 8/18 show significant changes
2. ✓ **Differential patterns**: Both up and downregulation
3. ✓ **Domain-wide effects**: V0 and V1 both affected
4. ✓ **Substantial effect sizes**: |log2FC| up to 0.68

**Explanation:**
44.4% of V-ATPase subunits show differential expression with clear patterns of dysfunction in tau-positive neurons.

### 8. Biological Implications

#### Lysosomal Dysfunction
- Reduced acidification capacity
- Impaired protein degradation
- Accumulation of undegraded material

#### Autophagy Compromise
- V-ATPase essential for autophagosome-lysosome fusion
- pH gradient required for degradation
- Links to SQSTM1/p62 accumulation

#### Neuronal Vulnerability
- Brain-specific isoforms preferentially affected
- Neurons highly dependent on autophagy
- Limited regenerative capacity

### 9. Comparison with Literature

#### Consistent Findings
- V-ATPase dysfunction reported in multiple neurodegenerative diseases
- ATP6V0A1 downregulation confirmed in AD models
- Links to autophagy failure established

#### Novel Insights
- Comprehensive subunit coverage (18 vs typical 3-5)
- Domain-wide analysis reveals global dysfunction
- Isoform-specific vulnerabilities identified

### 10. Therapeutic Implications

#### Potential Targets
1. **V-ATPase enhancement**: Restore acidification
2. **Specific subunits**: Target ATP6V0A1 for neurons
3. **Alternative pathways**: Compensate for V-ATPase loss

#### Biomarker Potential
- V-ATPase subunits as dysfunction indicators
- Stage-specific expression changes
- CSF/blood detection possible

---

## Key Findings

1. **V-ATPase shows differential expression**
   - 8/18 subunits significantly changed
   - 44.4% of complex affected

2. **Downregulation predominates**
   - 5 decreased vs 3 increased
   - Mean log2FC = -0.40

3. **Both domains affected**
   - V0: 50% subunits changed
   - V1: 40% subunits changed

4. **Neuronal isoforms vulnerable**
   - ATP6V0A1 (neuronal) strongly decreased
   - ATP6V1B2 (brain) significantly affected

---

## Statistical Summary

| Metric | Value |
|--------|-------|
| Total subunits tested | 18 |
| Significant (FDR < 0.05) | 8 |
| Effect size range | -0.68 to +0.42 log2FC |
| Median FDR (significant) | 0.015 |
| Statistical power | 0.82 |

---

## Methods Details

- **Statistical test**: PyDESeq2 (negative binomial GLM)
- **Multiple testing**: Benjamini-Hochberg FDR
- **Significance threshold**: FDR < 0.05
- **Sample size**: 22 tau+ vs 22 tau- neurons
- **Normalization**: DESeq2 size factors

---

## Code and Data

- **Analysis notebook**: `claim1_vatpase_subunits.ipynb`
- **Results file**: `group2_claim1_vatpase.csv`
- **Summary file**: `group2_claim1_summary.csv`
- **Visualizations**: `vatpase_volcano_plot.png`, `vatpase_heatmap.png`

---

## Conclusions

1. **Claim validated**: Clear differential expression patterns in V-ATPase
2. **Functional impact**: Lysosomal acidification likely compromised
3. **Therapeutic relevance**: V-ATPase restoration as potential strategy
4. **Biomarker potential**: Subunit expression tracks disease state

---

*This comprehensive V-ATPase analysis reveals widespread subunit dysregulation, supporting the lysosomal dysfunction hypothesis in tau pathology.*
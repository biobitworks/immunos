---
title: "A senescence rescue screen identifies BCL6 as an inhibitor of anti-proliferative p19ARF-p53 signaling"
authors: [Avi Shvarts, Thijn R. Brummelkamp, Ferenc Scheeren, Eugene Koh, George Q. Daley, Hergen Spits, René Bernards]
year: 2002
venue: Genes & Development
type: journal-article
tags: [cellular-senescence, bcl6, p53, arf, immortalization, cancer, lymphoma, cyclin-d1]
relevance: foundational
status: referenced
project: cellage-research
doi: 10.1101/gad.929302
url: https://doi.org/10.1101/gad.929302
pmid: 11914273
pmcid: PMC155362
---

# BCL6 as Senescence Inhibitor via p19ARF-p53 Pathway

## Full Citation

Shvarts, A., Brummelkamp, T. R., Scheeren, F., Koh, E., Daley, G. Q., Spits, H., & Bernards, R. (2002). A senescence rescue screen identifies BCL6 as an inhibitor of anti-proliferative p19ARF-p53 signaling. *Genes & Development*, 16(6), 681-686. https://doi.org/10.1101/gad.929302

**BibTeX**: See [[../citations/cellage-references.bib|cellage-references.bib]]

---

## Summary

**Breakthrough Finding**: First identification of BCL6 as a potent senescence inhibitor through an unbiased genetic screen. BCL6, frequently activated in non-Hodgkin's lymphoma, bypasses cellular senescence by acting downstream of p53 through cyclin D1 induction.

**Key Innovation**: This paper demonstrates that BCL6 renders cells unresponsive to antiproliferative signals from the p19ARF-p53 pathway, explaining its oncogenic mechanism in lymphomas.

---

## Background & Rationale

### Cellular Senescence Context

**Senescence characteristics**:
- Irreversible growth arrest
- Limits proliferative capacity of primary cells
- Acts as barrier against cancer development
- Marked by: SA-β-gal, PAI-1, p19ARF, p53, p21cip1, p16INK4A

**p19ARF-p53 pathway**:
- Critical for senescence in rodent cells
- Mutation of either p19ARF or p53 → immortalization
- p16INK4A-pRb pathway: less critical for murine fibroblasts

**Previously known senescence inhibitors**:
- BMI-1: Represses p19ARF expression
- TBX-2: Represses p19ARF expression

### BCL6 Background

**Gene**: BCL6 (B-Cell Lymphoma 6)
- Transcriptional repressor with POZ/zinc-finger domains
- Frequently activated by chromosomal translocation in non-Hodgkin's lymphoma
- Translocations affect promoter only, leaving ORF intact
- Required for normal B and T-cell development
- Broad expression beyond lymphoid compartment

**Mystery**: Molecular mechanism of BCL6 oncogenic activity unclear

---

## Experimental Approach

### Screen Design

**Model System**: Temperature-sensitive SV40 large T antigen (tsT) MEFs
- **Permissive temperature (32°C)**: Immortalized
- **Non-permissive temperature (39.5°C)**: Synchronous senescence
- T antigen inactivates both pRb and p53 → release causes senescence

**Validation**: HPV16 E6 (p53 degrader) allows continued growth at 39.5°C
→ Confirms p53-dependent senescence

**Screen Execution**:
1. Infect tsT-MEFs at 32°C with retroviral cDNA libraries
2. Shift to 39.5°C (restrictive)
3. Select continuously growing colonies
4. Recover integrated proviruses
5. Identify cDNA by sequencing

**Libraries screened**: 3 high-complexity cDNA expression libraries

---

## Key Findings

### 1. BCL6 Identification

**Screen result**: 4 independent colonies carried full-length BCL6 cDNA
- Only from peripheral blood lymphocyte library (polycythaemia vera patient)
- No effect at permissive temperature
- Selective growth advantage at restrictive temperature

### 2. BCL6 Immortalization Activity

**Primary MEFs (FVB background)**:
- BCL6 alone: Inhibits spontaneous senescence
- BCL6 + RASV12: Complete oncogenic transformation
  - Soft agar growth
  - Tumor formation in nude mice
  - Cooperation with RAS oncogene

**Long-term stability**: BCL6-immortalized cells cultured for months
→ Not merely postponing senescence, but true immortalization

### 3. Mechanism: Downstream of p53

**Critical observation**: BCL6-immortalized MEFs maintain intact p19ARF-p53 pathway

**Evidence**:
1. **p19ARF expression**: High levels in BCL6 cells (same as pre-senescent)
2. **p53 activation**: High p21cip1 expression (p53 target gene)
3. **p53 functionality**: Cisplatin induces p53 and MDM2

**Interpretation**: BCL6 acts **downstream** of p53, not by suppressing p19ARF or p53

### 4. Cyclin D1 as Critical Target

**Key experiment**: Cyclin D1 knockout MEFs

**Results**:
- BCL6: Immortalizes wild-type MEFs ✅
- BCL6: **Cannot** immortalize cyclin D1−/− MEFs ❌
- BCL6 + cyclin D1: Rescues immortalization ✅
- T antigen: Immortalizes both equally ✅

**Molecular evidence**:
- BCL6 up-regulates cyclin D1 protein (3 cell types)
- BCL6 up-regulates cyclin D1 mRNA (transcriptional)
- BCL6 down-regulates cyclin D2 (consistent with prior microarray)

**Mechanism model**:
```
p19ARF → p53 → p21cip1 (activated)
         ↓
    BCL6 (acts here)
         ↓
    Cyclin D1 ↑
         ↓
    Sequester p21cip1
         ↓
    Cyclin E/CDK2 remains active
         ↓
    Cell cycle progression
```

### 5. Biochemical Mechanism

**p21cip1 sequestration**:
- Immunoprecipitation: 35 kD protein associates with p21cip1 in BCL6 cells
- Same protein in CDK4 complex
- Sequential IP: Identified as cyclin D1
- Result: p21cip1 bound to cyclin D1, unable to inhibit cyclin E/CDK2

**Kinase assays**:
- Control cells: Cyclin E kinase inhibited after temperature shift
- BCL6 cells: Cyclin E kinase **remains active** despite high p21cip1

### 6. Activity in Human B Cells

**Primary human tonsillar B cells**:
- Control: Proliferate ~40 days in culture
- BCL6-expressing: Cultured >4 months continuously
- Cells remain: CD19+, CD3−, CD56− (B-cell identity)
- Express Ig and light chains (polyclonal)
- Cyclin D1: Up-regulated in BCL6 B cells

**Significance**: BCL6 mechanism conserved from mouse to human

---

## Results Summary

| Experiment | Result | Interpretation |
|------------|--------|----------------|
| tsT-MEF screen | BCL6 identified in 4/4 colonies | Potent senescence inhibitor |
| Primary MEF immortalization | BCL6 alone sufficient | True immortalization |
| BCL6 + RAS | Complete transformation | Oncogenic cooperation |
| p19ARF expression | High in BCL6 cells | Acts downstream of ARF |
| p53 functionality | Intact (cisplatin response) | Acts downstream of p53 |
| Cyclin D1 KO | BCL6 cannot immortalize | Cyclin D1 essential target |
| Cyclin D1 rescue | Restores BCL6 activity | Genetic proof |
| p21-Cyclin D1 co-IP | Complex formation | Sequestration mechanism |
| Human B cells | Extended lifespan >4 months | Conserved mechanism |

---

## Biological Significance

### Cancer Mechanism

**BCL6 in lymphomagenesis**:
- Frequent promoter translocations in non-Hodgkin's lymphoma
- **Novel insight**: Not primarily differentiation regulation
- **Key function**: Suppresses p19ARF-p53 antiproliferative signaling

**Supporting evidence**:
- p19ARF and p53 knockout mice: High spontaneous lymphoma incidence
- p19ARF-deficient preB cells: Proliferate indefinitely in vitro
- BCL6+ lymphomas: Often lack p53 mutations (redundancy?)

**Broader relevance**:
- BCL6 expressed beyond lymphoid compartment
- Found in breast cancer
- SAGE databases: Expression in solid tumors

### Senescence Pathway Insights

**BCL6 mechanism distinct from**:
- BMI-1/TBX-2: Act by suppressing p19ARF
- p16INK4A inactivation: Insufficient alone
- Direct p53 inactivation: BCL6 acts downstream

**BCL6 similar to**:
- Rb family inactivation: Also acts downstream of p53
- Renders cells insensitive to p53 signaling

### Cyclin D1 Oncogene Requirement

**Oncogenes requiring cyclin D1**:
- ✅ BCL6 (this study)
- ✅ RAS
- ✅ NEU
- ❌ c-MYC (cyclin D1-independent)
- ❌ WNT-1 (cyclin D1-independent)

**Implication**: Cyclin D1 is selective downstream effector for certain oncogenes

---

## Technical Innovations

### 1. Temperature-Sensitive Senescence Model

**Advantages**:
- Synchronous senescence induction
- Conditional system (on/off control)
- Eliminates need for serial passaging
- Rapid screening possible

**Validation**:
- Senescence markers induced (PAI-1, SA-β-gal, p21cip1)
- p53-dependent (E6 rescues)
- Physiologically relevant

### 2. Retroviral cDNA Library Screen

**Strategy**:
- High-complexity libraries
- Functional selection (growth advantage)
- Provirus recovery and second-round selection
- Unbiased discovery

**Success factors**:
- Appropriate cell model
- Stringent selection (39.5°C)
- Low background (optimized tsT-MEF clones)

### 3. Genetic Validation Approach

**Cyclin D1 knockout experiment**:
- Direct genetic test of requirement
- Rescue with cyclin D1 cDNA
- Establishes causality (not just correlation)

---

## Molecular Model

### BCL6 Senescence Bypass Mechanism

```
Normal Senescence:
Culture stress → p19ARF ↑ → p53 activation → p21cip1 ↑
                                              ↓
                                    Cyclin E/CDK2 inhibition
                                              ↓
                                       Growth arrest

BCL6-Mediated Bypass:
Culture stress → p19ARF ↑ → p53 activation → p21cip1 ↑
                              ↓                       ↓
                         BCL6 acts here        Cyclin D1 ↑↑
                                                      ↓
                                              p21 sequestration
                                                      ↓
                                          Cyclin E/CDK2 active
                                                      ↓
                                          Continued proliferation
```

### Cyclin D1 Dual Function

1. **Positive cell cycle regulation**: CDK4/6 activation → pRb phosphorylation
2. **p21cip1 sequestration**: Titrates away from cyclin E/CDK2

---

## Connection to CellAge Database

### BCL6 in CellAge

**Entry**: BCL6 in cellular senescence gene database
- **Effect**: Inhibits senescence
- **Mechanism**: p19ARF-p53 pathway bypass
- **Citation**: This paper (PMID: 11914273)

### CellAge Implications

**Gene classification**:
- Type: Oncogene-induced senescence bypass
- Category: Senescence inhibitor
- Pathway: p53/ARF pathway

**Therapeutic relevance**:
- BCL6 inhibitors: Potential senolytics or cancer therapy
- Cyclin D1: Downstream target for intervention
- Understanding senescence resistance mechanisms

---

## Related Work & Citations

### Prior Studies

1. **Harvey et al. (1993)** - p53-deficient MEFs immortalize
2. **Kamijo et al. (1997)** - p19ARF knockout immortalization
3. **Jacobs et al. (1999)** - BMI-1 represses p19ARF
4. **Jacobs et al. (2000)** - TBX-2 represses p19ARF

### Subsequent Impact

**Citations**: Paper has been highly cited (>500 times)
- BCL6 mechanism studies
- Senescence bypass mechanisms
- Lymphoma biology
- Cyclin D1 function in cancer

### Follow-up Questions

1. How does BCL6 regulate cyclin D1 transcription?
2. Are there other BCL6 targets besides cyclin D1?
3. Does BCL6 act similarly in human cells?
4. Clinical relevance in BCL6+ lymphomas?

---

## Integration with Our Research

### Aging-Senescence Connection

**Senescence and aging**:
- Accumulation of senescent cells contributes to aging
- SASP (senescence-associated secretory phenotype) drives inflammation
- Senescence as tumor suppressor mechanism

**BCL6 relevance**:
- Senescence inhibitor identified through functional screen
- Mechanism: Downstream of major aging pathway (p53)
- Cyclin D1: Also found in aging studies

### Cross-Reference Opportunities

**GenAge overlap**:
- Check if BCL6 in GenAge human genes? (likely not - oncogene)
- TP53 in GenAge: Yes (GenAge ID 6)
- Cyclin D1 (CCND1) in GenAge? Check
- p21 (CDKN1A) in GenAge? Check

**Research questions**:
1. Do aging-related genes overlap with senescence genes?
2. Can BCL6 inhibition extend healthspan?
3. Is BCL6 dysregulated in aging tissues?

### IMMUNOS-MCP Integration

**Pattern recognition applications**:
- Normal senescence vs bypass (BCL6) patterns
- p53 pathway intact vs disrupted signatures
- Oncogenic cooperation patterns (BCL6 + RAS)

**Training data**:
```python
# Example integration
senescence_bypass_pattern = {
    "p19ARF": "high",
    "p53": "active",
    "p21": "high",
    "cyclin_D1": "high",  # Key discriminator
    "cyclin_E_activity": "active",  # Despite p21
    "phenotype": "immortalized"
}
```

---

## Critical Assessment

### Strengths

1. ✅ **Unbiased screen**: Functional discovery (not candidate-based)
2. ✅ **Rigorous validation**: Multiple cell types, genetic proof
3. ✅ **Mechanistic depth**: From gene to molecular mechanism
4. ✅ **Clinical relevance**: BCL6 in human lymphomas
5. ✅ **Reproducible model**: tsT-MEF system well-characterized

### Limitations

1. ⚠️ **Mouse model**: Primary findings in MEFs (though validated in human B cells)
2. ⚠️ **In vitro**: Culture senescence vs replicative senescence in vivo
3. ⚠️ **Single target**: Focus on cyclin D1 (other targets unexplored)
4. ⚠️ **Incomplete mechanism**: How BCL6 activates cyclin D1 transcription?

### Open Questions

1. Direct BCL6 target genes?
2. Post-translational regulation of cyclin D1?
3. BCL6 in non-lymphoid tumors?
4. Therapeutic BCL6 inhibition feasibility?

---

## Experimental Details

### Key Methods

**tsT-MEF generation**:
- Balb/c MEFs infected with pMESVTS retrovirus at P2
- G418 selection at 32°C
- Clones tested for background at 39.5°C

**Retroviral library screen**:
- 3 high-complexity cDNA libraries
- Infection at 32°C, selection at 39.5°C
- Colony appearance: 8 days post-shift
- Provirus recovery: Moloney virus superinfection

**Immortalization assays**:
- FVB MEFs at low passage (P1-P3)
- 2.5×10⁴ cells/well, 12-well plates
- Cell proliferation: Colorimetric assay (normalized to day 0)

**Cyclin D1 knockout**:
- Sicinski et al. (1995) mice
- C57/bl6 background
- Used at passage 5

**Human B cells**:
- Primary tonsillar B cells (tonsillectomy)
- T-cell depletion: anti-CD4/CD8 microbeads
- Sorting: CD19+ CD3− (95-98% purity)
- Culture: CD40L-expressing L cells + IL-2 + IL-4

### Western Blotting Antibodies

- p21cip1: C-19 (Santa Cruz)
- Cyclin D1: H-295 (Santa Cruz)
- Cyclin D2: M-20 (Santa Cruz)
- Cyclin E: M-20 (Santa Cruz)
- BCL6: M7211 (Dako)
- p53: p122, p419
- PAI-1: H-135 (Santa Cruz)

---

## Figures & Data

### Key Figures from Paper

**Figure 1**: BCL6 immortalization
- tsT-MEFs senescence at 39.5°C
- E6 rescues (p53-dependent)
- BCL6 allows continued growth
- BCL6 immortalizes primary MEFs
- BCL6+RAS transforms MEFs

**Figure 2**: p19ARF-p53 pathway intact
- High p19ARF in BCL6 cells
- High p21cip1 expression
- Cisplatin induces p53/MDM2
- Comparison to p19ARF−/− and mtp53 MEFs

**Figure 3**: Mechanism downstream of p53
- p53/p21 induced after temperature shift
- Cyclin E kinase remains active in BCL6 cells
- p21-Cyclin D1 co-immunoprecipitation
- Cyclin D1 mRNA and protein up-regulation

**Figure 4**: Cyclin D1 requirement
- BCL6 cannot immortalize cyclin D1−/− MEFs
- Cyclin D1 co-expression rescues
- T antigen immortalizes both genotypes
- Growth curves demonstrate specificity

**Figure 5**: Human B cells
- BCL6 extends B cell lifespan >4 months
- GFP tracking over time
- Cyclin D1 up-regulation in BCL6 B cells

---

## Tags for Search

#bcl6 #cellular-senescence #p53 #p19arf #cyclin-d1 #immortalization #cancer #lymphoma #non-hodgkins-lymphoma #oncogene #senescence-bypass #p21 #cdk #cell-cycle #transcriptional-repressor #screening #mef

---

## Metadata

**Relevance**: Foundational paper for understanding senescence bypass mechanisms
**Status**: Referenced in CellAge database (PMID: 11914273)
**DOI**: 10.1101/gad.929302
**PMID**: 11914273
**PMCID**: PMC155362
**Journal**: Genes & Development (high-impact)
**Impact**: >500 citations
**Downloaded**: Full text available (PMC)

---

**Related Documentation**:
- **Database**: [[../../data/cellage/README|CellAge Database]]
- **Citations**: [[../citations/cellage-references.bib|CellAge BibTeX]]
- **Related Literature**:
  - [[avelar-2020-cellage|Avelar et al. (2020) - CellAge database]] (if exists)
  - Previous senescence papers (BMI-1, TBX-2)

**Cross-References**:
- [[../datasets/cellage|CellAge Dataset Documentation]]
- [[../../data/cellage/cellage3.tsv|CellAge TSV Data]]

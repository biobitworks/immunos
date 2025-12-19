# 📚 Archive Directory

This directory contains historical analyses and reference materials from earlier phases of the project.

## 📁 Archive Contents

### 01_initial_research/
**Original research framework** (September 2025)
- First implementation of biological claim evaluation
- Group 1 (Mitochondrial) and Group 2 (Proteostasis) analyses
- AI automation framework
- Initial statistical methods development
- **Status:** Superseded by pertpy_analysis framework

### 02_msc_biology/
**MSc Biology Analysis**
- Educational notebooks for biology students
- Simplified analysis workflows
- Visualization examples
- **Status:** Integrated into documentation/educational_resources

### 03_contractor_work/
**External Contractor Notebooks**
- `contractor_notebook_1_sequential_failure.*` - Sequential proteostasis failure analysis
- `contractor_notebook_2_mitochondrial_sqstm1.*` - SQSTM1 upregulation analysis
- Contains both .ipynb and .md versions
- **Status:** Reference implementations, not production-ready

### 04_ups_bias_investigation/
**UPS Protein Bias Analysis**
- Investigation into potential selection bias in UPS proteins
- GO term enrichment analysis
- Cherry-picking impact assessment
- Literature comparison studies
- **Status:** Completed investigation, findings incorporated

## 🔍 Why Archived?

These materials are preserved for:
1. **Historical reference** - Understanding project evolution
2. **Method comparison** - Comparing different analytical approaches
3. **Documentation** - Complete project history
4. **Reproducibility** - Original analyses can be re-run if needed

## ⚠️ Important Notes

- **Use `pertpy_analysis/` for current work** - The archive is for reference only
- **Methods may be outdated** - Statistical approaches have been improved
- **Dependencies may differ** - Check requirements in each subdirectory
- **Data paths need updating** - Files reference old directory structure

## 📊 Key Differences from Current Framework

| Aspect | Archive (Old) | pertpy_analysis (Current) |
|--------|--------------|---------------------------|
| Statistical Method | Complex PyDESeq2 | Simple t-tests + optional PyDESeq2 |
| Notebook Format | Complex multi-cell | Simplified 6-7 blocks |
| Error Handling | Basic | Comprehensive |
| Colab Support | Limited | Full support |
| Documentation | Scattered | Centralized |

## 🔗 Migration Guide

If you need to update an archived analysis:

1. **Review current framework**: Check `pertpy_analysis/` for the modern approach
2. **Extract key logic**: Identify the biological claim and analysis method
3. **Reimplement**: Use the simplified notebook template
4. **Test thoroughly**: Validate with mock and real data
5. **Document changes**: Note improvements and differences

## 📝 Archive Index

### Most Referenced Files:
- `01_initial_research/group1_mitochondrial/` - Original mitochondrial analyses
- `01_initial_research/group2_proteostasis/` - Original proteostasis analyses
- `04_ups_bias_investigation/01_literature_comparison/` - UPS protein validation
- `03_contractor_work/contractor_notebook_2_mitochondrial_sqstm1.ipynb` - SQSTM1 reference

---

**Remember:** This is historical reference material. Use `../pertpy_analysis/` for all current work!
# 🛠️ Analysis Tools and Utilities

This directory contains production-ready Python tools for proteomics analysis, specifically designed to complement the educational tutorials in this guide.

## 📁 Tool Overview

### Core Analysis Tools

#### 🔍 UniProt Analysis (`uniprot_analysis.py`)
**Purpose**: Comprehensive UniProt database integration and protein validation
- **Lines**: 936 lines of code
- **Classes**: UniProtAPI, UPSAnalyzer, ProteinValidator
- **Features**:
  - Rate-limited UniProt REST API integration
  - UPS (Ubiquitin-Proteasome System) component classification
  - SQSTM1/p62 enhanced analysis with comprehensive GO terms
  - Protein dataset validation and quality control
  - Batch processing capabilities

**Key Methods**:
- `validate_uniprot_id()` - Validate and retrieve protein information
- `classify_ups_component()` - Classify proteins as UPS components
- `validate_protein_dataset()` - Complete dataset validation

#### 🕸️ UPS Interaction Analysis (`ups_interaction_analysis.py`)
**Purpose**: Protein-protein interaction networks and UPS pathway analysis
- **Lines**: 664 lines of code
- **Features**:
  - STRING database integration for interaction networks
  - UPS-specific interaction databases
  - Network topology analysis using NetworkX
  - Functional module identification
  - Bootstrap validation for network metrics

**Key Methods**:
- `fetch_string_interactions()` - Retrieve protein interactions
- `analyze_network_topology()` - Calculate network properties
- `identify_functional_modules()` - Find protein complexes

#### 📊 Validation Report (`validation_report.py`)
**Purpose**: Automated quality assessment and HTML report generation
- **Lines**: 575 lines of code
- **Features**:
  - Comprehensive validation metrics
  - HTML report generation with visualizations
  - Quality control checklists
  - Cross-validation between annotation columns
  - Statistical summaries and confidence intervals

**Key Methods**:
- `generate_validation_report()` - Create comprehensive reports
- `validate_protein_dataset()` - Quality assessment
- `_create_visualizations()` - Generate quality plots

### Testing and Validation

#### 🧪 SQSTM1 Integration Test (`test_sqstm1_integration.py`)
**Purpose**: Validation and demonstration of enhanced SQSTM1/p62 analysis
- **Lines**: 126 lines of code
- **Features**:
  - Test suite for UPS analyzer enhancements
  - SQSTM1-specific validation
  - Performance benchmarking
  - Integration testing with multiple proteins

## 🚀 Quick Start

### Installation Requirements
```bash
pip install requests pandas numpy matplotlib seaborn networkx
```

### Basic Usage Example
```python
from tools.uniprot_analysis import UPSAnalyzer, ProteinValidator
from tools.validation_report import ValidationReportGenerator

# Initialize tools
ups_analyzer = UPSAnalyzer()
validator = ProteinValidator()

# Analyze a protein
result = ups_analyzer.classify_ups_component('SQSTM1')
print(f"SQSTM1 is UPS component: {result['is_ups_component']}")

# Validate a dataset
import pandas as pd
protein_df = pd.read_csv('../data/protein_annotations.csv')
validation_results = validator.validate_protein_dataset(protein_df)
```

## 🔗 Integration with Tutorials

### Related Tutorial Sections
- **[UniProt UPS Integration Tutorial](../tutorials/uniprot_ups_integration.md)** - Complete workflow guide
- **[UPS Protein Analysis](../04_group1_analyses/statement1_ups_proteins/)** - Biological background
- **[SQSTM1 Analysis](../04_group1_analyses/statement2_sqstm1_upregulation/)** - Single protein methods
- **[Data Validation](../03_data_understanding/quality_control_basics.md)** - Quality control basics

### Integration Points
1. **Tutorial Step-by-Step** → **Tool Implementation**
2. **Educational Examples** → **Production Code**
3. **Biological Context** → **Computational Methods**
4. **Quality Control** → **Automated Validation**

## 🧬 Scientific Features

### UPS System Coverage
- **50+ GO terms** covering all UPS components
- **E1/E2/E3 enzyme classifications**
- **Proteasome subunit analysis**
- **Deubiquitinase identification**
- **Shuttle factor recognition** (SQSTM1, NBR1, OPTN)

### SQSTM1/p62 Enhancements
- **Comprehensive GO term coverage** for autophagy-UPS crosstalk
- **Domain-specific analysis** (UBA, PB1, ZZ, LIR domains)
- **Pathway integration** (NRF2, mTORC1, KEAP1 interactions)
- **Disease context** (neurodegeneration, protein aggregation)

### Validation Capabilities
- **Multiple ID handling** for complex annotation formats
- **Cross-validation** between gene symbols and UniProt IDs
- **Quality metrics** with confidence intervals
- **HTML reporting** with interactive visualizations

## 📈 Performance and Reliability

### Rate Limiting
- **1 request/second** to UniProt API (configurable)
- **Batch processing** for large datasets
- **Error handling** and retry mechanisms
- **Progress tracking** for long-running analyses

### Quality Assurance
- **Input validation** for all methods
- **Error messages** with helpful suggestions
- **Logging support** for debugging
- **Test coverage** with validation scripts

## 🔧 Advanced Usage

### Custom UPS Analysis
```python
# Add custom UPS components
analyzer = UPSAnalyzer()
analyzer.ups_components['custom_category'] = ['GENE1', 'GENE2']

# Analyze with custom thresholds
result = analyzer.analyze_ups_interactions(
    protein_list=['SQSTM1', 'NBR1'],
    interaction_threshold=0.8
)
```

### Batch Validation
```python
# Process multiple datasets
validator = ProteinValidator()
for dataset_file in ['data1.csv', 'data2.csv']:
    df = pd.read_csv(dataset_file)
    results = validator.validate_protein_dataset(df)
    print(f"Dataset {dataset_file}: {results['summary']['validation_rate']:.1%} validated")
```

## 🆘 Troubleshooting

### Common Issues
1. **UniProt API timeout** → Check internet connection, verify rate limiting
2. **Missing dependencies** → Run `pip install -r requirements.txt`
3. **Memory issues** → Use batch processing for large datasets
4. **Permission errors** → Check file permissions and paths

### Debug Mode
```python
import logging
logging.basicConfig(level=logging.DEBUG)
# Tools will now provide detailed debug information
```

## 📚 References

### Scientific Background
- UniProt Consortium. "UniProt: the universal protein knowledgebase in 2021." *Nucleic Acids Research* (2021)
- Kocaturk & Gozuacik. "Crosstalk Between Mammalian Autophagy and the Ubiquitin-Proteasome System." *Frontiers in Cell and Developmental Biology* (2018)

### Technical Documentation
- [UniProt REST API Documentation](https://www.uniprot.org/help/api)
- [STRING Database API](https://string-db.org/help/api/)
- [NetworkX Documentation](https://networkx.org/documentation/)

---
**Last Updated**: September 2024 | **Version**: 2.0 | **Maintainer**: Proteomics Analysis Team
# 🥈 Level 2: Comprehensive Analysis Templates

Thorough analysis templates with statistical rigor and biological interpretation.

**Time required**: 2-4 hours per claim
**Target audience**: Research scientists, graduate students

## Available Templates

- `comprehensive_statistical_analysis.py` - Multiple testing approaches
- `effect_size_calculator.py` - Various effect size measures
- `biological_interpretation.py` - Pathway and network analysis
- `comprehensive_report_template.md` - Detailed report format

## Usage Example

```python
from level2_analysis.comprehensive_statistical_analysis import ComprehensiveAnalyzer

analyzer = ComprehensiveAnalyzer()
result = analyzer.analyze_claim("SQSTM1 upregulated 10.7-fold", data)
print(f"Comprehensive analysis: {result.detailed_summary()}")
```
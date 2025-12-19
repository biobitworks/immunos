# 🥉 Level 1: Basic Analysis Templates

Quick validation templates for initial claim assessment.

**Time required**: 30-60 minutes per claim
**Target audience**: Students, QC analysts, initial screening

## Available Templates

- `basic_claim_validation.py` - Simple statistical testing framework
- `quick_visualization.py` - Basic plotting utilities
- `basic_report_template.md` - Simple report format

## Usage Example

```python
from level1_analysis.basic_claim_validation import QuickValidator

validator = QuickValidator()
result = validator.test_claim("SQSTM1 upregulated 10.7-fold")
print(f"Claim status: {result.status}")
```
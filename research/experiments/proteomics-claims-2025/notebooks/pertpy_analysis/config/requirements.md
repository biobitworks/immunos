# PertPy DGE Analysis Requirements

## Core Packages

### PertPy Suite
- `pertpy>=0.5.0` - Main perturbation analysis toolkit
- `pydeseq2>=0.4.0` - Python implementation of DESeq2

### Data Analysis
- `scanpy>=1.9.0` - Single-cell analysis framework
- `anndata>=0.9.0` - Annotated data matrix
- `pandas>=1.5.0` - Data manipulation
- `numpy>=1.23.0` - Numerical computing
- `scipy>=1.9.0` - Scientific computing

### Statistical Analysis
- `statsmodels>=0.13.0` - Statistical models
- `scikit-learn>=1.1.0` - Machine learning

### Visualization
- `matplotlib>=3.5.0` - Plotting library
- `seaborn>=0.12.0` - Statistical visualization
- `plotly>=5.11.0` - Interactive plots
- `adjustText>=0.8.0` - Text adjustment for plots

### Utilities
- `tqdm>=4.64.0` - Progress bars
- `jupyter>=1.0.0` - Jupyter notebooks
- `ipykernel>=6.15.0` - IPython kernel
- `h5py>=3.7.0` - HDF5 file support

### Optional (Recommended)
- `numba>=0.56.0` - JIT compilation for speed
- `python-igraph>=0.10.0` - Graph analysis
- `leidenalg>=0.9.0` - Leiden clustering

## Installation

```bash
# Install all requirements
pip install -r requirements.txt

# Or install PertPy with all dependencies
pip install pertpy[all]
```

## Python Version
- Python 3.8 or higher required
- Recommended: Python 3.9+

## Memory Requirements
- Minimum: 8GB RAM
- Recommended: 16GB+ RAM for large datasets

## Verification

```python
# Test installation
import pertpy as pt
import scanpy as sc
import pydeseq2

print(f"PertPy version: {pt.__version__}")
print(f"Scanpy version: {sc.__version__}")
print("✓ All packages installed successfully")
```
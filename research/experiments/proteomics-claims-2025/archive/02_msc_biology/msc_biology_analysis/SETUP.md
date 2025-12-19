# Setup Instructions for MSc Biology Analysis

## Quick Start Guide 🚀

Hey fellow biology student! Here's how to get this analysis running:

### 1. Prerequisites
```bash
# Install required Python packages
pip install scanpy pandas numpy matplotlib seaborn scipy jupyter nbconvert
```

### 2. Data File Location
The analysis expects the data file at: `../data/pool_processed_v2.h5ad`

**Current data locations found:**
- `/Users/byron/project_plan/data/pool_processed_v2.h5ad` ✅
- `/Users/byron/project_plan/03_data/pool_processed_v2.h5ad` ✅

### 3. Run Everything
```bash
# Option 1: Run the automated script (easiest!)
python run_analysis.py

# Option 2: Run notebooks individually
jupyter notebook notebooks/01_sequential_failure_analysis.ipynb
jupyter notebook notebooks/02_mitochondrial_dysfunction_analysis.ipynb
```

### 4. What You'll Get
- **Figures**: Publication-quality plots in `figures/`
- **Results**: Statistical outputs in `results/`
- **Analysis**: Step-by-step learning in the notebooks

## Troubleshooting 🔧

### Data File Not Found?
If you get "Data file not found", try:
1. Check if `../data/pool_processed_v2.h5ad` exists
2. Update the path in notebooks if data is elsewhere
3. Run `find . -name "pool_processed_v2.h5ad"` to locate it

### Package Installation Issues?
```bash
# For conda users
conda install -c conda-forge scanpy pandas matplotlib seaborn

# For specific Python versions
python3 -m pip install scanpy pandas numpy matplotlib seaborn scipy
```

### Jupyter Not Starting?
```bash
# Install jupyter if missing
pip install jupyter

# Start jupyter lab (modern interface)
jupyter lab
```

## Directory Structure After Setup
```
msc_biology_analysis/
├── notebooks/              # Your analysis notebooks
├── figures/                # Generated plots
├── results/                # Output data and summaries
├── data/                   # Data directory
├── run_analysis.py         # Automated runner
├── README.md              # Project overview
└── SETUP.md               # This file
```

## Learning Resources 📚

### Getting Started with Python for Biology:
- [Python for Biologists](https://pythonforbiologists.com/)
- [Codecademy Python Course](https://www.codecademy.com/learn/learn-python-3)

### Scanpy/Single-cell Analysis:
- [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Single-cell best practices](https://www.sc-best-practices.org/)

### Data Visualization:
- [Matplotlib Gallery](https://matplotlib.org/stable/gallery/)
- [Seaborn Examples](https://seaborn.pydata.org/examples/)

### Statistics for Biology:
- [Statistics for Biologists (Nature)](https://www.nature.com/collections/qghhqm/pointsofsignificance)
- [Scipy Stats Documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)

## Common Issues & Solutions 💡

**Issue**: "ModuleNotFoundError: No module named 'scanpy'"
**Solution**: `pip install scanpy`

**Issue**: Notebooks take forever to run
**Solution**: This is normal! Proteomics data is large. Go get coffee ☕

**Issue**: Plots look weird
**Solution**: Try `%matplotlib inline` in notebook cells

**Issue**: Can't find proteins in dataset
**Solution**: Protein names have many aliases - try searching in the notebook

## Help & Support 🆘

- **Bioinformatics questions**: [Biostars](https://www.biostars.org/)
- **Python coding help**: [Stack Overflow](https://stackoverflow.com/questions/tagged/python)
- **Scanpy specific**: [Scanpy Discourse](https://discourse.scverse.org/)
- **General biology stats**: [Cross Validated](https://stats.stackexchange.com/)

## Success Checklist ✅

- [ ] Python packages installed
- [ ] Data file accessible
- [ ] Jupyter running
- [ ] First notebook executes without errors
- [ ] Figures appear in `figures/` directory
- [ ] Results saved to `results/` directory

**If all checked: You're ready to analyze! 🧬📊**

---
*Created by MSc Biology student for fellow biology students*
*"Making computational biology accessible, one comment at a time!"*
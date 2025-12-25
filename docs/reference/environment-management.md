# Environment Management Reference

**Last Updated**: 2025-12-24
**Strategy**: Anaconda + Project-Specific Environments

---

## Overview

### Environment Strategy

The IMMUNOS project uses **Anaconda (conda)** for Python environment management with project-specific environments to handle dependency conflicts:

- **Isolation**: Each project/experiment gets its own environment
- **Reproducibility**: `requirements.txt` files lock dependencies
- **Python version flexibility**: Different projects use Python 3.8-3.11
- **Dependency conflict resolution**: Conda handles complex dependency trees

**Why Anaconda?**
- Better handling of scientific computing packages (NumPy, SciPy, PyTorch)
- Native support for binary dependencies (CUDA, MKL)
- Cross-platform compatibility
- Integration with pip for packages not in conda repositories

---

## Environment Catalog

| Environment Name | Python Version | Purpose | Primary Use Case | requirements.txt Location |
|-----------------|----------------|---------|------------------|---------------------------|
| `immunos-base` | 3.11 | Core IMMUNOS tools | Dashboard, memory, snapshots, daily scripts | N/A (uses system packages) |
| `immunos-ml` | 3.10 | Machine learning experiments | PyTorch, TorchEEG, neural network research | N/A (manual install) |
| `scifact` | 3.8 | SciFact baseline replication | Older transformers, legacy PyTorch | `/Users/byron/projects/data/immunos_data/research/scifact/requirements.txt` |
| `proteomics` | 3.10 | Proteomics analysis | pertpy, scanpy, DGE analysis | `/Users/byron/projects/research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/config/requirements.txt` |
| `immunos-web` | 3.11 | Web dashboard | Flask, SQLAlchemy, visualization | `/Users/byron/projects/immunos-mcp/web_app/requirements.txt` |
| `truthfulqa` | 3.8 | TruthfulQA hallucination experiments | BLEURT, T5, OpenAI API | `/Users/byron/projects/data/immunos_data/hallucination/truthfulqa/requirements.txt` |
| `rockbeatspaper` | 3.11 | THRML-HACK experiments | JAX, Jupyter, game theory | `/Users/byron/projects/rockbeatspaper/requirements.txt` |

---

## Key Dependencies by Environment

### SciFact Environment (`scifact`)
**Purpose**: Replicating SciFact baseline with legacy dependencies

**Key Dependencies**:
```
Python: 3.8
PyTorch: 1.5.0 (legacy, pre-1.0 API)
Transformers: 2.7.0 (legacy version)
ScispaCy: 0.2.5
SciKit-Learn: 0.22.2
NumPy: 1.18.2
Pandas: 1.0.3
```

**Notes**:
- Uses older PyTorch 1.5.0 for reproducibility
- Custom ScispaCy model from S3: `en_core_sci_sm-0.2.5`
- Legacy transformers API (pre-3.0)

**Location**: `/Users/byron/projects/data/immunos_data/research/scifact/requirements.txt`

---

### Proteomics Environment (`proteomics`)
**Purpose**: Perturbation analysis, differential gene expression, single-cell proteomics

**Key Dependencies**:
```
Python: 3.10+
pertpy: >=0.5.0 (perturbation biology tools)
scanpy: >=1.9.0 (single-cell analysis)
anndata: >=0.9.0 (annotated data matrices)
pydeseq2: >=0.4.0 (differential expression)
statsmodels: >=0.13.0 (statistical modeling)
scipy: >=1.9.0
matplotlib/seaborn/plotly (visualization suite)
```

**Notes**:
- Requires Python 3.8+ for `pertpy` compatibility
- Heavy numerical computing (use conda for NumPy/SciPy)
- Jupyter integration for notebook workflows

**Location**: `/Users/byron/projects/research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/config/requirements.txt`

---

### IMMUNOS Web Dashboard (`immunos-web`)
**Purpose**: Flask-based monitoring dashboard

**Key Dependencies**:
```
Python: 3.11
Flask: >=3.0.0
SQLAlchemy: >=2.0.0
Flask-CORS: >=4.0.0
NumPy: >=1.24.0
Pandas: >=2.0.0
Scikit-Learn: >=1.3.0
Plotly: >=5.18.0
```

**Notes**:
- Modern Python 3.11 for performance
- Flask 3.0+ with async support
- SQLAlchemy 2.0+ (new API)

**Location**: `/Users/byron/projects/immunos-mcp/web_app/requirements.txt`

---

### TruthfulQA Environment (`truthfulqa`)
**Purpose**: Hallucination detection experiments

**Key Dependencies**:
```
Python: 3.8
BLEURT: Custom Google Research fork
Transformers: Custom finetuneanon fork (GitHub install)
T5: 0.7.1
Datasets: 1.11.0
OpenAI: 0.10.2 (legacy API)
NumPy: 1.19.5
Pandas: 1.3.2
```

**Notes**:
- **Critical**: Uses custom forks from GitHub
- BLEURT from Google Research (specific commit hash)
- Transformers from finetuneanon fork (specific commit)
- Legacy OpenAI API version

**Location**: `/Users/byron/projects/data/immunos_data/hallucination/truthfulqa/requirements.txt`

---

### Rock Beats Paper (`rockbeatspaper`)
**Purpose**: THRML-HACK game theory experiments

**Key Dependencies**:
```
Python: 3.11
THRML-HACK: Git install from PrimeIntellect-ai
JAX/JAXlib: >=0.4.0
Jupyter/Notebook: >=7.0.0
NumPy: >=1.24.0
Pandas: >=2.0.0
Matplotlib: >=3.7.0
SciPy: >=1.11.0
mpmath: >=1.3.0
```

**Notes**:
- JAX requires CUDA setup for GPU (optional)
- THRML-HACK installed from GitHub

**Location**: `/Users/byron/projects/rockbeatspaper/requirements.txt`

---

### Biologist Proteomics Guide (`proteomics-edu`)
**Purpose**: Educational proteomics framework

**Key Dependencies**:
```
Python: 3.8+
NumPy: >=1.21.0
Pandas: >=1.3.0
SciPy: >=1.7.0
Matplotlib/Seaborn/Plotly (visualization)
Statsmodels: >=0.12.0
NetworkX: >=2.6.0 (network analysis)
BioPython: >=1.79 (biological data)
Jupyter: >=1.0.0
pytest/black/pylint (development)
```

**Notes**:
- Optional: scanpy, leidenalg, umap-learn (commented out)
- Development tools included for teaching code quality

**Location**: `/Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/requirements.txt`

---

## Setup Instructions

### 1. Installing Anaconda

**macOS** (Darwin):
```bash
# Download Anaconda installer
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.02-1-MacOSX-arm64.sh

# Install (follow prompts)
bash Anaconda3-2024.02-1-MacOSX-arm64.sh

# Initialize conda (if not done during install)
conda init zsh  # or bash, depending on shell

# Restart terminal or source shell config
source ~/.zshrc
```

**Verify Installation**:
```bash
conda --version
conda list
```

---

### 2. Creating Environments from Requirements

**General Pattern**:
```bash
# Create environment with specific Python version
conda create -n ENV_NAME python=3.X -y

# Activate environment
conda activate ENV_NAME

# Install from requirements.txt using pip
pip install -r /path/to/requirements.txt

# Verify installation
conda list
```

**Example: SciFact Environment**:
```bash
# Create Python 3.8 environment
conda create -n scifact python=3.8 -y

# Activate
conda activate scifact

# Install dependencies
pip install -r /Users/byron/projects/data/immunos_data/research/scifact/requirements.txt

# Verify PyTorch version
python -c "import torch; print(torch.__version__)"
# Expected: 1.5.0
```

**Example: Proteomics Environment**:
```bash
# Create Python 3.10 environment
conda create -n proteomics python=3.10 -y

# Activate
conda activate proteomics

# Install scientific stack via conda (recommended)
conda install -c conda-forge numpy scipy pandas scikit-learn matplotlib seaborn -y

# Install specialized packages via pip
pip install -r /Users/byron/projects/research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/config/requirements.txt

# Verify pertpy
python -c "import pertpy; print(pertpy.__version__)"
```

**Example: IMMUNOS Web Dashboard**:
```bash
# Create Python 3.11 environment
conda create -n immunos-web python=3.11 -y

# Activate
conda activate immunos-web

# Install dependencies
pip install -r /Users/byron/projects/immunos-mcp/web_app/requirements.txt

# Verify Flask
python -c "import flask; print(flask.__version__)"
```

---

### 3. Activating and Deactivating Environments

**Activate**:
```bash
conda activate ENV_NAME
```

**Deactivate**:
```bash
conda deactivate
```

**Check Current Environment**:
```bash
conda info --envs
# Active environment marked with *
```

**Switch Between Environments**:
```bash
conda activate scifact
# Work on SciFact replication...

conda activate proteomics
# Switch to proteomics analysis...

conda deactivate
# Return to base environment
```

---

### 4. Updating Environments

**Update All Packages in Environment**:
```bash
conda activate ENV_NAME
conda update --all -y
```

**Update Specific Package**:
```bash
conda activate ENV_NAME
conda update PACKAGE_NAME -y
# or
pip install --upgrade PACKAGE_NAME
```

**Sync with Updated requirements.txt**:
```bash
conda activate ENV_NAME
pip install -r /path/to/requirements.txt --upgrade
```

**Export Current Environment**:
```bash
conda activate ENV_NAME

# Conda-style export (includes all dependencies)
conda env export > environment_backup.yml

# Pip-style export (only pip packages)
pip freeze > requirements_backup.txt
```

---

## Project-to-Environment Mapping

| Project/Directory | Environment | Activation Command |
|-------------------|-------------|-------------------|
| `/Users/byron/projects/scripts/immunos_*.py` | `immunos-base` or `base` | `conda activate immunos-base` |
| `/Users/byron/projects/data/immunos_data/research/scifact/` | `scifact` | `conda activate scifact` |
| `/Users/byron/projects/research/experiments/proteomics-claims-2025/` | `proteomics` | `conda activate proteomics` |
| `/Users/byron/projects/immunos-mcp/web_app/` | `immunos-web` | `conda activate immunos-web` |
| `/Users/byron/projects/data/immunos_data/hallucination/truthfulqa/` | `truthfulqa` | `conda activate truthfulqa` |
| `/Users/byron/projects/rockbeatspaper/` | `rockbeatspaper` | `conda activate rockbeatspaper` |
| Machine learning experiments | `immunos-ml` | `conda activate immunos-ml` |

**Quick Activation Pattern**:
```bash
# Navigate to project
cd /Users/byron/projects/data/immunos_data/research/scifact

# Activate corresponding environment
conda activate scifact

# Run project scripts
python claim_verification/whatever.py
```

---

## Commands Reference

### Environment Management

**List All Environments**:
```bash
conda env list
# or
conda info --envs
```

**Create New Environment**:
```bash
# With specific Python version
conda create -n ENV_NAME python=3.X -y

# From environment.yml file
conda env create -f environment.yml

# Clone existing environment
conda create -n NEW_ENV --clone EXISTING_ENV
```

**Remove Environment**:
```bash
conda env remove -n ENV_NAME
```

**Activate/Deactivate**:
```bash
conda activate ENV_NAME
conda deactivate
```

---

### Package Management

**Install Packages via Conda**:
```bash
# Single package
conda install PACKAGE_NAME -y

# Multiple packages
conda install numpy scipy pandas -y

# From specific channel
conda install -c conda-forge PACKAGE_NAME -y

# Specific version
conda install numpy=1.24.0 -y
```

**Install Packages via Pip** (within conda env):
```bash
# Activate environment first
conda activate ENV_NAME

# Install with pip
pip install PACKAGE_NAME

# From requirements.txt
pip install -r requirements.txt

# Specific version
pip install numpy==1.24.0
```

**List Installed Packages**:
```bash
conda list
# or
pip list
```

**Search for Package**:
```bash
conda search PACKAGE_NAME
```

**Update Packages**:
```bash
# Update all in environment
conda update --all -y

# Update specific package
conda update PACKAGE_NAME -y

# Update via pip
pip install --upgrade PACKAGE_NAME
```

---

### Export and Backup

**Export Conda Environment**:
```bash
# Full environment (conda + pip)
conda env export > environment.yml

# Cross-platform (no builds)
conda env export --no-builds > environment.yml

# Only explicitly installed packages
conda env export --from-history > environment.yml
```

**Export Pip Requirements**:
```bash
pip freeze > requirements.txt

# Only packages installed via pip (not conda)
pip list --format=freeze > requirements_pip_only.txt
```

**Restore from Backup**:
```bash
# From conda environment.yml
conda env create -f environment.yml

# From pip requirements.txt
conda create -n ENV_NAME python=3.X -y
conda activate ENV_NAME
pip install -r requirements.txt
```

---

### Environment Information

**Show Environment Details**:
```bash
conda info

# Show environment location
conda info --envs

# Show package channels
conda config --show channels
```

**Check Python Version in Environment**:
```bash
conda activate ENV_NAME
python --version
```

---

## Troubleshooting Common Issues

### 1. Conda Not Found in PATH

**Problem**: `conda: command not found` after installation

**Solution**:
```bash
# Initialize conda for your shell
~/anaconda3/bin/conda init zsh  # or bash

# Restart terminal or source config
source ~/.zshrc

# Verify
conda --version
```

---

### 2. Dependency Conflicts During Install

**Problem**: `pip` fails with dependency version conflicts

**Solution A**: Install scientific packages via conda first
```bash
conda activate ENV_NAME

# Install numerical packages via conda (better dependency resolution)
conda install -c conda-forge numpy scipy pandas scikit-learn -y

# Then install specialized packages via pip
pip install pertpy scanpy  # etc.
```

**Solution B**: Use `--no-deps` for specific packages
```bash
pip install PACKAGE_NAME --no-deps
# Then manually install its dependencies
```

**Solution C**: Create fresh environment with exact versions
```bash
conda create -n ENV_NAME python=3.X numpy=1.Y.Z scipy=1.A.B -y
```

---

### 3. Old PyTorch/CUDA Version Conflicts

**Problem**: SciFact environment needs PyTorch 1.5.0, but CUDA is incompatible

**Solution**:
```bash
# Install CPU-only version for old PyTorch
conda activate scifact
pip install torch==1.5.0+cpu -f https://download.pytorch.org/whl/torch_stable.html
```

**For Modern PyTorch with GPU**:
```bash
conda activate immunos-ml
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y
```

---

### 4. GitHub Install Failures (TruthfulQA)

**Problem**: Custom GitHub packages fail to install from requirements.txt

**Solution**: Install them separately
```bash
conda activate truthfulqa

# Install BLEURT from GitHub
pip install https://github.com/google-research/bleurt/archive/b610120347ef22b494b6d69b4316e303f5932516.zip

# Install custom transformers fork
pip install https://github.com/finetuneanon/transformers/archive/1dd21b074751593bc65a992f90953c974e329511.zip

# Then install remaining packages
pip install -r requirements.txt
```

---

### 5. ScispaCy Model Download Fails

**Problem**: ScispaCy model URL unreachable or SSL error

**Solution**:
```bash
conda activate scifact

# Download model separately
wget https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.5/en_core_sci_sm-0.2.5.tar.gz

# Install locally
pip install en_core_sci_sm-0.2.5.tar.gz

# Or use spacy download
pip install spacy==2.2.4
python -m spacy download en_core_sci_sm
```

---

### 6. Jupyter Kernel Not Available

**Problem**: Jupyter can't find conda environment as kernel

**Solution**:
```bash
conda activate ENV_NAME

# Install ipykernel
conda install ipykernel -y

# Register environment as Jupyter kernel
python -m ipykernel install --user --name=ENV_NAME --display-name "Python (ENV_NAME)"

# Verify
jupyter kernelspec list
```

---

### 7. Environment Takes Too Much Disk Space

**Problem**: Multiple environments consuming GBs of space

**Solution A**: Clean conda cache
```bash
conda clean --all -y
```

**Solution B**: Remove unused environments
```bash
conda env list
conda env remove -n UNUSED_ENV
```

**Solution C**: Use hardlinks for packages (conda default)
```bash
# Check if already enabled
conda info

# Enable if not (usually default)
conda config --set hardlinks true
```

---

### 8. Permission Errors During Install

**Problem**: `PermissionError` when installing packages

**Solution**:
```bash
# Install to user directory (avoid system-wide)
pip install --user PACKAGE_NAME

# Or ensure conda environment is in user directory
conda create -n ENV_NAME python=3.X --prefix ~/conda_envs/ENV_NAME
```

---

### 9. Requirements.txt Has Incompatible Versions

**Problem**: `requirements.txt` specifies versions incompatible with Python or OS

**Solution**: Create modified requirements
```bash
# Copy requirements
cp requirements.txt requirements_modified.txt

# Edit to loosen version constraints
# Change: numpy==1.18.2
# To: numpy>=1.18.0

# Install modified version
pip install -r requirements_modified.txt
```

---

### 10. Multiple Python Installations Conflict

**Problem**: System Python, Homebrew Python, and Anaconda Python interfere

**Solution**: Always use conda environments
```bash
# Check which Python is active
which python
# Should show: ~/anaconda3/envs/ENV_NAME/bin/python

# If using base conda instead
conda activate base
which python
# Should show: ~/anaconda3/bin/python

# Never use system Python for projects
# Always: conda activate ENV_NAME first
```

---

## Best Practices

### 1. Always Use Environments
- Never install project dependencies in base conda
- Create dedicated environment per project
- Use `conda activate` before running project scripts

### 2. Document Dependencies
- Keep `requirements.txt` updated
- Export environment after major changes
- Document Python version in project README

### 3. Install Order Matters
- Install scientific stack via **conda** (NumPy, SciPy, Pandas, scikit-learn)
- Install specialized packages via **pip** (pertpy, custom forks)
- Conda handles binary dependencies better

### 4. Version Pinning
- Pin exact versions in `requirements.txt` for reproducibility
- Use `>=` for flexibility when exact version doesn't matter
- Document why specific versions are required (e.g., "PyTorch 1.5.0 for SciFact compatibility")

### 5. Regular Maintenance
- Run `conda clean --all` periodically
- Remove unused environments
- Update packages in active environments quarterly

### 6. Backup Critical Environments
```bash
conda activate CRITICAL_ENV
conda env export > CRITICAL_ENV_backup_$(date +%Y%m%d).yml
```

---

## Quick Start Checklist

**Setting Up New Project Environment**:
- [ ] Create conda environment: `conda create -n PROJECT python=3.X -y`
- [ ] Activate: `conda activate PROJECT`
- [ ] Install scientific stack: `conda install -c conda-forge numpy scipy pandas -y`
- [ ] Install project dependencies: `pip install -r requirements.txt`
- [ ] Register Jupyter kernel: `python -m ipykernel install --user --name=PROJECT`
- [ ] Test imports: `python -c "import MODULE; print(MODULE.__version__)"`
- [ ] Export for backup: `conda env export > environment_$(date +%Y%m%d).yml`

**Daily Workflow**:
- [ ] Navigate to project directory
- [ ] Activate environment: `conda activate PROJECT`
- [ ] Run scripts
- [ ] Deactivate when done: `conda deactivate`

---

## Related Documentation

- **IMMUNOS System**: `/Users/byron/projects/CLAUDE.md`
- **Proteomics Analysis Guide**: `/Users/byron/projects/research/experiments/proteomics-claims-2025/CLAUDE.md`
- **SciFact Replication Notes**: `/Users/byron/projects/docs/reference/scifact-baseline-replication.md`
- **Conda Documentation**: https://docs.conda.io/
- **Pip Documentation**: https://pip.pypa.io/

---

**Maintained by**: Byron (with Claude Sonnet 4.5)
**Environment Count**: 7+ active environments
**Total Disk Usage**: Run `du -sh ~/anaconda3/envs/` to check

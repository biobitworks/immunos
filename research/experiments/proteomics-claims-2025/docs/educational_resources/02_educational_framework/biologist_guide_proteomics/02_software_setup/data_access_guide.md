# 📁 Complete Data Access Guide: Local vs Cloud Analysis

## 🎯 Why This Guide Matters

**The Problem**: Many biologists get confused about where their data is and how to access it when switching between local computers and cloud platforms like Google Colab.

**The Solution**: This guide provides clear, step-by-step instructions for accessing your proteomics data in any environment, with specific examples for our analyses.

---

## 🗂️ Understanding Data Locations

### Local Analysis (Jupyter on Your Computer)
```
Your Computer
├── Documents/
│   └── proteomics_analysis/
│       ├── data/
│       │   ├── pool_processed_v2.h5ad    ← Your data is HERE
│       │   └── metadata.csv
│       └── notebooks/
│           └── analysis.ipynb            ← Your code runs HERE
```
**Data Access**: Direct file system access, very fast

### Cloud Analysis (Google Colab)
```
Google's Servers
├── /content/                             ← Temporary storage (deleted!)
│   └── your_notebook.ipynb              ← Your code runs HERE
└── /content/drive/MyDrive/               ← Your Google Drive
    └── proteomics_data/
        ├── pool_processed_v2.h5ad        ← Your data is HERE
        └── metadata.csv
```
**Data Access**: Through Google Drive, requires mounting

---

## 📥 Getting Your Data Ready

### Step 1: Download the Proteomics Dataset

#### Option A: Direct Download (If Available Online)
```python
# Example download script (adapt URL as needed)
import requests
import os

def download_file(url, local_filename):
    """Download large file with progress"""
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        downloaded = 0

        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)

                # Progress indicator
                if total_size > 0:
                    percent = (downloaded / total_size) * 100
                    print(f"\rDownloading: {percent:.1f}%", end="")

    print(f"\nDownload complete: {local_filename}")

# Download the dataset
url = "https://your-data-source.com/pool_processed_v2.h5ad"
download_file(url, "pool_processed_v2.h5ad")
```

#### Option B: Contact Data Provider
If the dataset isn't publicly available:
1. **Contact the research group** that generated the data
2. **Request access** to the proteomics dataset
3. **Download to your computer** for local analysis
4. **Upload to Google Drive** for cloud analysis

### Step 2: Organize Your Data

#### Create Consistent Structure (Important!)
Whether using local or cloud, always use the same folder structure:

```
proteomics_analysis/
├── data/
│   ├── raw/
│   │   ├── pool_processed_v2.h5ad        # Main dataset (large file)
│   │   └── metadata.csv                  # Sample information
│   ├── processed/
│   │   └── [analysis outputs]
│   └── annotations/
│       └── protein_annotations.csv       # Protein functional info
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_ups_analysis.ipynb
│   └── 03_sqstm1_analysis.ipynb
├── results/
│   ├── tables/
│   └── statistics/
└── plots/
    ├── exploratory/
    └── publication/
```

---

## 💻 Local Data Access (Jupyter on Your Computer)

### Setup: Create Project Directory

#### Windows
```bash
# Open Command Prompt or PowerShell
mkdir C:\Users\YourName\Documents\proteomics_analysis
cd C:\Users\YourName\Documents\proteomics_analysis
mkdir data notebooks results plots
mkdir data\raw data\processed data\annotations
```

#### Mac/Linux
```bash
# Open Terminal
mkdir ~/Documents/proteomics_analysis
cd ~/Documents/proteomics_analysis
mkdir -p data/{raw,processed,annotations} notebooks results plots
```

### File Placement
1. **Place downloaded data** in `data/raw/` folder
2. **Start Jupyter** from the project root directory
3. **Create notebooks** in `notebooks/` folder

### Code Examples for Local Access

#### Loading Data (Relative Paths)
```python
# From a notebook in notebooks/ folder
import scanpy as sc
import pandas as pd
import numpy as np

# Load main dataset (adjust path from notebooks/ to data/)
print("Loading proteomics dataset...")
adata = sc.read_h5ad('../data/raw/pool_processed_v2.h5ad')
print(f"✅ Dataset loaded: {adata.shape}")

# Load metadata
metadata = pd.read_csv('../data/raw/metadata.csv')
print(f"✅ Metadata loaded: {metadata.shape}")

# Load protein annotations (if available)
try:
    protein_info = pd.read_csv('../data/annotations/protein_annotations.csv')
    print(f"✅ Protein annotations loaded: {protein_info.shape}")
except FileNotFoundError:
    print("⚠️ Protein annotations not found - continuing without")
```

#### Checking File Locations
```python
import os

# Check current directory
print("Current working directory:", os.getcwd())

# List files in data directory
print("\nFiles in data/raw/:")
try:
    files = os.listdir('../data/raw/')
    for file in files:
        size = os.path.getsize(f'../data/raw/{file}')
        print(f"  {file} ({size/1024/1024:.1f} MB)")
except FileNotFoundError:
    print("  Directory not found - check your file structure")

# Verify data file exists
data_path = '../data/raw/pool_processed_v2.h5ad'
if os.path.exists(data_path):
    print(f"✅ Main dataset found at: {data_path}")
else:
    print(f"❌ Main dataset NOT found at: {data_path}")
    print("Check file location and name")
```

#### Saving Results Locally
```python
# Save analysis results
import pandas as pd

# Example: Save UPS analysis results
ups_results = pd.DataFrame({
    'protein': ['PSMA1', 'PSMA2', 'PSMB1'],
    'log2_fold_change': [0.2, -0.1, 0.5],
    'p_value': [0.03, 0.4, 0.001]
})

# Save to results directory
ups_results.to_csv('../results/tables/ups_analysis_results.csv', index=False)
print("✅ Results saved to ../results/tables/ups_analysis_results.csv")

# Save plots
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.scatter(ups_results['log2_fold_change'], -np.log10(ups_results['p_value']))
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.title('UPS Protein Volcano Plot')

# Save plot
plt.savefig('../plots/publication/ups_volcano_plot.png', dpi=300, bbox_inches='tight')
plt.show()
print("✅ Plot saved to ../plots/publication/ups_volcano_plot.png")
```

---

## ☁️ Google Colab Data Access

### Setup: Upload Data to Google Drive

#### Method 1: Through Google Drive Web Interface
1. **Go to**: https://drive.google.com/
2. **Create folder**: "proteomics_analysis"
3. **Create subfolders**: data, notebooks, results, plots
4. **Upload your data files** to data/ folder
   - Drag and drop `pool_processed_v2.h5ad`
   - Upload `metadata.csv` and other files

#### Method 2: Using Google Drive Desktop App
1. **Install Google Drive** on your computer
2. **Sync proteomics_analysis folder** to Drive
3. **Files automatically available** in Colab

### Code Examples for Colab Access

#### Essential Setup (Run This First!)
```python
# MUST RUN FIRST: Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

# Verify mount successful
import os
print("Google Drive contents:")
os.listdir('/content/drive/MyDrive')
```

#### Navigate to Your Data
```python
# Set up paths (adjust to your folder structure)
PROJECT_PATH = '/content/drive/MyDrive/proteomics_analysis'
DATA_PATH = f'{PROJECT_PATH}/data/raw'
RESULTS_PATH = f'{PROJECT_PATH}/results'

# Check if project folder exists
if os.path.exists(PROJECT_PATH):
    print(f"✅ Project folder found: {PROJECT_PATH}")
    print("Data files:")
    try:
        for file in os.listdir(DATA_PATH):
            size = os.path.getsize(f'{DATA_PATH}/{file}')
            print(f"  {file} ({size/1024/1024:.1f} MB)")
    except FileNotFoundError:
        print("  Data folder not found - check folder structure")
else:
    print(f"❌ Project folder not found: {PROJECT_PATH}")
    print("Available folders in Drive:")
    print(os.listdir('/content/drive/MyDrive'))
```

#### Loading Data from Google Drive
```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load main dataset
print("Loading proteomics dataset from Google Drive...")
try:
    adata = sc.read_h5ad(f'{DATA_PATH}/pool_processed_v2.h5ad')
    print(f"✅ Dataset loaded: {adata.shape}")
    print(f"Memory usage: {adata.X.nbytes / 1024 / 1024:.1f} MB")
except FileNotFoundError:
    print("❌ Dataset file not found")
    print("Check file name and location in Google Drive")
except Exception as e:
    print(f"❌ Error loading dataset: {e}")

# Load metadata
try:
    metadata = pd.read_csv(f'{DATA_PATH}/metadata.csv')
    print(f"✅ Metadata loaded: {metadata.shape}")
except FileNotFoundError:
    print("⚠️ Metadata file not found - check if uploaded")

# Quick data exploration
if 'adata' in locals():
    print(f"\n📊 Dataset overview:")
    print(f"  Observations (samples): {adata.n_obs}")
    print(f"  Variables (proteins): {adata.n_vars}")
    print(f"  Available metadata: {list(adata.obs.columns)}")
```

#### Handling Large Files in Colab
```python
# Check available memory
!cat /proc/meminfo | grep MemAvailable

# For very large datasets, work in chunks
def load_data_efficiently():
    """Load large datasets efficiently in Colab"""
    try:
        # Try loading full dataset
        adata = sc.read_h5ad(f'{DATA_PATH}/pool_processed_v2.h5ad')
        print(f"✅ Full dataset loaded: {adata.shape}")
        return adata

    except MemoryError:
        print("⚠️ Memory error - trying efficient loading")
        # Alternative approaches for large data

        # Option 1: Load only specific proteins
        adata = sc.read_h5ad(f'{DATA_PATH}/pool_processed_v2.h5ad')
        # Select subset of proteins for initial exploration
        ups_proteins = ['PSMA1', 'PSMA2', 'PSMB1', 'SQSTM1', 'VDAC1']
        available_ups = [p for p in ups_proteins if p in adata.var_names]
        adata_subset = adata[:, available_ups]
        print(f"✅ Subset loaded: {adata_subset.shape}")
        return adata_subset

# Load data efficiently
adata = load_data_efficiently()
```

#### Saving Results in Colab
```python
# Create results directories if they don't exist
os.makedirs(f'{RESULTS_PATH}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_PATH}/plots', exist_ok=True)

# Save analysis results to Drive
results_df = pd.DataFrame({
    'protein': ['SQSTM1', 'VDAC1', 'PSMA1'],
    'fold_change': [10.7, 1.2, 0.8],
    'p_value': [0.001, 0.3, 0.6]
})

# Save to Google Drive (persistent)
results_df.to_csv(f'{RESULTS_PATH}/tables/protein_analysis_results.csv', index=False)
print("✅ Results saved to Google Drive")

# Also save to temporary Colab storage (faster access)
results_df.to_csv('/content/temp_results.csv', index=False)
print("✅ Temporary copy saved for current session")
```

### Colab-Specific Tips

#### Managing Session Timeouts
```python
# Save critical variables before long operations
import pickle

# Save important data structures
with open(f'{RESULTS_PATH}/analysis_state.pkl', 'wb') as f:
    pickle.dump({
        'adata': adata,
        'results': results_df,
        'parameters': analysis_parameters
    }, f)

print("✅ Analysis state saved to Google Drive")

# Load saved state in new session
with open(f'{RESULTS_PATH}/analysis_state.pkl', 'rb') as f:
    saved_state = pickle.load(f)
    adata = saved_state['adata']
    results_df = saved_state['results']

print("✅ Previous analysis state restored")
```

#### Download Results from Colab
```python
from google.colab import files

# Download specific files to your computer
files.download('/content/drive/MyDrive/proteomics_analysis/results/tables/final_results.csv')

# Or zip multiple files
!zip -r analysis_results.zip /content/drive/MyDrive/proteomics_analysis/results/
files.download('analysis_results.zip')
```

---

## 🔄 Switching Between Environments

### Making Code Work in Both Environments

#### Universal Path Detection
```python
import os

def get_data_path():
    """Automatically detect if running locally or in Colab"""
    if 'COLAB_GPU' in os.environ:
        # Running in Google Colab
        print("🌐 Detected Google Colab environment")
        base_path = '/content/drive/MyDrive/proteomics_analysis'

        # Check if Drive is mounted
        if not os.path.exists('/content/drive'):
            print("❌ Google Drive not mounted!")
            print("Run: from google.colab import drive; drive.mount('/content/drive')")
            return None

    else:
        # Running locally
        print("💻 Detected local environment")
        base_path = '../..'  # Adjust based on notebook location

    data_path = f'{base_path}/data/raw'

    if os.path.exists(data_path):
        print(f"✅ Data path found: {data_path}")
        return data_path
    else:
        print(f"❌ Data path not found: {data_path}")
        return None

# Use in your analysis
DATA_PATH = get_data_path()
if DATA_PATH:
    adata = sc.read_h5ad(f'{DATA_PATH}/pool_processed_v2.h5ad')
```

#### Universal Package Installation
```python
import subprocess
import sys

def install_if_missing(packages):
    """Install packages if not available"""
    for package in packages:
        try:
            __import__(package)
            print(f"✅ {package} already available")
        except ImportError:
            print(f"📦 Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Required packages for our analysis
required_packages = ['scanpy', 'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels']
install_if_missing(required_packages)

print("✅ All packages ready!")
```

### Syncing Work Between Environments

#### Export from Local to Colab
```python
# On local computer: Export notebook to Drive
# 1. Save notebook in Google Drive folder
# 2. Open in Colab from Drive
# 3. Data already available in same Drive folder
```

#### Export from Colab to Local
```python
# Download notebook from Colab
# File → Download → .ipynb
# Place in your local notebooks/ folder
# Data paths may need adjustment
```

---

## 🚨 Common Problems and Solutions

### Problem 1: "File Not Found" Errors

#### Diagnosis
```python
import os

# Check current directory
print("Current directory:", os.getcwd())

# List all files in current directory
print("\nFiles in current directory:")
for item in os.listdir('.'):
    print(f"  {item}")

# Check if specific file exists
file_path = 'pool_processed_v2.h5ad'  # Adjust as needed
if os.path.exists(file_path):
    print(f"✅ Found: {file_path}")
else:
    print(f"❌ Not found: {file_path}")

    # Search for similar files
    all_files = os.listdir('.')
    h5ad_files = [f for f in all_files if f.endswith('.h5ad')]
    if h5ad_files:
        print("Available .h5ad files:")
        for f in h5ad_files:
            print(f"  {f}")
```

#### Solutions
```python
# Solution 1: Use correct relative paths
# From notebooks/ folder:
correct_path = '../data/raw/pool_processed_v2.h5ad'

# Solution 2: Use absolute paths (Colab)
colab_path = '/content/drive/MyDrive/proteomics_analysis/data/raw/pool_processed_v2.h5ad'

# Solution 3: Find file automatically
import glob
h5ad_files = glob.glob('**/*.h5ad', recursive=True)
if h5ad_files:
    print("Found .h5ad files:")
    for f in h5ad_files:
        print(f"  {f}")
    data_file = h5ad_files[0]  # Use first found
else:
    print("No .h5ad files found in directory tree")
```

### Problem 2: Memory Errors with Large Files

#### For Local Jupyter
```python
# Check available memory
import psutil
memory_info = psutil.virtual_memory()
print(f"Available memory: {memory_info.available / 1024**3:.1f} GB")

# Load data efficiently
def load_large_dataset(filepath, max_memory_gb=4):
    """Load dataset with memory management"""
    file_size = os.path.getsize(filepath) / 1024**3
    print(f"File size: {file_size:.1f} GB")

    if file_size > max_memory_gb:
        print("File too large, loading subset...")
        # Implement chunked loading or subset selection
        adata = sc.read_h5ad(filepath)
        # Select subset of samples or features
        return adata[:1000, :1000]  # Example subset
    else:
        return sc.read_h5ad(filepath)
```

#### For Google Colab
```python
# Upgrade to Colab Pro for more memory
# Or work with data subsets

# Clear memory between analyses
import gc

def clear_memory():
    """Clear Python memory"""
    gc.collect()
    print("Memory cleared")

# Use after large operations
clear_memory()
```

### Problem 3: Package Version Conflicts

#### Check Versions
```python
import sys
import numpy as np
import pandas as pd
import scanpy as sc

print(f"Python version: {sys.version}")
print(f"NumPy version: {np.__version__}")
print(f"Pandas version: {pd.__version__}")
print(f"Scanpy version: {sc.__version__}")
```

#### Fix Version Issues
```bash
# Local Jupyter: Update packages
pip install --upgrade scanpy pandas numpy matplotlib

# Colab: Usually packages are up to date
# If needed:
!pip install --upgrade scanpy
```

---

## ✅ Final Checklist

### Before Starting Analysis

#### Local Jupyter Setup
- [ ] **Project directory created** with proper structure
- [ ] **Data files downloaded** and placed in data/raw/
- [ ] **Jupyter started** from project root directory
- [ ] **Test data loading** with simple script
- [ ] **Required packages installed** and working

#### Google Colab Setup
- [ ] **Google Drive folder created** (proteomics_analysis)
- [ ] **Data files uploaded** to Drive/proteomics_analysis/data/raw/
- [ ] **Drive mounted** in Colab session
- [ ] **Test data loading** from Drive paths
- [ ] **Packages installed** (usually automatic in Colab)

#### Universal Checks
- [ ] **Data file verified**: Can load pool_processed_v2.h5ad
- [ ] **Metadata available**: Sample information accessible
- [ ] **Paths working**: Code finds files without errors
- [ ] **Memory sufficient**: Dataset fits in available RAM
- [ ] **Backup strategy**: Know how to save work

---

**You're now equipped to handle data access in any environment!** Whether working locally for maximum control or in the cloud for easy collaboration, you can confidently load and analyze your proteomics data.

*Next: Understanding Your Proteomics Dataset - now that you can access the data, let's explore what's in it!*

*Remember: Organized data management prevents 90% of analysis headaches!* 📁🔬
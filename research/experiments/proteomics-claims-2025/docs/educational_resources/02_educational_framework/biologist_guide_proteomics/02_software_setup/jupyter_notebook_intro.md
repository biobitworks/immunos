# 📓 Complete Jupyter Notebook Guide for Biologists

## 🎯 What You'll Learn

By the end of this guide, you'll be able to:
- ✅ **Choose between local Jupyter vs Google Colab** based on your needs
- ✅ **Start and navigate Jupyter notebooks** in both environments
- ✅ **Upload and access data files** locally and in the cloud
- ✅ **Write and run Python code** for proteomics analysis
- ✅ **Create publication-quality reports** mixing code, text, and plots
- ✅ **Share your work** with collaborators effectively

---

## 📚 What is a Jupyter Notebook?

### Think of it as a Lab Notebook for Code

#### Traditional Lab Notebook
```
Date: March 15, 2024
Experiment: Western blot for SQSTM1 protein
Protocol: [detailed steps]
Results: [photos, measurements]
Analysis: [calculations, interpretations]
Conclusions: [what this means]
```

#### Jupyter Notebook (Digital Lab Notebook)
```
# March 15, 2024 - SQSTM1 Proteomics Analysis

## Background
SQSTM1 is an autophagy receptor protein...

## Code to Load Data
[Python code that runs immediately]

## Results
[Automatically generated plots and statistics]

## Analysis and Conclusions
[Your interpretations and next steps]
```

### Key Advantages for Biologists
- **Mix explanations with code**: Document your thinking process
- **Immediate feedback**: See results as you go
- **Reproducible**: Others can run your exact analysis
- **Shareable**: Send complete analysis to collaborators
- **Publication-ready**: Generate figures and reports

---

## 🔄 Local Jupyter vs Google Colab: Complete Comparison

### Local Jupyter Notebooks

#### ✅ Advantages
- **Your data stays private**: No uploading to cloud servers
- **Unlimited storage**: Use your computer's full hard drive
- **No internet required**: Work offline once set up
- **Full control**: Install any packages you want
- **Faster file access**: Data on your local machine
- **No time limits**: Run analyses for days if needed

#### ❌ Disadvantages
- **Installation required**: Need to set up Python environment
- **Limited by your hardware**: Your computer's RAM and CPU
- **No collaboration**: Harder to share with colleagues
- **Backup responsibility**: You must save and backup work
- **Software maintenance**: Update packages yourself

#### 💡 Best For
- **Sensitive data**: Patient data, proprietary datasets
- **Large datasets**: Multi-gigabyte proteomics files
- **Long analyses**: Complex statistical modeling
- **Regular use**: Daily bioinformatics work

### Google Colab

#### ✅ Advantages
- **No installation**: Works in any web browser
- **Free powerful computers**: Access to GPUs and high RAM
- **Easy sharing**: Share links like Google Docs
- **Automatic saving**: Work saved to Google Drive
- **Collaboration**: Multiple people can edit together
- **No maintenance**: Google handles all software updates

#### ❌ Disadvantages
- **Data privacy concerns**: Files uploaded to Google servers
- **Limited storage**: 15GB free Google Drive storage
- **Internet required**: Need stable connection
- **Session timeouts**: Long analyses may be interrupted
- **Package limitations**: Some specialized tools not available

#### 💡 Best For
- **Learning**: Perfect for tutorials and education
- **Small datasets**: <1GB proteomics files
- **Collaboration**: Working with remote colleagues
- **Occasional use**: Exploring computational biology
- **Teaching**: Classroom and workshop settings

---

## 🚀 Getting Started with Local Jupyter

### Prerequisites
- Python installed (Anaconda recommended - see installation guide)
- Basic familiarity with file navigation

### Step 1: Start Jupyter Notebook

#### Method A: From Anaconda Navigator (Easiest)
1. **Open Anaconda Navigator**
   - Windows: Start Menu → Anaconda3 → Anaconda Navigator
   - Mac: Applications → Anaconda Navigator
   - Linux: Type `anaconda-navigator` in terminal

2. **Launch Jupyter**
   - Click "Launch" button under Jupyter Notebook
   - Your web browser will open automatically
   - You'll see the Jupyter file browser

#### Method B: From Command Line/Terminal
```bash
# Open terminal (Mac/Linux) or Anaconda Prompt (Windows)
jupyter notebook

# Browser opens automatically to: http://localhost:8888
```

### Step 2: Navigate to Your Project Folder

#### Set Up Project Directory First
```bash
# Create organized project structure
mkdir ~/proteomics_analysis
cd ~/proteomics_analysis
mkdir data notebooks results plots

# Start Jupyter from project directory
jupyter notebook
```

#### In Jupyter File Browser
1. **Navigate to your project folder**
   - Click folders to navigate
   - Use "Upload" button to add data files
   - Create new folders with "New" → "Folder"

2. **Create your first notebook**
   - Click "New" → "Python 3"
   - Rename: Click "Untitled" at top → "SQSTM1_Analysis"

### Step 3: Upload and Access Local Data Files

#### Organize Your Data
```
proteomics_analysis/
├── data/
│   ├── pool_processed_v2.h5ad     # Main dataset
│   ├── metadata.csv               # Sample information
│   └── protein_annotations.csv    # Protein details
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_ups_analysis.ipynb
│   └── 03_sqstm1_analysis.ipynb
├── results/
└── plots/
```

#### Upload Data via Jupyter Interface
1. **Click "Upload" button** in Jupyter file browser
2. **Select your data files**
3. **Click blue "Upload" button** to confirm
4. **Files appear in current directory**

#### Access Data in Your Code
```python
# Load data from relative paths
import scanpy as sc
import pandas as pd

# Load main proteomics dataset
adata = sc.read_h5ad('../data/pool_processed_v2.h5ad')

# Load metadata
metadata = pd.read_csv('../data/metadata.csv')

print(f"Dataset loaded: {adata.shape}")
print(f"Metadata loaded: {metadata.shape}")
```

---

## 🌐 Getting Started with Google Colab

### Step 1: Access Google Colab

#### Open Colab
1. **Go to**: https://colab.research.google.com/
2. **Sign in** with your Google account
3. **You'll see the Colab welcome page**

#### Alternative: Direct from Google Drive
1. **Go to**: https://drive.google.com/
2. **Click "New"** → **"More"** → **"Google Colaboratory"**
3. **If not visible**: Click "Connect more apps" and search "Colab"

### Step 2: Create Your First Notebook

#### Start New Notebook
1. **Click "New Notebook"** or **File → New Notebook**
2. **Rename notebook**: Click "Untitled0.ipynb" at top
3. **Type**: "SQSTM1_Proteomics_Analysis"

#### Connect to Google Drive (Essential!)
```python
# ALWAYS run this first in Colab
from google.colab import drive
drive.mount('/content/drive')

# You'll see authorization prompt - click link and authorize
print("Google Drive mounted successfully!")
```

### Step 3: Upload and Access Data in Colab

#### Method 1: Upload Files Directly (Small Files <25MB)
```python
from google.colab import files

# Upload file from your computer
uploaded = files.upload()

# File appears in /content/ directory
import pandas as pd
for filename in uploaded.keys():
    print(f"Uploaded: {filename}")
```

#### Method 2: Via Google Drive (Recommended for Large Files)

##### Upload to Drive First
1. **Go to**: https://drive.google.com/
2. **Create folder**: "proteomics_data"
3. **Upload your files**: pool_processed_v2.h5ad, metadata.csv, etc.

##### Access from Colab
```python
# Mount Google Drive (run once per session)
from google.colab import drive
drive.mount('/content/drive')

# Navigate to your data
import os
os.listdir('/content/drive/MyDrive/proteomics_data')

# Load your data
import scanpy as sc
adata = sc.read_h5ad('/content/drive/MyDrive/proteomics_data/pool_processed_v2.h5ad')
print(f"Dataset loaded: {adata.shape}")
```

#### Method 3: Download from Web (If Data is Online)
```python
# Download directly from URL
!wget https://example.com/pool_processed_v2.h5ad

# Or using Python
import requests
url = "https://example.com/your_data.csv"
response = requests.get(url)
with open('data.csv', 'wb') as f:
    f.write(response.content)
```

### Step 4: Manage Colab Limitations

#### Storage Management
```python
# Check available space
!df -h

# Clean up temporary files
!rm -rf /content/sample_data  # Remove example files
```

#### Handle Session Timeouts
```python
# Save results frequently
import pickle

# Save important variables
with open('/content/drive/MyDrive/results.pkl', 'wb') as f:
    pickle.dump(results_dataframe, f)

# Load results in new session
with open('/content/drive/MyDrive/results.pkl', 'rb') as f:
    results_dataframe = pickle.load(f)
```

#### Keep Session Alive (For Long Analyses)
```javascript
// Run this in browser console to prevent disconnection
function ClickConnect(){
    console.log("Working");
    document.querySelector("colab-connect-button").click()
}
setInterval(ClickConnect,60000)
```

---

## 📝 Jupyter Notebook Basics

### Understanding Cells

#### Cell Types
1. **Code cells**: Contain Python code that executes
2. **Markdown cells**: Contain formatted text, explanations
3. **Raw cells**: Plain text (rarely used)

#### Working with Cells
```python
# This is a code cell
print("Hello, proteomics!")

# You can have multiple lines
import pandas as pd
import numpy as np
```

```markdown
# This is a markdown cell
## You can format text with headers

- Lists
- **Bold text**
- *Italic text*

You can include equations: $\alpha = 0.05$
```

#### Essential Keyboard Shortcuts
- **Shift + Enter**: Run current cell and move to next
- **Ctrl + Enter**: Run current cell and stay
- **A**: Insert cell above
- **B**: Insert cell below
- **DD**: Delete cell
- **M**: Convert to markdown
- **Y**: Convert to code
- **Esc**: Command mode
- **Enter**: Edit mode

### Running Code

#### Execute Cells
```python
# Cell 1: Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

print("Libraries imported!")
```

```python
# Cell 2: Load data (depends on Cell 1)
data = pd.read_csv('your_data.csv')
print(f"Data shape: {data.shape}")
```

```python
# Cell 3: Analyze data (depends on Cell 2)
print(data.describe())
```

#### Important: Cell Execution Order
- **Cells execute in order you run them**, not written order
- **Variables persist** between cells
- **Always run cells from top to bottom** when starting
- **Use "Restart & Run All"** to ensure reproducibility

### Creating Publication-Quality Reports

#### Structure Your Notebook
```markdown
# Proteomics Analysis Report
**Author**: Your Name
**Date**: March 2024
**Dataset**: pool_processed_v2.h5ad

## Abstract
Brief summary of analysis and key findings...

## Introduction
Biological background and objectives...

## Methods
```

```python
# Detailed analysis code with comments
import scanpy as sc

# Load and explore data
adata = sc.read_h5ad('data/pool_processed_v2.h5ad')
print(f"Dataset dimensions: {adata.shape}")

# Quality control
sc.pl.highest_expr_genes(adata, n_top=20)
```

```markdown
## Results
Key findings from the analysis:

1. **SQSTM1 Expression**: Shows 10.7-fold upregulation
2. **Statistical Significance**: p < 0.001 after FDR correction
3. **Biological Interpretation**: Indicates autophagy dysfunction

## Discussion
These results suggest...

## Conclusions
```

#### Professional Formatting Tips
```python
# Clean, well-commented code
def analyze_protein_expression(adata, protein_name, group_col):
    """
    Analyze differential expression of single protein

    Parameters:
    -----------
    adata : AnnData object
        Proteomics dataset
    protein_name : str
        Name of protein to analyze
    group_col : str
        Column name for grouping (e.g., 'tau_status')

    Returns:
    --------
    dict : Analysis results
    """
    # Implementation here
    pass

# Use descriptive variable names
tau_positive_samples = adata[adata.obs['tau_status'] == 'positive']
tau_negative_samples = adata[adata.obs['tau_status'] == 'negative']

# Clear output formatting
print(f"{'Protein':<12} {'Fold Change':<12} {'P-value':<10}")
print("-" * 35)
for protein, fc, pval in results:
    print(f"{protein:<12} {fc:<12.2f} {pval:<10.3e}")
```

---

## 📊 Data Management Best Practices

### Local Jupyter File Organization

#### Recommended Structure
```
your_proteomics_project/
├── data/
│   ├── raw/                    # Original, unmodified data
│   │   └── pool_processed_v2.h5ad
│   ├── processed/              # Cleaned, filtered data
│   │   └── filtered_proteins.h5ad
│   └── annotations/            # Metadata, protein info
│       └── protein_functions.csv
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_quality_control.ipynb
│   ├── 03_ups_analysis.ipynb
│   ├── 04_sqstm1_analysis.ipynb
│   └── 05_sliding_window.ipynb
├── scripts/                    # Reusable Python functions
│   ├── analysis_functions.py
│   └── plotting_functions.py
├── results/
│   ├── tables/
│   ├── statistics/
│   └── processed_data/
├── plots/
│   ├── exploratory/
│   ├── publication/
│   └── supplementary/
└── reports/
    ├── interim_reports/
    └── final_report.ipynb
```

#### File Naming Conventions
```python
# Good file names (descriptive, sorted)
01_data_exploration_2024_03_15.ipynb
02_ups_protein_analysis_v2.ipynb
03_sqstm1_bootstrap_analysis.ipynb

# Bad file names (vague, unsorted)
analysis.ipynb
test.ipynb
untitled.ipynb
```

### Google Colab File Management

#### Google Drive Organization
```
My Drive/
├── proteomics_project/
│   ├── data/
│   │   ├── pool_processed_v2.h5ad
│   │   └── metadata.csv
│   ├── notebooks/
│   │   ├── 01_data_exploration.ipynb
│   │   └── 02_analysis.ipynb
│   ├── results/
│   │   ├── ups_analysis_results.csv
│   │   └── sqstm1_results.csv
│   └── plots/
│       ├── volcano_plot.png
│       └── correlation_heatmap.png
```

#### Saving and Loading Between Sessions
```python
# Save important results to Drive
import pickle
import pandas as pd

# Save DataFrames
results_df.to_csv('/content/drive/MyDrive/proteomics_project/results/ups_results.csv')

# Save complex objects
with open('/content/drive/MyDrive/proteomics_project/results/analysis_object.pkl', 'wb') as f:
    pickle.dump(analysis_results, f)

# Load in new session
results_df = pd.read_csv('/content/drive/MyDrive/proteomics_project/results/ups_results.csv')

with open('/content/drive/MyDrive/proteomics_project/results/analysis_object.pkl', 'rb') as f:
    analysis_results = pickle.load(f)
```

---

## 🤝 Sharing and Collaboration

### Sharing Local Jupyter Notebooks

#### Export Options
```python
# From Jupyter interface:
# File → Download as → HTML (.html)      # For viewing
# File → Download as → PDF via LaTeX     # For reports
# File → Download as → Python (.py)      # Code only
```

#### Version Control with Git (Advanced)
```bash
# Initialize git repository
git init
git add notebooks/
git commit -m "Initial analysis"

# Share on GitHub
git remote add origin https://github.com/username/proteomics-analysis
git push -u origin main
```

### Sharing Google Colab Notebooks

#### Share Like Google Docs
1. **Click "Share" button** (top right)
2. **Set permissions**:
   - **Viewer**: Can see but not edit
   - **Commenter**: Can add comments
   - **Editor**: Can modify code

3. **Copy link** and send to collaborators

#### Export from Colab
```python
# Download to your computer
# File → Download → .ipynb (notebook format)
# File → Download → .py (Python script)

# Save to GitHub directly
# File → Save a copy in GitHub
```

### Collaboration Best Practices

#### Working Together Effectively
1. **Clear cell documentation**: Explain what each section does
2. **Consistent naming**: Use same variable names across notebooks
3. **Version comments**: Add dates and changes at top
4. **Separate concerns**: One notebook per major analysis
5. **Test on clean environment**: "Restart & Run All" before sharing

---

## 🔧 Troubleshooting Common Issues

### Local Jupyter Problems

#### Jupyter Won't Start
```bash
# Check if Jupyter is installed
jupyter --version

# If not installed:
pip install jupyter

# If still problems, try:
conda install jupyter
```

#### Can't Find Data Files
```python
# Check current directory
import os
print("Current directory:", os.getcwd())
print("Files here:", os.listdir('.'))

# Use absolute paths if needed
data_path = "/full/path/to/your/data/file.h5ad"
```

#### Package Import Errors
```bash
# Install missing packages
pip install scanpy pandas numpy matplotlib seaborn

# Or with conda:
conda install -c conda-forge scanpy pandas numpy matplotlib seaborn
```

### Google Colab Problems

#### Drive Mount Fails
```python
# Try remounting
from google.colab import drive
drive.mount('/content/drive', force_remount=True)

# Check authorization
import os
os.listdir('/content/drive')
```

#### Session Disconnects
```python
# Save work frequently
results.to_csv('/content/drive/MyDrive/backup_results.csv')

# Use try-except for long operations
try:
    long_running_analysis()
except Exception as e:
    print(f"Error occurred: {e}")
    # Save partial results
```

#### Memory Errors
```python
# Check memory usage
!cat /proc/meminfo

# Use data in chunks
chunk_size = 1000
for i in range(0, len(data), chunk_size):
    chunk = data[i:i+chunk_size]
    # Process chunk
```

### General Debugging Tips

#### Code Not Working
```python
# Add print statements
print("Starting analysis...")
print(f"Data shape: {data.shape}")
print("Analysis complete!")

# Check data types
print(data.dtypes)
print(data.head())

# Use try-except blocks
try:
    result = risky_operation()
except Exception as e:
    print(f"Error: {e}")
    print("Check your data format")
```

---

## 🎯 Quick Start Checklists

### ✅ Local Jupyter Setup Checklist
- [ ] **Anaconda installed** and working
- [ ] **Project directory created** with organized structure
- [ ] **Jupyter notebook launched** from project directory
- [ ] **Data files uploaded** to data/ folder
- [ ] **Required packages installed** (scanpy, pandas, etc.)
- [ ] **First notebook created** and running

### ✅ Google Colab Setup Checklist
- [ ] **Google account** signed in
- [ ] **Colab accessed** at colab.research.google.com
- [ ] **Google Drive mounted** with drive.mount()
- [ ] **Data uploaded** to Google Drive
- [ ] **File paths working** from /content/drive/MyDrive/
- [ ] **Required packages installed** (usually pre-installed)

### ✅ First Analysis Checklist
- [ ] **Libraries imported** successfully
- [ ] **Data loaded** without errors
- [ ] **Basic exploration** completed (shape, head, describe)
- [ ] **Simple plot created** and displayed
- [ ] **Results saved** to appropriate location
- [ ] **Notebook documented** with markdown explanations

---

## 🚀 Next Steps

### Immediate Actions
1. **Choose your environment**: Local Jupyter vs Google Colab
2. **Set up organized project structure**
3. **Download and upload the proteomics dataset**
4. **Create your first analysis notebook**
5. **Run through basic data exploration**

### Continue Learning
- **Data Understanding**: Move to exploring your proteomics dataset
- **Statistical Analysis**: Begin with UPS protein analysis tutorial
- **Advanced Methods**: Progress to SQSTM1 and sliding window analyses

### Build Good Habits
- **Document everything**: Mix code with explanations
- **Save frequently**: Prevent losing work
- **Organize systematically**: Consistent file and folder naming
- **Test thoroughly**: Run "Restart & Run All" before sharing

---

**You're now ready to dive into computational proteomics analysis!** Whether you choose local Jupyter or Google Colab, you have the foundation to perform sophisticated biological data analysis.

*Next: [Understanding Your Proteomics Dataset](../03_data_understanding/dataset_overview.md)*

*Remember: The best way to learn Jupyter notebooks is by using them for real analysis - start with small examples and build up complexity!* 📓🔬
# 🐍 Python Installation Guide for Biologists

## 🎯 Why Python for Biology?

Python has become the **standard language** for biological data analysis because:
- **Easy to learn**: Syntax is close to natural language
- **Powerful libraries**: Specialized tools for every type of biological analysis
- **Free and open**: No expensive licenses needed
- **Huge community**: Millions of users worldwide, lots of help available
- **Reproducible**: Share exact methods with colleagues

**Don't worry if you've never programmed before!** This guide assumes zero experience.

---

## 💻 Installation Options (Choose One)

### Option 1: Anaconda (Recommended for Beginners) ⭐

**What is Anaconda?**
- Complete Python environment designed for data science
- Includes Python + all the libraries we need
- Includes Jupyter Notebook (interactive analysis environment)
- Easy installation and management

**Why choose Anaconda?**
- ✅ Everything included in one download
- ✅ No complex dependency management
- ✅ Works the same on Windows, Mac, Linux
- ✅ Easy to update and maintain

### Option 2: Standard Python + pip (For Some Experience)

**What is this approach?**
- Install basic Python, then add libraries as needed
- More flexible but requires more technical knowledge
- What many programmers use

**Why choose this?**
- ✅ Smaller initial download
- ✅ More control over what's installed
- ✅ Closer to how programming is typically done

### Option 3: Cloud-Based (No Installation Required)

**What are cloud options?**
- Run Python in your web browser
- No installation on your computer
- Good for trying things out

**Options include:**
- **Google Colab**: Free, powerful, integrates with Google Drive
- **Binder**: Run Jupyter notebooks from GitHub
- **GitHub Codespaces**: Full development environment in browser

---

## 🛠️ Step-by-Step Installation

### Method 1: Anaconda Installation (Recommended)

#### Step 1: Download Anaconda

1. **Go to**: https://www.anaconda.com/products/distribution
2. **Click**: "Download" button
3. **Choose your operating system**:
   - Windows: Download the .exe file
   - Mac: Download the .pkg file
   - Linux: Download the .sh file

**File sizes**: ~500MB (it's big because it includes everything!)

#### Step 2: Install Anaconda

**Windows**:
1. Double-click the downloaded .exe file
2. Follow the installation wizard
3. **Important**: Check "Add Anaconda to my PATH environment variable" (if offered)
4. Choose "Install for: Just Me" unless you need system-wide installation

**Mac**:
1. Double-click the downloaded .pkg file
2. Follow the installation wizard
3. Default settings are usually fine

**Linux**:
1. Open terminal
2. Navigate to download folder: `cd ~/Downloads`
3. Run: `bash Anaconda3-2023.XX-Linux-x86_64.sh` (replace XX with version)
4. Follow prompts, accept license, choose installation location

#### Step 3: Verify Installation

**Open Anaconda Navigator**:
- Windows: Start Menu → Anaconda3 → Anaconda Navigator
- Mac: Applications → Anaconda Navigator
- Linux: Type `anaconda-navigator` in terminal

**You should see**: A window with various application tiles (Jupyter Notebook, Spyder, etc.)

**If Navigator doesn't open**:
- Try opening "Anaconda Prompt" (Windows) or Terminal (Mac/Linux)
- Type: `conda list`
- You should see a list of installed packages

#### Step 4: Test Python

1. **Open Anaconda Prompt** (Windows) or Terminal (Mac/Linux)
2. **Type**: `python`
3. **You should see**:
   ```
   Python 3.11.X |Anaconda, Inc.| (default, ...)
   [GCC ...] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>>
   ```
4. **Type**: `2 + 2` and press Enter
5. **You should see**: `4`
6. **Type**: `exit()` to quit Python

**Congratulations!** Python is working.

### Method 2: Standard Python Installation

#### Step 1: Download Python

1. **Go to**: https://www.python.org/downloads/
2. **Download** the latest Python 3.x version (currently 3.11 or 3.12)
3. **Avoid** Python 2.x (it's outdated)

#### Step 2: Install Python

**Windows**:
1. Run the downloaded installer
2. **IMPORTANT**: Check "Add Python to PATH"
3. Choose "Install Now"

**Mac**:
1. Run the downloaded installer
2. Follow the installation wizard
3. Default settings are usually fine

**Linux** (Ubuntu/Debian):
```bash
sudo apt update
sudo apt install python3 python3-pip
```

#### Step 3: Install Required Packages

**Open terminal/command prompt and run**:
```bash
pip install pandas numpy scipy matplotlib seaborn scikit-learn jupyter scanpy statsmodels tqdm
```

**This might take 5-10 minutes** to download and install everything.

### Method 3: Cloud-Based Setup

#### Google Colab (Easiest)

1. **Go to**: https://colab.research.google.com/
2. **Sign in** with Google account
3. **Click**: "New Notebook"
4. **Test it**: Type `2 + 2` in a cell and press Shift+Enter

**Advantages**:
- ✅ No installation required
- ✅ Free access to powerful computers
- ✅ Most packages pre-installed

**Disadvantages**:
- ❌ Requires internet connection
- ❌ Files stored in Google Drive
- ❌ Session timeouts after inactivity

---

## 📦 Required Packages Explained

### Core Data Science Stack
- **pandas**: Excel-like data manipulation
- **numpy**: Mathematical operations on arrays
- **scipy**: Advanced statistics and scientific computing
- **matplotlib**: Basic plotting and visualization
- **seaborn**: Beautiful statistical visualizations

### Specialized for Our Analysis
- **scanpy**: Single-cell analysis (works great for proteomics too)
- **statsmodels**: Advanced statistical modeling
- **scikit-learn**: Machine learning algorithms
- **tqdm**: Progress bars (makes waiting more pleasant!)

### How to Install Missing Packages

**If using Anaconda**:
```bash
conda install package_name
```

**If using pip**:
```bash
pip install package_name
```

**In Jupyter Notebook**:
```python
!pip install package_name
```

---

## 🚀 Setting Up Your Development Environment

### Option 1: Jupyter Notebook (Recommended for Learning)

**What is Jupyter Notebook?**
- Interactive environment where you mix code, text, and plots
- Perfect for data analysis and learning
- Industry standard for exploratory analysis

**How to start Jupyter**:

**From Anaconda Navigator**:
1. Open Anaconda Navigator
2. Click "Launch" under Jupyter Notebook
3. Your web browser will open with Jupyter

**From command line**:
1. Open terminal/command prompt
2. Type: `jupyter notebook`
3. Browser should open automatically

**Create your first notebook**:
1. Click "New" → "Python 3"
2. Type: `print("Hello, proteomics!")`
3. Press Shift+Enter to run
4. You should see: `Hello, proteomics!`

### Option 2: Jupyter Lab (More Advanced Interface)

**What is Jupyter Lab?**
- Next generation of Jupyter Notebook
- More features, better file management
- Can work with multiple notebooks simultaneously

**How to start**:
```bash
jupyter lab
```

### Option 3: VS Code (For More Programming)

**What is VS Code?**
- Full-featured code editor
- Great for writing longer scripts
- Excellent Python support

**Installation**:
1. Download from: https://code.visualstudio.com/
2. Install Python extension
3. Configure Python interpreter

---

## 🧪 Testing Your Installation

### Basic Functionality Test

**Create a new Jupyter notebook and run each cell**:

```python
# Cell 1: Test basic Python
print("Python is working!")
result = 2 + 2
print(f"2 + 2 = {result}")
```

```python
# Cell 2: Test pandas (data manipulation)
import pandas as pd
data = {'protein': ['SQSTM1', 'VDAC1', 'TAU'],
        'expression': [10.5, 5.2, 8.1]}
df = pd.DataFrame(data)
print(df)
```

```python
# Cell 3: Test plotting
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x)
plt.plot(x, y)
plt.title("Test Plot - Sine Wave")
plt.show()
```

```python
# Cell 4: Test scientific computing
from scipy import stats
import numpy as np

# Generate some fake biological data
group1 = np.random.normal(5, 1, 20)  # Control group
group2 = np.random.normal(7, 1, 20)  # Treatment group

# Run a t-test
t_stat, p_value = stats.ttest_ind(group1, group2)
print(f"T-statistic: {t_stat:.3f}")
print(f"P-value: {p_value:.3f}")
```

```python
# Cell 5: Test scanpy (for proteomics analysis)
import scanpy as sc
print("Scanpy version:", sc.__version__)
print("All packages working correctly!")
```

**If all cells run without errors, your installation is successful!**

### Troubleshooting Common Issues

#### Issue 1: "Module not found" error
```
ModuleNotFoundError: No module named 'pandas'
```
**Solution**:
- Package not installed: `pip install pandas`
- Wrong Python environment: Make sure you're using the same Python where you installed packages

#### Issue 2: Jupyter won't start
```
Command 'jupyter' not found
```
**Solutions**:
- **If using Anaconda**: Use Anaconda Navigator instead
- **If using pip**: `pip install jupyter`
- **Path issue**: Restart terminal/command prompt

#### Issue 3: Permission errors (Windows)
```
Permission denied
```
**Solutions**:
- Run command prompt as administrator
- Install packages with `--user` flag: `pip install --user pandas`

#### Issue 4: SSL/Certificate errors
```
SSL: CERTIFICATE_VERIFY_FAILED
```
**Solutions**:
- Update certificates: On Mac, run `Install Certificates.command` in Python folder
- Use conda instead: `conda install package_name`

---

## 📁 Organizing Your Project

### Create Your Project Structure

**Create these folders on your computer**:
```
My_Proteomics_Analysis/
├── data/                    # Store datasets here
├── notebooks/               # Jupyter notebooks
├── scripts/                 # Python scripts (.py files)
├── results/                 # Analysis outputs
├── plots/                   # Generated figures
└── documentation/           # Your notes and reports
```

**In Jupyter, navigate to this folder before starting work.**

### Best Practices

#### File Naming
- ✅ Use descriptive names: `sqstm1_analysis_2024_01_15.ipynb`
- ✅ Include dates: Helps track versions
- ✅ Use underscores: Not spaces
- ❌ Avoid: `analysis.ipynb`, `temp.ipynb`, `untitled.ipynb`

#### Code Organization
- **One analysis per notebook**: Don't try to do everything in one file
- **Clear section headers**: Use markdown cells to explain what each section does
- **Save frequently**: Ctrl+S (Windows/Linux) or Cmd+S (Mac)

#### Version Control
- **Save different versions**: `analysis_v1.ipynb`, `analysis_v2.ipynb`
- **Backup important work**: Copy to cloud storage
- **Document changes**: Note what you changed between versions

---

## 🔧 Advanced Setup (Optional)

### Virtual Environments

**What are virtual environments?**
- Isolated Python installations for different projects
- Prevents package conflicts
- Professional best practice

**Create a virtual environment**:
```bash
# Using conda (if you have Anaconda)
conda create -n proteomics python=3.11
conda activate proteomics

# Using venv (if you have standard Python)
python -m venv proteomics_env
# Activate: proteomics_env\Scripts\activate (Windows)
# Activate: source proteomics_env/bin/activate (Mac/Linux)
```

### Package Management

**Keep track of your packages**:
```bash
# Save current packages
pip freeze > requirements.txt

# Install from requirements (on another computer)
pip install -r requirements.txt
```

### IDE Configuration

**If using VS Code**:
1. Install Python extension
2. Set Python interpreter: Ctrl+Shift+P → "Python: Select Interpreter"
3. Enable Jupyter support: Install Jupyter extension

---

## 🆘 Getting Help

### When Installation Goes Wrong

#### First Steps
1. **Restart your computer**: Fixes many path/environment issues
2. **Try a different method**: If Anaconda fails, try standard Python
3. **Check system requirements**: Make sure your OS is supported

#### Where to Get Help
1. **Anaconda documentation**: https://docs.anaconda.com/
2. **Python.org**: https://www.python.org/about/help/
3. **Stack Overflow**: Search for your exact error message
4. **Reddit r/Python**: Helpful community for beginners

#### Common Solutions
- **Use conda instead of pip**: `conda install package_name`
- **Update everything**: `conda update --all`
- **Fresh installation**: Uninstall and reinstall

### Learning Resources

#### Interactive Learning
- **Codecademy Python Course**: Hands-on learning
- **DataCamp**: Data science focus
- **Kaggle Learn**: Free micro-courses

#### Documentation
- **Official Python Tutorial**: https://docs.python.org/3/tutorial/
- **Pandas documentation**: https://pandas.pydata.org/docs/
- **Matplotlib tutorials**: https://matplotlib.org/stable/tutorials/index.html

---

## ✅ Installation Complete Checklist

Before moving to the next section, ensure:

- [ ] Python is installed and working
- [ ] All required packages are installed
- [ ] Jupyter Notebook/Lab starts successfully
- [ ] Test code runs without errors
- [ ] Project folder structure is created
- [ ] You know how to create and run a simple notebook

**If all items are checked, you're ready for data analysis!**

---

## 🚀 Next Steps

### Immediate Next Steps
1. **Explore Jupyter**: Create a few test notebooks
2. **Practice basics**: Try the examples from our test code
3. **Download our data**: Get the proteomics dataset
4. **Read data guide**: Understand what we'll be analyzing

### Learning Path
- **Beginners**: Start with basic Python tutorials
- **Some experience**: Jump into data manipulation with pandas
- **Ready to analyze**: Move on to data understanding section

---

## 💡 Pro Tips for Biologists

### Making Python Less Intimidating
1. **Start small**: Don't try to understand everything at once
2. **Focus on concepts**: The syntax will become familiar with practice
3. **Use it for real work**: Apply Python to your actual research questions
4. **Join communities**: Connect with other scientists using Python

### Python for Biological Thinking
- **Think of code as protocols**: Step-by-step instructions
- **Variables are like samples**: You label them and work with them
- **Functions are like techniques**: Reusable procedures
- **Packages are like equipment**: Specialized tools for specific tasks

### Building Confidence
- **Everyone makes errors**: Even experts debug code regularly
- **Google is your friend**: Most errors have been solved before
- **Practice regularly**: 15 minutes daily beats 2 hours weekly
- **Teach others**: Explaining concepts solidifies your understanding

---

**Congratulations on setting up your computational biology toolkit!** 🎉

You now have all the software needed to perform sophisticated proteomic analyses. The next step is understanding your data.

*Next: [Understanding Your Data](../03_data_understanding/dataset_overview.md)*

*Remember: The goal isn't to become a programmer, but to use programming tools to answer biological questions better!* 🧬🐍
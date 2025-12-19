# 🚀 Complete Setup Guide - Multiple Skill Levels

This guide provides setup instructions tailored to different experience levels, ensuring everyone can successfully use this proteomics analysis framework.

---

## 🎯 Choose Your Skill Level

### 🌱 **Beginner** (Never programmed before)
- Complete step-by-step instructions
- No assumptions about prior knowledge
- Extra troubleshooting and verification steps
- **Time needed**: 2-3 hours

### 🌿 **Intermediate** (Some programming experience)
- Streamlined instructions with key steps
- Focus on proteomics-specific setup
- **Time needed**: 30-60 minutes

### 🌳 **Advanced** (Experienced developer)
- Quick setup commands
- Customization options
- **Time needed**: 10-15 minutes

---

## 🌱 BEGINNER SETUP

### Step 1: Check Your Computer
First, let's see what you're working with:

**Windows:**
1. Press `Windows + R` keys together
2. Type `cmd` and press Enter
3. In the black window that opens, type: `python --version`
4. If you see "Python 3.x.x", you have Python! If not, continue to Step 2.

**Mac:**
1. Press `Cmd + Space` and type "Terminal"
2. Press Enter to open Terminal
3. Type: `python3 --version`
4. If you see "Python 3.x.x", you have Python! If not, continue to Step 2.

**Linux:**
1. Open Terminal (usually Ctrl+Alt+T)
2. Type: `python3 --version`
3. If you see "Python 3.x.x", you have Python! If not, continue to Step 2.

### Step 2: Install Python (if needed)

**Windows:**
1. Go to [python.org/downloads](https://python.org/downloads)
2. Click "Download Python 3.x.x" (latest version)
3. Run the downloaded file
4. ⚠️ **IMPORTANT**: Check "Add Python to PATH" box!
5. Click "Install Now"
6. Restart your computer when done

**Mac:**
1. Go to [python.org/downloads](https://python.org/downloads)
2. Click "Download Python 3.x.x" (latest version)
3. Run the downloaded .pkg file
4. Follow the installation wizard

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

### Step 3: Verify Installation
Open Terminal/Command Prompt again and type:
```bash
python3 --version
pip3 --version
```
You should see version numbers for both.

### Step 4: Download This Project
1. Go to the project's main page
2. Click the green "Code" button
3. Click "Download ZIP"
4. Unzip the file to your Desktop
5. Rename the folder to "proteomics_analysis"

### Step 5: Open Terminal in Project Folder

**Windows:**
1. Open File Explorer
2. Navigate to Desktop → proteomics_analysis
3. In the address bar, type `cmd` and press Enter

**Mac:**
1. Open Finder
2. Navigate to Desktop → proteomics_analysis
3. Right-click in the folder and select "New Terminal at Folder"

**Linux:**
1. Open file manager
2. Navigate to Desktop → proteomics_analysis
3. Right-click and select "Open in Terminal"

### Step 6: Install Required Packages
Copy and paste this command (it might take 5-10 minutes):
```bash
pip3 install -r requirements.txt
```

If you see errors, try:
```bash
python3 -m pip install -r requirements.txt
```

### Step 7: Test Everything Works
```bash
python3 tools/test_sqstm1_integration.py
```

If you see test results with ✅ symbols, congratulations! You're ready to start learning.

### 🆘 Beginner Troubleshooting
- **"Command not found"**: Python wasn't added to PATH. Reinstall Python and check the PATH box.
- **"Permission denied"**: Add `sudo` before the command on Mac/Linux.
- **"Package not found"**: Your internet connection might be slow. Try again.

---

## 🌿 INTERMEDIATE SETUP

### Quick Prerequisites Check
```bash
python3 --version  # Should be 3.7+
pip3 --version     # Should be 20.0+
git --version      # Optional but recommended
```

### Clone or Download Project
```bash
# Option 1: Git clone (recommended)
git clone [repository-url]
cd biologist_guide_proteomics

# Option 2: Download and extract ZIP
# Then navigate to the folder
```

### Virtual Environment (Recommended)
```bash
python3 -m venv proteomics_env
source proteomics_env/bin/activate  # Mac/Linux
# OR
proteomics_env\Scripts\activate     # Windows

pip install --upgrade pip
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Verify Installation
```bash
python -c "import pandas, numpy, matplotlib; print('Core packages OK')"
python tools/test_sqstm1_integration.py
```

### Optional: Jupyter Setup
```bash
pip install jupyter jupyterlab
jupyter lab --generate-config
```

---

## 🌳 ADVANCED SETUP

### One-Command Setup
```bash
git clone [repository-url] && cd biologist_guide_proteomics
python3 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip && pip install -r requirements.txt
python tools/test_sqstm1_integration.py
```

### Development Setup
```bash
# Install development dependencies
pip install -e .
pip install pytest black pylint mypy pre-commit

# Setup pre-commit hooks
pre-commit install

# Run tests
pytest tools/ -v
```

### Docker Setup (Optional)
```dockerfile
FROM python:3.9-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
CMD ["python", "tools/test_sqstm1_integration.py"]
```

### Custom Configuration
```python
# config.py
import os

# API rate limits
UNIPROT_RATE_LIMIT = float(os.getenv('UNIPROT_RATE_LIMIT', 1.0))
STRING_API_KEY = os.getenv('STRING_API_KEY', None)

# Visualization settings
FIGURE_DPI = int(os.getenv('FIGURE_DPI', 300))
PLOT_STYLE = os.getenv('PLOT_STYLE', 'seaborn-whitegrid')
```

---

## 📚 Next Steps by Skill Level

### 🌱 Beginners: Start Here
1. [README_START_HERE.md](00_getting_started/README_START_HERE.md)
2. [What is Computational Biology?](00_getting_started/what_is_computational_biology.md)
3. [Biology Primer](01_background_knowledge/biology_primer.md)
4. [Statistics for Biologists](01_background_knowledge/statistics_for_biologists.md)

### 🌿 Intermediate: Jump to Analysis
1. [Dataset Overview](03_data_understanding/dataset_overview.md)
2. [UPS Protein Analysis](04_group1_analyses/statement1_ups_proteins/)
3. [SQSTM1 Analysis](04_group1_analyses/statement2_sqstm1_upregulation/)
4. [Practice Exercises](07_practice_exercises/)

### 🌳 Advanced: Dive Deep
1. [Tools Documentation](tools/README.md)
2. [UniProt UPS Integration](tutorials/uniprot_ups_integration.md)
3. [Proteome-wide Analysis](05_group2_analyses/)
4. [Custom Analysis Development](06_resources_and_support/comprehensive_resource_guide.md)

---

## 🔧 Platform-Specific Notes

### Windows Specifics
- Use `python` instead of `python3` after installation
- Use `pip` instead of `pip3`
- Use backslashes `\` in file paths
- Consider Windows Subsystem for Linux (WSL) for advanced usage

### Mac Specifics
- Homebrew alternative: `brew install python3`
- May need Xcode command line tools: `xcode-select --install`
- Use `python3` and `pip3` explicitly

### Linux Specifics
- Package manager installations available
- May need development headers: `sudo apt install python3-dev`
- Virtual environments recommended for system separation

---

## 🆘 Universal Troubleshooting

### Common Issues and Solutions

#### Package Installation Fails
```bash
# Try upgrading pip first
python3 -m pip install --upgrade pip

# Install packages one by one to isolate issues
pip install pandas
pip install numpy
# etc.
```

#### Import Errors
```bash
# Check Python path
python3 -c "import sys; print(sys.path)"

# Verify package installation
pip list | grep pandas
```

#### Memory Issues
```bash
# Install packages with no cache
pip install --no-cache-dir -r requirements.txt

# For large datasets, increase virtual memory
```

#### Permission Issues
```bash
# Use user installation
pip install --user -r requirements.txt

# On Mac/Linux, use sudo only if necessary
sudo pip3 install -r requirements.txt
```

### Getting Help
1. Check [Troubleshooting Guide](06_resources_and_support/troubleshooting_guide.md)
2. Review error messages carefully
3. Search for error messages online
4. Post issues with full error messages and system information

---

## ✅ Verification Checklist

### All Skill Levels
- [ ] Python 3.7+ installed and accessible
- [ ] pip/pip3 working and updated
- [ ] All requirements.txt packages installed successfully
- [ ] Test script runs without errors
- [ ] Can import key packages (pandas, numpy, matplotlib)

### Intermediate+
- [ ] Virtual environment activated
- [ ] Git repository cloned (if using git)
- [ ] Jupyter/JupyterLab working (if needed)

### Advanced
- [ ] Development environment configured
- [ ] Pre-commit hooks installed (if contributing)
- [ ] Custom configuration set up
- [ ] Docker image built (if using Docker)

---

**🎉 Congratulations!** You're now ready to start your computational proteomics journey. Choose your learning path based on your skill level and begin exploring the fascinating world of protein analysis!
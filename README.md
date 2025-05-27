# BEC_3D Fitting Framework

This repository contains code and scripts for fitting 3D Bose-Einstein correlation (BEC) functions to experimental or simulated data using ROOT and Python.

## Overview

- The main fitting logic is implemented in [`fit.cxx`](fit.cxx), which performs the 3D fit and outputs results and plots.
- The Python script [`do_fit.py`](do_fit.py) serves as the main interface for running fits, handling arguments, and preparing input/output.
- The process can be automated for multiple configurations using the shell script [`do_shortfit.sh`](do_shortfit.sh).

## Structure

- **fit.cxx**: Main C++ code for fitting, histogram handling, and output generation.
- **do_fit.py**: Python script to run the fit, parse arguments, and manage input/output files.
- **do_shortfit.sh**: Bash script to automate running [`do_fit.py`](do_fit.py) with various arguments for batch processing.
- **input/**: Directory for input configuration files (e.g., YAML).
- **output/**: Directory for output files (ROOT, CSV, plots, etc.).
- **plot_summary/**: Additional plotting and utility scripts.

## Usage

### 1. Running a Single Fit

You can run a fit using the Python script:

```sh
python do_fit.py [arguments]
```

This script will prepare the necessary arguments and call the compiled `fit` executable (from `fit.cxx`).

**Example:**
```sh
python do_fit.py --file-data input/data.root --file-mc input/mc.root --hist-data myhist --hist-mc myhist_mc --q-min 0.02 --q-max 1.2 --c2-index 1
```

### 2. Automating Multiple Fits

To automate fits over multiple configurations, use the shell script:

```sh
./do_shortfit.sh [arguments]
```

This script will loop over various settings and call `do_fit.py` as needed.

### 3. Compiling the C++ Code

Make sure you have ROOT installed. Then compile with:

```sh
make
```

This will produce the `fit` executable.

## Arguments

Both `do_fit.py` and `fit.cxx` accept a variety of arguments for specifying input files, histogram names, fit ranges, rejection regions, and fit options. See the help output for each script for details:

```sh
python do_fit.py --help
./fit --help
```

## Input/Output

- **Input:** Data and MC ROOT files, YAML configuration files.
- **Output:** ROOT files, CSV tables, LaTeX tables, and plots (PDF/PNG) in the `output/` directory.

## Automation

- Use `do_shortfit.sh` for batch processing.
- Edit or provide input YAML files in the `input/` directory as needed.

## Requirements

- [ROOT](https://root.cern/) (for C++ code and histogramming)
- Python 3 (for scripting and automation)
- `pyyaml` (for YAML parsing in Python scripts)

## License

This project is intended for academic and research use. Please contact the author for other uses.

---

**Author:** happyguy54  
**Contact:** [GitHub Profile](https://github.com/happyguy54)

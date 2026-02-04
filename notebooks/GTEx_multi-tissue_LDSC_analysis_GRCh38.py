import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def __():

    import marimo as mo
    import urllib.request
    import os
    import subprocess
    import pandas as pd
    import tarfile
    from pathlib import Path
    import glob
    import json
    return mo, urllib, os, subprocess, pd, tarfile, Path, glob, json


@app.cell
def __(mo):
    mo.md("""
# LDSC Multi-Tissue Analysis – GTEx Gene Sets

This notebook performs tissue-specific heritability analysis using LDSC (LD Score Regression)  
with GTEx tissue-specific gene sets. The analysis identifies which tissues contribute most  
to the heritability of complex traits.

All analyses are performed using **GRCh38 / hg38** coordinates.

## Analysis Overview:
1. Process tissue-specific top gene lists
2. Create genomic annotations for each tissue
3. Calculate tissue-specific LD scores
4. Run partitioned heritability analysis
""")
    return


@app.cell
def __(mo):
    mo.md("""
## 0. Setup: Download and configure LDSC

This cell will:
1. Create a Python 2.7 conda environment for LDSC
2. Install LDSC and its dependencies
3. Install BEDTools (required for annotation generation)
""")
    return


@app.cell
def __(Path, subprocess, os, json):
    TOOLS_DIR = Path("tools")
    LDSC_DIR = TOOLS_DIR / "ldsc"
    TOOLS_DIR.mkdir(exist_ok=True)
    
 
    env_check = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True,
        text=True
    )
    ldsc_env_exists = "ldsc27" in env_check.stdout
    
    if not ldsc_env_exists:
        print("Creating Python 2.7 conda environment for LDSC...")
        subprocess.run([
            "conda", "create", "-n", "ldsc27", 
            "python=2.7", "-y"
        ], check=True)
    
  
    conda_prefix = subprocess.run(
        ["conda", "env", "list", "--json"],
        capture_output=True,
        text=True,
        check=True
    )
    envs = json.loads(conda_prefix.stdout)["envs"]
    ldsc27_path = [e for e in envs if "ldsc27" in e][0]
    
   
    bedtools_check = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "bedtools"), "--version"],
        capture_output=True,
        text=True
    )
    
    if bedtools_check.returncode != 0:
        print("Installing BEDTools in ldsc27 environment...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27",
            "-c", "bioconda", "bedtools", "-y"
        ], check=True)
    else:
        print(f" BEDTools already installed: {bedtools_check.stdout.strip()}")
    
 
    if not LDSC_DIR.exists():
        print("Cloning LDSC repository...")
        subprocess.run([
            "git", "clone",
            "https://github.com/bulik/ldsc.git",
            str(LDSC_DIR)
        ], check=True)
    
    check_numpy = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "python"), "-c", "import numpy"],
        capture_output=True
    )
    
    if check_numpy.returncode != 0:
        print("Installing LDSC dependencies in ldsc27 environment...")
        print("Installing OpenSSL 1.0...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "openssl=1.0.2", "-c", "conda-forge"
        ], check=True)
        
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "numpy", "scipy", "pandas", "bitarray", "-c", "conda-forge"
        ], check=True)
        
        print("Installing pybedtools and pysam...")
        subprocess.run([
            "conda", "install", "-n", "ldsc27", "-y",
            "pybedtools", "pysam=0.15.3", "-c", "bioconda", "-c", "conda-forge"
        ], check=True)
        
        print("Dependencies installed!")
    
    print("\nVerifying BEDTools is accessible to pybedtools...")
    pybedtools_check = subprocess.run(
        [os.path.join(ldsc27_path, "bin", "python"), "-c", 
         "from pybedtools import BedTool; import pybedtools.helpers as helpers; print('BEDTools path:', helpers.get_bedtools_path())"],
        capture_output=True,
        text=True
    )
    
    if pybedtools_check.returncode != 0:
        print("WARNING: pybedtools can't find BEDTools!")
        print("Error:", pybedtools_check.stderr)
    else:
        print("pybedtools can access BEDTools")
        print(pybedtools_check.stdout)
    
    ldsc_script = LDSC_DIR / "ldsc.py"
    subprocess.run(["chmod", "+x", str(ldsc_script)], check=True)
    
    python27_path = os.path.join(ldsc27_path, "bin", "python")
    
    print(f"\n LDSC environment ready!")
    print(f"Python 2.7 path: {python27_path}")
    print(f"LDSC path: {ldsc_script}")
    
    return (ldsc_script, python27_path, ldsc27_path)



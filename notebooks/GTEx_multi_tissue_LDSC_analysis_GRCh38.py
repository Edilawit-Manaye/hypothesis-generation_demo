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

    @app.cell
def __(mo):
    mo.md("## 1. Download GTEx gene sets, reference panels, and GWAS summary statistics")
    return


@app.cell
def __(os, urllib):

    os.makedirs("data/gtex/tissue_gene_sets", exist_ok=True)
    os.makedirs("data/reference", exist_ok=True)
    os.makedirs("data/gwas", exist_ok=True)
    os.makedirs("data/gtex/annot_files", exist_ok=True)
    os.makedirs("data/gtex/ldscores", exist_ok=True)
    
    print(" Data directories created")
    return


@app.cell
def __(mo):

    mo.md("""
### GTEx Tissue Gene Sets

**Important**: Place your GTEx tissue-specific gene lists in `data/gtex/tissue_gene_sets/`

Expected file format: `{tissue_name}_top10.txt` (or your chosen suffix)

Each file should contain:
- One gene per line (gene symbols or Ensembl IDs)
- Top genes expressed in that tissue
- Example tissues: Brain_Cortex, Liver, Heart, etc.

You can also modify the code below to specify custom file paths.
""")
    return


@app.cell
def __(os, urllib):

    if not os.path.exists("data/reference/GRCh38.tgz"):
        print("\nDownloading GRCh38 reference with baseline LD scores (this may take several minutes)...")
        urllib.request.urlretrieve(
            "https://zenodo.org/records/10515792/files/GRCh38.tgz?download=1",
            "data/reference/GRCh38.tgz"
        )
        print(" GRCh38 reference downloaded")
    else:
        print("\n GRCh38 reference already downloaded")
    
    if not os.path.exists("data/gwas/AD_bellenguez_2022_hg38.tsv.gz"):
        print("Downloading example GWAS summary statistics (Alzheimer's disease)...")
        print("You can replace this with your trait of interest.")
        urllib.request.urlretrieve(
            "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz",
            "data/gwas/AD_bellenguez_2022_hg38.tsv.gz"
        )
        print(" GWAS data downloaded")
    else:
        print(" GWAS data already downloaded")
    
    print("\n All downloads complete!")
    return


@app.cell
def __(mo):

    mo.md("## 2. Extract reference LD panels and baseline LD scores")
    return


@app.cell  
def __(tarfile, os):
    print("Checking GRCh38.tgz contents...")
    
    if os.path.exists("data/reference/GRCh38.tgz"):
        with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
            members = _tar.getmembers()
            print(f"Archive contains {len(members)} items")
            print("\nTop-level structure:")
            seen = set()
            for m in members[:50]:  
                parts = m.name.split('/')
                if len(parts) > 1:
                    top = parts[0] + "/" + parts[1]
                    if top not in seen:
                        print(f"  {top}")
                        seen.add(top)
    return


@app.cell
def __(tarfile, os):

    print("Extracting GRCh38 reference files...")
    
    if not os.path.exists("data/reference/GRCh38"):
        try:
            print("  Verifying GRCh38.tgz integrity...")
            with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
                _tar.getmembers()
            print("  File integrity verified")
        except (EOFError, tarfile.ReadError) as e:
            print(f"\n  Error: GRCh38.tgz is corrupted!")
            print(f"  Please delete it and re-run: rm data/reference/GRCh38.tgz")
            raise
        
        print("  Extracting (this may take a few minutes)...")
        with tarfile.open("data/reference/GRCh38.tgz", "r:gz") as _tar:
            _tar.extractall("data/reference")
        print("GRCh38 reference extracted successfully")
    else:
        print("GRCh38 directory exists")
    
  
    nested_files = [
        ("data/reference/GRCh38/baselineLD_v2.2.tgz", "data/reference", "baselineLD_v2.2"),
        ("data/reference/GRCh38/plink_files.tgz", "data/reference/GRCh38", "1000G.EUR.hg38.1.bed"),
        ("data/reference/GRCh38/weights.tgz", "data/reference/GRCh38", "weights")
    ]
    
    for tar_file, extract_to, check_file in nested_files:
        check_path = os.path.join(extract_to, check_file)
        if os.path.exists(check_path):
            print(f"  {os.path.basename(tar_file)} already extracted")
            continue
            
        if os.path.exists(tar_file):
            print(f"  Extracting {os.path.basename(tar_file)}...")
            with tarfile.open(tar_file, "r:gz") as _tar:
                _tar.extractall(extract_to)
            print(f" {os.path.basename(tar_file)} extracted")
            
            if not os.path.exists(check_path):
                print(f" Warning: Expected file {check_path} not found after extraction")
        else:
            print(f"  Warning: {tar_file} not found")
    
    critical_file = "data/reference/GRCh38/1000G.EUR.hg38.1.bim"
    if os.path.exists(critical_file):
        print(f"\n All reference files ready! Verified: {critical_file}")
    else:
        print(f"\n ERROR: Critical file missing: {critical_file}")
        print("Checking what files exist in data/reference/GRCh38/:")
        if os.path.exists("data/reference/GRCh38"):
            files = os.listdir("data/reference/GRCh38")
            print(f"  Found {len(files)} files/directories")
            for _f in sorted(files)[:10]:  
                print(f"    - {_f}")
        else:
            print("  Directory doesn't exist!")
    
    return

@app.cell
def __(mo):
    """
    COMMIT: Add GWAS munging section header
    """
    mo.md("## 3. Munge GWAS summary statistics")
    return


@app.cell
def __(subprocess, os, python27_path):
    os.makedirs("data/munged", exist_ok=True)

    munged_file = "data/munged/AD_bellenguez_2022_hg38_munged.sumstats.gz"
    
    if os.path.exists(munged_file):
        print(" Munged GWAS file already exists, skipping munging step")
    else:
        print("\n" + "="*60)
        print("STEP 3: Munging GWAS summary statistics")
        print("="*60)
        print("\nNote: Modify column names below to match your GWAS format!")
        
        subprocess.run([
            python27_path, "tools/ldsc/munge_sumstats.py",
            "--sumstats", "data/gwas/AD_bellenguez_2022_hg38.tsv.gz",
            "--out", "data/munged/AD_bellenguez_2022_hg38_munged",
            "--a1", "effect_allele",
            "--a2", "other_allele",
            "--p", "p_value",
            "--snp", "variant_id",
            "--N-col", "n_total"
        ], check=True)
        
        print("\n GWAS munging complete")
    return
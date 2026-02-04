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
import pandas as pd
import os

base_dir = "/mnt/hdd_1/rediet/hypothesis-generation-demo"

# Check first 2 chromosomes to understand the structure
for chr_num in [1, 2]:
    print(f"=== Chromosome {chr_num} ===")
    
    # Read your ATAC annotation file
    atac_file = f"{base_dir}/ldsc_ph/ATAC_annot/chr{chr_num}.annot.gz"
    atac_df = pd.read_csv(atac_file, sep='\t')
    
    print(f"ATAC annotation file:")
    print(f"  Columns: {atac_df.columns.tolist()}")
    print(f"  Shape: {atac_df.shape}")
    print(f"  First few rows:")
    print(atac_df.head(3))
    print()
    
    # Read the corresponding 1000G bim file
    bim_file = f"{base_dir}/finemapping_susie_abdu_v2/finemapping/data/1000Genomes_phase3/plink_format_b37/EUR/EUR.{chr_num}.1000Gp3.20130502.bim"
    bim_df = pd.read_csv(bim_file, sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    
    print(f"1000G BIM file:")
    print(f"  Shape: {bim_df.shape}")
    print(f"  First few rows:")
    print(bim_df.head(3))
    print()
    
    # Check SNP overlap
    atac_snps = set(atac_df['SNP'])
    bim_snps = set(bim_df['SNP'])
    common_snps = atac_snps.intersection(bim_snps)
    
    print(f"SNP overlap analysis:")
    print(f"  ATAC SNPs: {len(atac_snps)}")
    print(f"  1000G SNPs: {len(bim_snps)}")
    print(f"  Common SNPs: {len(common_snps)}")
    print(f"  ATAC SNPs missing from 1000G: {len(atac_snps - bim_snps)}")
    print(f"  1000G SNPs missing from ATAC: {len(bim_snps - atac_snps)}")
    print()
    
    # Check if ATAC file has the right columns for LDSC
    required_cols = ['CHR', 'BP', 'SNP', 'CM']
    has_required = all(col in atac_df.columns for col in required_cols)
    print(f"Has required LDSC columns {required_cols}: {has_required}")
    
    if has_required:
        # Check if SNPs are in same order
        first_10_atac = atac_df['SNP'].head(10).tolist()
        first_10_bim = bim_df['SNP'].head(10).tolist()
        print(f"First 10 SNPs match: {first_10_atac == first_10_bim}")
        
        if first_10_atac != first_10_bim:
            print(f"  ATAC SNPs: {first_10_atac}")
            print(f"  BIM SNPs:  {first_10_bim}")
    
    print("="*50)
    print()
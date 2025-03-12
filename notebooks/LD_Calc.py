import os
import subprocess
import pandas as pd
from cyvcf2 import VCF, Writer

<<<<<<< Updated upstream
POPULATION = "EUR"
SAMPLE_PANEL_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR = "../data/susie"
SNP_LIST_FILE = "../data/susie/snplist/chr_sig_locus.snplist.txt"
=======
CHROMOSOME = "1"
POPULATION = "EUR"
#VCF_URL = f"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{CHROMOSOME}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
VCF_URL=F"ALL.chr{CHROMOSOME}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
SAMPLE_PANEL_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR = "../data/susie/ALL_chr/plink_b_files/new_chr_eur_38"
SNP_LIST_FILE = "../data/susie/chr/snplist/chr_sig_locus_38.snplist.txt"
>>>>>>> Stashed changes

def run_command(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {cmd}")
        print(result.stderr)
        exit(1)
    return result

def prepare_ld_matrix():
<<<<<<< Updated upstream
    print("Running plink...")
    
    
 
    vcf_dir = os.path.join(OUTPUT_DIR, "vcf")
    updated_vcf_dir = os.path.join(OUTPUT_DIR, "updated_vcf")
    plink_binary_dir = os.path.join(OUTPUT_DIR, "plink_binary")
    

    os.makedirs(vcf_dir, exist_ok=True)
    os.makedirs(updated_vcf_dir, exist_ok=True)
    os.makedirs(plink_binary_dir, exist_ok=True)
    
=======
    print("running plink")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
>>>>>>> Stashed changes
    if not os.path.exists("integrated_call_samples_v3.20130502.ALL.panel"):
        run_command(f"wget {SAMPLE_PANEL_URL}")

    panel = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", sep="\t")
    eur_samples = panel[panel["super_pop"] == POPULATION]["sample"].tolist()

    with open(f"{OUTPUT_DIR}/eur_samples.txt", "w") as f:
        f.write("\n".join([f"{s}\t{s}" for s in eur_samples]))

<<<<<<< Updated upstream
    for chrom in range(15, 16):  
        CHROMOSOME = str(chrom)
        VCF_URL = f"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{CHROMOSOME}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
       
        vcf_file = os.path.join(vcf_dir, f"ALL.chr{CHROMOSOME}.vcf.gz")
        if not os.path.exists(vcf_file):
            run_command(f"wget {VCF_URL} -O {vcf_file}")
        
        updated_vcf_file = os.path.join(updated_vcf_dir, f"ALL.chr{CHROMOSOME}.updated.vcf.gz")
        print(f"Processing chromosome {CHROMOSOME}...")

        vcf = VCF(vcf_file)
        writer = Writer(updated_vcf_file, vcf)
        for variant in vcf:
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT[0]
            variant.ID = f"{chrom}:{pos}:{ref}:{alt}"
            writer.write_record(variant)
        writer.close()

        plink_prefix = os.path.join(plink_binary_dir, f"chr{CHROMOSOME}_eur")
        if not os.path.exists(f"{plink_prefix}.bed"):
            run_command(
                f"plink --vcf {updated_vcf_file} "
                f"--keep {OUTPUT_DIR}/eur_samples.txt "
                f"--make-bed --out {plink_prefix} "
            )
        
 
        filtered_prefix = f"{plink_binary_dir}/chr{CHROMOSOME}_eur_filtered"
        run_command(
            f"plink --bfile {plink_prefix} "
            f"--extract {SNP_LIST_FILE} "
            f"--make-bed --out {filtered_prefix}"
        )
=======
    # vcf_file = f"{OUTPUT_DIR}/ALL.chr{CHROMOSOME}.vcf.gz"
    # if not os.path.exists(vcf_file):
    #     run_command(f"wget {VCF_URL} -O {vcf_file}")
    
    # vcf = VCF(vcf_file)
    # output_vcf_file = f"{OUTPUT_DIR}/ALL.chr{CHROMOSOME}.updated.vcf.gz"
    # writer = Writer(output_vcf_file, vcf)

    # for variant in vcf:
    #     chrom = variant.CHROM
    #     pos = variant.POS
    #     ref=variant.REF
    #     alt=variant.ALT[0]
    #     variant.ID =  f"{chrom}:{pos}:{ref}:{alt}"
    #     writer.write_record(variant)

    # writer.close()

    output_vcf_file='../data/susie/ALL_chr/chr_eur_38/ALL.chr22.updated.vcf.gz'
    print("output vcf file")
    
    plink_prefix = f"{OUTPUT_DIR}/chr22_eur"
    if not os.path.exists(f"{plink_prefix}.bed"):
        run_command(
            f"plink --vcf {output_vcf_file} "
            f"--keep {OUTPUT_DIR}/eur_samples.txt "
            f"--make-bed --out {plink_prefix} "
            
        )
       #zgrep "##reference" ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    # filtered_prefix = f"{OUTPUT_DIR}/chr{CHROMOSOME}_eur_filtered"
    # run_command(
    #     f"plink --bfile {plink_prefix} "
    #     f"--extract {SNP_LIST_FILE} "
    #     f"--make-bed --out {filtered_prefix}"
    # )
    # run_command(
    # f"plink --bfile {plink_prefix} "
    # f"--keep-allele-order "
    # f"--r square "
    # f"--extract {SNP_LIST_FILE} "
    # f"--out ../data/ld/sig_locus_mt"
    # )

    # run_command(
    # f"plink --bfile {plink_prefix} "
    # f"--keep-allele-order "
    # f"--r2 square "
    # f"--extract {SNP_LIST_FILE} "
    # f"--out ../data/ld/sig_locus_mt_r2"
    # )
>>>>>>> Stashed changes

        # # Generate LD matrices
        # run_command(
        #     f"plink --bfile {filtered_prefix} "
        #     "--r square "
        #     f"--out {OUTPUT_DIR}/chr{CHROMOSOME}_ld_matrix"
        # )

        # print(f"LD matrix for chromosome {CHROMOSOME} saved to: {OUTPUT_DIR}/chr{CHROMOSOME}_ld_matrix.ld")

if __name__ == "__main__":
    prepare_ld_matrix()

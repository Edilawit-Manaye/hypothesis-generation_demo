import os
import subprocess
import pandas as pd
from cyvcf2 import VCF, Writer

POPULATION = "EUR"
SAMPLE_PANEL_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR = "../data/susie"
SNP_LIST_FILE = "../data/susie/snplist/chr_sig_locus.snplist.txt"
# SNP_LIST_FILE="../data/susie/UK_Biobank/COL/significant_snp_lists.txt"
#SNP_LIST_FILE="../data/susie/UK_Biobank/BMI/significant_snplist.txt"


def run_command(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {cmd}")
        print(result.stderr)
        exit(1)
    return result

def prepare_ld_matrix():
    print("Running plink...")
    
    vcf_dir = os.path.join(OUTPUT_DIR, "vcf")
    updated_vcf_dir = os.path.join(OUTPUT_DIR, "updated_vcf")
    plink_binary_dir = os.path.join(OUTPUT_DIR, "new_plink_binary")
    
    os.makedirs(vcf_dir, exist_ok=True)
    os.makedirs(updated_vcf_dir, exist_ok=True)
    os.makedirs(plink_binary_dir, exist_ok=True)

    if not os.path.exists("integrated_call_samples_v3.20130502.ALL.panel"):
        run_command(f"wget {SAMPLE_PANEL_URL}")

    panel = pd.read_csv("integrated_call_samples_v3.20130502.ALL.panel", sep="\t")
    eur_samples = panel[panel["super_pop"] == POPULATION]["sample"].tolist()

    with open(f"{OUTPUT_DIR}/eur_samples.txt", "w") as f:
        f.write("\n".join([f"{s}\t{s}" for s in eur_samples]))

    for chrom in range(1, 23):  
        CHROMOSOME = str(chrom)
        VCF_URL = f"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{CHROMOSOME}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
        vcf_file = os.path.join(vcf_dir, f"ALL.chr{CHROMOSOME}.vcf.gz")
        if not os.path.exists(vcf_file):
            run_command(f"wget {VCF_URL} -O {vcf_file}")
        
        updated_vcf_file =  vcf_file = os.path.join(updated_vcf_dir, f"ALL.chr{CHROMOSOME}.updated.vcf.gz")
        print(f"Processing chromosome {CHROMOSOME}...")

        # vcf = VCF(vcf_file)
        # writer = Writer(updated_vcf_file, vcf)
        # for variant in vcf:
        #     chrom = variant.CHROM
        #     pos = variant.POS
        #     ref = variant.REF
        #     alt = variant.ALT[0]
        #     variant.ID = f"{chrom}:{pos}:{ref}:{alt}"
        #     writer.write_record(variant)
        # writer.close()

        plink_prefix = os.path.join(plink_binary_dir, f"chr{CHROMOSOME}_eur")
        print("plink_prefix")
        if not os.path.exists(f"{plink_prefix}.bed"):
                    run_command(
                        f"plink --vcf {updated_vcf_file} "
                        f"--keep {OUTPUT_DIR}/eur_samples.txt "
                        f"--keep-allele-order "
                        f"--make-bed --out {plink_prefix}"
                    )
                    
        filtered_prefix = f"{plink_binary_dir}/chr{CHROMOSOME}_eur_filtered"

        run_command(
                    f"plink --bfile {plink_prefix} "
                    f"--extract {SNP_LIST_FILE} "
                    f"--keep-allele-order "
                    f"--make-bed --out {filtered_prefix}"
                )

       
    # for chr_file in os.listdir("../data/susie/ALL_chr/"):
    #         chr_num = chr_file.split("_")[0]
    #         ld_output_file = f"../data/susie/ALL_chr/ld/test_sig_locus_mt_{chr_num}.ld"
    #         r2_output_file = f"../data/susie/ALL_chr/ld/test_sig_locus_mt_r2_{chr_num}.ld"

    #         if not os.path.exists(ld_output_file):
    #             run_command(
    #                 f"plink --bfile {filtered_prefix} "
    #                 f"--keep-allele-order --r square "
    #                 f"--extract ../data/susie/ALL_chr/{chr_file} "
    #                 f"--out ../data/susie/ALL_chr/ld/test_sig_locus_mt_{chr_num}"
    #             )

    #         if not os.path.exists(r2_output_file):
    #             run_command(
    #                 f"plink --bfile {filtered_prefix} "
    #                 f"--keep-allele-order --r2 square "
    #                 f"--extract ../data/susie/ALL_chr/{chr_file} "
    #                 f"--out ../data/susie/ALL_chr/ld/test_sig_locus_mt_r2_{chr_num}"
    #             )
        
if __name__ == "__main__":
    prepare_ld_matrix()
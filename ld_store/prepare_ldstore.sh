#!/bin/bash
set -euo pipefail

PLINK_PREFIX="/mnt/hdd_1/rediet/hypothesis-generation-demo/finemapping_susie_abdu_v2/data/1000Genomes_phase3/plink_format_b37/EUR/EUR.16.1000Gp3.20130502"       # prefix for .bed/.bim/.fam (without extension)
GWAS_FILE="/mnt/hdd_1/rediet/hypothesis-generation-demo/21001_munged.gwas.imputed_v3.both_sexes.tsv"      # your GWAS summary stats
SAMPLES=359983         
OUT_PREFIX="BGEN_data"         
LDSTORE="./ldstore_v2.0_x86_64" 


echo "Step 1. Convert PLINK binary to BGEN..."
plink --bfile "${PLINK_PREFIX}" \
       --make-bgen \
       --out "${OUT_PREFIX}"

echo "Step 2. Make .z file from GWAS summary stats..."
python3 <<EOF
import pandas as pd
gwas = pd.read_csv("${GWAS_FILE}", sep="\t")
z = gwas[["SNP", "CHR", "BP", "A1", "A2"]]
z.to_csv("${OUT_PREFIX}.z", sep=" ", index=False, header=False)
print("Saved SNP info to ${OUT_PREFIX}.z")
EOF

echo "Step 3. Create master file..."
cat <<EOT > master.txt
z;bgen;bgi;bcor;ld;n_samples
${OUT_PREFIX}.z;${OUT_PREFIX}.bgen;${OUT_PREFIX}.bgen.bgi;${OUT_PREFIX}.bcor;${OUT_PREFIX}.ld;${SAMPLES}
EOT

echo "Step 4. Run LDstore to compute LD..."
${LDSTORE} --in-files master.txt --write-bcor --read-only-bgen

echo "Step 5. (Optional) Convert .bcor to human-readable .ld..."
${LDSTORE} --bcor-to-text --bcor-file ${OUT_PREFIX}.bcor --ld-file ${OUT_PREFIX}.ld

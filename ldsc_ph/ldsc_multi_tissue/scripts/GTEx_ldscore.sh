#!/bin/bash


LDSC_PATH="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/ldsc_repo/ldsc.py"
GENE_SETS="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/tissue_gene_sets/*_top10.txt"
BFILE_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/GRCh38/plink_files/plink_files"
ANNOT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/annot_files"
OUT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/ldscores"


mkdir -p "$OUT_DIR"

python "$LDSC_PATH" \
    --l2 \
    --bfile "${BFILE_DIR}/1000G.EUR.hg38.${chr}" \
    --ld-wind-cm 1 \
    --annot "${ANNOT_DIR}/${t_name}.${chr}.annot.gz" \
    --thin-annot \
    --out "${OUT_DIR}/${t_name}.${chr}"
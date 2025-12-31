#!/bin/bash

LDSC_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/ldsc_repo"
BFILE_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/data/ldsc_plink/plink_format_b38"
ANNOT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch_38_annotaions"
OUTPUT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch38_ld_scores_updated"
SNPLIST_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch38_ld_scores"

CELL="Ast"

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)


mkdir -p ${OUTPUT_DIR}

echo "Processing Ast cell type for all chromosomes sequentially..."
echo "=========================================================="


for CHR in "${CHROMOSOMES[@]}"; do
    echo ""
    echo "Starting Ast chr${CHR} at $(date '+%H:%M:%S')..."
    echo "---------------------------------------------------"
    

    if [ ! -f "${BFILE_DIR}/EUR.chr${CHR}.GRCh38.hm3.bim" ]; then
        echo "ERROR: Filtered bfile not found: EUR.chr${CHR}.GRCh38.hm3"
        echo "Skipping chr${CHR}..."
        continue
    fi
    
  
    if [ ! -f "${ANNOT_DIR}/Ast.chr${CHR}.annot.gz" ]; then
        echo "ERROR: Annotation file not found: Ast.chr${CHR}.annot.gz"
        echo "Skipping chr${CHR}..."
        continue
    fi
    
    if [ ! -f "${SNPLIST_DIR}/hm3.chr${CHR}.snplist.filtered" ]; then
        echo "ERROR: SNP list not found: hm3.chr${CHR}.snplist.filtered"
        echo "Skipping chr${CHR}..."
        continue
    fi
    

    cd ${LDSC_DIR}
    
    python ./ldsc.py \
        --print-snps ${SNPLIST_DIR}/hm3.chr${CHR}.snplist.filtered \
        --ld-wind-cm 1.0 \
        --out ${OUTPUT_DIR}/Ast.chr${CHR} \
        --bfile ${BFILE_DIR}/EUR.chr${CHR}.GRCh38.hm3 \
        --thin-annot \
        --yes-really \
        --annot ${ANNOT_DIR}/Ast.chr${CHR}.annot.gz \
        --l2
    

    if [ $? -eq 0 ]; then
        echo "Ast chr${CHR} completed successfully"
    else
        echo "Ast chr${CHR} failed"
    fi
    
    echo "Finished Ast chr${CHR} at $(date '+%H:%M:%S')"
    echo ""
done

echo "=========================================================="
echo "All chromosomes processed for Ast cell type!"
echo "Output files in: ${OUTPUT_DIR}/Ast.chr*.l2.ldscore.gz"
echo "=========================================================="
#!/bin/bash


cd /mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/


MAKE_ANNOT_SCRIPT="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/ldsc_repo/make_annot.py"
GENE_COORD_FILE="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/data/gtex/gencode.v26.GRCh38.symbols.txt"
BIM_BASE_PATH="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/GRCh38/plink_files/plink_files"


mkdir -p annot_files


for list in tissue_gene_sets/*_top10.txt; do
  
    t_name=$(basename "$list" _top10.txt)
    echo "----------------------------------------------------------------"
    echo "Processing tissue: $t_name"
    echo "----------------------------------------------------------------"
    

    for chr in {1..22}; do
        echo "Creating annotation for Chromosome: $chr"
        
        python "$MAKE_ANNOT_SCRIPT" \
            --gene-set-file "$list" \
            --gene-coord-file "$GENE_COORD_FILE" \
            --windowsize 100000 \
            --bimfile "${BIM_BASE_PATH}/1000G.EUR.hg38.${chr}.bim" \
            --annot-file annot_files/"${t_name}.${chr}.annot.gz"
    done
done

echo "Annotation generation complete for all tissues."
#!/usr/bin/env bash
SUMSTATS_FILE=$1       
OUT_DIR=$2             
REF_LD_PREFIX=$3       
WEIGHTS_PREFIX=$4      
FREQ_PREFIX=$5        
OUTPUT_NAME=$6         


if [[ -z "$SUMSTATS_FILE" || -z "$OUT_DIR" || -z "$REF_LD_PREFIX" || -z "$WEIGHTS_PREFIX" || -z "$FREQ_PREFIX" || -z "$OUTPUT_NAME" ]]; then
  echo "Usage: $0 <sumstats_file> <output_dir> <ref_ld_chr_prefix> <weights_chr_prefix> <frq_chr_prefix> <output_name>"
  exit 1
fi


mkdir -p "$OUT_DIR"
OUTPUT_PATH="${OUT_DIR}/${OUTPUT_NAME}"


echo "Running LDSC heritability estimation..."

python ldsc.py \
  --h2 "$SUMSTATS_FILE" \
  --ref-ld-chr "$REF_LD_PREFIX" \
  --w-ld-chr "$WEIGHTS_PREFIX" \
  --overlap-annot \
  --frqfile-chr "$FREQ_PREFIX" \
  --out "$OUTPUT_PATH"

# === Check status ===
if [ $? -eq 0 ]; then
  echo "LDSC finished successfully. Output: ${OUTPUT_PATH}.results"
else
  echo "LDSC failed. Check log output."
  exit 1
fi




#  python2 ldsc.py --h2-cts ../../data/heritablity/UKBB_BMI.sumstats.gz --ref-ld-chr ../../data/heritablity/1000G_Phase3_baselineLD_ldscores/baselineLD. --out ../../data/heritablity/BMI_Multi_tissue_results --ref-ld-chr-cts ../../data/heritablity/Multi_tissue_gene_expr.ldcts  --w-ld-chr ../../data/heritablity/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. 
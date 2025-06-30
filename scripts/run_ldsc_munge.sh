#!/usr/bin/env bash

# This script munges a GWAS summary stats file using LDSC's munge_sumstats.py

SUMSTATS_FILE=$1             
OUT_DIR=$2                   
MERGE_ALLELES=$3             
NCORES=${4:-1}              


if [[ -z "$SUMSTATS_FILE" || -z "$OUT_DIR" || -z "$MERGE_ALLELES" ]]; then
  echo "Usage: $0 <sumstats_file> <output_dir> <merge_alleles_file> [ncores]"
  exit 1
fi


mkdir -p "$OUT_DIR"
BASENAME=$(basename "$SUMSTATS_FILE" .sumstats.gz)
OUT_PREFIX="$OUT_DIR/munged_${BASENAME}"


echo "Running munge_sumstats.py..."
python2 munge_sumstats.py \
  --sumstats "$SUMSTATS_FILE" \
  --out "$OUT_PREFIX" \
  --merge-alleles "$MERGE_ALLELES"


if [ $? -eq 0 ]; then
  echo "✅ Munge complete. Output written to ${OUT_PREFIX}.sumstats.gz"
else
  echo "❌ Error during munging."
  exit 1
fi

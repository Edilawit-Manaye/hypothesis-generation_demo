#!/usr/bin/env bash

PLINK2_PREFIX=$1
KEEP_FILE=$2
OUTPUT_DIR=$3
CORES=${4:-8}

if [[ -z "$PLINK2_PREFIX" || -z "$KEEP_FILE" || -z "$OUTPUT_DIR" ]]; then
  echo "Usage: $0 <plink2_prefix> <keep_file> <output_dir> <num_cores>"
  echo "Example: $0 /mnt/hdd_1/rediet/hypothesis-generation-demo/data/test_plink_files/all_hg38 EUR.id /mnt/hdd_1/rediet/ldsc_plink 8"
  exit 1
fi


if [[ ! -f "${PLINK2_PREFIX}.pgen" && ! -f "${PLINK2_PREFIX}.pgen.zst" ]]; then
  echo "ERROR: PLINK2 pgen file not found: ${PLINK2_PREFIX}.pgen or ${PLINK2_PREFIX}.pgen.zst"
  exit 1
fi

if [[ ! -f "$KEEP_FILE" ]]; then
  echo "ERROR: Keep file not found: $KEEP_FILE"
  exit 1
fi

PLINK_DIR="$OUTPUT_DIR/plink_format_b38"
mkdir -p "$PLINK_DIR"

echo "========================================"
echo "Converting PLINK2 to PLINK1 format"
echo "========================================"
echo "Input: $PLINK2_PREFIX"
echo "Keep file: $KEEP_FILE"
echo "Output directory: $PLINK_DIR"
echo "Cores: $CORES"
echo ""

if [[ -f "${PLINK2_PREFIX}.pgen.zst" ]]; then
  echo "Decompressing .pgen.zst file..."
  if [[ ! -f "${PLINK2_PREFIX}.pgen" ]]; then
    plink2 --zst-decompress "${PLINK2_PREFIX}.pgen.zst" > "${PLINK2_PREFIX}.pgen"
  fi
fi

if [[ -f "${PLINK2_PREFIX}.pvar.zst" ]]; then
  echo "Decompressing .pvar.zst file..."
  if [[ ! -f "${PLINK2_PREFIX}.pvar" ]]; then
    zstd -d "${PLINK2_PREFIX}.pvar.zst" -o "${PLINK2_PREFIX}.pvar"
  fi
fi

echo ""

process_chr() {
  chr=$1
  out="$PLINK_DIR/EUR.chr${chr}.GRCh38"
  
  echo "Processing chr$chr ..."
  
  plink2 \
    --pfile "$PLINK2_PREFIX" \
    --keep "$KEEP_FILE" \
    --chr $chr \
    --make-bed \
    --out "$out" 2>&1 | sed "s/^/[chr$chr] /"
  
  if [[ $? -eq 0 && -f "${out}.bed" ]]; then
    echo "✓ chr$chr complete"
    return 0
  else
    echo "✗ chr$chr FAILED" >&2
    return 1
  fi
}

export -f process_chr
export PLINK_DIR PLINK2_PREFIX KEEP_FILE

echo "Starting parallel processing..."
echo ""

parallel -j $CORES process_chr ::: {1..22}

echo ""
echo "========================================"
echo "Processing complete!"
echo "========================================"
echo ""

echo "Summary of output files:"
total=0
success=0
total_size=0

for chr in {1..22}; do
  bed="$PLINK_DIR/EUR.chr${chr}.GRCh38.bed"
  if [[ -s "$bed" ]]; then
    size=$(du -h "$bed" | cut -f1)
    variants=$(wc -l < "$PLINK_DIR/EUR.chr${chr}.GRCh38.bim")
    samples=$(wc -l < "$PLINK_DIR/EUR.chr${chr}.GRCh38.fam")
    echo "  chr$chr: ✓ ($size, $variants variants, $samples samples)"
    ((success++))
  else
    echo "  chr$chr: ✗ MISSING or EMPTY"
  fi
  ((total++))
done

echo ""
echo "Successful: $success/$total chromosomes"
echo ""
echo "Files are ready for LDSC in: $PLINK_DIR"
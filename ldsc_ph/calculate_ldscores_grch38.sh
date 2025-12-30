#!/bin/bash

set -euo pipefail


ANNOT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch_38_annotaions"
PLINK_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/data/ldsc_plink/plink_format_b38"
LDSCORE_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch38_ld_scores_updated"
LDSC_SCRIPT="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/ldsc_repo/ldsc.py"


HM3_SNPS="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/w_hm3.snplist"
CELL_TYPES="Ast Ex In Microglia Oligo OPC PerEndo"


mkdir -p "$LDSCORE_DIR"

echo "Starting LD score calculation (single-chromosome mode)"
echo "Using global HM3 SNP list: $HM3_SNPS"
echo ""


if [[ ! -f "$HM3_SNPS" ]]; then
    echo "ERROR: HapMap3 SNP list not found:"
    echo "  $HM3_SNPS"
    exit 1
fi

for CELL in $CELL_TYPES; do
    echo "========================================"
    echo "Processing cell type: $CELL"
    echo "========================================"

    for CHR in {21..22}; do
        ANNOT_FILE="$ANNOT_DIR/${CELL}.chr${CHR}.annot.gz"
        PLINK_PREFIX="$PLINK_DIR/EUR.chr${CHR}.GRCh38"
        OUT_PREFIX="$LDSCORE_DIR/${CELL}.chr${CHR}"

        if [[ ! -f "$ANNOT_FILE" ]]; then
            echo "  [SKIP] Missing annotation: chr${CHR}"
            continue
        fi

        if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
            echo "  [SKIP] Missing PLINK files: chr${CHR}"
            continue
        fi

        echo "  → Computing LD scores for chr${CHR}"

        python "$LDSC_SCRIPT" \
            --l2 \
            --bfile "$PLINK_PREFIX" \
            --ld-wind-cm 1 \
            --annot "$ANNOT_FILE" \
            --thin-annot \
            --out "$OUT_PREFIX" \
            --print-snps /mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/grch38_ld_scores/hm3.chr${CHR}.snplist.final \
            --yes-really


        if [[ -f "${OUT_PREFIX}.l2.ldscore.gz" ]]; then
            echo "    ✓ chr${CHR} completed"
        else
            echo "    ✗ chr${CHR} FAILED"
        fi
    done

    echo ""
done

echo "========================================"
echo "SUMMARY"
echo "========================================"
ls -lh "$LDSCORE_DIR"/*.l2.ldscore.gz 2>/dev/null || echo "No LD scores generated."

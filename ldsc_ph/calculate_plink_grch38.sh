
OUTPUT_DIR=$1
VCF_DIR=$2
POP_FILE=$3       
CORES=${4:-8}    

if [[ -z "$OUTPUT_DIR" || -z "$VCF_DIR" || -z "$POP_FILE" ]]; then
  echo "Usage: $0 <output_dir> <vcf_dir> <pop_list_file> <num_cores>"
  exit 1
fi

PLINK_DIR="$OUTPUT_DIR/plink_format_b38"
mkdir -p "$PLINK_DIR"


POPS=($(cat $POP_FILE))

process_chr() {
  pop=$1
  chr=$2
  outdir="$PLINK_DIR/$pop"
  mkdir -p "$outdir"

  keep_file="$outdir/${pop}.id"
  out="$outdir/${pop}.chr${chr}.GRCh38"

  vcf="$VCF_DIR/homo_sapiens-chr${chr}.vcf.gz"

  if [[ ! -s "$vcf" ]]; then
    echo "WARN: missing VCF for chr$chr: $vcf" >&2
    return
  fi

  if [[ ! -s "$keep_file" ]]; then
    echo "ERROR: keep file not found: $keep_file" >&2
    return
  fi

  echo "Processing pop=$pop chr=$chr ..."
  plink2 \
    --vcf "$vcf" \
    --keep "$keep_file" \
    --make-bed \
    --out "$out"
}

export -f process_chr
export PLINK_DIR VCF_DIR

for chr in {1..22} X Y MT; do
  parallel -j $CORES process_chr {1} $chr ::: "${POPS[@]}"
done

echo "All chromosomes processed for all populations."

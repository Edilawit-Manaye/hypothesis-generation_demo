
SUMSTATS="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/data/gwas/21001_raw.h.tsv.gz"
LDSCORE_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/ldscores"
BASELINE_PATH="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/GRCh38/baseline_v1.2/baseline."
WEIGHTS_PATH="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/GRCh38/weights/weights_hm3_no_hla/weights."
FREQ_PATH="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/GRCh38/1000G_frq/1000G.EUR.hg38."
OUTPUT_DIR="/mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/results/gtex/final_results"

mkdir -p $OUTPUT_DIR

python /mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc/ldsc_repo/munge_sumstats.py \
    --sumstats /mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/data/gwas/21001_raw.h.tsv.gz \
    --out /mnt/hdd_1/rediet/hypothesis-generation-demo/ldsc_ph/ldsc_multi_tissue/data/gwas/21001_munged \
    --snp rsid \
    --a1 effect_allele \
    --a2 other_allele \
    --p p_value \
    --signed-sumstats beta,0 \
    --N 461460

for ld_file in ${LDSCORE_DIR}/*.22.l2.ldscore.gz; do
    t_name=$(basename "$ld_file" .22.l2.ldscore.gz)
    

    python $LDSC_EXECUTABLE \
        --h2 $SUMSTATS \
        --ref-ld-chr "${BASELINE_PATH},${LDSCORE_DIR}/${t_name}." \
        --w-ld-chr $WEIGHTS_PATH \
        --overlap-annot \
        --frqfile-chr $FREQ_PATH \
        --out "${OUTPUT_DIR}/${t_name}_enrichment" \
        --print-coefficients
done
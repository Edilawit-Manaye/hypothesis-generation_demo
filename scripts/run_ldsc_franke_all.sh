
SUMSTATS=../../data/heritablity/munged_body_BMIz.sumstats.gz
WEIGHTS=../../data/heritablity/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=../../data/heritablity/1000G_Phase3_frq/1000G.EUR.QC.
OUTDIR=../../data/heritablity/franke_ldsc_results/

for i in {1..53}
do
  PREFIX=Franke.$i
  echo "Running LDSC for $PREFIX..."

  python ldsc.py \
    --h2 $SUMSTATS \
    --ref-ld-chr ../../data/heritablity/Multi_tissue_gene_expr_1000Gv3_ldscores/${PREFIX}. \
    --w-ld-chr $WEIGHTS \
    --frqfile-chr $FRQ \
    --overlap-annot \
    --out ${OUTDIR}BMI_${PREFIX}_results
done


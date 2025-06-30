
echo -e "Tissue_ID\tEnrichment\tEnrichment_p" > gtex_bmi_enrichment_summary.tsv

for f in BMI_GTEx.*_results.results; do
    tissue=$(echo $f | sed -E 's/BMI_GTEx\.([^.]+)_results\.results/\1/')
    enrichment=$(awk 'NR > 1 {print $5}' $f)
    pval=$(awk 'NR > 1 {print $7}' $f)
    echo -e "GTEx.${tissue}\t${enrichment}\t${pval}" >> gtex_bmi_enrichment_summary.tsv
done

sort -k2,2nr gtex_bmi_enrichment_summary.tsv > sorted_gtex_bmi_enrichment.tsv

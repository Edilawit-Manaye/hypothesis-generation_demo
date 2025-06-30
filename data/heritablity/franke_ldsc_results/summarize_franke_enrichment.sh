
echo -e "Tissue_ID\tEnrichment\tEnrichment_p" > franke_bmi_enrichment_summary.tsv

for f in BMI_Franke.*_results.results; do
    tissue=$(echo $f | sed -E 's/BMI_Franke\.([^.]+)_results\.results/\1/')
    enrichment=$(awk 'NR > 1 {print $5}' $f)
    pval=$(awk 'NR > 1 {print $7}' $f)
    echo -e "Franke.${tissue}\t${enrichment}\t${pval}" >> franke_bmi_enrichment_summary.tsv
done


sort -k2,2nr franke_bmi_enrichment_summary.tsv > sorted_franke_bmi_enrichment.tsv

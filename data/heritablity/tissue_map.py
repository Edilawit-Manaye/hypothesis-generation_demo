import pandas as pd


gtex_df = pd.read_csv("gtex_ldsc_results/sorted_gtex_bmi_enrichment.tsv", sep="\t")
franke_df = pd.read_csv("franke_ldsc_results/sorted_franke_bmi_enrichment.tsv", sep="\t")

ldcts_file = "Multi_tissue_gene_expr.ldcts"
mapping = {}

with open(ldcts_file, "r") as f:
    for line in f:
        if line.strip():
            tissue_name, file_paths = line.strip().split("\t")
            prefixes = file_paths.split(",")
            for prefix in prefixes:
                prefix_clean = prefix.strip().split("/")[-1].rstrip(".")
                mapping[prefix_clean] = tissue_name


gtex_df["Tissue_Name"] = gtex_df["Tissue_ID"].map(mapping)
franke_df["Tissue_Name"] = franke_df["Tissue_ID"].map(mapping)


gtex_df.to_csv("gtex_bmi_tissue_map.tsv", sep="\t", index=False)
franke_df.to_csv("franke_bmi_tissue_map.tsv", sep="\t", index=False)

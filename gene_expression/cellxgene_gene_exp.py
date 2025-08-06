import os
import json
import cellxgene_census
import pandas as pd
import numpy as np
import pickle
import os
from scipy.stats import pearsonr

class CellxgeneMock:
    def get_coexpression_matrix(self, gene, tissue, cell_type, k=500):
        with cellxgene_census.open_soma() as census:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                # obs_value_filter=f"cell_type == '{cell_type}'",
                obs_value_filter = f"tissue == '{tissue}'",
                obs_column_names=["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]

            )

            if 'feature_id' in adata.var.columns:
                print("feature id is found inside the data")
                adata.var_names = adata.var['feature_id']
            else:
                print("Gene names column 'feature_id' not found in var DataFrame")

            gene_expression_sum = np.array((adata.X > 0).sum(axis=0)).flatten()
            adata_filtered = adata[:, gene_expression_sum > 0]
            genes = adata_filtered.var['feature_id']
            df_expression = pd.DataFrame(adata_filtered.X.toarray(), columns=genes)

            if gene in df_expression.columns:
                non_zero_samples = df_expression[df_expression[gene] > 0]
                total_samples = df_expression.shape[0]
                non_zero_sample_count = non_zero_samples.shape[0]
                non_zero_percentage = (non_zero_sample_count / total_samples) * 100

                print(f"Total samples: {total_samples}")
                print(f"Samples with non-zero expression for '{gene}': {non_zero_sample_count} ({non_zero_percentage:.2f}%)")

                correlations = {}
                for g in non_zero_samples.columns:
                    if g != gene:
                        corr, p_value = pearsonr(non_zero_samples[gene], non_zero_samples[g])
                        if p_value < 0.05:
                            correlations[g] = corr

                sorted_correlations = sorted(correlations.items(), key=lambda x: x[1], reverse=True)
                top_positive = sorted_correlations[:k]
                top_negative = sorted_correlations[-k:]
                with open(f"top_positive_{tissue_type}.txt", "w") as f: f.writelines([f"{gene}\t{corr:.4f}\n" for gene, corr in top_positive])



                return top_positive, top_negative, genes
            else:
                print(f"Gene of interest '{gene}' not found in the dataset.")
                return [], [], []
gene_of_interest = 'ENSG00000140718'
# gene_of_interest = 'ENSG00000177508'  # IRX3
cell_type = 'preadipocyte'

with open("gtex_UBERONID_cellxgene_results.json", "r") as f:
    data = json.load(f)

mock = CellxgeneMock()

for key, value in data.items():
    tissue_type = value["cellxgene_descendant_ontology_name"]
    print(f"\nProcessing tissue: {tissue_type}")
    top_positive, top_negative, all_genes = mock.get_coexpression_matrix(
        gene=gene_of_interest,
        tissue=tissue_type,
        cell_type=cell_type
    )



ensembl_to_hgnc_map = pickle.load(open("../data/ensembl_to_hgnc.pkl", "rb"))
top_positive_hgnc = [(ensembl_to_hgnc_map.get(gene, gene), corr) for gene, corr in top_positive]
top_negative_hgnc = [(ensembl_to_hgnc_map.get(gene, gene), corr) for gene, corr in top_negative]
all_genes_hgnc = [ensembl_to_hgnc_map.get(gene, gene) for gene in all_genes]
with open("top_positive_hgnc.txt", "w") as f: f.writelines([f"{gene}\t{corr:.4f}\n" for gene, corr in top_positive_hgnc])

import gseapy as gp 

library = "GO_Biological_Process_2023"
organism = "Human"

res = gp.enrichr(gene_list=[gene[0] for gene in top_positive_hgnc],
                                gene_sets=library,
                                background=all_genes_hgnc,
                                organism=organism,
                                outdir=None).results
res.drop("Gene_set", axis=1, inplace=True)
res.insert(1, "ID", res["Term"].apply(
    lambda x: x.split("(")[1].split(")")[0]))
res["Term"] = res["Term"].apply(lambda x: x.split("(")[0])
res = res[res["Adjusted P-value"] < 0.05]

print(res)


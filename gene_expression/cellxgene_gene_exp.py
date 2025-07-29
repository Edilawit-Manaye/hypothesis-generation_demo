import os
import json
import cellxgene_census
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

class CellxgeneMock:
    def __init__(self):
        self.ensembl_hgnc_map = {}

    def download_expression_matrix(self, tissue=None, cell_type=None, save_path="expression_matrix.h5ad"):
        with cellxgene_census.open_soma() as census:
            if cell_type:
                obs_filter = f"cell_type == '{cell_type}'"
            elif tissue:
                print("it is a tissue....", tissue)
                obs_filter = f"tissue == '{tissue}'"
            else:
                raise ValueError("Provide either a tissue or a cell_type.")

            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                obs_value_filter=obs_filter,
                obs_column_names=["tissue", "cell_type"]
            )
            if 'feature_id' in adata.var.columns:  
                adata.var_names = adata.var['feature_id']
            else:
                print("Gene names column 'feature_id' not found in var DataFrame")

    
            gene_expression_sum = np.array((adata.X > 0).sum(axis=0)).flatten()
            adata_filtered = adata[:, gene_expression_sum > 0]
            genes = adata_filtered.var['feature_id']
            df_expression = pd.DataFrame(adata_filtered.X.toarray(), columns=genes)
            gene="IRX3"
            if gene in df_expression.columns:
                non_zero_samples = df_expression[df_expression[gene] > 0]
            else:
                print(f"Gene of interest '{gene}' not found in the dataset.")
                return [], []
            # adata.write(save_path)
            # print(adata)
            # print(f"Saved: {save_path} ({adata.n_obs} cells × {adata.n_vars} genes)")

        return save_path


os.makedirs("results", exist_ok=True)

with open("gtex_UBERONID_cellxgene_results.json", "r") as f:
    data = json.load(f)

mock = CellxgeneMock()

for key, value in data.items():
    tissue_type = value["cellxgene_descendant_ontology_name"]
    safe_name = tissue_type.replace(" ", "_").replace("/", "_")
    filename = os.path.join("results", f"brain_{safe_name}_expr.h5ad")
    mock.download_expression_matrix(
        tissue=tissue_type,
        save_path=filename
    )

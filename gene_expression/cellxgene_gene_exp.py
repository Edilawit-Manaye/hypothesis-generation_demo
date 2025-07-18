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
                print("it is a tissue....")
                obs_filter = f"tissue == '{tissue}'"
            else:
                raise ValueError("Provide either a tissue or a cell_type.")

            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                obs_value_filter=obs_filter,
                obs_column_names=["tissue", "cell_type"]
            )

            adata.write(save_path)
            print(adata)
            print(f"Saved: {save_path} ({adata.n_obs} cells × {adata.n_vars} genes)")

        return save_path

library = "GO_Biological_Process_2023"
organism = "Human"
tissue_type = "dorsolateral prefrontal cortex"

mock = CellxgeneMock()
mock.download_expression_matrix(
    tissue=tissue_type,
    save_path=f"brain_{tissue_type.replace(' ', '_')}_expr.h5ad"
)

import cellxgene_census
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

class CellxgeneMock:
    def __init__(self):
        self.ensembl_hgnc_map = {}

    def download_expression_matrix(self, tissue=None, cell_type=None, save_path="expression_matrix.h5ad",  k=500):
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
                 
            if 'feature_id' in adata.var.columns: 
                adata.var_names = adata.var['feature_id']
            else:
                print("Gene names column 'feature_id' not found in var DataFrame")
            gene_expression_sum = np.array((adata.X > 0).sum(axis=0)).flatten()
            adata_filtered = adata[:, gene_expression_sum > 0]
            genes = adata_filtered.var['feature_id']
            
            gene="ENSG00000177508"
            df_expression = pd.DataFrame(adata_filtered.X.toarray(), columns=genes)
            if gene in df_expression.columns:
                # Filter out samples where the gene of interest is not expressed (expression value = 0)
                non_zero_samples = df_expression[df_expression[gene] > 0]
                
                # Calculate percentage of samples with non-zero expression for the gene of interest
                total_samples = df_expression.shape[0]
                non_zero_sample_count = non_zero_samples.shape[0]
                non_zero_percentage = (non_zero_sample_count / total_samples) * 100

                print(f"Total samples: {total_samples}")
                print(f"Samples with non-zero expression for '{gene}': {non_zero_sample_count} ({non_zero_percentage:.2f}%)")

                # Calculate Pearson correlation coefficients and p-values
                correlations = {}
                for g in non_zero_samples.columns:
                    if g != gene:
                        corr, p_value = pearsonr(non_zero_samples[gene], non_zero_samples[g])
                        if p_value < 0.05:
                            correlations[g] = corr
                            
                # Sort correlations
                sorted_correlations = sorted(correlations.items(), key=lambda x: x[1], reverse=True)
                top_positive = sorted_correlations[:k]
                top_negative = sorted_correlations[-k:]

                return top_positive, top_negative, genes
            else:
                print(f"Gene of interest '{gene}' not found in the dataset.")
                return [], []
        

        return save_path

library = "GO_Biological_Process_2023"
organism = "Human"
tissue_type = "putamen"

mock = CellxgeneMock()
mock.download_expression_matrix(
    tissue=tissue_type,
    save_path=f"brain_{tissue_type.replace(' ', '_')}_expr.h5ad"
)

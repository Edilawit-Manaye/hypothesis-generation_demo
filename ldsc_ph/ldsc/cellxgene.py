import cellxgene_census
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import gseapy as gp

class cellxgene_mock:

    def __init__(self):
        self.ensembl_hgnc_map = {}  

    def get_coexpression_matrix(self, gene, tissue, cell_type, k=500):
        with cellxgene_census.open_soma() as census:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                obs_value_filter=f"cell_type == '{cell_type}'",  
                column_names={"obs": ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]},
            )

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

                return top_positive, top_negative, list(genes)

            else:
                print(f"Gene of interest '{gene}' not found in the dataset.")
                return [], [], []

    def run(self):
        library = "GO_Biological_Process_2023"
        organism = "Human"
        gene_of_interest = 'ENSG00000177508'
        tissue_type = 'Adipose'
        cell_type = 'preadipocyte'

        top_positive, top_negative, all_genes = self.get_coexpression_matrix(gene_of_interest, tissue_type, cell_type, k=500)
        print("top_positive", top_positive[:10])

        top_positive_hgnc = [(self.ensembl_hgnc_map.get(gene, gene), corr) for gene, corr in top_positive]
        with open("top_positive_hgnc.tsv", "w") as f:
            f.write("Gene\tCorrelation\n")  # header
            for gene, corr in top_positive_hgnc[:10]:
              f.write(f"{gene}\t{corr}\n")

        print("top_positive_hgnc", top_positive_hgnc[:10])
        top_negative_hgnc = [(self.ensembl_hgnc_map.get(gene, gene), corr) for gene, corr in top_negative]
        all_genes_hgnc = [self.ensembl_hgnc_map.get(gene, gene) for gene in all_genes]
        

        print(all_genes, "gene_list")

        res = gp.enrichr(
            gene_list=[gene[0] for gene in top_positive_hgnc],
            gene_sets=library,
            background=all_genes_hgnc,
            organism=organism,
            outdir=None
        ).results

        res.to_csv("res_results_test_00.csv")
mock = cellxgene_mock()
mock.run()

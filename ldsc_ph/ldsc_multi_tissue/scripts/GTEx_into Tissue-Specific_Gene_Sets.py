import pandas as pd
import numpy as np
from scipy import stats


tpm = pd.read_csv('../data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2, index_col=1)
meta = pd.read_csv('../data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')


tpm.drop(columns='Name', inplace=True)

sample_to_tissue = meta.set_index('SAMPID')['SMTSD'].to_dict()
tissues = meta['SMTSD'].unique()

for tissue in tissues:
    print(f"Processing {tissue}...")
    in_cols = [c for c in tpm.columns if sample_to_tissue.get(c) == tissue]
    out_cols = [c for c in tpm.columns if sample_to_tissue.get(c) != tissue and sample_to_tissue.get(c) is not None]
    
    if len(in_cols) < 5: continue 
    
    in_data = np.log2(tpm[in_cols] + 1)
    out_data = np.log2(tpm[out_cols] + 1)
    

    t_stat, _ = stats.ttest_ind(in_data, out_data, axis=1, equal_var=False)
    
 
    tpm[f'{tissue}_tstat'] = t_stat
    top_threshold = tpm[f'{tissue}_tstat'].quantile(0.9)
    top_genes = tpm[tpm[f'{tissue}_tstat'] >= top_threshold].index
    
    with open(f"{tissue.replace(' ', '_')}_top10.txt", 'w') as f:
        for gene in top_genes:
            f.write(gene.split('.')[0] + '\n')
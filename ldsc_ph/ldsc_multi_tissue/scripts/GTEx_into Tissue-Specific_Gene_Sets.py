import pandas as pd
import numpy as np
from scipy import stats
import os


TPM_FILE = '../data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
META_FILE = '../data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
CHUNK_SIZE = 5000  
OUTPUT_DIR = '../results/gtex/tissue_gene_sets'

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print("Loading metadata...")
meta = pd.read_csv(META_FILE, sep='\t')
sample_to_tissue = meta.set_index('SAMPID')['SMTSD'].to_dict()
tissues = [t for t in meta['SMTSD'].unique() if pd.notna(t)]


tissue_stats = {t: {'n': 0, 'sum': None, 'sum_sq': None} for t in tissues}


print("Processing TPM matrix in chunks...")
reader = pd.read_csv(TPM_FILE, sep='\t', skiprows=2, chunksize=CHUNK_SIZE, index_col=1)

all_t_stats = []
gene_ids = []

for i, chunk in enumerate(reader):
    print(f"  Processing chunk {i+1}...")
    chunk = chunk.drop(columns=['Name'])
    gene_ids.extend(chunk.index.tolist())
    
 
    data = np.log2(chunk + 1)
    
    chunk_t_stats = {}
    for tissue in tissues:
        in_cols = [c for c in data.columns if sample_to_tissue.get(c) == tissue]
        out_cols = [c for c in data.columns if sample_to_tissue.get(c) != tissue]
        
        if len(in_cols) < 5: continue
        
     
        t_stat, _ = stats.ttest_ind(data[in_cols], data[out_cols], axis=1, equal_var=False)
        chunk_t_stats[tissue] = t_stat
    
    all_t_stats.append(pd.DataFrame(chunk_t_stats))

print("Ranking genes...")
full_t_df = pd.concat(all_t_stats)
full_t_df.index = gene_ids

for tissue in full_t_df.columns:
    threshold = full_t_df[tissue].quantile(0.9)
    top_genes = full_t_df[full_t_df[tissue] >= threshold].index
    
  
    clean_name = tissue.replace(' ', '_').replace('-', '_').replace('(', '').replace(')', '')
    with open(f"{OUTPUT_DIR}/{clean_name}_top10.txt", 'w') as f:
        for gene in top_genes:
            f.write(gene.split('.')[0] + '\n')

print("Done! Files are in:", OUTPUT_DIR)
import pickle

with open("gtex_tissues_to_ontology_map.pkl", "rb") as f:
    gtex_map = pickle.load(f)
import pandas as pd

df = pd.DataFrame.from_dict(gtex_map, orient='index')
df.index.name = 'GTEx_tissue'
df.reset_index(inplace=True)
df.to_csv("gtex_tissue_mappings.tsv", sep='\t', index=False)

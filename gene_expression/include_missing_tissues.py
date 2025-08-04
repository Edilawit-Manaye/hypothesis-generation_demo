import pandas as pd


df = pd.read_csv("gtex_tissue_mappings.tsv", sep="\t", header=None, names=["GTEx_tissue", "UBERON_ID"])
new_row = pd.DataFrame([["Bladder", "UBERON_0001255"], ["Blood", "UBERON_0013756"], ["Cervix_Ectocervix", "UBERON_0012249"],["Cervix_Endocervix", "UBERON_0000458"],["Fallopian_Tube","UBERON_0003889"]], columns=["GTEx_tissue", "UBERON_ID"])
df = pd.concat([df, new_row], ignore_index=True)


df.to_csv("gtex_tissue_mappings_updated.tsv", sep="\t", index=False, header=False)

import pandas as pd
from pathlib import Path


DATA_DIR = Path("data/association_by_datasource_direct")  # your parquet folder
OUTPUT_FILE = Path("kg_triples.tsv")


print("Loading Parquet files from:", DATA_DIR)
df = pd.read_parquet(DATA_DIR)
print(f"Loaded {len(df)} rows with columns: {df.columns.tolist()}")


with open(OUTPUT_FILE, "w") as out:
    for _, row in df.iterrows():
        gene = row.get('targetId')
        disease = row.get('diseaseId')
        datatype = row.get('datatypeId')

        if gene and disease:
            predicate = f"associated_with_{datatype}"
            out.write(f"{gene}\t{predicate}\t{disease}\n")

print(f"Done! Triples saved to: {OUTPUT_FILE}")

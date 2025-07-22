import cellxgene_census
import pandas as pd
import json

with cellxgene_census.open_soma() as census:
    obs_df = census["census_data"]["homo_sapiens"]["obs"]
    tissue_data = obs_df.read(column_names=["tissue"]).concat().to_pandas()
    unique_tissues = sorted(tissue_data["tissue"].dropna().unique().tolist())


with open("all_unique_tissues.json", "w") as f:
    json.dump(unique_tissues, f, indent=4)

print(f"Found {len(unique_tissues)} tissues.")

import pandas as pd
from pathlib import Path

data_dir = Path("/mnt/hdd_1/rediet/hypothesis-generation-demo/kg_rule_mining/data/association_by_datasource_direct")


df = pd.read_parquet(data_dir)
print(df.head())
print(df.columns)

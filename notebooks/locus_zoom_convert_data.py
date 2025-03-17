import pandas as pd
import json


csv_file_path = "chr16_all_region_snps.csv"  
df = pd.read_csv(csv_file_path)


json_data = {
    "data": {
        "beta": df["beta"].tolist(),
        "chromosome": df["CHR"].tolist(),
        "log_pvalue": df["log_pvalue"].tolist(),
        "position": df["POS"].tolist(),
        "ref_allele": df["minor_allele"].tolist(),
        "ref_allele_freq": df["minor_AF"].tolist(),  
        "variant": df["SNPID"].tolist()
    },
    "lastPage": None
}


json_file_path = "chr16_all_region_snps.json"
with open(json_file_path, "w") as json_file:
    json.dump(json_data, json_file, indent=4)

print(f"JSON file saved as {json_file_path}")

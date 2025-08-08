# marimo notebook
import marimo

__generated_with = "0.1.0"
app = marimo.App()

# ╔════ Cell 1: Imports
@app.cell
def __():
    import requests
    import json
    import os
    from pronto import Ontology
    import pandas as pd
    
    return requests, json, os, Ontology, pd

# ╔════ Cell 2: Utility Functions
@app.cell
def __():
    def download_file(url, filename):
        if os.path.exists(filename):
            print(f"{filename} already exists, skipping download")
            return
        print(f"Downloading {filename}...")
        response = requests.get(url)
        response.raise_for_status()
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {filename}")

    def download_json_file(url, filename):
        if os.path.exists(filename):
            print(f"{filename} already exists, loading from local file")
            with open(filename, 'r') as f:
                data = json.load(f)
            return data
        print(f"Downloading {filename}...")
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Downloaded {filename}")
        return data

    return download_file, download_json_file

# ╔════ Cell 3: Mapping and Ontology Functions
@app.cell
def __():
    def create_cellxgene_mapping_from_tissue_descendants(tissue_descendants_data):
        cellxgene_uberon_map = {}
        for parent_uberon_id in tissue_descendants_data.keys():
            cellxgene_uberon_map[parent_uberon_id] = parent_uberon_id
        for parent_uberon_id, descendants in tissue_descendants_data.items():
            if isinstance(descendants, list):
                for descendant_id in descendants:
                    cellxgene_uberon_map[descendant_id] = parent_uberon_id
        print(f"Created mapping for {len(cellxgene_uberon_map)} UBERON IDs (including descendants)")
        return cellxgene_uberon_map

    def get_tissue_name_from_ontology(uberon_id, ontology):
        print("onology name from get tissue name from ontology:", ontology)
        if not ontology:
            return None
        try:
            term = ontology[uberon_id]
            return term.name
        except KeyError:
            return None
        except Exception as e:
            print(f"Error getting tissue name for {uberon_id}: {e}")
            return None

    return create_cellxgene_mapping_from_tissue_descendants, get_tissue_name_from_ontology



if __name__ == "__main__":
    app.run()
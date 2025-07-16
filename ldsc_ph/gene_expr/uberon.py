import requests
import json
import os
from pronto import Ontology
import pandas as pd

def download_file(url, filename):
    """Download file only if it doesn't exist locally"""
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
    """Download JSON file only if it doesn't exist locally"""
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

def create_cellxgene_mapping_from_tissue_descendants(tissue_descendants_data):
    """Create mapping from tissue descendants data - includes both parent tissues and their descendants"""
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
    """Get the human-readable tissue name from UBERON ID using ontology"""
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


    

if __name__ == "__main__":
    main()
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


def main():
    uberon_url = "http://purl.obolibrary.org/obo/uberon.owl"
    uberon_filename = "uberon.owl"
    tissue_descendants_url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene-ontology-guide/latest/ontology-assets/tissue_descendants.json"
    tissue_descendants_filename = "tissue_descendants.json"

    download_file(uberon_url, uberon_filename)
    tissue_descendants_data = download_json_file(tissue_descendants_url, tissue_descendants_filename)


    gtex_uberon_mapping = {
        'Brain_Frontal_Cortex_(BA9)': 'UBERON:0009834',
        'Brain_Putamen_(basal_ganglia)':'UBERON:0001874',
        'Brain_Caudate_(basal_ganglia)':'UBERON:0001873',
        'Brain_Nucleus_accumbens_(basal_ganglia)':'UBERON:0001882',
        'Brain_Anterior_cingulate_cortex_BA24': 'UBERON:0009835',
        'Brain_Cerebellar_Hemisphere':'UBERON:0002037',
        'Brain_Cerebellum':'UBERON:0002037'


    }

    cellxgene_uberon_mapping = create_cellxgene_mapping_from_tissue_descendants(tissue_descendants_data)

    uberon_ontology = None
    if os.path.exists(uberon_filename):
        try:
            print("Loading UBERON ontology... This may take a moment.")
            uberon_ontology = Ontology(uberon_filename)
            print("UBERON ontology loaded successfully.")
        except Exception as e:
            print(f"Error loading UBERON ontology: {e}")
            print("Continuing without ontology - only direct matches will be found.")
    else:
        print(f"UBERON ontology file '{uberon_filename}' not found. Only direct matches will be found.")

    results = {}

    for gtex_tissue_to_map in [
        'Brain_Frontal_Cortex_(BA9)',
        'Brain_Putamen_(basal_ganglia)',
        'Brain_Caudate_(basal_ganglia)',
        'Brain_Nucleus_accumbens_(basal_ganglia)',
        'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Cerebellar_Hemisphere',
        'Brain_Cerebellum'
    ]:
       
        gtex_uberon_id = gtex_uberon_mapping.get(gtex_tissue_to_map)
        gtex_ontology_name = get_tissue_name_from_ontology(gtex_uberon_id, uberon_ontology) if gtex_uberon_id else None


        mapped_parent_id, match_type, notes = map_gtex_to_cellxgene_tissue(
            gtex_tissue_to_map,
            gtex_uberon_mapping,
            cellxgene_uberon_mapping,
            uberon_ontology
        )
        parent_ontology_name = get_tissue_name_from_ontology(mapped_parent_id, uberon_ontology) if mapped_parent_id else None


        if mapped_parent_id and gtex_uberon_id and mapped_parent_id != gtex_uberon_id:
            descendant_id = gtex_uberon_id
            descendant_name = gtex_ontology_name
        else:
            descendant_id = None
            descendant_name = None

        results[gtex_tissue_to_map] = {
            "gtex_tissue_name": gtex_tissue_to_map,
            "gtex_uberon_id": gtex_uberon_id,
            "gtex_ontology_name": gtex_ontology_name,
            "cellxgene_parent_uberon_id": mapped_parent_id,
            "cellxgene_parent_ontology_name": parent_ontology_name,
            "cellxgene_descendant_uberon_id": descendant_id,
            "cellxgene_descendant_ontology_name": descendant_name,
            "match_type": match_type,
            "notes": notes
        }

  
        print(f"\nResult for '{gtex_tissue_to_map}':")
        print(json.dumps(results[gtex_tissue_to_map], indent=4))

    with open("gtex_cellxgene_detailed_results.json", "w") as f:
        json.dump(results, f, indent=4)

    print("Detailed results (with parent & descendant IDs) saved to gtex_cellxgene_detailed_results.json")


    

if __name__ == "__main__":
    main()
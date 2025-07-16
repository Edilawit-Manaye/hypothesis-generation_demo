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

def map_gtex_to_cellxgene_tissue(gtex_tissue_name, gtex_uberon_map, cellxgene_uberon_map, ontology):
    """Map GTEx tissue to CellxGene tissue using UBERON ontology"""
    print(f"\n--- Mapping GTEx: '{gtex_tissue_name}' ---")

    gtex_uberon_id = gtex_uberon_map.get(gtex_tissue_name)
    if not gtex_uberon_id:
        return None, "no_direct_uberon_found", f"No direct UBERON ID found for GTEx tissue '{gtex_tissue_name}' in mapping."

    print(f"GTEx UBERON ID: {gtex_uberon_id}")
    print(f"Checking for {gtex_uberon_id} in tissue descendants...")
    

    if gtex_uberon_id in cellxgene_uberon_map:
        mapped_parent = cellxgene_uberon_map[gtex_uberon_id]
        if mapped_parent == gtex_uberon_id:
            print(f"✓ Direct match found: {gtex_uberon_id} exists as parent tissue")
            tissue_name = get_tissue_name_from_ontology(gtex_uberon_id, ontology)
            notes = f"Direct UBERON ID match found: {gtex_uberon_id}"
            if tissue_name:
                notes += f" ({tissue_name})"
            return gtex_uberon_id, "direct", notes
        else:
            print(f"✓ Descendant match found: {gtex_uberon_id} is a descendant of {mapped_parent}")
            gtex_tissue_name = get_tissue_name_from_ontology(gtex_uberon_id, ontology)
            parent_tissue_name = get_tissue_name_from_ontology(mapped_parent, ontology)
            notes = f"GTEx tissue '{gtex_uberon_id}'"
            if gtex_tissue_name:
                notes += f" ({gtex_tissue_name})"
            notes += f" is a descendant of CellxGene tissue '{mapped_parent}'"
            if parent_tissue_name:
                notes += f" ({parent_tissue_name})"
            return mapped_parent, "descendant", notes
    else:
        print(f"✗ No direct match: {gtex_uberon_id} not found in tissue descendants")


    if ontology:
        try:
            gtex_term = ontology[gtex_uberon_id]
            print(f"GTEx UBERON Term Name: {gtex_term.name}")
            gtex_ancestor_ids = set()
            try:
                for ancestor_term in gtex_term.superclasses(with_self=True):
                    gtex_ancestor_ids.add(str(ancestor_term.id))
            except AttributeError:
                gtex_ancestor_ids.add(str(gtex_term.id))
                print("Warning: Could not retrieve ancestors, using only the term itself")

            print(f"Found {len(gtex_ancestor_ids)} ancestor terms")

            best_match = None
            match_level = float('inf')
            best_notes = ""

            for cellxgene_uberon_id in cellxgene_uberon_map.keys():
                if cellxgene_uberon_id in gtex_ancestor_ids:
                    try:
                        cellxgene_term = ontology[cellxgene_uberon_id]
                        distance = 0
                        current_term = gtex_term
                        
                        if cellxgene_uberon_id == gtex_uberon_id:
                            distance = 0
                        else:
                            distance = 1
                        
                        if distance < match_level:
                            best_match = cellxgene_uberon_id
                            match_level = distance
                            best_notes = f"Broader match found: GTEx UBERON ID '{gtex_uberon_id}' ({gtex_term.name}) is related to CellxGene UBERON ID '{cellxgene_uberon_id}' ({cellxgene_term.name}). Distance: {distance}"

                    except KeyError:
                        continue

            if best_match:
                return best_match, "broader_match_found", best_notes
            else:
                return None, "no_cellxgene_match", f"No suitable CellxGene Census tissue (direct or broader) found for GTEx UBERON ID '{gtex_uberon_id}'."

        except KeyError:
            return None, "no_uberon_in_ontology", f"UBERON ID '{gtex_uberon_id}' not found in the loaded ontology."
        except Exception as e:
            return None, "ontology_error", f"Error during ontology traversal: {e}"
    else:
        return None, "ontology_not_loaded", "UBERON ontology not loaded, cannot perform hierarchical mapping."


    

if __name__ == "__main__":
    main()
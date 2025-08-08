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

# ╔════ Cell 4: GTEx to CellxGene Mapping Logic
@app.cell
def __():
    def map_gtex_to_cellxgene_tissue(gtex_tissue_name, gtex_uberon_map, cellxgene_uberon_map, ontology):
        print(f"\n--- Mapping GTEx: '{gtex_tissue_name}' ---")
        gtex_uberon_id = gtex_uberon_map.get(gtex_tissue_name)
        if not gtex_uberon_id:
            return None, "no_direct_uberon_found", f"No direct UBERON ID found for GTEx tissue '{gtex_tissue_name}' in mapping."

        print(f"GTEx UBERON ID: {gtex_uberon_id}")
        print(f"Checking for {gtex_uberon_id} in tissue descendants...")
        if gtex_uberon_id in cellxgene_uberon_map:
            mapped_parent = cellxgene_uberon_map[gtex_uberon_id]
            if mapped_parent == gtex_uberon_id:
                print(f"Direct match found: {gtex_uberon_id} exists as parent tissue")
                tissue_name = get_tissue_name_from_ontology(gtex_uberon_id, ontology)
                notes = f"Direct UBERON ID match found: {gtex_uberon_id}"
                if tissue_name:
                    notes += f" ({tissue_name})"
                return gtex_uberon_id, "direct", notes
            else:
                print(f"Descendant match found: {gtex_uberon_id} is a descendant of {mapped_parent}")
                gtex_tissue_name = get_tissue_name_from_ontology(gtex_uberon_id, ontology)
                parent_tissue_name = get_tissue_name_from_ontology(mapped_parent, ontology)
                notes = f"GTEx tissue '{gtex_uberon_id}'"
                if gtex_tissue_name:
                    notes += f" ({gtex_tissue_name})"
                notes += f" is a descendant of CellxGene tissue '{mapped_parent}'"
                if parent_tissue_name:
                    notes += f" ({parent_tissue_name})"
                return mapped_parent, "descendant", notes

        print(f"No direct match: {gtex_uberon_id} not found in tissue descendants")
        if ontology:
            try:
                gtex_term = ontology[gtex_uberon_id]
                gtex_ancestor_ids = {str(term.id) for term in gtex_term.superclasses(with_self=True)}
                best_match = None
                match_level = float('inf')
                best_notes = ""
                for cellxgene_uberon_id in cellxgene_uberon_map.keys():
                    if cellxgene_uberon_id in gtex_ancestor_ids:
                        try:
                            cellxgene_term = ontology[cellxgene_uberon_id]
                            distance = 0 if cellxgene_uberon_id == gtex_uberon_id else 1
                            if distance < match_level:
                                best_match = cellxgene_uberon_id
                                match_level = distance
                                best_notes = (
                                    f"Broader match found: GTEx UBERON ID '{gtex_uberon_id}' "
                                    f"({gtex_term.name}) is related to CellxGene UBERON ID "
                                    f"'{cellxgene_uberon_id}' ({cellxgene_term.name}). Distance: {distance}"
                                )
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

    return map_gtex_to_cellxgene_tissue

# ╔════ Cell 5: Inputs, Mapping, Results
@app.cell
def __(
    download_file, download_json_file,
    create_cellxgene_mapping_from_tissue_descendants,
    get_tissue_name_from_ontology,
    map_gtex_to_cellxgene_tissue
):
    uberon_url = "http://purl.obolibrary.org/obo/uberon.owl"
    uberon_filename = "uberon.owl"
    tissue_descendants_url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene-ontology-guide/latest/ontology-assets/tissue_descendants.json"
    tissue_descendants_filename = "tissue_descendants.json"

    download_file(uberon_url, uberon_filename)
    tissue_descendants_data = download_json_file(tissue_descendants_url, tissue_descendants_filename)

    gtex_uberon_mapping = {
        'Adipose_Subcutaneous': 'UBERON:0002190',
        'Fallopian_Tube': 'UBERON:0003889',
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

    results = {}
    for gtex_tissue_to_map in gtex_uberon_mapping:
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

    with open("all_gtex_cellxgene_detailed_results.json", "w") as f:
        json.dump(results, f, indent=4)

    print("Detailed results saved to all_gtex_cellxgene_detailed_results.json")
    return results

if __name__ == "__main__":
    app.run()
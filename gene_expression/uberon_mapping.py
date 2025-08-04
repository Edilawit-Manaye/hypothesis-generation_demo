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
    print("onology name from get tissue name from ontology:", ontology)
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
    else:
        print(f"No direct match: {gtex_uberon_id} not found in tissue descendants")


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

def main():
    uberon_url = "http://purl.obolibrary.org/obo/uberon.owl"
    uberon_filename = "uberon.owl"
    tissue_descendants_url = "https://raw.githubusercontent.com/chanzuckerberg/cellxgene-ontology-guide/latest/ontology-assets/tissue_descendants.json"
    tissue_descendants_filename = "tissue_descendants.json"

    download_file(uberon_url, uberon_filename)
    tissue_descendants_data = download_json_file(tissue_descendants_url, tissue_descendants_filename)


    gtex_uberon_mapping = {
        'Adipose_Subcutaneous': 'UBERON:0002190',
        'Adipose_Visceral_Omentum': 'UBERON:0010414',
        'Adrenal_Gland': 'UBERON:0002369',
        'Artery_Aorta': 'UBERON:0001496',
        'Artery_Coronary': 'UBERON:0001621',
        'Artery_Tibial': 'UBERON:0007610',
        'Bladder': 'UBERON:0001255',
        'Brain_Amygdala': 'UBERON:0001876',
        'Brain_Anterior_cingulate_cortex_BA24': 'UBERON:0009835',
        'Brain_Caudate_basal_ganglia': 'UBERON:0001873',
        'Brain_Cerebellar_Hemisphere': 'UBERON:0002037',
        'Brain_Cerebellum': 'UBERON:0002037',
        'Brain_Cortex': 'UBERON:0001870',
        'Brain_Frontal_Cortex_BA9': 'UBERON:0009834',
        'Brain_Hippocampus': 'UBERON:0001954',
        'Brain_Hypothalamus': 'UBERON:0001898',
        'Brain_Nucleus_accumbens_basal_ganglia': 'UBERON:0001882',
        'Brain_Putamen_basal_ganglia': 'UBERON:0001874',
        'Brain_Spinal_cord_cervical_c-1': 'UBERON:0006469',
        'Brain_Substantia_nigra': 'UBERON:0002038',
        'Breast_Mammary_Tissue': 'UBERON:0008367',
        'Cells_Cultured_fibroblasts': 'EFO:0002009',
        'Cells_EBV-transformed_lymphocytes': 'EFO:0000572',
        'Colon_Sigmoid': 'UBERON:0001159',
        'Colon_Transverse': 'UBERON:0001157',
        'Esophagus_Gastroesophageal_Junction': 'UBERON:0004550',
        'Esophagus_Mucosa': 'UBERON:0006920',
        'Esophagus_Muscularis': 'UBERON:0004648',
        'Heart_Atrial_Appendage': 'UBERON:0006631',
        'Heart_Left_Ventricle': 'UBERON:0006566',
        'Kidney_Cortex': 'UBERON:0001225',
        'Liver': 'UBERON:0001114',
        'Lung': 'UBERON:0008952',
        'Minor_Salivary_Gland': 'UBERON:0006330',
        'Muscle_Skeletal': 'UBERON:0011907',
        'Nerve_Tibial': 'UBERON:0001323',
        'Ovary': 'UBERON:0002119',
        'Pancreas': 'UBERON:0001150',
        'Pituitary': 'UBERON:0000007',
        'Prostate': 'UBERON:0002367',
        'Skin_Not_Sun_Exposed_Suprapubic': 'UBERON:0001416',
        'Skin_Sun_Exposed_Lower_leg': 'UBERON:0001511',
        'Small_Intestine_Terminal_Ileum': 'UBERON:0001211',
        'Spleen': 'UBERON:0002106',
        'Stomach': 'UBERON:0000945',
        'Testis': 'UBERON:0000473',
        'Thyroid': 'UBERON:0002046',
        'Uterus': 'UBERON:0000995',
        'Vagina': 'UBERON:0000996',
        'Whole_Blood': 'UBERON:0013756',
        'Blood': 'UBERON:0013756',
        'Cervix_Ectocervix': 'UBERON:0012249',
        'Cervix_Endocervix': 'UBERON:0000458',
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
            print("Continuing without ontology - only direct matches will be found.")
    else:
        print(f"UBERON ontology file '{uberon_filename}' not found. Only direct matches will be found.")

    results = {}

    for gtex_tissue_to_map in [

            'Adipose_Subcutaneous',
            'Adipose_Visceral_Omentum',
            'Adrenal_Gland',
            'Artery_Aorta',
            'Artery_Coronary',
            'Artery_Tibial',
            'Bladder',
            'Brain_Amygdala',
            'Brain_Anterior_cingulate_cortex_BA24',
            'Brain_Caudate_basal_ganglia',
            'Brain_Cerebellar_Hemisphere',
            'Brain_Cerebellum',
            'Brain_Cortex',
            'Brain_Frontal_Cortex_BA9',
            'Brain_Hippocampus',
            'Brain_Hypothalamus',
            'Brain_Nucleus_accumbens_basal_ganglia',
            'Brain_Putamen_basal_ganglia',
            'Brain_Spinal_cord_cervical_c-1',
            'Brain_Substantia_nigra',
            'Breast_Mammary_Tissue',
            'Cells_Cultured_fibroblasts',
            'Cells_EBV-transformed_lymphocytes',
            'Colon_Sigmoid',
            'Colon_Transverse',
            'Esophagus_Gastroesophageal_Junction',
            'Esophagus_Mucosa',
            'Esophagus_Muscularis',
            'Heart_Atrial_Appendage',
            'Heart_Left_Ventricle',
            'Kidney_Cortex',
            'Liver',
            'Lung',
            'Minor_Salivary_Gland',
            'Muscle_Skeletal',
            'Nerve_Tibial',
            'Ovary',
            'Pancreas',
            'Pituitary',
            'Prostate',
            'Skin_Not_Sun_Exposed_Suprapubic',
            'Skin_Sun_Exposed_Lower_leg',
            'Small_Intestine_Terminal_Ileum',
            'Spleen',
            'Stomach',
            'Testis',
            'Thyroid',
            'Uterus',
            'Vagina',
            'Whole_Blood',
            'Blood',
            'Cervix_Ectocervix',
            'Cervix_Endocervix',
            'Fallopian_Tube'


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

    with open("all_gtex_cellxgene_detailed_results.json", "w") as f:
        json.dump(results, f, indent=4)

    print("Detailed results (with parent & descendant IDs) saved to gtex_cellxgene_detailed_results.json")


    

if __name__ == "__main__":
    main()
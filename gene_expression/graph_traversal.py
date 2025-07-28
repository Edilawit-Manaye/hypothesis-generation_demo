import collections
from owlready2 import * 

def get_all_ancestors(uberon_id, uberon_hierarchy, visited=None):
    """
    Recursively finds all ancestors of a given UBERON ID in the hierarchy.

    Args:
        uberon_id (str): The UBERON ID to find ancestors for.
        uberon_hierarchy (dict): A dictionary representing UBERON parent-child relationships.
                                 Keys are UBERON IDs, values are lists of direct parent UBERON IDs.
        visited (set): Internal set to keep track of visited nodes to prevent infinite loops.

    Returns:
        set: A set of all ancestor UBERON IDs (including the ID itself).
    """
    if visited is None:
        visited = set()
    
    ancestors = {uberon_id}
    
    if uberon_id in visited:
        return ancestors
    
    visited.add(uberon_id)

    direct_parents = uberon_hierarchy.get(uberon_id, [])
    for parent_id in direct_parents:
        ancestors.update(get_all_ancestors(parent_id, uberon_hierarchy, visited))
    
    return ancestors

def load_uberon_hierarchy_from_owl(owl_file_path):
    """
    Loads the UBERON ontology from an OWL file and builds a parent-child hierarchy dictionary.
    Focuses on 'is_a' relationships for direct parent mapping.

    Args:
        owl_file_path (str): The path to the uberon.owl file.

    Returns:
        dict: A dictionary where keys are UBERON IDs (child) and values are lists of
              their direct parent UBERON IDs.
    """
    print(f"Loading UBERON ontology from: {owl_file_path}...")
    try:
        onto = get_ontology(f"file://{owl_file_path}").load()
        print("UBERON ontology loaded successfully.")
    except Exception as e:
        print(f"Error loading ontology: {e}")
        print("Please ensure 'uberon.owl' is in the correct path and is a valid OWL file.")
        return {}

    uberon_hierarchy = {}
    for cls in onto.classes():
        iri = cls.iri
        uberon_id = None
        if "http://purl.obolibrary.org/obo/UBERON_" in iri:
            uberon_id = iri.split("http://purl.obolibrary.org/obo/UBERON_")[-1]
            uberon_id = f"UBERON:{uberon_id}" 
        elif "UBERON:" in iri: 
            uberon_id = iri.split("UBERON:")[-1]
            uberon_id = f"UBERON:{uberon_id}"

        if not uberon_id:
            continue 
        parents = []
        for parent_cls in cls.is_a:
            if isinstance(parent_cls, ThingClass): 
                parent_iri = parent_cls.iri
                parent_uberon_id = None
                if "http://purl.obolibrary.org/obo/UBERON_" in parent_iri:
                    parent_uberon_id = parent_iri.split("http://purl.obolibrary.org/obo/UBERON_")[-1]
                    parent_uberon_id = f"UBERON:{parent_uberon_id}"
                elif "UBERON:" in parent_iri:
                    parent_uberon_id = parent_iri.split("UBERON:")[-1]
                    parent_uberon_id = f"UBERON:{parent_uberon_id}"
                
                if parent_uberon_id:
                    parents.append(parent_uberon_id)
        
        if parents:
            uberon_hierarchy[uberon_id] = parents
    
    print(f"Built hierarchy with {len(uberon_hierarchy)} entries.")
    return uberon_hierarchy

def map_gtex_to_cellxgene(gtex_tissues, cellxgene_tissues, uberon_hierarchy):
    """
    Maps GTEx tissues to CellxGene Census tissues using UBERON IDs and hierarchy.

    Args:
        gtex_tissues (list): A list of dictionaries, each representing a GTEx tissue.
                             Expected keys: 'gtex_tissue_name', 'gtex_uberon_id', 'gtex_ontology_name'.
        cellxgene_tissues (list): A list of dictionaries, each representing a CellxGene tissue.
                                  Expected keys: 'cellxgene_tissue_name', 'cellxgene_uberon_id',
                                  'cellxgene_ontology_name', 'cellxgene_general_uberon_id',
                                  'cellxgene_general_ontology_name'.
        uberon_hierarchy (dict): A dictionary representing UBERON parent-child relationships.
                                 Keys are UBERON IDs, values are lists of direct parent UBERON IDs.

    Returns:
        dict: A dictionary where keys are GTEx tissue names and values are dictionaries
              containing mapping results.
    """
    results = collections.OrderedDict()

    cellxgene_lookup = {}
    for cx_tissue in cellxgene_tissues:
        cx_uberon_id = cx_tissue['cellxgene_uberon_id']
        cx_general_uberon_id = cx_tissue['cellxgene_general_uberon_id']

        if cx_uberon_id:
            if cx_uberon_id not in cellxgene_lookup:
                cellxgene_lookup[cx_uberon_id] = []
            cellxgene_lookup[cx_uberon_id].append(cx_tissue)
        
        if cx_general_uberon_id and cx_general_uberon_id != cx_uberon_id: 
            if cx_general_uberon_id not in cellxgene_lookup:
                cellxgene_lookup[cx_general_uberon_id] = []
            cellxgene_lookup[cx_general_uberon_id].append(cx_tissue)

    for gtex_tissue in gtex_tissues:
        gtex_name = gtex_tissue['gtex_tissue_name']
        gtex_uberon = gtex_tissue['gtex_uberon_id']
        gtex_ontology_name = gtex_tissue['gtex_ontology_name']

        match_found = False
        match_details = {
            "gtex_tissue_name": gtex_name,
            "gtex_uberon_id": gtex_uberon,
            "gtex_ontology_name": gtex_ontology_name,
            "cellxgene_parent_uberon_id": None,
            "cellxgene_parent_ontology_name": None,
            "cellxgene_descendant_uberon_id": None,
            "cellxgene_descendant_ontology_name": None,
            "match_type": "no_cellxgene_match",
            "notes": f"No suitable CellxGene Census tissue (direct or broader) found for GTEx UBERON ID '{gtex_uberon}'."
        }

        if not gtex_uberon:
            match_details["notes"] = "GTEx UBERON ID is missing."
            results[gtex_name] = match_details
            continue
        gtex_uberon_ancestors = get_all_ancestors(gtex_uberon, uberon_hierarchy)
        
        for ancestor_id in gtex_uberon_ancestors:
            if ancestor_id in cellxgene_lookup:
                matched_cx_tissues = cellxgene_lookup[ancestor_id]
            
                for matched_cx_tissue in matched_cx_tissues:
                    cx_uberon = matched_cx_tissue['cellxgene_uberon_id']
                    cx_general_uberon = matched_cx_tissue['cellxgene_general_uberon_id']
                    cx_name = matched_cx_tissue['cellxgene_ontology_name']
                    cx_general_name = matched_cx_tissue['cellxgene_general_ontology_name']

                    if gtex_uberon == cx_uberon:
                        match_details.update({
                            "cellxgene_parent_uberon_id": None,
                            "cellxgene_parent_ontology_name": None,
                            "cellxgene_descendant_uberon_id": cx_uberon,
                            "cellxgene_descendant_ontology_name": cx_name,
                            "match_type": "direct",
                            "notes": f"GTEx tissue '{gtex_uberon}' ({gtex_ontology_name}) directly matches CellxGene tissue '{cx_uberon}' ({cx_name})."
                        })
                        match_found = True
                        break 
                    elif gtex_uberon == cx_general_uberon:
                        match_details.update({
                            "cellxgene_parent_uberon_id": None,
                            "cellxgene_parent_ontology_name": None,
                            "cellxgene_descendant_uberon_id": cx_general_uberon,
                            "cellxgene_descendant_ontology_name": cx_general_name,
                            "match_type": "direct_general",
                            "notes": f"GTEx tissue '{gtex_uberon}' ({gtex_ontology_name}) directly matches CellxGene general tissue '{cx_general_uberon}' ({cx_general_name})."
                        })
                        match_found = True
                        break
                    elif ancestor_id == cx_uberon:
                        match_details.update({
                            "cellxgene_parent_uberon_id": cx_uberon,
                            "cellxgene_parent_ontology_name": cx_name,
                            "cellxgene_descendant_uberon_id": gtex_uberon,
                            "cellxgene_descendant_ontology_name": gtex_ontology_name,
                            "match_type": "descendant",
                            "notes": f"GTEx tissue '{gtex_uberon}' ({gtex_ontology_name}) is a descendant of CellxGene tissue '{cx_uberon}' ({cx_name})."
                        })
                        match_found = True
                        break
                    elif ancestor_id == cx_general_uberon:
                        match_details.update({
                            "cellxgene_parent_uberon_id": cx_general_uberon,
                            "cellxgene_parent_ontology_name": cx_general_name,
                            "cellxgene_descendant_uberon_id": gtex_uberon,
                            "cellxgene_descendant_ontology_name": gtex_ontology_name,
                            "match_type": "descendant_general",
                            "notes": f"GTEx tissue '{gtex_uberon}' ({gtex_ontology_name}) is a descendant of CellxGene general tissue '{cx_general_uberon}' ({cx_general_name})."
                        })
                        match_found = True
                        break
            if match_found:
                break 

        results[gtex_name] = match_details
    return results

gtex_tissues_data = [
    {"gtex_tissue_name": "Cervix_Endocervix", "gtex_uberon_id": "UBERON:0000458", "gtex_ontology_name": "endocervix"},
    {"gtex_tissue_name": "Fallopian_Tube", "gtex_uberon_id": "UBERON:0003889", "gtex_ontology_name": "fallopian tube"},
    {"gtex_tissue_name": "Heart_Atrial_Appendage", "gtex_uberon_id": "UBERON:0002078", "gtex_ontology_name": "heart atrial appendage"},
    {"gtex_tissue_name": "Lung", "gtex_uberon_id": "UBERON:0002107", "gtex_ontology_name": "lung"},
    {"gtex_tissue_name": "Brain", "gtex_uberon_id": "UBERON:0000955", "gtex_ontology_name": "brain"},
    {"gtex_tissue_name": "Liver", "gtex_uberon_id": "UBERON:0002108", "gtex_ontology_name": "liver"},
    {"gtex_tissue_name": "Skin_Sun_Exposed_Lower_leg", "gtex_uberon_id": "UBERON:0000009", "gtex_ontology_name": "skin"}, # Example of broader term
    {"gtex_tissue_name": "Adipose_Subcutaneous", "gtex_uberon_id": "UBERON:0001013", "gtex_ontology_name": "subcutaneous adipose tissue"},
    {"gtex_tissue_name": "Blood_Vessel_Aorta", "gtex_uberon_id": "UBERON:0000465", "gtex_ontology_name": "aorta"},
    {"gtex_tissue_name": "Pancreas", "gtex_uberon_id": "UBERON:0001264", "gtex_ontology_name": "pancreas"},
    {"gtex_tissue_name": "Spleen", "gtex_uberon_id": "UBERON:0002106", "gtex_ontology_name": "spleen"},
    {"gtex_tissue_name": "Kidney_Cortex", "gtex_uberon_id": "UBERON:0001225", "gtex_ontology_name": "kidney cortex"},
    {"gtex_tissue_name": "Testis", "gtex_uberon_id": "UBERON:0000002", "gtex_ontology_name": "testis"}, # Example of a term that might be too general or incorrect UBERON ID for demo
    {"gtex_tissue_name": "Unknown_Tissue", "gtex_uberon_id": "UBERON:9999999", "gtex_ontology_name": "unknown tissue"} # Example of non-existent UBERON ID
]


cellxgene_tissues_data = [
    {"cellxgene_tissue_name": "cervix", "cellxgene_uberon_id": "UBERON:0000002", "cellxgene_ontology_name": "cervix", "cellxgene_general_uberon_id": "UBERON:0000990", "cellxgene_general_ontology_name": "reproductive system"},
    {"cellxgene_tissue_name": "fallopian tube", "cellxgene_uberon_id": "UBERON:0003889", "cellxgene_ontology_name": "fallopian tube", "cellxgene_general_uberon_id": "UBERON:0000990", "cellxgene_general_ontology_name": "reproductive system"},
    {"cellxgene_tissue_name": "heart", "cellxgene_uberon_id": "UBERON:0000948", "cellxgene_ontology_name": "heart", "cellxgene_general_uberon_id": "UBERON:0001009", "cellxgene_general_ontology_name": "cardiovascular system"},
    {"cellxgene_tissue_name": "lung", "cellxgene_uberon_id": "UBERON:0002107", "cellxgene_ontology_name": "lung", "cellxgene_general_uberon_id": "UBERON:0002107", "cellxgene_general_ontology_name": "lung"}, # General same as specific
    {"cellxgene_tissue_name": "brain", "cellxgene_uberon_id": "UBERON:0000955", "cellxgene_ontology_name": "brain", "cellxgene_general_uberon_id": "UBERON:0000955", "cellxgene_general_ontology_name": "brain"},
    {"cellxgene_tissue_name": "liver", "cellxgene_uberon_id": "UBERON:0002108", "cellxgene_ontology_name": "liver", "cellxgene_general_uberon_id": "UBERON:0002108", "cellxgene_general_ontology_name": "liver"},
    {"cellxgene_tissue_name": "skin", "cellxgene_uberon_id": "UBERON:0000009", "cellxgene_ontology_name": "skin", "cellxgene_general_uberon_id": "UBERON:0000009", "cellxgene_general_ontology_name": "skin"},
    {"cellxgene_tissue_name": "adipose tissue", "cellxgene_uberon_id": "UBERON:0001013", "cellxgene_ontology_name": "adipose tissue", "cellxgene_general_uberon_id": "UBERON:0001013", "cellxgene_general_ontology_name": "adipose tissue"},
    {"cellxgene_tissue_name": "aorta", "cellxgene_uberon_id": "UBERON:0000465", "cellxgene_ontology_name": "aorta", "cellxgene_general_uberon_id": "UBERON:0001009", "cellxgene_general_ontology_name": "cardiovascular system"},
    {"cellxgene_tissue_name": "pancreas", "cellxgene_uberon_id": "UBERON:0001264", "cellxgene_ontology_name": "pancreas", "cellxgene_general_uberon_id": "UBERON:0001264", "cellxgene_general_ontology_name": "pancreas"},
    {"cellxgene_tissue_name": "spleen", "cellxgene_uberon_id": "UBERON:0002106", "cellxgene_ontology_name": "spleen", "cellxgene_general_uberon_id": "UBERON:0002106", "cellxgene_general_ontology_name": "spleen"},
    {"cellxgene_tissue_name": "kidney", "cellxgene_uberon_id": "UBERON:0002113", "cellxgene_ontology_name": "kidney", "cellxgene_general_uberon_id": "UBERON:0002113", "cellxgene_general_ontology_name": "kidney"},
    {"cellxgene_tissue_name": "reproductive system", "cellxgene_uberon_id": "UBERON:0000990", "cellxgene_ontology_name": "reproductive system", "cellxgene_general_uberon_id": "UBERON:0000990", "cellxgene_general_ontology_name": "reproductive system"},
    {"cellxgene_tissue_name": "cardiovascular system", "cellxgene_uberon_id": "UBERON:0001009", "cellxgene_ontology_name": "cardiovascular system", "cellxgene_general_uberon_id": "UBERON:0001009", "cellxgene_general_ontology_name": "cardiovascular system"}
]


if __name__ == "__main__":
    uberon_owl_file_path = "uberon.owl" 
    uberon_hierarchy = load_uberon_hierarchy_from_owl(uberon_owl_file_path)

    if not uberon_hierarchy:
        print("Could not load UBERON hierarchy. Exiting.")
    else:
        mapped_results = map_gtex_to_cellxgene(gtex_tissues_data, cellxgene_tissues_data, uberon_hierarchy)

        import json
        print("\n--- Mapping Results ---")
        print(json.dumps(mapped_results, indent=4))

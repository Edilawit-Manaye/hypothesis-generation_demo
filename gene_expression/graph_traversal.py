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

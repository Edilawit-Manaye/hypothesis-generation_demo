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

import json
from owlready2 import *

def load_uberon_hierarchy_from_owl(owl_file_path):
    print(f"Loading UBERON ontology from: {owl_file_path}")
    onto = get_ontology(f"file://{owl_file_path}").load()
    print("Ontology loaded.")

    uberon_hierarchy = {}

    for cls in onto.classes():
        iri = cls.iri
        if "UBERON_" in iri:
            uberon_id = f"UBERON:{iri.split('UBERON_')[-1]}"
        elif "UBERON:" in iri:
            uberon_id = f"UBERON:{iri.split('UBERON:')[-1]}"
        else:
            continue

        parents = []
        for parent in cls.is_a:
            if isinstance(parent, ThingClass):
                parent_iri = parent.iri
                if "UBERON_" in parent_iri:
                    parent_id = f"UBERON:{parent_iri.split('UBERON_')[-1]}"
                elif "UBERON:" in parent_iri:
                    parent_id = f"UBERON:{parent_iri.split('UBERON:')[-1]}"
                else:
                    continue
                parents.append(parent_id)

        if parents:
            uberon_hierarchy[uberon_id] = parents

    return uberon_hierarchy

def extract_unmapped_gtex_ids(json_path):
    with open(json_path) as f:
        data = json.load(f)

    unmapped = {}
    for tissue, record in data.items():
        if record["match_type"] == "no_cellxgene_match":
            unmapped[tissue] = record["gtex_uberon_id"]
    return unmapped

def load_tissue_descendants(json_path):
    with open(json_path) as f:
        return json.load(f)

def trace_and_save_parents(unmapped_ids, uberon_hierarchy, tissue_descendants, output_path, trace_log_path):
    parent_results = {}
    trace_log = {}

    def get_all_ancestors_with_trace(uberon_id, visited=None):
        if visited is None:
            visited = set()
        ancestors = [uberon_id] 
        if uberon_id in visited:
            return ancestors
        visited.add(uberon_id)

        direct_parents = uberon_hierarchy.get(uberon_id, [])
        trace_log[uberon_id] = direct_parents
        for parent_id in direct_parents:
            ancestors += [a for a in get_all_ancestors_with_trace(parent_id, visited) if a not in ancestors]
        return ancestors

    for tissue_name, uberon_id in unmapped_ids.items():
        print(f"\nTracing ancestors for: {tissue_name} ({uberon_id})")
        ancestors = get_all_ancestors_with_trace(uberon_id)
        found_in_descendants = []

        for a in ancestors:
            if a in tissue_descendants:
                print(f"Match found: ancestor {a} is in tissue_descendants")
                found_in_descendants.append(a)
            else:
                print(f"No match: ancestor {a} is NOT in tissue_descendants")


        found_flag = False
        for ancestor in ancestors:
            if ancestor in tissue_descendants:
                print(f"First match found in tissue_descendants for {tissue_name}: {ancestor}")
                found_flag = True
                break
        if not found_flag:
            print(f"No ancestor found in tissue_descendants for {tissue_name}")

        parent_results[tissue_name] = {
            "gtex_uberon_id": uberon_id,
            "ancestor_uberon_ids": ancestors,
            "ancestors_in_tissue_descendants": found_in_descendants
        }

    with open(output_path, "w") as f:
        json.dump(parent_results, f, indent=2)
    print(f"\n Saved ancestor mappings to: {output_path}")

    with open(trace_log_path, "w") as f:
        json.dump(trace_log, f, indent=2)
    print(f"Saved parent trace log to: {trace_log_path}")


gtex_json_path = "all_gtex_cellxgene_detailed_results.json"
uberon_owl_path = "uberon.owl"
tissue_descendants_path = "tissue_descendants.json"

ancestor_output_path = "unmapped_ancestors.json"
trace_log_output_path = "parent_lookup_traces.json"

unmapped_gtex_ids = extract_unmapped_gtex_ids(gtex_json_path)
uberon_hierarchy = load_uberon_hierarchy_from_owl(uberon_owl_path)
tissue_descendants = load_tissue_descendants(tissue_descendants_path)

trace_and_save_parents(
    unmapped_ids=unmapped_gtex_ids,
    uberon_hierarchy=uberon_hierarchy,
    tissue_descendants=tissue_descendants,
    output_path=ancestor_output_path,
    trace_log_path=trace_log_output_path
)

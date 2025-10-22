import re
import os

def parse_prolog_to_amie(prolog_file, output_file="output/pl_to_amie.tsv"):
    """Directly parse Prolog file to AMIE format without loading entire KB"""
    print(f"Parsing {prolog_file} to AMIE format...")
    
    os.makedirs("data", exist_ok=True)
    
    patterns = {
        'eqtl_association': re.compile(r'eqtl_association\(snp\(([^)]+)\), gene\(([^)]+)\)\)\.'),
        'maf': re.compile(r'maf\(eqtl_association\(snp\(([^)]+)\), gene\(([^)]+)\)\), ([^)]+)\)\.'),
        'slope': re.compile(r'slope\(eqtl_association\(snp\(([^)]+)\), gene\(([^)]+)\)\), ([^)]+)\)\.'),
        'p_value': re.compile(r'p_value\(eqtl_association\(snp\(([^)]+)\), gene\(([^)]+)\)\), ([^)]+)\)\.'),
        'biological_context': re.compile(r'biological_context\(eqtl_association\(snp\(([^)]+)\), gene\(([^)]+)\)\), uberon\(([^)]+)\)\)\.')
    }
    
    facts = []
    
    with open(prolog_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('%'):
                continue
                
            for pred_name, pattern in patterns.items():
                match = pattern.match(line)
                if match:
                    if pred_name == 'eqtl_association':
                        facts.append(f"{pred_name}\t{match.group(1)}\t{match.group(2)}")
                    else:
                        facts.append(f"{pred_name}\t{match.group(1)}\t{match.group(2)}\t{match.group(3)}")
                    break
    

    with open(output_file, 'w') as f:
        for fact in facts:
            f.write(fact + '\n')
    
    print(f"✓ Parsed {len(facts)} facts to {output_file}")
    
    print("\nSample facts:")
    for fact in facts[:10]:
        print(f"  {fact}")
    
    return output_file

if __name__ == "__main__":
    prolog_file = "/mnt/hdd_1/abdu/prolog_out_v4/gtex/eqtl/edges.pl"
    output_file = parse_prolog_to_amie(prolog_file)
    
    if output_file:
        print(f"\nAMIE format file ready: {output_file}")
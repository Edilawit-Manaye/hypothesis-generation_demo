import re
import sys
import os
from datetime import datetime

IN = "output/mined_rules.rules"
OUT = "output/mined_rules.lpad"
min_pca = 0.70
min_support = 100

p_pca = re.compile(r"pca\s*=\s*([0-9]*\.?[0-9]+)", re.I)
p_support = re.compile(r"(?:support|supp)\s*=\s*(\d+)", re.I)
p_rule = re.compile(r"^(.*?)\s*<=\s*(.*?)\s*(?:#.*)?$")

def create_output_directory():
    output_dir = os.path.dirname(OUT)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}")
            return False
    return True

def normalize_pred_atom(atom):
    atom = atom.strip()
    atom = re.sub(r'\s+', ' ', atom)
    atom = re.sub(r'\s*,\s*', ', ', atom)
    return atom

def validate_rule_syntax(head, body):
    if not head or not body:
        return False
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*\s*\([^)]+\)$', head):
        return False
    body_predicates = [p.strip() for p in body.split(',')]
    for pred in body_predicates:
        if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*\s*\([^)]+\)$', pred):
            return False
    return True

def safe_float_conversion(value, default=None):
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

def safe_int_conversion(value, default=0):
    try:
        return int(value)
    except (ValueError, TypeError):
        return default

try:
    if not create_output_directory():
        sys.exit(1)

    total_rules = 0
    converted_rules = 0
    skipped_rules = 0

    with open(IN, 'r') as f, open(OUT, 'w') as out:
        out.write(f"% LPAD rules converted from {IN}\n")
        out.write(f"% Conversion date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write(f"% Thresholds: min_pca={min_pca}, min_support={min_support}\n")
        out.write("% Format: head:probability :- body.\n\n")
        
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            total_rules += 1
            m = p_rule.match(line)
            if not m:
                skipped_rules += 1
                continue
            
            head = normalize_pred_atom(m.group(1))
            body = normalize_pred_atom(m.group(2))
            
            if not validate_rule_syntax(head, body):
                skipped_rules += 1
                continue
            
            pca_m = p_pca.search(line)
            sup_m = p_support.search(line)
            
            pca = safe_float_conversion(pca_m.group(1) if pca_m else None)
            sup = safe_int_conversion(sup_m.group(1) if sup_m else 0)
            
            if pca is None:
                conf = re.search(r"conf\s*=\s*([0-9]*\.?[0-9]+)", line, re.I)
                if conf:
                    pca = safe_float_conversion(conf.group(1))
            
            if pca is None or pca < min_pca or sup < min_support:
                skipped_rules += 1
                continue
            
            lpad_line = f"{head}:{pca:.4f} :- {body}.\n"
            out.write(f"% line={line_num} support={sup} pca={pca:.4f}\n")
            out.write(lpad_line)
            out.write("\n")
            converted_rules += 1
            
            if converted_rules % 50 == 0:
                print(f"Processed {converted_rules} rules...")

    print(f"Conversion completed successfully!")
    print(f"Total rules processed: {total_rules}")
    print(f"Rules converted: {converted_rules}")
    print(f"Rules skipped: {skipped_rules}")
    print(f"Output written to: {OUT}")

except FileNotFoundError:
    print(f"Error: Input file '{IN}' not found.")
    sys.exit(1)
except PermissionError:
    print(f"Error: Permission denied when accessing files.")
    sys.exit(1)
except Exception as e:
    print(f"Unexpected error during conversion: {e}")
    sys.exit(1)

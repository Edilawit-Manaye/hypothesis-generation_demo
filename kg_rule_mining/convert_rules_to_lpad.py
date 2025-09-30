
import re

IN = "mined_rules.rules"
OUT = "mined_rules.lpad"
min_pca = 0.70
min_support = 100

p_pca = re.compile(r"pca\s*=\s*([0-9]*\.?[0-9]+)", re.I)
p_support = re.compile(r"(?:support|supp)\s*=\s*(\d+)", re.I)
p_rule = re.compile(r"^(.*?)\s*<=\s*(.*?)\s*(?:#.*)?$")

def normalize_pred_atom(atom):
    # basic normalization: strip spaces; keep AMIE variable names (usually uppercase).
    return atom.strip()

with open(IN) as f, open(OUT, "w") as out:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"): continue
        m = p_rule.match(line)
        if not m:
            continue
        head = normalize_pred_atom(m.group(1))
        body = normalize_pred_atom(m.group(2))
        pca_m = p_pca.search(line)
        sup_m = p_support.search(line)
        pca = float(pca_m.group(1)) if pca_m else None
        sup = int(sup_m.group(1)) if sup_m else 0
        if pca is None:
            # fallback: try 'conf' or 'confPCA' labels, else skip
            conf = re.search(r"conf\s*=\s*([0-9]*\.?[0-9]+)", line, re.I)
            if conf:
                pca = float(conf.group(1))
        if pca is None:
            continue
        if pca < min_pca or sup < min_support:
            continue
        # produce LPAD line
        # AMIE rule head/body already looks like: predicate(X,Y)  and atoms separated by ','
        lpad_line = f"{head}:{pca} :- {body}.\n"
        out.write(f"% support={sup} pca={pca}\n")   # provenance comment
        out.write(lpad_line)
print("Wrote", OUT)

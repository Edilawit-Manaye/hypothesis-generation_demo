import os
import subprocess

def export_prolog_to_amie(prolog_file, output_file="data/kg_amie.tsv"):
    """Export Prolog facts to AMIE format"""
    print(f"Exporting {prolog_file} to AMIE format...")
    
    os.makedirs("data", exist_ok=True)
    
    export_script = f"""
:- consult('{prolog_file}').

main :-
    open('{output_file}', write, Stream),
    export_facts(Stream),
    close(Stream),
    halt.

export_facts(Stream) :-
    % Export eQTL associations
    (   eqtl_association(snp(SNP), gene(Gene)),
        format(Stream, 'eqtl_association\\t~w\\t~w\\n', [SNP, Gene]),
        fail
    ;   true
    ),
    
    % Export MAF facts
    (   maf(eqtl_association(snp(SNP), gene(Gene)), MAF),
        format(Stream, 'maf\\t~w\\t~w\\t~w\\n', [SNP, Gene, MAF]),
        fail
    ;   true
    ),
    
    % Export slope facts
    (   slope(eqtl_association(snp(SNP), gene(Gene)), Slope),
        format(Stream, 'slope\\t~w\\t~w\\t~w\\n', [SNP, Gene, Slope]),
        fail
    ;   true
    ),
    
    % Export p-value facts
    (   p_value(eqtl_association(snp(SNP), gene(Gene)), PValue),
        format(Stream, 'p_value\\t~w\\t~w\\t~w\\n', [SNP, Gene, PValue]),
        fail
    ;   true
    ),
    
    % Export biological context
    (   biological_context(eqtl_association(snp(SNP), gene(Gene)), uberon(Context)),
        format(Stream, 'biological_context\\t~w\\t~w\\t~w\\n', [SNP, Gene, Context]),
        fail
    ;   true
    ).

:- initialization(main).
"""
    
    with open("export_temp.pl", "w") as f:
        f.write(export_script)
    
    try:
        result = subprocess.run(["swipl", "-q", "-f", "export_temp.pl"], 
                              capture_output=True, text=True, timeout=60)
        
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                lines = f.readlines()
            print(f"✓ Exported {len(lines)} facts to {output_file}")
            return output_file
        else:
            print("✗ Export failed")
            return None
            
    except Exception as e:
        print(f"Error: {e}")
        return None
    finally:
        if os.path.exists("export_temp.pl"):
            os.remove("export_temp.pl")

if __name__ == "__main__":
    prolog_file = "/mnt/hdd_1/abdu/prolog_out_v4/gtex/eqtl/edges.pl"
    output_file = export_prolog_to_amie(prolog_file)
    
    if output_file:
        print(f"\nAMIE format file ready: {output_file}")
        print("\nFirst 10 lines:")
        with open(output_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 10:
                    print(f"  {line.strip()}")
% load_kg_with_rules.pl (modified)
:- style_check(-discontiguous).
:- multifile gene/1.
:- multifile snp/1.
:- multifile enhancer/1.
:- multifile associated_with/2.
:- multifile regulates/2.
:- multifile regulatory_effect/2.
:- multifile eqtl_association/2.

% File paths for knowledge graph files
user:file_search_path(gencode_gene, '/mnt/hdd_1/abdu/prolog_out_v2/gencode/gene').
user:file_search_path(dbsnp_nodes, '/mnt/hdd_1/abdu/prolog_out_v2/dbsnp/nodes.pl').
user:file_search_path(dbsnp_edges, '/mnt/hdd_1/abdu/prolog_out_v2/dbsnp/edges.pl').
user:file_search_path(enhancer_nodes, '/mnt/hdd_1/abdu/prolog_out_v2/enhancer_atlas/nodes.pl').
user:file_search_path(enhancer_edges, '/mnt/hdd_1/abdu/prolog_out_v2/enhancer_atlas/edges.pl').
user:file_search_path(dbsuper_edges, '/mnt/hdd_1/abdu/prolog_out_v2/dbsuper/edges.pl').
user:file_search_path(relevant_gene_kb, '/mnt/hdd_1/abdu/prolog_out_v2/metta_out_v5').  % for regulatory_effect / eqtl_association

% Load a list of files safely with timing
load_with_time(Files, Name) :-
    format("Loading ~w...~n", [Name]),
    maplist(consult, Files),
    format("Loaded ~w!~n", [Name]).

% Load minimal knowledge graph
load_minimal_kg :-
    load_with_time([gencode_gene('nodes.pl')], "genes"),
    load_with_time([dbsnp_nodes], "snps"),
    load_with_time([enhancer_nodes, enhancer_edges], "enhancers"),
    load_with_time([dbsuper_edges], "associations").

% Load relevant gene facts safely
load_relevant_gene_facts :-
    load_with_time([relevant_gene_kb('regulatory_effect.pl'),
                    relevant_gene_kb('eqtl_association.pl')], "relevant_gene_facts").

% Load KG + mined rules
load_kg_with_rules :-
    format("Loading knowledge graph...~n"),
    load_minimal_kg,
    load_relevant_gene_facts,
    format("Loading mined rules...~n"),
    consult('output/mined_rules.lpad'),
    format("Knowledge graph and rules loaded successfully!~n").

% Dump all triples to TSV: predicate(subject, object)
dump_triples(File) :-
    open(File, write, Stream),
    forall(
        (clause(Pred(X,Y), _)),
        format(Stream, "~w\t~w\t~w~n", [Pred, X, Y])
    ),
    close(Stream),
    format("Dumped all triples to ~w~n", [File]).

% Convenience: load KG and dump to a specific TSV
load_and_dump(File) :-
    load_kg_with_rules,
    dump_triples(File).

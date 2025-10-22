% load_small_relevant.pl
:- style_check(-discontiguous).
:- multifile regulatory_effect/2.
:- multifile eqtl_association/2.
:- multifile activity_by_contact/2.
:- multifile relevant_gene/2.

% File paths (adjust to your directory)
user:file_search_path(relevant_gene_kb, '/mnt/hdd_1/abdu/prolog_out_v2/metta_out_v5_small').

% Load a list of files safely with timing
load_with_time(Files, Name) :-
    format("Loading ~w...~n", [Name]),
    maplist(consult, Files),
    format("Loaded ~w!~n", [Name]).

% Load only the small subsets
load_relevant_gene_small :-
    load_with_time([relevant_gene_kb('regulatory_effect_small.pl'),
                    relevant_gene_kb('eqtl_association_small.pl'),
                    relevant_gene_kb('activity_by_contact_small.pl'),
                    relevant_gene_kb('relevant_gene_small.pl')],
                   "small relevant gene facts").

% Dump all triples to TSV
dump_triples(File) :-
    open(File, write, Stream),
    forall(
        (clause(Pred(X,Y), _)),
        format(Stream, "~w\t~w\t~w~n", [Pred, X, Y])
    ),
    close(Stream),
    format("Dumped all triples to ~w~n", [File]).

% Convenience
load_and_dump(File) :-
    load_relevant_gene_small,
    dump_triples(File).

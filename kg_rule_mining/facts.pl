% SNPs and genes
snp(s1).
snp(s2).
gene(gA).
gene(gB).

% Enhancers
enhancer(e1).
enhancer(e2).

% Evidence from eQTL database
eqtl_association(snp(s1), gene(gA)).
eqtl_association(snp(s2), gene(gB)).

% Activity-by-contact data
activity_by_contact(snp(s1), gene(gB)).

% Regulatory region relationships
in_regulatory_region(snp(s1), enhancer(e1)).
in_regulatory_region(snp(s2), enhancer(e2)).

% Enhancer-to-gene links
associated_with(enhancer(e1), gene(gA)).
associated_with(enhancer(e2), gene(gB)).

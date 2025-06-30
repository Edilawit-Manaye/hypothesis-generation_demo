library(echolocatoR)
FTO <- echolocatoR::finemap_loci(
    fullSS_path = "/mnt/hdd_1/miked/echoTest/fullSS.tsv.gz",
    loci = "FTO",
    LD_reference = "1KGphase3",
    superpopulation = "EUR"
)

LD_matrix <- echolocatoR::get_LD(
    locus_dir = tempdir(),
    snp_list = FTO$SNP,
    LD_reference = "1KGphase3",
    superpopulation = "EUR"
)
write.table(LD_matrix, "data/FTO_LD_matrix.tsv", sep = "\t")

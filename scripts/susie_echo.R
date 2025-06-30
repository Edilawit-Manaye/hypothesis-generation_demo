library(echolocatoR)
library(readr)



# topSNPs <- read_tsv("data/susie/gwas/formatted/recho_mapped_topSNPs_formatted.txt")

# names(topSNPs)


# fullSS_path < -read_tsv("data/susie/gwas/formatted/21001_formatted.tsv.gz")
# names(fullSS)

topSNPs <- read_tsv("data/susie/Recho/chr_16_mapped_munged.txt")

# names(topSNPs)


# topSNPs <- read_tsv("data/susie/gwas/formatted/recho_mapped_topSNPs_formatted.txt")
# names(topSNPs)
# topSNPs_clean <- topSNPs |> dplyr::select(-matches("(?i)^locus$"))

# # topSNPS <- read_tsv("data/susie/gwas/formatted/mike_topSNPs_mapped.txt")
# # names(topSNPs)
# topSNPs <- read_tsv("data/susie/Recho/newer_cojo_mapped_topSNPs_full.txt")


# topSNPs <- topSNPs |> dplyr::select(-matches("^LOCUS$"))


colnames(topSNPs)


FTO <- echolocatoR::finemap_loci(
    fullSS_path = "data/susie/mismatch_test/harmonized_munged_gwas_data.txt",
    topSNPs = topSNPs,
    LD_reference = "1KGphase3",
    superpopulation = "EUR",
    dataset_name = "21001bothSexes",
    fullSS_genome_build = "hg19",
    loci = c("FTO"),
    bp_distance = 600000,
    finemap_methods = c("ABF", "SUSIE", "FINEMAP"),
    munged = TRUE,
    verbose = TRUE
)




# SNP	CHR	POS	A1	A2	A1_INPUTTED	Freq	Effect	StdErr	P	N	INFO	MAF	tstat	leadSNP	SUSIE.CS	SUSIE.PP	Support	Consensus_SNP	mean.PP	mean.CS
# rs10521305	16	53908484	T	C	T	0.940083	0.0124529	0.00404203	0.0023	457824	1	0.059917	3.08085293775652	FALSE	3	0.991830450743225	1	FALSE	0.991830450743225	1

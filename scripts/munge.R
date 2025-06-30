library(MungeSumstats)
library(readr)
library(dplyr)


fullSS <- read_csv("data/susie/Recho/mapped_topSNPs.txt")
head(fullSS)


fullSS <- fullSS %>%
    rename(
        SNP = rsid,
        CHR = chromosome,
        BP = position,
        SE = SE,
        P = P,
        N = N
    )


stopifnot(length(names(fullSS)) == length(unique(names(fullSS))))

munge_path <- "data/susie/munged/mapped_topSNPs.txt"
dir.create(dirname(munge_path), showWarnings = FALSE, recursive = TRUE)
dir.create("data/susie/munged/logs", showWarnings = FALSE, recursive = TRUE)

formatted_sumstats <- MungeSumstats::format_sumstats(
    path = fullSS,
    ref_genome = "GRCh37",
    dbSNP = 144,
    return_data = TRUE,
    save_path = munge_path,
    bi_allelic_filter = TRUE,
    allele_flip_check = TRUE,
    log_folder = "data/susie/munged/logs"
)

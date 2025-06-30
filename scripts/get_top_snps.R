library(echolocatoR)
library(readr)
library(dplyr)
library(vautils)
library(data.table)

topSNPs <- read_tsv("data/susie/cojo/all_chr/all_chr_cojo.jma.cojo_Recho_formatted.tsv")

top_snps <- topSNPs %>%
    dplyr::rename(rsid = SNP, chromosome = CHR, position = POS)

mapped_genes <- vautils::find_nearest_gene(as.data.frame(top_snps), build = "hg19", collapse = FALSE, snp = "rsid", flanking = 1000)

mapped_genes <- mapped_genes %>%
    mutate(distance = recode(distance, "intergenic" = "0")) %>%
    mutate(distance = abs(as.numeric(distance))) %>%
    arrange(distance) %>%
    group_by(rsid) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    rename(gene_name = GENE)

final_df <- left_join(top_snps, mapped_genes, by = c("rsid", "chromosome", "position")) %>%
    rename(Effect = BETA)

fwrite(final_df, "data/susie/Recho/mapped_topSNPs_full.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

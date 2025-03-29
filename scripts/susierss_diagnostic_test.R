library(susieR)
library(curl)
library(readr)


data1 <- read_csv('./scripts/filtered_chr16_pos53828066_snps_ld_test.csv')

library(readr)


Rin <- read_delim('./scripts/test_sig_locus_mt_r2.ld', delim = "\t", col_names = FALSE)


z_scores <- data1$beta / data1$se

lambda = estimate_s_rss(z_scores, Rin, n=359983)
lambda
# [1] 6.474096e-09

condz_in = kriging_rss(z_scores, Rin, n=n)
condz_in$plot

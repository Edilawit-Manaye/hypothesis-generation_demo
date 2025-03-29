library(susieR)
library(curl)
library(readr)

data1 <- read_csv('./scripts/filtered_chr16_pos53828066_snps_ld_test.csv')

X <- as.matrix(data1[, c("ytx", "beta", "se", "tstat", "P", "POS")])
Y <- data1$beta  

n <- nrow(data1)
b <- Y  

sumstats <- univariate_regression(X, Y)

valid_indices <- sumstats$sebetahat > 0
z_scores <- sumstats$betahat[valid_indices] / sumstats$sebetahat[valid_indices]
b_filtered <- b[valid_indices]
susie_plot(z_scores, y = "z", b = b_filtered)


summary(z_scores)

Rin <- cor(X)
attr(Rin, "eigen") <- eigen(Rin, symmetric = TRUE)

susie_plot(z_scores, y = "z", b = b)

library(readr)


data("N3finemapping")

n = nrow(N3finemapping$X)
b = N3finemapping$true_coef[,1]
sumstats <- univariate_regression(N3finemapping$X, N3finemapping$Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
Rin = cor(N3finemapping$X)
attr(Rin, "eigen") = eigen(Rin, symmetric = TRUE)
susie_plot(z_scores, y = "z", b=b)

lambda = estimate_s_rss(z_scores, Rin, n=n)
lambda
condz_in = kriging_rss(z_scores, Rin, n=n)
condz_in$plot

fit <- susie_rss(z_scores, Rin, n=n, estimate_residual_variance = TRUE)
susie_plot(fit,y = "PIP", b=b)

data_file <- tempfile(fileext = ".RData")
data_url <- paste0("https://raw.githubusercontent.com/stephenslab/susieR/",
                   "master/inst/datafiles/SummaryConsistency1k.RData")
curl_download(data_url,data_file)
load(data_file)
zflip = SummaryConsistency$z
ld = SummaryConsistency$ldref
n=10000
b = numeric(length(zflip))
b[SummaryConsistency$signal_id] = zflip[SummaryConsistency$signal_id]
plot(zflip, pch = 16, col = "#767676", main = "Marginal Associations", 
     xlab="SNP", ylab = "z scores")
points(SummaryConsistency$signal_id, zflip[SummaryConsistency$signal_id], col=2, pch=16)
points(SummaryConsistency$flip_id, zflip[SummaryConsistency$flip_id], col=7, pch=16)

fit = susie_rss(zflip, ld, n=n)
susie_plot(fit, y='PIP', b=b)
points(SummaryConsistency$flip_id, fit$pip[SummaryConsistency$flip_id], col=7, pch=16)

lambda = estimate_s_rss(zflip, ld, n=n)
lambda


condz = kriging_rss(zflip, ld, n=n)
condz$plot
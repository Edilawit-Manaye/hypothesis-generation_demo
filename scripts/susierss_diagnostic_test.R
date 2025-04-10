library(susieR)
library(curl)
library(readr)

data1 <- read_csv("./scripts/filtered_chr16_pos53828066_snps_ld_test.csv")
Rin <- read_delim("data/susie/test_ld_nan/test_gwas_sig_locus_mt_r2.ld", delim = "\t", col_names = FALSE)
R <- as.matrix(read.table("data/susie/test_ld_nan/test_gwas_sig_locus_mt_r2.ld", header = FALSE))

z_scores <- data1$beta / data1$se
R <- R[1:length(z_scores), 1:length(z_scores)]

print("finished calculating z-score")

n <- 359983
lambda <- estimate_s_rss(z_scores, R, n = 359983)
print(lambda)
#  0.8391805

print("running kriging_rss")
condz_in <- kriging_rss(z_scores, R, n = n)

print("going to plot...")
condz_in$plot

expected_z <- condz_in$conditional_dist$condmean
observed_z <- z_scores

print("Expected Z values:")
print(head(expected_z))

print("Observed Z values:")
print(head(observed_z))

if (length(expected_z) > 0) {
    residuals <- abs(expected_z - observed_z)


    threshold <- quantile(residuals, 0.95, na.rm = TRUE)
    filtered_observed_z <- observed_z[residuals <= threshold]
    filtered_expected_z <- expected_z[residuals <= threshold]

    if (length(filtered_observed_z) > 0 && length(filtered_expected_z) > 0) {
        plot(filtered_expected_z, filtered_observed_z,
            main = "Filtered Plot",
            xlab = "Expected Z", ylab = "Observed Z", col = "black"
        )
        abline(0, 1, col = "black")
    } else {
        print("No data points left after filtering.")
    }
} else {
    print("Expected Z values are empty or invalid.")
}


filtered_data <- data1[residuals <= threshold, ]
write.csv(filtered_data, "scripts/filtered_chr16_snps_residuals95.txt", row.names = FALSE)

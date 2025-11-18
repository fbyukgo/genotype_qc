args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Not enough arguments provided. Usage: Rscript plot_het.R <prefix> <maf_label> <Cohort> <out_dir>")
}

prefix <- args[1]
maf <- args[2]
Cohort <- args[3]
plot_dir <- args[4]
report_dir <- args[4]

# Read the .het file
het_file <- paste0(prefix, ".het")
het <- read.table(het_file, header=TRUE)

# Compute heterozygosity rate
# Formula: (Number of Non-Missing Genotypes - Observed Homozygotes) / Number of Non-Missing Genotypes
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."

# Calculate mean and SD
het_mean <- mean(het$HET_RATE)
het_sd <- sd(het$HET_RATE)
lower_bound <- het_mean - 3 * het_sd
upper_bound <- het_mean + 3 * het_sd

# --- PLOTTING ---
png_filename <- paste0(Cohort, "_HetCheck_", maf, ".png")
png_path <- file.path(plot_dir, png_filename)

png(filename=png_path, width=1800, height=1200, res=300)

hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", 
     main=paste0("Heterozygosity Rate for: ", maf),
     col="lightblue", border="black", breaks=30, cex.lab=0.8, cex.main=1.2)

abline(v = het_mean, col="blue", lwd=3, lty=2)  # Mean
abline(v = lower_bound, col="red", lwd=3, lty=2)  # Lower threshold
abline(v = upper_bound, col="red", lwd=3, lty=2)  # Upper threshold

legend("topright", legend=c("Mean", "Outlier Thresholds"), 
       col=c("blue", "red"), lwd=3, lty=2, cex=0.5)
dev.off()

# --- REPORTING ---
# Identify outliers
het_fail <- subset(het, (het$HET_RATE < lower_bound) | (het$HET_RATE > upper_bound))
het_fail$HET_DST <- (het_fail$HET_RATE - het_mean) / het_sd

het_number <- nrow(het_fail)
print(paste("Number of heterozygosity outliers:", het_number))

# Construct filename and path
# Added .txt extension for consistency
report_filename <- paste0(het_number, "_", Cohort, "_HetCheck_fails_", maf, ".txt")
report_path <- file.path(report_dir, report_filename)

write.table(het_fail, report_path, row.names=FALSE, quote=FALSE, sep="\t")

args = commandArgs(trailingOnly=TRUE)

prefix <- args[1]
maf = args[2]
Cohort <- args[3]
plot_dir = args[4]  
report_dir = args[4]

het <- read.table(paste(prefix, ".het", sep=""), head=TRUE)

# Compute heterozygosity rate
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."

# Calculate mean and SD
het_mean <- mean(het$HET_RATE)
het_sd <- sd(het$HET_RATE)
lower_bound <- het_mean - 3 * het_sd
upper_bound <- het_mean + 3 * het_sd

# High-quality PNG output
png(filename=paste0(plot_dir, Cohort,"_HetCheck_", maf, ".png"), width=1800, height=1200, res=300)

hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main=paste0("Heterozygosity Rate for:", maf), 
     col="lightblue", border="black", breaks=30, cex.lab=0.8, cex.main=1.2)

abline(v = het_mean, col="blue", lwd=3, lty=2)  # Mean
abline(v = lower_bound, col="red", lwd=3, lty=2)  # Lower threshold
abline(v = upper_bound, col="red", lwd=3, lty=2)  # Upper threshold

legend("topright", legend=c("Mean", "Outlier Thresholds"), col=c("blue", "red"), lwd=3, lty=2, cex=0.5)
dev.off()

# Identify outliers
het_fail = subset(het, (het$HET_RATE < lower_bound) | (het$HET_RATE > upper_bound))
het_fail$HET_DST = (het_fail$HET_RATE - het_mean) / het_sd

het_number <- dim(het_fail)[1]
write.table(het_fail, paste(report_dir, het_number, "_", Cohort, "_HetCheck_fails_", maf, sep=""), row.names=FALSE)


args = commandArgs(trailingOnly=TRUE)

prefix <- args[1]
maf = args[2]
Cohort <- "Test"
plot_dir = "plots/"
report_dir = "reports/"

het <- read.table(paste(prefix, ".het", sep=""), head=TRUE)
png(filename=paste0(plot_dir, Cohort,"_HetCheck_", maf, ".png"))
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()

het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));

het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
het_number <- dim(het_fail)[1]
write.table(het_fail, paste(report_dir, het_number, "_", Cohort, "_HetCheck_fails_", maf  , sep=""), row.names=FALSE)

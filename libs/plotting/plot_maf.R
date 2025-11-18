
args = commandArgs(trailingOnly=TRUE)

file <- args[1]

Cohort <- "Test"
#file <- "/ictstr01/groups/itg/teams/kim-hellmuth/users/furkan.bueyuekgoel/qtl/geno/data/qc_out/PLINK_M00982_plate1-4-updated"

maf_freq <- read.table(paste(file, ".frq", sep=""), header =TRUE, as.is=T)

plot_dir = "plots/"
report_dir = "reports/"


png(filename=paste0(plot_dir, Cohort,"_MafDist.png"))

hist(maf_freq[,5], freq=TRUE, col="blue", border ="black", main = "MAF distribution", xlab="MAF", ylab="Frequency", breaks=100)
abline(v=0.01, lwd=2, col="firebrick", lty=2)
dev.off()


args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
Cohort = args[2]

#setwd("/ictstr01/groups/itg/teams/kim-hellmuth/users/furkan.bueyuekgoel/qtl/geno/data")
imiss<-read.table(paste(prefix,".imiss",sep=""),header=TRUE,fill=TRUE)

plot_dir = args[3]
report_dir = args[3]
missigness_threshold = as.numeric(args[4])
#Samples to be removed at callrate threshold:
cr <- which(imiss$F_MISS > missigness_threshold)  #98% call rate
print(length(cr))
CR <- imiss[cr,]
crnumber <- dim(CR)[1]
write.table(CR, paste(report_dir, crnumber, "_", Cohort, "_", "CallRate_fails", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


#1
png(filename=paste0(plot_dir, Cohort,"_CallRate_ToRemove_1.png"))
hist(imiss$F_MISS, freq=TRUE, col="blue", border ="black", main = "Sample Call Rate", sub = Cohort, xlab="F_MISS", ylab="Frequency", ylim = c(0, 10), breaks=100)
abline(v=0.02, lwd=2, col="firebrick", lty=2)
dev.off()
#2
png(filename=paste0(plot_dir, Cohort,"_CallRate_ToRemove_2.png"))
plot(sort(imiss$F_MISS), pch=20,main = "Sample Call Rate", xlab=paste0(Cohort," samples"), ylab="F_MISS",)
abline(h=0.02, lwd=2, col="firebrick", lty=2)
dev.off()
#3
png(filename=paste0(plot_dir, Cohort,"_CallRate_ToRemove_3.png"))
plot(y=rnorm(nrow(imiss)), x=imiss$F_MISS, pch=20, main = "Sample Call Rate", sub = Cohort, xlab="F_MISS", ylab=paste0(Cohort," samples"))
abline(v=0.02, lwd=2, col="firebrick", lty=2)
dev.off()


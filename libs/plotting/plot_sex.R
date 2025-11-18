
args = commandArgs(trailingOnly=TRUE)

prefix_sex <- args[1]
prefix_xchr <- args[2]

Cohort <- args[3]
sexcheck<-read.table(paste(prefix_sex,".sexcheck",sep=""),header=TRUE,fill=TRUE)
xchr_imiss<-read.table(paste(prefix_xchr,".imiss",sep=""),header=TRUE,fill=TRUE)

sexcheck_imiss<-data.frame(FID=sexcheck$FID,IID=sexcheck$IID,PEDSEX=sexcheck$PEDSEX, SNPSEX=sexcheck$SNPSEX, STATUS=sexcheck$STATUS, F_inbreed=sexcheck$F, F_MISS=xchr_imiss$F_MISS)

plot_dir = args[4]
report_dir = args[4]


png(filename=paste0(plot_dir, Cohort,"_SexCheck.png"))
plot(sexcheck_imiss$F_inbreed, sexcheck_imiss$F_MISS, col="grey", main="Sex check", sub= Cohort, xlab="X chr inbreeding (homozygosity) estimate F", ylab="Proportion of missing SNPs for the X chr")
temp <- subset(sexcheck_imiss, sexcheck_imiss$PEDSEX=="1") #1=males
points(temp$F_inbreed, temp$F_MISS, col="Blue")
temp <- subset(sexcheck_imiss, sexcheck_imiss$PEDSEX=="2") #2=females
points(temp$F_inbreed, temp$F_MISS, col="Red")
temp <- subset(sexcheck_imiss, sexcheck_imiss$STATUS=="PROBLEM") #STATUS
points(temp$F_inbreed, temp$F_MISS, col="Yellow", pch=16,cex=0.8)
abline(v=0.8,col=2,lty=3)
abline(v=0.2,col=2,lty=3)
abline(h=0.02,col=2,lty=3)
abline(h=0.05,col=2,lty=3)
legend("topright", c("Male PEDSEX","Female PEDSEX", "Problem Status"), fill=c("Blue","Red", "Yellow"))
dev.off()



sc <- which(sexcheck_imiss$STATUS=="PROBLEM")
print(length(sc))
SC <- sexcheck_imiss[sc,c(2,3,4,5,6,7)]
scnumber <- dim(SC)[1]
write.table(SC, paste(report_dir, scnumber, "_", Cohort,  "_SexCheck_fails", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

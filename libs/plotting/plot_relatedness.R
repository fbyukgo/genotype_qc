
args = commandArgs(trailingOnly=TRUE)

prefix <- args[1]
Cohort <- args[2]
plot_dir = args[3]  
report_dir = plot_dir  



png(filename=paste0(plot_dir, Cohort,"_RelatedCheck", ".png"))

genome<-read.table(paste(prefix, ".genome", sep =""), sep = "\t", header=TRUE, fill=TRUE)
sorted_PI_HAT <- sort(genome$PI_HAT)
plot(sorted_PI_HAT, main= Cohort, xlab= "rank order of PI_HAT", ylab = "PI_HAT")
abline(h=0.2,col=2,lty=3)

dev.off() 



#subsetting out the duplicates (PH_HAT > 0.9)
wPiHatmore0.9 <- which(genome$PI_HAT > 0.9)
print(length(wPiHatmore0.9))
print(head(wPiHatmore0.9))
PiHatmore0.9 <- genome[wPiHatmore0.9,]
dupnumber <- dim(PiHatmore0.9)[1]

#With the duplicates excluded (PI_HAT>0.9) what is the range for the remaining comparisons?
PiHatless0.9 <- genome[-wPiHatmore0.9,]
print(max(PiHatless0.9$PI_HAT))
print(min(PiHatless0.9$PI_HAT))

#Numbering the duplicate pairs and writing the file out for excel:
PiHatmore0.9$Dup <- (1:nrow(PiHatmore0.9))
write.table(PiHatmore0.9, paste(report_dir, dupnumber, "_", Cohort,  "_PI_HATmore0.9", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


#subsetting out the duplicates (PH_HAT > 0.2)
wPiHatmore0.2 <- which(genome$PI_HAT > 0.2)
print(length(wPiHatmore0.2))
print(head(wPiHatmore0.2))
PiHatmore0.2 <- genome[wPiHatmore0.2,]
dupnumber <- dim(PiHatmore0.2)[1]

PiHatless0.2 <- genome[-wPiHatmore0.2,]
print(max(PiHatless0.2$PI_HAT))
print(min(PiHatless0.2$PI_HAT))

PiHatmore0.2$Dup <- (1:nrow(PiHatmore0.2))
write.table(PiHatmore0.2, paste(report_dir, dupnumber, "_", Cohort,  "_PI_HATmore0.2", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Not enough arguments provided. Usage: Rscript plot_sex.R <prefix_sex> <prefix_xchr> <Cohort> <out_dir>")
}

prefix_sex <- args[1]
prefix_xchr <- args[2]
Cohort <- args[3]
plot_dir <- args[4]
report_dir <- args[4]

# Read the input files
# Use paste0 for cleaner filename construction
sexcheck_file <- paste0(prefix_sex, ".sexcheck")
xchr_file <- paste0(prefix_xchr, ".imiss")

sexcheck <- read.table(sexcheck_file, header=TRUE, fill=TRUE)
xchr_imiss <- read.table(xchr_file, header=TRUE, fill=TRUE)

# Combine relevant data
# Note: This assumes rows are in the exact same order. 
# If Plink outputs matched files, this is fine.
sexcheck_imiss <- data.frame(
  FID = sexcheck$FID,
  IID = sexcheck$IID,
  PEDSEX = sexcheck$PEDSEX, 
  SNPSEX = sexcheck$SNPSEX, 
  STATUS = sexcheck$STATUS, 
  F_inbreed = sexcheck$F, 
  F_MISS = xchr_imiss$F_MISS
)

# --- PLOTTING ---
png_filename <- paste0(Cohort, "_SexCheck.png")
# Use file.path to join directory and filename safely
png_path <- file.path(plot_dir, png_filename)

png(filename = png_path)

# Base plot
plot(sexcheck_imiss$F_inbreed, sexcheck_imiss$F_MISS, col="grey", 
     main="Sex check", sub=Cohort, 
     xlab="X chr inbreeding (homozygosity) estimate F", 
     ylab="Proportion of missing SNPs for the X chr")

# Color coding
# 1 = Males (Blue)
temp_m <- subset(sexcheck_imiss, sexcheck_imiss$PEDSEX == "1")
points(temp_m$F_inbreed, temp_m$F_MISS, col="Blue")

# 2 = Females (Red)
temp_f <- subset(sexcheck_imiss, sexcheck_imiss$PEDSEX == "2")
points(temp_f$F_inbreed, temp_f$F_MISS, col="Red")

# PROBLEM Status (Yellow)
temp_p <- subset(sexcheck_imiss, sexcheck_imiss$STATUS == "PROBLEM")
points(temp_p$F_inbreed, temp_p$F_MISS, col="Yellow", pch=16, cex=0.8)

# Threshold lines
abline(v=0.8, col=2, lty=3)
abline(v=0.2, col=2, lty=3)
abline(h=0.02, col=2, lty=3)
abline(h=0.05, col=2, lty=3)

legend("topright", 
       c("Male PEDSEX", "Female PEDSEX", "Problem Status"), 
       fill=c("Blue", "Red", "Yellow"))

dev.off()


# --- WRITING REPORT ---
sc <- which(sexcheck_imiss$STATUS == "PROBLEM")
print(paste("Number of sex check failures:", length(sc)))

# Subset columns: IID(2), PEDSEX(3), SNPSEX(4), STATUS(5), F_inbreed(6), F_MISS(7)
SC <- sexcheck_imiss[sc, c(2,3,4,5,6,7)]
scnumber <- nrow(SC)

# Construct filename and path
report_filename <- paste0(scnumber, "_", Cohort, "_SexCheck_fails.txt")
report_path <- file.path(report_dir, report_filename)

write.table(SC, report_path, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

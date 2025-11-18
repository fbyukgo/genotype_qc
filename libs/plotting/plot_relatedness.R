args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Not enough arguments provided. Usage: Rscript plot_relatedness.R <prefix> <Cohort> <out_dir>")
}

prefix <- args[1]
Cohort <- args[2]
plot_dir <- args[3]
report_dir <- args[3]

# Read the .genome file
# Assuming the input file is tab-delimited based on your previous script, 
# but usually .genome is space-delimited. Plink 1.9 .genome can be variable.
# Your previous script used sep="\t", so keeping that, but typical default is whitespace.
genome_file <- paste0(prefix, ".genome")
genome <- read.table(genome_file, header=TRUE, fill=TRUE) 

# If the file failed to read correctly (e.g. if it wasn't tab sep), try default whitespace
if (ncol(genome) == 1) {
  genome <- read.table(genome_file, header=TRUE, fill=TRUE)
}


# --- PLOTTING ---
png_filename <- paste0(Cohort, "_RelatedCheck.png")
png_path <- file.path(plot_dir, png_filename)

png(filename = png_path)

sorted_PI_HAT <- sort(genome$PI_HAT)
plot(sorted_PI_HAT, main=Cohort, xlab="Rank order of PI_HAT", ylab="PI_HAT")
abline(h=0.2, col=2, lty=3) # Threshold line

dev.off()


# --- REPORT 1: DUPLICATES (PI_HAT > 0.9) ---
# Subsetting out the duplicates
wPiHatmore0.9 <- which(genome$PI_HAT > 0.9)
print(paste("Number of duplicates (>0.9):", length(wPiHatmore0.9)))

if (length(wPiHatmore0.9) > 0) {
  PiHatmore0.9 <- genome[wPiHatmore0.9,]
  dupnumber <- nrow(PiHatmore0.9)
  
  # Range for remaining
  PiHatless0.9 <- genome[-wPiHatmore0.9,]
  print(paste("Max PI_HAT (remaining):", max(PiHatless0.9$PI_HAT)))
  
  # Numbering the duplicate pairs
  PiHatmore0.9$Dup <- 1:nrow(PiHatmore0.9)
  
  # Write table
  report_filename1 <- paste0(dupnumber, "_", Cohort, "_PI_HATmore0.9.txt")
  report_path1 <- file.path(report_dir, report_filename1)
  
  write.table(PiHatmore0.9, report_path1, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
} else {
  print("No duplicates found (>0.9).")
}


# --- REPORT 2: RELATED (PI_HAT > 0.2) ---
# Subsetting out the related individuals
wPiHatmore0.2 <- which(genome$PI_HAT > 0.2)
print(paste("Number of related pairs (>0.2):", length(wPiHatmore0.2)))

if (length(wPiHatmore0.2) > 0) {
  PiHatmore0.2 <- genome[wPiHatmore0.2,]
  relnumber <- nrow(PiHatmore0.2)
  
  # Range for remaining
  PiHatless0.2 <- genome[-wPiHatmore0.2,]
  if (nrow(PiHatless0.2) > 0) {
      print(paste("Max PI_HAT (remaining):", max(PiHatless0.2$PI_HAT)))
  }
  
  # Numbering the pairs
  PiHatmore0.2$Dup <- 1:nrow(PiHatmore0.2)
  
  # Write table
  report_filename2 <- paste0(relnumber, "_", Cohort, "_PI_HATmore0.2.txt")
  report_path2 <- file.path(report_dir, report_filename2)
  
  write.table(PiHatmore0.2, report_path2, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
} else {
    print("No related pairs found (>0.2).")
}

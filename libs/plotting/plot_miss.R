args = commandArgs(trailingOnly=TRUE)

# Check if enough arguments are provided
if (length(args) < 4) {
  stop("Not enough arguments provided. Usage: Rscript plot_miss.R <prefix> <Cohort> <out_dir> <threshold>")
}

prefix <- args[1]
Cohort <- args[2]
plot_dir <- args[3]
report_dir <- args[3] # Assuming report and plot dir are the same based on your input
missigness_threshold <- as.numeric(args[4])

# Read the .imiss file
# used paste0 as it is cleaner than paste(..., sep="")
imiss_file <- paste0(prefix, ".imiss")
imiss <- read.table(imiss_file, header=TRUE, fill=TRUE)

# Samples to be removed at callrate threshold:
cr <- which(imiss$F_MISS > missigness_threshold)
print(paste("Number of samples failing call rate:", length(cr)))

CR <- imiss[cr,]
crnumber <- nrow(CR) # nrow is cleaner than dim()[1]

# --- WRITE REPORT ---
# Construct filename first
report_filename <- paste0(crnumber, "_", Cohort, "_CallRate_fails.txt")
# Construct full path using file.path()
report_path <- file.path(report_dir, report_filename)

write.table(CR, report_path, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# --- PLOT 1: Histogram ---
png_filename1 <- paste0(Cohort, "_CallRate_ToRemove_1.png")
png(filename = file.path(plot_dir, png_filename1))

hist(imiss$F_MISS, freq=TRUE, col="blue", border ="black", 
     main = "Sample Call Rate", sub = Cohort, 
     xlab="F_MISS", ylab="Frequency", ylim = c(0, 10), breaks=100)

# Changed hardcoded 0.02 to the variable missigness_threshold
abline(v=missigness_threshold, lwd=2, col="firebrick", lty=2) 
dev.off()

# --- PLOT 2: Sorted Plot ---
png_filename2 <- paste0(Cohort, "_CallRate_ToRemove_2.png")
png(filename = file.path(plot_dir, png_filename2))

plot(sort(imiss$F_MISS), pch=20, 
     main = "Sample Call Rate", 
     xlab=paste0(Cohort," samples"), ylab="F_MISS")

# Changed hardcoded 0.02 to the variable missigness_threshold
abline(h=missigness_threshold, lwd=2, col="firebrick", lty=2)
dev.off()

# --- PLOT 3: Random Scatter ---
png_filename3 <- paste0(Cohort, "_CallRate_ToRemove_3.png")
png(filename = file.path(plot_dir, png_filename3))

plot(y=rnorm(nrow(imiss)), x=imiss$F_MISS, pch=20, 
     main = "Sample Call Rate", sub = Cohort, 
     xlab="F_MISS", ylab=paste0(Cohort," samples"))

# Changed hardcoded 0.02 to the variable missigness_threshold
abline(v=missigness_threshold, lwd=2, col="firebrick", lty=2)
dev.off()

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
	  stop("Not enough arguments. Usage: Rscript plot_mds_2.R <prefix> <Cohort> <out_dir> <panel_file> <expected_ancestry>")
}


prefix <- args[1]
Cohort <- args[2]
plot_dir <- args[3]
report_dir <- args[3]
population_panel_merged <- args[4]
expected_ancestry <- args[5] # e.g., "EUR", "AFR", "ASN"

# Fix Paths
mds_file <- paste0(prefix, ".mds")
data <- read.table(mds_file, header=TRUE)

race <- read.table(file=population_panel_merged, header=TRUE)
datafile <- merge(data, race, by=c("IID","FID"))

# -------------------------------------------------------
# 🧮 DYNAMIC OUTLIER DETECTION (Mean +/- 3 SD)
# -------------------------------------------------------

# 1. Extract the Reference Data for the expected ancestry
ref_data <- subset(datafile, datafile$ANC == expected_ancestry)

if (nrow(ref_data) == 0) {
	  stop(paste("Error: No samples found in panel matching expected ancestry:", expected_ancestry))
}

# 2. Calculate Mean and SD for Component 1 (C1) and Component 2 (C2)
# Assuming C1 is col 4 and C2 is col 5 based on your original script
c1_mean <- mean(ref_data[,4])
c1_sd   <- sd(ref_data[,4])
c2_mean <- mean(ref_data[,5])
c2_sd   <- sd(ref_data[,5])

# 3. Define Boundaries (3 Standard Deviations)
c1_lower <- c1_mean - (3 * c1_sd)
c1_upper <- c1_mean + (3 * c1_sd)
c2_lower <- c2_mean - (3 * c2_sd)
c2_upper <- c2_mean + (3 * c2_sd)

print(paste("Defining", expected_ancestry, "Cluster:"))
print(paste("C1 Range:", round(c1_lower,4), "to", round(c1_upper,4)))
print(paste("C2 Range:", round(c2_lower,4), "to", round(c2_upper,4)))

# 4. Identify "OWN" outliers falling outside these boundaries
# Logic: Is C1 outside range? OR Is C2 outside range?
highlight_own <- subset(datafile, datafile$ANC == "OWN" & 
			                        (datafile[,4] < c1_lower | datafile[,4] > c1_upper | 
						                          datafile[,5] < c2_lower | datafile[,5] > c2_upper))

outlier_count <- nrow(highlight_own)
print(paste("Detected", outlier_count, "ancestry outliers based on 3-SD distance from", expected_ancestry))

# Write the outliers to file
report_filename <- paste0(outlier_count, "_", Cohort, "_ethnicity_outliers.tsv")
write.table(highlight_own, file = file.path(report_dir, report_filename), 
	                sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------
# 🎨 PLOTTING
# -------------------------------------------------------

png_filename <- paste0(Cohort, "_MDSPlot.png")
png(filename=file.path(plot_dir, png_filename), width=2000, height=2000, res=300)

# Dynamic limits to ensure all points are seen, or keep your hardcoded ones if preferred
# Using your defaults: xlim=c(-0.1,0.2), ylim=c(-0.15,0.1)
plot(NULL, xlim=c(-0.1,0.2), ylim=c(-0.15,0.1),
          xlab="MDS Component 1", ylab="MDS Component 2", main=paste(Cohort, "Ancestry Analysis"))

colors <- c("EUR"="green", "ASN"="red", "AMR"=470, "AFR"="blue", "OWN"="black")

# Plot all points by population
# Note: datafile[,14] is assumed to be the 'ANC' column based on your script
for (pop in names(colors)) {
	  subset_data <- datafile[datafile$ANC == pop,]
  pch_val <- ifelse(pop == "OWN", 3, 1)
    points(subset_data[,4], subset_data[,5], pch=pch_val, cex=0.7, col=colors[pop])
}

# Draw the "Safe Zone" box for the expected ancestry
rect(c1_lower, c2_lower, c1_upper, c2_upper, border="orange", lty=2, lwd=2)

# Label outliers
if (nrow(highlight_own) > 0) {
	  text(highlight_own[,4], highlight_own[,5], labels=highlight_own$IID,
	              pos=4, cex=0.3, col="black")
}

# Legend
legend("topright", pch=c(1,1,1,1,3), legend=names(colors),
              col=colors, bty="o", cex=1)

# Add a subtitle about the box
mtext(paste("Orange Box: +/- 3SD from", expected_ancestry, "mean"), side=3, cex=0.6)

dev.off()

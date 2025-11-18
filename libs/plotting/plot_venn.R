library(VennDiagram)



cr_file = Sys.glob("reports/*_Test_CallRate_fails")
cr_df = read.table(cr_file, sep="\t", header=TRUE)
cr_samples = cr_df$IID
print(cr_samples)

sex_file = Sys.glob("reports/*_Test_SexCheck_fails")
sex_df = read.table(sex_file, sep="\t", header=TRUE)
sex_samples = sex_df$IID
print(sex_samples)

het_common_file = Sys.glob("reports/*_Test_HetCheck_fails_mafgte0.01")
het_common_df = read.table(het_common_file, sep=" ", header=TRUE)
het_common_samples = het_common_df$IID
print(het_common_samples)

het_rare_file = Sys.glob("reports/*_Test_HetCheck_fails_mafless0.01")
het_rare_df = read.table(het_rare_file, sep=" ", header=TRUE)
het_rare_samples = het_rare_df$IID
print(het_rare_samples)


rel_file = Sys.glob("reports/*_Test_PI_HATmore0.2")
rel_df = read.table(rel_file, sep="\t", header=TRUE)
rel_samples = rel_df$IID1
print(rel_samples)

eth_file = Sys.glob("reports/ethnicity_outliers.tsv")
eth_df = read.table(eth_file, sep="\t", header=TRUE)
eth_samples = eth_df$IID
print(eth_samples)


sample_sets <- list(CallRate = cr_samples, SexDiscrepancies = sex_samples, HeterozygosityCheckMafCommon = het_common_samples,
                    HeterozygosityCheckMafRare = het_rare_samples, "DuplicatedSamples" = rel_samples, EthnicityOutliers = eth_samples)

print(sample_sets)
# Calculate the number of unique samples
unique_samples <- unique(unlist(sample_sets))
num_unique_samples <- length(unique_samples)

# Print the number of unique samples
cat("Total unique samples across all analyses:", num_unique_samples, "\n")


venn.plot <- venn.diagram(
  x = sample_sets,
  filename = "plots/SampleQCVenn.png",
  imagetype = "png",
  fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#0072B2"),  # Softer colors
  alpha = 0.4,       # More transparency for better visualization
  cex = 1.2,         # Size of numbers inside the diagram
  cat.cex = 0.8,     # Reduce category label size so they fit
  cat.dist = 0.1,    # Adjust distance of labels from the circles
  fontface = "bold",
  margin = 0.15      # Increase margin to prevent label cropping
)


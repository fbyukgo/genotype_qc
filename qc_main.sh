#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

#-------------------------------------------------------------------------------
# USER CONFIGURATION
#-------------------------------------------------------------------------------
# This script assumes it is in the root of your project directory, with
# 'aux_data/' and 'libs/' as subdirectories.
#
# /your_project_root/
# ├── qc_main.sh      <-- This script
# ├── aux_data/
# │   ├── 1Kg.bed
# │   ├── 1Kg.bim
# │   ├── 1Kg.fam
# │   ├── Populationfile-1kg.txt
# │   └── complex_regions.txt
# └── libs/
#     ├── HRC-1000G-check-bim.pl
#     └── plotting/
#         └── ... (all your R scripts)
#-------------------------------------------------------------------------------

# --- 1. Define Relative Paths ---
# Get the absolute path of the directory where this script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Define paths to auxiliary data and libraries relative to the script
AUX_DATA_DIR="${SCRIPT_DIR}/aux_data"
LIBS_DIR="${SCRIPT_DIR}/libs"
PLOTTING_LIBS_DIR="${LIBS_DIR}/plotting"

# --- 2. Define Paths to Reference Files (in aux_data) ---
# This is the only section you might need to edit if you rename files.
ONEKG_BFILE="${AUX_DATA_DIR}/1Kg"
ONEKG_POPULATION_FILE="${AUX_DATA_DIR}/Populationfile-1kg.txt"
COMPLEX_REGIONS_FILE="${AUX_DATA_DIR}/complex_regions.txt"
# Will's strand check script
STRAND_CHECK_SCRIPT="${LIBS_DIR}/HRC-1000G-check-bim.pl"

# --- 3. Define QC Parameters ---
# These are standard QC values. Edit them here if you need to.
PREFILTER_GENO=0.1     # Initial pre-filter SNP missingness
PREFILTER_MIND=0.1     # Initial pre-filter sample missingness
HWE_THRESHOLD=1e-6     # HWE p-value threshold for variant removal
LD_WINDOW=50           # LD pruning window size (in SNPs)
LD_STEP=5              # LD pruning step size
LD_R2=0.2              # LD pruning r^2 threshold

#-------------------------------------------------------------------------------
# SCRIPT ARGUMENTS
#-------------------------------------------------------------------------------

if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <bfile_in> <qc_out_dir> <report_plot_dir> <analysis_name> <maf_threshold> <missigness_threshold> <expected_anc> <genome_build>"
    echo "Example: $0 /path/to/mydata /qc/output /qc/reports 'MyCohort' 0.01 0.02 'EUR' 'b37'"
    exit 1
fi

bfile_in=$1           # Input bfile (e.g., /path/to/mydata)
qc_out_dir=$2         # Output directory for QC'd files
report_plot_dir=$3    # Directory for plots and reports
analysis_name=$4      # Name of the analysis (for plots)
maf_threshold=$5      # Final MAF threshold (e.g., 0.01)
missigness_threshold=$6 # Final *SNP* missingness threshold (e.g., 0.02)
expected_anc=$7 #expected ancestry to determine outliers based on 
genome_build=$8 #genome build, either b37 or b38
# Get just the filename (e.g., "mydata")
bfile_prefix=$(basename $bfile_in)

# Create output directories if they don't exist
mkdir -p $qc_out_dir
mkdir -p $report_plot_dir

echo "Starting QC process for: ${analysis_name}"
echo "Input file: ${bfile_in}"
echo "QC output directory: ${qc_out_dir}"
echo "Report directory: ${report_plot_dir}"
echo "Final MAF threshold: ${maf_threshold}"
echo "Final SNP missingness threshold: ${missigness_threshold}"

#-------------------------------------------------------------------------------
# Step 0. Pre-processing
#-------------------------------------------------------------------------------
# Comments mention 0.1, 0.2 were done, but 0.3 is the first plink command.
# Running 0.3 as the first step.

# 0.3. Prefilter exclude samples and SNPs with call rate >10% missing
echo "Step 0.3: Pre-filtering with --geno ${PREFILTER_GENO} and --mind ${PREFILTER_MIND}..."
plink --bfile ${bfile_in} --allow-no-sex --geno ${PREFILTER_GENO} --mind ${PREFILTER_MIND} --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss0.1_wopseudo_rsid_split

# splitting X pseudoautosomal region as XY (25)
echo "Step 0.3: Splitting pseudoautosomal region..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss0.1_wopseudo_rsid_split --split-x $genome_build no-fail --make-bed --out ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1

# convert variant IDS from rsids to a more standard format
echo "Step 0.3: Standardizing variant IDs..."
plink2 --bfile ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1 --set-all-var-ids @:#[$genome_build]\$r,\$a --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss0.1

# This is the base file for the rest of the QC
BASE_BFILE="${qc_out_dir}/${bfile_prefix}_miss0.1"


#-------------------------------------------------------------------------------
# Sample level QC
#-------------------------------------------------------------------------------

# Step 1  Individuals with outlying missing genotype
echo "Step 1: Calculating genome-wide missingness..."
plink --bfile ${BASE_BFILE} --missing --out ${qc_out_dir}/${bfile_prefix}_missing

echo "Step 1: Plotting missingness..."
Rscript ${PLOTTING_LIBS_DIR}/plot_miss.R ${qc_out_dir}/${bfile_prefix}_missing $analysis_name $report_plot_dir $missigness_threshold

# Step 2 Individuals with discordant sex information
echo "Step 2: Performing sex check..."
plink --bfile ${BASE_BFILE} --check-sex --out ${qc_out_dir}/${bfile_prefix}_sexcheck

plink --bfile ${BASE_BFILE} --chr 23 --make-bed --out ${qc_out_dir}/${bfile_prefix}_xchr
plink --bfile ${qc_out_dir}/${bfile_prefix}_xchr --missing --out ${qc_out_dir}/${bfile_prefix}_xchr-missing

echo "Step 2: Plotting sex discrepancies..."
Rscript ${PLOTTING_LIBS_DIR}/plot_sex.R ${qc_out_dir}/${bfile_prefix}_sexcheck ${qc_out_dir}/${bfile_prefix}_xchr-missing $analysis_name $report_plot_dir

# Step 3: Individuals with outlying heterozygosity rate
echo "Step 3: Calculating heterozygosity..."
# generate autosomal common plink files
plink --bfile ${BASE_BFILE} --autosome --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22

# Note: Using the *user-provided* MAF threshold here
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22 --maf ${maf_threshold} --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22 --exclude ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}.bim --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}

# exclude high ld regions for common (>= ${maf_threshold}) snps
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold} --exclude ${COMPLEX_REGIONS_FILE} --range --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning.prune.in --het --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned

echo "Step 3: Plotting heterozygosity (MAF >= ${maf_threshold})..."
Rscript ${PLOTTING_LIBS_DIR}/plot_het_2.R ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned mafgte${maf_threshold} $analysis_name $report_plot_dir

# then for MAF < ${maf_threshold} ones
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold} --exclude ${COMPLEX_REGIONS_FILE} --range --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_pruning
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_pruning.prune.in --het --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_ldpruned

echo "Step 3: Plotting heterozygosity (MAF < ${maf_threshold})..."
Rscript ${PLOTTING_LIBS_DIR}/plot_het_2.R ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_ldpruned mafless${maf_threshold} $analysis_name $report_plot_dir

# This LD-pruned file is used for relatedness and PCA
LD_PRUNED_FILE="${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned"
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning.prune.in --make-bed --out ${LD_PRUNED_FILE}


# Step 4: Duplicated or related individuals
echo "Step 4: Calculating relatedness (IBD)..."
plink --bfile ${LD_PRUNED_FILE} --genome --out ${qc_out_dir}/${bfile_prefix}_ldpruned_genome
# make the file tab delimited
awk -v OFS='\t' '{$1=$1}1' ${qc_out_dir}/${bfile_prefix}_ldpruned_genome.genome > ${qc_out_dir}/${bfile_prefix}_ldpruned_genome_tab.genome

echo "Step 4: Plotting relatedness..."
Rscript ${PLOTTING_LIBS_DIR}/plot_relatedness.R ${qc_out_dir}/${bfile_prefix}_ldpruned_genome_tab $analysis_name $report_plot_dir


# Step 5: Detecting ancestry outliers
echo "Step 5: Preparing population panel files for PCA..."

# Define output paths for generated panel files
PANEL_1KG_OUT="${AUX_DATA_DIR}/population_panel_1KG.txt"
PANEL_OUR_OUT="${qc_out_dir}/population_panel_our.txt"
PANEL_MERGED_OUT="${qc_out_dir}/population_panel_merged.txt"

# for our samples
fam_to_modify=${bfile_in}.fam
cat $fam_to_modify | awk '{print $1, $2, "OWN"}' > ${PANEL_OUR_OUT}

# merge panel files
cat ${PANEL_OUR_OUT} ${PANEL_1KG_OUT} | sed -e '1i\FID IID ANC' > ${PANEL_MERGED_OUT}

echo "Step 5: Processing 1000 Genomes reference data..."
onekg_prefix=$(basename $ONEKG_BFILE)

# QC the 1KG file
plink --bfile $ONEKG_BFILE --set-missing-var-ids @:#[$genome_build]\$1,\$2 --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid
plink --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid --geno ${PREFILTER_GENO} --mind ${PREFILTER_MIND} --allow-no-sex --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1
plink --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1 --maf ${maf_threshold} --allow-no-sex --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1_maf${maf_threshold}

# set variant id format for 1kg file
plink2 --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1_maf${maf_threshold} --set-all-var-ids @:#[$genome_build]\$r,\$a --make-bed --out ${qc_out_dir}/${onekg_prefix}_idcorrected

# Extract intersecting variants
plink --bfile ${qc_out_dir}/${onekg_prefix}_idcorrected --extract ${LD_PRUNED_FILE}.bim --recode --make-bed --out ${qc_out_dir}/${onekg_prefix}_intersect
plink --bfile ${LD_PRUNED_FILE} --extract ${qc_out_dir}/${onekg_prefix}_intersect.bim --recode --make-bed --out ${qc_out_dir}/${bfile_prefix}_ldpruned_intersect

# merge the 2 files
echo "Step 5: Merging cohort data with 1000 Genomes data..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_ldpruned_intersect -bmerge ${qc_out_dir}/${onekg_prefix}_intersect --make-bed --out ${qc_out_dir}/${bfile_prefix}_1KG_merged

# Run MDS
echo "Step 5: Running MDS analysis..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_1KG_merged --genome --out ${qc_out_dir}/${bfile_prefix}_1KG_merged_genome
plink --bfile ${qc_out_dir}/${bfile_prefix}_1KG_merged --read-genome ${qc_out_dir}/${bfile_prefix}_1KG_merged_genome.genome --cluster --mds-plot 10 --out ${qc_out_dir}/${bfile_prefix}_1KG_merged_mds

# plot the MDS file
echo "Step 5: Plotting MDS..."
Rscript ${PLOTTING_LIBS_DIR}/plot_mds_2.R ${qc_out_dir}/${bfile_prefix}_1KG_merged_mds $analysis_name $report_plot_dir ${PANEL_MERGED_OUT} ${expected_anc}


#-------------------------------------------------------------------------------
# Variant level QC
#-------------------------------------------------------------------------------
# Apply final filters to the main file (not the LD-pruned one)

# Step 6 Variants with outlying missing genotype
echo "Step 6: Filtering variants with missingness > ${missigness_threshold}..."
plink --bfile ${BASE_BFILE} --geno ${missigness_threshold} --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}

# Step 7 Variants out of Hardy-Weinberg Equilibrium
echo "Step 7: Filtering variants out of HWE (p < ${HWE_THRESHOLD})..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold} --hwe ${HWE_THRESHOLD} --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD}

# Step 8 Variants with low minor allele frequency
echo "Step 8: Filtering variants with MAF < ${maf_threshold}..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD} --maf ${maf_threshold} --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD}_mafgte${maf_threshold}


#-------------------------------------------------------------------------------
# QC Complete
#-------------------------------------------------------------------------------

FINAL_BFILE_NAME="${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD}_mafgte${maf_threshold}"
FINAL_BFILE_PATH="${qc_out_dir}/${FINAL_BFILE_NAME}"

echo -e "\n-----------------------------------------------------------------"
echo -e "QC is done."
echo -e "Final variant-cleaned plink files are: \n  ${FINAL_BFILE_PATH}.bed\n  ${FINAL_BFILE_PATH}.bim\n  ${FINAL_BFILE_PATH}.fam"
echo -e "\nUse the reports in ${report_plot_dir} to create a list of samples to remove (e.g., 'samples_to_remove.txt')."
echo -e "Then run the following command to get the final sample-cleaned files:"
echo -e "\n  plink --bfile ${FINAL_BFILE_PATH} --remove samples_to_remove.txt --make-bed --out ${FINAL_BFILE_PATH}_samplecleaned"
echo -e "-----------------------------------------------------------------\n"


echo "Deleting intermediate files to save space..."

# Find and delete all files in the output directory EXCEPT the final cleaned files
find $qc_out_dir -type f -not -name "${FINAL_BFILE_NAME}.*" -delete

echo "Done."

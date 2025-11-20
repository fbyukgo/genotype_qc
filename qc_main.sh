#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

#-------------------------------------------------------------------------------
# USER CONFIGURATION
#-------------------------------------------------------------------------------
# This script assumes it is in the root of your project directory.
#
# /your_project_root/
# ├── qc_main.sh      <-- This script
# ├── aux_data/       <-- Reference files
# └── libs/           <-- Scripts
#-------------------------------------------------------------------------------

# --- 1. Define Relative Paths ---
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
AUX_DATA_DIR="${SCRIPT_DIR}/aux_data"
LIBS_DIR="${SCRIPT_DIR}/libs"
PLOTTING_LIBS_DIR="${LIBS_DIR}/plotting"

# --- 2. Define Paths to Reference Files ---
ONEKG_BFILE="${AUX_DATA_DIR}/1Kg"
ONEKG_POPULATION_FILE="${AUX_DATA_DIR}/Populationfile-1kg.txt"
COMPLEX_REGIONS_FILE="${AUX_DATA_DIR}/complex_regions.txt"

# --- 3. Define QC Parameters ---
PREFILTER_GENO=0.1
PREFILTER_MIND=0.1
HWE_THRESHOLD=1e-6
LD_WINDOW=50
LD_STEP=5
LD_R2=0.2

#-------------------------------------------------------------------------------
# SCRIPT ARGUMENTS
#-------------------------------------------------------------------------------

if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <bfile_in> <qc_out_dir> <report_plot_dir> <analysis_name> <maf_threshold> <missigness_threshold> <expected_anc> <genome_build>"
    exit 1
fi

bfile_in=$1
qc_out_dir=$2
report_plot_dir=$3
analysis_name=$4
maf_threshold=$5
missigness_threshold=$6
expected_anc=$7
genome_build=$8
bfile_prefix=$(basename $bfile_in)

mkdir -p $qc_out_dir
mkdir -p $report_plot_dir

echo "Starting QC process for: ${analysis_name}"

#-------------------------------------------------------------------------------
# Step 0. Pre-processing
#-------------------------------------------------------------------------------

# 0.3. Prefilter
echo "Step 0.3: Pre-filtering with --geno ${PREFILTER_GENO} and --mind ${PREFILTER_MIND}..."
plink --bfile ${bfile_in} --allow-no-sex --geno ${PREFILTER_GENO} --mind ${PREFILTER_MIND} --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss0.1_wopseudo_rsid_split

# splitting X pseudoautosomal region
echo "Step 0.3: Splitting pseudoautosomal region..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss0.1_wopseudo_rsid_split --split-x $genome_build no-fail --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1

# convert variant IDS (Using AWK for strict column control)
# UPDATED: Explicitly setting ID = Chr:Pos_$5_$6 (ignoring Ref/Alt logic).
# We copy .bed/.fam and rewrite .bim to ensure exact column usage.
echo "Step 0.3: Standardizing variant IDs (awk based)..."
cp ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1.bed ${qc_out_dir}/${bfile_prefix}_miss0.1.bed
cp ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1.fam ${qc_out_dir}/${bfile_prefix}_miss0.1.fam
# $1=Chr, $4=Pos, $5=A1, $6=A2. We modify $2 (ID)
awk '{$2=$1":"$4"_"$5"_"$6; print $0}' ${qc_out_dir}/${bfile_prefix}_rsid_miss0.1.bim > ${qc_out_dir}/${bfile_prefix}_miss0.1.bim

BASE_BFILE="${qc_out_dir}/${bfile_prefix}_miss0.1"

#-------------------------------------------------------------------------------
# Sample level QC
#-------------------------------------------------------------------------------

# Step 1: Missingness
echo "Step 1: Calculating genome-wide missingness..."
plink --bfile ${BASE_BFILE} --missing --out ${qc_out_dir}/${bfile_prefix}_missing

echo "Step 1: Plotting missingness..."
Rscript ${PLOTTING_LIBS_DIR}/plot_miss.R ${qc_out_dir}/${bfile_prefix}_missing $analysis_name $report_plot_dir $missigness_threshold

# Step 2: Sex Check
echo "Step 2: Performing sex check..."
plink --bfile ${BASE_BFILE} --check-sex --out ${qc_out_dir}/${bfile_prefix}_sexcheck

# Extract Chr23 for missingness check
plink --bfile ${BASE_BFILE} --chr 23 --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_xchr
plink --bfile ${qc_out_dir}/${bfile_prefix}_xchr --missing --out ${qc_out_dir}/${bfile_prefix}_xchr-missing

echo "Step 2: Plotting sex discrepancies..."
Rscript ${PLOTTING_LIBS_DIR}/plot_sex.R ${qc_out_dir}/${bfile_prefix}_sexcheck ${qc_out_dir}/${bfile_prefix}_xchr-missing $analysis_name $report_plot_dir

# Step 3: Heterozygosity
echo "Step 3: Calculating heterozygosity..."
# generate autosomal common plink files
plink --bfile ${BASE_BFILE} --autosome --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22

# MAF Filter
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22 --maf ${maf_threshold} --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22 --exclude ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}.bim --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}

# Exclude high LD regions
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold} --exclude ${COMPLEX_REGIONS_FILE} --range --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning.prune.in --het --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned

echo "Step 3: Plotting heterozygosity (MAF >= ${maf_threshold})..."
Rscript ${PLOTTING_LIBS_DIR}/plot_het_2.R ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned mafgte${maf_threshold} $analysis_name $report_plot_dir

# Same for MAF < threshold
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold} --exclude ${COMPLEX_REGIONS_FILE} --range --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_pruning
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_pruning.prune.in --het --out ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_ldpruned

echo "Step 3: Plotting heterozygosity (MAF < ${maf_threshold})..."
Rscript ${PLOTTING_LIBS_DIR}/plot_het_2.R ${qc_out_dir}/${bfile_prefix}_chr1-22_mafless${maf_threshold}_nocr_ldpruned mafless${maf_threshold} $analysis_name $report_plot_dir

# Create LD Pruned File
LD_PRUNED_FILE="${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_ldpruned"
plink --bfile ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr --extract ${qc_out_dir}/${bfile_prefix}_chr1-22_mafgte${maf_threshold}_nocr_pruning.prune.in --keep-allele-order --make-bed --out ${LD_PRUNED_FILE}

# Step 4: Relatedness
echo "Step 4: Calculating relatedness (IBD)..."
plink --bfile ${LD_PRUNED_FILE} --genome --out ${qc_out_dir}/${bfile_prefix}_ldpruned_genome
awk -v OFS='\t' '{$1=$1}1' ${qc_out_dir}/${bfile_prefix}_ldpruned_genome.genome > ${qc_out_dir}/${bfile_prefix}_ldpruned_genome_tab.genome

echo "Step 4: Plotting relatedness..."
Rscript ${PLOTTING_LIBS_DIR}/plot_relatedness.R ${qc_out_dir}/${bfile_prefix}_ldpruned_genome_tab $analysis_name $report_plot_dir

# Step 5: Ancestry / PCA
echo "Step 5: Preparing population panel files for PCA..."

PANEL_1KG_OUT="${AUX_DATA_DIR}/population_panel_1KG.txt"
PANEL_OUR_OUT="${qc_out_dir}/population_panel_our.txt"
PANEL_MERGED_OUT="${qc_out_dir}/population_panel_merged.txt"

fam_to_modify=${bfile_in}.fam
cat $fam_to_modify | awk '{print $1, $2, "OWN"}' > ${PANEL_OUR_OUT}
cat ${PANEL_OUR_OUT} ${PANEL_1KG_OUT} | sed -e '1i\FID IID ANC' > ${PANEL_MERGED_OUT}

echo "Step 5: Processing 1000 Genomes reference data..."
onekg_prefix=$(basename $ONEKG_BFILE)

# QC 1KG
# UPDATED: Using PLINK2 to infer Ref ($r) and Alt ($a) for naming.
# CHANGED: Switched to --set-missing-var-ids so it only fills IDs that are missing (.).
# Existing IDs (e.g. rsIDs) will be preserved.
plink2 --bfile $ONEKG_BFILE --set-missing-var-ids @:#_\$r_\$a --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid

plink --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid --geno ${PREFILTER_GENO} --mind ${PREFILTER_MIND} --allow-no-sex --keep-allele-order --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1
plink --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1 --maf ${maf_threshold} --allow-no-sex --keep-allele-order --make-bed --out ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1_maf${maf_threshold}

# Extract intersecting variants
plink --bfile ${qc_out_dir}/${onekg_prefix}_nomisingid_miss0.1_mind0.1_maf${maf_threshold} --extract ${LD_PRUNED_FILE}.bim --keep-allele-order --recode --make-bed --out ${qc_out_dir}/${onekg_prefix}_intersect
plink --bfile ${LD_PRUNED_FILE} --extract ${qc_out_dir}/${onekg_prefix}_intersect.bim --keep-allele-order --recode --make-bed --out ${qc_out_dir}/${bfile_prefix}_ldpruned_intersect

# Merge
echo "Step 5: Merging cohort data with 1000 Genomes data..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_ldpruned_intersect -bmerge ${qc_out_dir}/${onekg_prefix}_intersect --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_1KG_merged

# Run MDS
echo "Step 5: Running MDS analysis..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_1KG_merged --genome --out ${qc_out_dir}/${bfile_prefix}_1KG_merged_genome
plink --bfile ${qc_out_dir}/${bfile_prefix}_1KG_merged --read-genome ${qc_out_dir}/${bfile_prefix}_1KG_merged_genome.genome --cluster --mds-plot 10 --out ${qc_out_dir}/${bfile_prefix}_1KG_merged_mds

echo "Step 5: Plotting MDS..."
Rscript ${PLOTTING_LIBS_DIR}/plot_mds_2.R ${qc_out_dir}/${bfile_prefix}_1KG_merged_mds $analysis_name $report_plot_dir ${PANEL_MERGED_OUT} ${expected_anc}

#-------------------------------------------------------------------------------
# Variant level QC
#-------------------------------------------------------------------------------

# Step 6: Geno Filter
echo "Step 6: Filtering variants with missingness > ${missigness_threshold}..."
plink --bfile ${BASE_BFILE} --geno ${missigness_threshold} --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}

# Step 7: HWE
echo "Step 7: Filtering variants out of HWE (p < ${HWE_THRESHOLD})..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold} --hwe ${HWE_THRESHOLD} --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD}

# Step 8: MAF
echo "Step 8: Filtering variants with MAF < ${maf_threshold}..."
plink --bfile ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD} --maf ${maf_threshold} --keep-allele-order --make-bed --out ${qc_out_dir}/${bfile_prefix}_miss${missigness_threshold}_hwe${HWE_THRESHOLD}_mafgte${maf_threshold}

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
echo -e "\n  plink --bfile ${FINAL_BFILE_PATH} --remove samples_to_remove.txt --keep-allele-order --make-bed --out ${FINAL_BFILE_PATH}_samplecleaned"
echo -e "-----------------------------------------------------------------\n"

echo "Deleting intermediate files to save space..."
find $qc_out_dir -type f -not -name "${FINAL_BFILE_NAME}.*" -delete
echo "Done."

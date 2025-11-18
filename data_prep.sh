#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

#-------------------------------------------------------------------------------
# USER CONFIGURATION
#-------------------------------------------------------------------------------
# This script assumes it is in the root of your project directory, with
# 'aux_data/' and 'libs/' as subdirectories.
#-------------------------------------------------------------------------------

# --- 1. Define Relative Paths ---
# Get the absolute path of the directory where this script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
start_dir=${PWD} # Save the directory where the user *ran* the script

# Define paths to auxiliary data and libraries relative to the script
AUX_DATA_DIR="${SCRIPT_DIR}/aux_data"
LIBS_DIR="${SCRIPT_DIR}/libs"

# --- 2. Define Paths to Reference Files ---
STRAND_CHECK_SCRIPT="${LIBS_DIR}/HRC-1000G-check-bim.pl"
HRC_REFERENCE_GZ="${AUX_DATA_DIR}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"

#-------------------------------------------------------------------------------
# SCRIPT ARGUMENTS
#-------------------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file_prefix> <data_directory>"
    echo "Example: $0 mycohort /path/to/my/data"
    echo "  (This will look for /path/to/my/data/mycohort.ped and .map)"
    exit 1
fi

file_prefix=$1 # Prefix of input ped/map files (e.g., "mycohort")
data_dir=$2    # Directory for plink files (e.g., "/path/to/data")

# Full path to the file prefix
bfile_path="${data_dir}/${file_prefix}"

echo "Starting data preparation..."
echo "Data directory: ${data_dir}"
echo "File prefix: ${file_prefix}"
echo "Script directory: ${SCRIPT_DIR}"

#-------------------------------------------------------------------------------
# SCRIPT EXECUTION
#-------------------------------------------------------------------------------

# Step 1: Convert ped/map to bed/fam/bim files if not already done
if [ -f "${bfile_path}.bed" ]; then
    echo "bfile already exists, skipping the conversion step."
else
    echo "Converting ped/map to bed/fam/bim files..."
    plink --file ${bfile_path} --make-bed --out ${bfile_path}
fi

# Step 2: Create a frequency file
echo "Creating frequency file (.frq)..."
plink --freq --bfile ${bfile_path} --out ${bfile_path}
echo "FINISHED bfile conversion and frequency calculation."

# Step 3: Run Will's pre-qc script (HRC-1000G-check-bim.pl)
# This script needs to be in the same directory as the data, so we copy it.
echo "Preparing for strand check..."
cp "${STRAND_CHECK_SCRIPT}" "${data_dir}/"

# Unzip the reference file *to* the data directory, not in aux_data
hrc_ref_unzipped_name=$(basename "${HRC_REFERENCE_GZ}" .gz)
hrc_ref_unzipped_path="${data_dir}/${hrc_ref_unzipped_name}"
echo "Unzipping HRC reference to ${hrc_ref_unzipped_path}..."
gunzip -c "${HRC_REFERENCE_GZ}" > "${hrc_ref_unzipped_path}"

# Enter the data directory to run the perl script
# The perl script generates Run-plink.sh in the *current* directory
echo "ENTERING data directory: ${data_dir}"
cd "${data_dir}"

echo "RUNNING strand check QC script (HRC-1000G-check-bim.pl)..."
# Note: Now that we are in data_dir, we use relative names
perl HRC-1000G-check-bim.pl -b "${file_prefix}.bim" -f "${file_prefix}.frq" -r "${hrc_ref_unzipped_name}" -h

echo "Running the generated 'Run-plink.sh' script..."
sh Run-plink.sh
echo "QC script is finished."

# Step 4: Clean up intermediate files
echo "Cleaning up temporary files..."
rm HRC-1000G-check-bim.pl
rm "${hrc_ref_unzipped_name}"
rm Run-plink.sh
# The perl script also creates other logs, like Exclude-*.txt, Force-Allele1-*.txt, etc.
# We can remove them to keep the directory clean.
rm -f ${file_prefix}-updated-chr*

# Return to the original starting directory
cd "${start_dir}"
echo "Returned to ${start_dir}"

#-------------------------------------------------------------------------------
# PREP COMPLETE
#-------------------------------------------------------------------------------
# The 'Run-plink.sh' script typically creates files with '-updated' suffix
FINAL_BFILE="${data_dir}/${file_prefix}-updated"

echo -e "\n-----------------------------------------------------------------"
echo -e "Data preparation is done."
echo -e "The output file for the main QC script is:"
echo -e "  ${FINAL_BFILE}"
echo -e "\nUse this path as the <bfile_in> argument for qc_main.sh"
echo -e "-----------------------------------------------------------------"

# Genotype QC Pipeline

A standardized BASH and R pipeline for performing Quality Control (QC) on human genotype data. This workflow starts from raw PED/MAP files, prepares them using the Michigan Imputation Server pre-phasing checks, and then runs a comprehensive QC analysis (sample and variant QC) using PLINK.

## 1. Installation & Setup

### 1.1. Get the Pipeline

First, clone this repository to your local machine:

```bash
git clone https://github.com/fbyukgo/genotype_qc.git
cd genotype_qc
```

### 1.2. Create the Conda Environment

This pipeline uses plink, plink2, bcftools, and R. All dependencies are defined in the `environment.yml` file.

Recommended (using Mamba):

```bash
mamba env create -f environment.yml
```

If you use Conda:

```bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate genotype_qc
```

You must activate this environment every time you run the pipeline.

### 1.3. Add Reference Data


This pipeline requires a set of reference data files (`aux_data/`). These files are hosted on GitHub Releases to keep the main repository lightweight.

**You must download and extract these files before running the pipeline.**

**Important note: the reference files in this directory are all in genome build GRCh37**

You can do this in two steps from your Linux shell.

1.  **Download the Data:**
    First, you need the URL for the release file.
    * Go to the [Releases page](https://github.com/fbyukgo/genotype_qc/releases) for this repository.
    * Find the latest release, right-click on the `aux_data_v1.tar.gz` asset, and select "Copy Link Address".

    Now, use `wget` in your terminal to download it. (Paste the URL you just copied).

    ```bash 
    wget https://github.com/fbyukgo/genotype_qc/releases/download/1.0/aux_data_v1.tar.gz
    ```

2.  **Unzip the Data:**
    Once the download is finished, run this command to unzip the file. This will create the `aux_data/` directory.

    ```bash
    # This will create the 'aux_data/' directory
    tar -xzvf aux_data_v1.tar.gz
    # we can delete the gzipped archive
    rm aux_data_v1.tar.gz 
    ```

After these two steps, you will have the `aux_data/` folder, and the pipeline is ready to run.


- aux_data/1Kg.bed  
- aux_data/1Kg.bim  
- aux_data/1Kg.fam  
- aux_data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz (originally from ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz) 
- aux_data/Populationfile-1kg.txt  
- aux_data/complex_regions.txt  

---

## 2. Running the QC Pipeline

### Step 1: Prepare Your Input Data

Place your raw `.ped` and `.map` genotype files in a directory like:

```
/path/to/my/data/
```

### Step 2: Run Data Preparation (`data_prep.sh`)

Usage:

```bash
./data_prep.sh <file_prefix> <data_directory>
```

Example:

```bash
./data_prep.sh mycohort /path/to/my/data
```

This generates:

```
/path/to/my/data/mycohort-updated
```

### Step 3: Run the Main QC (`qc_main.sh`)

Usage:

```bash
./qc_main.sh <bfile_in> \
 <qc_out_dir> \
 <report_plot_dir> \
 <analysis_name> \
 <maf_threshold> \
 <missingness_threshold> \
 <expected_anc> \
 <genome_build>
```

Arguments:

- `<bfile_in>`: path to updated bfile  
- `<qc_out_dir>`: PLINK output directory  
- `<report_plot_dir>`: plots + reports  
- `<analysis_name>`:str: label for plots i.e "MyAnalysis" 
- `<maf_threshold>`: e.g., 0.01  
- `<missingness_threshold>`: e.g., 0.02  
- `<expected_anc>`: e.g. "EUR", options: str --> ["EUR", "ASN", "AMR", "AFR"] 
- `<genome_build>`:str: either "b37" or "b38" 

**when using b38, you also need to use b38 reference files (will be added to the releases soon)** 

Example:

```bash
./qc_main.sh /path/to/my/data/mycohort-updated \
 ./qc_output \
 ./qc_reports \
 "MyCohort" \
 0.01 \
 0.02 \
 "EUR" \
 "b37"

```

---

## 3. Understanding the Output

### QC Reports & Plots

Saved in:

```
./qc_reports
```

Check for:

- High missingness  
- Sex errors  
- Heterozygosity outliers  
- Relatedness (PI_HAT > 0.1875)  
- PCA/MDS outliers  

### Final Cleaned Data

Stored in:

```
./qc_output
```

Example:

```
mycohort-updated_miss0.02_hwe1e-6_mafgte0.01.{bed,bim,fam}
```

---

## Final Step: Manual Sample Removal

Create:

```bash
$ cat samples_to_remove.txt
FID_101 IID_101
FID_202 IID_202
```

Run:

```bash
plink --bfile ./qc_output/mycohort-updated_miss0.02_hwe1e-6_mafgte0.01 \
      --remove samples_to_remove.txt \
      --make-bed \
      --out ./qc_output/mycohort_FINAL_CLEANED
```

Final dataset:

```
mycohort_FINAL_CLEANED.bed/bim/fam
```

---

## Project Structure

```
/genotype_qc/
├── qc_main.sh
├── data_prep.sh
├── environment.yml
├── README.md
│
├── aux_data/
│   ├── 1Kg.bed
│   ├── 1Kg.bim
│   ├── 1Kg.fam
│   ├── HRC.r1-1.GRCh37...gz
│   ├── Populationfile-1kg.txt
│   └── complex_regions.txt
│
└── libs/
    ├── HRC-1000G-check-bim.pl
    └── plotting/
        ├── plot_het_2.R
        ├── plot_mds_2.R
        ├── plot_miss.R
        ├── plot_relatedness.R
        └── plot_sex.R
```
```


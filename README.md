# Genotype QC Pipeline

A standardized BASH and R pipeline for performing Quality Control (QC) on human genotype data. This workflow starts from raw PED/MAP files, prepares them using the Michigan Imputation Server pre-phasing checks, and then runs a comprehensive QC analysis (sample and variant QC) using PLINK.

## 1. Installation & Setup

### 1.1. Get the Pipeline

First, clone this repository to your local machine:

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
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

This pipeline requires several external reference files. Most are included in this repository in the `aux_data/` directory (i.e 1Kg ref files, population file for this and complex regions all in GRCh37 build).

However, one file, `HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz`, must be downloaded manually due to its size.

**Action Required:**

1.  **Download the HRC file:**
    * Go to the Michigan Imputation Server tools page: [https://imputationserver.readthedocs.io/en/latest/prepare-your-data/](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)
    * Find the "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" file.
    * Download it directly using `wget`:
    ```bash
    wget [https://www.well.ox.ac.uk/~wrayner/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz](https://www.well.ox.ac.uk/~wrayner/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz)
    ```

2.  **Move the file:**
    * Place the downloaded `HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz` file directly into your `aux_data/` directory.

After this step, your `aux_data/` directory should contain all the required files and is ready to use.
Place all required files inside the `aux_data/` directory:

- aux_data/1Kg.bed  
- aux_data/1Kg.bim  
- aux_data/1Kg.fam  
- aux_data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz  
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
             <missingness_threshold>
```

Arguments:

- `<bfile_in>`: path to updated bfile  
- `<qc_out_dir>`: PLINK output directory  
- `<report_plot_dir>`: plots + reports  
- `<analysis_name>`: label for plots  
- `<maf_threshold>`: e.g., 0.01  
- `<missingness_threshold>`: e.g., 0.02  

Example:

```bash
./qc_main.sh /path/to/my/data/mycohort-updated \
             ./qc_output \
             ./qc_reports \
             "MyCohort" \
             0.01 \
             0.02
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

```
samples_to_remove.txt
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
/your-repo-name/
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


configfile: "config.yaml"

# --- Helper Variables ---
PREFIX = config["analysis_name"]
BASE_NAME = config["bfile_in"].split("/")[-1]

# Construct the final filename based on params to match your script logic
FINAL_NAME = f"{BASE_NAME}_miss{config['params']['final_miss']}_hwe{config['params']['hwe']}_mafgte{config['params']['final_maf']}"

# --- Target Rule ---
rule all:
    input:
        # 1. The Final PLINK files (Using double brackets {{ }} to escape f-string)
        expand(f"{config['out_dir']}/{FINAL_NAME}.{{ext}}", ext=['bed', 'bim', 'fam']),

        # 2. Plots (Sample QC)
        f"{config['plot_dir']}/{PREFIX}_missingness_plot.done",
        f"{config['plot_dir']}/{PREFIX}_sexcheck_plot.done",
        f"{config['plot_dir']}/{PREFIX}_het_plot.done",
        f"{config['plot_dir']}/{PREFIX}_relatedness_plot.done",
        f"{config['plot_dir']}/{PREFIX}_MDS_plot.done",

        # 3. Cleanup Trigger (Note the comma on the line above!)
        f"{config['out_dir']}/cleanup.done"
#-------------------------------------------------------------------------------
# Step 0: Pre-processing
#-------------------------------------------------------------------------------

rule preprocess:
    input:
        bed = config["bfile_in"] + ".bed",
        bim = config["bfile_in"] + ".bim",
        fam = config["bfile_in"] + ".fam"
    output:
        bed = temp(f"{config['out_dir']}/{BASE_NAME}_miss0.1.bed"),
        bim = temp(f"{config['out_dir']}/{BASE_NAME}_miss0.1.bim"),
        fam = temp(f"{config['out_dir']}/{BASE_NAME}_miss0.1.fam")
    params:
        geno = config["params"]["pre_geno"],
        mind = config["params"]["pre_mind"],
        build = config["genome_build"],
        out_step1 = f"{config['out_dir']}/{BASE_NAME}_miss0.1_wopseudo_rsid_split",
        out_step2 = f"{config['out_dir']}/{BASE_NAME}_rsid_miss0.1",
        out_final = f"{config['out_dir']}/{BASE_NAME}_miss0.1"
    shell:
        """
        # 1. Filter missingness
        plink --bfile {config[bfile_in]} --allow-no-sex --geno {params.geno} --mind {params.mind} --make-bed --out {params.out_step1}
        
        # 2. Split X
        plink --bfile {params.out_step1} --split-x {params.build} no-fail --make-bed --out {params.out_step2}
        
        # 3. Standardize IDs
        # FIX: We use single quotes around the ID format string so Bash doesn't mess up the $r or []
        plink2 --bfile {params.out_step2} --set-all-var-ids '@:#[{params.build}]$r,$a' --make-bed --out {params.out_final}
        
        # Clean intermediates
        rm {params.out_step1}.* {params.out_step2}.*
        """


#-------------------------------------------------------------------------------
# Step 1: Missingness
#-------------------------------------------------------------------------------
rule check_missingness:
    input:
        bed = rules.preprocess.output.bed,
        bim = rules.preprocess.output.bim,
        fam = rules.preprocess.output.fam
    output:
        imiss = temp(f"{config['out_dir']}/{BASE_NAME}_missing.imiss"),
        lmiss = temp(f"{config['out_dir']}/{BASE_NAME}_missing.lmiss"),
        plot_flag = f"{config['plot_dir']}/{PREFIX}_missingness_plot.done"
    params:
        out_prefix = f"{config['out_dir']}/{BASE_NAME}_missing"
    shell:
        """
        plink --bfile {config[out_dir]}/{BASE_NAME}_miss0.1 --missing --out {params.out_prefix}
        
        Rscript {config[plot_script_dir]}/plot_miss.R {params.out_prefix} {PREFIX} {config[plot_dir]} {config[params][final_miss]}
        touch {output.plot_flag}
        """

#-------------------------------------------------------------------------------
# Step 2: Sex Check
#-------------------------------------------------------------------------------
rule check_sex:
    input:
        bed = rules.preprocess.output.bed,
        bim = rules.preprocess.output.bim,
        fam = rules.preprocess.output.fam
    output:
        sexcheck = temp(f"{config['out_dir']}/{BASE_NAME}_sexcheck.sexcheck"),
        plot_flag = f"{config['plot_dir']}/{PREFIX}_sexcheck_plot.done"
    params:
        out_sex = f"{config['out_dir']}/{BASE_NAME}_sexcheck",
        out_x = f"{config['out_dir']}/{BASE_NAME}_xchr",
        out_x_miss = f"{config['out_dir']}/{BASE_NAME}_xchr-missing"
    shell:
        """
        plink --bfile {config[out_dir]}/{BASE_NAME}_miss0.1 --check-sex --out {params.out_sex}
        
        # X-chr specific missingness
        plink --bfile {config[out_dir]}/{BASE_NAME}_miss0.1 --chr 23 --make-bed --out {params.out_x}
        plink --bfile {params.out_x} --missing --out {params.out_x_miss}
        
        Rscript {config[plot_script_dir]}/plot_sex.R {params.out_sex} {params.out_x_miss} {PREFIX} {config[plot_dir]}
        
        # Clean up X chr intermediate
        rm {params.out_x}.*
        touch {output.plot_flag}
        """

#-------------------------------------------------------------------------------
# Step 3: Heterozygosity & LD Pruning
#-------------------------------------------------------------------------------
# This rule creates the LD pruned file used for Relatedness and Ancestry
rule ld_pruning_and_het:
    input:
        bed = rules.preprocess.output.bed,
        bim = rules.preprocess.output.bim,
        fam = rules.preprocess.output.fam
    output:
        pruned_bed = f"{config['out_dir']}/{BASE_NAME}_ldpruned.bed",
        pruned_bim = f"{config['out_dir']}/{BASE_NAME}_ldpruned.bim",
        pruned_fam = f"{config['out_dir']}/{BASE_NAME}_ldpruned.fam",
        plot_flag = f"{config['plot_dir']}/{PREFIX}_het_plot.done"
    # ... params and shell remain same ...
    params:
        maf = config["params"]["final_maf"],
        complex = config["complex_regions"],
        base_prefix = f"{config['out_dir']}/{BASE_NAME}_miss0.1",
        out_prefix = f"{config['out_dir']}/{BASE_NAME}",
        ld_win = config["params"]["ld_window"],
        ld_step = config["params"]["ld_step"],
        ld_r2 = config["params"]["ld_r2"]
    shell:
        """
        # 1. Extract Autosomes
        plink --bfile {params.base_prefix} --autosome --make-bed --out {params.out_prefix}_chr1-22
        
        # 2. Split by MAF (High vs Low)
        plink --bfile {params.out_prefix}_chr1-22 --maf {params.maf} --make-bed --out {params.out_prefix}_highMAF
        plink --bfile {params.out_prefix}_chr1-22 --exclude {params.out_prefix}_highMAF.bim --make-bed --out {params.out_prefix}_lowMAF
        
        # 3. Prune and Het (High MAF)
        plink --bfile {params.out_prefix}_highMAF --exclude {params.complex} --range --make-bed --out {params.out_prefix}_highMAF_nocr
        plink --bfile {params.out_prefix}_highMAF_nocr --indep-pairwise {params.ld_win} {params.ld_step} {params.ld_r2} --out {params.out_prefix}_highMAF_pruning
        plink --bfile {params.out_prefix}_highMAF_nocr --extract {params.out_prefix}_highMAF_pruning.prune.in --het --out {params.out_prefix}_highMAF_het
        
        # 4. Prune and Het (Low MAF)
        plink --bfile {params.out_prefix}_lowMAF --exclude {params.complex} --range --make-bed --out {params.out_prefix}_lowMAF_nocr
        plink --bfile {params.out_prefix}_lowMAF_nocr --indep-pairwise {params.ld_win} {params.ld_step} {params.ld_r2} --out {params.out_prefix}_lowMAF_pruning
        plink --bfile {params.out_prefix}_lowMAF_nocr --extract {params.out_prefix}_lowMAF_pruning.prune.in --het --out {params.out_prefix}_lowMAF_het
        
        # 5. Create the final LD pruned file for next steps (using High MAF set)
        plink --bfile {params.out_prefix}_highMAF_nocr --extract {params.out_prefix}_highMAF_pruning.prune.in --make-bed --out {config[out_dir]}/{BASE_NAME}_ldpruned

        # 6. Plotting
        Rscript {config[plot_script_dir]}/plot_het_2.R {params.out_prefix}_highMAF_het mafgte{params.maf} {PREFIX} {config[plot_dir]}
        Rscript {config[plot_script_dir]}/plot_het_2.R {params.out_prefix}_lowMAF_het mafless{params.maf} {PREFIX} {config[plot_dir]}
        
        touch {output.plot_flag}
        """

#-------------------------------------------------------------------------------
# Step 4: Relatedness
#-------------------------------------------------------------------------------
rule check_relatedness:
    input:
        bed = rules.ld_pruning_and_het.output.pruned_bed,
        bim = rules.ld_pruning_and_het.output.pruned_bim,
        fam = rules.ld_pruning_and_het.output.pruned_fam
    output:
        genome = temp(f"{config['out_dir']}/{BASE_NAME}_ldpruned_genome.genome"),
        plot_flag = f"{config['plot_dir']}/{PREFIX}_relatedness_plot.done"
    params:
        out_prefix = f"{config['out_dir']}/{BASE_NAME}_ldpruned_genome"
    shell:
        """
        plink --bfile {config[out_dir]}/{BASE_NAME}_ldpruned --genome --out {params.out_prefix}
        
        # Tab convert
        awk -v OFS='\t' '{{$1=$1}}1' {params.out_prefix}.genome > {params.out_prefix}_tab.genome
        
        Rscript {config[plot_script_dir]}/plot_relatedness.R {params.out_prefix}_tab {PREFIX} {config[plot_dir]}
        touch {output.plot_flag}
        """

#-------------------------------------------------------------------------------
# Step 5: Ancestry (MDS)
#-------------------------------------------------------------------------------

rule prep_1kg:
    input:
        bed = config["onekg_bfile"] + ".bed"
    output:
        # REMOVED temp() wrappers here so files persist
        bed = f"{config['out_dir']}/1kg_ref.bed",
        bim = f"{config['out_dir']}/1kg_ref.bim",
        fam = f"{config['out_dir']}/1kg_ref.fam"
    params:
        build = config["genome_build"],
        out_base = f"{config['out_dir']}/1kg_ref_base",
        out_final = f"{config['out_dir']}/1kg_ref",
        geno = config["params"]["pre_geno"],
        mind = config["params"]["pre_mind"],
        maf = config["params"]["final_maf"]
    shell:
        """
        # Clean 1KG data
        plink --bfile {config[onekg_bfile]} --set-missing-var-ids '@:#[{params.build}]$1,$2' --make-bed --out {params.out_base}
        
        plink --bfile {params.out_base} --geno {params.geno} --mind {params.mind} --maf {params.maf} --allow-no-sex --make-bed --out {params.out_base}_clean
        
        # Standardize IDs to match our data
        plink2 --bfile {params.out_base}_clean --set-all-var-ids '@:#[{params.build}]$r,$a' --make-bed --out {params.out_final}
        
        rm {params.out_base}*
        """

rule ancestry_mds:
    input:
        # FIXED: Explicitly request all 3 LD-pruned files so Snakemake doesn't delete them
        pruned_bed = rules.ld_pruning_and_het.output.pruned_bed,
        pruned_bim = rules.ld_pruning_and_het.output.pruned_bim,
        pruned_fam = rules.ld_pruning_and_het.output.pruned_fam,
        
        # Reference files
        ref_bed = rules.prep_1kg.output.bed,
        ref_bim = rules.prep_1kg.output.bim,
        ref_fam = rules.prep_1kg.output.fam,
        pop_file = config["onekg_pop_file"]
    output:
        mds = f"{config['out_dir']}/{BASE_NAME}_1KG_merged_mds.mds",
        plot_flag = f"{config['plot_dir']}/{PREFIX}_MDS_plot.done"
    params:
        out_dir = config["out_dir"],
        ref_prefix = f"{config['out_dir']}/1kg_ref",
        our_prefix = f"{config['out_dir']}/{BASE_NAME}_ldpruned",
        merged_prefix = f"{config['out_dir']}/{BASE_NAME}_1KG_merged",
        panel_merged = f"{config['out_dir']}/population_panel_merged.txt",
        anc = config["params"]["expected_ancestry"]
    shell:
        """
        # 1. Prepare Panel File
        cat {config[bfile_in]}.fam | awk '{{print $1, $2, "OWN"}}' > {params.out_dir}/population_panel_our.txt
        cat {params.out_dir}/population_panel_our.txt {input.pop_file} | sed -e '1i\FID IID ANC' > {params.panel_merged}

        # 2. Intersect
        plink --bfile {params.ref_prefix} --extract {params.our_prefix}.bim --recode --make-bed --out {params.ref_prefix}_intersect
        plink --bfile {params.our_prefix} --extract {params.ref_prefix}_intersect.bim --recode --make-bed --out {params.our_prefix}_intersect
        
        # 3. Merge (With Retry Logic)
        if ! plink --bfile {params.our_prefix}_intersect -bmerge {params.ref_prefix}_intersect --make-bed --out {params.merged_prefix}; then
            echo "Merge failed. Checking for .missnp file..."
            if [ -f "{params.merged_prefix}.missnp" ]; then
                echo "Removing incompatible variants and retrying..."
                plink --bfile {params.our_prefix}_intersect --exclude {params.merged_prefix}.missnp --make-bed --out {params.our_prefix}_intersect_retry
                plink --bfile {params.ref_prefix}_intersect --exclude {params.merged_prefix}.missnp --make-bed --out {params.ref_prefix}_intersect_retry
                plink --bfile {params.our_prefix}_intersect_retry -bmerge {params.ref_prefix}_intersect_retry --make-bed --out {params.merged_prefix}
            else
                echo "Merge failed with no .missnp file. Exiting."
                exit 1
            fi
        fi
        
        # 4. MDS
        plink --bfile {params.merged_prefix} --genome --out {params.merged_prefix}_genome
        plink --bfile {params.merged_prefix} --read-genome {params.merged_prefix}_genome.genome --cluster --mds-plot 10 --out {params.merged_prefix}_mds
        
        # 5. Plot
        Rscript libs/plotting/plot_mds_2.R {params.merged_prefix}_mds {PREFIX} {config[plot_dir]} {params.panel_merged} {params.anc}
        touch {output.plot_flag}
        """



#-------------------------------------------------------------------------------
# Step 6-8: Final Variant QC (On Original Data, not LD Pruned)
#-------------------------------------------------------------------------------
rule final_variant_qc:
    input:
        bed = rules.preprocess.output.bed,
        bim = rules.preprocess.output.bim,
        fam = rules.preprocess.output.fam
    output:
        bed = f"{config['out_dir']}/{FINAL_NAME}.bed",
        bim = f"{config['out_dir']}/{FINAL_NAME}.bim",
        fam = f"{config['out_dir']}/{FINAL_NAME}.fam"
    params:
        base = f"{config['out_dir']}/{BASE_NAME}_miss0.1",
        out_final = f"{config['out_dir']}/{FINAL_NAME}",
        miss = config["params"]["final_miss"],
        hwe = config["params"]["hwe"],
        maf = config["params"]["final_maf"]
    shell:
        """
        plink --bfile {params.base} --geno {params.miss} --hwe {params.hwe} --maf {params.maf} --make-bed --out {params.out_final}
        """
#-------------------------------------------------------------------------------
# Step 9: Cleanup
#-------------------------------------------------------------------------------
rule cleanup:
    input:
        # We use the Bed file and Plot flags to ensure this runs LAST
        bed = rules.final_variant_qc.output.bed,
        mds_plot = rules.ancestry_mds.output.plot_flag,
        rel_plot = rules.check_relatedness.output.plot_flag
    output:
        touch(f"{config['out_dir']}/cleanup.done")
    params:
        dir = config["out_dir"],
        final_name = FINAL_NAME,
        ref_name = "1kg_ref"
    shell:
        """
        # Find all files in output dir
        # Exclude final files, reference files, and the cleanup flag
        # Delete everything else
        find {params.dir} -type f \
            -not -name "{params.final_name}.*" \
            -not -name "cleanup.done" \
            -delete
        """

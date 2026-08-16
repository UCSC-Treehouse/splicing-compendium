# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on

# Run from project root: snakemake --snakefile filter_min_junc.smk --profile pheonix-profile

import os
import pandas as pd
from datetime import datetime

configfile: "config/test_minreads_junctions.yaml"

# read in configfile values
sample_table = pd.read_table(config["sample_sheet"])
SAMPLES = sample_table["sample"].tolist()
GROUPS = sample_table["group"].tolist()
REF_GTF = config["reference_gtf"]

# make a dictionary to map samples to groups so merged GTF output won't need group wildcards
sample_to_group = dict(zip(SAMPLES, GROUPS))

pathvars:
    test_results_dir = f"results/test_minreads_psi_calc"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<test_results_dir>/psi"

rule filter_junctions:
    input: "<test_results_dir>/merged_junctions.bed"
    output: "<test_results_dir>/minread_filtered_junc/filtered_merged_junctions.bed"
    threads: 8
    resources:
        mem_mb = 2200000,
        runtime = 720
    shell:
        """
        # remove rows corresponding to junctions where all samples have reads < 20
        # make output tab-separated
        awk 'BEGIN{FS=OFS="\t"}
        # keep header row
        NR==1 {print; next}
        {
            # instantiate number of rows to keep
            keep=0
            # loop through column 5-last column and check if count > 20
            for (i=5; i<=NF; i++) {
                if ($i+0 >= 20) { keep=1; break }
            }
            # if one column has count of at least 20, keep the pos_id row
            if (keep) print
        }' {input} > {output}
        """

rule calculate_psi:
    input:
        events_dir = "<test_results_dir>/events",
        merged_junctions = "<test_results_dir>/minread_filtered_junc/filtered_merged_junctions.bed",
        merged_gtf = "<test_results_dir>/merged_gtf.gtf",
    output:
        shiba_psi_out = directory("<test_results_dir>/psi")
    params:
        shiba_scripts = config["shiba_scripts_path"],
        min_reads = config["shiba_min_reads"]
    priority: 1
    threads: 8
    resources:
        mem_mb = 2200000,
        runtime = 20160
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/psi.py -m {params.min_reads} -p {threads} -v --onlypsi {input.merged_junctions} {input.events_dir} {output.shiba_psi_out}

        # zip PSI results
        pigz -p {threads} \
        {output.shiba_psi_out}/*.txt

        """

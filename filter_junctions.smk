# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on

# Run from project root: snakemake --snakefile filter_junctions.smk --profile pheonix-profile

import os
import pandas as pd
from datetime import datetime

configfile: "config/test_split_junctions.yaml"

# read in configfile values
sample_table = pd.read_table(config["sample_sheet"])
SAMPLES = sample_table["sample"].tolist()
GROUPS = sample_table["group"].tolist()
REF_GTF = config["reference_gtf"]

# make a dictionary to map samples to groups so merged GTF output won't need group wildcards
sample_to_group = dict(zip(SAMPLES, GROUPS))

pathvars:
    merged_shiba_results = f"results/merged_shiba/splice_compendium_v1",
    test_results_dir = f"results/test_split_psi_calc"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<test_results_dir>/psi"

rule filter_junctions:
    input: "<merged_shiba_results>/merged_junctions.bed"
    output: "<test_results_dir>/junctions_by_chromosome/chr1_merged_junctions.bed"
    threads: 8
    resources:
        mem_mb = 5000,
        runtime = 60
    shell:
        """
        # filter for chr1 value under chromosome column (first column)
        awk -F'\t' 'NR==1 || $1=="chr1"' {input} > {output}
        """

rule calculate_psi:
    input:
        events_dir = "<merged_shiba_results>/events",
        merged_junctions = "<test_results_dir>/junctions_by_chromosome/chr1_merged_junctions.bed",
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf",
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

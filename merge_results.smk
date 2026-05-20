# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on
# Run from project root: snakemake --snakefile merge_results.smk -j 15

import os
import pandas as pd
from datetime import datetime

configfile: "config/merge_shiba_config.yaml"

GENOME_ID = config["genome_id"]
SAMPLES = pd.read_table(config["sample_sheet"])["samples"].tolist()
GROUPS = pd.read_table(config["sample_sheet"])["group"].tolist()
VERSION = config["version"]

# these pathvars are from our "separate shiba runs" snakemake in main
# the main snakefile will also need to be changed to reflect how we handle sample groups
pathvars:
    merged_shiba_results = os.path.join("results", "merged_shiba", VERSION),
    references = "references"

# path to group shiba results
GROUP_SHIBA_RESULTS = expand(
    os.path.join("results", "{group}", "shiba"),
    group=GROUPS
)

# path to junction.bed files
JUNCTION_BEDS = expand(
    os.path.join(
        "{group_shiba_result}",
        "{sample}/junctions/junctions.bed"
    ),
    group_shiba_result=GROUP_SHIBA_RESULTS,
    sample=SAMPLES
)

# path to sample gtfs from separate shiba runs
SAMPLE_GTFS = expand(
    os.path.join(
        "{group_shiba_result}",
        "{sample}/annotation/assembled_annotation.gtf.gz"
    ),
    group_shiba_result=GROUP_SHIBA_RESULTS,
    sample=SAMPLES
)

# path to unzipped sample gtfs
UNZIPPED_GTFs = expand(
    os.path.join(
        "{group_shiba_result}",
        "{sample}/annotation/assembled_annotation.gtf"
    ),
    group_shiba_result=GROUP_SHIBA_RESULTS,
    sample=SAMPLES
)

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed",

rule merge_junctions:
    input: JUNCTION_BEDS
    output: "<merged_shiba_results>/merged_junctions.bed"
    # compute joined junctions string prior to passing into join script
    params:
        junctions=lambda wildcards, input: ",".join(input)
    priority: 1
    threads: 4
    shell:
        """
        Rscript scripts/03-merge_separate_junctions.R --junctions={params.junctions} --output={output}
        """

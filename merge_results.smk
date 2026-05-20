# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all samples with separate Shiba results is given to the config file to identify files to run the workflow on
# Otherwise, this workflow will operate on all samples in the shiba results directory
# Run from project root: snakemake --snakefile merge_results.smk -j 15

import os
import pandas as pd
from datetime import datetime

configfile: "config/merge_shiba_config.yaml"

GENOME_ID = config["genome_id"]
SAMPLE_GROUP = config["sample_group"]
if config.get("sample_sheet"):
    SAMPLES = pd.read_table(config["sample_sheet"])["samples"].tolist()
else: # this part needs to be changed as we sometimes need to junctions and sometimes the GTFs - maybe make sample sheet required?
    SAMPLES, = glob_wildcards(os.path.join("results", SAMPLE_GROUP, "shiba", "{sample}", "annotation", "assembled_annotation.gtf.gz"))

# a timestamp might be fragile if I run parts of the pipeline over several days
VERSION = "test" # probably move this to configfile

# these pathvars are from our "separate shiba runs" snakemake in main
pathvars:
    data = os.path.join("data", SAMPLE_GROUP),
    reports = os.path.join("reports", SAMPLE_GROUP),
    results = os.path.join("results", SAMPLE_GROUP),
    shiba_results = "<results>/shiba",
    merged_shiba_results = os.path.join("<shiba_results>/merged_results", VERSION),
    logs = os.path.join("logs", SAMPLE_GROUP),
    references = "references"

# path to junction.bed files
JUNCTION_BEDS = expand(
    os.path.join(
        "<shiba_results>",
        "{sample}/junctions/junctions.bed"
    ),
    sample=SAMPLES
)

# path to sample gtfs from separate shiba runs
SAMPLE_GTFS = expand(
    os.path.join(
        "<shiba_results>",
        "{sample}/annotation/assembled_annotation.gtf.gz"
    ),
    sample=SAMPLES
)

# path to unzipped sample gtfs
UNZIPPED_GTFs = expand(
    os.path.join(
        "<shiba_results>",
        "{sample}/annotation/assembled_annotation.gtf"
    ),
    sample=SAMPLES
)

# path to junctions manifest file to pass into junction merging rule
# How do I get it to grab both TARGET and GTEx sample groups?
# Another rule to merge manifest files from multiple groups together?
# What's best for if we want this script to be usable for splice compendium updates? e.g. adding new datasets to compendium
# Example: Adding CBTN + GTEx brain
JCN_MANIFEST = "<merged_shiba_results>/junction_manifest.tsv"

# path to GTF manifest file
GTF_MANIFEST = "<merged_shiba_results>/gtf_manifest.tsv"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed",
        merged_gtf = "<merged_shiba_results>/merged_gtf.bed"

# Don't do manifest, just pass in joined samples in the other example in the PR
# update bed joining code in the other branch too
rule merge_junctions:
    input: JUNCTION_BEDS
    output: "<merged_shiba_results>/merged_junctions.bed"
    priority: 1
    threads: 4
    shell:
        """
        Rscript scripts/03-merge_separate_junctions.R --junctions={",".join(input)} --output={output}
        """

rule unzip_gtfs:
    input:
        SAMPLE_GTFS
    output:
        temp(UNZIPPED_GTFs)
    priority: 1
    threads: 4
    shell:
        """
        # keep zipped input so we don't need to zip the GTFs again
        gunzip {input} --keep
        """

# don't give this a manifest, pass it a list of GTF paths (avoid copying files across shared filesystems)
rule merge_gtfs:
    input:
        sample_gtfs = UNZIPPED_GTFs,
        reference_gtf = f"<references>/{GENOME_ID}.annotation.gtf",
    output:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 4
    shell:
        """
        # make list of all input files and print it into manifest one line at a time
        echo "{'\n'.join(input.manifest)}" > gtf_manifest.txt
        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} gtf_manifest.txt
        """

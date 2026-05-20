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
else:
    SAMPLES, = glob_wildcards(os.path.join("results", SAMPLE_GROUP, "shiba", "{sample}", "annotation", "assembled_annotation.gtf.gz"))

# a timestamp might be fragile if I run parts of the pipeline over several days
VERSION = "test"

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

# path to junctions manifest file to pass into junction merging rule
# How do I get it to grab both TARGET and GTEx sample groups?
# Another rule to merge manifest files from multiple groups together?
# What's best for if we want this script to be usable for splice compendium updates? e.g. adding new datasets to compendium
# Example: Adding CBTN + GTEx brain
JCN_MANIFEST = "<merged_shiba_results>/junction_manifest.tsv"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<merged_shiba_results>/merged_junctions.bed"

rule make_junction_manifest:
    input:
        JUNCTION_BEDS
    output:
        JCN_MANIFEST
    run:
        rows = []
        for f in input:
            # Extract sample name from path
            # results/<group>/shiba/<sample>/junctions/junctions.bed
            sample = f.split(os.sep)[3]

            rows.append({
                "sample": sample,
                "junction_bed": f
            })

        # write dataframe into a tsv to pass into merge script
        df = pd.DataFrame(rows)
        df.to_csv(output[0], sep="\t", index=False)

rule merge_junctions:
    input: JCN_MANIFEST
    output: "<merged_shiba_results>/merged_junctions.bed"
    priority: 1
    threads: 4
    shell:
        """
        Rscript scripts/03-merge_separate_junctions.R --junctions={input} --output={output}
        """

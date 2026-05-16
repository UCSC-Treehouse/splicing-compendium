# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all samples with separate Shiba results is given to the config file to identify files to run the workflow on
# Otherwise, this workflow will operate on all samples in the shiba results directory
# Run from project root: snakemake --snakefile merge_results.smk -j 15

import os
import pandas as pd

configfile: "config/merge_shiba_config.yaml"

GENOME_ID = config["genome_id"]
SAMPLE_GROUP = config["sample_group"]
if config.get("sample_sheet"):
    SAMPLES = pd.read_table(config["sample_sheet"])["samples"].tolist()
else:
    SAMPLES, = glob_wildcards(os.path.join("results", SAMPLE_GROUP, "shiba", "{sample}", "annotation", "assembled_annotation.gtf.gz"))

pathvars:
    data = os.path.join("data", SAMPLE_GROUP),
    reports = os.path.join("reports", SAMPLE_GROUP),
    merged_results = os.path.join("results", "merged_results"),
    shiba_results = os.path.join(results, "shiba"),
    logs = os.path.join("logs", SAMPLE_GROUP),
    references = "references"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        expand("<results>/shiba/{sample}/", sample = SAMPLES)

# Fifth rule: merge junction bedfiles
rule merge_junctions:
    input: "<shiba_results>/{sample}/junctions/junctions.bed"
    output: "results/merged_shiba/timestamp/merged_junctions.bed"
    priority: 1
    threads: 4
    shell:
    """
    Rscript 03-merge_separate_junctions.R --junctions={input} --output={output}
    """

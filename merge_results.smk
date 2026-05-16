# snakefile for merging shiba separate results
# Build reference indices first if necessary: snakemake --snakefile build_references.smk --cores 15
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
    results = os.path.join("results", SAMPLE_GROUP),
    shiba_results = os.path.join(results, "shiba"),
    logs = os.path.join("logs", SAMPLE_GROUP),
    references = "references"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        expand("<results>/shiba/{sample}/", sample = SAMPLES)

# Fifth rule: merge GTFs from separate shiba runs
rule merge_separate_gtfs:
    input:
        reference_gtf = f"<references>/{GENOME_ID}.annotation.gtf"
        sample_gtf = "<shiba_results>/{sample}/annotation/assembled_annotation.gtf.gz"
    output:
        merged_gtf = "<results>/shiba/merged_gtf.gtf"
    priority: 1
    log:
        "<log>/merge_gtf/merge_gtf.log"
    threads: 7
    resources:
        mem_mb = 60000
    shell:
    """
    # create output directory for merged gtf
    mkdir -p "<results>/shiba/merged"

    # unzip gtfs so stringtie can operate on them
    gunzip {input.sample_gtf}

    # make list of GTFs to merge
    echo 'sample\tbam_path\tgroup\ttechnology' > {params.experiment_table}

    # merge GTFs in sample sheet with stringtie
    stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} {input.gtf_list}

				"stringtie",
				"-p", str(num_processors),
				"-G", reference_gtf,
				"-o", sample_gtf

    """

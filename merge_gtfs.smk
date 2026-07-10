# Snakefile for generating and merging GTFs of a defined set of samples
# This snakefile uses a sample sheet (experiment.tsv) containing sample accessions and bam paths to identify files to act on
# To generate test files to assess differences caused by merging the reference in multiple times, 
# this workflow will run on experiment.tsvs of one sample at a time to generate a GTF.
# Then, the GTFs generated and merged with references separately will be merged together by passing merge_gtf_list.tsv into the merge_gtf rule.

# Run from project root: snakemake --snakefile merge_gtfs.smk --profile pheonix-profile

import pandas as pd
import os
from datetime import datetime

configfile: "config/merge_shiba_config.yaml"
merge_list: "exploration/merge_gtf_list.tsv"

# read in configfile values
SAMPLES = pd.read_table(config["sample_sheet"])["sample"].tolist()
GROUPS = pd.read_table(config["sample_sheet"])["group"].tolist()

pathvars:
    merged_shiba_results = f"results/merged_gtf_tests/{config["version"]}"

# path to sample gtfs from separate shiba runs
SAMPLE_GTFS = expand(
    "results/{group}/shiba/{sample}/annotation/assembled_annotation.gtf.gz",
    zip,
    group=GROUPS,
    sample=SAMPLES
)
SAMPLE_BAMS = expand(
    "data/{group}/{sample}/Aligned.sorted.out.bam",
    group=GROUPS,
    sample=SAMPLES
)

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<merged_shiba_results>/psi"

rule bam2gtf:
    input: 
        # sample bams
        sample_bam = placeholder
    output:
        temp("annotation/{sample}.gtf")
    threads:
        20
    shell:
        """
        stringtie -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {params.longread_option} \
        {input.bam} >& {log}
        """

rule merge_gtfs:
    input:
        # sample_gtfs are zipped
        sample_gtfs = SAMPLE_GTFS,
        reference_gtf = config["reference_gtf"],
        merge_list = pd.read_table(merge_list)
    output:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 4
    shell:
        """
        # merge gtfs with stringtie for splice analysis
        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} {input.merge_list}
        """
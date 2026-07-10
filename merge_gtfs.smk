# Snakefile for generating and merging GTFs of a defined set of samples
# This snakefile uses a sample sheet (experiment.tsv) containing sample accessions and bam paths to identify files to act on
# To generate test files to assess differences caused by merging the reference in multiple times,
# this workflow will run on experiment.tsvs of one sample at a time to generate a GTF.
# Then, the GTFs generated and merged with references separately will be merged together by passing merge_gtf_list.tsv into the merge_gtf rule.

# Run from project root: snakemake --snakefile merge_gtfs.smk --profile pheonix-profile

import pandas as pd
import os
from datetime import datetime

configfile: "config/merge_gtf_test_config.yaml"
merge_list: "exploration/merge_gtf_list.tsv"

# read in configfile values
SAMPLES = pd.read_table(config["sample_sheet"])["sample"].tolist()
GROUPS = pd.read_table(config["sample_sheet"])["group"].tolist()
REF_GTF = config["reference_gtf"]

# make a dictionary to map samples to groups so merged GTF output won't need group wildcards
sample_to_group = dict(zip(SAMPLES, GROUPS))

pathvars:
    merged_gtf_results = f"results/merged_gtf_tests/{config["version"]}"

# path to sample bams from separate shiba runs
SAMPLE_BAMS = expand(
    "data/{group}/star-output/{sample}/Aligned.sortedByCoord.out.bam",
    group=GROUPS,
    sample=SAMPLES
)

SAMPLE_GTFS =  expand(
    "<merged_gtf_results>/{sample}.gtf",
    sample = SAMPLES
)

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<merged_gtf_results>/merged_gtf.gtf"

rule bam2gtf:
    input:
        ref_gtf = REF_GTF,
        bam = lambda wildcard: (
            f"data/{sample_to_group[wildcard.sample]}"
            f"/star-output/{wildcard.sample}/Aligned.sortedByCoord.out.bam"
        )
    output: temp("<merged_gtf_results>/{sample}.gtf")
    threads: 1
    log: "<merged_gtf_results>/logs/{sample}_bam2gtf.log"
    shell:
        """
        stringtie -p {threads} -G {input.ref_gtf} -o {output} {input.bam} >& {log}
        """

rule merge_gtfs:
    input:
        reference_gtf = config["reference_gtf"],
        sample_gtfs = SAMPLE_GTFS
    output: "<merged_gtf_results>/merged_gtf.gtf"
    priority: 1
    threads: 1
    log: "<merged_gtf_results>/logs/merge_gtfs.log"
    shell:
        """
        # instantiate manifest as a temp file
        manifest=$(mktemp)

        # add gtf to manifest for merging
        for gtf in {input.sample_gtfs}; do
            # Add to the manifest
            echo "$gtf" >> $manifest
        done

        # merge gtfs with stringtie for splice analysis
        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output} $manifest >& {log}
        """

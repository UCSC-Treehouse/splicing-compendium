# Snakefile for generating and merging GTFs of a defined set of samples
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on
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
        reference_gtf = config["reference_gtf"]
    output:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 4
    shell:
        """
        # instantiate manifest as a temp file
        manifest=$(mktemp)

        # unzip input gtfs for stringtie
        for gz_gtf in {input.sample_gtfs}; do
            # create the uncompressed file path
            gtf="${{gz_gtf%.gz}}"
            # decompress (leaving the original) explicitly - decompressed files are in the same directory as the compressed files
            gunzip -c "$gz_gtf" > "$gtf"
            # Add to the manifest
            echo "$gtf" >> $manifest
        done

        # merge gtfs with stringtie for splice analysis
        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} $manifest

        # delete gtf files using the manifest
        xargs rm < $manifest

        # remove temporary manifest
        rm $manifest
        """
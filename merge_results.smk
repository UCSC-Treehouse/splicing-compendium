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
    merged_shiba_results = f"results/merged_shiba/{VERSION}",
    references = "references"

# path to junction.bed files
# need to zip paths so group/sample pairs are matched rowwise
JUNCTION_BEDS = expand(
    "results/{group}/shiba/{sample}/junctions/junctions.bed",
    zip,
    group=GROUPS,
    sample=SAMPLES
)

# path to sample gtfs from separate shiba runs
SAMPLE_GTFS = expand(
    "results/{group}/shiba/{sample}/annotation/assembled_annotation.gtf.gz",
    zip,
    group=GROUPS,
    sample=SAMPLES
)

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed",
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"

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

rule merge_gtfs:
    input:
        # sample_gtfs are zipped
        sample_gtfs = SAMPLE_GTFS,
        reference_gtf = f"<references>/{GENOME_ID}.annotation.gtf"
    output:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 4
    shell:
        """
        # instantiate manifest as a temp file
        manifest=$(mktemp)

        # make list of all input files and print it into manifest one line at a time
        echo "{params.gtf_manifest}" > $manifest

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

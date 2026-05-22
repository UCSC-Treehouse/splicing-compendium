# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on
# Run from project root: snakemake --snakefile merge_results.smk -j 15

import os
import pandas as pd
from datetime import datetime

configfile: "config/merge_shiba_config.yaml"

# read in configfile values
GENOME_ID = config["genome_id"]
SAMPLES = pd.read_table(config["sample_sheet"])["samples"].tolist()
GROUPS = pd.read_table(config["sample_sheet"])["group"].tolist()
VERSION = config["version"]
SHIBA_SCRIPTS = config["shiba_scripts_path"]
REFERENCE_GTF = config["reference_gtf"]
MIN_READS = config["min_reads"]

# the main snakefile will also need to be changed to reflect how we handle sample groups
pathvars:
    merged_shiba_results = f"results/merged_shiba/{VERSION}"

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
        shiba_out = "<merged_shiba_results>/psi"

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
        reference_gtf = REFERENCE_GTF
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

rule gtf_to_events:
    input:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf",
        reference_gtf = REFERENCE_GTF
    output:
        shiba_out = directory("<merged_shiba_results>/events")
    params:
        shiba_scripts = SHIBA_SCRIPTS
    priority: 1
    threads: 10 # too many? I want it to run fast
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/gtf2event.py -i {input.merged_gtf} -r {input.reference_gtf} -o {output.shiba_out} -p {threads} -v
        """

rule calculate_psi:
    input:
        events_dir = directory("<merged_shiba_results>/events"),
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed"
    output:
        shiba_out = directory("<merged_shiba_results>/psi")
    params:
        shiba_scripts = SHIBA_SCRIPTS,
        min_reads = MIN_READS
    priority: 1
    threads: 15
    shell:
        """
        python ${{CONDA_PREFIX}}/{params.shiba_scripts}/psi.py -m {params.min_reads} -p {threads} -v --onlypsi {input.merged_junctions} {input.events_dir} {output.shiba_out}
        """

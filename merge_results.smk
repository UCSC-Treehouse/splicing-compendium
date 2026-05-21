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
SHIBA_SCRIPTS = config["shiba_scripts_path"]

# the main snakefile will also need to be changed to reflect how we handle sample groups
pathvars:
    merged_shiba_results = f"results/merged_shiba/{VERSION}",
    references = "references",

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

# path to unzipped gtfs from separate shiba runs
UNZIPPED_GTFS = expand(
    "results/{group}/shiba/{sample}/annotation/assembled_annotation.gtf",
    zip,
    group=GROUPS,
    sample=SAMPLES
)

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        shiba_out = "<merged_shiba_results>/events"

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
    params:
        # compute unzipped GTF paths string to print into manifest file for stringtie merge
        gtf_manifest = lambda wildcards, input: "\n".join(
            path.removesuffix(".gz") for path in input.sample_gtfs)
    output:
        # mark unzipped gtfs as temp so they are deleted once merging is complete
        unzipped_gtfs = temp(UNZIPPED_GTFS),
        gtf_manifest_file = "<merged_shiba_results>/gtf_manifest.txt",
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 4
    shell:
        """
        # make list of all input files and print it into manifest one line at a time
        echo "{params.gtf_manifest}" > {output.gtf_manifest_file}

        # unzip input gtfs for stringtie
        gunzip -k {input.sample_gtfs}

        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output.merged_gtf} {output.gtf_manifest_file}
        """

rule gtf_to_events:
    input:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf",
        reference_gtf = f"<references>/{GENOME_ID}.annotation.gtf"
    output:
        shiba_out = directory("<merged_shiba_results>/events")
    params:
        shiba_scripts = SHIBA_SCRIPTS
    priority: 1
    threads: 10 # too many? I want it to run fast
    shell:
        """
        python ${{CONDA_PREFIX}}/{params.shiba_scripts}/gtf2event.py -i {input.merged_gtf} -r {input.reference_gtf} -o {output.shiba_out} -p {threads} -v
        """

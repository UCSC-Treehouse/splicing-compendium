# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on
# Run from project root: snakemake --snakefile merge_results.smk -j 15

import os
import pandas as pd
from datetime import datetime

configfile: "config/merge_shiba_config.yaml"

# read in configfile values
SAMPLES = pd.read_table(config["sample_sheet"])["samples"].tolist()
GROUPS = pd.read_table(config["sample_sheet"])["group"].tolist()

# the main snakefile will also need to be changed to reflect how we handle sample groups
pathvars:
    merged_shiba_results = f"results/merged_shiba/{config["version"]}"

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
        shiba_psi_out = directory("<merged_shiba_results>/psi")

rule merge_junctions:
    input: JUNCTION_BEDS
    output: "<merged_shiba_results>/merged_junctions.bed"
    # compute joined junctions string prior to passing into join script
    params:
        # use comma-separated list of inputs for R script
        junctions=lambda wildcards, input: ",".join(input)
    priority: 1
    threads: 4
    shell:
        """
        tempdir=$(mktemp -d)
        for file in {input}; do
            # Remove duplicated fields in bedfiles separated by ";" inside the tab-delimited bedfile
            # look for "value;value" inside tab-delimited fields and replace these entries with "value"
            sed -E ':a; s/(^|\t)([^\t;]*);\2(\t|$)/\1\2\3/g; ta' {input} > $tempdir/$(basename $file)
        done

        rm -rf $tempdir

        Rscript scripts/03-merge_separate_junctions.R --junctions={params.junctions} --output={output}
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

rule gtf_to_events:
    input:
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf",
        reference_gtf = config["reference_gtf"]
    output:
        shiba_out = directory("<merged_shiba_results>/events")
    params:
        shiba_scripts = config["shiba_scripts_path"]
    priority: 1
    threads: 10 # too many? I want it to run fast
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/gtf2event.py -i {input.merged_gtf} -r {input.reference_gtf} -o {output.shiba_out} -p {threads} -v
        """

rule calculate_psi:
    input:
        events_dir = "<merged_shiba_results>/events",
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed",
        merged_gtf = "<merged_shiba_results>/merged_gtf.gtf",
    output:
        shiba_psi_out = directory("<merged_shiba_results>/psi")
    params:
        shiba_scripts = config["shiba_scripts_path"],
        min_reads = config["shiba_min_reads"]
    priority: 1
    threads: 15
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/psi.py -m {params.min_reads} -p {threads} -v --onlypsi {input.merged_junctions} {input.events_dir} {output.shiba_psi_out}
        """

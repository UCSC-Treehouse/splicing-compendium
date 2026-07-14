# Snakefile for generating and merging GTFs from one sample at a time with the reference GTF, for a defined set of samples
# This snakefile uses a sample sheet containing sample accessions and bam paths to identify files to act on
# To generate test files to assess differences caused by merging the reference in multiple times,
# this workflow will run on experiment.tsvs of one sample at a time to generate a GTF.
# Then, the GTFs generated for each sample awill be merged with the reference GTF, creating one merged GTF per sample.
# Finally, the merged GTFs of all samples will be merged together with the reference again to mimic the workflow in merge_results.smk.

# Run from project root: snakemake --snakefile merge_gtf_multiple_refs.smk --profile pheonix-profile

import pandas as pd
import os
from datetime import datetime

configfile: "config/merge_gtf_test_config.yaml"

# read in configfile values
sample_table = pd.read_table(config["sample_sheet"])
SAMPLES = sample_table["sample"].tolist()
GROUPS = sample_table["group"].tolist()
REF_GTF = config["reference_gtf"]

# make a dictionary to map samples to groups so merged GTF output won't need group wildcards
sample_to_group = dict(zip(SAMPLES, GROUPS))

pathvars:
    merged_gtf_results = f"results/multiple_refs_merged_gtf_tests/{config["version"]}"


rule all:
    input:
        "<merged_gtf_results>/merged_gtf.gtf"

rule bam2gtf:
    input:
        ref_gtf = REF_GTF,
        bam = lambda wildcard: (
            os.path.join(
                "data",
                sample_to_group[wildcard.sample],
                "star-output",
                wildcard.sample,
                "Aligned.sortedByCoord.out.bam"
            )
        )
    output: temp("<merged_gtf_results>/pre-merge/{sample}.gtf")
    threads: 8
    log: "logs/{sample}_bam2gtf.log"
    shell:
        """
        stringtie -p {threads} -G {input.ref_gtf} -o {output} {input.bam} >& {log}
        """

rule first_ref_merge:
    input:
        reference_gtf = config["reference_gtf"],
        sample_gtf = "<merged_gtf_results>/pre-merge/{sample}.gtf"
    output: temp("<merged_gtf_results>/first_ref_merge/{sample}.gtf")
    threads: 8
    log: "logs/{sample}_first_ref_merge.log"
    shell:
        """
        # instantiate manifest as a temp file
        manifest=$(mktemp)

        # Add gtf path to the manifest
        echo "{input.sample_gtf}" > $manifest

        # merge gtfs with stringtie for splice analysis
        stringtie --merge -p {threads} -G {input.reference_gtf} -o {output} $manifest

        # remove temporary manifest
        rm $manifest
        """

rule merge_all_gtfs:
    input:
        reference_gtf = config["reference_gtf"],
        sample_gtfs = expand(
            "<merged_gtf_results>/first_ref_merge/{sample}.gtf",
            sample=SAMPLES
            )
    output: "<merged_gtf_results>/merged_gtf.gtf"
    priority: 1
    threads: 8
    log: "logs/merge_gtfs.log"
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

        # remove temporary manifest
        rm $manifest
        """

# Snakefile for merging GTF and junction counts bed files produced by separate shiba runs and running Shiba on the merged files
# A sample sheet of all sample accessions with separate Shiba results, and their group, is given to the config file to identify files to run the workflow on

# Run from project root: snakemake --snakefile merge_results.smk --profile pheonix-profile

import os
import pandas as pd
from datetime import datetime

configfile: "config/compendium_v1_merge_config.yaml"

# read in configfile values
sample_table = pd.read_table(config["sample_sheet"])
SAMPLES = sample_table["sample"].tolist()
GROUPS = sample_table["group"].tolist()
REF_GTF = config["reference_gtf"]

# make a dictionary to map samples to groups so merged GTF output won't need group wildcards
sample_to_group = dict(zip(SAMPLES, GROUPS))

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

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        "<merged_shiba_results>/merged_persample_psi"

rule bam2gtf:
    input:
        ref_gtf = REF_GTF,
        bam = lambda wildcards: (
            os.path.join(
                "data",
                sample_to_group[wildcards.sample],
                "star-output",
                wildcards.sample,
                "Aligned.sortedByCoord.out.bam"
            )
        )
    output: "<merged_shiba_results>/pre-merge/{sample}.gtf"
    threads: 8
    resources:
        mem_mb = 5500,
        runtime = 360
    log: "logs/{sample}_bam2gtf.log"
    shell:
        """
        stringtie -v -p {threads} -G {input.ref_gtf} -o {output} {input.bam} > {log} 2>&1
        """

rule merge_gtfs:
    input:
        reference_gtf = config["reference_gtf"],
        sample_gtfs = expand(
            "<merged_shiba_results>/pre-merge/{sample}.gtf",
            sample = SAMPLES
            )
    output: "<merged_shiba_results>/merged_gtf.gtf"
    priority: 1
    threads: 8
    resources:
        mem_mb = 200000,
        runtime = 360
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
        stringtie -v --merge -p {threads} -G {input.reference_gtf} -o {output} $manifest > {log} 2>&1

        # remove temp manifest file
        rm $manifest
        """

rule merge_junctions:
    input: JUNCTION_BEDS
    output: "<merged_shiba_results>/merged_junctions.bed"
    priority: 10
    threads: 15
    resources:
        mem_mb = 2200000,
        runtime = 720
    shell:
        """
        tempdir=$(mktemp -d)
        # initialize a number for identifying the temp junction files
        n=1

        for file in {input}; do
            # Remove duplicated fields in bedfiles separated by ";" inside the tab-delimited bedfile
            awk 'BEGIN{{FS=OFS="\t"}} {{ # set tab as delimiter
              for (i = 1; i <= NF; i++) {{
                  if ($i ~ /;/) {{ # only process fields containing ";"
                      n = split($i, parts, ";") # split into parts on semicolons
                      new = parts[1]
                      for (j = 2; j <= n; j++) {{
                          if (parts[j] != parts[j-1]) {{ # skip consecutive duplicates
                              new = new ";" parts[j] # reassemble if values are different
                          }}
                      }}
                      $i = new
                  }}
              }}
              print
          }}' $file > $tempdir/$n.bed

          ((n ++))
        done

        Rscript scripts/03-merge_separate_junctions.R --junctions=$tempdir --output={output}

        # remove tempdir of deduplicated junctions
        rm -rf $tempdir
        """


rule merge_junctions_persample:
    input: JUNCTION_BEDS
    # output is one bedfile per sample, all in a single directory
    output:
        bedfile = "<merged_shiba_results>/merged_junctions/all_junctions.bed",
        counts = expand("<merged_shiba_results>/merged_junctions/{sample}_junction_counts.tsv", sample = SAMPLES)
    params:
        output_dir = subpath(output.bedfile, parent=True)
    priority: 10
    threads: 15
    resources:
        mem_mb = 1600000,
        runtime = 400
    shell:
        """
        tempdir=$(mktemp -d)
        # initialize a number for identifying the temp junction files
        n=1

        for file in {input}; do
            # Remove duplicated fields in bedfiles separated by ";" inside the tab-delimited bedfile
            awk 'BEGIN{{FS=OFS="\t"}} {{ # set tab as delimiter
              for (i = 1; i <= NF; i++) {{
                  if ($i ~ /;/) {{ # only process fields containing ";"
                      n = split($i, parts, ";") # split into parts on semicolons
                      new = parts[1]
                      for (j = 2; j <= n; j++) {{
                          if (parts[j] != parts[j-1]) {{ # skip consecutive duplicates
                              new = new ";" parts[j] # reassemble if values are different
                          }}
                      }}
                      $i = new
                  }}
              }}
              print
          }}' $file > $tempdir/$n.bed

          ((n ++))
        done

        Rscript scripts/03-merge_separate_junctions.R --junctions=$tempdir --output_dir {params.output_dir}

        # remove tempdir of deduplicated junctions
        rm -rf $tempdir
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
    threads: 10
    resources:
        mem_mb = 60000,
        runtime = 360
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/gtf2event.py -i {input.merged_gtf} -r {input.reference_gtf} -o {output.shiba_out} -p {threads} -v
        """

rule calculate_sample_psi:
    input:
        junctions = "<merged_shiba_results>/merged_junctions/all_junctions.bed",
        junction_sample_count = "<merged_shiba_results>/merged_junctions/{sample}_junction_counts.tsv",
        events_dir = "<merged_shiba_results>/events",
    output:
        shiba_psi_out = directory("<merged_shiba_results>/sample_psi/{sample}")
    params:
        shiba_scripts = config["shiba_scripts_path"],
        min_reads = config["shiba_min_reads"]
    priority: 1
    threads: 8
    resources:
        mem_mb = 60000,
        runtime = 360
    shell:
        """
        tempfile=$(mktemp)
        trap "rm -f $tempfile" EXIT

        paste -d "\t" {input.junctions} {input.junction_sample_count} > $tempfile
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/psi.py \
            -m {params.min_reads} \
            -p {threads} \
            -v \
            --onlypsi \
            $tempfile {input.events_dir} \
            {output.shiba_psi_out}
        """

rule merge_persample_psi:
    input:
        persample_psi_dir = expand("<merged_shiba_results>/sample_psi/{sample}", sample = SAMPLES)
    output:
        directory("<merged_shiba_results>/merged_persample_psi")
    params:
        sample_sheet = config["sample_sheet"],
        version = config["version"]
    priority: 1
    threads: 15 # not currently in use
    resources:
        mem_mb = 20000,
        runtime = 480
    shell:
        """
        Rscript scripts/04-merge_sample_psi.R --sample_sheet={params.sample_sheet} --version_dir={params.version} --output_dir={output}
        """

rule calculate_merged_psi:
    input:
        events_dir = "<merged_shiba_results>/events",
        merged_junctions = "<merged_shiba_results>/merged_junctions.bed",
    output:
        shiba_psi_out = directory("<merged_shiba_results>/psi")
    params:
        shiba_scripts = config["shiba_scripts_path"],
        min_reads = config["shiba_min_reads"]
    priority: 1
    threads: 8
    resources:
        mem_mb = 2200000,
        runtime = 20160
    shell:
        """
        python ${{CONDA_PREFIX:-.}}/{params.shiba_scripts}/psi.py -m {params.min_reads} -p {threads} -v --onlypsi {input.merged_junctions} {input.events_dir} {output.shiba_psi_out}

        # zip PSI results
        pigz -p {threads} \
        {output.shiba_psi_out}/*.txt

        """

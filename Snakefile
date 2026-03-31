# snakefile for compendium shiba run
# Build reference indexes first if necessary: snakemake --snakefile build_references.smk --cores 15
# Usage: snakemake --cores 15

import os

configfile: "config.yaml"


GENOME_ID = config["genome_id"]
SAMPLE_GROUP = config["sample_group"]
if config.get("samples"):
    SAMPLES = config["samples"]
else:
    SAMPLES, = glob_wildcards(os.path.join("data", SAMPLE_GROUP, "fastq", "{sample}_1.fastq.gz"))

pathvars:
    data = os.path.join("data", SAMPLE_GROUP),
    reports = os.path.join("reports", SAMPLE_GROUP),
    results = os.path.join("results", SAMPLE_GROUP),
    logs = os.path.join("logs", SAMPLE_GROUP),
    references = "references"

# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        expand("<results>/shiba/{sample}/", sample = SAMPLES)

# first rule: trim reads with fastp
rule trim_reads:
    input:
        fastq1 = "<data>/fastq/{sample}_1.fastq.gz",
        fastq2 = "<data>/fastq/{sample}_2.fastq.gz"
    output:
        fastq1 = temp("<data>/trimmed/{sample}_1.fastq.gz"),
        fastq2 = temp("<data>/trimmed/{sample}_2.fastq.gz"),
        fastp_html = "<reports>/{sample}_fastp.html",
        fastp_json = "<reports>/{sample}_fastp.json"
    log:
        "<logs>/fastp/{sample}-fastp.log"
    threads: 7
    resources:
      mem_mb = 20000
    shell:
        """
        fastp \
            --in1 {input.fastq1} \
            --in2 {input.fastq2} \
            --out1 {output.fastq1} \
            --out2 {output.fastq2} \
            -h {output.fastp_html} \
            -j {output.fastp_json} \
            --thread {threads} \
            > {log} 2>&1
        """

# Second rule: Align reads
rule align_reads:
    input:
        fastq1 = "<data>/trimmed/{sample}_1.fastq.gz",
        fastq2 = "<data>/trimmed/{sample}_2.fastq.gz",
        index = f"<references>/star/{GENOME_ID}"
    output:
        bam = "<data>/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        sj = "<data>/star-output/{sample}/SJ.out.tab"
    log:
        "<logs>/star/{sample}-star.log"
    threads: 7
    resources:
      mem_mb = 60000
    params:
        star_dir = lambda wildcards, output: os.path.dirname(output.bam),
        samtools_memory = lambda wildcards, resources, threads: int(resources.mem_mb / threads * 0.7)  # use 70% of available memory for samtools sorting
    shell:
        """
        STAR \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --readFilesCommand zcat \
            --outFileNamePrefix "{params.star_dir}/" \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --twopassMode Basic \
            --quantMode GeneCounts \
            --bamRemoveDuplicatesType UniqueIdentical \
            > {log} 2>&1

        # sort separately to save memory
        samtools sort \
          -@ {threads} \
          -m {params.samtools_memory}M \
          -o {output.bam} \
          {params.star_dir}/Aligned.out.bam

        rm {params.star_dir}/Aligned.out.bam

        rm -rf "{params.star_dir}/_STARpass1"
        rm -rf "{params.star_dir}/_STARgenome"
      """

# Third rule: index alignments
rule index_bams:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    threads: 4
    shell:
        "samtools index -@ {threads} {input}"

# Fourth rule: run shiba on a single sample
rule run_shiba:
    input:
        bam = "<data>/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "<data>/star-output/{sample}/Aligned.sortedByCoord.out.bam.bai",
        gtf = f"<references>/{GENOME_ID}.annotation.gtf",
        config_template = "templates/shiba_config_template.yaml"
    output:
        shiba_out = directory("<results>/shiba/{sample}")
    log:
        "<logs>/shiba/{sample}-shiba.log"
    params:
        experiment_table = lambda wildcards, output: os.path.join(output.shiba_out, "experiment.tsv"),
        config_file = lambda wildcards, output: os.path.join(output.shiba_out, "shiba_config.yaml"),
        sample_group = SAMPLE_GROUP
    threads: 4
    resources:
      mem_mb = 32000
    shell:
        """
        mkdir -p {output.shiba_out}

        # Create a one sample experiment.tsv for Shiba run
        echo 'sample\tbam_path\tgroup\ttechnology' > {params.experiment_table}
        echo '{wildcards.sample}\t{input.bam}\t{params.sample_group}\tshort' >> {params.experiment_table}

        # create config file for each experiment.tsv generated
        cp {input.config_template} {params.config_file}
        echo 'workdir: {output.shiba_out}' >> {params.config_file}
        echo "experiment_table: {params.experiment_table}" >> {params.config_file}
        echo 'gtf: {input.gtf}' >> {params.config_file}

        # Run Shiba
        shiba.py -p {threads} {params.config_file} &> {log}

        # zip splicing results to save space
        pigz -p {threads} \
         {output.shiba_out}/annotation/*.gtf \
         {output.shiba_out}/events/*.txt \
         {output.shiba_out}/results/splicing/*.txt \
         {output.shiba_out}/results/expression/*.txt

        """

# snakefile for compendium shiba run
import os

# define sample paths
# a symbolic link is used to refer to the data directory /private/groups/treehouse/working-projects/celiang/bulk as "shiba-run-data"
# symlink is created with ln -s /private/groups/treehouse/working-projects/celiang/bulk shiba-run-data

# Usage example for testing one job at a time: snakemake -p -j 1

pathvars:
    data_dir = "shiba-run-data",
    reports_dir = "reports",
    results_dir = "results",
    references_dir = "references",
    logs = "logs"

# replace SAMPLE_GROUP with "target" or "gtex" depending on group of files to be preprocessed
SAMPLE_GROUP = "gtex"
GENOME_BASE = "gencode.v47.primary_assembly"
# SAMPLE, = glob_wildcards(os.path.join("shiba-run-data", SAMPLE_GROUP, "fastq", "{sample}_1.fastq.gz"))
SAMPLES = ["SRR601500"]


# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        expand("<results_dir>/{sample_group}_{sample}_shiba/", sample_group = SAMPLE_GROUP, sample = SAMPLES)

# rule to temporarily unzip a file (which will be deleted when done with it)
rule unzip_file:
    input:
        "{file}.gz"
    output:
        temp("{file}")
    localrule: True
    shell:
        "gunzip -c {input} > {output}"

rule star_index:
    input:
        genome_fasta = "<references_dir>/GRCh38.primary_assembly.genome.fa",
        genome_gtf = f"<references_dir>/{GENOME_BASE}.annotation.gtf"
    output:
        index_dir = f"<references_dir>/{GENOME_BASE}-STAR-database"
    log:
        "<logs>/star_index.log"
    threads: 16
    shell:
        """
        STAR \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.genome_fasta} \
            --sjdbGTFfile {input.genome_gtf} \
            --runMode genomeGenerate \
            --runThreadN {threads}
        """

# first rule: trim reads with fastp
rule trim_reads:
    input:
        fastq1 = "<data_dir>/{sample_group}/fastq/{sample}_1.fastq.gz",
        fastq2 = "<data_dir>/{sample_group}/fastq/{sample}_2.fastq.gz"
    output:
        fastq1 = temp("<data_dir>/{sample_group}/trimmed/{sample}_1.fastq.gz"),
        fastq2 = temp("<data_dir>/{sample_group}/trimmed/{sample}_2.fastq.gz"),
        fastp_html = "<reports_dir>/{sample_group}/{sample}_fastp.html",
        fastp_json = "<reports_dir>/{sample_group}/{sample}_fastp.json"
    log:
        "<logs>/{sample_group}/trimming-{sample}.log"
    threads: 8
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
            # output errors to snakemake logfile
            >> {log} 2>&1
        """

# Second rule: Align reads
rule align_reads:
    input:
        fastq1 = "<data_dir>/{sample_group}/trimmed/{sample}_1.fastq.gz",
        fastq2 = "<data_dir>/{sample_group}/trimmed/{sample}_2.fastq.gz",
        index = f"<references_dir>/{GENOME_BASE}-STAR-database"
    output:
        bam = "<data_dir>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        sj = "<data_dir>/{sample_group}/star-output/{sample}/SJ.out.tab"
    log:
        "<logs>/{sample_group}/star-alignment-{sample}.log"
    threads: 16
    resources:
      mem_mb=48000
    shell:
        """
        STAR \
            --readFilesCommand zcat \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --outFileNamePrefix "<data_dir>/{wildcards.sample_group}/star-output/{wildcards.sample}/" \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --quantMode GeneCounts \
            --bamRemoveDuplicatesType UniqueIdentical \
            --limitBAMsortRAM 60500000000 \
            --outBAMsortingBinsN 200 \
            >> {log} 2>&1

        rm -rf "<data_dir>/{wildcards.sample_group}/star-output/{wildcards.sample}/_STARpass1"
        rm -rf "<data_dir>/{wildcards.sample_group}/star-output/{wildcards.sample}/_STARgenome"
      """

# Third rule: index alignments
rule index_bams:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"

# Fourth rule: run shiba on a single sample
rule run_shiba:
    input:
        bam = "<data_dir>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "<data_dir>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam.bai",
        config_template = "shiba_config_template.yaml"
    output:
        experiment_file = "{sample_group}_{sample}_pilot_shiba_experiment.tsv",
        config_file = "{sample_group}_{sample}_shiba_config.yaml",
        shiba_output = directory("results/{sample_group}_{sample}_shiba/")
    log:
        "<logs>/{sample_group}/shiba-run-{sample}.log"
    threads: 1
    shell:
        """
        # Create experiment.tsv for Shiba run
        Rscript generate_experiment_file.R --sample={wildcards.sample} --bam={input.bam} --output={output.experiment_file} --group={wildcards.sample_group}

        # create config file for each experiment.tsv generated
        cp {input.config_template} {output.config_file}
        echo 'workdir:' >> {output.config_file}
        echo '  {output.shiba_output}' >> {output.config_file}
        echo 'experiment_table:' >> {output.config_file}
        echo '  {output.experiment_file}' >> {output.config_file}

        # Run Shiba
        shiba.py -p {threads} {output.config_file} &> {log}

        # zip splicing results to save space
        pigz {output.shiba_output}/*.txt
        """

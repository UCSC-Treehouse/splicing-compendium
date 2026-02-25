# snakefile for compendium shiba run
import os

# define sample paths

# Usage example for testing one job at a time: snakemake -p -j 1

pathvars:
    data = "shiba-run-data",
    reports = "reports",
    results = "results",
    references = "references",
    logs = "logs"

# replace SAMPLE_GROUP with "target" or "gtex" depending on group of files to be preprocessed
SAMPLE_GROUP = "gtex"
GENOME_BASE = "gencode.v47.primary_assembly"
# SAMPLE, = glob_wildcards(os.path.join("shiba-run-data", SAMPLE_GROUP, "fastq", "{sample}_1.fastq.gz"))
SAMPLES = ["ERR15741037"]


# create all rule with expanded wildcards because cannot run target rules with wildcards
rule all:
    input:
        expand("<results>/{sample_group}/{sample}_shiba/", sample_group = SAMPLE_GROUP, sample = SAMPLES)

# rule to temporarily unzip fasta or gtf files (which will be deleted when done)
rule unzip_refs:
    input:
        "{file}.gz"
    output:
        temp("{file}")
    wildcard_constraints:
        file = r".+\.(fa|fasta|gtf)"
    localrule: True
    shell:
        "gunzip -c {input} > {output}"

rule star_index:
    input:
        genome_fasta = "<references>/GRCh38.primary_assembly.genome.fa",
        genome_gtf = f"<references>/{GENOME_BASE}.annotation.gtf"
    output:
        index_dir = directory(f"<references>/star/{GENOME_BASE}")
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
        fastq1 = "<data>/{sample_group}/fastq/{sample}_1.fastq.gz",
        fastq2 = "<data>/{sample_group}/fastq/{sample}_2.fastq.gz"
    output:
        fastq1 = temp("<data>/{sample_group}/trimmed/{sample}_1.fastq.gz"),
        fastq2 = temp("<data>/{sample_group}/trimmed/{sample}_2.fastq.gz"),
        fastp_html = "<reports>/{sample_group}/{sample}_fastp.html",
        fastp_json = "<reports>/{sample_group}/{sample}_fastp.json"
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
        fastq1 = "<data>/{sample_group}/trimmed/{sample}_1.fastq.gz",
        fastq2 = "<data>/{sample_group}/trimmed/{sample}_2.fastq.gz",
        index = f"<references>/star/{GENOME_BASE}"
    output:
        bam = "<data>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        sj = "<data>/{sample_group}/star-output/{sample}/SJ.out.tab"
    log:
        "<logs>/{sample_group}/star-alignment-{sample}.log"
    threads: 16
    resources:
      mem_mb=48000
    shell:
        """
        star_dir=$(dirname "{output.bam}")
        STAR \
            --readFilesIn <(gunzip -c {input.fastq1}) <(gunzip -c {input.fastq2}) \
            --outFileNamePrefix "${{star_dir}}/" \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --quantMode GeneCounts \
            --bamRemoveDuplicatesType UniqueIdentical \
            --limitBAMsortRAM 60500000000 \
            --outBAMsortingBinsN 200 \
            >> {log} 2>&1

        rm -rf "${{star_dir}}/_STARpass1"
        rm -rf "${{star_dir}}/_STARgenome"
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
        bam = "<data>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "<data>/{sample_group}/star-output/{sample}/Aligned.sortedByCoord.out.bam.bai",
        config_template = "shiba_config_template.yaml"
    output:
        shiba_out = directory("<results>/{sample_group}/{sample}_shiba/")
    log:
        "<logs>/{sample_group}/shiba-run-{sample}.log"
    threads: 1
    shell:
        """
        # Create a one sample experiment.tsv for Shiba run
        experiment_file={output.shiba_out}/experiment.tsv
        echo 'sample\tbam_path\tgroup\ttechnology' > ${{experiment_file}}
        echo '{wildcards.sample}\t{input.bam}\t{wildcards.sample_group}\tshort' >> ${{experiment_file}}

        # create config file for each experiment.tsv generated
        config_file={output.shiba_out}/shiba_config.yaml
        cp {input.config_template} ${{config_file}}
        echo 'workdir: {output.shiba_out}' >> ${{config_file}}
        echo "experiment_table: ${{experiment_file}}" >> ${{config_file}}

        # Run Shiba
        shiba.py -p {threads} ${{config_file}} &> {log}

        # zip splicing results to save space
        pigz {output.shiba_out}/splicing/*.txt
        pigz {output.shiba_out}/expression/*.txt
        """

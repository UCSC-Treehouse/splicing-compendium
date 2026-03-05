# Snakemake workflow for downloading reference files and building genome indexes.
# Run from project root: snakemake --snakefile build_references.smk -j 16

configfile: "config.yaml"

pathvars:
    references = "references",
    logs = "logs"


rule all:
    input:
        index_dir = f"<references>/star/{config['genome_id']}",
        genome_gtf = f"<references>/{config['genome_id']}.annotation.gtf",
        genome_fasta = f"<references>/{config['genome_assembly']}.fa.gz"


rule download_fasta:
    output:
        f"<references>/{config['genome_assembly']}.fa.gz"
    localrule: True
    shell:
        "curl -fsSL '{config[genome_fasta_url]}' -o {output}"


rule download_gtf:
    output:
        f"<references>/{config['genome_id']}.annotation.gtf.gz"
    localrule: True
    shell:
        "curl -fsSL '{config[genome_gtf_url]}' -o {output}"


rule unzip_files:
    input:
        "{file}.gz"
    output:
        temp("{file}")
    wildcard_constraints:
        file = r".*\.(fa|fasta|gtf)$"
    localrule: True
    shell:
        "gunzip -c {input} > {output}"


rule star_index:
    input:
        genome_fasta = f"<references>/{config['genome_assembly']}.fa",
        genome_gtf = f"<references>/{config['genome_id']}.annotation.gtf"
    output:
        index_dir = directory(f"<references>/star/{config['genome_id']}")
    log:
        f"<logs>/star-index/{config['genome_id']}.log"
    threads: 16
    resources:
      mem_mb = 48000
    shell:
        """
        STAR \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.genome_fasta} \
            --sjdbGTFfile {input.genome_gtf} \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            > {log} 2>&1

        """

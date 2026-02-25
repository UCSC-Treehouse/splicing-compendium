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


rule download_genome_fasta:
    output:
        f"<references>/{config['genome_id']}.fa.gz"
    log:
        "<logs>/download_genome_fasta.log"
    localrule: True
    shell:
        "curl -fsSL '{config[genome_fasta_url]}' -o {output} >> {log} 2>&1"


rule download_genome_gtf:
    output:
        f"<references>/{config['genome_id']}.annotation.gtf.gz"
    log:
        "<logs>/download_genome_gtf.log"
    localrule: True
    shell:
        "curl -fsSL '{config[genome_gtf_url]}' -o {output} >> {log} 2>&1"


rule unzip_file:
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
        genome_fasta = f"<references>/{config['genome_id']}.fa",
        genome_gtf = f"<references>/{config['genome_id']}.annotation.gtf"
    output:
        index_dir = directory(f"<references>/star/{config['genome_id']}")
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

# Treehouse Splice Compendium repo

The code in this repository uses [Shiba (Kubota 2025)](https://github.com/Sika-Zheng-Lab/Shiba) to generate transcriptome-wide alternative splicing profiles for pediatric cancer research.
The work is described in our manuscript on biorxiv at LINK HERE.
Splice data (matrix of PSI values for each sample, accompanied by chromosome coordinates, annotation status of splice event, splice event type, gene ID, and gene name) quantified for the samples analyzed in this compendium can be downloaded at [10.5281/zenodo.22737043
](https://zenodo.org/records/22737043).
Clinical and tissue type metadata for samples present in this resource can be found at the Zenodo link, or at [metadata/combined_compendium_metadata.tsv](https://github.com/UCSC-Treehouse/splicing-compendium/blob/main/metadata/combined_compendium_metadata.tsv).

The datasets used in this project are under controlled access on dbGaP. They consist of:

- [GTEx](https://www.gtexportal.org/home/aboutAdultGtex): RNA sequences from non-diseased tissue sites collected from adult patients, ranging from 20-79 years of age.
- [TARGET](https://www.cancer.gov/ccg/research/genome-sequencing/target): RNA sequences from pediatric cancer samples spanning 6 cancer types 
(Acute Lymphoblastic Leukemia, Acute Myeloid Leukemia, Wilms Tumor, Clear Cell Sarcoma, Rhabdoid Tumor, and Neuroblastoma). 
Osteosarcoma samples are excluded for now, due to their complex transcriptomes and the relative rarity of bone RNA-seq samples.

## Setting up the environment

The environment for this project is managed with Pixi.
After cloning the repo, the environment can be set up with the following commands:

Install Pixi if it is not already installed:

```
curl -fsSL https://pixi.sh/install.sh | sh
https://pixi.prefix.dev/latest/installation/

```

Then, call the pixi environment to start downloading dependencies

```
pixi shell
```

All commands in this readme are intended to be run from within the pixi shell environment.

Note that currently the only supported platform is Linux, so this command must be run on a Linux machine.

## Filtering accessions for download

Data filtering for GTEx and TARGET accession download was performed with `notebooks/data_filtering/filter_gtex_sra.qmd` and `notebooks/data_filtering/filter_target_sra.qmd`.
Running these notebooks produce `gtex_accessions.txt` and `target_accessions.txt` files in `metadata/filter_target_gtex`, 
which are used to download the raw sequence files processed in this compendium. 
These notebooks can be run interactively within an Rstudio session, or with the following command:

```
quarto render filter_gtex_sra.qmd
```

## Data download

Sequence file download is performed with `scripts/01-fastq-download.sh`, in batches defined by the accessions files generated above.

General usage: 

```
scripts/01-fastq-download.sh [dataset] [batch number]
```

example command for downloading the first 1.2TB batch of target: 

```
scripts/01-fastq-download.sh target 1
```

## Workflow configuration

Most configuration for this project is managed with the `config.yaml` and `compendium_v1_merge_config.yaml` files within the `config/` directory.
These files contain information about the reference genome and annotation to be used for STAR indexing and alignment, as well as sample group information for runs of the workflows.

## Download reference files and generate STAR index for file processing

Downloading reference files and generating the STAR index can be done with the Snakemake workflow in `build_references.smk`:

```
snakemake --cores 16 -s build_references.smk
```

This will create a `references` directory in the project directory, which will contain the reference genome fasta and annotation gtf files, as well as the STAR index.
Other index files may be added later.

## Filtering metadata for accession numbers to download and downloading files :

`pixi run scripts/01-target-subset-download.sh [sample group]`

Replace [sample group] with the dataset to be downloaded (e.g. `target`, `gtex`)

## Running the main Snakemake workflow

<img width="2126" height="1646" alt="biorxiv_fullpage_pipeline_workflow_only" src="https://github.com/user-attachments/assets/22a09d4d-7d36-4921-bbe2-027743c71cd9" />

Input data files for the main snakemake workflow consist of raw sequence `.fastq` files, which are expected to be found in the `data` directory of this project.
Output files consist of Shiba results generated for each individual sample, and will be written to the `results` directory.
Data files are organized by the sample group (e.g. TARGET, GTEx) in the `data` directory, and the sample group is specified in the `config.yaml` file.
You may want to use symbolic links to point to the location of data files on your machine, if they are not stored in the project directory.

For example, a symbolic link is used in this workflow could be created with

```
ln -s /mnt/bulk data
```

The Snakemake workflow specified in `Snakefile` is used to trim fastq sequence files, 
align and index trimmed files, and run the Shiba splicing workflow to produce outputs for each individual sample.
Outputs of this workflow consist of:
- Aligned reads for each sample
- Per-sample Shiba output. 
Although this workflow produces all Shiba output (junction counts, GTF, PSI tables), 
only the `junctions.bed` file is used in the downstream workflow.


After updating the `config.yaml` file with the desired sample groups, reference genome information, and sample sheet, the main Snakemake workflow can be run using the default config file with:

```
snakemake --cores 16
```

To run a different sample group, you can either modify the config file or specify the sample group on the command line with:

```
snakemake --cores 16 --config sample_group={sample group}
```

By default, the workflow will run on all samples found in the `data` directory for the specified sample group.
To run on a subset of samples, you can either specify the desired samples in the `config.yaml` file or on the command line with:

```
snakemake --cores 16 --config sample_group={sample group} samples="['sample1','sample2',...]"
```

## Running the merge_results workflow

Although Shiba can be run on multiple samples using `Snakefile` to produce PSI tables, 
this method does not take advantage of the strategy Shiba employs to quantify unannotated splice events using a GTF constructed from input .bam files.
To maximize the repertoire of transcripts annotated in the GTF and therefore the possible splice events quantified by Shiba for a given sample,
`merge_results.smk` is run.

Inputs to this workflow consist of:
- Aligned, sorted, and indexed `.bam`s generated from `Snakefile`
- Per-sample junction bedfiles generated from Shiba via `Snakefile`

Outputs consist of:
- Merged GTF created from all samples' `.bam` files and the reference GTF specified by `compendium_v1_merge_config.yaml`
- Merged junction bedfile for all samples
- Shiba-generated PSI tables produced from the above two files

The workflow is run via the following command:

```
snakemake --snakefile merge_results.smk --profile pheonix-profile
```

Note that `pheonix-profile` points to a cluster-specific config file within the `pheonix-profile/` directory and may need to be reconfigured depending on the user's cluster.

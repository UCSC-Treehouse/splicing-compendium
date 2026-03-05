# Shiba run on bulk pediatric cancer RNA-seq samples

To determine the existence of interesting aberrant splicing events in pediatric cancer samples, we are interested in running Shiba on a dataset of pediatric tumor bulk RNA-seq samples.
This project currently contains scripts to filter metadata from publicly available datasets.

The datasets used in this project are under controlled access on dbGaP. They consist of:

- [GTEx](https://www.gtexportal.org/home/aboutAdultGtex): RNA sequences from 54 non-diseased tissue sites collected from adult patients, ranging from 20-79 years of age.
- [TARGET](https://www.cancer.gov/ccg/research/genome-sequencing/target): RNA sequences from pediatric cancer samples spanning 6 cancer types (Acute Lymphoblastic Leukemia, Acute Myeloid Leukemia, Wilms Tumor, Clear Cell Sarcoma, Rhabdoid Tumor, Neuroblastoma, and Osteosarcoma). Osteosarcoma samples are excluded for now, due to their complex transcriptomes.

To run these scripts, access to datasets must have been granted in dbGaP.
The repository key file (`.ngc`) must be downloaded from dbGaP and the path to that file specified in `01-target-subset-download.sh`.
Since the dataset consists of large files, the location of the user's cache may need to be changed to a directory with more storage.
This can be done by following the instructions outlined [here](https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration).

## Setting up the environment

The environment for this project is managed with Pixi and can be set up with the following command:

```
pixi install
```

Note that currently the only supported platform is Linux, so this command must be run on a Linux machine (e.g. OpenStack instance).

## Configuration

Most configuration for this project is managed with the `config.yaml` file.
This file contains information about the reference genome and annotation to be used for STAR indexing and alignment, as well as sample group information for runs of the workflow.

Note that input data files are expected to be found in the `data` directory of this project, and output files will be written to the `results` directory.
Note that data files are organized by the sample group (e.g. TARGET, GTEx) in the `data` directory, and the sample group is specified in the `config.yaml` file.
You may want to use symbolic links to point to the location of data files on your machine, if they are not stored in the project directory.

For example, a symbolic link is used in this workflow could be created with

```
ln -s /mnt/bulk data
```

## Download reference files and generate STAR index for file processing

Downloading reference files and generating the STAR index can be done with the Snakemake workflow in `build_references.smk`:

```
pixi run snakemake --cores 16 -s build_references.smk
```

This will create a `references` directory in the project directory, which will contain the reference genome fasta and annotation gtf files, as well as the STAR index.
Other index files may be added later.

## Filtering metadata for accession numbers to download and downloading files :

`pixi run scripts/01-target-subset-download.sh [sample group]`

Replace [sample group] with the dataset to be downloaded (e.g. `target`, `gtex`)


## Running the main Snakemake workflow

After updating the `config.yaml` file with the desired sample groups and reference genome information, the main Snakemake workflow can be run using the default config file with:

```
pixi run snakemake --cores 16
```

To run a different sample group, you can either modify the config file or specify the sample group on the command line with:

```
pixi run snakemake --cores 16 --config sample_group={sample group}
```

By default, the workflow will run on all samples found in the `data` directory for the specified sample group.
To run on a subset of samples, you can either specify the desired samples in the `config.yaml` file or on the command line with:

```
pixi run snakemake --cores 16 --config sample_group={sample group} samples="['sample1','sample2',...]"
```

# Transferring files between OpenStack instances

Due to limitations on allotted storage per OpenStack instance, some files need to be offloaded to a separate instance to make room for both raw sequence files and alignments.
This is managed with `transfer-and-offload-target.sh`.

Subsets of TARGET and GTEx were downloaded on openstack instance ubuntu@10.50.100.19

To make room for alignments of the sequence files, TARGET `.fastq`s were transferred to ubuntu@10.50.100.47 using the following command:

From within ubuntu@10.50.100.19:

`bash transfer-and-offload-target.sh transfer`

Then, from within ubuntu@10.50.100.47, transferred files were checked:

`bash transfer-and-offload-target.sh checksums`

Finally, from within ubuntu@10.50.100.19, the transferred files were offloaded using the checksum results:

`bash transfer-and-offload-target.sh offload /mnt/splicing-project/data/logs/2025-10-01T04:12:24_transfer_checksum.txt`

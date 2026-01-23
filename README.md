# Shiba run on bulk pediatric cancer RNA-seq samples

To determine the existence of interesting aberrant splicing events in pediatric cancer samples, we are interested in running Shiba on a dataset of pediatric tumor bulk RNA-seq samples.
This project currently contains scripts to filter metadata from publicly availible datasets.

The datasets used in this project are under controlled access on dbGaP. They consist of:

- [GTEx](https://www.gtexportal.org/home/aboutAdultGtex): RNA sequences from 54 non-diseased tissue sites collected from adult patients, ranging from 20-79 years of age.
- [TARGET](https://www.cancer.gov/ccg/research/genome-sequencing/target): RNA sequences from pediatric cancer samples spanning 6 cancer types (Acute Lymphoblastic Leukemia, Acute Myeloid Leukemia, Wilms Tumor, Clear Cell Sarcoma, Rhabdoid Tumor, Neuroblastoma, and Osteosarcoma). Osteosarcoma samples are excluded for now, due to their complex transcriptomes.

To run these scripts, access to datasets must have been granted in dbGaP.
The repository key file (`.ngc`) must be downloaded from dbGaP and the path to that file specified in `01-target-subset-download.sh`.
Since the dataset consists of large files, the location of the user's cache may need to be changed to a directory with more storage.
This can be done by following the instructions outlined [here](https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration).

## Filtering metadata for accession numbers to download and downloading files :

`pixi run scripts/01-target-subset-download.sh [sample group]`

Replace [sample group] with the dataset to be downloaded (e.g. target, gtex)

## Running Snakemake workflow for analyzing pediatric cancer samples with Shiba

Activate environment

```
pixi shell
```

Run Snakemake workflow
```
snakemake -p --cores 8
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

# Running Shiba on TARGET subset of 88 samples to obtain splice results for clustering
Eventually, this step will be incorporated into the Snakemake pipeline.
For the time being, the Shiba run on the TARGET subset was run using:

`bash 02-target-subset-shiba-run.sh`

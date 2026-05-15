### merge junctions.bed files produced by separate shiba runs ###

# define sample group
group <- "target"

## Directories and files ##
# define the data directories
# sample sheet dir
config_dir <- file.path("config")
# shiba results dir
shiba_dir <- file.path("results", group, "shiba")
# output dir for merged bedfile
out_dir <- file.path("results/merged_shiba")

# define files
# define output file
out_junctions <- file.path(out_dir, "merged_junctions.bed")

# define input files
# sample sheet of samples to merge
sample_sheet_file <- file.path(config_dir, "two_sample_merge_samples.tsv")

## make directories if they dont exist ##
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## Read in files ##
sample_sheet <- readr::read_tsv(sample_sheet_file, col_names = FALSE)
# convert df to list of samples
samples <- sample_sheet$X1

# make sample paths to the junctions.bed files for separate splice runs
junction_paths <- file.path(
  shiba_dir,
  samples,
  "junctions",
  "junctions.bed"
)
# name the separate junction file paths
names(junction_paths) <- names(samples)


## read in files and merge##
# read in junctions.bed files created from separate Shiba runs
merged_junctions <- purrr::map(junction_paths, \(file) {
    read.table(file,
               header = TRUE,
               sep="\t",
               stringsAsFactors=FALSE,
               quote="",
               # make sure columns are all the same class for merging
               colClasses = "character")

  }) |>
    # merge junctions tables from multiple samples
    # the resulting table separates junction counts from each sample by columns with the sample ID
    purrr::reduce(\(x, y) dplyr::full_join(x, y, by = c("ID", "start", "end", "chr")))

## Save merged junction counts as output
readr::write_tsv(merged_junctions, out_junctions)

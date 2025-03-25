## This script takes in a directory full of
## CZID output files, where each file corresponds to the
## analysis done on a single sample.
## This script then outputs one file per taxonomic level,
## which combines information across all samples for the given
## taxonomic level

## input format: 1 CZID file per sample, each row is a pathogen identified
## from the sample, columns contain information about that pathogen,
## including nr_count and nt_count

## output format: 1 file per taxonomic level, where each row is a pathogen,
## and each column is a sample. first row has sample names, and first
## column has the pathogen names.

library(dplyr)
library(tidyr)
library(ggplot2)

# file path to directory containing input files
inputDir <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/"
outputDir <- "./"
taxLevel <- "genus"    # can do species or genus level

if(taxLevel == "genus") {
  taxNum <- 2
  taxSym <- "g"
} else if(taxLevel == "species") {
  taxNum <- 1
  taxSym <- "s"
} else {
  print("invalid taxonomic level")
  stop()
}
outputFileNR <- paste0(outputDir, "CZID_nr_", taxLevel, "_counts.csv")
outputFileNT <- paste0(outputDir, "CZID_nt_", taxLevel, "_counts.csv")

if (dir.exists(inputDir)) {
  files <- list.files(path = inputDir)
} else {
  print("input directory does not exist")
}

# keep only the files that end in 'taxon_report.csv'
indices <- grep("taxon_report\\.csv$", files)
filtered_files <- files[indices]
file_paths <- paste0(inputDir, filtered_files)

# extract the pathogens present and nr_counts and nt_counts for each sample into a list
extract_data <- function(file) {
  sampleName <- paste0("Sample_", sub("_taxon_report.csv", "", tail(strsplit(file, "/")[[1]] , 1) ) )
  data <- read.csv(file)
  taxFilteredData <- data[data$'tax_level' == taxNum, ]
  taxFilteredData <- taxFilteredData[, c("category", "name", "nt_count", "nr_count")]
  taxFilteredData$"Sample" <- sampleName
  return(taxFilteredData)
}
extracted_data <- lapply(file_paths, extract_data)
combined_data <- do.call(rbind, extracted_data)
combined_data[is.na(combined_data)] <- 0

# if there is a duplicate pathogen-sample combo (i.e. multiple species in a genus), sum them
combined_data <- combined_data %>%
  group_by(name, Sample) %>%
  summarize(
    nr_count = sum(nr_count),
    nt_count = sum(nt_count),
    category = category[1]
  )
  
# separate into 2 dataframes with new name (combined from category and name),
# and either nr_counts or nt_counts as data, and the appropriate row and column headers
pathogenNames <- paste0("k__", combined_data$category, ";-", taxSym, "__", combined_data$name)
nr_data <- data.frame(Pathogens = pathogenNames, Sample = combined_data$Sample, count = combined_data$nr_count)
nt_data <- data.frame(Pathogens = pathogenNames, Sample = combined_data$Sample, count = combined_data$nt_count)

# condense into final format
nr_formatted <- nr_data %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)
nt_formatted <- nt_data %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)

# print 2 output files
write.csv(nr_formatted, outputFileNR, row.names = FALSE)
write.csv(nt_formatted, outputFileNT, row.names = FALSE)

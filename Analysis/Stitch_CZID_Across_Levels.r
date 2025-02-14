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
  speciesOnlyData <- data[data$'tax_level' == "1", ]
  speciesOnlyData <- speciesOnlyData[, c("category", "name", "nt_count", "nr_count")]
  speciesOnlyData$"Sample" <- sampleName
  return(speciesOnlyData)
}
extracted_data <- lapply(file_paths, extract_data)
combined_data <- do.call(rbind, extracted_data)
combined_data[is.na(combined_data)] <- 0
  
# separate into 2 dataframes with new name (combined from category and name),
# and either nr_counts or nt_counts as data, and the appropriate row and column headers
pathogenNames <- paste0("k__", combined_data$category, ";-s__", combined_data$name)
nr_data <- data.frame(species = pathogenNames, Sample = combined_data$Sample, count = combined_data$nr_count)
nt_data <- data.frame(species = pathogenNames, Sample = combined_data$Sample, count = combined_data$nt_count)

# make histograms showing counts for nr and nt to show equivalent comparisons
countsHistogram <- function(data) {
  plot <- ggplot(data, aes(x=count)) +
    geom_histogram(binwidth = 100000, fill = "blue") +
    labs(title = "Histogram of all counts from CZID",
         x = "counts of pathogen found",
         y = "number of sample-pathogen pairs with this count")
  
  return(plot)
}

nr_hist <- countsHistogram(nr_data)
nt_hist <- countsHistogram(nt_data)

# make CDF showing counts for nr and nt
countsCDF <- function( data ) {
  plot <- ggplot( data , aes(x=rowNum, y=count )) +
    geom_point() +
    labs(title = "CDF of pathogens found from Krystelle's microbiome data from CZID",
         x = "sample-pathogen combination",
         y = "counts found")
  
  return(plot)
}

nr_ordered <- nr_data[order(nr_data$count), ]
nr_ordered$rowNum <- as.numeric(rownames(nr_data))
nr_cdf_zeros <- countsCDF(nr_ordered)
nr_ordered_noZeros <- subset(nr_ordered, count !=0)
nr_cdf <- countsCDF(nr_ordered)

# condense into final format with cutoff informed by histogram
nr_formatted <- nr_data %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)
nt_formatted <- nt_data %>%
  pivot_wider(names_from = Sample,
              values_from = count,
              values_fill = 0)

# print 2 output files
write.csv(nr_formatted, "CZID_nt_species_counts.csv", row.names = FALSE)
write.csv(nt_formatted, "CZID_nr_species_counts.csv", row.names = FALSE)

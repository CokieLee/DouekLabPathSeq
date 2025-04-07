library(dplyr)
library(stringr)
library(tidyr)

## takes in path to csv file (format: rows are bacteria, columns are rownames + samples)
## outputs R dataframe in same format, simplifies and modifies names
## TODO: should maybe split this into two functions, load and name cleaner
loadCountsData <- function(countsPath, taxAbbrev) {
  df <- read.csv(countsPath) %>%
    filter(str_detect(taxID, paste0(taxAbbrev, "__\\S+") )) %>%
    mutate( taxID = sapply(taxID, (function(x) str_split_1(x, paste0(taxAbbrev, "__") )[2])) ) %>%
    mutate(across(-1, as.numeric)) %>%
    replace(is.na(.), 0) %>%
    rename_with(~ str_split(., "_", simplify = TRUE)[, 2], -1)
  
  return(df)
}

## transposes the data, maintaining the convention of no having colnames
## but no rownames (rownames are first column)
transposeDf <- function(countsDf) {
  pivot_longer(cols = -(taxID), names_to="SampleName", values_to="counts") %>%
  pivot_wider(names_from="taxID", values_from="counts")
}

## takes counts dataframe (in format output by loadCountsData())
## takes sample sheet dataframe (in format: c(col of sample names (named SampleName), col of treatment label (named Treat)))
## outputs counts dataframe, but with extra cols "Treat", containing treatment label
countsWithMetadata <- function(countsDf, sampleInfoDf) {
  df <- countsDf %>%
    transposeDF() %>%
    left_join(sampleInfoDf, by="SampleName") %>%
    select(SampleName, Treatment, everything())
    
  return(df)
}

filterMedian <- function() {
  pathseqThreshByMean <- meanPathseqCtrl + (3 * stdPathseqCtrl)
  medianPathseqCtrl <- median(as.numeric(pathseqCounts$Water))
  pathseqThreshByPercent <- quantile(pathseqCounts$Water, probs=c(0.99))[1] %>%
    unname()
  
  pathseq_filteredByPerc <- pathseqCounts %>%
    mutate(across(-1, ~ ifelse(. < pathseqThreshByPercent, 0, .)))
}

filterMean <- function() {
  meanPathseqCtrl <- mean(as.numeric(pathseqCounts$Water))
  stdPathseqCtrl <- sd(as.numeric(pathseqCounts$Water))
  
  pathseq_filteredByMean <- pathseqCounts %>%
    mutate(across(-1, ~ ifelse(. < pathseqThreshByMean, 0, .) ))
}

topOverallPathogens <- function() {
  
}

## test functions
inputPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/pathseqTaxonomicSummaries/RNA_genus_tpm.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"
taxAbbrev <- "g"

counts <- loadCountsData(inputPath, taxAbbrev)

sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(Sample.name, HIV.status)) %>%
  rename(!!!setNames(names(.), c("SampleName", "Treatment"))) %>%
  filter(Treatment != "")

countsMetaData <- countsWithMetadata(counts, sampleSheet)

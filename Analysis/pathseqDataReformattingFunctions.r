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
transposeForward <- function(countsDf) {
  df <- countsDf %>%
    pivot_longer(cols = -(taxID), names_to="SampleName", values_to="counts") %>%
    pivot_wider(names_from="taxID", values_from="counts")
  
  return(df)
}

## takes counts dataframe (in format output by loadCountsData())
## takes sample sheet dataframe (in format: c(col of sample names (named SampleName), col of treatment label (named Treat)))
## outputs counts dataframe, but with extra cols "Treat", containing treatment label
countsMergeMetadata <- function(countsDf, sampleInfoDf) {
  df <- countsDf %>%
    pivot_longer(cols = -(taxID), names_to="SampleName", values_to="counts") %>%
    pivot_wider(names_from="taxID", values_from="counts") %>%
    left_join(sampleInfoDf, by="SampleName") %>%
    select(SampleName, Treatment, everything())
  
  return(df)
}

countsUnmergeMetaData <- function(countsDf) {
  df <- countsDf %>%
    select(-Treatment) %>%
    pivot_longer(cols = -(SampleName), names_to="taxID", values_to="counts") %>%
    pivot_wider(names_from="SampleName", values_from="counts")
  
  return(df)
}

filterMedian <- function(countsDf, negControlCol, percentThreshold = 0.99) {
  pathseqThreshByPercent <- quantile(countsDf[[negControlCol]], probs=c(percentThreshold))[1] %>%
    unname()
  
  pathseq_filteredByPerc <- countsDf %>%
    mutate(across(-1, ~ ifelse(. < pathseqThreshByPercent, 0, .)))
  
  return(pathseq_filteredByPerc)
}

filterMean <- function(countsDf, negControlCol, stdThreshold = 3) {
  meanCtrl <- mean(as.numeric(countsDf[[negControlCol]] ))
  stdCtrl <- sd(as.numeric(countsDf[[negControlCol]] ))
  threshByMean <- (stdCtrl * stdThreshold) + meanCtrl
  
  filteredByMean <- countsDf %>%
    mutate(across(-1, ~ ifelse(. < threshByMean, 0, .) ))
  
  return(filteredByMean)
}

topOverallPathogens <- function() {
  
}

poolAllCounts <- function(countsDf, removeZerosNas, keepTopN) {
  df <- countsDf %>%
    pivot_longer(cols = !taxID,
                 names_to = "SampleName",
                 values_to = "count") %>%
    mutate(count = as.numeric(count)) %>%
    arrange(count)
  
  if (!missing(removeZerosNas) && (removeZerosNas == TRUE)) {
    df <- df %>%
      filter(count != 0)
  }
  if (!missing(keepTopN)) {
    df <- head(df, keepTopN)
  }
  
  df <- mutate(df, rowNum = as.numeric(rownames(df)) )
  
  return(df)
}

poolCountsAcrossSamples <- function(countsDf, removeZeros, keepTopN) {
  df <- countsDf %>%
    pivot_longer(cols = !taxID,
                 names_to = "SampleName",
                 values_to = "count") %>%
    mutate(count = as.numeric(count)) %>%
    group_by(taxID) %>%
    summarize(SummedCount = sum(count, na.rm = TRUE)) %>%
    arrange(desc(SummedCount))
  
  if (!missing(removeZeros) && (removeZeros == TRUE)) {
    df <- df %>%
      filter(SummedCount != 0)
  }
  if (!missing(keepTopN)) {
    df <- head(df, keepTopN)
  }
  
  return(df)
}

## test functions
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/pathseqTaxonomicSummaries/RNA_genus_tpm.csv"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/summarizedData/CZID_nr_genus_counts.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"
taxAbbrev <- "g"

counts <- loadCountsData(pathseqPath, taxAbbrev)
counts1 <- loadCountsData(czidPath, taxAbbrev)

sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(Sample.name, HIV.status)) %>%
  rename(!!!setNames(names(.), c("SampleName", "Treatment"))) %>%
  filter(Treatment != "")

countsMetaData <- countsMergeMetadata(counts, sampleSheet)

countTransposedBack <- countsUnmergeMetaData(countsMetaData)

transposedDf <- transposeForward(counts)

pooledDf <- poolCountsAcrossSamples(counts, TRUE, 10)

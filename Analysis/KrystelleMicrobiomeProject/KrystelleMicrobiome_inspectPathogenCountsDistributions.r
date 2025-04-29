library("tidyr")
library("dplyr")
library("stringr")
library("ggplot2")
source("/data/home/parkercol/PathSeq/Analysis/pathseqDataReformattingFunctions.r")
source("/data/home/parkercol/PathSeq/Analysis/pathseqAnalysisFunctions.r")

## Define variable paths
# import csv file with genus counts
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/pathseqTaxonomicSummaries/RNA_genus_tpm.csv"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/summarizedData/CZID_nr_genus_counts.csv"
# define paths for writing filtered files
pathseqFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/pathseq_genus_tpm_filtered.csv"
czidFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/summarizedData/czid_nr_genus_counts_filtered.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"

taxLevel = "genus"
taxAbbrev = 'g'
###############################################################################
## Define data

# define pathseq all
pathseqCounts <- loadCountsData(pathseqPath, taxAbbrev)

# define czid all
# TODO: czid all currently has duplicate pathogen names (combined here). investigate why.
czidCounts <- (loadCountsData(czidPath, taxAbbrev)) %>%
  group_by(taxID) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup()

# define sample sheet with treatment information
sampleSheet <-
  read.csv(sampleSheetPath, sep="\t")

treatmentSheet <-
  sampleSheet %>%
  subset(select=c(Sample.name, HIV.status)) %>%
  rename(!!!setNames(names(.), c("SampleName", "Treatment"))) %>%
  filter(Treatment != "")

# TODO: finish fixing filtering functions and calling
# define datasets filtered by everything below 3std above negative control mean
pathseqFilteredByStd <- filterMean(pathseqCounts, "Water", 3)
czidFilteredByStd <- filterMean(czidCounts, "Water", 3)
# define datasets filtered of everything below 99th percentile of negative control
pathseqFilteredByPerc <- filterMedian(pathseqCounts, "Water", percentThreshold = 0.99)
czidFilteredByPerc <- filterMedian(czidCounts, "Water", 0.99)

# define pathseq merged with treatment information
pathseqMergeTreat <- countsMergeMetadata(pathseqCounts, treatmentSheet)
czidMergeTreat <- countsMergeMetadata(czidCounts, treatmentSheet)

## define data groups by treatment labels
pathseqHIVPos <- pathseqMergeTreat %>% filter(Treatment == "Positive") %>%
  select(-Treatment)
pathseqHIVNeg <- pathseqMergeTreat %>% filter(Treatment == "Negative") %>%
  select(-Treatment)
pathseqHUU <- pathseqMergeTreat %>% filter(Treatment == "HUU") %>%
  select(-Treatment)
pathseqHEU <- pathseqMergeTreat %>% filter(Treatment == "HEU") %>%
  select(-Treatment)
pathseqAdult <- pathseqMergeTreat %>%
  filter(Treatment == "Positive" | Treatment == "Negative")
pathseqChild <- pathseqMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
pathseqChildMeanFiltered <- pathseqMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
pathseqChildMedianFiltered <- pathseqMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")

czidHIVPos <- czidMergeTreat %>% filter(Treatment == "Positive") %>%
  select(-Treatment)
czidHIVNeg <- czidMergeTreat %>% filter(Treatment == "Negative") %>%
  select(-Treatment)
czidHUU <- czidMergeTreat %>% filter(Treatment == "HUU") %>%
  select(-Treatment)
czidHEU <- czidMergeTreat %>% filter(Treatment == "HEU") %>%
  select(-Treatment)
czidAdult <- czidMergeTreat %>%
  filter(Treatment == "Positive" | Treatment == "Negative")
czidChild <- czidMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
czidChildMeanFiltered <- czidMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
czidChildMedianFiltered <- czidMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")

################################################################################
## look at counts distribution and compare CZID to pathseq, filtering to non-filtering

## get overall pathseq and czid results distribution before filtering
pathseqPooled <- poolAllCounts(pathseqCounts)
czidPooled <- poolAllCounts(czidCounts)
cdf_zeros <-
  basePlot("Distribution of pipeline results (All results, including 0 values for \"not found\")",
           "All pipeline results (one dot for each microbial genus found across all samples)", "Amount of microbe found (tpm)") +
  geom_point(data = pathseqPooled, aes(x=rowNum, y=count), color='green') +
  geom_point(data = czidPooled, aes(x=rowNum, y=count), color='orange') +
  scale_x_continuous(limits = c(0, 750000))
cdf_zeros

pathseqPooledPositive <- poolAllCounts(pathseqCounts, removeZerosNas = TRUE)
czidPooledPositive <- poolAllCounts(czidCounts, removeZerosNas = TRUE)
cdf <-
  basePlot("Distribution of pipeline results (Excluding 0 values for \"not found\")",
           "All pipeline results (one dot for each microbial genus found across all samples)", "Amount of microbe found (tpm)") +
  geom_point(data = pathseqPooledPositive, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czidPooledPositive, aes(x=rowNum, y=count), color='orange') +
  scale_x_continuous(limits = c(0, 750000))
cdf

# get overall pathseq and czid results distribution after mean filtering
pathseqPooledMeanFiltered <- poolAllCounts(pathseqFilteredByStd, removeZerosNas = TRUE)
czidPooledMeanFiltered <- poolAllCounts(czidFilteredByStd, removeZerosNas = TRUE)
cdf_meanFiltered <-
  basePlot("Distribution of pipeline results (Excluding all values below 3 std above negative control's mean)",
           "All pipeline results (one dot for each microbial genus found across all samples)", "Amount of microbe found (tpm)") +
  geom_point(data = pathseqPooledMeanFiltered, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czidPooledMeanFiltered, aes(x=rowNum, y=count), color='orange') +
  scale_x_continuous(limits = c(0, 7000))
cdf_meanFiltered

# get overall pathseq and czid results distribution after median filtering
pathseqMedianFiltered <- poolAllCounts(pathseqFilteredByPerc, removeZerosNas = TRUE)
czidMedianFiltered <- poolAllCounts(czidFilteredByPerc, removeZerosNas = TRUE)
cdf_medianFiltered <-
  basePlot("Distribution of pipeline results (Excluding all values below 99th percentile of negative control values)",
           "All pipeline results (one dot for each microbial genus found across all samples)", "Amount of microbe found (tpm)") +
  geom_point(data = pathseqMedianFiltered, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czidMedianFiltered, aes(x=rowNum, y=count), color='orange') +
  scale_x_continuous(limits = c(0, 7000))
cdf_medianFiltered

################################################################################
## Examine the overlap between adults and children
pathseqAdultPooled <- poolCountsAcrossSamples(countsUnmergeMetaData(pathseqAdult), removeZeros = TRUE)
pathseqChildPooled <- poolCountsAcrossSamples(countsUnmergeMetaData(pathseqChild), removeZeros = TRUE)

pathseqAdultTop10Microbes <- poolCountsAcrossSamples(countsUnmergeMetaData(pathseqAdult), removeZeros = TRUE, 10)
pathseqChildTop10Microbes <- poolCountsAcrossSamples(countsUnmergeMetaData(pathseqChild), removeZeros = TRUE, 10)

# make dataframe where for every overlapping bacteria between child and adult,
# have a column of child counts and adult counts (both summed across samples)
adultChildOverallOverlap <-
  inner_join(pathseqAdultPooled, pathseqChildPooled, by="taxID")
adultChildTop10Overlap <-
  inner_join(pathseqAdultTop10Microbes, pathseqChildTop10Microbes, by="taxID")

## function 
adultChildPairwiseOverlap <- function(inputCountsDf) {
  inputCountsDf %>%
  transposeForward() %>%
  filter(SampleName != "Positive" & SampleName != "Water") %>%
  mutate(Group = substr(SampleName, start=3, stop=5)) %>%
  select(SampleName, Group, everything()) %>%
  group_by(Group) %>%
  summarize(
    across(
      -c("SampleName"),
      ~ if (all(. != 0)) list(1) else list(0),
      .names="overlapping_{.col}"),
    .groups = "keep") %>%
  rowwise() %>%
  mutate(NumOverlapped = sum(c_across(everything()) == 1) ) %>%
  ungroup() %>%
  select(Group, NumOverlapped, everything())
}

## TODO: remove
czidadultChildPairwiseOverlap <-
  czidCounts %>%
    transposeForward() %>%
    filter(SampleName != "Positive" & SampleName != "Water") %>%
    mutate(Group = substr(SampleName, start=3, stop=5)) %>%
    select(SampleName, Group, everything()) %>%
    group_by(Group) %>%
    summarize(
      across(
        -c("SampleName"),
        ~ if (all(. != 0)) list(1) else list(0),
        .names="overlapping_{.col}"),
      .groups = "keep") %>%
    rowwise() %>%
    mutate(NumOverlapped = sum(c_across(everything()) == 1) ) %>%
    ungroup() %>%
    select(Group, NumOverlapped, everything())

## function
childMicrobeNum <- function(childCountsDf, overlapDf) {
  childCountsDf %>%
  select(-c(Treatment)) %>%
  rowwise() %>%
  mutate(SumBacteriaInChild = sum(c_across(-SampleName) != 0), na.rm = TRUE ) %>%
  ungroup() %>%
  mutate(Group = substr(SampleName, start=3, stop=5)) %>%
  left_join(select(overlapDf, c(Group, NumOverlapped)), by="Group") %>%
  mutate(PercentInherited = (NumOverlapped / SumBacteriaInChild)) %>%
  select(SampleName, Group, SumBacteriaInChild, NumOverlapped, PercentInherited, everything()) %>%
  arrange(PercentInherited)
}

## graphing function
childInheritenceCDF <- function(pathseqInheritance, czidInheritance, title) {
  basePlot(title,
            "Child", "Proportion of microbes that also existed in mother") +
  geom_point(data = mutate(pathseqInheritance, RowNum = as.numeric(row.names(pathseqInheritance))),
             aes(x=RowNum, y=PercentInherited), color='green') +
  geom_point(data = mutate(czidInheritance, RowNum = as.numeric(row.names(czidInheritance))),
             aes(x=RowNum, y=PercentInherited), color='orange')
}

## get mother-child overlap percentage (of total child microbes) for each child, unfiltered
# get a dataframe of how many microbes have non-zero values for both the mother and child in each pair
pathseqAdultChildPairwise <- adultChildPairwiseOverlap(pathseqCounts)
czidAdultChildPairwise <- adultChildPairwiseOverlap(czidCounts)
# get the number of non-zero pathogens for each child
pathseqChildMicrobeNum <- childMicrobeNum(pathseqChild, pathseqAdultChildPairwise)
czidChildMicrobeNum <- childMicrobeNum(czidChild, czidAdultChildPairwise)
# plot the percentage of microbes found in each child, which were also found in child's mother
childInheritenceUnfiltered <- childInheritenceCDF(pathseqChildMicrobeNum, czidChildMicrobeNum,
                                                  "Proportion of microbes in each child inherited from mother")
childInheritenceUnfiltered

## get mother-child overlap percentage (of total child microbes) for each child, filtered by mean
# get dataframes showing how many microbes have non-zero values for both the mother and child in each pair
pathseqPairwiseMeanFiltered <- adultChildPairwiseOverlap(pathseqFilteredByStd)
czidPairwiseMeanFiltered <- adultChildPairwiseOverlap(czidFilteredByStd)
# get the number of non-zero pathogens for each mean filtered child
pathseqChildMicrobeFilterMean <- childMicrobeNum(pathseqChildMeanFiltered, pathseqPairwiseMeanFiltered)
czidChildMicrobeFilterMean <- childMicrobeNum(czidChildMeanFiltered, czidPairwiseMeanFiltered)
# plot the percentage of microbes found in each child, which are also found in child's mother
childInheritenceMeanFilt <- childInheritenceCDF(pathseqChildMicrobeFilterMean, czidChildMicrobeFilterMean,
                                                "Proportion of microbes in each child inherited from mother\n(filtered vals under 3std above mean of negative controls)")
childInheritenceMeanFilt

## get mother-child overlap percentage (of total child microbes) for each child, filtered by mean
# get dataframes showing how many microbes have non-zero values for both the mother and child in each pair
pathseqPairwiseMedianFiltered <- adultChildPairwiseOverlap(pathseqFilteredByPerc)
czidPairwiseMedianFiltered <- adultChildPairwiseOverlap(czidFilteredByPerc)
# get the number of non-zero pathogens for each mean filtered child
pathseqChildMicrobeFilterMedian <- childMicrobeNum(pathseqChildMedianFiltered, pathseqPairwiseMedianFiltered)
czidChildMicrobeFilterMedian <- childMicrobeNum(czidChildMedianFiltered, czidPairwiseMedianFiltered)
# plot the percentage of microbes found in each child, which are also found in child's mother
childInheritenceMedianFilt <- childInheritenceCDF(pathseqChildMicrobeFilterMedian, czidChildMicrobeFilterMedian,
                                                "Proportion of microbes in each child inherited from mother\n(filtered under 99th percentile of negative controls)")
childInheritenceMedianFilt

################################################################################
## Examine the distributions of particular microbes

compareMicrobeCDFByTreat <- function(treat1Counts, treat2Counts, microbeName, title, labelA, labelB) {
  # Counts for treatment 1
  treat1Microbe <- treat1Counts %>%
    select(all_of(c("SampleName", microbeName))) %>%
    rename(count = microbeName) %>%
    mutate(count = as.numeric(count)) %>%
    arrange(count) %>%
    mutate(rowNum = sort(as.numeric(rownames(.)))) %>%
    mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) ) %>%
    mutate(percentile = scaledPercentile(rowNum)) %>%
    mutate(label = labelA)

  # Counts for treatment 2
  treat2Microbe <- treat2Counts %>%
    select(all_of(c("SampleName", microbeName))) %>%
    rename(count = microbeName) %>%
    mutate(count = as.numeric(count)) %>%
    arrange(count) %>%
    mutate(rowNum = sort(as.numeric(rownames(.)))) %>%
    mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) ) %>%
    mutate(percentile = scaledPercentile(rowNum)) %>%
    mutate(label = labelB)

  combinedDf <- rbind(treat1Microbe, treat2Microbe)

  CDF <-
    ggplot() +
    geom_point(data = combinedDf, aes(x=percentile, y=count, color=label)) +
    labs(x = "Each dot is an individual (arranged by percentile of microbe content)", y = "Transcripts per million found") +
    scale_color_manual(values = setNames( c("red", "blue"), c(as.character(labelA), as.character(labelB))) ) +
    labs(title = title) +
    theme_bw()
  return(CDF)
}

# define treatment labels for adults and children (for display on graphs, and for marking in dataframe)

## Examine Gardnerella (in adults) results from pathseq
pathseqAdultCDF_Gardnerella <-
  compareMicrobeCDFByTreat(pathseqHIVPos, pathseqHIVNeg, "Gardnerella",
                           "Distribution of Gardnerella reported in adults by Pathseq", "HIV Positive", "HIV Negative")
pathseqAdultCDF_Gardnerella
# Gardnerella HIV Pos vs Neg in log scale
pathseqAdultCDF_Gardnerella_log <- pathseqAdultCDF_Gardnerella + scale_y_log10() +
  labs(title = "Distribution of Gardnerella reported in adults by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqAdultCDF_Gardnerella_log
## Examine Gardnerella (in adults) results from pathseq
czidAdultCDF_Gardnerella <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Gardnerella",
                           "Distribution of Gardnerella reported in adults by CZID", "HIV Positive", "HIV negative")
czidAdultCDF_Gardnerella
# Gardnerella HIV Pos vs Neg in log scale
czidAdultCDF_Gardnerella_log <- czidAdultCDF_Gardnerella + scale_y_log10() +
  labs(title = "Distribution of Gardnerella reported in adults by CZID", y="Transcripts per million found (LOG SCALE**)")
czidAdultCDF_Gardnerella_log

## Examine Lactobacillus (in adults) results from pathseq
pathseqAdultCDF_Lactobacillus <-
  compareMicrobeCDFByTreat(pathseqHIVPos, pathseqHIVNeg, "Lactobacillus",
                           "Distribution of Lactobacillus reported in adults by Pathseq", "HIV Positive", "HIV Negative")
pathseqAdultCDF_Lactobacillus
# Lactobacillus HIV Pos vs Neg in log scale
pathseqAdultCDF_Lactobacillus_log <- pathseqAdultCDF_Lactobacillus + scale_y_log10() +
  labs(title = "Distribution of Lactobacillus reported in adults by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqAdultCDF_Lactobacillus_log
## Examine Lactobacillus (in adults) results from pathseq
czidAdultCDF_Lactobacillus <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Lactobacillus",
                           "Distribution of Lactobacillus reported in adults by CZID", "HIV Positive", "HIV negative")
czidAdultCDF_Lactobacillus
# Lactobacillus HIV Pos vs Neg in log scale
czidAdultCDF_Lactobacillus_log <- czidAdultCDF_Lactobacillus + scale_y_log10() +
  labs(title = "Distribution of Lactobacillus reported in adults by CZID", y="Transcripts per million found (LOG SCALE**)")
czidAdultCDF_Lactobacillus_log

## Examine streptococcus (in adults) results from pathseq
pathseqAdultCDF_Streptococcus <-
  compareMicrobeCDFByTreat(pathseqHIVPos, pathseqHIVNeg, "Streptococcus",
                           "Distribution of Streptococcus reported in adults by Pathseq", "HIV Positive", "HIV Negative")
pathseqAdultCDF_Streptococcus
# Streptococcus HIV Pos vs Neg in log scale
pathseqAdultCDF_Streptococcus_log <- pathseqAdultCDF_Streptococcus + scale_y_log10() +
  labs(title = "Distribution of Streptococcus reported in adults by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqAdultCDF_Streptococcus_log
## Examine Streptococcus (in adults) results from pathseq
czidAdultCDF_Streptococcus <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Streptococcus",
                           "Distribution of Streptococcus reported in adults by CZID", "HIV Positive", "HIV negative")
czidAdultCDF_Streptococcus
# Streptococcus HIV Pos vs Neg in log scale
czidAdultCDF_Streptococcus_log <- czidAdultCDF_Streptococcus + scale_y_log10() +
  labs(title = "Distribution of Streptococcus reported in adults by CZID", y="Transcripts per million found (LOG SCALE**)")
czidAdultCDF_Streptococcus_log
##############################################################################
## Chilren

## Examine Lactobacillus (in children) results from pathseq
pathseqChildCDF_Lactobacillus <-
  compareMicrobeCDFByTreat(pathseqHEU, pathseqHUU, "Lactobacillus",
                           "Distribution of Lactobacillus reported in babies by Pathseq", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
pathseqChildCDF_Lactobacillus
# Lactobacillus HIV Pos vs Neg in log scale
pathseqChildCDF_Lactobacillus_log <- pathseqChildCDF_Lactobacillus + scale_y_log10() +
  labs(title = "Distribution of Lactobacillus reported in babies by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqChildCDF_Lactobacillus_log
## Examine Lactobacillus (in Childs) results from pathseq
czidChildCDF_Lactobacillus <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Lactobacillus",
                           "Distribution of Lactobacillus reported in babies by CZID", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
czidChildCDF_Lactobacillus
# Lactobacillus HIV Pos vs Neg in log scale
czidChildCDF_Lactobacillus_log <- czidChildCDF_Lactobacillus + scale_y_log10() +
  labs(title = "Distribution of Lactobacillus reported in babies by CZID", y="Transcripts per million found (LOG SCALE**)")
czidChildCDF_Lactobacillus_log

## Examine Bifidobacterium (in children) results from pathseq
pathseqChildCDF_Bifidobacterium <-
  compareMicrobeCDFByTreat(pathseqHEU, pathseqHUU, "Bifidobacterium",
                           "Distribution of Bifidobacterium reported in babies by Pathseq", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
pathseqChildCDF_Bifidobacterium
# Bifidobacterium HIV Pos vs Neg in log scale
pathseqChildCDF_Bifidobacterium_log <- pathseqChildCDF_Bifidobacterium + scale_y_log10() +
  labs(title = "Distribution of Bifidobacterium reported in babies by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqChildCDF_Bifidobacterium_log
## Examine Bifidobacterium (in Childs) results from pathseq
czidChildCDF_Bifidobacterium <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Bifidobacterium",
                           "Distribution of Bifidobacterium reported in babies by CZID", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
czidChildCDF_Bifidobacterium
# Bifidobacterium HIV Pos vs Neg in log scale
czidChildCDF_Bifidobacterium_log <- czidChildCDF_Bifidobacterium + scale_y_log10() +
  labs(title = "Distribution of Bifidobacterium reported in babies by CZID", y="Transcripts per million found (LOG SCALE**)")
czidChildCDF_Bifidobacterium_log

## Examine blautia (in children) results from pathseq
pathseqChildCDF_Blautia <-
  compareMicrobeCDFByTreat(pathseqHEU, pathseqHUU, "Blautia",
                           "Distribution of Blautia reported in babies by Pathseq", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
pathseqChildCDF_Blautia
# Blautia HIV Pos vs Neg in log scale
pathseqChildCDF_Blautia_log <- pathseqChildCDF_Blautia + scale_y_log10() +
  labs(title = "Distribution of Blautia reported in babies by Pathseq", y="Transcripts per million found (LOG SCALE**)")
pathseqChildCDF_Blautia_log
## Examine Blautia (in Childs) results from pathseq
czidChildCDF_Blautia <-
  compareMicrobeCDFByTreat(czidHIVPos, czidHIVNeg, "Blautia",
                           "Distribution of Blautia reported in babies by CZID", "HIV Exposed\n(Uninfected, mothers on ART)", "HIV Unexposed")
czidChildCDF_Blautia
# Blautia HIV Pos vs Neg in log scale
czidChildCDF_Blautia_log <- czidChildCDF_Blautia + scale_y_log10() +
  labs(title = "Distribution of Blautia reported in babies by CZID", y="Transcripts per million found (LOG SCALE**)")
czidChildCDF_Blautia_log

##################################################################################
## write filtered data back to file
write.csv(pathseq_filteredByMean, pathseqFilteredPath, row.names = FALSE)
write.csv(czid_filteredByMean, czidFilteredPath, row.names = FALSE)

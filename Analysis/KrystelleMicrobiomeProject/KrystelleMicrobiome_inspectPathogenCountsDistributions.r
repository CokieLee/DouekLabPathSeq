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
pathseqAdult <- pathseqMergeTreat %>%
  filter(Treatment == "Positive" | Treatment == "Negative")
pathseqChild <- pathseqMergeTreat %>%
  filter(Treatment == "HEU" | Treatment == "HUU")

czidHIVPos <- czidMergeTreat %>% filter(Treatment == "Positive") %>%
  select(-Treatment)
czidHIVNeg <- czidMergeTreat %>% filter(Treatment == "Negative") %>%
  select(-Treatment)
czidAdult <- czidMergeTreat %>%
  filter(Treatment == "Positive" | Treatment == "Negative")
czidChild <- czidMergeTreat %>%
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
  geom_point(data = czidPooledMeanFiltered, aes(x=rowNum, y=count), color='orange')
cdf_meanFiltered

# get overall pathseq and czid results distribution after median filtering
pathseqMedianFiltered <- poolAllCounts(pathseqFilteredByPerc, removeZerosNas = TRUE)
czidMedianFiltered <- poolAllCounts(czidFilteredByPerc, removeZerosNas = TRUE)
cdf_medianFiltered <-
  basePlot("Distribution of pipeline results (Excluding all values below 99th percentile of negative control values)",
           "All pipeline results (one dot for each microbial genus found across all samples)", "Amount of microbe found (tpm)") +
  geom_point(data = pathseqMedianFiltered, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czidMedianFiltered, aes(x=rowNum, y=count), color='orange')
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

pathseqAdultChildPairwise <- adultChildPairwiseOverlap(pathseqCounts)
czidAdultChildPairwise <- adultChildPairwiseOverlap(czidCounts)

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

pathseqChildMicrobeNum <- childMicrobeNum(pathseqChild)
czidChildMicrobeNum <- childMicrobeNum(czidChild)

childInheritenceCDF <- function(pathseqInheritance, czidInheritance) {
  basePlot("Percentage of childrens' microbes existing in their mother",
            "Child", "Percentage of microbes that also existed in mother") +
  geom_point(data = mutate(pathseqInheritance, RowNum = as.numeric(row.names(pathseqInheritance))),
             aes(x=RowNum, y=PercentInherited), color='green') +
  geom_point(data = mutate(czidInheritance, RowNum = as.numeric(row.names(czidInheritance))),
             aes(x=RowNum, y=PercentInherited), color='orange')
}
childInheritenceUnfiltered <- childInheritenceCDF(pathseqChildMicrobeNum, czidChildMicrobeNum)

pathseqChildFilterMean <- countsMergeMetadata(pathseqFilteredByStd, treatmentSheet) %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
czidChildFilterMean <- countsMergeMetadata(czidFilteredByStd, treatmentSheet) %>%
  filter(Treatment == "HEU" | Treatment == "HUU")
pathseqChildMicrobeFilterMean <- childMicrobeNum(pathseqChildFilterMean)


################################################################################
## Examine the distributions of particular microbes

# pathseq plot HIV positive Gardnerella CDF
pathseqHIVPosGardnerella <- pathseqHIVPos %>%
  select(c(SampleName, `k__Bacteria;-p__Actinobacteria;-c__Actinomycetia;-o__Bifidobacteriales;-f__Bifidobacteriaceae;-g__Gardnerella`))
colnames(pathseqHIVPosGardnerella) <- c("SampleName", "count")
pathseqHIVPosGardnerella$count <- as.numeric(pathseqHIVPosGardnerella$count)
pathseqHIVPosGardnerella <- pathseqHIVPosGardnerella[order(pathseqHIVPosGardnerella$count), ]
pathseqHIVPosGardnerella$rowNum <- sort( as.numeric( rownames(pathseqHIVPosGardnerella) ) )

# pathseq plot HIV negative Gardnerella CDF
pathseqHIVNegGardnerella <- pathseqHIVNeg %>%
  select(c(SampleName, `k__Bacteria;-p__Actinobacteria;-c__Actinomycetia;-o__Bifidobacteriales;-f__Bifidobacteriaceae;-g__Gardnerella`))
colnames(pathseqHIVNegGardnerella) <- c("SampleName", "count")
pathseqHIVNegGardnerella$count <- as.numeric(pathseqHIVNegGardnerella$count)
pathseqHIVNegGardnerella <- pathseqHIVNegGardnerella[order(pathseqHIVNegGardnerella$count), ]
pathseqHIVNegGardnerella$rowNum <- sort( as.numeric( rownames(pathseqHIVNegGardnerella) ) )

# pathseq Gardnerella HIV Pos + Neg overlay
pathseqHIVPosGardnerella_format <- pathseqHIVPosGardnerella %>%
  mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) )
pathseqHIVNegGardnerella_format <- pathseqHIVNegGardnerella %>%
  mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) )
pathseqHIVAdult_CDF <-
  ggplot() +
  geom_point(data = cbind(pathseqHIVPosGardnerella_format, percentile = scaledPercentile(pathseqHIVPosGardnerella$rowNum)),
             aes(x=percentile, y=count, color="HIVpositive")) +
  geom_point(data = cbind(pathseqHIVNegGardnerella_format, percentile = scaledPercentile(pathseqHIVNegGardnerella$rowNum)),
             aes(x=percentile, y=count, color="HIVnegative")) +
  labs(x = "samples (mothers, by percentile)", y = "Pathseq bacterial counts (tpm)") +
  scale_color_manual(values = c("HIVpositive" = "red", "HIVnegative" = "blue")) +
  labs(title = "Distribution of Gardnerella found by PathSeq") +
  theme_bw()
pathseqHIVAdult_CDF
# Gardnerella HIV Pos + Neg log overlay
pathseqHIVAdultLog_CDF <- pathseqHIVAdult_CDF + scale_y_log10() +
  labs(title = "(Log Scale) Distribution of Gardnerella by Pathseq")
pathseqHIVAdultLog_CDF

# czid plot HIV positive Gardnerella CDF
czidHIVPosGardnerella <- czidHIVPos %>%
  select(c(SampleName, `k__bacteria;-g__Gardnerella`))
colnames(czidHIVPosGardnerella) <- c("SampleName", "count")
czidHIVPosGardnerella$count <- as.numeric(czidHIVPosGardnerella$count)
czidHIVPosGardnerella <- czidHIVPosGardnerella[order(czidHIVPosGardnerella$count), ]
czidHIVPosGardnerella$rowNum <- sort( as.numeric( rownames(czidHIVPosGardnerella) ) )

# czid plot HIV negative Gardnerella CDF
czidHIVNegGardnerella <- czidHIVNeg %>%
  select(c(SampleName, `k__bacteria;-g__Gardnerella`))
colnames(czidHIVNegGardnerella) <- c("SampleName", "count")
czidHIVNegGardnerella$count <- as.numeric(czidHIVNegGardnerella$count)
czidHIVNegGardnerella <- czidHIVNegGardnerella[order(czidHIVNegGardnerella$count), ]
czidHIVNegGardnerella$rowNum <- sort( as.numeric( rownames(czidHIVNegGardnerella) ) )

# czid Gardnerella HIV Pos + Neg overlay
czidHIVPosGardnerella_format <- czidHIVPosGardnerella %>%
  mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) )
czidHIVNegGardnerella_format <- czidHIVNegGardnerella %>%
  mutate(count = ifelse(count == 0, min(count[count > 0]) / 1000, count) )
czidHIVAdult_CDF <-
  ggplot() +
  geom_point(data = cbind(czidHIVPosGardnerella_format, percentile = scaledPercentile(czidHIVPosGardnerella$rowNum)),
             aes(x=percentile, y=count, color="HIVpositive")) +
  geom_point(data = cbind(czidHIVNegGardnerella_format, percentile = scaledPercentile(czidHIVNegGardnerella$rowNum)),
             aes(x=percentile, y=count, color="HIVnegative")) +
  labs(x = "samples (mothers, by percentile)", y = "counts found") +
  scale_color_manual(values = c("HIVpositive" = "red", "HIVnegative" = "blue")) +
  labs(title = "CDF of CZID counts of Gardnerella bacteria")
czidHIVAdult_CDF
# Gardnerella HIV Pos + Neg log overlay
czidHIVAdultLog_CDF <- czidHIVAdult_CDF + scale_y_log10() +
  labs(title = "CDF of CZID counts of Gardnerella bacteria (log scale)")
czidHIVAdultLog_CDF

##################################################################################
## write filtered data back to file
write.csv(pathseq_filteredByMean, pathseqFilteredPath, row.names = FALSE)
write.csv(czid_filteredByMean, czidFilteredPath, row.names = FALSE)

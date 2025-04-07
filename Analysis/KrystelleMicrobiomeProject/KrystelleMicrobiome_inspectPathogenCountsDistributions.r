library("tidyr")
library("dplyr")
library("stringr")
library("ggplot2")
source("/data/home/parkercol/PathSeq/Analysis/pathseqDataReformattingFunctions.r")

## Define variable paths

# import csv file with genus counts
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/pathseqTaxonomicSummaries/RNA_genus_tpm.csv"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/summarizedData/CZID_nr_genus_counts.csv"
# define paths for writing filtered files
pathseqFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/pathseq_genus_tpm_filtered.csv"
czidFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/summarizedData/czid_nr_genus_counts_filtered.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"

taxLevel = "genus"
###############################################################################
## Define data

# define pathseq all
pathseqCounts <- loadCountsData(pathseqPath)
czidCounts <- loadCountsData(czidPath)

# define sample sheet with treatment information
sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(SampleName, Treatment)) %>%
  filter(!str_detect(tolower(SampleName), tolower("positive")) &
           !str_detect(tolower(SampleName), tolower("water")))

# define datasets filtered by everything below 3std above negative control mean
pathseqMeanFilterInfo <- filterMean(pathseqCounts)
pathseqFilteredByStd <- pathseqMeanFilterInfo[3]
czidMeanFilterInfo <- filterMean(czidCounts)
czidFilteredByStd <- czidMeanFilterInfo[3]

# define datasets filtered of everything below 99th percentile of negative control
pathseqPercFilterInfo <- filterMedian(pathseqCounts)
pathseqFilteredByPerc <- pathseqPercFilterInfo[3]
czidPercFilterInfo <- filterMedian(czidCounts)
czidFilteredByPerc <- czidPercFilterInfo[3]

# define pathseq merged with treatment information
pathseqMergeTreat <- countsWithMetadata(pathseqCounts, sampleSheet)
czidMergeTreat <- countsWithMetadata(czidCounts, sampleSheet)

## define data groups by treatment labels
pathseqHIVPos <- pathseqMergeTreat %>% filter(Treatment == "Positive") %>%
  select(-Treatment)
pathseqHIVNeg <- pathseqMergeTreat %>% filter(Treatment == "Negative") %>%
  select(-Treatment)
pathseqChild <- pathseqMergeTreat %>% filter(Treatment == "HEU") %>%
  filter(Treatment == "HUU") %>%
  select(-Treatment)

czidHIVPos <- czidMergeTreat %>% filter(Treatment == "Positive") %>%
  select(-Treatment)
czidHIVNeg <- czidMergeTreat %>% filter(Treatment == "Negative") %>%
  select(-Treatment)

################################################################################
## look at counts distribution and compare CZID to pathseq, filtering to non-filtering

## get pooled pathseq and CZID plots before filtering
pathseqPooled <- pivot_longer(pathseqCounts,
                              cols = !genus,
                              names_to = "Sample",
                              values_to = "count")
czidPooled <- pivot_longer(czidCounts,
                           cols = !genus,
                           names_to = "Sample",
                           values_to = "count")

# plot pathseq HIV positive for filtering
pathseqHIVPosPooled <- pivot_longer(pathseqHIVPos,
                                    cols = starts_with("k_"),
                                    names_to = "Pathogen",
                                    values_to = "count")
pathseqHIVPosPooled$count <- as.numeric(pathseqHIVPosPooled$count)
pathseqHIVPosPooled <- pathseqHIVPosPooled[order(pathseqHIVPosPooled$count), ]
pathseqHIVPosPooled$rowNum <- as.numeric(rownames(pathseqHIVPosPooled))

pathseqPooled$count <- as.numeric(pathseqPooled$count)
pathseq_ordered <- pathseqPooled[order(pathseqPooled$count), ]
pathseq_ordered$rowNum <- as.numeric(rownames(pathseq_ordered)) %>% replace(is.na(.), 0)
czidPooled$count <- as.numeric(czidPooled$count)
czid_ordered <- czidPooled[order(czidPooled$count), ]
czid_ordered$rowNum <- as.numeric(rownames(czid_ordered)) %>% replace(is.na(.), 0)

cdf_zeros <-
  basePlot("CDF of pathogens found from Krystelle's microbiome data", "sample-pathogen combination", "counts found") +
  geom_point(data = pathseq_ordered, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czid_ordered, aes(x=rowNum, y=count), color='orange')
cdf_zeros

pathseq_ordered_noZeros <- subset(pathseq_ordered, count !=0)
czid_ordered_noZeros$rowNum <- as.numeric(rownames(czid_ordered_noZeros)) %>% replace(is.na(.), 0)
czid_ordered_noZeros <- subset(czid_ordered, count !=0)
czid_ordered_noZeros$rowNum <- as.numeric(rownames(czid_ordered_noZeros)) %>% replace(is.na(.), 0)

cdf <-
  basePlot("Pathogens found by Pathseq vs CZID", "sample-pathogen combination", "amount of pathogen found (tpm)") +
  geom_point(data = pathseq_ordered_noZeros, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czid_ordered_noZeros, aes(x=rowNum, y=count), color='orange')
cdf

################################################################################
## Examine the overlap between adults and children


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

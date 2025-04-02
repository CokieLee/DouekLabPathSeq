library("tidyr")
library("dplyr")
library("stringr")
library("ggplot2")

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
pathseqCounts <- read.csv(pathseqPath) %>%
  subset(genus != "taxonomy") %>%
  mutate(across(-1, as.numeric)) %>%
  replace(is.na(.), 0) %>%
  rename_with(~ str_split(., "_", simplify = TRUE)[, 2], -1)

# define czid all
czidCounts <- read.csv(czidPath) %>%
  mutate(across(-1, as.numeric)) %>%
  replace(is.na(.), 0) %>%
  rename_with(~str_split(., "_", simplify = TRUE)[, 2], -1) %>%
  rename(genus = Pathogens)

# define sample sheet with treatment information
sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(Sample.name, HIV.status)) %>%
  filter(!str_detect(tolower(Sample.name), tolower("positive")) &
           !str_detect(tolower(Sample.name), tolower("water")))

# define pathseq filtered by negative control
meanPathseqCtrl <- mean(as.numeric(pathseqCounts$Water))
stdPathseqCtrl <- sd(as.numeric(pathseqCounts$Water))
pathseqThreshByMean <- meanPathseqCtrl + (3 * stdPathseqCtrl)
medianPathseqCtrl <- median(as.numeric(pathseqCounts$Water))
pathseqThreshByPercent <- quantile(pathseqCounts$Water, probs=c(0.99))[1] %>%
  unname()

pathseq_filteredByMean <- pathseqCounts %>%
  mutate(across(-1, ~ ifelse(. < pathseqThreshByMean, 0, .) ))
pathseq_filteredByPerc <- pathseqCounts %>%
  mutate(across(-1, ~ ifelse(. < pathseqThreshByPercent, 0, .)))

# define czid filtered by negative control
meanCzidCtrl <- mean(czidCounts$Water)
stdCzidCtrl <- sd(czidCounts$Water)
czidThreshByMean <- meanCzidCtrl + (3 * stdCzidCtrl)
medianCzidCtrl <- median(as.numeric(czidCounts$Water))
czidThreshByPercent <- quantile(czidCounts$Water, probs=c(0.99))[1] %>%
  unname()

czid_filteredByMean <- czidCounts %>%
  mutate(across(-1, ~ ifelse(. < czidThreshByMean, 0, .)))
czid_filteredByPerc <- czidCounts %>%
  mutate(across(-1, ~ ifelse(. < czidThreshByPercent, 0, .)) )

# define pathseq merged with treatment information
# TODO: fix this into dplyr style
pathseqMergeTreat <- t(pathseq_filteredByPerc)
colnames(pathseqMergeTreat) <- pathseqMergeTreat[1, ]
pathseqMergeTreat <- tail(pathseqMergeTreat, -1)
pathseqMergeTreat <- cbind(Sample.name = rownames(pathseqMergeTreat), pathseqMergeTreat )
pathseqMergeTreat <- merge(pathseqMergeTreat, sampleSheet, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define czid merged with treatment information
czidMergeTreat <- t(czid_filteredByPerc)
colnames(czidMergeTreat) <- czidMergeTreat[1, ]
czidMergeTreat <- tail(czidMergeTreat, -1)
czidMergeTreat <- cbind(Sample.name = rownames(czidMergeTreat), czidMergeTreat)
czidMergeTreat <- merge(czidMergeTreat, sampleSheet, by = "Sample.name")

# define pathseq HIV+ set
pathseqHIVPos <- pathseqMergeTreat %>% filter(HIV.status == "Positive") %>%
  select(-HIV.status)

# define pathseq HIV- set
pathseqHIVNeg <- pathseqMergeTreat %>% filter(HIV.status == "Negative") %>%
  select(-HIV.status)

# define pathseq children set
pathseqChild <- pathseqMergeTreat %>% filter(HIV.status == "HEU") %>%
  filter(HIV.status == "HUU") %>%
  select(-HIV.status)

# define czid HIV+ set
czidHIVPos <- czidMergeTreat %>% filter(HIV.status == "Positive") %>%
  select(-HIV.status)

# define czid HIV- set
czidHIVNeg <- czidMergeTreat %>% filter(HIV.status == "Negative") %>%
  select(-HIV.status)

# define czid children set

###############################################################################
## Analysis functions

# make histograms showing counts for nr and nt to show equivalent comparisons
countsHistogram <- function(data) {
  plot <- ggplot(data, aes(x=count)) +
    geom_histogram(binwidth = 100000, fill = "blue") +
    labs(title = "Histogram of all counts",
         x = "counts of pathogen found",
         y = "number of sample-pathogen pairs with this count")
  
  return(plot)
}

# make CDF showing counts for nr and nt
countsCDF <- function( data ) {
  plot <- ggplot( data , aes(x=rowNum, y=count )) +
    geom_point() +
    labs(title = "CDF of pathogens found from Krystelle's microbiome data",
         x = "sample-pathogen combination",
         y = "counts found")
  
  return(plot)
}

scaledPercentile <- function(countsCol) {
  n <- length(countsCol)
  rank <- rank(countsCol)
  percentile <- (rank - 1) / (n - 1)
  scaled_percentile <- percentile * (n / (n + 1))
  
  return(scaled_percentile)
}

################################################################################
## format data before plotting

## plot pathseq and CZID cdf of counts before filtering
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
################################################################################
## define graphs
pathseqPooled$count <- as.numeric(pathseqPooled$count)
pathseq_ordered <- pathseqPooled[order(pathseqPooled$count), ]
pathseq_ordered$rowNum <- as.numeric(rownames(pathseq_ordered)) %>% replace(is.na(.), 0)
czidPooled$count <- as.numeric(czidPooled$count)
czid_ordered <- czidPooled[order(czidPooled$count), ]
czid_ordered$rowNum <- as.numeric(rownames(czid_ordered)) %>% replace(is.na(.), 0)

pathseq_cdf_zeros <-
  ggplot() +
  geom_point(data = pathseq_ordered, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czid_ordered, aes(x=rowNum, y=count), color='orange') +
  labs(title = "CDF of pathogens found from Krystelle's microbiome data",
       x = "sample-pathogen combination",
       y = "counts found") +
  theme_bw()
pathseq_cdf_zeros

pathseq_ordered_noZeros <- subset(pathseq_ordered, count !=0)
czid_ordered_noZeros$rowNum <- as.numeric(rownames(czid_ordered_noZeros)) %>% replace(is.na(.), 0)
czid_ordered_noZeros <- subset(czid_ordered, count !=0)
czid_ordered_noZeros$rowNum <- as.numeric(rownames(czid_ordered_noZeros)) %>% replace(is.na(.), 0)

pathseq_cdf <-
  ggplot() +
  geom_point(data = pathseq_ordered_noZeros, aes(x=rowNum, y=count ), color='green') +
  geom_point(data = czid_ordered_noZeros, aes(x=rowNum, y=count), color='orange') +
  labs(title = "Pathogens found by Pathseq vs CZID",
       x = "sample-pathogen combination",
       y = "amount of pathogen found (tpm)") +
  theme_bw()
pathseq_cdf



################################################################################
# plot pathseq graphs
pathseqPooled$count <- as.numeric(pathseqPooled$count)
pathseq_ordered <- pathseqPooled[order(pathseqPooled$count), ]
pathseq_ordered$rowNum <- as.numeric(rownames(pathseq_ordered)) %>% replace(is.na(.), 0)
pathseq_cdf_zeros <- countsCDF(pathseq_ordered)

pathseq_ordered_noZeros <- subset(pathseq_ordered, count !=0)
pathseq_cdf <- countsCDF(pathseq_ordered_noZeros)

# plot czid graphs
czidPooled$count <- as.numeric(czidPooled$count)
czid_ordered <- czidPooled[order(czidPooled$count), ]
czid_ordered$rowNum <- as.numeric(rownames(czid_ordered)) %>% replace(is.na(.), 0)
czid_cdf_zeros <- countsCDF(czid_ordered)

czid_ordered_noZeros <- subset(czid_ordered, count !=0)
czid_cdf <- countsCDF(czid_ordered_noZeros)

filteredType1 <- pathseq_filteredByMean
filteredType2 <- czid_filteredByMean

## plot filtered pathseq and CZID cdf of counts
pathseqPooledFiltered <- pivot_longer(filteredType1,
                              cols = !genus,
                              names_to = "Sample",
                              values_to = "count")
czidPooledFiltered <- pivot_longer(filteredType2,
                           cols = !genus,
                           names_to = "Sample",
                           values_to = "count")

pathseqPooledFiltered_ordered <- pathseqPooledFiltered[order(pathseqPooledFiltered$count), ]
pathseqPooledFiltered_ordered$rowNum <- as.numeric(rownames(pathseqPooledFiltered_ordered))
pathseqPooledFiltered_ordered <- subset(pathseqPooledFiltered_ordered, count != 0)
pathseqPooledFiltered_cdf <- countsCDF(pathseqPooledFiltered_ordered)
pathseqPooledFiltered_cdf

czidPooledFiltered_ordered <- czidPooledFiltered[order(czidPooledFiltered$count), ]
czidPooledFiltered_ordered$rowNum <- as.numeric(rownames(czidPooledFiltered_ordered))
czidPooledFiltered_ordered <- subset(czidPooledFiltered_ordered, count != 0)
czidPooledFiltered_cdf <- countsCDF(czidPooledFiltered_ordered)
czidPooledFiltered_cdf

# plot HIV positive CDF
pathseqHIVPosCDF <- countsCDF(pathseqHIVPosPooled)

# pathseq plot HIV positive Gardnerella CDF
pathseqHIVPosGardnerella <- pathseqHIVPos %>%
  select(c(Sample.name, `k__Bacteria;-p__Actinobacteria;-c__Actinomycetia;-o__Bifidobacteriales;-f__Bifidobacteriaceae;-g__Gardnerella`))
colnames(pathseqHIVPosGardnerella) <- c("Sample.name", "count")
pathseqHIVPosGardnerella$count <- as.numeric(pathseqHIVPosGardnerella$count)
pathseqHIVPosGardnerella <- pathseqHIVPosGardnerella[order(pathseqHIVPosGardnerella$count), ]
pathseqHIVPosGardnerella$rowNum <- sort( as.numeric( rownames(pathseqHIVPosGardnerella) ) )

# pathseq plot HIV negative Gardnerella CDF
pathseqHIVNegGardnerella <- pathseqHIVNeg %>%
  select(c(Sample.name, `k__Bacteria;-p__Actinobacteria;-c__Actinomycetia;-o__Bifidobacteriales;-f__Bifidobacteriaceae;-g__Gardnerella`))
colnames(pathseqHIVNegGardnerella) <- c("Sample.name", "count")
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
  select(c(Sample.name, `k__bacteria;-g__Gardnerella`))
colnames(czidHIVPosGardnerella) <- c("Sample.name", "count")
czidHIVPosGardnerella$count <- as.numeric(czidHIVPosGardnerella$count)
czidHIVPosGardnerella <- czidHIVPosGardnerella[order(czidHIVPosGardnerella$count), ]
czidHIVPosGardnerella$rowNum <- sort( as.numeric( rownames(czidHIVPosGardnerella) ) )

# czid plot HIV negative Gardnerella CDF
czidHIVNegGardnerella <- czidHIVNeg %>%
  select(c(Sample.name, `k__bacteria;-g__Gardnerella`))
colnames(czidHIVNegGardnerella) <- c("Sample.name", "count")
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

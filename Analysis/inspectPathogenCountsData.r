library("tidyr")

# import csv file with genus counts
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/diversityMetrics/RNA_genus_tpm.csv"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/diversityMetrics/CZID_nr_species_counts.csv"

# define paths for writing filtered files
pathseqFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/diversityMetrics/pathseq_genus_tpm_filtered.csv"
czidFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/diversityMetrics/czid_nr_genus_counts_filtered.csv"

PathogenCounts1 <- read.csv(pathseqPath) %>%
  subset(genus != "taxonomy")
PathogenCounts1 <- mutate(PathogenCounts1, across(-1, as.numeric))
PathogenCounts1 <- PathogenCounts1 %>% replace(is.na(.), 0)
PathogenCounts2 <- read.csv(czidPath)
PathogenCounts2 <- PathogenCounts2 %>% replace(is.na(.), 0)

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

## plot pathseq and CZID cdf of counts before filtering
pathseqPooled <- pivot_longer(PathogenCounts1,
                              cols = !genus,
                              names_to = "Sample",
                              values_to = "count")
czidPooled <- pivot_longer(PathogenCounts2,
                           cols = !genus,
                           names_to = "Sample",
                           values_to = "count")

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

## filter counts that are more than 2 standard deviations above negative control mean
meanPathseqCtrl <- mean(as.numeric(PathogenCounts1$Sample_Water_Ctrl_S61_L001_R1_001))
stdPathseqCtrl <- sd(as.numeric(PathogenCounts1$Sample_Water_Ctrl_S61_L001_R1_001))
pathseqThreshByMean <- meanPathseqCtrl + (3 * stdPathseqCtrl)

meanCzidCtrl <- mean(PathogenCounts2$Sample_Water_Ctrl_S61_681293)
stdCzidCtrl <- sd(PathogenCounts2$Sample_Water_Ctrl_S61_681293)
czidThreshByMean <- meanCzidCtrl + (3 * stdCzidCtrl)

medianPathseqCtrl <- median(as.numeric(PathogenCounts1$Sample_Water_Ctrl_S61_L001_R1_001))
pathseqThreshByPercent <- quantile(PathogenCounts1$Sample_Water_Ctrl_S61_L001_R1_001, probs=c(0.99))[1] %>%
  unname()

medianCzidCtrl <- median(as.numeric(PathogenCounts2$Sample_Water_Ctrl_S61_681293))
czidThreshByPercent <- quantile(PathogenCounts2$Sample_Water_Ctrl_S61_681293, probs=c(0.99))[1] %>%
  unname()

Path1_filteredByMean <- PathogenCounts1
Path1_filteredByMean[Path1_filteredByMean < pathseqThreshByMean] <- 0

Path1_filteredByPerc <- PathogenCounts1
Path1_filteredByPerc[Path1_filteredByPerc < pathseqThreshByPercent] <- 0

Path2_filteredByMean <- PathogenCounts2
Path2_filteredByMean[Path2_filteredByMean < czidThreshByMean] <- 0
Path2_filteredByPerc <- PathogenCounts2
Path2_filteredByPerc[Path2_filteredByPerc < pathseqThreshByPercent]

filteredType1 <- Path1_filteredByMean
filteredType2 <- Path2_filteredByMean

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

czidPooledFiltered_ordered <- czidPooledFiltered[order(czidPooledFiltered$count), ]
czidPooledFiltered_ordered$rowNum <- as.numeric(rownames(czidPooledFiltered_ordered))
czidPooledFiltered_ordered <- subset(czidPooledFiltered_ordered, count != 0)
czidPooledFiltered_cdf <- countsCDF(czidPooledFiltered_ordered)

## write filtered data back to file
write.csv(Path1_filteredByMean, pathseqFilteredPath, row.names = FALSE)
write.csv(Path2_filteredByMean, czidFilteredPath, row.names = FALSE)


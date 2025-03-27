library("dplyr")
library("ggplot2")

## Define variable paths
# import csv file with genus counts
pathseqFilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/pathseq_genus_tpm_filtered.csv"
pathseqUnfilteredPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/pathseqTaxonomicSummaries/RNA_genus_tpm.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"

# output paths
topGenusListPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10GenusList.csv"
topGenusInfoPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_filtered.csv"
topGenusInfoPathUnfiltered <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
topN = 10

###############################################################################
## Define data

# define sample sheet with treatment information
sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(Sample.name, HIV.status)) %>%
  filter(!str_detect(tolower(Sample.name), tolower("positive")) &
           !str_detect(tolower(Sample.name), tolower("water")))

# define pathseq all
pathseqCounts <- read.csv(pathseqFilteredPath) %>%
  subset(genus != "taxonomy") %>%
  mutate(across(-1, as.numeric)) %>%
  replace(is.na(.), 0) %>%
  mutate(genus = (sapply(str_split(genus, "_"), function(x) paste0("g_", tail(x, 1)) ) ) ) %>%
  filter(genus != "g_")

pathseqTransposed <- t(pathseqCounts)
colnames(pathseqTransposed) <- pathseqTransposed[1, ]
pathseqTransposed <- pathseqTransposed[-1, ]

pathseqMergeTreat <- pathseqTransposed %>%
  bind_cols(Sample.name = rownames(pathseqTransposed)) %>%
  merge(sampleSheet, by = "Sample.name") %>%
  mutate(across(-all_of(c("HIV.status", "Sample.name")), ~as.numeric(.)))

# define pathseq all for unfiltered data
pathseqUnfilteredCounts <- read.csv(pathseqUnfilteredPath) %>%
  subset(genus != "taxonomy") %>%
  mutate(across(-1, as.numeric)) %>%
  replace(is.na(.), 0) %>%
  rename_with(~str_split(., "_", simplify = TRUE)[, 2], -1) %>%
  mutate(genus = (sapply(str_split(genus, "_"), function(x) paste0("g_", tail(x, 1)) ) ) ) %>%
  filter(genus != "g_")

pathseqUnfilteredTransposed <- t(pathseqUnfilteredCounts)
colnames(pathseqUnfilteredTransposed) <- pathseqUnfilteredTransposed[1, ]
pathseqUnfilteredTransposed <- pathseqUnfilteredTransposed[-1, ]

pathseqUnfilteredMergeTreat <- pathseqUnfilteredTransposed %>%
  bind_cols(Sample.name = rownames(pathseqUnfilteredTransposed)) %>%
  merge(sampleSheet, by = "Sample.name") %>%
  mutate(across(-all_of(c("HIV.status", "Sample.name")), ~as.numeric(.)))

###############################################################################
## get highest counts genuses across samples

pathseqSummedCounts <- pathseqTransposed %>%
  as.data.frame() %>%
  mutate_all(as.numeric) %>%
  summarize(across(everything(), sum) ) %>%
  t() %>%
  data.frame() %>%
  rename_with(~ c("counts")) %>%
  arrange(desc(counts))

pathseqHighestGenuses <- head(pathseqSummedCounts, topN)
write.csv(pathseqHighestGenuses, topGenusListPath, row.names = TRUE)

topCountsAndInfo <- pathseqMergeTreat %>%
  select( "HIV.status", "Sample.name", any_of(rownames(pathseqHighestGenuses)) )
write.csv(topCountsAndInfo, topGenusInfoPath, row.names = FALSE)

topUnfilteredCountsAndInfo <- pathseqUnfilteredMergeTreat %>%
  select( "HIV.status", "Sample.name", any_of(rownames(pathseqHighestGenuses)) )
write.csv(topUnfilteredCountsAndInfo, topGenusInfoPathUnfiltered, row.names = FALSE)
###############################################################################
## define heatmap

pathseqLong <- pathseqMergeTreat %>%
  pivot_longer(cols = starts_with("g_"),
               names_to = "Pathogen",
               values_to = "counts")

pathseqLongTop <- pathseqLong %>%
  filter(Pathogen %in% rownames(pathseqHighestGenuses) )

pathseqHeatMap <- ggplot(pathseqLongTop, aes(x = Sample.name, y = Pathogen, fill = counts)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  facet_wrap(~ HIV.status, scales = "free_x")
pathseqHeatMap
  
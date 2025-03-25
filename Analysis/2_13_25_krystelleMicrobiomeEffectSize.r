library(stringr)
library(dplyr)
library(purrr)
library(vegan)

## read in sample sheet and pathseq and czid counts
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/diversityMetrics/pathseq_genus_tpm_filtered.csv"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/diversityMetrics/czid_nr_genus_counts_filtered.csv"
sampleSheetPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Moritz_Sample_sheet_infection_status.txt"

##############################################################################
## define data to run on
sampleSheet <-
  read.csv(sampleSheetPath, sep="\t") %>%
  subset(select=c(Sample.name, HIV.status))
sampleSheet <- sampleSheet[order(sampleSheet$Sample.name), ]
sampleSheet <- head(sampleSheet, -2)

## sorts sampleSheet values into infection vs non-infection
infectionSplitSheet <- sampleSheet
infectionSplitSheet$HIV.status <-
  sapply(infectionSplitSheet$HIV.status,
         function (x) if(x != "Positive") {"Uninfected"} else {"Positive"}
        )

adultPosNegSplit <- sampleSheet
adultPosNegSplit <-
  adultPosNegSplit %>%
  filter(HIV.status != "HUU") %>%
  filter(HIV.status != "HEU")
  
childExposureSplit <- sampleSheet
childExposureSplit <-
  childExposureSplit %>%
  filter(HIV.status != "Positive") %>%
  filter(HIV.status != "Negative")

# define datasets for each split and merge with correct sample sheet
pathseqCounts <-
  read.csv(pathseqPath) %>%
  t() %>%
  replace(is.na(.), 0) %>%
  as.data.frame()
colnames(pathseqCounts) <- pathseqCounts[1, ]
pathseqCounts <- pathseqCounts[-1,]
pathseqCounts$Sample.name <- unlist(
  sapply(rownames(pathseqCounts), function(x) str_split_1(x, "_")[2])
)
pathseqCounts <- pathseqCounts[order(pathseqCounts$Sample.name), ]
pathseqCounts <- head(pathseqCounts, -2)
pathseqBacteriaCols <- grepl("Bacteria", names(pathseqCounts)) | (names(pathseqCounts) == "Sample.name")
pathseqCounts <- pathseqCounts[, pathseqBacteriaCols]

czidCounts <-
  read.csv(czidPath) %>%
  t() %>%
  replace(is.na(.), 0) %>%
  as.data.frame()
colnames(czidCounts) <- czidCounts[1, ]
czidCounts <- czidCounts[-1, ]
czidCounts$Sample.name <- unlist(
  sapply(rownames(czidCounts), function(x) str_split_1(x, "_")[2])
)
czidCounts <- czidCounts[order(czidCounts$Sample.name), ]
czidCounts <- head(czidCounts, -2)
czidBacteriaCols <- grepl("bacteria", names(czidCounts)) | (names(czidCounts) == "Sample.name")
czidCounts <- czidCounts[, czidBacteriaCols]

# define pathseq positive vs uninfected (all samples) split
pathseqInfectionMerged <- merge(pathseqCounts, infectionSplitSheet, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define czid positive vs uninfected (all samples) split
czidInfectionMerged <- merge(czidCounts, infectionSplitSheet, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define czid counts with adult infection split
czidAdultMerged <- merge(czidCounts, adultPosNegSplit, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define pathseq counts with adult infection split
pathseqAdultMerged <- merge(pathseqCounts, adultPosNegSplit, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define czid counts with child exposure split
czidChildMerged <- merge(czidCounts, childExposureSplit, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# define pathseq counts with child exposure split
pathseqChildMerged <- merge(pathseqCounts, childExposureSplit, by = "Sample.name", all.x = FALSE, all.y = FALSE)

# get pathogen names lists of interest
pathseqPathNames <- names(pathseqCounts)[names(pathseqCounts) != "Sample.name"]
czidPathNames <- names(czidCounts)[names(czidCounts) != "Sample.name"]

##############################################################################
## define functions

## perform t-tests for each pathogen across all samples
# perform_t_test <- function(pathogen, treatment) {
#   return(t.test(pathogen ~ treatment))
# }
# extract_p_value <- function(t_test_result) {
#   return(t_test_result$p.value)
# }
# 
# t_test_results <-
#   map(
#     pathogenNames,
#     function(pathogen_name) {
#       path <- as.numeric(infect1_merged[[pathogen_name]])
#       perform_t_test(path, infect1_merged$HIV.status)
#     }
# )
# p_values_unadjusted <- map_dbl(t_test_results, extract_p_value)
# sigIndices_unadjusted <- p_values_unadjusted < 0.05
# sigPaths_unadjusted <- data.frame(species = pathogenNames[sigIndices_unadjusted],
#                                   p_values_unadjusted = p_values_unadjusted[sigIndices_unadjusted])
# 
# p_values <- p.adjust(p_values_unadjusted, method = "bonferroni")
# significantIndices <- p_values < 0.05
# significantPaths <- data.frame(species = pathogenNames[significantIndices],
#                                   p_vals = p_values[significantIndices])


# wilcox test function
runWilcoxTest <- function(data, pathogenNames) {
  # remove all data columns that are zero only (wilcox cant handle all 0 variance)
  char_cols <- c("Sample.name", "HIV.status")
  counts_wilcox <- data %>%
    mutate(across(-all_of(char_cols), ~ as.numeric(.))) %>%
    select(where(~ (is.numeric(.)) && (sum(.) != 0) ), all_of(char_cols) )
  
  # downselect pathogen names
  wilcox_pathNames <- intersect(pathogenNames, colnames(counts_wilcox))
  
  # wilcox test
  wilcox_results <-
    map(
      wilcox_pathNames,
      function(pathogen_name) {
        path <- counts_wilcox[, c(pathogen_name, "HIV.status")]
        wilcox.test(as.numeric(path[[pathogen_name]]) ~ path[["HIV.status"]] )
      }
    )
  
  wilcox_pVals <- map_dbl(wilcox_results, function(x) x$p.value)
  summaryDf <-
    data.frame(
    pathogenNames = wilcox_pathNames,
    pVals = wilcox_pVals)

  return(summaryDf[order(summaryDf$pVals), ])
}

getMeanDiff <- function(data, pathogenNames, treatA, treatB) {
  pathNames <- intersect(pathogenNames, colnames(data))
  
  meanDiff <-
    map(
      pathNames,
      function(pathName) {
        path <- data[, c(pathName, "HIV.status")]
        return(
          mean( as.numeric(filter(path, HIV.status == treatA)[[pathName]]) )
          - mean( as.numeric(filter(path, HIV.status == treatB)[[pathName]]) )
        )
      }
    )

  return(data.frame(pathogenNames = pathNames, meanDiff = unlist(meanDiff)))
}

getBrayCurtisDiff <- function(data, pathogenNames, treatA, treatB) {
  
}

#############################################################################
## Compute results

# Positive HIV vs all uninfected
pathseq_wilcox_pVals_a <- runWilcoxTest(pathseqInfectionMerged, pathseqPathNames)
pathseq_wilcox_effectSize_a <- getMeanDiff(pathseqInfectionMerged, pathseqPathNames, "Positive", "Uninfected")
pathseq_wilcox_bonferroni_a <- pathseq_wilcox_pVals_a %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
pathseq_wilcox_hochberg_a <- pathseq_wilcox_pVals_a %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
pathseq_wilcox_summarized_a <- data.frame(
  pathogenNames = pathseq_wilcox_pVals_a$pathogenNames,
  unadjusted = pathseq_wilcox_pVals_a$pVals,
  bonferroni = pathseq_wilcox_bonferroni_a$pVals,
  hochberg = pathseq_wilcox_hochberg_a$pVals) %>%
  merge(pathseq_wilcox_effectSize_a, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(pathseq_wilcox_summarized_a)

czid_wilcox_pVals_a <- runWilcoxTest(czidInfectionMerged, czidPathNames)
czid_wilcox_effectSize_a <- getMeanDiff(czidInfectionMerged, czidPathNames, "Positive", "Uninfected")
czid_wilcox_bonferroni_a <- czid_wilcox_pVals_a %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
czid_wilcox_hochberg_a <- czid_wilcox_pVals_a %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
czid_wilcox_summarized_a <- data.frame(
  pathogenNames = czid_wilcox_pVals_a$pathogenNames,
  unadjusted = czid_wilcox_pVals_a$pVals,
  bonferroni = czid_wilcox_bonferroni_a$pVals,
  hochberg = czid_wilcox_hochberg_a$pVals) %>%
  merge(czid_wilcox_effectSize_a, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(czid_wilcox_summarized_a)

# Adult disease status (pos/neg) wilcox test
pathseq_wilcox_pVals_b <- runWilcoxTest(pathseqAdultMerged, pathseqPathNames)
pathseq_wilcox_effectSize_b <- getMeanDiff(pathseqAdultMerged, pathseqPathNames, "Positive", "Negative")
pathseq_wilcox_bonferroni_b <- pathseq_wilcox_pVals_b %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
pathseq_wilcox_hochberg_b <- pathseq_wilcox_pVals_b %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
pathseq_wilcox_summarized_b <- data.frame(
  pathogenNames = pathseq_wilcox_pVals_b$pathogenNames,
  unadjusted = pathseq_wilcox_pVals_b$pVals,
  bonferroni = pathseq_wilcox_bonferroni_b$pVals,
  hochberg = pathseq_wilcox_hochberg_b$pVals) %>%
  merge(pathseq_wilcox_effectSize_b, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(pathseq_wilcox_summarized_b)

czid_wilcox_pVals_b <- runWilcoxTest(czidAdultMerged, czidPathNames)
czid_wilcox_effectSize_b <- getMeanDiff(czidAdultMerged, czidPathNames, "Positive", "Negative")
czid_wilcox_bonferroni_b <- czid_wilcox_pVals_b %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
czid_wilcox_hochberg_b <- czid_wilcox_pVals_b %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
czid_wilcox_summarized_b <- data.frame(
  pathogenNames = czid_wilcox_pVals_b$pathogenNames,
  unadjusted = czid_wilcox_pVals_b$pVals,
  bonferroni = czid_wilcox_bonferroni_b$pVals,
  hochberg = czid_wilcox_hochberg_b$pVals) %>%
  merge(czid_wilcox_effectSize_b, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(czid_wilcox_summarized_b)

# Child disease exposure status (HEU vs HUU) wilcox test
pathseq_wilcox_pVals_c <- runWilcoxTest(pathseqChildMerged, pathseqPathNames)
pathseq_wilcox_effectSize_c <- getMeanDiff(pathseqChildMerged, pathseqPathNames, "HEU", "HUU")
pathseq_wilcox_bonferroni_c <- pathseq_wilcox_pVals_c %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
pathseq_wilcox_hochberg_c <- pathseq_wilcox_pVals_c %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
pathseq_wilcox_summarized_c <- data.frame(
  pathogenNames = pathseq_wilcox_pVals_c$pathogenNames,
  unadjusted = pathseq_wilcox_pVals_c$pVals,
  bonferroni = pathseq_wilcox_bonferroni_c$pVals,
  hochberg = pathseq_wilcox_hochberg_c$pVals) %>%
  merge(pathseq_wilcox_effectSize_c, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(pathseq_wilcox_summarized_c)

czid_wilcox_pVals_c <- runWilcoxTest(czidChildMerged, czidPathNames)
czid_wilcox_effectSize_c <- getMeanDiff(czidChildMerged, czidPathNames, "HEU", "HUU")
czid_wilcox_bonferroni_c <- czid_wilcox_pVals_c %>% mutate(pVals = p.adjust(pVals, method = "bonferroni"))
czid_wilcox_hochberg_c <- czid_wilcox_pVals_c %>% mutate(pVals = p.adjust(pVals, method = "hochberg"))
czid_wilcox_summarized_c <- data.frame(
  pathogenNames = czid_wilcox_pVals_c$pathogenNames,
  unadjusted = czid_wilcox_pVals_c$pVals,
  bonferroni = czid_wilcox_bonferroni_c$pVals,
  hochberg = czid_wilcox_hochberg_c$pVals) %>%
  merge(czid_wilcox_effectSize_c, by = "pathogenNames", x.all = FALSE, y.all = FALSE) %>%
  arrange(unadjusted)
View(czid_wilcox_summarized_c)


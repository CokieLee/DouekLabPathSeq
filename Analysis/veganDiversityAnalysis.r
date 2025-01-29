library("vegan")
library("dplyr")
library("tidyr")
library("ggplot2")
library("purrr")

# import csv file with species counts
filePath1 <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/diversityMetrics/RNA_species_tpm.csv"
filePath2 <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/diversityMetrics/CZID_nr_species_counts.csv"
PathogenCounts1 <- read.csv(filePath1)
PathogenCounts1 <- PathogenCounts1 %>% replace(is.na(.), 0)
PathogenCounts2 <- read.csv(filePath2)
PathogenCounts2 <- PathogenCounts2 %>% replace(is.na(.), 0)

print(colnames(PathogenCounts1))
print(colnames(PathogenCounts2))

# file naming information
dateToday = "1_2_25"
dataInfo_1 = "Pathseq"
dataInfo_2 = "CZID"
projectLabel = "Krystelle_Microbiome"

## format data for vegan on separate input dataset
Path1_formatted <- data.frame(t( as.data.frame(PathogenCounts1)))
colnames(Path1_formatted) <- as.character(unlist(Path1_formatted[1, ]))
Path1_formatted <- Path1_formatted[-1,]
# alter the names of samples so they are unique to treatment method
rownames(Path1_formatted) <- paste(dataInfo_1, rownames(Path1_formatted), sep="_")

Path2_formatted <- data.frame(t( as.data.frame(PathogenCounts2)))
colnames(Path2_formatted) <- as.character(unlist(Path2_formatted[1, ]))
Path2_formatted <- Path2_formatted[-1,]
# alter the names of samples so they are unique to treatment method
rownames(Path2_formatted) <- paste(dataInfo_2, rownames(Path2_formatted), sep="_")

## combine the two different pathogen counts into one dataframe
Path1_formatted$row_name <- rownames(Path1_formatted)
Path1_formatted <- mutate_all(Path1_formatted,
                              function(x) as.numeric(as.character(x))) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))
Path2_formatted$row_name <- rownames(Path2_formatted)
Path2_formatted <- mutate_all(Path2_formatted,
                              function(x) as.numeric(as.character(x))) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))

## get species count (richness)
richness1 <- specnumber(Path1_formatted)
richness1 <- data.frame(Lab=dataInfo_1, Richness=unlist(richness1))
richness2 <- specnumber(Path2_formatted)
richness2 <- data.frame(Lab=dataInfo_2, Richness=unlist(richness2))
richness <- rbind(richness1, richness2)

Title = paste0(dateToday, "_PathogenRichness_", projectLabel)

richBoxes <-
    ggplot(
        data = richness,
        aes(x=Lab, y=Richness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=Title) +
        theme_bw()
ggsave(paste0(Title, ".png"))

## get vegan diversity scores
divScores1 <- diversity(Path1_formatted, index = "shannon")
divScores1 <- data.frame(Lab=dataInfo_1, Diversity=unlist(divScores1))
divScores2 <- diversity(Path2_formatted, index = "shannon")
divScores2 <- data.frame(Lab=dataInfo_2, Diversity=unlist(divScores2))
diversity <- rbind(divScores1, divScores2)

Title = paste0(dateToday, "_PathogenShannonDiversity_", projectLabel)

divBoxes <-
    ggplot(
          data = diversity,
          aes(x=Lab, y=Diversity, color=Lab),
          position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=Title) +
        theme_bw()
ggsave(paste0(Title, ".png"))

## get evenness score
evenScores1 <- divScores1$Diversity / log( richness1$Richness)
evenScores1 <- data.frame(Lab=dataInfo_1, Evenness=unlist(evenScores1))
evenScores2 <- divScores2$Diversity / log(richness2$Richness)
evenScores2 <- data.frame(Lab=dataInfo_2, Evenness=unlist(evenScores2))
Evenness <- rbind(evenScores1, evenScores2)

Title = paste(dateToday, "_PathogenEvenness_", projectLabel)

evenBoxes <-
    ggplot(
        data = Evenness,
        aes(x=Lab, y=Evenness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=Title) +
        theme_bw()
ggsave(paste0(Title, ".png"))

## append summary stats to main file
# PathogenCountsSummarized <- cbind(Sample=rownames(Path1_formatted),
#                                   "Treat1 Richness"=richness1$Richness, "Treat2 Richness"=richness2$Richness)
# write.csv(PathogenCountsSummarized, "PathogenCountsSummarized.csv")

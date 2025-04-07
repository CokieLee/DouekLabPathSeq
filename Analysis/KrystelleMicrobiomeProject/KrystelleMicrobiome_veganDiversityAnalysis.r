library("vegan")
library("dplyr")
library("tidyr")
library("ggplot2")
library("purrr")

# import csv file with species counts
pathseqPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/diversityMetrics/"
# filePath <- "/data/vrc_his/douek_lab/cokie/TestFiles/"
czidPath <- "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/CZID_Results/diversityMetrics/"

PathogenCounts1 <- read.csv(paste0(pathseqPath, "RNA_species_tpm.csv") )
PathogenCounts1 <- PathogenCounts1 %>% replace(is.na(.), 0)
PathogenCounts2 <- read.csv(paste0(czidPath, "CZID_nr_species_counts.csv") )
PathogenCounts2 <- PathogenCounts2 %>% replace(is.na(.), 0)

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
Path1_formatted$label <- dataInfo_1
Path2_formatted$row_name <- rownames(Path2_formatted)
Path2_formatted$label <- dataInfo_2
allData_formatted <- full_join(Path1_formatted, Path2_formatted)
rownames(allData_formatted) <- allData_formatted$row_name

# change all values to numeric values & remove rownames column
## TODO: figure out why NAs are being introduced by coersion
allData_formatted <- mutate_all(allData_formatted[, -1],
                                function(x) (as.numeric(as.character(x))) ) %>%
  mutate(across(everything(), ~tidyr::replace_na(., 0)))
# create a list corresponding to the order of samples in allData_formatted
# which says whether each sample is from analysis group 1 or 2
originLabList <- allData_formatted$label
# remove the rownames and labels columns
allData_formatted <- subset(allData_formatted, select=-c(row_name, label))

## get species count (richness)
richness <- specnumber(allData_formatted)
richFrame <- tibble(
    Path1 = richness[originLabList == dataInfo_1],
    Path2 = richness[originLabList == dataInfo_2]
)
rm()
ggplot(richFrame, aes(x=dataInfo_1, y=dataInfo_2)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, projectLabel, "_pathogen_richness_")) +
    theme_bw()
ggsave(paste(dateToday, "_RichnessScatter_", projectLabel, ".png"))

richBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Richness=unlist(richness)),
        aes(x=Lab, y=Richness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_PathogenCounts pathogen richness_", projectLabel)) +
        theme_bw()
ggsave(paste(dateToday, "_RichnessBoxPlots_", projectLabel, ".png"))

## get vegan diversity scores
divScores <- diversity(allData_formatted, index = "shannon")
divFrame <- tibble(
    Path1 = divScores[originLabList == dataInfo_1],
    Path2 = divScores[originLabList == dataInfo_2]
)
divScatter <- ggplot(divFrame, aes(x=dataInfo_1, y=dataInfo_2)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, "_PathogenCounts pathogen diversity_", projectLabel)) +
    theme_bw()
ggsave(paste(dateToday, "_DiversityScatter_", projectLabel, ".png"))

divBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Diversity=unlist(divScores)),
        aes(x=Lab, y=Diversity, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_PathogenCounts pathogen diversity_", projectLabel)) +
        theme_bw()
ggsave(paste(dateToday, "_DiversityBoxPlots_", projectLabel, ".png"))

## get evenness score
evenScores <- divScores / log(richness)
evenFrame <- tibble(
    Path1 = evenScores[originLabList == dataInfo_1],
    Path2 = evenScores[originLabList == dataInfo_2]
)
evenScatter <- ggplot(evenFrame, aes(x=dataInfo_1, y=dataInfo_2)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, "_PathogenCounts pathogen evenness_", projectLabel)) +
    theme_bw()
ggsave(paste(dateToday, "_EvennessScatter_", projectLabel, ".png"))

evenBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Evenness=unlist(evenScores)),
        aes(x=Lab, y=Evenness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_PathogenCounts pathogen evenness_", projectLabel)) +
        theme_bw()
ggsave(paste(dateToday, "_EvennessBoxPlots_", projectLabel, ".png"))

## append summary stats to main file
PathogenCountsSummarized <- cbind(Lab=unlist(originLabList), richness=richness, diversity=divScores, evenness=evenScores, PathogenCounts_formatted)
write.csv(PathogenCountsSummarized, paste0(projectLabel, "PathogenCountsSummarized.csv"))

# return to old working directory
setwd(old_wd)



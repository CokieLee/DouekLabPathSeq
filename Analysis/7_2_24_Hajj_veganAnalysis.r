library("vegan")
library("dplyr")
library("tidyr")
library("ggplot2")
library("purrr")

# import csv file with species counts
old_wd <- getwd()
filePath <- "C:\\Users\\cokie\\OneDrive - Parker Tong LLP\\PersonalFiles\\Work\\DouekLab\\PathSeq\\HajjComparisons\\"
setwd(filePath)
hajj <- read.csv("RNA_species_tpm.csv")
hajj_sheba <- read.csv("Filtered_Combined_output_Hajj_Sheba-Sequencing.tsv", sep='\t')
hajj_douek <- read.csv("Filtered_Combined_output_Seq_Core.tsv", sep='\t')
hajj_sheba <- hajj_sheba %>% replace(is.na(.), 0)
hajj_douek <- hajj_douek %>% replace(is.na(.), 0)

# file naming information
dateToday = "7_17_24"
dataInfo = "CZID"

## format data for vegan on combined dataset
# hajj <- as.data.frame(hajj)
# hajj_formatted <- data.frame(t(hajj))
# colnames(hajj_formatted) <- as.character(unlist(hajj_formatted[1, ]))
# hajj_ <- head(hajj_formatted[-1, ], -2)
# hajj_formatted <- mutate_all(hajj_formatted, function(x) (as.numeric(as.character(x))) )
# ## create an element for each sample saying which lab it originated from
# originLabList <- map(
#     rownames(hajj_formatted),
#     function(x) if(length(strsplit(x, "_")[[1]]) == 3) {"Sheba"} else {"Douek"}
# )

## format data for vegan on separate input dataset
sheba_formatted <- data.frame(t( as.data.frame(hajj_sheba)[, -1] ))
colnames(sheba_formatted) <- as.character(unlist(sheba_formatted[1, ]))
sheba_formatted <- sheba_formatted[-1, ]
# douk lab samples
douek_formatted <- data.frame(t( as.data.frame(hajj_douek)[, -1] ))
colnames(douek_formatted) <- as.character(unlist(douek_formatted[1, ]))
douek_formatted <- head(douek_formatted[-1,], 8)
rownames(douek_formatted) <- gsub("X", "Douek", rownames(douek_formatted))
## combine the two labs
douek_formatted$row_name <- rownames(douek_formatted)
sheba_formatted$row_name <- rownames(sheba_formatted)
hajj_formatted <- full_join(sheba_formatted, douek_formatted)
rownames(hajj_formatted) <- hajj_formatted$row_name
hajj_formatted <- mutate_all(hajj_formatted[, -1], function(x) (as.numeric(as.character(x))) ) %>%
    mutate(across(everything(), ~tidyr::replace_na(., 0)))
## create an element for each sample saying which lab it originated from
originLabList <- map(
    rownames(hajj_formatted),
    function(x) if(length(strsplit(x, "_")[[1]]) == 2) {"Douek"} else {"Sheba"}
)

## get species count (richness)
richness <- specnumber(hajj_formatted)
richFrame <- tibble(
    Douek = richness[originLabList == "Douek"],
    Sheba = richness[originLabList == "Sheba"]
)
richScatter <- ggplot(richFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, "_Hajj pathogen richness_", dataInfo)) +
    theme_bw()
ggsave(paste(dateToday, "_RichnessScatter_", dataInfo, ".png"))

richBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Richness=unlist(richness)),
        aes(x=Lab, y=Richness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_Hajj pathogen richness_", dataInfo)) +
        theme_bw()
ggsave(paste(dateToday, "_RichnessBoxPlots_", dataInfo, ".png"))

## get vegan diversity scores
divScores <- diversity(hajj_formatted, index = "shannon")
divFrame <- tibble(
    Douek = divScores[originLabList == "Douek"],
    Sheba = divScores[originLabList == "Sheba"]
)
divScatter <- ggplot(divFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, "_Hajj pathogen diversity_", dataInfo)) +
    theme_bw()
ggsave(paste(dateToday, "_DiversityScatter_", dataInfo, ".png"))

divBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Diversity=unlist(divScores)),
        aes(x=Lab, y=Diversity, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_Hajj pathogen diversity_", dataInfo)) +
        theme_bw()
ggsave(paste(dateToday, "_DiversityBoxPlots_", dataInfo, ".png"))

## get evenness score
evenScores <- divScores / log(richness)
evenFrame <- tibble(
    Douek = evenScores[originLabList == "Douek"],
    Sheba = evenScores[originLabList == "Sheba"]
)
evenScatter <- ggplot(evenFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title=paste(dateToday, "_Hajj pathogen evenness_", dataInfo)) +
    theme_bw()
ggsave(paste(dateToday, "_EvennessScatter_", dataInfo, ".png"))

evenBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Evenness=unlist(evenScores)),
        aes(x=Lab, y=Evenness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title=paste(dateToday, "_Hajj pathogen evenness_", dataInfo)) +
        theme_bw()
ggsave(paste(dateToday, "_EvennessBoxPlots_", dataInfo, ".png"))

## append summary stats to main file
hajjSummarized <- cbind(Lab=unlist(originLabList), richness=richness, diversity=divScores, evenness=evenScores, hajj_formatted)
write.csv(hajjSummarized, "hajjSummarized.csv")

# return to old working directory
setwd(old_wd)



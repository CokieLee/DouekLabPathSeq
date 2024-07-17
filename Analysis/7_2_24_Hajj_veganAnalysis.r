library("vegan")
library("dplyr")
library("ggplot2")
library("purrr")

# import csv file with species counts
old_wd <- getwd()
filePath <- "C:\\Users\\cokie\\OneDrive - Parker Tong LLP\\PersonalFiles\\Work\\DouekLab\\PathSeq\\HajjComparisons\\"
setwd(filePath)
hajj <- read.csv("RNA_species_tpm.csv")

# format data for vegan
hajj <- as.data.frame(hajj)
hajj_t <- data.frame(t(hajj))
colnames(hajj_t) <- as.character(unlist(hajj_t[1, ]))
hajj_t <- head(hajj_t[-1, ], -2)
hajj_t <- mutate_all(hajj_t, function(x) (as.numeric(as.character(x))) )
## create an element for each sample saying which lab it originated from
originLabList <- map(
    rownames(hajj_t),
    function(x) if(length(strsplit(x, "_")[[1]]) == 3) {"Sheba"} else {"Douek"}
)

## get species count (richness)
richness <- specnumber(hajj_t)
richFrame <- tibble(
    Douek = richness[originLabList == "Douek"],
    Sheba = richness[originLabList == "Sheba"]
)
richScatter <- ggplot(richFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title="Richness Sheba vs Douek") +
    theme_bw()
ggsave("RichnessScatter.png")

richBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Richness=unlist(richness)),
        aes(x=Lab, y=Richness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title="Hajj pathogen richness") +
        theme_bw()
ggsave("RichnessBoxPlots.png")

## get vegan diversity scores
divScores <- diversity(hajj_t, index = "shannon")
divFrame <- tibble(
    Douek = divScores[originLabList == "Douek"],
    Sheba = divScores[originLabList == "Sheba"]
)
divScatter <- ggplot(divFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title="Diversity Scores Sheba vs Douek") +
    theme_bw()
ggsave("DiversityScatter.png")

divBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Diversity=unlist(divScores)),
        aes(x=Lab, y=Diversity, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title="Hajj pathogen shannon diversity score") +
        theme_bw()
ggsave("DiversityBoxPlots.png")

## get evenness score
evenScores <- divScores / log(richness)
evenFrame <- tibble(
    Douek = evenScores[originLabList == "Douek"],
    Sheba = evenScores[originLabList == "Sheba"]
)
evenScatter <- ggplot(evenFrame, aes(x=Douek, y=Sheba)) +
    geom_point(color = "blue", size=3) +
    labs(title="Evenness Scores Sheba vs Douek") +
    theme_bw()
ggsave("EvennessScatter.png")

evenBoxes <-
    ggplot(
        data = data.frame(Lab=unlist(originLabList), Evenness=unlist(evenScores)),
        aes(x=Lab, y=Evenness, color=Lab),
        position = position_dodge(width=0.5)
        ) +
        geom_boxplot() +
        geom_point() +
        labs(title="Hajj pathogen pielou evenness score") +
        theme_bw()
ggsave("EvennessBoxPlots.png")

## append summary stats to main file
hajjSummarized <- cbind(Lab=unlist(originLabList), richness=richness, diversity=divScores, evenness=evenScores, hajj_t)
write.csv(hajjSummarized, "hajjSummarized.csv")

# return to old working directory
setwd(old_wd)

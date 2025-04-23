import sys
sys.path.append("/data/home/parkercol/PathSeq/Analysis/")
import maxSeparatingLine

import pandas as pd
import statsmodels.stats.multitest as smt

## adult pathseq unfiltered paths
adultPathseqUniltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment.csv"
CDFOutAdultPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propOutPathseqUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
## adult pathseq filtered paths
adultPathseqMeanFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutAdultPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
adultPathseqMedianFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutAdultPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
## adult czid unfiltered paths
adultCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment.csv"
CDFOutAdultCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propOutCZIDUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
## adult pathseq filtered paths
adultPathseqMeanFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutAdultPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
adultPathseqMedianFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutAdultPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"

## child pathseq unfiltered paths
childPathseqUniltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutChildPathseqUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Lactobacillus_childUnfiltered_threshold.png"
propOutChildPathseqUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Blautia_childUnfiltered_proportions.png"
## child pathseq filtered paths
childPathseqMeanFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutChildPathseqMeanFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
childPathseqMedianFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutChildPathseqMedianFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
## child czid unfiltered paths
childCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment.csv"
CDFOutChildCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propOutChildCZIDUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
## child pathseq filtered paths
childCZIDMeanFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutChildCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"
childPathseqMedianFiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/pathseqResults/summaryFormatData/Top10Genus_CountsAndTreatment_unfiltered.csv"
CDFOutChildCZIDUnfiltered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_threshold.png"
propPlotPathUnfiltertered = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Analysis/Gardnerella_unfiltered_proportions.png"

countsPath = childPathseqUniltered
plotOutPath = CDFOutChildPathseqUnfiltered
print(plotOutPath)
adultTreatA = "Positive"
adultTreatB = "Negative"
childTreatB = "HUU"
childTreatA = "HEU"
microbeLabel = "g_Blautia"

## read in counts table (of top 30 bacteria, with HIV status information)
counts = pd.read_csv(countsPath)
gardnerellaCounts = counts[ [microbeLabel, 'HIV.status'] ]
gardnerella_A = gardnerellaCounts[gardnerellaCounts['HIV.status'] == childTreatA][microbeLabel].to_list()
gardnerella_B = gardnerellaCounts[gardnerellaCounts['HIV.status'] == childTreatB][microbeLabel].to_list()

## get best threshold for separating positive and negative
gardnerellaThresh = maxSeparatingLine.bestThreshold(gardnerella_A, gardnerella_B, True,
                                                    "Cumulative distributions of "+ microbeLabel +" found by Pathseq")
print("Threshold: " + str(gardnerellaThresh))

## find the proportions of HIV + and - above and below the threshold
p_value = maxSeparatingLine.propSignificance(gardnerella_A, gardnerella_B, gardnerellaThresh)
p_value_bonferroni = maxSeparatingLine.bonferroni_correction(p_value, 10)

print(microbeLabel + " raw p-value: " + str(p_value))
print(microbeLabel + " bonferroni p-value: " + str(p_value_bonferroni))

## do significance test for all bacteria for children (based on HIV exposure status)
raw_p_values = []
for col in counts.columns:
    if col.startswith("g_"):
        bacteriaCounts = counts[ [col, 'HIV.status'] ]
        readsA = bacteriaCounts[bacteriaCounts['HIV.status'] == childTreatA][col].to_list()
        readsB = bacteriaCounts[bacteriaCounts['HIV.status'] == childTreatB][col].to_list()
        
        if (any(readsA) or any(readsB)):
            threshold = maxSeparatingLine.bestThreshold(readsA, readsB)
            p_value = maxSeparatingLine.propSignificance(readsA, readsB, threshold)

            raw_p_values.append(p_value)

reject, pvals_benjamini, _, _ = smt.multipletests(raw_p_values, alpha=0.05, method='fdr_bh')

## plot the proportions as barcharts with confidence intervals
maxSeparatingLine.visualizeProportions(
    gardnerella_A, gardnerella_B,
    gardnerellaThresh,
    "Proportion of samples with "+ microbeLabel +" above threshold " + str(gardnerellaThresh), childTreatA, childTreatB,
    propOutChildPathseqUnfiltertered)

print(raw_p_values)
print(pvals_benjamini)

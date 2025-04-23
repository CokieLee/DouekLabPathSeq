import numpy as np
from numpy import interp
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import statsmodels.stats.multitest as smt

def plotCDFThreshold(threshold, readsA, quantA, readsB, quantB, 
                     labelA, labelB, plotOutPath, title="Cumulative distribution of pathogen counts",):
    plt.figure()
    plt.scatter(quantA, readsA, color='red', label=labelA)
    plt.scatter(quantB, readsB, color='blue', label=labelB)
    plt.axhline(y=threshold, color='black', label="threshold of max separation (" + str(threshold) + ")")
    plt.title(title)
    plt.ylabel("pathogen counts (tpm) detected")
    plt.xlabel("Percentile of sample's pathogen counts")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotOutPath)

## function to draw maximally separating line
def bestThreshold(readsA, readsB, plot=False, plotTitle="Cumulative distribution of pathogen counts"):
    readsA = sorted(readsA)
    readsB = sorted(readsB)
    nA = len(readsA)
    nB = len(readsB)

    quantA = tuple( (k + (1/2)) / nA for k in range(nA) )
    quantB = tuple( (k + (1/2)) / nB for k in range(nB) )

    quantB_at_A = interp(readsA, readsB, quantB)
    quantA_at_B = interp(readsB, readsA, quantA)

    sep_at_A = tuple( ( abs(quantA[k] - quantB_at_A[k]), readsA[k]) for k in range(nA) if readsA[k] > 0)
    sep_at_B = tuple( ( abs(quantB[k] - quantA_at_B[k]), readsB[k]) for k in range(nB) if readsB[k] > 0)

    allSep = sep_at_A + sep_at_B
    maxSepDist, maxSepRead = max(allSep)

    if plot:
        plotCDFThreshold(maxSepRead, readsA, quantA, readsB, quantB, "HIV positive", "HIV negative", plotTitle)

    return( maxSepRead )

def getProportions(readsA, readsB, threshold):
    nA = len(readsA)
    A_aboveThresh = [ n for n in readsA if n > threshold ]
    nB = len(readsB)
    B_aboveThresh = [ n for n in readsB if n > threshold ]

    A_propAbove = len(A_aboveThresh) / nA
    B_propAbove = len(B_aboveThresh) / nB

    return nA, nB, len(A_aboveThresh), len(B_aboveThresh), A_propAbove, B_propAbove


def propSignificance(readsA, readsB, threshold):
    nA, nB, A_aboveThresh, B_aboveThresh, A_propAbove, B_propAbove = getProportions(readsA, readsB, threshold)

    ## do significance test for proportions
    p_pooled = (A_aboveThresh + B_aboveThresh) / (nA + nB)
    se_pooled = np.sqrt(p_pooled * (1 - p_pooled) * (1/nA + 1/nB))
    z_stat = (A_propAbove - B_propAbove) / se_pooled
    p_value = 2 * (1 - norm.cdf(abs(z_stat)))

    return p_value

def bonferroni_correction(p_value, num_comparisons, alpha=0.05):
    corrected_p_value = p_value * num_comparisons

    return corrected_p_value

def visualizeProportions(readsA, readsB, threshold, title="Proportion above threshold", A_label="group A", B_label="group B", pltPath="./proportionsBarChar.png"):
    nA, nB, A_aboveThresh, B_aboveThresh, A_propAbove, B_propAbove = getProportions(readsA, readsB, threshold)

    groups = [A_label, B_label]
    proportions = [A_propAbove, B_propAbove]

    confidence_level = 0.95
    z = norm.ppf( 1 - (1 - confidence_level) )

    se1 = np.sqrt(A_propAbove * (1 - B_propAbove) / nA)
    se2 = np.sqrt(B_propAbove * (1 - B_propAbove) / nB)

    ci1_lower = A_propAbove - z * se1
    ci1_upper = A_propAbove + z * se1
    ci2_lower = B_propAbove - z * se2
    ci2_upper = B_propAbove - z * se2

    ci_lower = [ci1_lower, ci2_lower]
    ci_upper = [ci1_upper, ci2_upper]
    ci_errors = [[abs(proportions[i] - ci_lower[i]), abs(ci_upper[i] - proportions[i])] for i in range(2)]

    plt.figure()
    plt.bar(groups, proportions, color=['red', 'blue'], yerr=np.array(ci_errors).T, capsize=5)
    plt.ylabel("Proportion above threshold (95% CI)")
    plt.title(title)
    
    plt.savefig(pltPath)
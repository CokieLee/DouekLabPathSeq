import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import pandas as pd
import os

## plot the size of the Input files
infoFile = open("/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/InputFilesInfo.txt")
outputDir = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Profiling/"

####################################################################
## This code is for producing plots which visualize various performance
# aspects of the pathseq snakemake pipeline, across muliple input sample runs.
# It produces visualizations of the memory and compute time performance for each rule,
# and visualizes the relationship between performance and traits of the input data of each rule
# (including input size, and in some cases, input structure).

## Assumptions about input data and directory structures

# Snakemake log files
# Sacct file

####################################################################


sizes = []

for line in infoFile.readlines()[1:]:
    size = int(line.split()[4]) / 1000000000    # size as gigabytes
    sizes.append(size)

counts, edges = np.histogram(sizes, bins=50)
plt.stairs(counts, edges)

plt.xlabel("Input file size (in Gb)")
plt.ylabel("Number of files")
plt.title("Krystelle's microbiome data (1Gb ~= 15M reads)")

plt.savefig("InputSizeHistogram.jpg")

####################################################################
# concatenate sacct files
logsPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/"
bowtieSacctPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/bowtieUnmaskedSacct.txt"
sacctPath3 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_3_25_sacct_formatted.txt"
sacctPath5 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_5_25_sacct_formatted.txt"
sacctPath11 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_11_25_sacct_formatted.txt"
sacctPath12 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_12_25_sacct_formatted.txt"
sacctPath16 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_16_25_sacct_formatted.txt"

column_widths3 = [13, 11, 11, 11, 11, 9]
column_widths5 = [13, 11, 11, 11]

sacctBowtie = pd.read_fwf(bowtieSacctPath, widths=column_widths5).drop(index=0)
sacct3 = pd.read_fwf(sacctPath3, widths=column_widths3).drop(index=0)
sacct5 = pd.read_fwf(sacctPath5, widths=column_widths5).drop(index=0)
sacct11 = pd.read_fwf(sacctPath11, widths=column_widths3).drop(index=0)
sacct12 = pd.read_fwf(sacctPath12, widths=column_widths3).drop(index=0)
sacct16 = pd.read_fwf(sacctPath16, widths=column_widths3).drop(index=0)

sacctMerge = pd.concat([sacctBowtie, sacct3, sacct5, sacct11, sacct12, sacct16], ignore_index=True)
sacctMerge['JobIDBase'] = sacctMerge['JobID'].apply(lambda x: x.split(".")[0])

# get list of successful and relevant job IDs to filter sacct list
filenames = [ f for f in os.listdir(logsPath) if ( os.path.isfile(os.path.join(logsPath, f)) and f.endswith(".out") ) ]
jobIDs = [(f.split(".")[0]).split("-")[1] for f in filenames]

uniqueSacct = sacctMerge[ ( sacctMerge['JobIDBase'].isin(jobIDs) ) ]

uniqueSacct.to_csv(os.path.join(logsPath, "2_3_25_sacctSummary.txt"), sep="\t", index=False)

sacctPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_3_25_sacctSummary.txt"
sacct = pd.read_csv(sacctPath, sep='\t')

##################################################################
## function for getting the full path of the input file used by a snakemake rule
# given a snakemake rule output log
def getInputFile(filename, inputNumber):
    with open(filename, "r") as file:
        line = file.readline()
        while line:
            if line.startswith("    input:"):
                inputPath = (line.split(",")[inputNumber]).strip()
                return(inputPath)
            else:
                line = file.readline()

## function for finding the CPU time (in hours) for a batch job, given the job's jobID, and a sacct file.
#   The sacct file is a dataframe, which much contain at least the columns "JobID" (a column with a number,
#   and an optional suffix. the ".bat+" suffix indicates entry for the batch job (rather than the job launching overhead)),
#   "JobIDBase" (job id number, without any suffix), and "CPUTime" (time spent computing for the job).
def getCPUTime(sacctSummary, jobID):
    filterForJobs = sacctSummary[sacctSummary['JobIDBase'] == jobID]
    batchEntry = filterForJobs[filterForJobs['JobID'].str.endswith(".bat+")]
    if len(batchEntry) == 0:
        cpuTime = 0
        return(cpuTime)
    else:
        cpuTime = batchEntry.iloc[0, 2].split("-")
        if (len(cpuTime) > 1):
            cpuTimePT = datetime.strptime(cpuTime[1], '%H:%M:%S')
            days = float(cpuTime[0])
        elif (len(cpuTime) == 1):
            cpuTimePT = datetime.strptime(cpuTime[0], '%H:%M:%S')
            days = 0
        cpuSeconds = float(cpuTimePT.second) + float(60*cpuTimePT.minute) + float(3600*cpuTimePT.hour) + float(3600 * 24 * days)
        cpuHours = (float(cpuSeconds) / float(3600))
        return(cpuHours)

## function for finding the maxRSS (max memory usage) for a batch job, given the job's jobID, and a sacct file.
#   The sacct file is a dataframe, which much contain at least the columns "JobID" (a column with a number,
#   and an optional suffix. the ".bat+" suffix indicates entry for the batch job (rather than the job launching overhead)),
#   "JobIDBase" (job id number, without any suffix), and "MaxRSS" (maximum memory used by the job entry).
def getMaxRSS(sacctSummary, jobID):
    filterForJobs = sacctSummary[sacctSummary['JobIDBase'] == jobID]
    batchEntry = filterForJobs[filterForJobs['JobID'].str.endswith(".bat+")]
    if len(batchEntry) == 0:
        maxRSS = 0
    else:
        maxRSS = (float("nan") if isinstance(batchEntry.iloc[0, 3], float) else (batchEntry.iloc[0, 3])[:-1])
    return(maxRSS)

## function for gathering resource usage information from sacct file for particular job name. 
# INPUT: sacct (df with columns: JobID, JobName, CPUTime, MaxRSS, State, NTasks, JobIDBase),
#   logsPath (path string to directory containing directories (named for each job/rule name),
#       containing snakemake-style rule log files for each job, and whose filenames are of the form [rule name].[job ID].out),
#   ruleName (job name from sacct's "JobName" column)
# OUTPUT: csv file of resource usage information for one particular job name
#   (columns: jobID, inputSize, CPUTime, maxRSS)
def compileRuleInfo(sacct, logsPath, ruleName):
    ruleOuts = [ f for f in os.listdir(os.path.join(logsPath, ruleName)) if f.endswith(".out")]
    ruleJobIDs = [ (f.split("-")[-1]).split(".")[0] for f in ruleOuts ]
    ruleInfo = pd.DataFrame({'jobID': ruleJobIDs})
    ruleInput = [ getInputFile( os.path.join(logsPath, ruleName, f), 1 ) for f in ruleOuts ]

    # get file size for each input file, populate into a column of dataframe
    ruleInfo['inputSize'] = [ os.path.getsize(f) for f in ruleInput ]
    ruleInfo['inputSize'] = ruleInfo['inputSize']
    ruleInfo['CPUTime'] = [ ( getCPUTime(sacct, float(id)) ) for id in ruleInfo['jobID'] ]
    ruleInfo['maxRSS'] = [ ( getMaxRSS(sacct, float(id)) ) for id in ruleInfo['jobID'] ]

    return(ruleInfo)

## function to make 2 scatter plots: one of file's size (x-axis) against the max memory usage
#   of a program which used the file as an input. one of file's size against program time.
# Input: dataframe with at least columns "inputSize", "CPUTime", and "maxRSS"
def plotTimeAndMem(ruleInfo, ruleName):
    ruleInfoTime = ruleInfo.dropna(subset=['CPUTime'])
    TimeFig, TimeAx = plt.subplots()
    TimeAx.scatter(ruleInfoTime['inputSize'].astype(int), ruleInfoTime['CPUTime'].astype(float))
    TimeAx.set_title(ruleName + " on Skyline (89 samples from Krystelle's Microbiome data)")
    TimeAx.set_xlabel("Input file size (1GB ~= 15M reads)")
    TimeAx.set_ylabel("CPUTime in hours")
    TimeAx.set_xticks(np.arange(0, ( (7 + 1) + 1000000000), 1000000000))
    TimeFig.savefig(os.path.join(outputDir, ruleName + "_inputSize_vs_CPUTime.png"))

    # # plot scatter plot of input size against memory usage
    ruleInfoMem = ruleInfo.dropna(subset=['maxRSS'])
    MemFig, MemAx = plt.subplots()
    MemAx.scatter(ruleInfoMem['inputSize'].astype(int), ruleInfoMem['maxRSS'].astype(float))
    MemAx.set_title(ruleName + " on Skyline (89 samples from Krystelle's Microbiome data)")
    MemAx.set_xlabel("Input file size (1GB ~= 15M reads)")
    MemAx.set_ylabel("Max memory used in kb")
    MemFig.savefig(os.path.join(outputDir, ruleName + "_inputSize_vs_maxRSS.png"))

## function which takes in a list of unsorted numbers (input sizes for a particular rule),
# and produces a distribution plot (scatter plot of list, sorted).
# ruleName is a string, used to title the plot
def plotCDF(unsortedNumbers, ruleName):
    sortedList = np.sort(unsortedNumbers)
    n = len(sortedList)

    fig, ax = plt.subplots()
    ax.scatter(range(n), sortedList)
    ax.set_title(ruleName + " input size cdf")
    ax.set_xlabel("input size percentile")
    ax.set_ylabel("input size in GB")
    fig.savefig(os.path.join(outputDir, ruleName + "_inputSize_CDF.png"))

###################################################################
## plot the time taken and memory used for bowtieUnmasked against input file size
bowtieUnmaskedInfo = compileRuleInfo(sacct, logsPath, "bowtieUnmasked")
plotTimeAndMem(bowtieUnmaskedInfo, "bowtieUnmasked")
# plot CDF of all input files for bowtieUnmasked
plotCDF(bowtieUnmaskedInfo['inputSize'].to_numpy(), "bowtieUnmasked")

## plot the time taken for star against input file size
starInfo = compileRuleInfo(sacct, logsPath, "star")
plotTimeAndMem(starInfo, "star")
# plot CDF of all input files for star
plotCDF(starInfo['inputSize'].to_numpy(), "star")

## plot the time taken for bowtiePrimate against input file size
bowtiePrimateInfo = compileRuleInfo(sacct, logsPath, "bowtiePrimate")
plotTimeAndMem(bowtiePrimateInfo, "bowtiePrimate")
# plot CDF of all input files for bowtiePrimate
plotCDF(bowtiePrimateInfo['inputSize'].to_numpy(), "bowtiePrimate")

## plot the time taken for trinity against post-filtering size, and contig size
trinityInfo = compileRuleInfo(sacct, logsPath, "trinity")
plotTimeAndMem(trinityInfo, "trinity")
# plot CDF of all input files for trinity
plotCDF(trinityInfo['inputSize'].to_numpy(), "trinity")
# plot contig size CDF for 2 lowest, 2 highest, and 2 median time trinity run

# add averge and largest contig size to trinity info, and plot against this


## plot the time taken for filterHost against input size
filterHostInfo = compileRuleInfo(sacct, logsPath, "filterHost")
plotTimeAndMem(filterHostInfo, "filterHost")
# TODO: plot the time taken for filterHost against contig size of inputs
# plot CDF of all input files for filterHost
plotCDF(filterHostInfo['inputSize'].to_numpy(), "filterHost")

## plot the time taken for kaiju against input size and contig size
## TODO: modify the kaiju rule so that the .fasta input comes immediately after the script input

## plot the time taken for buildSalmon against input size and contig size

## plot the time taken for salmonQuant against input size and contig size

## plot the time taken for mergeTaxAndQuant against input size and contig size
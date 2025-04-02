import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import pandas as pd
import os

## plot the size of the Input files
infoFile = open("/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/InputFilesInfo.txt")
outputDir = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/Visualizations/Profiling/"

sizes = []

for line in infoFile.readlines()[1:]:
    size = int(line.split()[4]) / 1000000000    # size as gigabytes
    sizes.append(size)

counts, edges = np.histogram(sizes, bins=50)
plt.stairs(counts, edges)

plt.xlabel("Input file size (in Gb)")
plt.ylabel("Number of files")
plt.title("Microbiome sequencing data (1Gb ~= 15M reads)")

plt.savefig(os.path.join(outputDir, "InputSizeHistogram.jpg") )
####################################################################
## plot input size of bowtieUnmasked

# sacctFile = open("/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/bowtieUnmaskedSacct.txt")

# times = []
# mems = []

# for line in sacctFile.readlines()[2:]:
#     if "batch" in line:
#         entries = line.split()

#         if len(entries) >= 4:
#             time = entries[2]
#             mem = entries[3]
#             times.append(time)
#             mems.append(mem)

# def clockToSeconds(clockString):
#     totalSecs = 0
#     h,m,s = clockString.split(":")
#     totalSecs = totalSecs + (int(h) * 60 * 60)
#     totalSecs = totalSecs + (int(m) * 60)
#     totalSecs = totalSecs + int(s)

#     return totalSecs

# timeInSecs = list(map(clockToSeconds, times))

# counts, edges = np.histogram(timeInSecs, bins=100)
# plt.stairs(counts, edges)

# plt.xlabel("Time to run bowtieUnmasked.sh (in Secs)")
# plt.ylabel("Number of files")
# plt.title("Krystelle's microbiome data")

# plt.savefig("BowtieUnmaskedRunTimeHistogram.jpg")

##################################################################
# concatenate sacct files
# logsPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/"
# bowtieSacctPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/bowtieUnmaskedSacct.txt"
# sacctPath3 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_3_25_sacct_formatted.txt"
# sacctPath5 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_5_25_sacct_formatted.txt"
# sacctPath11 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_11_25_sacct_formatted.txt"
# sacctPath12 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_12_25_sacct_formatted.txt"
# sacctPath16 = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_16_25_sacct_formatted.txt"

# column_widths3 = [13, 11, 11, 11, 11, 9]
# column_widths5 = [13, 11, 11, 11]

# sacctBowtie = pd.read_fwf(bowtieSacctPath, widths=column_widths5).drop(index=0)
# sacct3 = pd.read_fwf(sacctPath3, widths=column_widths3).drop(index=0)
# sacct5 = pd.read_fwf(sacctPath5, widths=column_widths5).drop(index=0)
# sacct11 = pd.read_fwf(sacctPath11, widths=column_widths3).drop(index=0)
# sacct12 = pd.read_fwf(sacctPath12, widths=column_widths3).drop(index=0)
# sacct16 = pd.read_fwf(sacctPath16, widths=column_widths3).drop(index=0)

# sacctMerge = pd.concat([sacctBowtie, sacct3, sacct5, sacct11, sacct12, sacct16], ignore_index=True)
# sacctMerge['JobIDBase'] = sacctMerge['JobID'].apply(lambda x: x.split(".")[0])

# # get list of successful and relevant job IDs to filter sacct list
# filenames = [ f for f in os.listdir(logsPath) if ( os.path.isfile(os.path.join(logsPath, f)) and f.endswith(".out") ) ]
# jobIDs = [(f.split(".")[0]).split("-")[1] for f in filenames]

# uniqueSacct = sacctMerge[ ( sacctMerge['JobIDBase'].isin(jobIDs) ) ]

# uniqueSacct.to_csv(os.path.join(logsPath, "2_3_25_sacctSummary.txt"), sep="\t", index=False)

# sacctPath = "/data/vrc_his/douek_lab/projects/PathSeq/Krystelle/logs_successful_2_3_25/2_3_25_sacctSummary.txt"
# sacct = pd.read_csv(sacctPath, sep='\t')

# ##################################################################
# def getInputFile(filename, inputNumber):
#     with open(filename, "r") as file:
#         line = file.readline()
#         while line:
#             if line.startswith("    input:"):
#                 inputPath = (line.split(",")[inputNumber]).strip()
#                 return(inputPath)
#             else:
#                 line = file.readline()

# # lookup the cputime and memory usage for batch job with each ID, store into 2 more columns
# def getCPUTime(sacctSummary, jobID):
#     filterForJobs = sacctSummary[sacctSummary['JobIDBase'] == jobID]
#     batchEntry = filterForJobs[filterForJobs['JobID'].str.endswith(".bat+")]
#     if len(batchEntry) == 0:
#         cpuTime = 0
#         return(cpuTime)
#     else:
#         cpuTime = batchEntry.iloc[0, 2].split("-")
#         if (len(cpuTime) > 1):
#             cpuTimePT = datetime.strptime(cpuTime[1], '%H:%M:%S')
#             days = float(cpuTime[0])
#         elif (len(cpuTime) == 1):
#             cpuTimePT = datetime.strptime(cpuTime[0], '%H:%M:%S')
#             days = 0
#         cpuSeconds = float(cpuTimePT.second) + float(60*cpuTimePT.minute) + float(3600*cpuTimePT.hour) + float(3600 * 24 * days)
#         cpuHours = (float(cpuSeconds) / float(3600))
#         return(cpuHours)

# def getMaxRSS(sacctSummary, jobID):
#     filterForJobs = sacctSummary[sacctSummary['JobIDBase'] == jobID]
#     batchEntry = filterForJobs[filterForJobs['JobID'].str.endswith(".bat+")]
#     if len(batchEntry) == 0:
#         maxRSS = 0
#     else:
#         maxRSS = (float("nan") if isinstance(batchEntry.iloc[0, 3], float) else (batchEntry.iloc[0, 3])[:-1])
#     return(maxRSS)

# def compileRuleInfo(sacct, logsPath, ruleName):
#     ruleOut = [ f for f in os.listdir(os.path.join(logsPath, ruleName)) if f.endswith(".out")]
#     ruleJobID = [ (f.split("-")[-1]).split(".")[0] for f in ruleOut ]
#     ruleInfo = pd.DataFrame({'jobID': ruleJobID})
#     ruleInput = [ getInputFile( os.path.join(logsPath, ruleName, f), 1 ) for f in ruleOut ]

#     # get file size for each input file, populate into a column of dataframe
#     ruleInfo['inputSize'] = [ os.path.getsize(f) for f in ruleInput ]
#     ruleInfo['inputSize'] = ruleInfo['inputSize']
#     ruleInfo['CPUTime'] = [ ( getCPUTime(sacct, float(id)) ) for id in ruleInfo['jobID'] ]
#     ruleInfo['maxRSS'] = [ ( getMaxRSS(sacct, float(id)) ) for id in ruleInfo['jobID'] ]

#     return(ruleInfo)

# def plotTimeAndMem(ruleInfo, ruleName):
#     ruleInfoTime = ruleInfo.dropna(subset=['CPUTime'])
#     TimeFig, TimeAx = plt.subplots()
#     TimeAx.scatter(ruleInfoTime['inputSize'].astype(int), ruleInfoTime['CPUTime'].astype(float))
#     TimeAx.set_title(ruleName + " on Skyline (89 samples from Krystelle's Microbiome data)")
#     TimeAx.set_xlabel("Input file size (1GB ~= 15M reads)")
#     TimeAx.set_ylabel("CPUTime in hours")
#     TimeAx.set_xticks(np.arange(0, ( (7 + 1) + 1000000000), 1000000000))
#     TimeFig.savefig(os.path.join(outputDir, ruleName + "_inputSize_vs_CPUTime.png"))

#     # # plot scatter plot of input size against memory usage
#     ruleInfoMem = ruleInfo.dropna(subset=['maxRSS'])
#     MemFig, MemAx = plt.subplots()
#     MemAx.scatter(ruleInfoMem['inputSize'].astype(int), ruleInfoMem['maxRSS'].astype(float))
#     MemAx.set_title(ruleName + " on Skyline (89 samples from Krystelle's Microbiome data)")
#     MemAx.set_xlabel("Input file size (1GB ~= 15M reads)")
#     MemAx.set_ylabel("Max memory used in kb")
#     MemFig.savefig(os.path.join(outputDir, ruleName + "_inputSize_vs_maxRSS.png"))

# def plotCDF(unsortedNumbers, ruleName):
#     sortedList = np.sort(unsortedNumbers)
#     n = len(sortedList)

#     fig, ax = plt.subplots()
#     ax.scatter(range(n), sortedList)
#     ax.set_title(ruleName + " input size cdf")
#     ax.set_xlabel("input size percentile")
#     ax.set_ylabel("input size in GB")
#     fig.savefig(os.path.join(outputDir, ruleName + "_inputSize_CDF.png"))

# ###################################################################
# ## plot the time taken and memory used for bowtieUnmasked against input file size
# bowtieUnmaskedInfo = compileRuleInfo(sacct, logsPath, "bowtieUnmasked")
# plotTimeAndMem(bowtieUnmaskedInfo, "bowtieUnmasked")
# # plot CDF of all input files for bowtieUnmasked
# plotCDF(bowtieUnmaskedInfo['inputSize'].to_numpy(), "bowtieUnmasked")

# ## plot the time taken for star against input file size
# starInfo = compileRuleInfo(sacct, logsPath, "star")
# plotTimeAndMem(starInfo, "star")
# # plot CDF of all input files for star
# plotCDF(starInfo['inputSize'].to_numpy(), "star")

# ## plot the time taken for bowtiePrimate against input file size
# bowtiePrimateInfo = compileRuleInfo(sacct, logsPath, "bowtiePrimate")
# plotTimeAndMem(bowtiePrimateInfo, "bowtiePrimate")
# # plot CDF of all input files for bowtiePrimate
# plotCDF(bowtiePrimateInfo['inputSize'].to_numpy(), "bowtiePrimate")

# ## plot the time taken for trinity against post-filtering size, and contig size
# trinityInfo = compileRuleInfo(sacct, logsPath, "trinity")
# plotTimeAndMem(trinityInfo, "trinity")
# # plot CDF of all input files for trinity
# plotCDF(trinityInfo['inputSize'].to_numpy(), "trinity")
# # plot contig size CDF for 2 lowest, 2 highest, and 2 median time trinity run

# # add averge and largest contig size to trinity info, and plot against this


# ## plot the time taken for filterHost against input size and contig size

# ## plot the time taken for kaiju against input size and contig size

# ## plot the time taken for buildSalmon against input size and contig size

# ## plot the time taken for salmonQuant against input size and contig size

# ## plot the time taken for mergeTaxAndQuant against input size and contig size
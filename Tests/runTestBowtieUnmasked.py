import sys, getopt

def compareBowtieOutputs(bowtie1, bowtie2):
    bowtie1File = open(bowtie1, "r")
    bowtie2File = open(bowtie2, "r")

    bowtie1Lines = bowtie1File.readlines()
    bowtie2Lines = bowtie2File.readlines()
    for line1, line2 in zip(bowtie1Lines, bowtie2Lines):
        bowtie1Val = line1.split()[0]
        bowtie2Val = line2.split()[0]
        if bowtie1Val != bowtie2Val:
            return (False, line1 + " != " + line2)

    return (True, None)

def main(argv):
    truthDir = None
    testDir = None

    opts, arg = getopt.getopt(argv, "h:t:i:", ["truthDir=", "testDir="])
    for opt, arg in opts:
        if opt == '-h':
            print('Program runs the bowtieUnmasked step of pathseq pipeline, and then tests that the output is correct')
        if opt == '-t':
            truthDir = arg
            print('Ground truth file path is ', truthDir)
        if opt == '-i':
            testDir = arg
            print('Output to be tested is at file path ', testDir)

    if truthDir == None:
        print('Error: must provide path to output of ground truth bowtieUnmasked alignment')
        exit(1)
    if testDir == None:
        print('Error: must provide path to output of test of bowtieUnmasked alignment')
        exit(1)

    if truthDir.endswith('/'):
        truthDir = truthDir[:-1]
    if testDir.endswith('/'):
        testDir = testDir[:-1]

    sampleNameTruth = (truthDir.split("/")[-1]).split(".")[0]
    sampleNameTest = (truthDir.split("/")[-1]).split(".")[0]

    if sampleNameTruth != sampleNameTest:
        print("Error: mismatching sampleOutput folders")
        exit(1)
    sampleNameBase = sampleNameTruth[len("Sample_"):]

    # get the alignment outputs for the first and second alignment for test and truth
    truthBowtieStats = []
    testBowtieStats = []

    truthBowtieStats.append(truthDir + "/ERCC_alignment_rates/ERCC_alignment_rate_" + sampleNameBase + ".txt")
    truthBowtieStats.append(truthDir + "/Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome/genome_alignment_rate_" + sampleNameBase + ".txt")

    testBowtieStats.append(testDir + "/ERCC_alignment_rates/ERCC_alignment_rate_" + sampleNameBase + ".txt")
    testBowtieStats.append(testDir + "/Generated_Data_2nd_Bowtie_Alignment_Unmasked_Genome/genome_alignment_rate_" + sampleNameBase + ".txt")

    if len(truthBowtieStats) != len(testBowtieStats):
        print("Error: differing number of sample outputs for test and ground truth")
        exit(1)

    for truth, test in zip(truthBowtieStats, testBowtieStats):
        compareResults = compareBowtieOutputs(test, truth)
        if compareResults[0]:
            print("bowtieUnmasked first alignment correct for ", sampleNameTest)
        else:
            print("bowtieUnmasked first alignment has mismatch for ", sampleNameTest, ": ", compareResults[1])
    

if __name__ == "__main__":
    main(sys.argv[1:])

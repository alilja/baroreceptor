import csv

def baroAnalysis(fileName, clusterWidth = 3, lag = 1):
    f = open(fileName,"r")
    reader = csv.reader(f)

    HR = [0]
    systolic = [0]

    lineNum = 0
    for line in reader:
        if(lineNum > 1): #skip the headers
            if(line[5] != HR[-1]):
                HR.append(line[5])
            if(line[4] != systolic[-1]):
                systolic.append(line[4])

        lineNum += 1
    f.close()

    clusters = []
    testNum = 0
    currentClusterStart = 0
    currentCluserEnd = 0

    for i in range(0, len(systolic) - clusterWidth):
        if(i <= len(HR) - (lag + clusterWidth)):


            currentSysCluster = systolic[i:i+clusterWidth]
            currentHRCluster = HR[i+lag:i+lag+clusterWidth]
            if(descending(currentHRCluster) and descending(currentSysCluster)):
                #both descending
            if(ascending(currentHRCluster) and ascending(currentSysCluster)):
                #both ascending





baroAnalysis("testData.csv")
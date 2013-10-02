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
                HR.append(float(line[5]))
            if(line[4] != systolic[-1]):
                systolic.append(float(line[2]))

        lineNum += 1
    f.close()

    clusters = []
    currentSysDirection = ""
    currentHRDirection = ""

    currentSysRun = [0]
    currentHRRun = [0]

    for i in range(0, len(systolic)):
        if(i <= len(HR) - (lag + clusterWidth)):
            if(currentSysDirection == "+" or currentSysDirection == ""):
                if(systolic[i] >= currentSysRun[-1]):
                    currentSysRun.append(systolic[i]) #might have to grab beginning and ending i instead to make it play nice with HR
                    currentSysDirection = "+"
                else: #end of run
                    if(len(currentSysRun) >= clusterWidth):
                        clusters.append(currentSysRun)
                    currentSysRun = [0]
                    currentSysDirection = ""

            if(currentSysDirection == "-" or currentSysDirection == ""):
                if(systolic[i] <= currentSysRun[-1]):
                    currentSysRun.append(systolic[i])
                    currentSysDirection = "-"
                else: #end of run
                    if(len(currentSysRun) >= clusterWidth):
                        clusters.append(currentSysRun)
                    currentSysRun = [0]
                    currentSysDirection = ""

    print(clusters[0:10])









baroAnalysis("testData.csv")
import csv

#readCSVFile: str, int, int, int --> tuple-of-list-of-num
def readCSVFile(fileName, headerLength = 1, HRChannel = 42, SBPChannel = 40):
    f = open(fileName,"r")
    reader = csv.reader(f,delimiter="\t")

    HR = [0]
    SBP = [0]

    HRIndex = 0
    SBPIndex = 0

    lineNum = 0

    runs = []
    currentRun = {"HR":[],"SBP":[]}
    direction = 0

    for line in reader:
        if("CH"+str(HRChannel) in line):
            HRIndex = line.index("CH"+str(HRChannel))
            SBPIndex = line.index("CH"+str(SBPChannel))
            print("HR index: "+str(HRIndex))
            print("SBPIndex: "+str(SBPIndex))

        if(lineNum > headerLength): #skip the headers
            if(float(line[HRIndex]) != HR[-1]):
                HR.append(float(line[HRIndex]))
            if(float(line[SBPIndex]) != SBP[-1]):
                SBP.append(float(line[SBPIndex]))

        lineNum += 1
    f.close()
    print("Finished analyzing file \""+fileName+"\"")
    return (SBP, HR)  

## Combine these two functions so that all the data processing is done online ##
## Going to need two items at a time. Good luck & godspeed. ##

#findMatchingRuns: list-of-num, list-of-num, num, num --> list-of-list-of-num
def findMatchingRuns(SBP, HR, clusterWidth = 3, lag = 0):
    """Takes two lists of numbers and determines when they are both moving
    in the same direction; that is, when they are both increasing at the same
    time or both decreasing at the same time. Example:

        a = [0,1,2, 3,4,5,1,0,9, 8]
        b = [1,4,9,13,0,3,2,0,11,14]

    The function will find the indices 0:3, 4:5, 5:7, and 7:8. It will also
    determine the length and direction. Output is in the form of lists:

        [runStart, runEnd, runLength, runDirection]

    Direction is either +1 or -1.

    clusterWidth is the minimum length of each run. lag is the difference
    in offset between the second list and the first list."""

    print(len(HR))
    print(len(SBP))

    runs = []
    currentRun = {}
    direction = 0
    for i in range(1, len(SBP)):
        SBPDiff = SBP[i] - SBP[i - 1]
        HRDiff = HR[i + lag] - HR[i + lag - 1]

        if(HRDiff > 0 and SBPDiff > 0):
            if(direction >= 0):
                currentRun["HR"].append(HR[i])
                currentRun["SBP"].append(SBP[i])
                direction = 1
            else:
                runs.append(currentRun)
                currentRun = {"HR":[HR[i]], "SBP":[SBP[i]]}
                direction = 1
        elif(HRDiff < 0 and SBPDiff < 0):
            if(direction <= 0):
                currentRun["HR"].append(HR[i])
                currentRun["SBP"].append(SBP[i])
                direction = -1
            else:
                #print(currentRun)
                runs.append(currentRun)
                currentRun = {"HR":[HR[i]], "SBP":[SBP[i]]}
                direction = -1
        elif(HRDiff == 0 and SBPDiff == 0):
            currentRun["HR"].append(HR[i])
            currentRun["SBP"].append(SBP[i])
        else:
            #print(currentRun)
            runs.append(currentRun)
            currentRun = {"HR":[HR[i]], "SBP":[SBP[i]]}

    return [r for r in runs[1:-1] if len(r["SBP"]) >= clusterWidth]



data = readCSVFile("davis cold pressor0000.csv", 32)
runs = findMatchingRuns(data[0],data[1], 2, 1)

print(runs[0:10])

for run in runs:
    runStart = run[0]
    runEnd = run[1]
    runLength = run[2]
    runDirection = run[3]

    for i in range(runStart - 3, runEnd + 3):
        out = ""
        if(i == runStart or i == runEnd):
            out += "> "
        out += str(data[0][i]) + "    " + str(data[1][i])
        #print(out)

    #print("\n\n")
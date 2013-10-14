import csv

def pearsonR(x, y):
    # Assume len(x) == len(y)
    n = len(x)
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))
    psum = sum(map(lambda x, y: x * y, x, y))
    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
    if den == 0: return 0
    return num / den

#readCSVFile: str, int, int, int --> tuple-of-list-of-num
def readCSVFile(fileName, headerLength = 1, RRChannel = 42, SBPChannel = 40):
    f = open(fileName,"r")
    reader = csv.reader(f,delimiter="\t")

    RR = [0]
    SBP = [0]

    RRIndex = 0
    SBPIndex = 0

    lineNum = 0

    runs = []
    currentRun = {"RR":[],"SBP":[]}
    direction = 0

    for line in reader:
        if("CH"+str(RRChannel) in line):
            RRIndex = line.index("CH"+str(RRChannel))
            SBPIndex = line.index("CH"+str(SBPChannel))
            print("RRIndex: "+str(RRIndex))
            print("SBPIndex: "+str(SBPIndex))

        if(lineNum > headerLength): #skip the headers
            if(float(line[RRIndex]) != RR[-1]):
                RR.append(float(line[RRIndex]))
            if(float(line[SBPIndex]) != SBP[-1]):
                SBP.append(float(line[SBPIndex]))

        lineNum += 1
    f.close()
    print("Finished analyzing file \""+fileName+"\"")
    return (SBP, RR)  

## Combine these two functions so that all the data processing is done online ##
## Going to need two items at a time. Good luck & godspeed. ##

#findMatchingRuns: list-of-num, list-of-num, num, num --> list-of-list-of-num
def findMatchingRuns(SBP, RR, clusterWidth = 3, lag = 0):
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

    print("RR length: "+str(len(RR)))
    print("SBP length: "+str(len(SBP)))

    runs = []
    currentRun = {"RR":[],"SBP":[]}
    direction = 0
    for i in range(1, len(SBP)):
        SBPDiff = SBP[i] - SBP[i - 1]
        RRDiff = RR[i + lag] - RR[i + lag - 1]

        if(RRDiff > 0 and SBPDiff > 0):
            if(direction >= 0):
                currentRun["RR"].append(RR[i])
                currentRun["SBP"].append(SBP[i])
                direction = 1
            else:
                runs.append(currentRun)
                currentRun = {"RR":[RR[i]], "SBP":[SBP[i]]}
                direction = 1
        elif(RRDiff < 0 and SBPDiff < 0):
            if(direction <= 0):
                currentRun["RR"].append(RR[i])
                currentRun["SBP"].append(SBP[i])
                direction = -1
            else:
                #print(currentRun)
                runs.append(currentRun)
                currentRun = {"RR":[RR[i]], "SBP":[SBP[i]]}
                direction = -1
        elif(RRDiff == 0 and SBPDiff == 0):
            currentRun["RR"].append(RR[i])
            currentRun["SBP"].append(SBP[i])
        else:
            #print(currentRun)
            runs.append(currentRun)
            currentRun = {"RR":[RR[i]], "SBP":[SBP[i]]}
    return [r for r in runs[1:-1] if len(r["SBP"]) >= clusterWidth]

#findCorrelatedRuns: list-of-dict-of-list-of-num, num --> list-of-dict-of-list-of-num
def findCorrelatedRuns(runs, minCorrelation = 0.85):
    return [run for run in runs if pearsonR(run["SBP"], run["RR"]) > minCorrelation]


data = readCSVFile("davis cold pressor0000.csv", 32)

f = open("davisColdPressorRR.csv")
reader = csv.reader(f)
lineNum = 0
davisRR = []
for line in reader:
    if(lineNum > 3):
        davisRR.append(float(line[0]))
    lineNum += 1
f.close()

runs = findMatchingRuns(data[0],davisRR, 3, 1)
correlatedRuns = findCorrelatedRuns(runs, .95)

print(correlatedRuns)

f = open("CSVOutput.csv","w")
output = "SBP, RR\n"
for run in correlatedRuns:
    zippedPairs = zip(run["SBP"], run["RR"])
    for pair in zippedPairs:
        output += str(pair[0]) + ", " + str(pair[1]) + "\n"
    output += ",\n"
f.write(output)
f.close()
#!/usr/bin/python3
import csv
import sys
import getopt
import os

_verbose = False
_fileName = ""
_headerLength = 33
_RR = "CH42"
_SBP = "CH40"
_ECG = "CH14"
_filter = 1.5
_pearson = 0.85
_width = 3
_lag = 0

try:
    opts, args = getopt.getopt(sys.argv[1:],"hvi:d:r:s:e:f:p:w:l:",["input=",
                        "header=","rrchannel=","sbpchannel=","ecgchannel=",
                        "ecgfilter=","pearsonr=","clusterwidth=","lag="])
except getopt.GetoptError:
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print("""analysis.py -i --input        <input file> 
          -d --header       <header length> 
          -r --hrchannel    <hr channel> 
          -s --sbpchannel   <sbp channel> 
          -e --ecgchannel   <ecg channel>
          -f --ecgfilter    <high pass filter>
          -p --pearsonr     <minimum run correlation>
          -w --clusterwidth <minimum run n>
          -l --lag          <hr offset from sbp>""")
        sys.exit()
    elif opt == '-v':
        _verbose = True
    elif opt in ("-i", "--input"):
        _fileName = arg
    elif opt in ("-d", "--header"):
        _headerLength = int(arg)
    elif opt in ("-r", "--hrchannel"):
        _RR = arg
    elif opt in ("-s", "--sbpchannel"):
        _SBP = arg
    elif opt in ("-e", "--ecgchannel"):
        _ECG = arg
    elif opt in ("-f", "--ecgfilter"):
        _filter = float(arg)
    elif opt in ("-p", "--pearsonr"):
        _pearson = float(arg)
    elif opt in ("-w", "--clusterwidth"):
        _width = int(arg)
    elif opt in ("-l", "--lag"):
        _lag = int(arg)
    
def pearsonR(x, y):
    # Assume len(x) == len(y)
    n = len(x)
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))
    psum = sum(map(lambda x, y: x * y, x, y))
    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * 
            (sum_y_sq - pow(sum_y, 2) / n), 0.5)
    if den == 0: return 0
    return num / den

# channel 5 is NIBP, 40 is SBP, 14 is ECG and 42 is HR.

#readCSVFile: str, int, int, int --> tuple-of-list-of-num
def readCSVFile(fileName, headerLength = 1, RRChannel = "CH42", 
            SBPChannel = "CH40", ECGChannel = "CH14", ECGFilter = 1.5):
    f = open(fileName,"r")
    print("Opening \"%s\""%fileName)
    reader = csv.reader(f,delimiter="\t")

    RR = [0]
    SBP = [0]

    RRIndex = 0
    SBPIndex = 0
    ECGIndex = 0

    movingAverageList = []

    ECGGrabLine = -1
    grabNewLine = True

    for lineNum, line in enumerate(reader):
        if(RRChannel in line):
            RRIndex = line.index(RRChannel)
            SBPIndex = line.index(SBPChannel)
            ECGIndex = line.index(ECGChannel)
            if(_verbose):
                print("ECGIndex: "+str(ECGIndex))
                print("RRIndex: "+str(RRIndex))
                print("SBPIndex: "+str(SBPIndex))

        if(lineNum > headerLength): #skip the headers
            if(float(line[ECGIndex]) >= ECGFilter): #filter out anything lower than the spike height
                if(grabNewLine):      #make sure we haven't already grabbed a number
                    
                    ECGGrabLine = lineNum + round(float(line[RRIndex]))*2
                    grabNewLine = False

            if(lineNum == ECGGrabLine):
                if(_verbose):
                    print("Grabbed @ "+str(lineNum))
                RR.append(float(line[RRIndex])/60)
                SBP.append(float(line[SBPIndex]))
                grabNewLine = True

    f.close()
    if(_verbose):
        print("SBP: "+str(RR))
    print("Closing \"%s\""%fileName)
    return (SBP, RR)  

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

    runDirection is either +1 or -1.

    clusterWidth is the minimum length of each run. lag is the difference
    in offset between the second list and the first list."""
    if(_verbose):
        print("RR length: "+str(len(RR)))
        print("SBP length: "+str(len(SBP)))

    runs = []
    currentRun = {"RR":[],"SBP":[]}
    direction = 0
   
    iterList = enumerate(zip(SBP, RR))
    next(iterList)

    for i, (SBPEntry, RREntry) in iterList:
        SBPDiff = SBPEntry - SBP[i - 1]
        RRDiff = RR[i + lag] - RR[i + lag - 1]
        if(RRDiff > 0 and SBPDiff > 0):
            if(direction >= 0):
                currentRun["RR"].append(RREntry)
                currentRun["SBP"].append(SBPEntry)
                direction = 1
            else:
                runs.append(currentRun)
                currentRun = {"RR":[RREntry], "SBP":[SBPEntry]}
                direction = 1
        elif(RRDiff < 0 and SBPDiff < 0):
            if(direction <= 0):
                currentRun["RR"].append(RREntry)
                currentRun["SBP"].append(SBPEntry)
                direction = -1
            else:
                #print(currentRun)
                runs.append(currentRun)
                currentRun = {"RR":[RREntry], "SBP":[SBPEntry]}
                direction = -1
        elif(RRDiff == 0 and SBPDiff == 0):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        else:
            #print(currentRun)
            runs.append(currentRun)
            currentRun = {"RR":[RREntry], "SBP":[SBPEntry]}
    return [r for r in runs if len(r["SBP"]) >= clusterWidth]

#findCorrelatedRuns: list-of-dict-of-list-of-num, num --> list-of-dict-of-list-of-num
def findCorrelatedRuns(runs, minCorrelation = 0.75):
    return [run for run in runs if pearsonR(run["SBP"], run["RR"]) > minCorrelation]

path = os.path.dirname(os.path.realpath(__file__)) + "/%s/"%_fileName
for file in os.listdir(path):
    currentFile = os.path.join(path, file)

    data = readCSVFile(currentFile, _headerLength, _RR, _SBP, _ECG, _filter)
    runs = findMatchingRuns(data[0],data[1], _width, _lag)
    correlatedRuns = findCorrelatedRuns(runs, _pearson)

    output = ["SBP, RR"]
    stuff = list(zip(data[0], data[1]))
    sbpAverage = 0
    rrAverage = 0
    for sbp, rr in stuff:
        output.append(",".join([str(sbp),str(rr)]))
        sbpAverage += sbp
        rrAverage  += rr
    count = len(stuff)
    output.append(",".join([str(sbpAverage/count), str(rrAverage/count)]))
    f = open(currentFile[:-4]+"_raw.csv","w")
    f.write("\n".join(output))
    f.close()

    if(_verbose):
        print("Correlated runs: "+str(correlatedRuns))

    f = open(currentFile[:-4]+"_correlated.csv","w")
    output = ["SBP, RR, SPRS"]
    for run in runs:
        zippedPairs = zip(run["SBP"], run["RR"])
        for sbp, rr in zippedPairs:
            output.append(",".join((str(sbp),str(rr), str(rr/sbp))))
    f.write("\n".join(output))
    f.close()
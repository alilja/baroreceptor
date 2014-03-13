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
_debug = False

try:
    opts, args = getopt.getopt(sys.argv[1:],"hvi:d:r:s:e:f:p:w:l:",["input=",
                        "header=","rrchannel=","sbpchannel=","ecgchannel=",
                        "ecgfilter=","pearsonr=","clusterwidth=","lag="])
except getopt.GetoptError:
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print("""analysis.py -i --input        <input file> 
          -h --header       <header length> 
          -r --hrchannel    <hr channel> 
          -s --sbpchannel   <sbp channel> 
          -e --ecgchannel   <ecg channel>
          -f --ecgfilter    <high pass filter>
          -p --pearsonr     <minimum run correlation>
          -w --clusterwidth <minimum run n>
          -l --lag          <hr offset from sbp>
          -d --divider      <divider character>
          --debug""")
        sys.exit()
    elif opt == '-v':
        _verbose = True
    elif opt in ("-i", "--input"):
        _fileName = arg
    elif opt in ("-h", "--header"):
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
    elif opt in ("-d", "--divider"):
        _lag = str(arg)
    elif opt in "--debug":
        _debug = True

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

# channel 5 is NIBP, 40 is SBP, 14 is ECG and 42 is HR.

#readCSVFile: str, int, int, int --> tuple-of-list-of-num
def readCSVFile(fileName, headerLength = 1, RRChannel = "CH42", 
            SBPChannel = "CH40", ECGChannel = "CH14", ECGFilter = 1.5):
    f = open(fileName,"r")
    print("Opening \"%s\""%fileName)
    reader = csv.reader(f,delimiter=",")

    RR = [0]
    SBP = [0]

    RRIndex = 0
    SBPIndex = 0
    ECGIndex = 0

    ECGGrabLine = -1
    grabNewLine = True

    header_end = 0
    header_read = False

    for lineNum, line in enumerate(reader):
        #if(lineNum < 50): print(line)
        if(RRChannel in line):
            if(_verbose): print("Grabbing header information @ line %s" % lineNum)
            header_end = lineNum + 1
            header_read = True
            RRIndex = line.index(RRChannel)
            SBPIndex = line.index(SBPChannel)
            ECGIndex = line.index(ECGChannel)
            if(_verbose):
                print("ECGIndex: "+str(ECGIndex))
                print("RRIndex: "+str(RRIndex))
                print("SBPIndex: "+str(SBPIndex))

        if(lineNum > header_end and header_read): #skip the headers
            if(float(line[ECGIndex]) >= ECGFilter): #filter out anything lower than the spike height
                #print("Filter exceeded")
                if(grabNewLine):      #make sure we haven't already grabbed a number
                    ECGGrabLine = lineNum + round(float(line[RRIndex]))*2
                    grabNewLine = False

            if(lineNum == ECGGrabLine):
                if(_debug): print("Grabbed @ "+str(lineNum))
                RR.append(60000/float(line[RRIndex])) #60000/HR = RR
                SBP.append(float(line[SBPIndex]))
                grabNewLine = True

    f.close()
    if(_verbose):p rint("SBP: "+str(RR))
    print("Closing \"%s\""%fileName)
    return (SBP, RR)  

#findMatchingRuns: list-of-num, list-of-num, num, num --> list-of-list-of-num
def findMatchingRuns(SBP, RR, clusterWidth = 3, lag = 0):
    if(_verbose):
        print("RR length: "+str(len(RR)))
        print("SBP length: "+str(len(SBP)))

    runs = []
    currentRun = {"RR":[],"SBP":[]}
    
    processing_list = zip(SBP, RR)
    list_length = len(processing_list)

    if(_verbose):
        print("Beginning positive checks.")

    iterList = enumerate(processing_list)
    next(iterList)

    for i, (SBPEntry, RREntry) in iterList:
        RREntry = RR[i + lag]

        if(lag > 0 and i + lag + 1 == list_length):
            print("Lag exceeding length of list; not processing final %s items" % lag)
            break

        prevRR  = RR[i - 1]
        prevSBP = SBP[i - 1]

        if(prevRR < RREntry and prevSBP < SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR == RREntry and prevSBP == SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR < RREntry and prevSBP == SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR == RREntry and prevSBP < SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)        
        else:
            runs.append(currentRun)
            currentRun = {"RR":[RREntry],"SBP":[SBPEntry]}

    if(_verbose):
        print("Beginning negative checks.")

    iterList = enumerate(processing_list)
    next(iterList)

    for i, (SBPEntry, RREntry) in iterList:
        RREntry = RR[i + lag]

        if(lag > 0 and i + lag + 1 == list_length):
            print("Lag exceeding length of list; not processing final %s items" % lag)
            break

        prevRR  = RR[i - 1]
        prevSBP = SBP[i - 1]

        if(prevRR > RREntry and prevSBP > SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR == RREntry and prevSBP == SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR > RREntry and prevSBP == SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)
        elif(prevRR == RREntry and prevSBP > SBPEntry):
            currentRun["RR"].append(RREntry)
            currentRun["SBP"].append(SBPEntry)        
        else:
            runs.append(currentRun)
            currentRun = {"RR":[RREntry],"SBP":[SBPEntry]}

    return [r for r in runs if len(r["RR"]) >= clusterWidth]

#findCorrelatedRuns: list-of-dict-of-list-of-num, num --> list-of-dict-of-list-of-num
def findCorrelatedRuns(runs, minCorrelation = 0.75):
    output = []
    for run in runs:
        # check to see if all of one column is identical
        sb_entries = run["SBP"]
        rr_entries = run["RR"]

        if(sum(sb_entries) - len(sb_entries)*sb_entries[0] == 0 or sum(rr_entries) - len(rr_entries)*rr_entries[0] == 0):
            if(_verbose):
                print("Skipping entry starting with SBP = {0}, RR = {1} because one column sums to zero.".format(sb_entries[0],rr_entries[0]))
        else:
            r = pearsonR(run["SBP"], run["RR"])
            if(r >= minCorrelation or r <= -minCorrelation):
                output.append(run)
    return output

print("Reading input file...")
data = readCSVFile(_fileName, _headerLength, _RR, _SBP, _ECG, _filter)

output = ["SBP, RR"]
stuff = zip(data[0], data[1])
for a, b in stuff:
    output.append(",".join([str(a),str(b)]))
f = open(_fileName[:-4]+"_raw.csv","w")
f.write("\n".join(output))
f.close()

print("Looking for matching runs...")
matchingRuns = findMatchingRuns(data[0],data[1], _width, _lag)
print("> Found %s matching, non-correlated runs." % len(matchingRuns))
print("Calculating correlation...")
correlatedRuns = findCorrelatedRuns(matchingRuns, _pearson)

print("> Found %s correlated runs." % len(correlatedRuns))
if(_verbose):
    print("Correlated runs: %s" % correlatedRuns)

output_file = _fileName[:-4]+"_correlated.csv"

f = open(output_file,"w")
output = ["SBP, RR"]
for run in matchingRuns:
    zippedPairs = zip(run["SBP"], run["RR"])
    for pair in zippedPairs:
        output.append(str(pair[0]) + ", " + str(pair[1]))
    output.append(",\n")
f.write("\n".join(output))
f.close()
print("Output written to %s." % output_file)
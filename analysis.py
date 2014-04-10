#!/usr/bin/python
import csv
import sys
import getopt
from collections import namedtuple

_verbose = False
_fileName = ""
_headerLength = 33
_RR = ""
_SBP = ""
_ECG = ""
_filter = 1.5
_pearson = 0.85
_width = 3
_lag = 0

_debug = False


## TODO improve docstrings (by which I mean add them)
## TODO move command parser to its own function
## TODO remove doubling in find_matches and make it run with positive/negative numbers
## TODO write some unit tests

BiometricData = namedtuple("BiometricData","hr_index nibp_index ecg_index data")

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
          -g --debug""")
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

def find_header(data, flag):
    for i, line in enumerate(data):
        for item in line:
            if(flag in item):
                return i
    else:
        return -1

# channel 5 is NIBP, 40 is SBP, 14 is ECG, and 42 is HR.
## TODO Remove hr_channel if we definitely don't need it ##
def process_csv(file, hr_channel = "CH42", nibp_channel = "CH5", ecg_channel = "CH14"):
    try:
        f = open(file, "r")
    except IOError:
        print "Could not read \"%s\""%file

    csv_f = csv.reader(f,delimiter="\t")
    search_list = [line for line in csv_f]

    header_location = find_header(search_list, "CH")
    header = search_list[header_location]
    search_list = search_list[header_location+2:-1]
    if(_verbose): print("Found header @ %s"%header_location)

    for i, line in enumerate(search_list):
        # Convert items in list to floats
        for j, item in enumerate(line):
            try:
                if(item.strip()):
                    line[j] = float(item) 
                else:
                    # strip out whitespace at the end
                    line = line[0:j]
            except ValueError as detail:
                print("ValueError on post-header line {0}, {1}: {2}. Interpretor: {3}".format(i, j, item, detail))
                raise
        search_list[i] = line
    f.close()
    return BiometricData(hr_index   = header.index(hr_channel),
                         nibp_index = header.index(nibp_channel),
                         ecg_index  = header.index(ecg_channel),
                         data = search_list)

# expects the BiometricData namedtuple
def find_spikes(bio_data, threshold = 1.5):
    SBP = [-1] # start with something in the list so we can just append new RR intervals
    RR  = []

    spike_a_location = 0
    look_for_b_spike = False
    nibp_search = []

    # basically the idea is to find the first spike in the ECG (just look for it to exceed a threshold)
    # then find the next spike
    # use the indexes of those two as a window to look in the NIBP to find the max
    # that max is the SBP for that interval
    # then use the index of that max to find the HR for the PREVIOUS SBP
    for line_num, line in enumerate(bio_data.data):
        if(line[bio_data.ecg_index] >= threshold):
            if(look_for_b_spike == False):
                if(_debug): print("Found spike A @ %s" % line_num)

                spike_a_location = line_num
                look_for_b_spike = True
                nibp_search.append(line[bio_data.nibp_index])

            elif(look_for_b_spike == True and line_num > spike_a_location + 50):
                if(_debug): print("Found spike B @ %s \n" % line_num)

                sbp = max(nibp_search)
                rr = line_num - spike_a_location

                SBP.append(sbp)
                RR.append(rr)

                nibp_search = []
                look_for_b_spike = False

        if(look_for_b_spike): nibp_search.append(line[bio_data.nibp_index])
    if(SBP == [-1] or not RR):
        raise RuntimeError("SBP or RR appear empty; try changing the ECG threshold.")
    return (SBP, RR)    

def find_matches(SBP, RR, clusterWidth = 3, lag = 0):
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

#correlate_runs: list-of-dict-of-list-of-num, num --> list-of-dict-of-list-of-num
def correlate_runs(runs, minCorrelation = 0.75):
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

def main():
    # Read input, find spikes
    print("Reading input file...")
    bio_data = process_csv(_fileName)
    if(_debug): print("Processed CSV. bio_data: %s" % bio_data)
    data = find_spikes(bio_data, _filter)

    # Produce _raw file
    output = ["SBP, RR"]
    stuff = zip(data[0], data[1])
    for a, b in stuff:
        output.append(",".join([str(a),str(b)]))
    f = open(_fileName[:-4]+"_raw.csv","w")
    f.write("\n".join(output))
    f.close()

    # Match
    print("Looking for matching runs...")
    matchingRuns = find_matches(data[0],data[1], _width, _lag)
    print("Found %s matching, non-correlated runs." % len(matchingRuns))

    # Correlate
    print("Calculating correlation...")
    correlatedRuns = correlate_runs(matchingRuns, _pearson)
    print("Found %s correlated runs." % len(correlatedRuns))
    if(_verbose):
        print("Correlated runs: %s" % correlatedRuns)

    # Write correlated file
    f = open(_fileName[:-4]+"_correlated.csv","w")
    output = ["SBP, RR"]
    for run in matchingRuns:
        zippedPairs = zip(run["SBP"], run["RR"])
        for pair in zippedPairs:
            output.append(str(pair[0]) + ", " + str(pair[1]))
        output.append(",\n")
    f.write("\n".join(output))
    f.close()
    print("Output written to %s." % _fileName[:-4]+"_correlated.csv")

if __name__ == "__main__":
    main()
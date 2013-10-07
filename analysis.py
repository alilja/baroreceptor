import csv

#readCSVFile: str --> tuple-of-list-of-num
def readCSVFile(fileName):
    f = open(fileName,"r")
    reader = csv.reader(f)

    HR = [0]
    systolic = [0]

    lineNum = 0
    for line in reader:
        if(lineNum > 1): #skip the headers
            if(float(line[5]) != HR[-1]):
                HR.append(float(line[5]))
            if(float(line[4]) != systolic[-1]):
                systolic.append(float(line[4]))

        lineNum += 1
    f.close()

    return (systolic, HR)  

#findMatchingRuns: list-of-num, list-of-num, num, num --> list-of-list-of-num
def findMatchingRuns(systolic, HR, clusterWidth = 3, lag = 0):
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

    runs = []
    count = 0
    direction = 0
    for i in range(1, len(systolic)):
        sysDiff = systolic[i] - systolic[i - 1]
        HRDiff = HR[i + lag] - HR[i + lag - 1]

        if(HRDiff > 0 and sysDiff > 0):
            if(direction >= 0):
                count += 1
                direction = 1
            else:
                runs.append([i - count - 1, i - 1, count, direction])
                count = 0
                direction = 1
                count = 1
        elif(HRDiff < 0 and sysDiff < 0):
            if(direction <= 0):
                count += 1
                direction = -1
            else:
                runs.append([i - count - 1, i - 1, count, direction])
                count = 0
                direction = -1
                count = 1
        elif(HRDiff == 0 and sysDiff == 0):
            count +=1 
        else:
            runs.append([i - count - 1, i - 1, count, direction])
            count = 0



    return [r for r in runs if r[2] >= clusterWidth]

data = readCSVFile("testData.csv")
runs = findMatchingRuns(data[0],data[1], 2, 1)

print(runs)

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
        print(out)

    print("\n\n")
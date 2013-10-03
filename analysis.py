import csv

def get_monotonic_subsequences(data, min_length=3):
    direction = data[1] - data[0]  # determine direction of initial subsequence
    subsequences = []
    cur_seq = [data[0]]
    for i in range(1, len(data)):
        if direction > 0:
            if data[i] >= data[i - 1]:
                cur_seq.append((i,data[i]))
            else:
                subsequences.append(cur_seq)
                cur_seq = [(i, data[i])]
                direction = -1
        else:
            if data[i] <= data[i - 1]:
                cur_seq.append((i,data[i]))
            else:
                subsequences.append(cur_seq)
                cur_seq = [(i, data[i])]
                direction = -1

    subsequences.append(cur_seq)
    return [x for x in subsequences if len(x) >= min_length]

def baroAnalysis(fileName, clusterWidth = 3, lag = 1):
    f = open(fileName,"r")
    reader = csv.reader(f)

    HR = [0]
    systolic = [0]

    lineNum = 0
    for line in reader:
        if(lineNum > 1): #skip the headers
            if(line[5] != str(HR[-1])):
                HR.append(float(line[5]))
            if(line[4] != str(systolic[-1])):
                systolic.append(float(line[4]))

        lineNum += 1
    f.close()  

    print(len(systolic))
    print(len(HR))

    count = 0

    for i in range(0, len(HR)-1):
        sysDiff = systolic[i+1] - systolic[i]
        HRDiff = HR[i+1] - HR[i]

        if(sysDiff * HRDiff > 0): #moving in the same direction
            print("same d")















baroAnalysis("testData.csv")
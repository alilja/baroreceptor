import csv

def baroAnalysis(fileName, clusterWidth = 3, lag = 1):
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

    print(len(systolic))
    print(len(HR))

    runs = []
    count = 0
    direction = 0
    for i in range(1, len(systolic)):
        sysDiff = systolic[i] - systolic[i - 1]
        HRDiff = HR[i] - HR[i - 1]

        if(HRDiff > 0 and sysDiff > 0):
            if(direction >= 0):
                count += 1
                direction = 1
            else:
                if(count >= 0):
                    runs.append([i - count - 1, i - 1, count, direction])
                    count = 0
                direction = 1
                count = 1
        elif(HRDiff < 0 and sysDiff < 0):
            if(direction <= 0):
                count += 1
                direction = -1
            else:
                if(count >= 0):
                    runs.append([i - count - 1, i - 1, count, direction])
                    count = 0
                direction = -1
                count = 1
        elif(HRDiff == 0 and sysDiff == 0):
            count +=1 
        else:
            if(count >= 0):
                runs.append([i - count - 1, i - 1, count, direction])
                count = 0
        print(str(i)+": "+str(direction))



    return [r for r in runs if r[2] >= clusterWidth]

print(baroAnalysis("testData.csv"))
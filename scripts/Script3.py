import pandas as pd

def getPos(subDf):
    list_ = []
    for index, row in subDf.iterrows():
        start = row['Start']
        end = row['End']
        list_.append([start, end])
    return list_
def checkDist(positions, minDist):
    passCount = 0
    for el in range(0, len(positions), 1):
        if el < 2:
            end1 = positions[el][1]
            start2 = positions[el + 1][0]

            if (start2 - end1) >= minDist:
                passCount += 1

        else:
            end1 = positions[el][1]
            start2 = positions[0][0]
            if abs((start2 - end1)) >= minDist:
                passCount += 1

    if passCount == 3:
        return True
    else:
        return False
def distanceChecker():
    genes = open('').read().splitlines()
    locationsDF = pd.read_csv('', sep='\t')
    minDist = 100000
    count = 0
    save = []
    for comb in genes:
        print(comb)
        splits = comb.split('_')
        subDf = locationsDF[locationsDF['Gene'].isin(splits)]
        positions = getPos(subDf)
        distCheck = checkDist(positions, minDist)
        print(count)
        save.append([comb, distCheck])
        count += 1

    with open('', 'w') as out:
        for i in save:
            for x in i:
                out.write(str(x) + '\t')
            out.write('\n')

if __name__ == "__main__":
    distanceChecker()




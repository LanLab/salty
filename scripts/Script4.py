import glob
import pandas as pd
import pandas as pd
import multiprocessing

def countRedundantGenesMain():

    #cluster sizes
    clusterSizesIn = open('').read().splitlines()
    clusterSizesDict = {}
    for line in clusterSizesIn[1:]:
        col = line.split(',')
        clusterSizesDict[col[0]]= col[1]

    #relational information
    combFolderIn = open('').read().splitlines()
    combFolderDict = {}
    for lines in combFolderIn:
        comb = lines.split('\t')[0].strip()
        folder = lines.split('\t')[1].strip()
        combFolderDict[comb] = folder

    manager = multiprocessing.Manager()
    summaryList = manager.list()

    combs = open('').read().splitlines()
    ins = [(comb, combFolderDict, clusterSizesDict, summaryList) for comb in combs]
    run_multiprocessing(countRedundantGenes, ins, 1)

    with open ('', 'w') as out:
        out.write(f'Combination\tNon-Redundant-Lineages\tIsolates-In-Non-Redundant-Lineages\tSingleton-Non-Redundant-Lineages\n')
        for element in summaryList:
            for i in element:
                out.write(str(i) + '\t')
            out.write('\n')

def countRedundantGenes(comb, combFolderDict, clusterSizesDict, summaryList):
    folder = combFolderDict[comb]
    input = open(
        f'XX{folder}/XX/{comb}XX').read().splitlines()
    nonRedunLin = 0
    isolateInNonRedunLin = 0
    for line in input[1:]:
        cols = line.split(',')
        alleleCount = 0
        for allele in cols[1:]:
            if allele != '':
                alleleCount += 1

        if alleleCount == 1:
            nonRedunLin += 1
            lineage = cols[0]
            linStrainSize = int(clusterSizesDict[lineage])
            isolateInNonRedunLin = isolateInNonRedunLin + linStrainSize
            percen = round((isolateInNonRedunLin/10000*100),2)

    summaryList.append([comb, nonRedunLin, isolateInNonRedunLin, percen])
    print(len(summaryList), comb, nonRedunLin, isolateInNonRedunLin, percen) #add out dict

def run_multiprocessing(func, i, n_processors):
    with multiprocessing.Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)

if __name__ == "__main__":
    countRedundantGenesMain()




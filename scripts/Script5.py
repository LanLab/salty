import pandas as pd
import glob
from matplotlib import pyplot as plt
import numpy as np

def resProfCleaner(resProf):
    out = []
    for i in resProf:
        if ';' in i:
            cols1 = i.split(';')
            for x in cols1:
                    out.append(x)
        else:
            out.append(i)
    print(out)
    out = list(set(sorted(out)))
    print(out)

    return out

def getAMTPattern():
    AMRclassesIn = ''
    abricateFolder = ''
    SALTyLineagesIn = ''
    SALTyLineages = pd.read_csv(SALTyLineagesIn, sep='\t')
    AMRclasses = open(AMRclassesIn).read().splitlines()
    subsetGenomes = SALTyLineages[(SALTyLineages['clonal complex'] == 'CC1') & (SALTyLineages['Lineage'] == '20')]['Genome'].tolist()

    prevComb = {}
    compList = []
    for filename in glob.iglob(f'{abricateFolder}/*_out.csv'):
        genome = filename.split('/')[-1].split('_')[0]
        if genome in subsetGenomes:
            subDF = pd.read_csv(filename)
            resProf = sorted(subDF['RESISTANCE'].tolist())
            cleanedResProf = resProfCleaner(resProf)
            combProf = '|'.join(cleanedResProf)
            if combProf not in prevComb.keys():
                prevComb[combProf] = 1
            else:
                prevComb[combProf] = prevComb[combProf] + 1

            sublist = []
            for i in AMRclasses:
                if i in cleanedResProf:
                    sublist.append(1)
                else:
                    sublist.append(0)
            compList.append(sublist)

    completeMatrix = pd.DataFrame(columns=AMRclasses, data=compList)
    completeMatrix.to_csv('')

    freqMatrix = pd.DataFrame(index=prevComb.keys(), columns=AMRclasses + ['Isolate Count', 'AMR Count'])
    print(f'Writing Out Combs For: {len(prevComb.keys())}')
    for key, values in prevComb.items():
        els = key.split('|')
        for el in els:
            freqMatrix.loc[key, el] = 'True'
        freqMatrix.loc[key, 'Isolate Count'] = values
        freqMatrix.loc[key, 'AMR Count'] = len(els)
    freqMatrix.to_csv('')
    subfreqMatrix = freqMatrix.sort_values('Isolate Count', ascending=False)
    subfreqMatrix = subfreqMatrix[subfreqMatrix['Isolate Count'] >= 10]
    subfreqMatrix = subfreqMatrix.dropna(axis=1, how='all')
    subfreqMatrix = subfreqMatrix.fillna('False')
    subfreqMatrix.to_csv('')

if __name__ == "__main__":
    getAMTPattern()

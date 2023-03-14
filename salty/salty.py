import argparse
import os, shutil
from os.path import join
import subprocess
import glob
import pandas as pd
from time import sleep as sl
import collections
import multiprocessing
import time
import textwrap
import sys

#caller
def caller(path, args, start_time_ongoing):
    start_time_analysis = time.time() ##FIX
    alleles = {'Lineage':'-','SACOL1908': '-', 'SACOL0451': '-', 'SACOL2725': '-'}
    parseInput(path, args) #alter DB path option
    alleles = filtCalledAlleles(alleles, args)
    alleles = getLineageFromAllele(alleles, args)
    alleles = checkFailedLineage(alleles, path, args)
    generateReport(alleles, args)
    timer(os.environ['accession'], start_time_ongoing, start_time_analysis)
def parseInput(path, args):
    if path[0] == 'assembly':
        os.environ['accession'] = getAccession(path[1])
        mkdirOutput(args)
        outPath = args.output_folder + '/' + os.environ['accession'] + '/kma_' + os.environ['accession']
        p = subprocess.Popen(['kma','-i', path[1], '-t_db', args.kma_index, '-o', outPath, '-t', '1'])
        (output, err) = p.communicate()
        p_status = p.wait()

    elif path[0] == 'pairedEndReadForward':
        os.environ['accession'] = getAccession(path[1])
        mkdirOutput(args)
        outPath = args.output_folder + '/' + os.environ['accession'] + '/kma_' + os.environ['accession']
        inputSTR = f"""kma -ipe {path[1].replace('_1.','_*.')} -t_db {args.kma_index} -o {outPath} -t 1"""
        p = subprocess.Popen(inputSTR, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
def filtCalledAlleles(alleles, args):
    resPath = args.output_folder + '/' + os.environ['accession'] + '/kma_' + os.environ['accession'] + '.res'
    df = pd.read_csv(resPath ,sep='\t')
    df = df[(df['Template_Identity'] == 100.00) & (df['Template_Coverage'] == 100.00) & (df['Query_Identity'] == 100.00)
        & (df['Query_Coverage'] == 100.00)]
    for index, row in df.iterrows():
        gene = row['#Template'].split(':')[0]
        allele = int(row['#Template'].split(':')[1])
        alleles[gene] = allele

        print('Passed: \t gene:' + str(gene) + '\t allele:' + str(allele))

    return alleles
def getLineageFromAllele(alleles, args):
    lineageAlleles = pd.read_csv(args.lineages)
    filtLineageAllelesDF = filtLineageAlleles(alleles, lineageAlleles)

    if filtLineageAllelesDF.shape[0] == 1:
        lineage = filtLineageAllelesDF['Lineage'].values[0]
        alleles['Lineage'] = lineage
        return alleles

    elif filtLineageAllelesDF.shape[0] == 0:
        alleles['Lineage'] = 'No lineages association.'
        return alleles

    elif filtLineageAllelesDF.shape[0] > 1:
        alleles['Lineage'] = 'Mulitple Lineage Found'
        filtLineageAllelesDF.to_csv(args.output_folder + "/" + os.environ['accession'] + "/" + os.environ['accession'] + '_multipleLineageAlleles.csv', index=False)
        return alleles

def checkFailedLineage(alleles, path, args):
    if args.mlstPrediction:
        if alleles['Lineage'] == 'No lineages association.':
            base = os.path.join(os.path.dirname(os.path.realpath(__file__)))
            MLSTin = base + '/resources/MLSTtoSaLTy.csv'
            MLSTtoSaLTyDf = pd.read_csv(MLSTin)
            MLSTtype = getMLSTtype(path, args)
            if MLSTtype.isdigit():
                MLSTtype = int(MLSTtype)
                associatedSaLTy = MLSTtoSaLTyDf[MLSTtoSaLTyDf['7-Gene-MLST'] == MLSTtype]['SaLTy-Lineage'].values[0]
                alleles['Lineage'] = f'*{associatedSaLTy}'
                return alleles
            else:
                os.environ['accession'] = getAccession(path[1])
                print(f""""{os.environ['accession']} is untypable with MLST. Therefore, lineage can not be predicted with SaLTy or MLST.""")
                alleles['Lineage'] = '*No lineages association.'
                return alleles
        else:
            return alleles
    else:
        return alleles

def getMLSTtype(path, args):
    os.environ['accession'] = getAccession(path[1])
    MLSTOutput = subprocess.check_output(['mlst', '--quiet','--scheme', 'saureus', '--threads', str(args.threads), path[1]])
    type = MLSTOutput.decode('utf-8').split('\t')[2]
    return type

def filtLineageAlleles(alleles, lineageAllelesDF):
    filtLineageAllelesDF = pd.DataFrame(columns=alleles.keys())
    for gene in alleles.keys():
        if gene != 'Lineage':
            filtAllele = lineageAllelesDF[lineageAllelesDF[gene] == alleles[gene]]
            if filtAllele.shape[0] > 0:
                filtAlleleIndex = filtAllele.index.values[0]
                if filtAlleleIndex not in list(filtLineageAllelesDF.index):
                    filtLineageAllelesDF = filtLineageAllelesDF.append(filtAllele, sort=True)
    return filtLineageAllelesDF

#aux functions
def getAccession(path):
    return path.split('/')[-1].split('.')[0].split('_')[0]
def generateReport(alleles, args):
    print(os.environ['accession'] + ': writing output.')
    print(alleles)
    outMeta = getOutMeta(args)
    with open(args.output_folder + "/" + os.environ['accession'] + "/" + os.environ['accession'] + "_lineage." + outMeta[0], 'w') as out:
        out.write("Genome")
        for header in alleles.keys():
            out.write(outMeta[1] + str(header))
        out.write('\n')
        out.write(str(os.environ['accession']))
        for allele in alleles.values():
            out.write(outMeta[1] + str(allele))
        out.write('\n')
def generateSummary(args):
    #TODO   change to share a managers list between multiple processing. avoid reading in the end.
    print('Generating Summary Report. Writing to ' + args.output_folder)
    saveList = []
    header=False
    for report in glob.iglob(args.output_folder + '/*/*_lineage.*'):
        infile = open(report).read().splitlines()
        if header == False:
            saveList.append(infile[0])
            header = True
        saveList.append(infile[-1])

    meta = getOutMeta(args)
    with open(args.output_folder + '/summaryReport.' + meta[0], 'w') as out:
        for line in saveList:
            out.write(line + '\n')
def getOutMeta(args):
    if args.csv_format:
        outMeta = ('csv', ',')
        return outMeta
    else:
        outMeta = ('txt', '\t')
        return outMeta
def timer(genome, start_time_ongoing, start_time_analysis):
    # Total time elapsed
    # since the timer started
    totaltime = round((time.time() - start_time_ongoing), 4)
    isolateTime = round((time.time() - start_time_analysis), 4)

    # # Printing the lap number,
    # # lap-time and total time
    # print('TESTINGLINE',totaltime, isolateTime, genome)

#I/O
def argsParser():
    parser = argparse.ArgumentParser(prog='PROG')

    general = parser.add_argument_group('GENERAL')
    general.add_argument('-t','--threads', default=1, type=int, help='Number of threads (speeds up parsing raw reads).')
    general.add_argument('-f','--force', default=False, action='store_true',  help='Overwite existing output folder.')
    general.add_argument('--report', default=False, action='store_true',  help='Only generate summary report from previous SALTy outputs.')
    general.add_argument('-v','--version', default=False, action='store_true')
    general.add_argument("--check", action='store_true', help="check dependencies are installed")

    inputs = parser.add_argument_group('INPUT')
    inputs.add_argument('-i','--input_folder', help='Folder of genomes (*.fasta or *.fna) and/or pair end reads (each accession must have *_1.fastq.qz and *_2.fastq.')

    output = parser.add_argument_group('OUTPUT')
    output.add_argument('-o','--output_folder', help='Output Folder to save result.', default='stdout')
    output.add_argument('-c','--csv_format', action='store_true', help='Output file in csv format.')
    output.add_argument('-s','--summary', action='store_true', help='Concatenate all output assignments into single file.')

    paths = parser.add_argument_group('DATABASE & PROGRAM Paths')
    base = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    paths.add_argument('-l','--lineages', default=base + '/resources/alleles/alleles.csv', help='Path to specific alleles for each lineage.')
    paths.add_argument('-k','--kma_index', default=base + '/resources/kmaIndex/kmaIndex', help='Path to indexed KMA database.')
    paths.add_argument('-m','--mlstPrediction', action='store_true', default=True, help='Explained in ReadMe. Used as backup when lineage is unable to be called through SaLTy screening.')

    return(parser.parse_args())
def check_deps(checkonly, args):
    depslist = ["kma"]
    f = 0
    for dep in depslist:
        rc = subprocess.call(['which', dep], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if rc == 0:
            if checkonly:
                sys.stderr.write(f'{dep:<10}:{"installed":<10}\n')
        else:
            sys.stderr.write(f'{dep:<10}:{" Missing in path, Please install ":<10}{dep}\n')
            f += 1
    if f > 0:
        sys.stderr.write("One or more dependencies missing.'\n")
        sys.stderr.write("If KMA use: 'conda install -c bioconda kma'\n")
        sys.exit(1)
    else:
        if checkonly:
            sys.stderr.write("All dependencies are present.\n")
            sys.exit(0)
        else:
            return
def mkdirOutput(args):

    if os.path.exists(args.output_folder + '/' + os.environ['accession']):
        if args.force:
            shutil.rmtree(args.output_folder + '/' + os.environ['accession'])
            os.mkdir(args.output_folder + '/' + os.environ['accession'])
        else:
            print('Error: Directory Exists: ' + args.output_folder + '/' + os.environ['accession'])
            print('Manul kill required (ctrl + C on Linux)')
            sys.exit(0)
    else:
        os.mkdir(args.output_folder + '/' + os.environ['accession'])
def checkInputReads(reads_folder):
    files = list(glob.iglob(reads_folder + '/*'))
    count = 0
    for file in files:
        if 'fastq.gz' in file:
            count += 1

    if count == 2:
        return True
    else:
        return False
def collectGenomes(args):
    paths = []
    for fasta in glob.iglob(args.input_folder + '/*fasta'):
        paths.append(['assembly',fasta])
    for fna in glob.iglob(args.input_folder + '/*fna'):
        paths.append(['assembly',fna])
    for forwardRead in glob.iglob(args.input_folder + '/*_1.fastq.gz'):
        accesion = forwardRead.split('/')[-1].replace('_1.fastq.gz','')
        reverseRead = f"{args.input_folder}/{accesion}_2.fastq.gz"
        if os.path.isfile(reverseRead):
            paths.append(['pairedEndReadForward', forwardRead])
        else:
            print(f"Failed to find paired end reads for {accesion}")
    return paths

def run_multiprocessing(func, i, n_processors):
    with multiprocessing.Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)

#main
def main():

    args = argsParser()
    version = "1.0.4"
    start_time_ongoing = time.time()

    if args.check:
        check_deps(True, args)
        sys.exit(0)

    elif args.version:
        print(f"salty version: {version}")
        sys.exit(0)

    elif args.report:
        generateSummary(args)

    elif args.input_folder:
        inputPaths = collectGenomes(args)
        if (len(inputPaths) > 0):
            inps = [(path, args, start_time_ongoing) for path in inputPaths]
            run_multiprocessing(caller, inps, args.threads)
            generateSummary(args)

    else:
        print("Input folder doees not have assemblies and/or paired end reads.")
        exit()

if __name__ == '__main__':
    main()

import argparse
import os, shutil
from os.path import join
import subprocess
import glob
import pandas as pd
from time import sleep as sl
import collections
import multiprocessing

#testing
import time

def timer(genome, start_time_ongoing, start_time_analysis):
    # Total time elapsed
    # since the timer started
    totaltime = round((time.time() - start_time_ongoing), 4)
    isolateTime = round((time.time() - start_time_analysis), 4)

    # Printing the lap number,
    # lap-time and total time
    print('TESTINGLINE',totaltime, isolateTime, genome)

#caller
def caller(genome, args, start_time_ongoing):
    start_time_analysis = time.time() ##FIX

    args.input_genome = genome
    alleles = {'Lineage':'-','SACOL1908': '-', 'SACOL0451': '-', 'SACOL2725': '-'}
    parseInput(args) #alter DB path option
    alleles = filtCalledAlleles(alleles, args)
    alleles = getLineageFromAllele(alleles, args)
    generateReport(alleles, args)
    timer(genome, start_time_ongoing, start_time_analysis)
def parseInput(args):
    if args.input_genome:
        os.environ['accession'] = getAccession(args.input_genome)
        mkdirOutput(args)
        outPath = args.output_folder + '/' + os.environ['accession'] + '/kma_' + os.environ['accession']
        p = subprocess.Popen(['kma','-i', args.input_genome, '-t_db', args.kma_index, '-o', outPath, '-t', '1'])
        (output, err) = p.communicate()
        p_status = p.wait()

    elif args.reads_folder:
        os.environ['accession'] = getAccession(args.reads_folder.split(',')[0])
        mkdirOutput(args)
        outPath = args.output_folder + '/' + os.environ['accession'] + '/kma_' + os.environ['accession']
        inputSTR = f"""kma -ipe {args.reads_folder}/*.fastq.gz -t_db {args.kma_index} -o {outPath} -t 1"""
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
    if args.format_csv:
        outMeta = ('csv', ',')
        return outMeta
    else:
        outMeta = ('txt', '\t')
        return outMeta

#I/O
def argsParser():
    parser = argparse.ArgumentParser(prog='PROG')

    general = parser.add_argument_group('GENERAL')
    general.add_argument('-t','--threads', default=1, type=int, help='Number of threads (speeds up parsing raw reads).')
    general.add_argument('-f','--force', default=False, action='store_true',  help='Overwite existing output folder.')
    general.add_argument('--report', default=False, action='store_true',  help='Only generate summary report from previous SALTy outputs.')


    inputs = parser.add_argument_group('INPUT')
    # input = inputs.add_mutually_exclusive_group(required=True)
    input = inputs.add_mutually_exclusive_group()

    input.add_argument('-i','--genome_folder', help='Input folder with assembled DNA sequence file.')
    input.add_argument('-r','--reads_folder', help='Folder with forward and reverse raw reads (fastq.gz)')

    output = parser.add_argument_group('OUTPUT')
    output.add_argument('-o','--output_folder', help='Output Folder to save result.', default='stdout')
    output.add_argument('-csv','--format_csv', action='store_true', help='Output file in csv format.')
    output.add_argument('-s','--summary', action='store_true', help='Concatenate all output assignments into single file.')

    paths = parser.add_argument_group('DB PATHS')
    base = os.path.dirname(os.path.realpath(__file__))
    paths.add_argument('-l','--lineages', default=base + '/alleles/alleles.csv', help='Path to specific alleles for each lineage.')
    paths.add_argument('-k','--kma_index', default=base + '/kmaIndex/kmaIndex', help='Path to indexed KMA database.')

    # gitHub = parser.add_argument_group('PUBLIC ACCESS')
    # gitHub.add_argument('FIX ADD LANLAB LINK')

    return(parser.parse_args())
def mkdirOutput(args):

    if os.path.exists(args.output_folder + '/' + os.environ['accession']):
        if args.force:
            shutil.rmtree(args.output_folder + '/' + os.environ['accession'])
            os.mkdir(args.output_folder + '/' + os.environ['accession'])
        else:
            print('Error: Directory Exists: ' + args.output_folder + '/' + os.environ['accession'])
            exit()
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
def checkInputDatabases(args):
    kmaInput = list(set([x.split(':')[0] for x in open(args.kma_index + '.name').read().splitlines()]))
    lineageInput = open(args.lineages).read().splitlines()[0].split(',')[1:]

    if set(kmaInput) == set(lineageInput):
        return True
    else:
        return False
def collectGenomes(args):
    outList = []
    for type in ('*.fasta','*.fna'):
        outList.extend(glob.iglob(join(args.genome_folder + '/', type)))
    return outList
def collectRawReads(args):
    outList = []
    for forwardRead in glob.iglob(args.reads_folder + '/*_1.fastq*'):
        # genome = filename.split('/')[-1].split('_')[0]
        # if os.path.isfile(args.reads_folder + '/' + genome + '_2.fastq*'):
        reverseRead = forwardRead.replace('_1.','_2')
        paths = forwardRead + ',' + reverseRead
        outList.append(paths)
    return outList

def run_multiprocessing(func, i, n_processors):
    with multiprocessing.Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)

#main
def SALTy():
    import time ##FIX
    args = argsParser()
    start_time_ongoing = time.time()

    if not checkInputDatabases(args):
        print('Input KMA Index and Predefined Lineages do not contain same genes. Rebuild either.')
        exit()

    if args.report:
        generateSummary(args)

    else:
        if args.genome_folder:
            inputGenomePath = collectGenomes(args)
            if len(inputGenomePath) > 0:
                inps = [(x, args, start_time_ongoing) for x in inputGenomePath]
                run_multiprocessing(caller, inps, args.threads)
                generateSummary(args)

            else:
                print("Input folder doees not have assemblies with extension '.fna' or '.fasta'")
                exit()

        elif args.reads_folder: #only configured for paired ends raw reads
            inputRawReads = collectRawReads(args)
            if len(inputRawReads) > 0:
                inps = [(x, args) for x in inputRawReads]
                run_multiprocessing(caller, inps, args.threads)
                generateSummary(args)
            else:
                print('Input folder does not contain pair end raw reads.')
                exit()

if __name__ == '__main__':
    SALTy()

##TODO add in dependency checker
#Install Dependencies
#python3
#kma
#pandas
#np
#mlst
#sqlite
#glob


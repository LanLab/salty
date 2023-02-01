import pandas as pd
from itertools import combinations
import multiprocessing
from sys import argv
import os
from sklearn import metrics
import numpy as np
import sqlite3
import contextlib
from filelock import FileLock
import ctypes
from time import sleep as sl
import callSpecificAlleles

#creating databases
def createConnection(db_name):
    """ create a database connection to a SQLite database """
    conn = None
    conn = sqlite3.connect(f'{os.environ["outputfolder"]}/{db_name}.db')
    # conn = sqlite3.connect(':memory:')
    return conn
def createDatabases():

    createConnection('report_sqlite') #create report database

    createConnection('summary_sqlite') #create summary database
def createTables():
    #summary table
    sql_create_summary_table = """CREATE TABLE IF NOT EXISTS summary (
								  gene TEXT,
                                    accuracy REAL,
                                    weightedAvg REAL,
                                    macroAvgF1 REAL,
                                    macroAvgPrecision REAL,
                                    macroAvgRecall REAL,
                                    failedStrains INTEGER
                                    );"""
    conn = createConnection('summary_sqlite')
    with conn:
        c = conn.cursor()
        c.execute(sql_create_summary_table)

    #reports table
    sql_create_report_table = """CREATE TABLE IF NOT EXISTS reports (
                                    gene TEXT,
                                    cluster INTEGER,
                                    precision REAL,
                                    recall REAL,
                                    f1score REAL,
                                    support REAL 
                                    );"""
    reportconn = createConnection('report_sqlite')
    with reportconn:
        c = reportconn.cursor()
        c.execute(sql_create_report_table)
def getMacroStat(metric, statsDict, roundNum):

    tmpList = []
    for k,v in statsDict.items():
        if k.isnumeric() and int(k) != 100000:
            tmpList.append(v[metric])
    return round(np.average(tmpList), roundNum)

#checking DB out
def counterPrinter(gene, counterStuff):
    with counterStuff[1]:
        print(gene, counterStuff[0].value)
        counterStuff[0].value = counterStuff[0].value + 1

#writing to databases
def updateSummaryList(gene, reportDict, summaryList, failedStrains):
    accuracy = round(reportDict['accuracy'], 3)
    weightedAvg = round(reportDict['weighted avg']['f1-score'], 3)
    macroAvgF1 = round(reportDict['macro avg']['f1-score'], 3)
    macroAvgPrecision = round(getMacroStat('precision', reportDict, 3), 2)
    macroAvgRecall = round(getMacroStat('recall', reportDict, 3), 2)
    summaryList.append([gene, accuracy, weightedAvg, macroAvgF1, macroAvgPrecision, macroAvgRecall, failedStrains])
def updateSummaryDB(summaryList):
    with contextlib.closing(createConnection('summary_sqlite')) as conn:  # auto-closes
        with conn:  # auto-commits
            with contextlib.closing(conn.cursor()) as cursor:  # auto-closes
                ql_summary_headers = "INSERT INTO summary VALUES(?,?,?,?,?,?,?);"
                cursor.executemany(ql_summary_headers, summaryList)
    # print(f'Outputing {len(summaryList)} to summary db.')
    # print('\t'.join([str(x) for x in summaryList]))
    sl(1)

def updateReportList(gene, reportDict, reportList):
    for k, v in reportDict.items():
        tmpList = []
        if k.isnumeric() and int(k) != 100000:
            tmpList.append(gene)
            tmpList.append(k)

            for x, y in v.items():
                tmpList.append(round(y, 2))

        if len(tmpList) == 6:
            reportList.append(tmpList)
def updateReportDB(reportList):
    with contextlib.closing(createConnection('report_sqlite')) as connRep:  # auto-closes
        with connRep:  # auto-commits
            with contextlib.closing(connRep.cursor()) as cursorRep:  # auto-closes
                sql_reports_to_table = f"""INSERT INTO reports (gene, cluster, precision, recall, f1score, support) VALUES (?,?,?,?,?,?);"""
                cursorRep.executemany(sql_reports_to_table, reportList)
    sl(1)
    # with open('/srv/scratch/lanlab/liam/7_Saureus_lineages/2021-11-03-testingRandom5000v3/4_clusterSpecificAssemblies/testing/testreports.txt', 'w') as out:
    #     for i in reportList:
    #         for x in i:
    #             out.write(str(x) + '\t')
    #         out.write('\n')

    # print(f'Outputing {len(reportList)} to report db.')

def run_multiprocessing(func, i, n_processors):
    with multiprocessing.Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)
def main():
    alleles = pd.read_csv(argv[1], sep='\t', low_memory=False)
    # combs = [[x] for x in list(alleles.columns) if 'SACOL' in x]
    os.environ['outputfolder'] = argv[2] + '/' + argv[3] + '/'
    combs = [(x.split(',')) for x in open(argv[4], encoding='utf-8-sig').read().splitlines()]
    workers=int(argv[5])

	#create databases
    createDatabases()
    createTables()

    #create output variables
    manager = multiprocessing.Manager()
    summaryList = manager.list()
    reportList = manager.list()
    count = manager.Value(ctypes.c_ulonglong, 0)
    counter_lock = manager.Lock()

    # multiprocess
    inps = [(x, alleles, summaryList,reportList, (count, counter_lock)) for x in combs]
    run_multiprocessing(callSpecificAlleles.run, inps, workers)

    updateSummaryDB(summaryList)
    updateReportDB(reportList)

if __name__ == '__main__':
    main()


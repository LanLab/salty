from time import sleep as sl
import argparse
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
import pandas as pd
import time
import itertools
from datetime import datetime
from csv import reader
import sys
import os
from collections import OrderedDict
from collections import Counter

def unneg(a):
    if "-" in a:
        return int(a.split("_")[0][1:])
    else:
        return int(a)

def import_data(args):
    inprofiles = open(args.alleleprofiles, "r").read().splitlines()

    profs = {}
    st_to_strain = {}
    strain_to_st = {}

    for line in inprofiles[1:]:
        col = line.split("\t")
        z=0
        if col[0] not in profs and col[1] != "":
            if args.grapetree_format:
                noneg = [unneg(x) for x in col[2:]]
            else:
                noneg = [unneg(x) for x in col[3:]]

            for i in noneg:
                if i == 0:
                    z+=1
            zperc = (float(z)/float(len(noneg)))*100
            if zperc <= args.zeromaxperc and z <= args.maxzerocount:
                profs[col[0]] = noneg
                st_to_strain[col[1]] = [str(col[0])]
        else:
            st_to_strain[col[1]].append(str(col[0]))
        # print(col[1],z)
        strain_to_st[col[0]] = col[1]

    idlist = list(profs.keys())
    inprofs = [profs[x] for x in idlist]
    # print(st_to_strain)
    return idlist, profs, st_to_strain, strain_to_st

def mgt_dist_metric(a,b,args):
    match = 0.0
    missmatch = 0.0
    for i in range(len(a)):
        aAllele = a[i]
        bAllele = b[i]
        if aAllele == 0 or bAllele == 0:
            missmatch += args.zeroweight
        elif aAllele == bAllele:
            match += 1.0
        else:
            missmatch += 1.0
            if missmatch >= float(args.max_missmatch):
                return missmatch
    return missmatch

def run_dist(args,profs,idlist):
    # get existing distances if present, also get list of strain ids
    if args.distances:
        prin
        initiald = pd.read_csv(args.distances,sep="\t",index_col=0)
        header = list(initiald.head())
    else:
        header = []

    # make empty dataframe with all current strains
    # for size in range(100,1000,100):
    # print(len(idlist))

    start_time = time.time()
    newdf = pd.DataFrame(index=idlist,columns=idlist)
    pairs = itertools.combinations(idlist,r=2)
    tot = (len(idlist)*len(idlist)-len(idlist))/2
    frac = int(tot*0.01)
    c=0

    newids = len(list(set(idlist).difference(set(header))))
    if idlist == list(header):
        print("\nAll isolates found in input pairwise distances file")
        return initiald
    print("\nCalculating pairwise distances\nDone\t\texists\t\tnew\t\t% done\t\ttime")
    pre = 0
    new = 0
    for pair in pairs:
        id1 = pair[0]
        id2 = pair[1]
        if id1 in header and id2 in header:
            dist = int(initiald.at[int(id1),str(id2)])
            pre +=1

        else:
            ap1 = profs[id1]
            ap2 = profs[id2]

            dist = mgt_dist_metric(ap1,ap2,args)
            new+=1

        newdf.at[id1,id2] = dist
        newdf.at[id2,id1] = dist
        c+=1

        if c%frac == 0:
            print("{}\t\t{}\t\t{}\t\t{}%\t\t{:.3f} seconds".format(c,pre,new,int((c/tot)*100),time.time() - start_time), end="\r", flush=True)
            # pre = 0
            # new = 0
    print("{}\t\t{}\t\t{}\t\t{}%\t\t{:.3f} seconds".format(c, pre, new, 100, time.time() - start_time),
          end="\r", flush=True)
    for i in idlist:
        newdf.at[i,i] = 0

    print("\n")
    return newdf

def get_distances_frm_args(args):
    diststring = args.dist_limits
    dists = diststring.split(",")
    distances = []
    for i in dists:
        if "-" in i:
            n = i.split("-")
            nlist = list(range(int(n[0])+1,int(n[1])+2))
            # distance cutoffs seems to be non inclusive i.e. cutoff of 3 means max distance is 2
            # therefore need to add 1 to all values
        else:
            nlist = [int(i)+1]
        distances += nlist
    return distances

def run_agglom_cluster(args,idlist,distancedf,id_to_strain,strain_to_st):

    distances = get_distances_frm_args(args)

    clusterlists = {}

    preference = []
    # for id in idlist:
    #     print(id_to_strain)
    #     preference.append(len(id_to_strain[id]))
    print(len(idlist))
    start_time = time.time()
    print("\n\nCalculating cluster distance")

    max_all_merged = 0
    calced_distances = []
    for dist in distances:
        print("{}\t\t{:.3f} seconds".format(dist-1,time.time() - start_time), end="\r", flush=True)
        clusters = AgglomerativeClustering(n_clusters=None,distance_threshold=dist,affinity="precomputed", linkage="single").fit_predict(distancedf)
        # print(dir(clusters))
        print(dist, len(clusters))

        # labels = clusters.labels
        #chindexscore = metrics.calinski_harabasz_score(distancedf, clusters)
        noclusters = len(set(clusters))

        #print("ch index for cutoff {}: {}, per cluster: {}".format(dist-1,chindexscore,float(chindexscore)/float(noclusters)))
        clusterls = list(clusters)
        clusterlists[dist] = clusterls

        unique_clusters = len(set(clusterls))
        if unique_clusters <= 2:
            max_all_merged += 1
            calced_distances.append(dist)
            if max_all_merged == 10:
                print(f'All isolates clustered at level: {dist}')
                break
        else:
            calced_distances.append(dist)

        strain_to_cluster = {}

    distances = calced_distances
    for i in range(len(idlist)):
        id = idlist[i]
        strain=str(id)
        # for strain in id_to_strain[id]:
        if strain not in strain_to_cluster:
            strain_to_cluster[strain] = {}
        for d in distances:
            dreal = d-1
            strain_to_cluster[strain][dreal] = int(clusterlists[d][i])
        strain_to_cluster[strain][0] = strain_to_st[strain]
    print("\nclustering time", (" --- %s seconds ---" % (time.time() - start_time)))

    # for strain in list(sorted(strain_to_cluster.keys())):
    #     for distance in list(sorted(strain_to_cluster[strain].keys())):
    #         print(strain,distance,strain_to_cluster[strain][distance])

    return strain_to_cluster, distances

def writeout(args,strain_to_cluster, distances):

    outf = open(args.clusters_out,"w")
    outf.write("Isolate")
    for d in distances:
        outf.write("\tODC_"+str(d-1))
    outf.write("\n")

    for strain in strain_to_cluster:
        outf.write(strain)
        for distance in distances:
            dreal = distance-1
            clust = strain_to_cluster[strain][dreal]
            outf.write("\t"+str(clust))
        outf.write("\n")
    outf.close()

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-a", "--alleleprofiles", default='/Users/liamcheneyy/Desktop/heirCC/aps8.txt',
                        help="file containing allele profiles (output from get_ap9_for_strainls.py)")
    parser.add_argument("-l", "--dist_limits",
                        help="comma separated list of cluster cutoffs or range or both i.e 1,2,5 or 1-8 or 1,2,5-10",
                        default="0-2000")
    parser.add_argument("--clusters_out", default="/Users/liamcheneyy/Desktop/heirCC/species_clusters.txt")
    parser.add_argument("--outputPrefix", help="output path and prefix for output file generation", default="/Users/liamcheneyy/Desktop/heirCC/")
    parser.add_argument("--grapetree_format",
                        help="allele profiles are in grapetree format (no dST column and neg STs converted)",action = 'store_true')
    parser.add_argument("-m", "--max_missmatch",
                        help="maximum number of missmatches reported between 2 isolates", default=2001, type=int)
    parser.add_argument("-d", "--distances",
                        help="file containing distances corresponding to the alleleprofiles file (from previous run of this script if applicable)")
    parser.add_argument("-z", "--zeroweight",
                         help="fraction of an intact missmatch to assign to missing alleles", default=0.0, type=float)
    parser.add_argument("--zeromaxperc",
                        help="percentage of allele profile allowed to be 0", default=50, type=float)
    parser.add_argument("--maxzerocount",
                        help="number of 0 allele calls allowed in allele profile", default=50, type=float)

    args = parser.parse_args()


    return args

def main():
    args = parseargs()

    idlist, profs, st_to_strain, strain_to_st = import_data(args)

    distance_df = run_dist(args,profs,idlist)
    distance_df.to_csv(args.outputPrefix + '_distances.txt', sep='\t')
    strain_to_cluster,distances = run_agglom_cluster(args, idlist, distance_df, st_to_strain, strain_to_st)

    writeout(args,strain_to_cluster, distances)

if __name__ == '__main__':
    main()

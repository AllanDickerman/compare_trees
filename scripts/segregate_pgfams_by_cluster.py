import sys
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract PGFams to files based on cluster (column of tsv)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv", metavar='tsv file', type=str, help="tsv file with pgfams in first column, cluster in Nth col")
    parser.add_argument("--column", metavar='number', type=int, help="column to use (1-based)")
    parser.add_argument("--prefix", type=str, default="pgfam_cluster_", help="output to files: pgfam_cluster_#.txt")
    parser.add_argument("--max_cluster", type=int, default=10, help="max number of clusters to output")

    args = parser.parse_args()

    print(f"tsv file = {args.tsv}")
    print(f"column = {args.column}")
    print(f"prefix = {args.prefix}")

    pgfam_cluster = {}
    cluster_pgfam = {}
    with open(args.tsv) as F:
        for line in F:
            fields = line.split('\t')
            m = re.match(".*(PGF_\d+)", fields[0])
            if m:
                pgfam = m.group(1)
                cluster = int(fields[args.column-1])
                if cluster < 0:
                    continue # skip -1, which means unclassified
                pgfam_cluster[pgfam] = cluster
                #print(f"{m.group(1)}\t{cluster}")
                if cluster not in cluster_pgfam:
                    cluster_pgfam[cluster] = set()
                cluster_pgfam[cluster].add(pgfam)

    for cluster in sorted(cluster_pgfam)[:args.max_cluster]:
        cluster_file = args.prefix+"_"+str(cluster)+".txt"
        F = open(cluster_file, 'w')
        F.write("\n".join(cluster_pgfam[cluster])+"\n")

    
if __name__ == '__main__':
    #print("Here")
    exit(main())

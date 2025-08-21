
### Analyze all aligned pgfams individually
### Given n alignments for the same genomes,
### build trees from each pgfam separately
### using FastTree


import sys
import re
import os
import os.path
import glob
import argparse
import subprocess
import tempfile
import copy
from time import time, localtime, strftime                        
from Bio import SeqIO
from Bio import AlignIO

def genomeIdFromFigId(figId):
    m = re.match("fig\|(\d+\.\d+)", figId)
    if m:
        return m.group(1)
    return None

def readAlignment(file, format="fasta"):
    alignment = AlignIO.read(file, format)
    alignment.sort()
    return alignment

def concatenate_alignments(alignments):
# alignments is dictionary of Bio.multipleSeqAlignments
    id_set = set(seq.id for aln in alignments for seq in alignments[aln])
    totalLength = 0
    taxonIdList=sorted(id_set)
    alignmentIdList = sorted(alignments)
    seqDict = {}
    for seqid in id_set:
        seqDict[seqid] = ""
    for alignmentId in alignmentIdList:
        alignment = alignments[alignmentId]
        length = alignment.get_alignment_length()
        temp_id_set = set(id_set)
        for record in alignment:
            seqDict[record.id] += str(record.seq)
            temp_id_set.remove(record.id)
        for missing_id in temp_id_set:
            seqDict[missing_id] += '-' * length
    return seqDict

def build_gene_tree(alignment, alphabet='dna', program=None, protein_model=None, threads = 2):
    print(f"build_gene_tree() type = {program}")
    # Create a temporary directory
    original_dir = os.getcwd()
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
    tree_string = ''
    tree_likelihood = -1e6
    time = 1.0
    if program == 'raxml':
        alignment_path = os.path.abspath(alignment)
        with tempfile.TemporaryDirectory(prefix="assess_codon_tree") as tmpdirname:
            print('created temporary directory', tmpdirname)
            os.chdir(tmpdirname)
            # -f d new rapid hill-climbing
            # DEFAULT: ON This is the default RAxML tree search algorithm and is substantially faster than the original
            raxml_command = ['raxml', '-s', alignment_path, '-n', 'tree', '-f', 'd', '-p', '1234', '-T', str(threads)]
            if alphabet == 'dna':
                raxml_command.extend(['-m', 'GTRGAMMA'])
            elif program == 'protein':
                raxml_command.extend(['-m', 'PROTGAMMA'+protein_model])
            proc = subprocess.run(raxml_command)
            if os.path.exists("RAxML_result.tree"):
                tree_string = open("RAxML_result.tree").read()
            if os.path.exists("RAxML_info.tree"):
                F = open("RAxML_info.tree")
                for line in F:
                    m = re.match("Final .* of best tree (\S+)", line)
                    if m:
                        tree_likelihood = float(m.group(1))
                    m = re.match("Overall execution time: (\S+)", line)
                    if m:
                        time = float(m.group(1))
                F.close()
        os.chdir(original_dir)
    if program == 'fasttree':
        nt_flag = ''
        if alphabet == 'dna':
            nt_flag = '-nt'
        result = subprocess.run(f"FastTree {nt_flag} {alignment}", shell=True, capture_output=True, text=True)
        tree_string = result.stdout
        #print(f"tree = {tree_string}")
        #print(f"stats = {result.stderr}")
        m = re.search("LogLk = (\S+)", result.stderr)
        if m:
            tree_likelihood = m.group(1)
        m = re.search("Total time: (\S+)", result.stderr)
        if m:
            time = float(m.group(1))

    return [tree_string, tree_likelihood, time]

def main():
    parser = argparse.ArgumentParser(description="Build one tree per PGFam alignment", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--alphabet", metavar="dna/protein", type=str, help="dna or protein")
    parser.add_argument("--program", metavar="fasttree/raxml", type=str, help="raxml or fasttree")
    parser.add_argument("--proteinModel", metavar="model", type=str, help="protein substitution model")
    parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for tree building")
    parser.add_argument("--max_trees", metavar="max", type=int, default=0, help="stop after this many trees")

    args = parser.parse_args()

    pgfam_tree_file = open(f"pgfam_{args.alphabet}_trees.nex", 'w') 
    pgfam_tree_file.write(f"#NEXUS\n")
    #ptfam_tree_file.write("\nBegin Taxa;\n  Dimensions ntax={len(genomeIds)};\n")
    #pgfam_tree_file.write("  TaxLabels "+" ".join(sorted(genomeIds))+";\nEND;\n\n")
    pgfam_tree_file.write("Begin TREES;\n")

    al_score_file_name = f"pgfam_{args.alphabet}_tree_stats_file.txt"
    alignment_score_file = open(al_score_file_name, "w")
    alignment_score_file.write("pgfam\tnumpos\tlogL\tperSiteLL\ttime\n")
    alignment_files = glob.glob(f"{args.alphabet}_alignments/PGF*")
    i = 0
    for alignment_file in alignment_files:
        i += 1
        if args.max_trees and i >= args.max_trees:
            print(f" hit max_trees {args.max_trees}")
            break
        m = re.search("/(P.F_\d+)", alignment_file)
        pgfam = m.group(1)
        print(f" build tree for {pgfam}")

        (tree, likelihood, time) = build_gene_tree(alignment_file, program='fasttree', alphabet=args.alphabet, protein_model = args.proteinModel, threads = args.threads)

        pgfam_tree_file.write(f"Tree {pgfam}_{args.alphabet}_tree = {tree}\n")
        F = open(alignment_file)
        F.readline() # first def line
        aligned_seq1 = F.readline().rstrip()
        alignment_length = len(aligned_seq1)
        per_site_loglikelihood = float(likelihood) / alignment_length
        alignment_score_file.write(f"{pgfam}\t{alignment_length}\t{likelihood}\t{per_site_loglikelihood}\t{time}\n")
        alignment_score_file.flush()
        F.close()

    pgfam_tree_file.write("END;\n")

if __name__ == '__main__':
    exit(main())


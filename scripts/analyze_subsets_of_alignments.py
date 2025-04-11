
### Analyze subsets of alignments
### Given n alignments for the same genomes,
### for various values of m < n
### build trees from concatenations of m alignments
### using RAxML


import sys
import re
import os
import os.path
import glob
import argparse
import subprocess
import tempfile
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

def build_dna_tree(dna_alignments, threads = 2):
    print(f"build_dna_tree() num alignments passed = {len(dna_alignments)}")
    print(f"  ids={dna_alignments.keys()}")
    concatenated_alignment = concatenate_alignments(dna_alignments)
    # Create a temporary directory
    original_dir = os.getcwd()
    return_value = None
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
    with tempfile.TemporaryDirectory(prefix="assess_codon_tree") as tmpdirname:
        print('created temporary directory', tmpdirname)
        os.chdir(tmpdirname)
        F = open("dna_alignment.fa", 'w')
        for seq_id in sorted(concatenated_alignment):
            F.write(f">{seq_id}\n{concatenated_alignment[seq_id]}\n")
        F.close()
        # -f d new rapid hill-climbing
        # DEFAULT: ON This is the default RAxML tree search algorithm and is substantially faster than the original
        raxml_command = ['raxml', '-s', 'dna_alignment.fa', '-n', 'dna_tree', '-m', 'GTRGAMMA', '-f', 'd', '-p', '1234', '-T', str(threads)]
        proc = subprocess.run(raxml_command)
        if os.path.exists("RAxML_result.dna_tree"):
            return_value = open("RAxML_result.dna_tree").read()
        os.chdir(original_dir)
    return return_value

def build_protein_tree(protein_alignments, protein_model, threads=2):
    print(f"build_protein_tree() num alignments passed = {len(protein_alignments)}")
    print(f"  ids={protein_alignments.keys()}")
    concatenated_alignment = concatenate_alignments(protein_alignments)
    # Create a temporary directory
    original_dir = os.getcwd()
    return_value = None
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
    with tempfile.TemporaryDirectory(prefix="assess_codon_tree") as tmpdirname:
        print('created temporary directory', tmpdirname)
        os.chdir(tmpdirname)
        F = open("protein_alignment.fa", 'w')
        for seq_id in sorted(concatenated_alignment):
            F.write(f">{seq_id}\n{concatenated_alignment[seq_id]}\n")
        F.close()
        # -f d new rapid hill-climbing
        # DEFAULT: ON This is the default RAxML tree search algorithm and is substantially faster than the original
        raxml_command = ['raxml', '-s', 'protein_alignment.fa', '-n', 'protein_tree', '-m', protein_model, '-f', 'd', '-p', '1234', '-T', str(threads)]
        proc = subprocess.run(raxml_command)
        if os.path.exists("RAxML_result.protein_tree"):
            return_value = open("RAxML_result.protein_tree").read()
        os.chdir(original_dir)
    return return_value

def build_joint_tree(dna_alignments, protein_alignments, protein_model, threads=2):
    print(f"build_joint_tree() num alignments passed = {len(protein_alignments)}")
    print(f"  ids={protein_alignments.keys()}")
    concatenated_dna_alignment = concatenate_alignments(dna_alignments)
    concatenated_protein_alignment = concatenate_alignments(protein_alignments)
    # Create a temporary directory
    original_dir = os.getcwd()
    return_value = None
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
    with tempfile.TemporaryDirectory(prefix="assess_codon_tree") as tmpdirname:
        print('created temporary directory', tmpdirname)
        os.chdir(tmpdirname)
        total_dna_length = 0
        F = open("joint_alignment.fa", 'w')
        first = True
        for seq_id in sorted(concatenated_dna_alignment):
            F.write(f">{seq_id}\n{concatenated_dna_alignment[seq_id]}\n")
            if first:
                total_dna_length = len(concatenated_dna_alignment[seq_id])
                total_protein_length = len(concatenated_protein_alignment[seq_id])
                first = False
            F.write(concatenated_protein_alignment[seq_id]+"\n")
        F.close()
        if (total_protein_length * 3) != total_dna_length:
            raise Exception(f"Hey! total protein length = {total_protein_length} (x3={total_protein_length*3}, total dna length = {total_dna_length}")
        
        with open("partitions.txt", 'w') as PART:
            for i in range(1,4):
                PART.write(f"DNA, codon{i} = {i}-{total_dna_length}\\3\n")
            PART.write(f"{protein_model}, proteins = {total_dna_length+1}-{total_dna_length + total_protein_length}\n")

        # -f d new rapid hill-climbing
        # DEFAULT: ON This is the default RAxML tree search algorithm and is substantially faster than the original
        raxml_command = ['raxml', "-s", "joint_alignment.fa",   "-n", "joint_tree",   "-m",  "GTRGAMMA",   "-q",  "partitions.txt",  "-p", "12345", "-T", str(threads)]
        proc = subprocess.run(raxml_command)
        if os.path.exists("RAxML_result.joint_tree"):
            return_value = open("RAxML_result.joint_tree").read()
        os.chdir(original_dir)
    return return_value

def main():
    parser = argparse.ArgumentParser(description="Analye subsets of alignments, building trees", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dna", action='store_true', help="analyze DNA alignments")
    parser.add_argument("--protein", action='store_true', help="analyze DNA alignments")
    parser.add_argument("--joint", action='store_true', help="analyze merged DNA+protein alignments")
    parser.add_argument("--subsampleSize", "-s", metavar="#", type=int, help="number of genes per subsample")
    #parser.add_argument("--rateModel", metavar="model", type=str, choices = ['CAT', 'GAMMA'], default="GAMMA", help="variable rate category model CAT|GAMMA")
    parser.add_argument("--proteinModel", metavar="model", type=str, help="raxml protein substitution model")
    parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for raxml")

    args = parser.parse_args()
    starttime = time()
    genomeIds = set() # list of genome IDs for tree building (enforcing maxGenomesMissing and maxAllowedDupes)


    protein_alignment = {}
    for alignment_file in glob.glob("protein_alignments/*afa"):
        m = re.search("/(P.F_\d+)", alignment_file)
        if m:
            homology_group = m.group(1)
            al = readAlignment(alignment_file)
            protein_alignment[homology_group] = al
    print(f"num protein alignments is {len(protein_alignment)}")
    protein_list = sorted(protein_alignment.keys())

    dna_alignment = {}
    sample_id = None
    for alignment_file in glob.glob("dna_alignments/*afn"):
        m = re.search("/(P.F_\d+)", alignment_file)
        if m:
            homology_group = m.group(1)
            al = readAlignment(alignment_file)
            dna_alignment[homology_group] = al
    print(f"num dna alignments is {len(dna_alignment)}")
    dna_list = sorted(dna_alignment.keys())

    seq_id_set = set()
    for seq_record in dna_alignment[dna_list[0]]:
        seq_id_set.add(seq_record.id)
    print(f"num unique seq ids is {len(seq_id_set)}")

    for al in protein_alignment:
        num_found = 0
        for seq_record in protein_alignment[al]:
            if seq_record.id in seq_id_set:
                num_found += 1
        if num_found != len(seq_id_set):
            print(f"Hey! alignment {al} had only {num_found} seqs vs {len(seq_id_set)}")
        num_found = 0
        for seq_record in dna_alignment[al]:
            if seq_record.id in seq_id_set:
                num_found += 1
        if num_found != len(seq_id_set):
            print(f"Hey! alignment {al} had only {num_found} seqs vs {len(seq_id_set)}")

    print("ok")

    sample_size = args.subsampleSize
    if args.dna:
        dna_tree_file = open(f"dna_trees_sampleSize{sample_size}.nwk", 'w') 
    if args.protein:
        protein_tree_file = open(f"protein_trees_sampleSize{sample_size}.nwk", 'w') 
    if args.joint:
        joint_tree_file = open(f"joint_trees_sampleSize{sample_size}.nwk", 'w') 
    (dna_tree, protein_tree, joint_tree) = (None, None, None)
    print(f" generate {int(sample_size/len(dna_list))} trees")
    for sample_pos in range(sample_size-1, len(dna_list), sample_size):
        print(f" build tree {int(sample_pos/sample_size)} of {int(len(dna_list)/sample_size)}")
        dna_sample = {}
        protein_sample = {}
        for i in range(sample_pos, len(dna_list)):
            homology_group = dna_list[i]
            dna_sample[homology_group] = dna_alignment[homology_group]
            protein_sample[homology_group] = protein_alignment[homology_group]
        if args.dna:
            dna_tree = build_dna_tree(dna_sample, threads = args.threads)
            if dna_tree:
                dna_tree_file.write(dna_tree + "\n")
                dna_tree_file.flush()
        if args.protein:
            protein_tree = build_protein_tree(protein_sample, args.proteinModel, threads = args.threads)
            if protein_tree:
                protein_tree_file.write(protein_tree + "\n")
                protein_tree_file.flush()
        if args.joint:
            joint_tree = build_joint_tree(dna_sample, protein_sample, args.proteinModel, threads = args.threads)
            if joint_tree:
                joint_tree_file.write(joint_tree + "\n")
                joint_tree_file.flush()

if __name__ == '__main__':
    exit(main())


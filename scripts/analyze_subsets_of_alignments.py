
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

def build_tree(dna_alignments=None, protein_alignments=None, analysis_type=None, protein_model=None, threads = 2):
    print(f"build_tree() type = {analysis_type}")
    #print(f"  ids={dna_alignments.keys()}")
    # Create a temporary directory
    original_dir = os.getcwd()
    return_value = None
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
    with tempfile.TemporaryDirectory(prefix="assess_codon_tree") as tmpdirname:
        print('created temporary directory', tmpdirname)
        os.chdir(tmpdirname)
        # -f d new rapid hill-climbing
        # DEFAULT: ON This is the default RAxML tree search algorithm and is substantially faster than the original
        raxml_command = ['raxml', '-s', 'alignment.fa', '-n', 'tree', '-f', 'd', '-p', '1234', '-T', str(threads)]
        if analysis_type == 'dna':
            raxml_command.extend(['-m', 'GTRGAMMA'])
            concatenated_alignment = concatenate_alignments(dna_alignments)
        elif analysis_type == 'protein':
            raxml_command.extend(['-m', 'PROTGAMMA'+protein_model])
            concatenated_alignment = concatenate_alignments(protein_alignments)
        elif analysis_type == 'joint':
            raxml_command.extend(['-m', 'GTRGAMMA', '-q', 'partitions.txt'])
            concat_dna = concatenate_alignments(dna_alignments)
            concat_protein = concatenate_alignments(protein_alignments)
            concatenated_alignment = {}
            for seq_id in concat_dna:
                concatenated_alignment[seq_id] = concat_dna[seq_id]+concat_protein[seq_id]
            total_dna_length = len(concat_dna[seq_id])
            total_protein_length = len(concat_protein[seq_id])
            with open("partitions.txt", 'w') as PART:
                for i in range(1,4):
                    PART.write(f"DNA, codon{i} = {i}-{total_dna_length}\\3\n")
                PART.write(f"{protein_model}, proteins = {total_dna_length+1}-{total_dna_length + total_protein_length}\n")
        F = open("alignment.fa", 'w')
        for seq_id in sorted(concatenated_alignment):
            F.write(f">{seq_id}\n{concatenated_alignment[seq_id]}\n")
        F.close()
        proc = subprocess.run(raxml_command)
        (tree_string, tree_likelihood, raxml_time) = ['', 0, 0]
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
                    raxml_time = float(m.group(1))

        os.chdir(original_dir)
    return [tree_string, tree_likelihood, raxml_time]

def main():
    parser = argparse.ArgumentParser(description="Analye subsets of alignments, building trees", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dna", action='store_true', help="analyze DNA alignments")
    parser.add_argument("--protein", action='store_true', help="analyze DNA alignments")
    parser.add_argument("--joint", action='store_true', help="analyze merged DNA+protein alignments")
    parser.add_argument("--sampleSize", "-s", metavar="#", type=int, help="number of genes per sample")
    parser.add_argument("--proteinModel", metavar="model", type=str, help="raxml protein substitution model")
    parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for raxml")
    parser.add_argument("--alignment_score_file", type=str, help="file with scores generated by codon_tree.py")
    parser.add_argument("--max_trees", metavar="max", type=int, default=0, help="stop after this many trees")

    args = parser.parse_args()
    starttime = time()


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

    genomeIds = set()
    for seq_record in dna_alignment[dna_list[0]]:
        genomeIds.add(seq_record.id)
    print(f"num unique seq ids is {len(genomeIds)}")

    for al in protein_alignment:
        num_found = 0
        for seq_record in protein_alignment[al]:
            if seq_record.id in genomeIds:
                num_found += 1
        if num_found != len(genomeIds):
            print(f"Hey! alignment {al} had only {num_found} seqs vs {len(genomeIds)}")
        num_found = 0
        for seq_record in dna_alignment[al]:
            if seq_record.id in genomeIds:
                num_found += 1
        if num_found != len(genomeIds):
            print(f"Hey! alignment {al} had only {num_found} seqs vs {len(genomeIds)}")

    print("ok")

    alignment_score = {}
    alignment_score_headers = ['mean_squared_freq', 'num_pos', 'prop_gaps']
    pgfam_used = {}
    PGFam = ''
    with open(args.alignment_score_file) as F:
        header_list = F.readline().rstrip().split("\t") #PGFam   gaps    sum_squared_freq    mean_squared_freq   num_pos num_seqs    worstSeqScore   used_in_tree
        header_list.pop(0)
        for line in F:
            if not line.startswith('PGF'):
                continue
            fields = line.rstrip().split("\t")
            PGFam = fields.pop(0)
            alignment_score[PGFam] = {} 
            for index, header in enumerate(header_list):
                alignment_score[PGFam][header] = float(fields[index])
            alignment_score[PGFam]['prop_gaps'] = alignment_score[PGFam]['gaps']/(alignment_score[PGFam]['num_pos']*alignment_score[PGFam]['num_seqs'])
            pgfam_used[PGFam] = bool(fields[-1])
    print(f"Num alignment scores read = {len(alignment_score)}")

    sample_size = args.sampleSize
    if args.dna:
        dna_tree_file = open(f"dna_trees_sampleSize{sample_size}.nex", 'w') 
        dna_tree_file.write(f"#NEXUS\n\nBegin Taxa;\n  Dimensions ntax={len(genomeIds)};\n")
        dna_tree_file.write("  TaxLabels "+" ".join(sorted(genomeIds))+";\nEND;\n\nBegin TREES;\n")
    if args.protein:
        protein_tree_file = open(f"protein_trees_sampleSize{sample_size}.nex", 'w') 
        protein_tree_file.write(f"#NEXUS\n\nBegin Taxa;\n  Dimensions ntax={len(genomeIds)};\n")
        protein_tree_file.write("  TaxLabels "+" ".join(sorted(genomeIds))+";\nEND;\n\nBegin TREES;\n")
    if args.joint:
        joint_tree_file = open(f"joint_trees_sampleSize{sample_size}.nex", 'w') 
        joint_tree_file.write(f"#NEXUS\n\nBegin Taxa;\n  Dimensions ntax={len(genomeIds)};\n")
        joint_tree_file.write("  TaxLabels "+" ".join(sorted(genomeIds))+";\nEND;\n\nBegin TREES;\n")

    homologs_per_sample_file = open(f"homologs_sampleSize{sample_size}.txt", "w")

    al_score_file_name = f"alignment_sample_scores_sampleSize{sample_size}"
    if args.dna:
        al_score_file_name += "_dna"
    if args.protein:
        al_score_file_name += "_protein"
    if args.joint:
        al_score_file_name += "_joint"
    al_score_file_name += ".txt"
    alignment_sample_score_file = open(al_score_file_name, "w")
    (dna_tree, protein_tree, joint_tree) = (None, None, None)
    sample_scores = copy.deepcopy(alignment_score[PGFam])
    first_iteration = True
    print(f" generate {int(sample_size/len(dna_list))} trees")
    for sample_pos in range(0, len(dna_list)-sample_size+1, sample_size):
        if args.max_trees and (sample_pos/sample_size) >= args.max_trees:
            print(f" hit max_trees {args.max_trees}")
            break
        print(f" build tree {int(sample_pos/sample_size)} of {int(len(dna_list)/sample_size)}")
        dna_sample = {}
        protein_sample = {}
        homology_list = []
        for i in range(sample_size):
            homology_group = dna_list[sample_pos + i]
            dna_sample[homology_group] = dna_alignment[homology_group]
            protein_sample[homology_group] = protein_alignment[homology_group]
            homology_list.append(homology_group)
        sample_label = f"tree_{sample_size}genes_{sample_pos:04}"
        homologs_per_sample_file.write(sample_label+"\t"+"\t".join(homology_list)+"\n")

        all_pgfams_used = True
        for homolog in homology_list[1:]: # first element is already there
            for score in alignment_score_headers:
                sample_scores[score] += alignment_score[homolog][score]
            all_pgfams_used = all_pgfams_used and pgfam_used[homolog]
        alignment_score_list = []
        for score in alignment_score_headers:
            alignment_score_list.append(sample_scores[score] / sample_size)

        tree_score_headers = []
        tree_score_list = []

#def build_tree(dna_alignments=None, protein_alignments=None, analysis_type=None, protein_model=None, threads = 2):
        if args.dna:
            (dna_tree, dna_likelihood, dna_time) = build_tree(analysis_type='dna', dna_alignments = dna_sample, threads = args.threads)
            if dna_tree:
                dna_tree_file.write(f"Tree dna_{sample_label} = {dna_tree}\n")
                dna_tree_file.flush()
                tree_score_list.extend([dna_likelihood, dna_time])
                tree_score_headers.extend(['dna_likelihood', 'dna_time'])
        if args.protein:
            (protein_tree, protein_likelihood, protein_time) = build_tree(analysis_type='protein', protein_alignments = protein_sample, protein_model=args.proteinModel, threads = args.threads)
            if protein_tree:
                protein_tree_file.write(f"Tree protein_{sample_label} = {protein_tree}\n")
                protein_tree_file.flush()
                tree_score_list.extend([protein_likelihood, protein_time])
                tree_score_headers.extend(['protein_likelihood', 'protein_time'])
        if args.joint:
            (joint_tree, joint_likelihood, joint_time) = build_tree(analysis_type='joint', dna_alignments = dna_sample, protein_alignments = protein_sample, protein_model = args.proteinModel, threads = args.threads)
            if joint_tree:
                joint_tree_file.write(f"Tree joint_{sample_label} = {joint_tree}\n")
                joint_tree_file.flush()
                tree_score_list.extend([joint_likelihood, joint_time])
                tree_score_headers.extend(['joint_likelihood', 'joint_time'])

        if first_iteration:
            alignment_sample_score_file.write("sample\t"+ "\t".join(alignment_score_headers)+"\tall_used\t"+"\t".join(tree_score_headers)+"\n")
            first_iteration = False
        alignment_sample_score_file.write(sample_label)
        for s in alignment_score_list:
            alignment_sample_score_file.write(f"\t{s:.3}")
        alignment_sample_score_file.write(f"\t{all_pgfams_used}")
        for s in tree_score_list:
            alignment_sample_score_file.write(f"\t{s:.3}")
        alignment_sample_score_file.write("\n")

        for homolog in homology_list:
            if not pgfam_used[homolog]:
                alignment_sample_score_file.write(f"{homolog} not used!\n")
        alignment_sample_score_file.flush()

    if args.dna:
        dna_tree_file.write("END;")
    if args.protein:
        protein_tree_file.write("END;")
    if args.joint:
        joint_tree_file.write("END;")

if __name__ == '__main__':
    exit(main())


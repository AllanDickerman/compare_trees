import sys
import BitVector
import glob
import re

def read_fasta(fh):
    seqs={}
    seqid = None
    for line in fh:
        if line.startswith('>'):
            #print("line starts with "+line)
            seqid = line.rstrip()[1:].split(' ')[0]
            #print("Seq_id = "+seqid)
            seqs[seqid] = ''
        else:
            seqs[seqid] += line.rstrip()
    return seqs

def count_partitions(al, seqid_list):
    #seqid_list = sorted(al.keys())
    partition_counts = {}
    al_length = len(al[seqid_list[0]])
    print(f"count_partitions, al_length = {al_length}")
    for pos in range(al_length):
        bitvector_dict = {}
        for i, seqid in enumerate(seqid_list):
            state = al[seqid][pos]
            if state not in bitvector_dict:
                bitvector_dict[state] = BitVector.BitVector(size = len(seqid_list))
            bitvector_dict[state][i] = 1
        for state in bitvector_dict:
            bv = bitvector_dict[state]
            if not bv[0]:
                bv = ~bv
            bitstring = str(bv)
            if bitstring not in partition_counts:
                partition_counts[bitstring] = 0
            partition_counts[bitstring] += 1
    return(partition_counts, al_length)


def main():

    global_counts = {}
    alignment_counts = {}
    genome_list = []
    total_length = 0
    for alignment_file in glob.glob("dna_alignments/PGF*afn"):
        #if len(alignment_counts) > 400:
        #    break
        print(f"read {alignment_file}")
        m = re.match("dna_alignments/(PGF_\d+)", alignment_file)
        pgfam = m.group(1)
        alignment_counts[pgfam] = {}
        with open(alignment_file) as F:
            al = read_fasta(F)
            if not len(genome_list):
                for genome in sorted(al):
                    genome_list.append(genome)
            pc, alength = count_partitions(al, genome_list)
            total_length += alength
            for bs in sorted(pc, key=pc.get):
                if bs not in global_counts:
                    global_counts[bs] = 0
                global_counts[bs] += pc[bs]
                alignment_counts[pgfam][bs] = (pc[bs] / alength)

    for bs in global_counts:
        global_counts[bs] /= total_length # change to freqency

    
    F = open("alignment_partition_freq_matrix.txt", 'w')
    sorted_bs = sorted(global_counts, key=global_counts.get, reverse=True)
    alignment_score = {}
    # first row describes overall global pattern
    F.write('global')
    for bs in sorted_bs[:int(len(genome_list)*4)]:
        if global_counts[bs] <= 0.0001:
            break
        F.write(f"\t{global_counts[bs]:.4f}")
    F.write("\n")

    for pgfam in sorted(alignment_counts):
        score = 0
        F.write(pgfam)
        for bs in sorted_bs[:int(len(genome_list)*4)]:
            if global_counts[bs] <= 0.0001:
                break
            freq = 0
            if bs in alignment_counts[pgfam] and alignment_counts[pgfam][bs] > 0:
                freq = alignment_counts[pgfam][bs]
                score += (freq * global_counts[bs]) 
            F.write(f"\t{freq:.4f}")
        F.write("\n")
        alignment_score[pgfam] = score
    print("alignment partition freq matrix written to alignment_partition_freq_matrix.txt");
        

    F = open("alignments_scoreed_by_partition_freq_product.txt", 'w')
    for pgfam in sorted(alignment_counts):
        F.write(f"{pgfam}\t{alignment_score[pgfam]:.4f}\n")
    print("alignment scores written to alignments_scoreed_by_partition_freq_product.txt");
    return(0)

if __name__ == '__main__':
    exit(main())

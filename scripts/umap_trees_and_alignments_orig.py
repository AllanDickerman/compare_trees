import dendropy
import sys
import pandas
import glob
import re
import BitVector
import argparse

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

def count_alignment_partitions(al, seqid_list):
    #seqid_list = sorted(al.keys())
    partition_counts = {}
    al_length = len(al[seqid_list[0]])
    #print(f"count_partitions, al_length = {al_length}")
    for pos in range(al_length):
        bitvector_dict = {}
        for i, seqid in enumerate(seqid_list):
            state = al[seqid][pos]
            if state not in bitvector_dict:
                bitvector_dict[state] = BitVector.BitVector(size = len(seqid_list))
            bitvector_dict[state][i] = 1
        for state in bitvector_dict:
            bv = bitvector_dict[state]
            bitcount = bv.count_bits()
            if bitcount > len(seqid_list)/2:
                bv = ~bv
            bitstring = str(bv)
            if bitstring not in partition_counts:
                partition_counts[bitstring] = 0
            partition_counts[bitstring] += 1
    return(partition_counts, al_length)

def generate_umap_plot(data, n_neighbors=[10,30,100]):
    # yields 2-D embeddings of data vectors
    print("generate_umap_plot()")
    import umap.umap_ as umap
    umap_df = pandas.DataFrame(index = data.index)
    for n in n_neighbors:
        prefix = f"UMAP_N{n}_"
        reducer = umap.UMAP(n_neighbors = n)
        model = reducer.fit(data)
        u = model.transform(data)
        udf = pandas.DataFrame(u, columns=[prefix+"0", prefix+"1"], index=data.index)
        umap_df = pandas.concat([umap_df, udf], axis="columns")
        print(f"umap shape = {umap_df.shape}")
    return umap_df # 

def cluster_umap_plot(udf, min_samples=3):
    # groups gene families into clusters based on umap coordinates
    # run over range of eps and save in additional columns of the umap dataframe
    from sklearn.cluster import DBSCAN
    columns = udf.columns;
    umap_column_prefix = []
    for c in udf.columns:
        m = re.match("(UMAP_N\d+_)0", c)
        if m:
            umap_column_prefix.append(m.group(1))
    for prefix in umap_column_prefix:
        for eps in [0.1, 0.5, 1]:
            dbclust = DBSCAN(eps=eps, min_samples=min_samples).fit(udf[[prefix+'0', prefix+'1']])
            udf[prefix+"C"+str(eps)] = dbclust.labels_
    return


def main():
    parser = argparse.ArgumentParser(description="Analye subsets of alignments, building trees", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", metavar='files', type=str, nargs='*', help="nexus tree files")
    parser.add_argument("--alignments", metavar='dir', type=str, help="directory with alignments")
    parser.add_argument("--write_matrix", action='store_true', help="write partition matrix")

    args = parser.parse_args()

    print(f"tree files = {args.trees}")
    print(f"alignment dir = {args.alignments}")

    tree_bitstring_length = {}
    partition_count = {}
    global_tree_partition_length = {}
    taxa = dendropy.TaxonNamespace()
    tree_yielder = dendropy.Tree.yield_from_files(
            files=args.trees,
            schema='nexus',
            taxon_namespace=taxa,
            )
    for tree_idx, one_tree in enumerate(tree_yielder):
        print(one_tree.label+"    \t"+str(tree_idx))
        tree_label = one_tree.label.replace(' ', '_')
        m = re.search("(PGF.\d+)", one_tree.label)
        #if m:
        #    tree_label = m.group(1) #pgfam
        tree_bitstring_length[tree_label] = {}
        bem = one_tree.bipartition_edge_map
        for bipart in bem:
            if (bem[bipart].length is not None):
                bitstring = bipart.leafset_as_bitstring(reverse=True)
                tree_bitstring_length[tree_label][bitstring] = bem[bipart].length
                if not bitstring in partition_count:
                    partition_count[bitstring] = 0
                    global_tree_partition_length[bitstring] = 0
                partition_count[bitstring] += 1
                global_tree_partition_length[bitstring] += bem[bipart].length
    tree_bl_mean = 0
    for bs in global_tree_partition_length:
        global_tree_partition_length[bs] /= partition_count[bs]
        tree_bl_mean += global_tree_partition_length[bs]
    tree_bl_mean /= len(global_tree_partition_length)
    print(f"num tree partitions = {len(global_tree_partition_length)}, mean len={tree_bl_mean}")

    print(f"num trees = {len(tree_bitstring_length)}")
    taxon_list = []
    for i, taxon in enumerate(taxa):
        print(f"{i}\t{taxon.label}")
        taxon_list.append(taxon.label)
    min_bp_count = 0 #len(tree_partition_length)/2

    # now count partitions seen in alignments
    global_apart_counts = {}
    alignment_counts = {}
    alignment_length = {}
    total_length = 0
    alignment_type = '_alignment'
    m = re.search('([a-z]+_alignment)', args.alignments)
    if m:
        alignment_type = m.group(1)
    alignment_files = glob.glob(args.alignments+"/PGF*")
    for alignment_file in alignment_files:
        #print(f"read {alignment_file}")
        m = re.search("(PGF_\d+)", alignment_file)
        if m:
            pgfam = m.group(1)
        else:
            pgfam = alignment_file
        pgfam += alignment_type
        alignment_counts[pgfam] = {}
        with open(alignment_file) as F:
            al = read_fasta(F)
            pc, alength = count_alignment_partitions(al, taxon_list)
            total_length += alength
            alignment_length[pgfam] = alength
            #print(f" num partitions = {len(pc)}")
            for bs in sorted(pc, key=pc.get):
                if bs not in global_apart_counts:
                    global_apart_counts[bs] = 0
                global_apart_counts[bs] += pc[bs]
                alignment_counts[pgfam][bs] = (pc[bs] / alength)
    print(f"num bitstrings seen in alignments = {len(global_apart_counts)}")

    global_apart_freq = {}
    afreq_mean = 0
    for bs in global_apart_counts:
        global_apart_freq[bs] = global_apart_counts[bs] / total_length
        afreq_mean += global_apart_freq[bs]
    afreq_mean /= len(global_apart_freq)
    print(f"num global counts from alignments = {len(global_apart_counts)}, mean = {afreq_mean}")


    F = open("Alignment_vs_tree_partition_values.txt", 'w')
    i=0
    F.write("partition\talign\ttree\n")
    shared_bitstrings = []
    threshold_for_saving_alignment_partition = 0 #afreq_mean * 0.25
    threshold_for_saving_tree_partition = 0 #tree_bl_mean * 0.25
    for bs in sorted(global_apart_counts, key=global_apart_counts.get, reverse=True):
        i+=1
        tree_partition_length = 0
        if bs in global_tree_partition_length:
            tree_partition_length = global_tree_partition_length[bs]
            if (global_apart_freq[bs] > threshold_for_saving_alignment_partition) and (tree_partition_length > threshold_for_saving_tree_partition):
                shared_bitstrings.append(bs)
        F.write(f"{bs}\t{global_apart_freq[bs]}\t{tree_partition_length}\n")
    F.close()
    print(f"threshold for alignment parition = {threshold_for_saving_alignment_partition}")
    print(f"threshold for tree parition = {threshold_for_saving_tree_partition}")
    print(f"number of shared bitstrings = {len(shared_bitstrings)}")


    tree_vector_list = []
    tree_names = []
    for tree in tree_bitstring_length:
        tree_vector = []
        #print(tree)
        for bitstring in shared_bitstrings:
            bl = 0.0
            if bitstring in tree_bitstring_length[tree]:
                bl = tree_bitstring_length[tree][bitstring]
            tree_vector.append(bl)
        tree_vector_list.append(tree_vector)
        tree_names.append(tree.replace(' ', '_'))


    tree_v_df = pandas.DataFrame(tree_vector_list, index=tree_names, columns=shared_bitstrings)
    print(f"tree_v_df = {tree_v_df}")


    alignment_global_score = {}
    alignment_vector_list = []
    alignment_names = []
    average_alignment_score = 0
    max_alignment_score = 0
    for alignment in alignment_counts:
        alignment_names.append(alignment+"_dna_alignment")
        alignment_vector = []
        #print(alignment)
        score = 0
        for bitstring in shared_bitstrings:
            pf = 0.0
            if bitstring in alignment_counts[alignment]:
                pf = alignment_counts[alignment][bitstring]
                score += alignment_counts[alignment][bitstring] * global_apart_freq[bitstring]
            alignment_vector.append(pf)
        alignment_vector_list.append(alignment_vector)
        alignment_global_score[alignment] = score
        average_alignment_score += score
        if score > max_alignment_score:
            max_alignment_score = score

    average_alignment_score /= len(alignment_counts)
    print(f"average alignment score = {average_alignment_score}")
    print(f"max alignment score = {max_alignment_score}")
    alignment_global_score_file = "alignment_dot_product.txt"
    with open(alignment_global_score_file, 'w') as F:
        F.write("PGFam\tdot_product\tpropMax\tlength\n")
        for pgfam in alignment_counts:
            score = alignment_global_score[pgfam] 
            prop = score / max_alignment_score
            F.write(f"{pgfam}\t{score:0.4f}\t{prop:0.4f}\t{alignment_length[pgfam]}\n")
    print(f"alignment global scores written to {alignment_global_score_file}")
    #sys.exit()


    alignment_names.append("total_alignment")
    alignment_vector = []
    #print(alignment)
    for bitstring in shared_bitstrings:
        pf = 0.0
        if bitstring in global_apart_freq:
            pf = global_apart_freq[bitstring]
        alignment_vector.append(pf)
    alignment_vector_list.append(alignment_vector)


    alignment_v_df = pandas.DataFrame(alignment_vector_list, index=alignment_names, columns=shared_bitstrings)

    tree_alignment_partition_matrix = pandas.concat([tree_v_df, alignment_v_df])
    tree_alignment_partition_matrix.fillna(0, inplace=True)
    if args.write_matrix:
        tree_alignment_partition_matrix.to_csv("tree_alignment_partition_matrix.txt", sep="\t")
        F = open("ordered_genomes.txt", "w")
        F.write("genomes\n"+"\n".join(taxon_list)+"\n")

    joint_umap = generate_umap_plot(tree_alignment_partition_matrix, n_neighbors=[5,10,20])
    cluster_umap_plot(joint_umap)

    joint_umap.to_csv(f"tree_alignment_bipartition_umap.txt", sep="\t")
    print(f"umap data written to tree_alignment_bipartition_umap.txt")

if __name__ == '__main__':
    exit(main())

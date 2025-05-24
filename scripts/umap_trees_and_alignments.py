import dendropy
import sys
import pandas
import glob
import re
import BitVector
import argparse
import math

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

def score_tree_partitions(tree, normalize=False):
    bem = tree.bipartition_edge_map
    total_length = 0 # can use to normalize partition scores to sum to 1.0
    for bipart in bem:
        if (bem[bipart].length is not None):
            total_length += bem[bipart].length
    partition_length = {}
    for bipart in bem:
        if (bem[bipart].length is not None):
            bitstring = bipart.leafset_as_bitstring(reverse=True)
            length = bem[bipart].length
            if normalize:
                lenght /= total_length # convert to proportion of total tree length
            partition_length[bitstring] = length 
    return (partition_length, total_length)

def count_alignment_partitions(al, seqid_list, normalize=False):
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
            increment = 1
            if normalize:
                increment = 1/len(bitvector_dict)
            partition_counts[bitstring] += increment
    return(partition_counts, al_length)

def generate_umap_plot(data, n_neighbors=[10,30,100]):
    # yields 2-D embeddings of data vectors
    print("generate_umap_plot()")
    import umap.umap_ as umap
    umap_df = pandas.DataFrame(index = data.index)
    data.fillna(0, inplace=True) # avoid errors around NaNs
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
    parser.add_argument("--write_matrix", action='store_true', help="write partition matrix and separate column list")
    parser.add_argument("--normalize_alignment_scores", action='store_true', help="adjust for multiple states, alignment scores sum to 1.0")
    parser.add_argument("--normalize_tree_scores", action='store_true', help="divide tree scores by tree length, tree scores sum to 1.0")

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
        (tree_scores, tree_length) = score_tree_partitions(one_tree, normalize = args.normalize_tree_scores)
        tree_bitstring_length[tree_label] = tree_scores
        for bitstring in tree_scores:
            if not bitstring in partition_count:
                partition_count[bitstring] = 0
                global_tree_partition_length[bitstring] = 0
            partition_count[bitstring] += 1
            global_tree_partition_length[bitstring] += tree_scores[bitstring]

    tree_bl_mean = 0
    unnormalized_sum = 0
    for bs in global_tree_partition_length:
        #global_tree_partition_length[bs] /= partition_count[bs]
        # how to normalize??? above line or below
        unnormalized_sum += global_tree_partition_length[bs]
        global_tree_partition_length[bs] /= len(tree_bitstring_length)
        tree_bl_mean += global_tree_partition_length[bs]
    print(f"sum of tree partition freq = {tree_bl_mean}")
    print(f"unnormalized sum of global tree values = {unnormalized_sum}")
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
    monomorphic_bs = '0'*len(taxon_list)
    alignment_monomorphic_count = {}
    for alignment_file in alignment_files:
        #print(f"read {alignment_file}")
        m = re.search("(PGF_\d+)", alignment_file)
        if m:
            pgfam = m.group(1)
        else:
            pgfam = alignment_file
        pgfam += "_"+alignment_type
        alignment_counts[pgfam] = {}
        with open(alignment_file) as F:
            al = read_fasta(F)
            pc, alength = count_alignment_partitions(al, taxon_list, normalize=args.normalize_alignment_scores)
            total_length += alength
            alignment_length[pgfam] = alength
            alignment_monomorphic_count[pgfam] = pc[monomorphic_bs]
            #print(f" num partitions = {len(pc)}")
            for bs in pc:
                if bs not in global_apart_counts:
                    global_apart_counts[bs] = 0
                global_apart_counts[bs] += pc[bs]
                alignment_counts[pgfam][bs] = (pc[bs] / alength)
    print(f"num bitstrings seen in alignments = {len(global_apart_counts)}")

    global_apart_freq = {}
    afreq_mean = 0
    num_in_tree = 0
    sum_in_tree = 0
    for bs in global_apart_counts:
        global_apart_freq[bs] = global_apart_counts[bs] / total_length
        afreq_mean += global_apart_freq[bs]
        if bs in global_tree_partition_length:
            num_in_tree += 1
            sum_in_tree += global_tree_partition_length[bs]
    print(f"global alignment partitions in tree = {num_in_tree}")
    print(f"sum of tree partiton lengths in alignments = {sum_in_tree}")
    print(f"global alignment sum freq = {afreq_mean}")
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
    alignment_self_score = {}
    alignment_dist = {}
    alignment_max_dev = {}
    alignment_vector_list = []
    alignment_names = []
    average_alignment_score = 0
    max_alignment_score = [0,0,0]
    alignment_freq_sum = {}
    for alignment in alignment_counts:
        alignment_names.append(alignment)
        alignment_vector = []
        #print(alignment)
        score_list = [0,0,0]
        self_list = [0,0,0]
        dist_list = [0,0,0]
        alignment_sum = 0
        max_dev = 0.0
        #for bitstring in global_apart_freq: #shared_bitstrings:
        for bitstring in shared_bitstrings:
            global_freq = global_apart_freq[bitstring]
            numon = 0
            for i in bitstring:
                if i == '1':
                    numon += 1
            pf = 0.0
            if bitstring in alignment_counts[alignment]:
                pf = alignment_counts[alignment][bitstring]
                for i in range(3):
                    if numon >= i: #segregate scores for fixed states, autapomorphies, and synapomorphies
                        score_list[i] += (pf * global_freq)  # sum of cross products = dot product = scalar product
                        self_list[i] += (pf * pf)  # self dot product == sum of squares
                alignment_sum += pf
            alignment_vector.append(pf)
            dist = (pf - global_freq) * (pf - global_freq)
            for i in range(3):
                if numon >= i: #segregate scores for fixed states, autapomorphies, and synapomorphies
                    dist_list[i] += dist
            if (pf - global_freq) > max_dev:
                max_dev = (pf - global_freq)

        alignment_vector_list.append(alignment_vector)
        alignment_global_score[alignment] = score_list
        alignment_self_score[alignment] = self_list
        alignment_dist[alignment] = dist_list
        alignment_max_dev[alignment] = max_dev
        alignment_freq_sum[alignment] = alignment_sum
        average_alignment_score += score_list[0]
        for i in range(3):
            if score_list[i] > max_alignment_score[i]:
                max_alignment_score[i] = score_list[i]

    average_alignment_score /= len(alignment_counts)
    print(f"average alignment score = {average_alignment_score}")
    print(f"max alignment score = {max_alignment_score}")
    alignment_global_score_file = "alignment_dot_product.txt"
    with open(alignment_global_score_file, 'w') as F:
        F.write("PGFam\tdprod0\tpMax0\tdprod1\tpMax1\tdprod2\tpMax2\tmaxDev\tdist0\tdist1\tdist2\tlength\tal_sum\tnum_patterns\tprop_fixed\n")
        for pgfam in alignment_counts:
            F.write(pgfam)
            for i in range(3):
                score = alignment_global_score[pgfam][i] 
                prop = score / max_alignment_score[i]
                F.write(f"\t{score:0.5f}\t{prop:0.3f}")
            F.write(f"\t{alignment_max_dev[pgfam]:0.5}")
            for i in range(3):
                #score = 0
                #if alignment_self_score[pgfam][i]:
                score = alignment_dist[pgfam][i] #alignment_global_score[pgfam][i] / alignment_self_score[pgfam][i] 
                #    if score > 1.0:
                #        score = 1.0
                #prop = score / max_alignment_score[i]
                F.write(f"\t{score:0.5f}")
            propFixed = alignment_monomorphic_count[pgfam]/alignment_length[pgfam]
            F.write(f"\t{alignment_length[pgfam]}\t{alignment_freq_sum[pgfam]:.2f}\t{len(alignment_counts[pgfam])}\t{propFixed:.3f}\n")
    print(f"alignment global scores written to {alignment_global_score_file}")


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

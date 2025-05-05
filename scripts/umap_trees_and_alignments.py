import dendropy
import sys
import pandas
import glob
import re
import BitVector

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


#sampleSize = sys.argv[1]
#tree_files = glob.glob(f"*_trees_sampleSize{sampleSize}.nex")
#tree_files = glob.glob(f"*_trees_sampleSize*.nex")
tree_file = sys.argv[1]

tree_partition_length = {}
tree_bitstring_length = {}
partition_count = {}
global_tree_partition_length = {}
tree_types = set()
taxa = dendropy.TaxonNamespace()
tree_yielder = dendropy.Tree.yield_from_files(
        files=[tree_file],
        schema='nexus',
        taxon_namespace=taxa,
        )
for tree_idx, one_tree in enumerate(tree_yielder):
    bem = one_tree.bipartition_edge_map
    print(one_tree.label+"    \t"+str(tree_idx))
    tree_type = one_tree.label.split(' ')[2]
    #tree_label = one_tree.label.replace(' ', '_')
    m = re.search("(PGF.\d+)", one_tree.label)
    tree_label = m.group(1) #pgfam
    tree_partition_length[tree_label] = {}
    tree_bitstring_length[tree_label] = {}
    tree_types.add(tree_type)
    for bipart in bem:
        if (bem[bipart].length is not None):
            bitstring = bipart.leafset_as_bitstring(reverse=True)
            tree_bitstring_length[tree_label][bitstring] = bem[bipart].length
            if not bitstring in partition_count:
                partition_count[bitstring] = 0
                global_tree_partition_length[bitstring] = 0
            partition_count[bitstring] += 1
            global_tree_partition_length[bitstring] += bem[bipart].length
    ##if tree_idx > 100:
    #    break
for bs in global_tree_partition_length:
    global_tree_partition_length[bs] /= partition_count[bs]
print(f"tree types = {tree_types}")
print(f"num trees = {len(tree_partition_length)}")
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
al_directory = sys.argv[2]
alignment_files = glob.glob(al_directory+"/PGF*")
for alignment_file in alignment_files:
    print(f"read {alignment_file}")
    m = re.search("(PGF_\d+)", alignment_file)
    pgfam = m.group(1)
    alignment_counts[pgfam] = {}
    with open(alignment_file) as F:
        al = read_fasta(F)
        pc, alength = count_alignment_partitions(al, taxon_list)
        total_length += alength
        print(f" num partitions = {len(pc)}")
        for bs in sorted(pc, key=pc.get):
            if bs not in global_apart_counts:
                global_apart_counts[bs] = 0
            global_apart_counts[bs] += pc[bs]
            alignment_counts[pgfam][bs] = (pc[bs] / alength)

global_apart_freq = {}
for bs in global_apart_counts:
    global_apart_freq[bs] = global_apart_counts[bs] / total_length
print(f"num global counts from alignments = {len(global_apart_counts)}")


F = open("Alignment_vs_tree_partition_values.txt", 'w')
i=0
F.write("partition\talign\ttree\n")
shared_bitstrings = set()
threshold_for_saving_alignment_partition = 5e-6
threshold_for_saving_tree_partition = 1e-3
for bs in sorted(global_apart_counts, key=global_apart_counts.get, reverse=True):
    i+=1
    #if i > 50:
    #   break
    #bipartition_list.append(bp)
    #partition_taxon_list = []
    #for i, val in enumerate(bs): # converts bs to list
    #    if val == '1':
    #        partition_taxon_list.append(taxon_list[i])
    tree_partition_length = 0
    if bs in global_tree_partition_length:
        tree_partition_length = global_tree_partition_length[bs]
        if (global_apart_freq[bs] > threshold_for_saving_alignment_partition) and (tree_partition_length > threshold_for_saving_tree_partition):
            shared_bitstrings.add(bs)
    F.write(f"{bs}\t{global_apart_counts[bs]}\t{global_tree_partition_length}\n")
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
    tree_names.append(tree.replace(' ', '_')+"_dna_tree")


tree_v_df = pandas.DataFrame(tree_vector_list, index=tree_names)
print(f"tree_v_df = {tree_v_df}")


alignment_vector_list = []
alignment_names = []
for alignment in alignment_counts:
    alignment_names.append(alignment+"_dna_alignment")
    alignment_vector = []
    #print(alignment)
    for bitstring in shared_bitstrings:
        pf = 0.0
        if bitstring in alignment_counts[alignment]:
            pf = alignment_counts[alignment][bitstring]
        alignment_vector.append(pf)
    alignment_vector_list.append(alignment_vector)

alignment_v_df = pandas.DataFrame(alignment_vector_list, index=alignment_names)

#alignment_umap = generate_umap_plot(alignment_v_df, n_neighbors=[5,10,20])
#cluster_umap_plot(alignment_umap)
#tree_umap = generate_umap_plot(tree_v_df, n_neighbors=[5,10,20])
#cluster_umap_plot(tree_umap)

joint_umap = generate_umap_plot(pandas.concat([tree_v_df, alignment_v_df]), n_neighbors=[5,10,20])
cluster_umap_plot(joint_umap)

#sys.exit()

#print(f"tree umap = {tree_umap}")
#tree_umap.to_csv(f"tree_bipartition_umap.txt", sep="\t")
#print(f"umap data written to tree_bipartition_umap.txt")


#alignment_umap.to_csv(f"alignment_bipartition_umap.txt", sep="\t")
#print(f"umap data written to alignment_bipartition_umap.txt")

joint_umap.to_csv(f"joint_bipartition_umap.txt", sep="\t")
print(f"umap data written to joint_bipartition_umap.txt")

import dendropy
import sys
import pandas
import glob
import re


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
tree_files = sys.argv[1:]

tree_partition_length = {}
all_partitions = set()
tree_types = set()
taxa = dendropy.TaxonNamespace()
tree_yielder = dendropy.Tree.yield_from_files(
        files=tree_files,
        schema='nexus',
        taxon_namespace=taxa,
        )
for tree_idx, one_tree in enumerate(tree_yielder):
    bem = one_tree.bipartition_edge_map
    #print(one_tree.label+"    \t"+str(tree_idx))
    tree_type = one_tree.label.split(' ')[2]
    tree_label = one_tree.label.replace(' ', '_')
    tree_partition_length[tree_label] = {}
    tree_types.add(tree_type)
    for bipart in bem:
        if (bem[bipart].length is not None):
            tree_partition_length[tree_label][bipart] = bem[bipart].length
            all_partitions.add(bipart)
    ##if tree_idx > 100:
    #    break
print(f"tree types = {tree_types}")
print(f"num trees = {len(tree_partition_length)}")
bipartition_list = list(all_partitions)
tree_vector_list = []
tree_names = []
for tree in tree_partition_length:
    tree_names.append(tree)
    tree_vector = []
    #print(tree)
    for bipart in bipartition_list:
        bl = 0.0
        if bipart in tree_partition_length[tree]:
            bl = tree_partition_length[tree][bipart]
        tree_vector.append(bl)
    tree_vector_list.append(tree_vector)

tree_v_df = pandas.DataFrame(tree_vector_list, index=tree_names)

print(f"tree_v_df = {tree_v_df}")

#generate_umap_plot(data, n_neighbors=[10,30,100]):
tree_umap = generate_umap_plot(tree_v_df, n_neighbors=[5,10,20])
cluster_umap_plot(tree_umap)

tree_umap['avg_site_ll'] = 0

for alignment_score_file in glob.glob("alignment_sample_scores_sampleSize*.txt"):
    F = open(alignment_score_file)
    m = re.match("alignment_sample_scores_sampleSize(\d+).txt", alignment_score_file)
    sampleSize = m.group(1)
    dna_likelihood_column, protein_likelihood_column, joint_likelihood_column = None, None, None 
    headers = F.readline().rstrip().split("\t")
    for i, column in enumerate(headers):
        if (column == 'dna_likelihood') and ('dna' in tree_types):
            dna_likelihood_column = i;
            print(f"dna_likelihood_column = {dna_likelihood_column}")
        if (column == 'protein_likelihood') and ('protein' in tree_types):
            protein_likelihood_column = i;
            print(f"protein_likelihood_column = {protein_likelihood_column}")
        if (column == 'joint_likelihood') and ('joint' in tree_types):
            joint_likelihood_column = i;
            print(f"joint_likelihood_column = {joint_likelihood_column}")
    for line in F:
        fields = line.rstrip().split("\t")
        sample_no = fields[0] #.split("_")[1]
        num_pos = float(fields[2])
        if dna_likelihood_column and ((sample_no + "_dna") in tree_names):
            dna_avg_site_likelihood = float(fields[dna_likelihood_column])/(num_pos * 3)
            tree_umap.loc[sample_no+"_dna"]['avg_site_ll'] = dna_avg_site_likelihood
            #print(f"for sample {sample_no}, dna_likelihood_column = {dna_likelihood_column}, table cell = {fields[dna_likelihood_column]}, dna_avg_ll = {dna_avg_site_likelihood}")
        if protein_likelihood_column and ((sample_no + "_protein") in tree_names):
            protein_avg_site_likelihood = float(fields[protein_likelihood_column])/num_pos
            tree_umap.loc[sample_no+"_protein"]['avg_site_ll'] = protein_avg_site_likelihood
        if joint_likelihood_column and ((sample_no + "_joint") in tree_names):
            joint_avg_site_likelihood = float(fields[joint_likelihood_column])/(num_pos * 4)
            tree_umap.loc[sample_no+"_joint"]['avg_site_ll'] = joint_avg_site_likelihood

print(f"tree umap = {tree_umap}")
tree_umap.to_csv(f"all_trees_bipartition_umap.txt", sep="\t")

print(f"umap data written to all_trees_bipartition_umap.txt")

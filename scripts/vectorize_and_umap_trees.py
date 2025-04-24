import dendropy
import sys
import pandas
import umap.umap_ as umap
import glob


def generate_umap_plot(data, n_neighbors=[10,30,100]):
    # yields 2-D embeddings of data vectors
    print("generate_umap_plot()")
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

sampleSize = sys.argv[1]
tree_files = glob.glob(f"*_trees_sampleSize{sampleSize}.nex")

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
    tree_partition_length[one_tree.label] = {}
    bem = one_tree.bipartition_edge_map
    #print(one_tree.label+"    \t"+str(tree_idx))
    tree_type = one_tree.label.split('_')[0]
    tree_types.add(tree_type)
    for bipart in bem:
        if (bem[bipart].length is not None):
            tree_partition_length[one_tree.label][bipart] = bem[bipart].length
            all_partitions.add(bipart)
    if tree_idx > 100:
        break
bipartition_list = list(all_partitions)
tree_vector_list = []
tree_names = []
for tree in tree_partition_length:
    tree_names.append(tree)
    tree_vector = []
    print(tree)
    for bipart in bipartition_list:
        bl = 0.0
        if bipart in tree_partition_length[tree]:
            bl = tree_partition_length[tree][bipart]
        tree_vector.append(bl)
    tree_vector_list.append(tree_vector)

tree_v_df = pandas.DataFrame(tree_vector_list, index=tree_names)

print(f"tree_v_df = {tree_v_df}")

tree_umap = generate_umap_plot(tree_v_df)

tree_umap['avg_site_ll'] = 0

F = open(f"alignment_sample_scores_sampleSize{sampleSize}.txt")
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
    sample_no = fields[0].split("_")[1]
    num_pos = float(fields[2])
    if dna_likelihood_column:
        dna_avg_site_likelihood = float(fields[dna_likelihood_column])/(num_pos * 3)
        tree_umap.loc["dna sample "+sample_no]['avg_site_ll'] = dna_avg_site_likelihood
        print(f"for sample {sample_no}, dna_likelihood_column = {dna_likelihood_column}, table cell = {fields[dna_likelihood_column]}, dna_avg_ll = {dna_avg_site_likelihood}")
    if protein_likelihood_column:
        protein_avg_site_likelihood = float(fields[protein_likelihood_column])/num_pos
        tree_umap.loc["protein sample "+sample_no]['avg_site_ll'] = protein_avg_site_likelihood
    if joint_likelihood_column:
        joint_avg_site_likelihood = float(fields[joint_likelihood_column])/(num_pos * 4)
        tree_umap.loc["joint sample "+sample_no]['avg_site_ll'] = joint_avg_site_likelihood

print(f"tree umap = {tree_umap}")
tree_umap.to_csv(f"samplesize{sampleSize}_bipartition_umap.txt", sep="\t")

print(f"umap data written to samplesize{sampleSize}_bipartition_umap.txt")

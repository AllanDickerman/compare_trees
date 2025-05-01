
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


infile = sys.argv[1]
df = pandas.read_csv(infile, sep="\t", index_col=0, header=None)
umapdf = generate_umap_plot(df)
cluster_umap_plot(umapdf)

umapdf.to_csv("umap_out.txt", sep="\t")

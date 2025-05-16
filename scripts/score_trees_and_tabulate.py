import dendropy
import sys
import re
import glob
import argparse
import statistics
from dendropy.calculate import treecompare

def compare_using_bitmaps(t1, t2):
    bi1 = t1.bipartition_edge_map
    bi2 = t2.bipartition_edge_map
    #print(f"tree1 has {len(bi1)} edges.")

    common = {}
    t1_length = 0
    t2_length = 0
    t1_shared_length = 0
    t2_shared_length = 0
    t1_shared = []
    t2_shared = []
    for b1 in bi1:
        if True or (not b1.is_trivial()):
            if (bi1[b1].length is not None):
                t1_length += bi1[b1].length
            if (b1 in bi2):
                common[b1] = bi2[b1]
                if (bi2[b1].length is not None) and (bi1[b1].length is not None):
                    t1_shared_length += bi1[b1].length
                    t2_shared_length += bi2[b1].length
                    t1_shared.append(bi1[b1].length)
                    t2_shared.append(bi2[b1].length)
                    #print(f"shared {bi1[b1].length:.5}\t{bi2[b1].length:.5}\tsb:{b1}")

    for b2 in bi2:
        if (bi2[b2].length is not None):
            t2_length += bi2[b2].length

    shared_cor = None
    if len(t1_shared) > 2:
        shared_cor = statistics.correlation(t1_shared, t2_shared)
    retval = [(t1_length/t2_length), (t1_shared_length/(t1_length*2))+(t2_shared_length/(t2_length*2)), shared_cor, len(t1_shared)/len(bi1), len(bi1)-len(t1_shared)]
    return retval

def process_tree_file(ref_tree, treefile):
    F = open(treefile)
    tree_string_list = F.read().split(";")
    stats = {}
    stats["length_ratio"] = []
    stats["shared_length_prop"] = []
    stats["shared_length_cor"] = []
    stats["shared_clades_prop"] = []
    stats["diff_clades"] = []
    stats["robinsonfoulds"] = []
    stats["rf_distance"] = []
    stats["euclidean"] = []
    for line in tree_string_list:
        if not '(' in line:
            continue # blank lines
        line += ';' # this got stripped off
        #print("tree string="+line)
        tree = dendropy.Tree.get(data=line, schema="newick", taxon_namespace=ref_tree.taxon_namespace) 
        comp_stats = compare_using_bitmaps(ref_tree, tree)
        stats["length_ratio"].append(comp_stats[0])
        stats["shared_length_prop"].append(comp_stats[1])
        stats["shared_length_cor"].append(comp_stats[2])
        stats["shared_clades_prop"].append(comp_stats[3])
        stats["diff_clades"].append(comp_stats[4])
        stats["robinsonfoulds"].append(treecompare.symmetric_difference(ref_tree, tree))
        stats["rf_distance"].append(treecompare.weighted_robinson_foulds_distance(ref_tree, tree))
        stats["euclidean"].append(treecompare.euclidean_distance(ref_tree, tree))
    return stats

def compare_trees_to_ref(tree_file, ref_tree):
    tree_yielder = dendropy.Tree.yield_from_files(
            files=[tree_file],
            schema='nexus',
            taxon_namespace=ref_tree.taxon_namespace,
            )
    stats = {'rfsym':[], 'rfdist':[], 'euclid':[]}
    for tree_idx, one_tree in enumerate(tree_yielder):
        tree_label = one_tree.label.replace(' ', '_')
        stats['rfsym'].append(treecompare.symmetric_difference(ref_tree, one_tree))
        stats['rfdist'].append(treecompare.weighted_robinson_foulds_distance(ref_tree, one_tree))
        stats['euclid'].append(treecompare.euclidean_distance(ref_tree, one_tree))
    return stats  

def main():
    parser = argparse.ArgumentParser(description="Compare trees to reference tree", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", metavar='files', type=str, nargs='*', help="nexus tree files")
    parser.add_argument("--ref_tree", metavar='file', type=str, help="reference tree")

    args = parser.parse_args()

    print(f"tree files = {args.trees}")
    print(f"ref tree = {args.ref_tree}")
    
    genomes = dendropy.TaxonNamespace()
    ref_tree = dendropy.Tree.get(path=args.ref_tree, schema='nexus', taxon_namespace=genomes) 

    print("sample\trfsymm\trfdist\teuclid")
    for tree_file in args.trees:
        stats = compare_trees_to_ref(tree_file, ref_tree)
        m = re.match("dna_ags_decile_0(\d)", tree_file)
        sys.stdout.write(f"Decile {m.group(1)}")
        for s in ('rfsym', 'rfdist', 'euclid'):
            #print(f"stats[{s}] = {stats[s]}")
            mean = statistics.mean(stats[s])
            stdev = statistics.stdev(stats[s])
            sys.stdout.write(f"\t{mean:0.3f}+-{stdev:.3f}")
        sys.stdout.write("\n")

    return(0)

if __name__ == '__main__':
    #print("Here")
    exit(main())

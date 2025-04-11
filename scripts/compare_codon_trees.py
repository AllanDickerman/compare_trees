import sys
import glob
import statistics
import re
import shutil
import subprocess
import os
import os.path


def run_treedist(intree, intree2):
    print(f"run_treedist({intree}, {intree2}")
    command = 'treedist'
    if os.path.exists("outfile"):
        os.remove("outfile")
    prompt_string = f"{intree}\n2\nL\nF\nY\n{intree2}\n"
    #print(f"input to process will be:\n{prompt_string}")
    devnull = open(os.devnull, 'w')
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=devnull, text=True)
    proc.communicate(input=prompt_string)
    proc.wait()
    if not os.path.exists("outfile"):
        print("treedist failed")
        return None
    print("treedist generated outfile!")
    F = open("outfile")
    tree_scores = []
    for line in F:
        m = re.match("\s+(\d+)\s+\|\s+([\d\.]+)\s*$", line)
        if m:
            tree_scores.append(float(m.group(2)))
        #else:
        #    print("no match for "+line, end="")
    print(f"Num tree scores found = {len(tree_scores)}")
    return tree_scores


ref_trees = {}
ref_tree_files = glob.glob("*_ref.nwk")
for tf in ref_tree_files:
    print(f"read ref tree {tf}")
    x = tf.split("_")
    data = x[0]
    ref_trees[data] = tf

data_trees = {}
tree_files = glob.glob("*_trees_sampleSize*nwk")
for tf in tree_files:
    print(tf)
    m = re.match("(\w+)_trees_sampleSize(\d+).nwk", tf)
    data = m.group(1)
    size = int(m.group(2))
    print(f"  data\t{data}\tsize\t{size}\t{tf}")
    if data not in data_trees:
        data_trees[data] = {}
    data_trees[data][size] = tf

ref_vs_sample = {}
sample_sizes = set()
for ref_data in ref_trees:
    ref_vs_sample[ref_data] = {}
    for sample_data in data_trees:   # ['dna', 'protein', 'joint' ]:
        ref_vs_sample[ref_data][sample_data] = {}
        for sample_size in data_trees[sample_data]:
            sample_sizes.add(sample_size)
            scores = run_treedist(data_trees[sample_data][sample_size], ref_trees[ref_data])
            if len(scores) > 1:
                mean_score = statistics.mean(scores)
                ref_vs_sample[ref_data][sample_data][sample_size] = mean_score
            else:
                raise Exception(f"hey! scores returned from run_treedist was len({len(scores)} \n   ref={ref_trees[ref_data]} \n   sample={data_trees[sample_data][sample_size]}\n")

F = open("sample_trees_vs_reference_scores.txt", 'w')
F.write("Ref Data >")
for ref_data in ref_trees:
    for sample_data in data_trees:
        F.write("\t"+ref_data)
F.write("\n")
F.write("Sample Data >")
for ref_data in ref_trees:
    for sample_data in data_trees:
        F.write("\t"+sample_data)
F.write("\n")
for sample_size in sorted(sample_sizes):
    F.write(str(sample_size))
    for ref_data in ref_trees:
        for sample_data in data_trees:
            if (ref_data in ref_vs_sample) and (sample_data in ref_vs_sample[ref_data]) and (sample_size in ref_vs_sample[ref_data][sample_data]): 
                F.write(f"\t{ref_vs_sample[ref_data][sample_data][sample_size]:.4}")
            else:
                F.write("\tna")
    F.write("\n")
F.close()
print("wrote data to sample_trees_vs_reference_scores.txt")



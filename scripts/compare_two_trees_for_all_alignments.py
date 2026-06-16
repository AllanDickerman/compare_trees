import sys
import subprocess
import os
import shutil
import os.path
import glob
import argparse
import re
import tempfile

def raxml_sh_test(alignment, tree_0, tree_1, temp_dir = '.', threads = 2):
    m = re.search("(PGF_\d+)", alignment)
    pgfam = m.group(1)
    # Create a temporary directory
    print(f"run_sh_test({pgfam})")
    subprocess.run("rm RAxML_info.sh-test", shell=True) # avoid clashing file names
    ret_val = {}
    #run raxml Shimodaira-Hasegawa log-ratio test of best tree versus reference trees
    raxml_command = ['raxml', '-s', alignment, '-n', 'sh-test', '-f', 'h', '-p', '1234', '-T', str(threads), '-m', 'GTRGAMMA', '-t', tree_0, '-z', tree_1]
    print(f"run {raxml_command}")
    subprocess.run(raxml_command, capture_output=True)
    tree_likelihood = 0
    with open("RAxML_info.sh-test") as F:
        for line in F:
            m = re.match("Model optimization, best Tree: (\S+)", line)
            if m:
                #print("got tree likelihood: "+line)
                ret_val["tree_0_LL"] = float(m.group(1))
                #ret_val['best_likelihood'] = tree_likelihood
            m = re.match("Tree: (\S+) Likelihood: (\S+) D\(LH\): (\S+) SD: (\S+)", line)
            if m:
                index = int(m.group(1))
                print("got comparison: "+line)
                ret_val["tree_1_LL"] = float(m.group(2))
                ret_val["dLL"] = float(m.group(3))
                ret_val["SD"] = float(m.group(4))

    print(f"ret_val: {ret_val}")
    return ret_val  


if __name__ == '__main__':    
    parser = argparse.ArgumentParser(description="Test two trees against each PGFam alignment", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--alignment_dir", metavar="path", type=str, help="path where dna alignments are found")
    parser.add_argument("--tree_0", metavar="file", type=str, help="file with first tree (newick)")
    parser.add_argument("--tree_1", metavar="file", type=str, help="file with second tree (newick)")
    parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for tree building")
    parser.add_argument("--max_alignments", metavar="max", type=int, default=0, help="stop after this many trees")
    args = parser.parse_args()

    alignment_files = []
    for file in glob.glob(f"{args.alignment_dir}/PGF*afn"):
        alignment_files.append(os.path.abspath(file))

    #temp_dir = tempfile.TemporaryDirectory(prefix="test_trees_by_logRatio_", delete= not args.debug)
    temp_dir = "./tmp"
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        print('created temporary directory', temp_dir)

    shutil.copy(args.tree_0, temp_dir)
    tree_0_file = os.path.basename(args.tree_0)
    tree_0_name = tree_0_file[:-4]
    shutil.copy(args.tree_1, temp_dir)
    tree_1_file = os.path.basename(args.tree_1)
    tree_1_name = tree_1_file[:-4]

    original_dir = os.getcwd()
    os.chdir(temp_dir)

    i = 0
    data = {}
    for alignment_file in alignment_files:
        i += 1
        print(f"round {i}")
        if args.max_alignments and i > args.max_alignments:
            print(f" hit max_alignments {args.max_alignments}")
            break
        m = re.search("/(P.F_\d+)", alignment_file)
        pgfam = m.group(1)
        subprocess.run("rm RAxML*", shell=True) # start fresh
        data[pgfam] = raxml_sh_test(alignment_file, tree_0_file, tree_1_file, args.threads)

    os.chdir(original_dir)
    print("write results to file")
    out = open("LogRatio_test_results.txt", 'w')
    print(f"PGFam\t{tree_0_name}_LL\t{tree_1_name}_LL\tdLL\tSD\tz-score", file=out)
    for pgfam in data:
        print(f"{pgfam}", end='', file=out)
        print(f"\t{data[pgfam]['tree_0_LL']:.0f}", end='', file=out)
        print(f"\t{data[pgfam]['tree_1_LL']:.0f}", end='', file=out)
        print(f"\t{data[pgfam]['dLL']:.1f}", end='', file=out)
        print(f"\t{data[pgfam]['SD']:.1f}", end='', file=out)
        z_score = 0
        if data[pgfam]['SD'] > 1: 
            z_score = data[pgfam]['dLL']/data[pgfam]['SD']
        print(f"\t{z_score:.3f}", end='\n', file=out)
    print("", file=out)

    out.close()
    sys.exit(0)

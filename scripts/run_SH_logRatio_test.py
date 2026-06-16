import sys
import subprocess
import os
import os.path
import os.path
import shutil
import glob
import argparse
import re
import tempfile

def raxml_best_tree_dna(alignment, name="ml_tree", threads=2):
    # run regular raxml best tree search
    # -f d new rapid hill-climbing
    raxml_command = ['raxml', '-s', alignment, '-n', name, '-f', 'd', '-p', '1234', '-T', str(threads), '-m', 'GTRGAMMA']
    print(f"run {raxml_command}")
    proc = subprocess.run(raxml_command, capture_output=True)
    return f"RAxML_bestTree.{name}"

def raxml_constrained_tree(alignment, constraint_tree, name="constrained_tree", threads=2):
    # run 
    raxml_command = ['raxml', '-s', alignment, '-n', name, '-g', constraint_tree, '-p', '1234', '-T', str(threads), '-m', 'GTRGAMMA']
    print(f"run {raxml_command}")
    proc = subprocess.run(raxml_command, capture_output=True)
    #if os.path.exists("RAxML_bestTree.best_tree"):
    return f"RAxML_bestTree.{name}"

def raxml_sh_test(alignment, test_tree, alt_trees, alt_names, temp_dir = '.', threads = 2):
    m = re.search("(PGF_\d+)", alignment)
    pgfam = m.group(1)
    # Create a temporary directory
    print(f"run_sh_test({pgfam})")
    print(f"alt_names = {alt_names}")
    subprocess.run("rm RAxML_info.sh-test", shell=True) # avoid clashing file names
    ret_val = {}
    #print("now run ls:")
    #subprocess.run(['ls'])
    #print("now run ls RAxML*:")
    #run raxml Shimodaira-Hasegawa log-ratio test of best tree versus reference trees
    # use -f H to signal re-optimize branch lengths for specified alignment
    raxml_command = ['raxml', '-s', alignment, '-n', 'sh-test', '-f', 'h', '-p', '1234', '-T', str(threads), '-m', 'GTRGAMMA', '-t', test_tree, '-z', alt_trees]
    print(f"run {raxml_command}")
    subprocess.run(raxml_command, capture_output=True)
    tree_likelihood = 0
    with open("RAxML_info.sh-test") as F:
        for line in F:
            m = re.match("Model optimization, best Tree: (\S+)", line)
            if m:
                #print("got tree likelihood: "+line)
                tree_likelihood = float(m.group(1))
                #ret_val['best_likelihood'] = tree_likelihood
            m = re.match("Tree: (\S+) Likelihood: (\S+) D\(LH\): (\S+) SD: (\S+)", line)
            if m:
                index = int(m.group(1))
                alt_name = alt_names[index]
                print(f"index= {index}, alt name = {alt_names[index]}")
                print("got reftree "+line)
                ret_val[alt_name] = {}
                ret_val[alt_name]['best_likelihood'] = tree_likelihood
                ret_val[alt_name]["comp_likelihood"] = float(m.group(2))
                ret_val[alt_name]["DLH"] = float(m.group(3))
                ret_val[alt_name]["SD"] = float(m.group(4))
                ret_val[alt_name]["alpha"] = "NA"
                m2 = re.match(".*Significantly Worse:.*Yes \((\d+%)", line)
                if m2:
                    ret_val[alt_name]["alpha"] = m2.group(1)

    print(f"ret_val:")
    for an in alt_names:
        print(f"{an}:\t{ret_val[an]}")
    return ret_val  


if __name__ == '__main__':    
    parser = argparse.ArgumentParser(description="Test each PGFam alignment against reference trees", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--alignment_dir", metavar="path", type=str, help="path where dna alignments are found")
    parser.add_argument("--ref_trees", metavar="file", type=str, nargs='*', help="file with reference trees to compare to")
    parser.add_argument("--base_constraint", metavar="file", type=str, help="mildly constrained tree.")
    parser.add_argument("--strong_constraints", metavar="file(s)", type=str, nargs='*', help="reference trees to compare to")
    parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for tree building")
    parser.add_argument("--max_alignments", metavar="max", type=int, default=0, help="stop after this many trees")
    parser.add_argument("--debug", action="store_true", help="perserve temp dir")

    args = parser.parse_args()

    #temp_dir = tempfile.TemporaryDirectory(prefix="test_trees_by_logRatio_", delete= not args.debug)
    temp_dir = "./tmp"
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        print('created temporary directory', temp_dir)

    #ref_path = os.path.abspath(args.refTrees)
    alignment_files = []
    for file in glob.glob(f"{args.alignment_dir}/PGF*afn"):
        alignment_files.append(os.path.abspath(file))
    if args.base_constraint:
        constraint_trees = {}
        constraint_names = []
        base_constraint_name = os.path.basename(args.base_constraint)
        base_constraint_name = re.sub(".nwk", "", base_constraint_name)
        base_constraint_name = re.sub("constraint_?", "", base_constraint_name)
        constraint_trees[base_constraint_name] = os.path.abspath(args.base_constraint)
        constraint_names.append(base_constraint_name)
        stronger_trees = ['ml_tree', base_constraint_name]

        for file in args.strong_constraints:
            name = os.path.basename(file)
            name = name[:-4]
            name = re.sub("constraint_?", "", name)
            constraint_trees[name] = os.path.abspath(file)
            constraint_names.append(name)
    elif (args.ref_trees):
        stronger_trees = ['ml_tree']
        subprocess.run(f"cat {' '.join(args.ref_trees)} > {temp_dir}/ref_trees.nwk", shell=True)
        ref_names = []
        ref_trees = []
        for file in args.ref_trees:
            ref_trees.append(file)
            ref_name = os.path.basename(file)
            ref_name = re.sub(".nwk", "", ref_name)
            ref_names.append(ref_name)
            print(f"parsed {args.ref_trees} to names: {ref_names}")
        subprocess.run(f"rm {temp_dir}/ref_trees.nwk", shell=True)
        subprocess.run(f"cat {' '.join(ref_trees)} > {temp_dir}/ref_trees.nwk", shell=True)

    original_dir = os.getcwd()
    os.chdir(temp_dir)
    #tmpdirname = tempfile.mkdtemp(prefix="assess_codon_tree")
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
        data[pgfam] = {}
        subprocess.run("rm RAxML*", shell=True) # start fresh
        ml_tree = raxml_best_tree_dna(alignment_file)
        os.rename(ml_tree, f"{pgfam}_ml_tree.nwk")
        ml_tree = f"{pgfam}_ml_tree.nwk"

        if (args.base_constraint): # do constrained tree comparison
            print(f"create constrained trees")
            base_constraint_tree_file = raxml_constrained_tree(alignment_file, constraint_trees[base_constraint_name], name=base_constraint_name)
            base_strong_constraint_trees = "base_strong_constraint_trees.nwk"
            subprocess.run(f"cp {base_constraint_tree_file} {base_strong_constraint_trees}", shell=True)
            strong_constraint_trees = "strong_constraint_treees.nwk" # 
            subprocess.run(f"rm {strong_constraint_trees}", shell=True)
            for name in constraint_names[1:]: # the ones after the base constraint
                sctree = raxml_constrained_tree(alignment_file, constraint_trees[name], name=name)
                subprocess.run(f"grep Final {sctree}", shell=True)
                subprocess.run(f"cat {sctree} >> {strong_constraint_trees}", shell=True)
                subprocess.run(f"cat {sctree} >> {base_strong_constraint_trees}", shell=True)
            print(f"done writting constrained trees to {strong_constraint_trees} and {base_strong_constraint_trees}")
            subprocess.run(f"uniq {sctree} | wc", shell=True)
            
            ml_vs_const_sh = raxml_sh_test(alignment_file, ml_tree, base_strong_constraint_trees, constraint_names, threads = args.threads)
            #def raxml_sh_test(alignment, test_tree, alt_trees, alt_names, temp_dir = '.', threads = 2):
            data[pgfam]["ml_tree"] = ml_vs_const_sh
            base_vs_strong_sh = raxml_sh_test(alignment_file, base_constraint_tree_file, strong_constraint_trees, constraint_names[1:], threads = args.threads)
    #def raxml_sh_test(alignment, test_tree, alt_trees, alt_names, temp_dir = '.', threads = 2):
            data[pgfam][constraint_names[0]] = base_vs_strong_sh
        
        elif (args.ref_trees): # compare to fixed trees with fixed branch lengths
            ml_vs_ref_sh = raxml_sh_test(alignment_file, ml_tree, "ref_trees.nwk", ref_names, threads = args.threads)
            data[pgfam]["ml_tree"] = ml_vs_ref_sh

    os.chdir(original_dir)
    pgfam_sd = {}
    outfile = "SH_test_vs_constrained_trees.txt"
    if args.ref_trees: # not constrained trees, fixed trees
        outfile = "SH_test_vs_fixed_trees.txt"
    print(f"write results to {outfile}")
    out = open(outfile, 'w')
    print("PGFam\tstrong\tcomp_tree\tstrongLL\tcompLL\tSD\tDHL\talpha", file=out)
    for pgfam in data:
        pgfam_sd[pgfam] = 0
        i = 0
        for stronger_tree in stronger_trees:
            for comp_tree in data[pgfam][stronger_tree]:
                #print(f"data[{pgfam}][{comp_set}][{key}] = {data[pgfam][comp_set][key]}")
                print(f"{pgfam}\t{stronger_tree}\t{comp_tree}", end='', file=out)
                print(f"\t{data[pgfam][stronger_tree][comp_tree]['best_likelihood']:.0f}", end='', file=out)
                print(f"\t{data[pgfam][stronger_tree][comp_tree]['comp_likelihood']:.0f}", end='', file=out)
                print(f"\t{data[pgfam][stronger_tree][comp_tree]['SD']:.1f}", end='', file=out)
                print(f"\t{data[pgfam][stronger_tree][comp_tree]['DLH']:.1f}", end='', file=out)
                print(f"\t{data[pgfam][stronger_tree][comp_tree]['alpha']}", file=out)
                pgfam_sd[pgfam] += data[pgfam][stronger_tree][comp_tree]['SD']
                i+= 1
        pgfam_sd[pgfam] /= i
    out.close()

    outfile = re.sub(".txt", "_zscores.txt", outfile)
    print(f"write z-score results to {outfile}")
    out = open(outfile, 'w')
    print("PGFAM", end='', file=out)
    for stronger_tree in stronger_trees:
        for comp_tree in data[pgfam][stronger_tree]:
            print(f"\t{stronger_tree}.v.{comp_tree}", end='', file=out)
    print("", file=out)
    for pgfam in data:
        print(pgfam, end='', file=out)
        for stronger_tree in stronger_trees:
            for comp_tree in data[pgfam][stronger_tree]:
                #print(f"data[{pgfam}][{comp_set}][{key}] = {data[pgfam][comp_set][key]}")
                std_dev = data[pgfam][stronger_tree][comp_tree]['SD']
                if std_dev < pgfam_sd[pgfam]:
                    std_dev = pgfam_sd[pgfam] # use more conservative estimate
                zscore = data[pgfam][stronger_tree][comp_tree]['DLH']/std_dev
                if (abs(zscore) > 30) and (std_dev < 1):
                    zscore = -1
                print(f"\t{zscore:.2}", end='', file=out)
        print("", file=out)
    out.close()
    sys.exit(0)


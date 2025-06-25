
import sys
import re
import glob
import math
import phylocode
from random import sample
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio import BiopythonExperimentalWarning
#import warnings
#with warnings.catch_warnings():
#    warnings.simplefilter('ignore', BiopythonExperimentalWarning)

align_dir = "dna_alignments/"


d_calculator = DistanceCalculator('identity')
comp_dmat = {}
dmat_dict = {}
pgfam_alstats = {}
pgfam_list = []
perGenomeLogFreq = {}
genomeAlignmentCount = {}
for file in glob.glob(align_dir+"*"):
    print(f"file = {file}")
    m = re.search("(PGF_\d+)", file)
    if m:
        pgfam = m.group(1)
    else:
        continue #skip this non-pgfam file
    pgfam_list.append(pgfam)
    alignment = AlignIO.read(file, "fasta")
    stats = phylocode.calcAlignmentStats(alignment)
    pgfam_alstats[pgfam] = stats
    for genomeId in stats['perseq_meanlogfreq']:
        #genomeId = genomeIdFromFigId(figId)
        #print(f"looking at genome {genomeId}")
        if not genomeId in perGenomeLogFreq:
            perGenomeLogFreq[genomeId] = 0
            genomeAlignmentCount[genomeId] = 0
        perGenomeLogFreq[genomeId] += stats['perseq_meanlogfreq'][genomeId]
        genomeAlignmentCount[genomeId] += 1

    dmat = d_calculator.get_distance(alignment)
    dmat_dict[pgfam] = dmat

    for i, id1 in enumerate(dmat.names):
        for j, id2 in enumerate(dmat.names):
            if id1 == id2:
                break
            pair = tuple(sorted([id1, id2]))
            if pair not in comp_dmat:
                comp_dmat[pair] = {}
            comp_dmat[pair][pgfam] = dmat.matrix[i][j]

    if len(pgfam_list) > 200:
        break

for genomeId in perGenomeLogFreq:
    if genomeAlignmentCount[genomeId]:
        perGenomeLogFreq[genomeId] /= genomeAlignmentCount[genomeId]



print(f"number of pgfams processed = {len(dmat_dict)}")
print("calculate average distance matrix")
for pair in comp_dmat:
    dsum = 0
    num = 0
    for pgfam in pgfam_list:
        if pgfam in comp_dmat[pair]:
            dsum += comp_dmat[pair][pgfam]
            num += 1
    avgdist = dsum/num
    comp_dmat[pair]['avg'] = avgdist

pgfam_taxon_delta = {}
pgfam_delta = {}
pgfam_avg_dist = {}
for pgfam in pgfam_list:
    #sys.stdout.write("{pgfam}\n")
    pgfam_delta[pgfam] = 0
    pgfam_taxon_delta[pgfam] = {}
    pgfam_avg_dist[pgfam] = 0
    dmat = dmat_dict[pgfam]
    for id1 in dmat.names:
        pgfam_taxon_delta[pgfam][id1] = 0
    for i, id1 in enumerate(dmat.names):
        #sys.stdout.write(f"\n id1({i})")
        for j, id2 in enumerate(dmat.names):
            if id1 == id2:
                break
            #sys.stdout.write(f"  {j}")
            pair = tuple(sorted([id1, id2]))
            #if not pair in pgfam_delta:
            #    pgfam_delta[pgfam] = 0
            pgfam_avg_dist[pgfam] += dmat.matrix[i][j]
            delta = comp_dmat[pair]['avg'] - dmat.matrix[i][j]
            delta2 = (delta * delta)
            pgfam_taxon_delta[pgfam][id1] += delta2
            pgfam_taxon_delta[pgfam][id2] += delta2
            pgfam_delta[pgfam] += delta2
            #sys.stdout.write(f" {comp_dmat[pair][pgfam]},{dmat.matrix[i][j]},{delta}")
    #print(f"\tf:{pgfam_delta[pgfam]}")
    pgfam_delta[pgfam] = math.sqrt(pgfam_delta[pgfam])
    pgfam_delta[pgfam] /= (len(dmat.names)*(len(dmat.names)-1))
    pgfam_avg_dist[pgfam] /= (len(dmat.names)*(len(dmat.names)-1))
    for id1 in dmat.names:
        pgfam_taxon_delta[pgfam][id1] = math.sqrt(pgfam_taxon_delta[pgfam][id1])
        pgfam_taxon_delta[pgfam][id1] /= len(dmat.names)

pgfam_msf = {}
pgfam_worst_seq_score = {}
alignmentStatsFile = glob.glob("detail_files/*.homologAlignmentStats.txt")[0]
F = open(alignmentStatsFile)
F.readline()
for line in F:
    if line.startswith("PGF_"):
        fields = line.strip().split("\t")
        pgfam = fields[0]
        msf = float(fields[3])
        worst_seq_score = float(fields[6])
        pgfam_msf[pgfam] = msf
        pgfam_worst_seq_score[pgfam] = worst_seq_score

print("PGFam\tDeltaMS\tAvgDist\tWorstS\tWorstDelta\tWorstTax")
for pgfam in pgfam_list:
    dmat = dmat_dict[pgfam]
    worst_taxa = sorted(pgfam_taxon_delta[pgfam], key=pgfam_taxon_delta[pgfam].get, reverse=True)
    wt = worst_taxa[0]
    print(f"{pgfam}\t{pgfam_delta[pgfam]:.5f}\t{pgfam_avg_dist[pgfam]:.5f}\t{pgfam_worst_seq_score[pgfam]:.3f}\t{pgfam_taxon_delta[pgfam][wt]:.4f}\t{wt}")
    #print(f"{pgfam}\t{pgfam_delta[pgfam]:.5f}\t{pgfam_delta[pgfam] / pgfam_avg_dist[pgfam]:5f}\t{pgfam_worst_seq_score[pgfam]:.3f}")
    wts = "  worst taxa: "
    for t in worst_taxa[:4]:
        wts += f"\t{t}:{pgfam_taxon_delta[pgfam][t]:.4f}:{pgfam_taxon_delta[pgfam][t]/pgfam_avg_dist[pgfam]:.5f}"
    #print(wts)

F = open("test_al_dist_vs_log_freq.txt", 'w')
F.write("PGFam\tGenome\tDeltaDist\tlogFreqRat\tlogFRatE\n")
for pgfam in pgfam_list:
    for genome in perGenomeLogFreq:
        delta_dist = 0
        logFreqRatio = 0
        logFreqRatio = 0
        logFreqRatioE = 0
        if genome in pgfam_taxon_delta[pgfam]:
            delta_dist = pgfam_taxon_delta[pgfam][genome]
            pgfam_genome_logfreq = pgfam_alstats[pgfam]['perseq_meanlogfreq'][genome]
            genome_avg_logfreq = perGenomeLogFreq[genome]
            logFreqRatio = pgfam_genome_logfreq / genome_avg_logfreq
            logFreqRatioE = pgfam_genome_logfreq / (genome_avg_logfreq - 1e-1)
            #logFreqRatio = pgfam_alstats[pgfam]['perseq_meanlogfreq'][genome] / perGenomeLogFreq[pgfam][genome]
        F.write(f"{pgfam}\t{genome}\t{delta_dist:.5f}\t{logFreqRatio:.5f}\t{logFreqRatioE:.5f}\n")
F.close()
print("delta_dist vs logfreqRat data in test_al_dist_vs_log_freq.txt")


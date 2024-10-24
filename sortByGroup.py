import sys
from collections import OrderedDict

if(len(sys.argv) < 2):
    print("Usage: Python3 sortByGroup.py Variants.[groupname].coverage.Select.Seqs.Rows.Classified.tab.reconstructed")                                    
    exit()

# python3 sortByGroup.py Variants.ont.kg.hf1.hf2.1.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed

ifile = open(sys.argv[1], 'r')

output_files = {}
counts = {}
insdel_counts = {}

for i in range(1, 6):
    output_files[str(i)] = open(sys.argv[1] + ".group" + str(i) + ".only", "w")
    counts[str(i)] = 0
    insdel_counts[str(i)] = {}

for i in range(1, 6):
    for j in range(i + 1, 6):
        id = str(i) + str(j)
        output_files[id] = open(sys.argv[1] + ".group" + str(i) + ".group" + str(j), "w")
        counts[id] = 0
        insdel_counts[id] = {}

output_files["many"] = open(sys.argv[1] + ".morethan3groups", "w")
output_files["empty"] = open(sys.argv[1] + ".noHaps", "w")
counts["many"] = 0
counts["empty"] = 0
insdel_counts["many"] = {}
insdel_counts["empty"] = {}

def increment(dict, key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1



for full_line in ifile.readlines():
    line = full_line.split('\t')
    hap_seq = line[10].strip().split(' ')
    insdel_hap_seq = line[11].strip().split(' ')

    total_haps = 0

    groups = {}
    insdel_haps = {}

    for insdel_hap in insdel_hap_seq:
        increment(insdel_haps, insdel_hap)
        total_haps += 1
    
    for hap in hap_seq:
        group = hap[0]
        if group != "-" and group != "_":
            groups[group] = "present"
    id = ""
    for group in sorted(groups.keys()):
        id += group

    if len(id) > 2:
        id = "many"

    if len(id) == 0:
        id = "empty"

    if id in output_files:
        output_files[id].write(full_line)
        counts[id] += 1

        for insdel_hap in insdel_haps:
            increment(insdel_counts[id], insdel_hap)
            increment(insdel_counts[id], "total")
            
    else:
        print("unknown id: ", id)


GENOMIC_COVERAGE = 30
# it is about 30 per NOR region

stats_file = open(sys.argv[1] + ".stats.tsv", "w")
stats_file.write("Group(s)\tRead_Count\tGenomic_Coverage\tAvgNHap\tFeatureLen\tInsDelHapCounts\n")
for group in sorted(insdel_counts.keys()):
    num_reads = counts[group]
    if "total" not in insdel_counts[group]:
        insdel_counts[group]["total"] = 0
    avg_num_hap = 0
    if (num_reads != 0):
        avg_num_hap = insdel_counts[group]["total"] / counts[group]
    feature_len = avg_num_hap * ((counts[group] / GENOMIC_COVERAGE) - 1)

    stats_file.write(group.ljust(5) + "\t" + str(num_reads).ljust(5) + "\t" + str(GENOMIC_COVERAGE) + "\t" + str(avg_num_hap)[:7].ljust(7) + "\t" + str(feature_len)[:7].ljust(7) + "\tInsdelhaplotypes:")
    
    sorted_insdel_counts = OrderedDict(sorted(insdel_counts[group].items(), key = lambda item: item[1], reverse = True))

    for insdel_hap in sorted_insdel_counts:
        stats_file.write("\t" + insdel_hap + ":" + str(sorted_insdel_counts[insdel_hap]).ljust(4))

    stats_file.write("\n")

total = 0
for group in output_files:
    total += counts[group]

stats_file.write("TOTAL\t" + str(total))


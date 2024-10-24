import sys
from collections import OrderedDict
import re

# python3 NORclassifier.py /share/rwmstore/Assemblies_RWM/Lettuce/Lsat_V14/Reads_Ribo_NOR_CLC_Variants/Variants_Classification/Variants.Lsat_NOR_Segments_13864.X.Variants.Rows.Seqs /share/rwmstore/Assemblies_RWM/Lettuce/Lsat_V14/Reads_Ribo_NOR_CLC_Variants/Variants_Classification/Variants.Lsat_NOR_Segments_90.X.tab
# python3 NORclassifier.py ../alignments/Variants.hifi3.Select.Seqs.Rows Haplotypes.hifi3.125.tab



if(len(sys.argv) < 3):
    print("Usage: Python3 NORclassifier.py Variants.[readgroup].Select.Seqs.Rows Haplotypes.[readgroup].[count].tab Wildcards.tab")
    print(" Variants: nuelcotide info (trace) for new segments:")
    print("     [readname] [trace]")
    print(" Haplotypes: trace info for reference haplotypes:")
    print("     [haplotype] [trace]")
    print(" Wildcards: trace info for wildcard haplotypes:")
    print("     [haplotype] [trace]")
    exit()

read_variants = open(sys.argv[1],'r')
reference_variants = open(sys.argv[2],'r')
wildcard_variants = open(sys.argv[3], 'r')

output_file = open(sys.argv[1].split('/')[-1] + ".Classified.tab", 'w')
# output_file.write("Read_Name\tHaplotype\tTrace\n")

stats_file = open(sys.argv[1].split('/')[-1] + ".STATS.tab", 'w')
stats_file.write("Haplotype\tCount\tPercentage\n")

# 1. Load dictionary queried by haplotype into memory
reference_haplotype_names = OrderedDict({})     # trace      -> haplotype
wildcard_haplotype_names = OrderedDict({})      # regex trace-> haplotype
reference_haplotype_count = OrderedDict({})                  # haplotype  -> count
haplotype_name_len = 0
for line in reference_variants.readlines():
    line = line.split()
    reference_haplotype_names[line[1]] = line[0]
    reference_haplotype_count[line[0]] = 0
    reference_haplotype_count[line[0].lower()] = 0
    haplotype_name_len = len(line[0])

for line in wildcard_variants.readlines():
    line = line.split()
    wildcard_haplotype_names[line[1]] = line[0]
    reference_haplotype_count[line[0]] = 0

no_haplotype_name = "-"
while len(no_haplotype_name) < haplotype_name_len:
    no_haplotype_name += "-"
reference_haplotype_count[no_haplotype_name] = 0

# 2. Categorize reads
total_reads = 0
for line in read_variants.readlines():
    line = line.split()
    # default is uncategorized
    read_haplotype = no_haplotype_name
    read_trace = line[1]
    read_trace = read_trace.replace("N", ".")

    # 1. PERFECT MATCH
    # check for haplotype categorization
    for ref_trace in reference_haplotype_names:
        if re.fullmatch(read_trace, ref_trace):
            # read_haplotype = ref_trace
            read_haplotype = reference_haplotype_names[ref_trace]

    # if read_trace in reference_haplotype_names.keys():
        

    # 2. SOFT MATCH
    # parse each haplotype and check if it is close enough to it (i.e. differs by at most 1 nucleotide)
    if read_haplotype == no_haplotype_name:
        for ref_trace in reference_haplotype_names:
            mismatch = 0
            for (i, nucleotide) in enumerate(ref_trace):
                if (read_trace[i] != nucleotide):
                    mismatch += 1
                    if mismatch > 2:
                        break
            
            if (mismatch < 2):
                read_haplotype = reference_haplotype_names[ref_trace].lower()
                break
    
    # 3. WILDCARD MATCH
    # check if matches with a wildcard haplotype
    if read_haplotype == no_haplotype_name:
        for wildcard in wildcard_haplotype_names:
            if re.fullmatch(wildcard, read_trace):
                read_haplotype = wildcard_haplotype_names[wildcard]
    
    output_file.write(line[0] + '\t' + read_haplotype  + '\t' + read_trace.strip() + '\n')
    reference_haplotype_count[read_haplotype] += 1
    total_reads += 1

# 3. Write stats

all_haplotype_names = reference_haplotype_names.copy()
all_haplotype_names.update(wildcard_haplotype_names)

for haplotype in reference_haplotype_count:
    count = reference_haplotype_count[haplotype]
    percentage = 100 * count / total_reads
    stats_file.write(haplotype + "\t" + str(count) + "\t" + str(percentage)[0:10] + '\n')
# for haplotype in wildcard_haplotype_names:

# count = reference_haplotype_count[no_haplotype_name]
# percentage = 100 * count / total_reads
# stats_file.write(no_haplotype_name + "\t" + str(count) + "\t" + str(percentage)[0:10] + '\n')

stats_file.write("TOTAL\t" + str(total_reads) + "\t100")
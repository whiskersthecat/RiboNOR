import sys
from collections import OrderedDict

# python3 NORclassifier.py /share/rwmstore/Assemblies_RWM/Lettuce/Lsat_V14/Reads_Ribo_NOR_CLC_Variants/Variants_Classification/Variants.Lsat_NOR_Segments_13864.X.Variants.Rows.Seqs /share/rwmstore/Assemblies_RWM/Lettuce/Lsat_V14/Reads_Ribo_NOR_CLC_Variants/Variants_Classification/Variants.Lsat_NOR_Segments_90.X.tab
# python3 NORclassifier.py ../alignments/Variants.hifi3.Select.Seqs.Rows Haplotypes.hifi3.125.tab


if(len(sys.argv) < 2):
    print("Usage: Python3 NORclassifier.py Variants.[readgroup].Select.Seqs.Rows Haplotypes.[readgroup].[count].tab")
    print(" Variants: nuelcotide info (trace) for new segments:")
    print("     [readname] [trace]")
    print(" Haplotypes: trace info for reference haplotypes:")
    print("     [haplotype] [trace]")
    exit()

read_variants = open(sys.argv[1],'r')
reference_variants = open(sys.argv[2],'r')

output_file = open(sys.argv[1].split('/')[-1] + ".Classified.tab", 'w')
# output_file.write("Read_Name\tHaplotype\tTrace\n")

stats_file = open(sys.argv[1].split('/')[-1] + ".STATS.tab", 'w')
stats_file.write("Haplotype\tCount\tPercentage\n")

# 1. Load dictionary queried by haplotype into memory
reference_haplotype_names = OrderedDict({})
reference_haplotype_count = {}
haplotype_name_len = 0
for line in reference_variants.readlines():
    line = line.split()
    reference_haplotype_names[line[1]] = line[0]
    reference_haplotype_count[line[0]] = 0
    haplotype_name_len = len(line[0])
no_haplotype_name = "-"
while len(no_haplotype_name) < haplotype_name_len:
    no_haplotype_name += "-"
reference_haplotype_count[no_haplotype_name] = 0

# 2. Categorize reads
total_reads = 0
for line in read_variants.readlines():
    line = line.split()
    read_haplotype = no_haplotype_name
    if line[1] in reference_haplotype_names.keys():
        read_haplotype = reference_haplotype_names[line[1]]
    output_file.write(line[0] + '\t' + read_haplotype  + '\t' + line[1].strip() + '\n')
    reference_haplotype_count[read_haplotype] += 1
    total_reads += 1

# 3. Write stats

for haplotype in reference_haplotype_names.keys():
    count = reference_haplotype_count[reference_haplotype_names[haplotype]]
    percentage = 100 * count / total_reads
    stats_file.write(reference_haplotype_names[haplotype] + "\t" + str(count) + "\t" + str(percentage)[0:10] + '\n')
count = reference_haplotype_count[no_haplotype_name]
percentage = 100 * count / total_reads
stats_file.write(no_haplotype_name + "\t" + str(count) + "\t" + str(percentage)[0:10] + '\n')

stats_file.write("TOTAL\t" + str(total_reads) + "\t100")
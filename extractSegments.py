# extracts segments from input fasta file and blastN results

import sys
# from collections import OrderedDict

if(len(sys.argv) < 3):
    print("Usage: Python3 reads.fasta blastN.blastn")
    exit()

#keygene # python3 extractSegments.py Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta ./blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.l9k.full.Sorted

fasta_file = open(sys.argv[1], 'r')
blastn_file = open(sys.argv[2], 'r')

segment_file = open(sys.argv[1] + ".segmented", 'w')
undetected_file = open(sys.argv[1] + ".undetected", 'w')

all_file = open(sys.argv[1] + ".allregions", 'w')

raw_reads = {}
# 1. read in fasta file into dictionary in memory
print("Reading reads into memory")
while True:
    name = fasta_file.readline().split()
    seq = fasta_file.readline().strip()
    if not seq:
        break
    name = name[0][1:]
    raw_reads[name] = seq

# 2. extract segments from fasta file
print("Extracting segments")
read_ctr = 0
previous_read_id = ""
for result in blastn_file.readlines():
    result = result.split()
    read_id = result[1]
    if read_id != previous_read_id:
        read_ctr += 1
        segment_ctr = 0
        previous_end_coord = 0
    previous_read_id = read_id
    segment_ctr += 1

    start_coord = int(result[12])
    end_coord = int(result[13])

    segment_len = start_coord - previous_end_coord
    if segment_ctr > 1 and (segment_len) > 19: # undetected segment
        segment_id = read_id + "_" + str(read_ctr) + "_" + str(segment_ctr) + "\tS:" + str(start_coord) + "\tE:" + str(end_coord) + "\tL:" + str(segment_len)
        segment_seq = raw_reads[read_id][previous_end_coord : start_coord]
        
        undetected_file.write(">" + segment_id + "\n" + segment_seq + "\n")
        all_file.write(">" + segment_id + "\tC:UNDETECTED" + "\n" + segment_seq + "\n")
        segment_ctr += 1

    segment_len = end_coord - start_coord
    
    previous_end_coord = end_coord
    segment_id = read_id + "_" + str(read_ctr) + "_" + str(segment_ctr) + "\tS:" + str(start_coord) + "\tE:" + str(end_coord) + "\tL:" + str(segment_len)
    segment_seq = raw_reads[read_id][start_coord : end_coord]
    segment_file.write(">" + segment_id + "\n" + segment_seq + "\n")
    all_file.write(">" + segment_id + "\tC:RIBO" + "\n" + segment_seq + "\n")
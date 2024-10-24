import sys
import os
from collections import OrderedDict

if(len(sys.argv) < 2):
    print("Usage: Python3 findNORTelomericLinkerSeq.py reads.fasta")
    exit()

# python3 findNORTelomericLinkerSeq.py Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta.Telomeric
# python3 findNORTelomericLinkerSeq.py Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta.Telomeric

fasta_file = open(sys.argv[1], 'r')
output_fasta = open(sys.argv[1] + ".linkerseq", 'w')

total_num = 0

while True:
    name = fasta_file.readline().strip()
    seq = fasta_file.readline()
    if not seq:
        break
    last_TEL = -1
    try:
        last_TEL = seq.rindex("AAACCCTAAACCCT")
    except:
        pass
    
    first_NOR = -1
    try:
        first_NOR = seq.find("CTGGAAACGACTCAGTCGGAGGTAG")
    except:
        pass

   # assert (first_NOR - last_TEL) > 0, "not have sequences present in correct order"
    len = first_NOR - last_TEL
    if len < 7000 and last_TEL >= 0 and first_NOR >= 0 and len > 0:
        total_num +=1
        output_fasta.write(name + "_linkerseq_len:" + str(len) + "\n")
        output_fasta.write(seq[last_TEL : (first_NOR + 25)] + "\n")
    else:
        print("skipped sequence: " , name, " with length: ", len , " and telomere, NOR pos", last_TEL, first_NOR)
        
print("total sequences: ", total_num)

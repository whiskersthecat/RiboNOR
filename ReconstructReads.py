import sys

if(len(sys.argv) < 2):
    print("Usage: Python3 ReconstructReads.py Variants.[groupname].Select.Seqs.Rows.Classified.tab")
    exit()

class_file = open(sys.argv[1], 'r')
reconstructed_file = open(sys.argv[1] + ".reconstructed", 'w')

MISSING_SEG = "__"

cur_read_name = ""
prev_read_name = ""
hap_seq = ""
read_index = 0
for line in class_file.readlines():
    line = line.split()
    read = line[0].split('_')
    cur_read_name = ''.join(read[:-2])
    if cur_read_name != prev_read_name and prev_read_name != "":
        reconstructed_file.write(prev_read_name + "\t" + hap_seq + "\n")
        hap_seq = ""
        read_index = 0
    read_index += 1
    while read_index < int(read[-1]):
        read_index += 1
        hap_seq += MISSING_SEG + " "
    
    hap_seq += line[-2] + " "
    prev_read_name = cur_read_name

reconstructed_file.write(prev_read_name + "\t" + hap_seq + "\n")
    


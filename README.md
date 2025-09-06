# Note: the work here is being deprecated. For a contemporary repository (though not up to date), please see: https://github.com/whiskersthecat/tandem-repeat-assembly-workflow/tree/main



# PART A: Generate haplotypes using a large threshold for SNPs (2%)

### 1. align hifi reads to reference haplotype
../Recombination_Rate/minimap2 -t 20 --eqx -a ./reference/Lsat_NOR_Reference_13864.fasta ./reads/test.v3.hifi.fasta > ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam

### 2. variant calling with my program
python3 hic_variant_caller/hicvariantcaller.py ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam ./reference/Lsat_NOR_Reference_13864.fasta 0.02

### 3. run findvariants2
python3 ./find_variants/findvariants2.py ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam.VAR.0.02.tab ./reference/Lsat_NOR_Reference_13864.fasta

### analyze the bigtable file
cd alignments
### 4. cut and Transpose
bash x-cut-transpose-x.sh test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam.AllVariants.BIGtable.tab hifi4
### 5. 
bash x-sort-group-x.sh hifi4

### 6. add names to higher ones that have at least 14 reads
awk '$1 >= 14 {printf "NOR_HAPLOTYPE_%03dC%04d\t%s\n", NR, $1, $2 }' Variants.hifi4.Select.Seqs.Uniq.Count > Variants.hifi4.Select.Seqs.Uniq.Count.Named

### 7. convert to fasta
awk '{print ">" $1 "\n" $2}' Variants.hifi4.Select.Seqs.Uniq.Count.Named > Variants.hifi4.Select.Seqs.Uniq.Count.fasta

### 8. run the phylogenetic tree analysis
// input fasta to simple phylogeny
// download tree file and screenshot output
python3 phylogeny_id.py hifi3/tree.tree

### 9. sort the renamed tree
sort hifi4/tree.tree.RENAMED.tsv > hifi4/tree.tree.RENAMED.sorted.tsv

### 10. join the renamed table and the .Named table from step 2
join ./phylogeny/hifi4/tree.tree.RENAMED.sorted.tsv ./alignments/Variants.hifi4.Select.Seqs.Uniq.Count.Named | awk '{print $2"\t"$4}' > ./classify_variants/Haplotypes.hifi4.30.tab

### 10.5 generate the wildcard haplotypes 1X and 2X (then rename and move this file)
python3 find_variants/generateWildcardHaplotypes.py alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam.VAR.0.02.tab

### 11. TEST classify the haplotypes for the hifi4 reads
python3 NORclassifier2.py ../alignments/Variants.hifi4.Select.Seqs.Rows Haplotypes.hifi4.30.tab Haplotypes.hifi4.Wildcard.tab
###### the counts check out with the counts in the original names!

### 12. Reconstruct the original reads
python3 ReconstructReads.py Variants.hifi4.Select.Seqs.Rows.Classified.tab

# Generate haplotypes for indels using a large threshold

### 2.0 change hard-coded parameters for variant caller (small indels > 0.04 or larger indels > 0.01)
### 2. variant calling with my program
python3 hic_variant_caller/hicvariantcaller.py ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam ./reference/Lsat_NOR_Reference_13864.fasta 0.04
... rename the file to .INDEL

### 3. run findvariants2
python3 ./find_variants/findvariants2.py ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam.VAR.INDEL.0.04.tab ./reference/Lsat_NOR_Reference_13864.fasta

### analyze the bigtable file
cd alignments
### 4. cut and Transpose
bash x-cut-transpose-x.sh test.v3.hifi.mapping_to_Lsat_NOR_Reference_13864.sam.AllVariants.BIGtable.tab hifiInDel4
### 5. 
bash x-sort-group-x.sh hifiInDel4

### 6. add names to the first 99 (all have 16 or more reads)
awk 'NR == 1, NR == 99 {printf "%02d\t%s\n", NR, $2 }' Variants.hifiInDel4.Select.Seqs.Uniq.Count > Haplotypes.INDEL.hifi4.99.tab
... moved it to the classify variants directory

### 7. TEST classify the haplotypes for the hifi4 reads
python3 NORclassifier2.py ../alignments/Variants.hifiInDel4.Select.Seqs.Rows Haplotypes.INDEL.hifi4.99.tab Haplotypes.INDEL.hifi4.99.tab

### 8. Reconstruct the original reads
python3 ReconstructReads.py Variants.hifiInDel4.Select.Seqs.Rows.Classified.tab


# PART B CLASSIFY ALL READS
### -1. create segments from raw reads
... in directory /share/rwmstore/Assemblies_RWM/Lettuce/Lsat_V14/Ribo_NOR_Scan/BlastN_Filter
... BLAST information
cat X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.l9k | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.l9k.full.Sorted
cat X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_ULR_2X.l9k | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_ULR_2X.l9k.full.Sorted
cat X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_UC_HiFi.l9k | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_UC_HiFi.l9k.full.Sorted
cat X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_HiFi_1X_2024_Select.l9k | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_HiFi_1X_2024_Select.l9k.full.Sorted

awk -F'\t' '$6 >= 6000' Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.m7.e240 | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.l6k.full.Sorted
awk -F'\t' '$6 >= 6000' Lsat_NOR_Reference_13864.vs.Salinas_ULR_2X_2024_08.Select.Ribo_NOR.m7.e240 | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_ULR_2X.l6k.full.Sorted
awk -F'\t' '$6 >= 6000' Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_UC_HiFi_UniDir.m7.e240 | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_UC_HiFi.l6k.full.Sorted
awk -F'\t' '$6 >= 6000' Lsat_NOR_Reference_13864.vs.Salinas_HiFi_1X_2024_Select_Ribo_NOR.m7.e240 | sort -k17nr -k2,2 -k13,13n > /share/rwmwork/peter/Ribo_Variants/reads/blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_HiFi_1X_2024_Select.l6k.full.Sorted


python3 extractSegments.py ./full_reads/Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta ./blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_ULR_2X.l6k.full.Sorted
python3 extractSegments.py ./full_reads/Salinas_HiFi_1X_2024_Select_Ribo_NOR.XXX.Fasta ./blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Salinas_HiFi_1X_2024_Select.l6k.full.Sorted
python3 extractSegments.py ./full_reads/Ribo_NOR_Reads_UC_HiFi_UniDir.2024.06.12.Sorted.Fasta ./blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_UC_HiFi.l6k.full.Sorted
python3 extractSegments.py ./full_reads/Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta ./blastN/X-Filter-Lsat_NOR_Reference_13864.vs.Ribo_NOR_Reads_KG_ONT_UniDir.l6k.full.Sorted


... manually move those files to segmented_reads
... I will try to analyze all of the segments even if they are not NOR
... they will be classified as -- if any of the following are true:
...  1. the segment does not align

### 0a. merge ONT, Keygene, HiFi1, and HiFi2 Segments
cat Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta.segmented Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta.segmented Ribo_NOR_Reads_UC_HiFi_UniDir.2024.06.12.Sorted.Fasta.segmented Salinas_HiFi_1X_2024_Select_Ribo_NOR.XXX.Fasta.segmented > ONT.KG.HiFi1.HiFi2.Fasta.segmented
cat Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta.undetected Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta.undetected Ribo_NOR_Reads_UC_HiFi_UniDir.2024.06.12.Sorted.Fasta.undetected Salinas_HiFi_1X_2024_Select_Ribo_NOR.XXX.Fasta.undetected > ONT.KG.HiFi1.HiFi2.Fasta.undetected
cat Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta.allregions Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta.allregions Ribo_NOR_Reads_UC_HiFi_UniDir.2024.06.12.Sorted.Fasta.allregions Salinas_HiFi_1X_2024_Select_Ribo_NOR.XXX.Fasta.allregions > ONT.KG.HiFi1.HiFi2.Fasta.allregions

... 233016 total lines
### 0b. merge ONT, Keygene, Hifi1, and Hifi2 Full Reads
cat Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta Ribo_NOR_Reads_UC_HiFi_UniDir.2024.06.12.Sorted.Fasta Salinas_HiFi_1X_2024_Select_Ribo_NOR.XXX.Fasta > ONT.KG.HiFi1.HiFi2.Fasta
... 112014 total lines

### 1a. align new segments to reference
../Recombination_Rate/minimap2 -t 20 --eqx -a ./reference/Lsat_NOR_Reference_9847.fasta ./reads/segmented_reads/ONT.KG.HiFi1.HiFi2.Fasta.segmented > ./alignments/ONT.KG.HiFi1.HiFi2.segmented.sam
... 195.092 seconds real time
../Recombination_Rate/minimap2 -t 20 --eqx -a ./reference/Lsat_NOR_Reference_9847.fasta ./reads/segmented_reads/other/ONT.KG.HiFi1.HiFi2.Fasta.allregions > ./alignments/ONT.KG.HiFi1.HiFi2.allregions.sam


### 1b. align Full Reads to reference duplicated
awk 'NR == 1 {print $0} NR == 2 {for(i = 0; i < 40; i++) {printf($0)};}' Lsat_NOR_Reference_9847.fasta > Lsat_NOR_Reference_9847x40.fasta
../Recombination_Rate/minimap2 -t 20 --eqx -a ./reference/Lsat_NOR_Reference_9847x40.fasta ./reads/full_reads/ONT.KG.HiFi1.HiFi2.Fasta > ./alignments/ONT.KG.HiFi1.HiFi2.sam


### 2a. find SNV variants in segments
python3 ./find_variants/findvariants2.py ./alignments/ONT.KG.HiFi1.HiFi2.allregions.sam ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_9847.sam.VAR.snvmnv.0.02.tab ./reference/Lsat_NOR_Reference_9847.fasta
... [INFO] 151167 total reads
... rename this bigtable file and log file with SNV.MNV

### 2b. find INDEL variants in segments
python3 ./find_variants/findvariants2.py ./alignments/ONT.KG.HiFi1.HiFi2.allregions.sam ./alignments/test.v3.hifi.mapping_to_Lsat_NOR_Reference_9847.sam.VAR.insdel.0.04.tab ./reference/Lsat_NOR_Reference_9847.fasta
... rename this bigtable file and log file with INS.DEL

cd alignments
### 3. cut and Transpose
... trial 1 is with 9k plus BLAST filtering
... trial 2 is with 6k plus BLAST filtering
... trial 3 is with 6k plus BLAST filtering, working with all segments (including undetected ones)
... trial 4 is with 6k plus BLAST filtering, working with all segments (including undetected ones), fixed findvariants adding N tokens
bash x-cut-transpose-x.sh ONT.KG.HiFi1.HiFi2.allregions.sam.SNV.MNV.BIGtable.tab ont.kg.hf1.hf2.4.snvmnv
bash x-cut-transpose-x.sh ONT.KG.HiFi1.HiFi2.allregions.sam.INS.DEL.BIGtable.tab ont.kg.hf1.hf2.4.insdel

### 4. classify the snvmnv haplotypes and insdel haplotypes
cd ../classify_variants
python3 NORclassifier3.py ../alignments/Variants.ont.kg.hf1.hf2.4.snvmnv.Select.Seqs.Rows Haplotypes.hifi4.30.tab Haplotypes.hifi4.Wildcard.tab
python3 NORclassifier3.py ../alignments/Variants.ont.kg.hf1.hf2.4.insdel.Select.Seqs.Rows Haplotypes.INDEL.hifi4.99.tab Haplotypes.INDEL.hifi4.99.tab

... -- = segment did not match haplotype
... __ = segment did not align
### 5. reconstruct the original reads
python3 ReconstructReads.py Variants.ont.kg.hf1.hf2.4.snvmnv.Select.Seqs.Rows.Classified.tab
... 55829 reads
python3 ReconstructReads.py Variants.ont.kg.hf1.hf2.4.insdel.Select.Seqs.Rows.Classified.tab
... 55829 reads

### 6. generate and merge the tables of coverage information
... Ribo Telo Copia Gypsy Chlor NBS
paste <(cut -f 1-3 <(sort LST_01_Ribo.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort LST_01_Ribo.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort LST_00_Telo.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort RVT_Copia.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort RVT_Gypsy.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort Lsat_Chlp.vs.XRibo_Sal_ONT_2X.asm.coverage)) <(cut -f 7 <(sort NBS_Cons.vs.XRibo_Sal_ONT_2X.asm.coverage)) | sort -k4n >  /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_ONT_2X.asm.coverage
paste <(cut -f 1-3 <(sort LST_01_Ribo.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort LST_01_Ribo.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort LST_00_Telo.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort RVT_Copia.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort RVT_Gypsy.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort Lsat_Chlp.vs.XRibo_Sal_KG_ONT.asm.coverage)) <(cut -f 7 <(sort NBS_Cons.vs.XRibo_Sal_KG_ONT.asm.coverage)) | sort -k4n >  /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_KG_ONT.asm.coverage
paste <(cut -f 1-3 <(sort LST_01_Ribo.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort LST_01_Ribo.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort LST_00_Telo.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort RVT_Copia.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort RVT_Gypsy.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort Lsat_Chlp.vs.XRibo_Sal_HiFi23.asm.coverage)) <(cut -f 7 <(sort NBS_Cons.vs.XRibo_Sal_HiFi23.asm.coverage)) | sort -k4n >  /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi23.asm.coverage
paste <(cut -f 1-3 <(sort LST_01_Ribo.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort LST_01_Ribo.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort LST_00_Telo.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort RVT_Copia.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort RVT_Gypsy.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort Lsat_Chlp.vs.XRibo_Sal_HiFi24.asm.coverage)) <(cut -f 7 <(sort NBS_Cons.vs.XRibo_Sal_HiFi24.asm.coverage)) | sort -k4n >  /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi24.asm.coverage

awk '{printf "%s%s%f\n", $0, "\t", 1 - $9 - $8 - $7 - $6 - $5 - $4}' /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_ONT_2X.asm.coverage | sort -k10nr > /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_ONT_2X.asm.coverage.sum
awk '{printf "%s%s%f\n", $0, "\t", 1 - $9 - $8 - $7 - $6 - $5 - $4}' /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_KG_ONT.asm.coverage | sort -k10nr > /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_KG_ONT.asm.coverage.sum
awk '{printf "%s%s%f\n", $0, "\t", 1 - $9 - $8 - $7 - $6 - $5 - $4}' /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi23.asm.coverage | sort -k10nr | sed 's/_//g' > /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi23.asm.coverage.sum
awk '{printf "%s%s%f\n", $0, "\t", 1 - $9 - $8 - $7 - $6 - $5 - $4}' /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi24.asm.coverage | sort -k10nr | sed 's/_//g' > /share/rwmwork/peter/Ribo_Variants/classify_variants/AllComponents.vs.XRibo_Sal_HiFi24.asm.coverage.sum

... Telo Copia Gypsy Chlor NBS other
cat AllComponents.vs.XRibo_Sal_KG_ONT.asm.coverage.sum AllComponents.vs.XRibo_Sal_ONT_2X.asm.coverage.sum AllComponents.vs.XRibo_Sal_HiFi23.asm.coverage.sum AllComponents.vs.XRibo_Sal_HiFi24.asm.coverage.sum | awk '{printf ("%-36s\t%d", $1, $3); for(i = 5; i < 11; i++) { if ($i != 0) {printf("\t%03.1i", $i * 1000)} else {printf("\t---")}}; print("")}' | sort > AllComponents.vs.ONT.KG.HiFi1.HiFi2.coverage.sum
... 56007 lines

### 7. join the reconstructed tables with the coverages
... first sort the reconstructed files
sort Variants.ont.kg.hf1.hf2.4.snvmnv.Select.Seqs.Rows.Classified.tab.reconstructed | awk -F'\t' {'printf("\t%-36s\t%s\n", $1, $2)'} > Variants.ont.kg.hf1.hf2.4.snvmnv.Select.Seqs.Rows.Classified.tab.reconstructed.sorted
sort Variants.ont.kg.hf1.hf2.4.insdel.Select.Seqs.Rows.Classified.tab.reconstructed | awk -F'\t' {'printf("%-36s\t%s\n", $1, $2)'}    > Variants.ont.kg.hf1.hf2.4.insdel.Select.Seqs.Rows.Classified.tab.reconstructed.sorted
... note: the components file has more reads than the classified file

... join the three tables
... join the first two
join -t $'\t' -2 2 AllComponents.vs.ONT.KG.HiFi1.HiFi2.coverage.sum Variants.ont.kg.hf1.hf2.4.snvmnv.Select.Seqs.Rows.Classified.tab.reconstructed.sorted > Variants.ont.kg.hf1.hf2.4.snvmnv.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted
join -t $'\t' Variants.ont.kg.hf1.hf2.4.snvmnv.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted Variants.ont.kg.hf1.hf2.4.insdel.Select.Seqs.Rows.Classified.tab.reconstructed.sorted > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted

### 7c. rename the keygene reads
awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {  ids[(substr($1, 2, 36))]; next } { if ($1 in ids) $1 = $1"\tK"; print }' ./reads/full_reads/Ribo_NOR_Reads_KG_ONT_UniDir.2024.06.12.Sorted.Fasta ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted > ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed
awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {  ids[(substr($1, 2, 36))]; next } { if ($1 in ids) $1 = $1"\tO"; print }' ./reads/full_reads/Salinas_ULR_2X_2024_08.Select.Ribo_NOR.Three_Units_Plus.UniDir.XXX.Fasta ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed > ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed2


awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {  ids[(substr($1, 2, 16))]; next } { if ((substr($1, 2, 16) in ids)) $1 = $1"\th"; print }' ./classify_variants/AllComponents.vs.XRibo_Sal_HiFi23.asm.coverage.sum ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed2 > ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed3
awk -F'\t' 'BEGIN {OFS = "\t"} NR == FNR {  ids[(substr($1, 2, 16))]; next } { if ((substr($1, 2, 16) in ids)) $1 = $1"\tH"; print }' ./classify_variants/AllComponents.vs.XRibo_Sal_HiFi24.asm.coverage.sum ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed3 > ./classify_variants/Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4


### 7c. sort by the presence of telomeric sequence (field 3), if tied, length of read
sort -k4nr -k3nr Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted

### 7d. reformat the table so that the haplotype sequences have a cosntant distance seperation, and the length has a constant distance seperation
awk -F'\t' '{printf("%s\t%s\t%-6s\t", $1, $2, $3); for(i = 4; i < 11; i++) {printf("%s\t", $i)}; printf("%-78s\t%s\n", $11, $12)}' Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat

... set the tab size to 3 in settings!!!

### 7e. add the linkerinfo
python3 appendLinkerInfo.py reads_TATCAATTTAAAATTATTTATGAAACTTGTACAGCAAAA.fasta.sorted reads_TATCAATTTAAAAATATTTATGAAACTTGTACAGCAAAA.fasta.sorted Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat

### 8. split the reads by what groups of haplotypes they contain
python3 sortByGroup.py Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final
group 1 only
group 2 only
group 3 only
group 4 only
group 5 only
groups 1 and 2
groups 1 and 3
groups 1 and 4
groups 1 and 5
groups 2 and 3
groups 2 and 4
groups 2 and 5
groups 3 and 4
groups 3 and 5
groups 4 and 5
three or more total groups
no groups

### 8.b make subgroups
... for groups 1 and 2 reads, there are a few different types
... d1: take the ones that have 2. 1. 2. (1m haplotype usually)
grep -E "2. 1. 2. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.d1
grep -vE "2. 1. 2. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd1

... d2: the 1 to 2 transition, take the ones that have 1. 2. (and are not in d1)
grep -E "1. 2. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd1 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.d2
grep -vE "1. 2. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd1 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd2

... d3: the 2 to 1 transition, take the ones that have 2. 1. (and are not in d1 nor d2)
grep -E "2. 1. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd2 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.d3
grep -vE "2. 1. " Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd2 > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.group1.group2.nd3

### 8.c extract the genic ends
grep "L1" Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.linker1
grep "L2" Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final > Variants.ont.kg.hf1.hf2.4.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.sorted.renamed4.sorted.reformat.final.linker2



### 9. duplicate the files, then manually annotate with a:11 where a is the chromosome and 11 is the repeat position (starting with 0) (c belongs to both chromosomes)
... duplicated 4 only

### 10. concatenate back together the partitioned reads
cat <ls -1 | grep '.*group.*'> > Variants.ont.kg.hf1.hf2.1.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.labeled

### 11. generate two sam files from the sam file in step 0b:
... place the left alignment coordinate to 9487 * n of each read based on the repeat position assigned in step 9
... generate two sam files for both chromosomes (c labeled will be in both files)
|| TO COMPLETE
python3 NORconstructor.py Variants.ont.kg.hf1.hf2.1.snvmnv.insdel.coverage.Select.Seqs.Rows.Classified.tab.reconstructed.ALL.labeled

### 12. visualize on CLC with reference as repeated reference repeat, collapse the alignments into one contig for both choromosomes


### ... QC, realign the reads to this contig and check for constant coverage

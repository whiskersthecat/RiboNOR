if [ "$#" -ne 2 ]; then
    echo "Usage: bash x-cut-transpose-x.sh bigtablefile.tab nameofreadgroup"
    exit 2
fi
echo "cutting big table"
cut -f 1-7 --complement $1 > Variants.$2.Select.Seqs.Columns
echo "transposing table"
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]"\t"
        for(i=2; i<=NR; i++){
            str=str""a[i,j];
        }
        print str
    }
}' Variants.$2.Select.Seqs.Columns > Variants.$2.Select.Seqs.Rows

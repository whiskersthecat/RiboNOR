if [ "$#" -ne 1 ]; then
    echo "Usage: bash x-sort-group-x.sh nameofreadgroup"
    exit 2
fi
echo "sorting"
sort -k2 Variants.$1.Select.Seqs.Rows > Variants.$1.Select.Seqs.Sorted
echo "grouping and counting"
cut -f 2 Variants.$1.Select.Seqs.Sorted | sort | uniq -c | sort -k1,1nr > Variants.$1.Select.Seqs.Uniq.Count
#!/bin/bash

# Get function arguments
countsDIR=$1        # counts_dir: 

echo "Sort the count files by first column ..."
# sort the count files by first column
cd ${countsDIR}
mkdir unsorted && mv *.txt unsorted
for f in unsorted/*.txt
do
 of=$(basename $f)
 # skip the footer which starts with "__"
 #grep -vE "__no_feature|__ambiguous|__too_low_Qual|__not_aligned|__alignment_not_unique" ${f} | sort -k1,1 -u > ${of}
 grep -vE "__no_feature|__ambiguous|__too_low_Qual|__not_aligned|__alignment_not_unique" ${f} | sort -k1,1 -u > ${of}
 grep -E  "__no_feature|__ambiguous|__too_low_Qual|__not_aligned|__alignment_not_unique" ${f} >> ${of}
done
rm -rf unsorted
cd -

echo ""
echo "Check if all the files has same number of mirs"
# check if all the files has same number of mirs
find ${countsDIR} -name *_*.txt -type f  -exec wc -l {} \;
mkdir test_counts
for f in ${countsDIR}/*.txt
do
 b=$(basename ${f} .txt)
 cut -f1 ${f} > test_counts/${b}.txt
done
a1=`ls test_counts | head -1`
for f in test_counts/*.txt
do
 diff test_counts/${a1} $f
done
rm -rf test_counts


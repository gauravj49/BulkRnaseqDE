
# Setup relevant environment/system variables
Log2FC=0.25
padj=0.05
baseMean=1
rnaClass="gene"
sort=0
species="hg" # hg or "mm"
normCounts="TRUE"
filterGenes="TRUE"
effectsToControl=""
sampleNameColumn=""
saveRdata="FALSE"

# Project Specific variables
projprefix="user_${rnaClass}_projectname"
projdirprefix="user/rnaseq/projectname"

# Get DEseq2 input
jobdir="/path/to/bkRnaseqDE"
projdir="${jobdir}/output/${projdirprefix}"
DEseq2OT="${projdir}/01.1_${projprefix}_${rnaClass}_bm${baseMean}"
DEseq2IN="${DEseq2OT}/deseq2IN"
countsDIR="${jobdir}/output/${projdirprefix}/counts"
qlogdir="${DEseq2OT}/logs/bsub"
mkdir -p ${DEseq2IN} ${DEseq2OT} ${qlogdir}

# sort the count files by first column
if [ ${sort} -ne "0" ]; then
 bash scripts/sort_counts_dir.sh ${countsDIR}
fi

# Get the DEseq2IN files
cd ${countsDIR}
ls *_A_* > ${DEseq2IN}/DisTrt.txt
ls *_B_* > ${DEseq2IN}/DisUntrt.txt
ls *_C_* > ${DEseq2IN}/CltTrt.txt
ls *_D_* > ${DEseq2IN}/CltUntrt.txt
cd - 

# Check the number of lines in the DEseq2IN files
wc -l ${DEseq2IN}/*.txt

# Define comparisons
t=(DisTrt CltTrt)
c=(DisUntrt CltUntrt)

# Run RUVseq-DEseq2
 for ((i=0;i<${#t[@]};i++))
 do
    echo "-------------------------------------------------------------------"
    echo "sh scripts/rnaseq_run_deseq2_RUVseq.sh ${countsDIR} ${DEseq2IN}/${c[$i]}.txt ${DEseq2IN}/${t[$i]}.txt "${c[$i]},${t[$i]}"  ${DEseq2OT}/${t[$i]}_over_${c[$i]} ${Log2FC} ${padj} ${baseMean} -xc '' ${species} FALSE ${normCounts} ${filterGenes} ${saveRdata}"
    echo "-------------------------------------------------------------------"
    echo ""
    bsub -J ${projprefix}_${rnaClass}_${t[$i]}_over_${c[$i]} -o ${qlogdir}/${t[$i]}_over_${c[$i]}.out -e ${qlogdir}/${t[$i]}_over_${c[$i]}.err sh scripts/rnaseq_run_deseq2_RUVseq.sh ${countsDIR} ${DEseq2IN}/${c[$i]}.txt ${DEseq2IN}/${t[$i]}.txt "${c[$i]},${t[$i]}"  ${DEseq2OT}/${t[$i]}_over_${c[$i]} ${Log2FC} ${padj} ${baseMean} "-xc" "" ${species} "FALSE" ${normCounts} ${filterGenes} ${saveRdata}
 done




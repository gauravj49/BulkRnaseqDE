#!/bin/bash

# For DEseq2 analysis
# Get function arguments
ID=$1             # raw_counts_dir:
CL=$2             # Control sample files list
TR=$3             # Treatment sample files list
CV=$4             # Condition values <control,treatment>: "WTKO"
OF=$5             # Output file
LF=$6             # Log2 Fold Change treshold
PD=$7             # Padjusted threshold
BM=$8             # Basemean threshold
XC=${9:-""}       # Xaxis clustering flag for sample names
YC=${10:-""}      # Yaxis clustering flag for feature names
SP=${11:-"mm"}    # hg or mm ....Species for annotation file
SK=${12:-"FALSE"} # Use ERCC controls for normalization
NR=${13:-"FALSE"} # If set, also get normalized reads
FG=${14:-"TRUE"}  # If set FALSE, do not filter genes/smallRNAs during RUVseq
SV=${15:-"FALSE"} # If set, save all data in the logs folder. This can restore all variables for debugging purpose later


# Create a log file
mkdir -p $(dirname ${OF})"/logs"
LOG=$(dirname ${OF})"/logs/"$(basename ${OF} .txt).log

# Get base directories
d=$(dirname ${CL} )
c=$(basename ${CL} .txt)
t=$(basename ${TR} .txt)

# Run DEseq2 analysis
echo "-------------------------------------------------------------------"
echo "Rscript scripts/R_RUVseq_DESeq2_DEanalysis.R ${ID} ${CL} ${TR} ${CV} ${OF} ${SK} gene ${NR} ${FG} ${SV} 2>&1 | tee -a ${LOG}"
echo "-------------------------------------------------------------------"
Rscript scripts/R_RUVseq_DESeq2_DEanalysis.R ${ID} ${CL} ${TR} ${CV} ${OF} ${SK} gene ${NR} ${FG} ${SV} 2>&1 | tee -a ${LOG}

# Draw differential expression heatmaps
RF=$(dirname ${OF})"/results_DEseq2/"$(basename ${OF} .txt)_DE_RESULTS.txt   # Results file
HM=$(dirname ${OF})"/heatmaps_DEseq2/"$(basename ${OF} .txt)_DE_heatmaps.txt # Heatmap file

# Add gene annotation to results file
echo "-------------------------------------------------------------------"
echo "python scripts/add_gene_annotation2deseq2_results.py -af=input/annotation/bed/${SP}_geneSymbols.bed -if=${RF} 2>&1 | tee -a ${LOG}"
echo "-------------------------------------------------------------------"
python scripts/add_gene_annotation2deseq2_results.py -af=input/annotation/bed/${SP}_geneSymbols.bed -if=${RF} 2>&1 | tee -a ${LOG}

# Draw the heatmaps
GRF=$(dirname ${OF})"/results_DEseq2/"$(basename ${OF} .txt)_DE_RESULTS_geneSymbols.txt    # Results file with common gene names
echo "-------------------------------------------------------------------"
echo "python scripts/draw_mirs_heatmap.py -if=${GRF} -of=${HM} -lf=${LF} -pj=${PD} -bm=${BM} ${XC} ${YC} -rq 2>&1 | tee -a ${LOG}"
echo "-------------------------------------------------------------------"
python scripts/draw_mirs_heatmap.py -if=${GRF} -of=${HM} -lf=${LF} -pj=${PD} -bm=${BM} ${XC} ${YC} -rq 2>&1 | tee -a ${LOG}
mv $(dirname ${OF})"/heatmaps_DEseq2/"*.log $(dirname ${OF})"/logs/"

# Add gene annotation to filtered heatmap file
HMF=$(dirname ${HM})/$(basename ${HM} .txt)"_differentially_expressed.txt"
# Create summary files in excel
SRD=$(dirname ${OF})"/summary_excel"
mkdir -p ${SRD}

# Get filenames
SRF=${SRD}/$(basename ${RF} .txt)".xls"
SMF=${SRD}/$(basename ${HM} .txt)"_differentially_expressed.xls"

# Sort the filtered files on padj (smallest to largest)
head -1 ${RF}  > tmp | tail -n+2 ${RF}  | sort -k9,9n >> tmp && mv tmp ${RF}
head -1 ${HMF} > tmp | tail -n+2 ${HMF} | sort -k9,9n >> tmp && mv tmp ${HMF}

# Convert them to excel
ssconvert ${RF} ${SRF} 
ssconvert ${HMF} ${SMF}

# Move excel files to summary folder
mv $(dirname ${OF})"/heatmaps_DEseq2/"*.xls ${SRD}
mv $(dirname ${OF})"/results_DEseq2/"*.xls ${SRD}

# Move all the files without common gene names to ohneGeneNames folder
OGND=${SRD}/ohneGeneNames
mkdir -p ${OGND}
mv $(find ${SRD} \! -name "*geneSymbols*" -type f) ${OGND}

##### NOTES ####
# 1) on appending with tee
# source: https://askubuntu.com/questions/808539/how-to-append-tee-to-a-file-in-bash
## echo -e "First Line"  | tee ~/output.log
## echo -e "Second Line" | tee -a ~/output.logs
# -a modifier is for 'append', or add to the end. Without -a, the tee command overwrites the file

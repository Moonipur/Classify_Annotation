#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
PACKAGE="/mnt/hdd/public/share_resource/OS_WES"
cur_dir=`pwd`

Help() {
    echo "usage: CNV_annotation [-i|-o|-r|h]"
    echo "description: CNV_annotation is a software suite script that uses the ACE and ClassifyCNV packages"
    echo "             for annotating gains and losses of CNVs and predicting list of dosage-sensitive genes"
    echo "             that affect as pathogenic variant."
    echo "optional argruments:"
    echo "    -i       Input file path (BAM file)"
    echo "    -o       Output directory path (/path/of/output/directory/), Default is the currect directory"
    echo "    -r       Type of reference human genome (hg19/hg38)"
    echo "    -h       Show this help message and exit"
    echo
    echo "author's email:"
    echo "    songphon_sutthittha@cmu.ac.th"
    echo
    echo "** Please contact us if you have any questions or problems with this script."
    echo "------------------------------------------------------------------------------------------"
}

while getopts ":hi:o:r:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # default [Current directory]
         name=${OPTARG};;
      o) # default [Current directory]
         outdir=${OPTARG};;
      r) # default [Current directory]
         ref=${OPTARG};;
     \?) # Invalid option
         echo "Error: Unrecognized arguments"
         exit;;
   esac
done

NAME=`echo ${name}| cut -d. -f1`

if [ ! -z ${outdir} ]
then
    OUTDIR=${outdir}
else
    OUTDIR=`pwd`
fi

#Step 1
echo -e ">> 1ST: The ACE scipt is generating and reading using the ACE software"
conda activate Rdevtools
if [ "${ref}" = "hg19" ]
then
    echo "#!/usr/bin/env Rscript" > ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "library(QDNAseq)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "library(ACE)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "bins <- getBinAnnotations(binSize=100, genome='hg19')" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "readCounts <- binReadCounts(bins)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "readCountsFiltered <- estimateCorrection(readCountsFiltered)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "copyNumbers <- correctBins(readCountsFiltered)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "copyNumbersNormalized <- normalizeBins(copyNumbers)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun='sqrt')" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "object <- copyNumbersSegmented" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "template <- objectsampletotemplate(object,index=1)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "segmentdf <- getadjustedsegments(template, cellularity = 0.25)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "sqmodel <- squaremodel(template, prows = 150, ptop = 3.3, pbottom = 1.8, " >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "                       penalty = 0.5, penploidy = 0.5)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "models <- data.frame(segmentdf)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "write.table(models, file.path('${NAME}_segmentation.tsv'), quote = FALSE, " >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo "            sep = '\t', na = '', row.names = FALSE)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
    Rscript ${OUTDIR}/ACE_SCRIPT_hg19.R
elif [ "${ref}" = "hg38" ]
then
    echo "#!/usr/bin/env Rscript" > ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "library(QDNAseq)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "library(QDNAseq.hg38)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "library(ACE)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "bins <- getBinAnnotations(binSize=100, genome='hg38')" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "readCounts <- binReadCounts(bins)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "readCountsFiltered <- estimateCorrection(readCountsFiltered)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "copyNumbers <- correctBins(readCountsFiltered)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "copyNumbersNormalized <- normalizeBins(copyNumbers)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun='sqrt')" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "object <- copyNumbersSegmented" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "template <- objectsampletotemplate(object,index=1)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "segmentdf <- getadjustedsegments(template, cellularity = 0.25)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "sqmodel <- squaremodel(template, prows = 150, ptop = 3.3, pbottom = 1.8, " >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "                       penalty = 0.5, penploidy = 0.5)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "models <- data.frame(segmentdf)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "write.table(models, file.path('${NAME}_segmentation.tsv'), quote = FALSE, " >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo "            sep = '\t', na = '', row.names = FALSE)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    echo -e "" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
    Rscript ${OUTDIR}/ACE_SCRIPT_hg38.R
fi
conda deactivate

#Step 2
echo -e ">> 2ND: The BED file is converting from TSV file using Python module"
conda activate Bedtools
cut -f1,2,3,8 ${NAME}_segmentation.tsv | sed -z 's/\t/,/g' > ${NAME}_segmentation_cut.csv
python3 ${PACKAGE}/Package_Software/Classify_Annotation/TSV2BED.py ${NAME}_segmentation_cut.tsv ${OUTDIR} ${NAME}_segmentation.csv
sed -z 's/,/\t/g' > ${NAME}_segmentation.bed
rm ${NAME}_segmentation.csv

#Step 3
echo -e ">> 3RD: The list of dosage-sensitive genes and pathogenic variants are predicting using ClassifyCNV package"
python3 ${PACKAGE}/ClassifyCNV/ClassifyCNV.py --infile ${OUTDIR}/${NAME}_segmentation.bed --GenomeBuild ${ref} --outdir ${OUTDIR}
conda deactivate

echo -e ">> 4TH: All steps already finish and saved out at; ${OUTDIR}"

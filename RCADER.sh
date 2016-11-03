#!/bin/bash

echo ""
echo ""
echo "********** RCADER **********"
echo ""
echo ""

####################### define executables
FASTAtoRF="./bin/FASTAtoRF"
rndForest="./src/_R/_predict.RF.R"
RCADER="./bin/RCADER"

####################### identify the input arguments
jobid=$1
proteins=$2
peaks=$3
scores=$4

if [ "$jobid" = "" ]; then
	echo -e "\nUsage: bash RCADER.sh <jobID> <C2H2_ZFP.fasta> <ChIP_seq.fasta> <ChIP_seq_scores.txt>\n"
	exit
fi

echo "Job ID: "$jobid
echo "Input FASTA file for the target protein(s): "$proteins
echo "Input FASTA file for the peaks: "$peaks

if [ -e "$proteins" ]; then
	echo "Protein sequence file found."
else
	echo "ERROR: Protein sequence file was not found."
	exit
fi

if [ -e "$peaks" ]; then
	echo "Peak sequence file found."
else
	echo "ERROR: Peak sequence file not found."
	exit
fi

if [ -e "$scores" ]; then
	echo "Peak score file found."
else
	echo "ERROR: Peak score file not found."
	exit
fi


####################### define temporary path
tmp_folder="./tmp/"$jobid
mkdir -p $tmp_folder
RF_in=$tmp_folder"/_predict.in"
RF_out=$tmp_folder"/_predict.RF.out"

####################### define the output path
out_folder="./out/"$jobid
mkdir -p $out_folder
rm -f $out_folder/log.step1.txt
rm -f $out_folder/log.step2.txt
out_file=$out_folder"/results"

####################### convert the input FASTA file to a covariate matrix file for the RF script
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
	$FASTAtoRF -minl 2 -maxl 8 -span $i -fasta $proteins -out $RF_in.span$i >>$out_folder/log.step1.txt
	if [ $i == 2 ]; then
		cat $RF_in.span$i > $RF_in
	else
		cat $RF_in.span$i | sed 1d >> $RF_in
	fi
done

####################### run the RF script, and reformat it for the next step
echo ""
echo ""
echo "Generating seed motifs using recognition code ..."
Rscript $rndForest $jobid
sed 's/"//g' $RF_out > $out_file.RF_out.txt

####################### run the RCADER script
echo ""
echo ""
echo "Optimizing motifs ..."
$RCADER -rf $out_file.RF_out.txt -fasta $peaks -score $scores -out $out_file  >>$out_folder/log.step2.txt









#*****************************************************************************************
# The following lines check the input/output, and produce appropriate messages
# If no error was detected in either input or output, the info messages will be written in
# ./out/<jobID>/log.info.txt
# Otherwise, the error messages will be written in
# ./out/<jobID>/log.error.txt
#*****************************************************************************************



####################### identify the input arguments
err=""
info=""

####################### define log files
step1=$out_folder/log.step1.txt
step2=$out_folder/log.step2.txt
report=$out_folder/results.report.txt

####################### check if any of the B1H-RC motifs were optimized

optimized=`cat $report | cut -f3 | grep '1' | wc -l`

if [ "$optimized" -le 0 ]; then
	err="ERROR: None of the predicted motifs from the provided C2H2-ZF proteins were enriched in the ChIP-seq peaks.\n"
else
	info="The predicted motifs of $optimized possible C2H2-ZF arrays were enriched in the ChIP-seq peaks.\n"$info
fi


####################### check if the peak sequence file had any valid sequences

numPeaks=`cat $step2 | grep 'sequences had associated scores and were retained.' | head -n 1 | cut -d ' ' -f1`

if [ "$numPeaks" = "" ]; then
	err="ERROR: Not enough sequences with associated scores were found in the input files.\n"
elif [ "$numPeaks" = "ERROR:" ]; then
	err="ERROR: Not enough sequences with associated scores were found in the input files.\n"
elif [ "$numPeaks" -le 0 ]; then
	err="ERROR: Not enough sequences with associated scores were found in the input files.\n"
else
	info="$numPeaks sequences with associated scores were found in the input files.\n"$info
fi


####################### check if the C2H2-ZF sequences have had any ZF arrays

numArrays=`cat $step2 | grep 'motifs were read.' | head -n 1 | cut -d ' ' -f1`

if [ "$numArrays" = "" ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
elif [ "$numArrays" = "ERROR:" ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
elif [ "$numArrays" -le 0 ]; then
	err="ERROR: The input C2H2-ZF sequences must have at least two adjacent canonical C2H2-ZF domains.\n"
else
	info="$numArrays possible C2H2-ZF arrays were tested.\n"$info
fi


####################### check if the C2H2-ZF file had any valid sequences

numC2H2=`cat $step1 | grep 'sequences were read.' | head -n 1 | cut -d ' ' -f1`

if [ "$numC2H2" = "" ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
elif [ "$numC2H2" = "ERROR:" ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
elif [ "$numC2H2" -le 0 ]; then
	err="ERROR: No sequences were found in the input FASTA for C2H2-ZF proteins. Please check the input format.\n"
else
	info="$numC2H2 sequences were found in the input FASTA for C2H2-ZF proteins.\n"$info
fi


####################### write the appropriate messages to the output

if [ "$err" = "" ]; then
	echo -e -n $info > $out_folder/log.info.txt
else
	echo -e -n $err > $out_folder/log.error.txt
fi

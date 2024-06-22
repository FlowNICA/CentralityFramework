#!/bin/bash

#
#$ -wd /weekly/$USER/Soft/Centrality/TMP/
#$ -cwd
#$ -N run_centrality_fit
#$ -q all.q
#$ -l h_rt=8:00:00
#$ -l s_rt=8:00:00
#$ -t 1-200
#
#$ -o /weekly/$USER/Soft/Centrality/TMP/
#$ -e /weekly/$USER/Soft/Centrality/TMP/
#

export JOB_ID=${JOB_ID}
export TASK_ID=${SGE_TASK_ID}

START_DIR=$PWD

# Where main script & template config file are stored
MAIN_DIR=/weekly/parfenov/Soft/Centrality/scripts/new/7gev_reco_urqmd_nosecondary

# Creating output
COMMIT=new_recoUrQMD_STARlike_7gev_NoSecondary
OUT=/weekly/parfenov/Soft/Centrality/OUT/${COMMIT}/${JOB_ID}
OUT_FILE=$OUT/file
OUT_ROOT=$OUT_FILE/root
OUT_PDF=$OUT_FILE/pdf
OUT_LOG=$OUT/log

mkdir -p $OUT_ROOT/glauber_qa
mkdir -p $OUT_ROOT/fit
mkdir -p $OUT_ROOT/b_test
#mkdir -p $OUT_PDF
mkdir -p $OUT_LOG

LOG=$OUT_LOG/JOB_${JOB_ID}_${TASK_ID}.log
echo "Job output." &> $LOG
echo "-------------------------------" &>> $LOG

PARAMETER_LIST=$MAIN_DIR/parameter.list

# Read N-th line from the parameter list. N = Number of the job in job array
PARAMETER_LINE=`sed "${TASK_ID}q;d" ${PARAMETER_LIST}`
# PARAMETER_LINE=`sed "250q;d" ${PARAMETER_LIST}`


# Parsing parameters from the parameter list
PARAMETER_f0=${PARAMETER_LINE%%:*[0-9]}
TAIL=${PARAMETER_LINE#[0-9]*:}
PARAMETER_f1=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_k0=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_k1=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_mult0=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_mult1=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}

echo "Parameter list: ${PARAMETER_LINE}" &>>$LOG
echo "Reading:" &>>$LOG
echo " " &>>$LOG
echo "f0 = ${PARAMETER_f0}" &>>$LOG
echo "f1 = ${PARAMETER_f1}" &>>$LOG
echo "k0 = ${PARAMETER_k0}" &>>$LOG
echo "k1 = ${PARAMETER_k1}" &>>$LOG
echo "mult_min = ${PARAMETER_mult0}" &>>$LOG
echo "mult_max = ${PARAMETER_mult1}" &>>$LOG
echo "-------------------------------" &>> $LOG

# Creating tmp directory & log file
TMPALL=/weekly/parfenov/TMP/
TMP=$TMPALL/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $TMP

# Copy config file to TMP directory
cp ${MAIN_DIR}/config.txt.template ${TMP}/config.txt

# Replacing template fit params with specific variable values
sed -e "s|fminfmin|${PARAMETER_f0}|" -i ${TMP}/config.txt
sed -e "s|fmaxfmax|${PARAMETER_f1}|" -i ${TMP}/config.txt
sed -e "s|kminkmin|${PARAMETER_k0}|" -i ${TMP}/config.txt
sed -e "s|kmaxkmax|${PARAMETER_k1}|" -i ${TMP}/config.txt
sed -e "s|multminmultmin|${PARAMETER_mult0}|" -i ${TMP}/config.txt
sed -e "s|multmaxmultmax|${PARAMETER_mult1}|" -i ${TMP}/config.txt

cat ${TMP}/config.txt &>>$LOG
echo "-------------------------------" &>> $LOG

# Sourcing ROOT
source /weekly/parfenov/Soft/MPDROOT/mpdroot_140920/build/config.sh

# Where is CentralityFramework stored
CENTRALITY_FRAMEWORK_DIR=/weekly/parfenov/Soft/Centrality/Framework/centrality-master

# Compile binaries
cd ${TMP}/
cmake $CENTRALITY_FRAMEWORK_DIR/ &>>$LOG
make &>>$LOG
echo "-------------------------------" &>> $LOG

# Do main program
./glauber ./config.txt &>>$LOG
echo "-------------------------------" &>> $LOG

# Copy output files into output directory
mv ${TMP}/glauber_qa.root $OUT_ROOT/glauber_qa/glauber_qa_${JOB_ID}_${TASK_ID}.root
mv ${TMP}/fit*.root $OUT_ROOT/fit/fit_${JOB_ID}_${TASK_ID}.root
mv ${TMP}/b_test.root $OUT_ROOT/b_test/b_test_${JOB_ID}_${TASK_ID}.root

cd $START_DIR

# Delete temporary directory
rm -rf ${TMP}
echo "Done!" &>> $LOG

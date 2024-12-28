#!/bin/bash

#
#SBATCH -D /scratch1/$USER/Soft/Centrality/TMP/
#SBATCH -J run_centrality_fit
#SBATCH -p all.q
#SBATCH --time=8:00:00
#SBATCH -a 1-200
#
#SBATCH -o /scratch1/$USER/Soft/Centrality/TMP/slurm_%A_%a.out
#SBATCH -e /scratch1/$USER/Soft/Centrality/TMP/slurm_%A_%a.out
#

ls /eos/nica
sleep 20

ls /cvmfs/nica.jinr.ru/
sleep 20

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.12.24-1

export JOB_ID=$SLURM_ARRAY_JOB_ID
export TASK_ID=$SLURM_ARRAY_TASK_ID

# Input file with MC-Glauber tree & tree name
export inFileGlauberName=
export inTreeGlauberName=

# Input file with multiplicity histogram & histogram name
export inFileDataName=
export inHistDataName=

START_DIR=$PWD

# Where main script & template config file are stored
MAIN_DIR=$START_DIR

# Where is CentralityFramework stored
CENTRALITY_FRAMEWORK_DIR=${MAIN_DIR}/../Framework/McGlauber/centrality-master

# Setting up tmp directory & log file
TMPALL=${MAIN_DIR}/TMP
TMP=$TMPALL/TMP_${JOB_ID}_${TASK_ID}
mkdir -p ${TMP}

# Creating output
COMMIT=System_Energy_ShortInfo
OUT=${MAIN_DIR}/OUT/${COMMIT}/${JOB_ID}
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

PARAMETER_LIST=$MAIN_DIR/../scripts/template/parameter.list

# Read N-th line from the parameter list. N = Number of the job in job array
PARAMETER_LINE=`sed "${TASK_ID}q;d" ${PARAMETER_LIST}`

# Parsing parameters from the parameter list
PARAMETER_f0=${PARAMETER_LINE%%:*[0-9]}
TAIL=${PARAMETER_LINE#[0-9]*:}
PARAMETER_f1=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_k0=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_k1=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_p0=${TAIL%%:*[0-9]}
TAIL=${TAIL#[0-9]*:}
PARAMETER_p1=${TAIL%%:*[0-9]}
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
echo "p0 = ${PARAMETER_p0}" &>>$LOG
echo "p1 = ${PARAMETER_p1}" &>>$LOG
echo "mult_min = ${PARAMETER_mult0}" &>>$LOG
echo "mult_max = ${PARAMETER_mult1}" &>>$LOG
echo "-------------------------------" &>> $LOG

# Copy config file to TMP directory
cp ${MAIN_DIR}/../scripts/template/config.c.template ${TMP}/config.c

# Replacing template fit params with specific variable values
sed -e "s|fminfmin|${PARAMETER_f0}|" -i ${TMP}/config.c
sed -e "s|fmaxfmax|${PARAMETER_f1}|" -i ${TMP}/config.c
sed -e "s|kminkmin|${PARAMETER_k0}|" -i ${TMP}/config.c
sed -e "s|kmaxkmax|${PARAMETER_k1}|" -i ${TMP}/config.c
sed -e "s|pminpmin|${PARAMETER_p0}|" -i ${TMP}/config.c
sed -e "s|pmaxpmax|${PARAMETER_p1}|" -i ${TMP}/config.c
sed -e "s|multminmultmin|${PARAMETER_mult0}|" -i ${TMP}/config.c
sed -e "s|multmaxmultmax|${PARAMETER_mult1}|" -i ${TMP}/config.c
sed -e "s|Path-to-Glauber-file|${inFileGlauberName}|" -i ${TMP}/config.c
sed -e "s|Mc-Glauber-tree-name|${inTreeGlauberName}|" -i ${TMP}/config.c
sed -e "s|Path-to-RefMult-file|${inFileDataName}|" -i ${TMP}/config.c
sed -e "s|Mult-histo-name|${inHistDataName}|"      -i ${TMP}/config.c

cat ${TMP}/config.c &>>$LOG
echo "-------------------------------" &>> $LOG

# Compile binaries
cd ${TMP}/
cmake $CENTRALITY_FRAMEWORK_DIR/ -Duse_multithreading=OFF &>>$LOG
make &>>$LOG
echo "-------------------------------" &>> $LOG

# Do main program
./glauber ./config.c &>>$LOG
echo "-------------------------------" &>> $LOG

# Copy output files into output directory
mv ${TMP}/glauber_qa.root $OUT_ROOT/glauber_qa/glauber_qa_${JOB_ID}_${TASK_ID}.root
mv ${TMP}/fit*.root $OUT_ROOT/fit/fit_${JOB_ID}_${TASK_ID}.root
mv ${TMP}/b_test.root $OUT_ROOT/b_test/b_test_${JOB_ID}_${TASK_ID}.root

cd $START_DIR

# Delete temporary directory
rm -rf ${TMP}
echo "Done!" &>> $LOG

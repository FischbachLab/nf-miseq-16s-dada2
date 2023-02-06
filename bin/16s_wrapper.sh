#!/bin/bash -x

set -euoE pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
DEFAULT=16
coreNum=${coreNum:-$DEFAULT}

if [ -z "${coreNum:-}" ]; then
  echo "coreNum was not set"
fi

# s3 inputs from env variables
#S3INPUTPATH="${1}"
#S3OUTPUTPATH="${2}"
#DB="${3}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_TMP="${OUTPUTDIR}/tmp"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
#DB="${LOCAL}/silvaDB"
INPUT_DIR="${LOCAL_TMP}/data"

mkdir -p ${OUTPUTDIR} ${LOCAL_OUTPUT} ${LOG_DIR}
mkdir -p ${LOCAL_TMP} ${INPUT_DIR}
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

# download 16s database
#aws s3 cp --quiet s3://maf-users/Xiandong_Meng/16S/silvaDB/silva_nr99_v138_train_set.fa.gz ${DB}

# download samples
aws s3 sync --quiet $S3INPUTPATH ${INPUT_DIR}

# run DADA2 pipeline
#input parameter: ouptut dir & input fastq dir
#timem ${scriptFolder}/post_processing.sh ${LOCAL_OUTPUT} ${INPUT_DIR} ${DB} &>> ${LOG_DIR}/16s.log.txt
16S.R ${LOCAL_OUTPUT} ${INPUT_DIR} ${DB}
fix_summary.py ${LOCAL_OUTPUT}/Sample_stats.tsv ${LOCAL_OUTPUT}/ASVs_taxonomy_dada2.tsv ${LOCAL_OUTPUT}/ASVs_counts.tsv

mv ${LOCAL}/Rplots.pdf ${LOCAL_OUTPUT}
mv ${LOCAL}/DADA2_summary ${LOCAL_OUTPUT}

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet ${LOCAL_OUTPUT} ${S3OUTPUTPATH}
rm -rf ${OUTPUTDIR}

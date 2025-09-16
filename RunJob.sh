#!/bin/bash

HOMEDIR=/cephfs/user/s6subans/ChargedHiggsAna_Boost/Code

cd ${HOMEDIR}

source /etc/profile

setupATLAS
lsetup  "root 6.18.04-x86_64-centos7-gcc8-opt"


IN_PATH=${1}
FILE=${2}
WP=${3}
OUTDIR=${4}
Cluster=${5}
Process=${6}

echo "${IN_PATH}"
echo "${FILE}"
echo "${WP}"
echo "${OUTDIR}"


./execute_V2.exe ${IN_PATH} ${FILE} ${WP} ${OUTDIR}
echo "Write Output File to: ${OUTDIR} ....."


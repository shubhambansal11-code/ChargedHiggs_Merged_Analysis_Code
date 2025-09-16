#!/bin/bash

#Submit script for gbb Tuple Ana
echo "Submit jobs for Charged Higgs Analysis......"

LOG_FOLDER="/cephfs/user/s6subans/ChargedHiggsLog/"

echo "Logs in: ${LOG_FOLDER}"
echo "Output in: ${OUTPATH}"

OUTDIR="/cephfs/user/s6subans/ChargedHiggsBoost/"


WP=(
#"77p"
"77p"
#"85p"
#"60p"
)

export EOS_MGM_URL=root://eosuser.cern.ch
echo "reading files !!! "
Paths=(  
"HpWh_BoostSample/MC16a/"
"HpWh_BoostSample/MC16d/"
"HpWh_BoostSample/MC16e/"
)

File=(
hp800_AFII.root   
hp1200_AFII.root  
hp1600_AFII.root  
hp2000_AFII.root  
hp2500_AFII.root  
hp3000_AFII.root  
)

echo "sucessfuly opened tar files"
for path in "${Paths[@]}"
do
   for file in "${File[@]}"
   do
   condor_submit IN_PATH="${path}" FILE="${file}" WP="${WP}" OUTDIR="${OUTDIR}" /cephfs/user/s6subans/ChargedHiggsAna_Boost/Code/run_Hplus.sub
   done
done

echo "all done !!! "

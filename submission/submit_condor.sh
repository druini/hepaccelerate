#!/bin/bash
if [[ $USER == *"druini"* ]]; then
    export X509_USER_PROXY=/afs/cern.ch/user/d/druini/.x509up_proxy
    source /afs/cern.ch/work/d/druini/miniconda3/etc/profile.d/conda.sh
    export PATH=/afs/cern.ch/work/d/druini/miniconda3/bin:$PATH
    source activate hepaccelerate_cpu
    cd /afs/cern.ch/work/d/druini/public/hepaccelerate
else
    export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148
    source /afs/cern.ch/work/a/algomez/miniconda3/etc/profile.d/conda.sh
    export PATH=/afs/cern.ch/work/a/algomez/miniconda3/bin:$PATH
    source activate hepaccelerate_cpu
    cd /afs/cern.ch/work/a/algomez/ttH/CMSSW_10_6_5/src/TTH/Analyzer/hepaccelerate/
fi

echo $USER
echo "submitting sample $1 $2"

if [[ $1 == *"Run"* ]]; then
  folder=Nano25Oct2019
  sampleFile=$1
else
    if [[ $2 == *"2016"* ]]; then
        folder=RunIISummer16NanoAODv6
    elif [[ $2 == *"2017"* ]]; then
        folder=RunIIFall17NanoAODv6
    else
        folder=RunIIAutumn18NanoAODv6
    fi
    sampleFile=$1_$2
fi

time PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py \
  --filelist ${PWD}/datasets/$folder/$sampleFile.txt \
  --sample $1  \
  --year $2 \
  --outdir ${PWD}/results/v04/$2/ \
  --boosted \
  --categories all
  #--cache-location /eos/user/d/druini/cache/ \
  #--from-cache
  #--filelist /afs/cern.ch/work/d/druini/public/hepaccelerate/datasets/algomez/$1.txt \

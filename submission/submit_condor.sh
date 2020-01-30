#!/bin/bash
export X509_USER_PROXY=/afs/cern.ch/user/d/druini/.x509up_proxy
source /afs/cern.ch/work/d/druini/miniconda3/etc/profile.d/conda.sh
export PATH=/afs/cern.ch/work/d/druini/miniconda3/bin:$PATH
source activate hepaccelerate_cpu
#cd /afs/cern.ch/user/a/algomez/workingArea/ttH/hepaccelerate
cd /afs/cern.ch/work/d/druini/public/hepaccelerate
    
echo "submitting sample $1"

if [[ $1 == *"Run"* ]]; then
  folder=Nano25Oct2019
else
  folder=RunIIFall17NanoAODv5
fi

PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py \
  --filelist /afs/cern.ch/work/d/druini/public/hepaccelerate/datasets/$folder/$1.txt \
  --sample $1  \
  --outdir /afs/cern.ch/work/d/druini/public/hepaccelerate/results/tau21_045 \
  --boosted \
  --categories all \
  --cache-location /eos/user/d/druini/cache/ \
  --from-cache 
  #--filelist /afs/cern.ch/work/d/druini/public/hepaccelerate/datasets/algomez/$1.txt \

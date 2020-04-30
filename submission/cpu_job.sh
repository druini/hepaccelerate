#!/bin/bash

#SBATCH -p wn

#SBATCH --account=cn-test

#SBATCH --job-name=cpu_ana

#SBATCH -o slurm_output/slurm-%j.out

echo "submitting sample $1"

if [[ $1 == *"Run"* ]]; then
  folder=Nano25Oct2019
else
  folder=RunIIFall17NanoAODv5
fi

PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py \
  --filelist /afs/cern.ch/work/d/druini/public/hepaccelerate/datasets/$folder/$1.txt \
  --sample $1  \
  --outdir /afs/cern.ch/work/d/druini/public/hepaccelerate/results/Wmass_60 \
  --boosted \
  --categories all \
  --cache-location /eos/user/d/druini/cache/ \
  --from-cache

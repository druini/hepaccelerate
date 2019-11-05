#!/bin/bash

#SBATCH -p wn

#SBATCH --account=cn-test

#SBATCH --job-name=cpu_ana

#SBATCH -o slurm_output/slurm-%j.out

echo "submitting sample $1"

PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py \
  --filelist datasets/algomez/$1.txt \
  --sample $1  \
  --outdir tests/ \
  --boosted \
  --categories boosted_higgs_only boosted_higgs_or_W \
  --cache-location ~/eos/cache/ \
  --from-cache 

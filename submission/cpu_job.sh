#!/bin/bash

#SBATCH -p wn

#SBATCH --account=cn-test

#SBATCH --job-name=cpu_analysis                 

#SBATCH -o slurm_output/slurm-%j.out

echo "submitting sample $1"

PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py \
  --filelist datasets/algomez/$1.txt \
  --sample $1  \
  --outdir results/ \
  --boosted \
  --cache-location ~/eos/cache/ \
  --from-cache 

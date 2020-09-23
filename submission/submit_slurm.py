#!/usr/bin/env python

import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
parser.add_argument('--datasets', help='directory with list of inputs', type=str, default='datasets')
parser.add_argument('--samples', nargs='+', help='List of samples to process', type=str, default=None)
parser.add_argument('--files-per-job', action='store', help='Number of files to process per job', type=int, default=5)
parser.add_argument('--postproc', action='store_true', help='Flag for running on postprocessed datasets and include corrections')
parser.add_argument('-q','--quick', action='store_true', help='submit jobs to the quick queue')

parser.add_argument('--parameters', nargs='+', help='change default parameters, syntax: name value, eg --parameters met 40 bbtagging_algorithm btagDDBvL', default=None)
parser.add_argument('--from-cache', action='store_true', help='Load from cache (otherwise create it)')
parser.add_argument('--files-per-batch', action='store', help='Number of files to process per batch', type=int, default=1, required=False)
parser.add_argument('--nthreads', action='store', help='Number of CPU threads to use', type=int, default=4, required=False)
parser.add_argument('--cache-location', action='store', help='Path prefix for the cache, must be writable', type=str, default=os.path.join(os.getcwd(), 'cache'))
parser.add_argument('--outdir', action='store', help='directory to store outputs', type=str, default=os.path.join(os.getcwd(),'results'))
parser.add_argument('--categories', nargs='+', help='categories to be processed (default: all)', default="all")
parser.add_argument('--year', action='store', choices=['2016', '2017', '2018'], help='Year of data/MC samples', default='2017')
parser.add_argument('--version', action='store', help='tag added to the output directory', type=str, default='')
parser.add_argument('--dryrun', action='store_true', help='Option to only create job files, without submitting them')
args = parser.parse_args()

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.makedirs(dir)

def partitionFileList(filelist, chunkSize=1):
    sampleFileList = np.loadtxt(filelist, dtype=str, ndmin=1)
    return [sampleFileList[i:i+chunkSize] for i in range(0, len(sampleFileList), chunkSize)]   

import time
timestr = time.strftime("%Y%m%d-%H%M%S")
ver = '' if args.version=='' else f'{args.version}_'
job_directory = os.path.join(os.getcwd(),'logs_submission',f'{ver}{args.year}_{timestr}')

# Make top level directories
mkdir_p(job_directory)

if args.samples is None: 
  args.samples = [f'submission/allSamples_2018.txt']
if len(args.samples)==1 and args.samples[0].endswith('.txt'):
  samples = [l.strip() for l in open(args.samples[0]).readlines() if not l.startswith(('#','Single'))]
else:
  samples = args.samples

if not type(args.categories)==list:
    categories = [args.categories]
else:
    categories = args.categories

for s in samples:
    sample_directory = os.path.join(job_directory,s)
    mkdir_p(sample_directory)
    if 'Single' in s:
      #subdir = 'Nano02Apr2020'
      subdir = 'Nano25Oct2019'
    else:
      if args.year=='2016':
        subdir = 'RunIISummer16NanoAODv7'
      elif args.year=='2017':
        subdir = 'RunIIFall17NanoAODv7'
      else:
        subdir = 'RunIIAutumn18NanoAODv7'
    if args.postproc: subdir += 'PostProc'
    filelist = os.path.join(args.datasets,subdir,f'{s}_{args.year}.txt')
    
    if not os.path.isfile(filelist):
      print(f'cannot find {filelist}, skipping...\n')
      continue
    par_files = partitionFileList(filelist, args.files_per_job * args.files_per_batch)

    for njob,f in enumerate(par_files): 
        job_file = os.path.join(sample_directory,f'{njob}.job')

        with open(job_file, "w") as fh:
            fh.write("#!/bin/bash\n")
            fh.write(f"#SBATCH --job-name={s}_{njob}.job\n")
            fh.write("#SBATCH -p wn\n")
            fh.write(f"#SBATCH -o {sample_directory}/slurm-{njob}.out\n")
            fh.write(f"#SBATCH -e {sample_directory}/slurm-{njob}.err\n")

            #fh.write("mkdir /scratch/c/\n")
            fh.write(f'echo sample: {s}\n') 
            fh.write(f"time PYTHONPATH=hepaccelerate:coffea:. python {os.getcwd()}/run_analysis.py ")
            fh.write("--categories ")
            fh.write(' '.join(map(str, categories)))
            if args.parameters is not None:
              fh.write(f' --parameters ')
              fh.write(' '.join(map(str, args.parameters)))
            if args.version!='':
                fh.write(f' --version {args.version} ')
            if args.from_cache:
                fh.write(" --from-cache ")
            if args.postproc:
              fh.write(' --corrections ')
            fh.write(f" --boosted --sample {s} --files-per-batch {args.files_per_batch} --nthread {args.nthreads}  --year {args.year} --outtag _{njob} --outdir {os.path.join(args.outdir,args.year)} ")
            for fi in f:
              local = fi.replace('root://xrootd-cms.infn.it/','/pnfs/psi.ch/cms/trivcat')
              if os.path.isfile(local):
                fh.write(f'{local} ')
              else:
                fh.write(f'{fi.replace("xrootd-cms.infn.it","cms-xrd-global.cern.ch")} ')

        os.system(f"sbatch {'-p quick' if args.quick else ''} --mem=4000 {job_file}")

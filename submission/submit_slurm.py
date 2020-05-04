#!/usr/bin/env python

import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
parser.add_argument('--datasets', help='directory with list of inputs', type=str, default='datasets')
parser.add_argument('--samples', nargs='+', help='List of samples to process', type=str, default=None, required=True)
parser.add_argument('--files-per-job', action='store', help='Number of files to process per job', type=int, default=10, required=False)
parser.add_argument('--batchSystem', help='batch system to submit to', type=str, default='slurm_cpu', required=False)

parser.add_argument('--from-cache', action='store_true', help='Load from cache (otherwise create it)')
parser.add_argument('--files-per-batch', action='store', help='Number of files to process per batch', type=int, default=1, required=False)
parser.add_argument('--nthreads', action='store', help='Number of CPU threads to use', type=int, default=4, required=False)
parser.add_argument('--cache-location', action='store', help='Path prefix for the cache, must be writable', type=str, default=os.path.join(os.getcwd(), 'cache'))
parser.add_argument('--outdir', action='store', help='directory to store outputs', type=str, default=os.path.join(os.getcwd(),'results'))
parser.add_argument('--DNN', action='store', choices=['save-arrays','cmb_binary', 'cmb_multiclass', 'ffwd_binary', 'ffwd_multiclass',False, 'mass_fit'], help='options for DNN evaluation / preparation', default=False)
parser.add_argument('--categories', nargs='+', help='categories to be processed (default: all)', default="all")
parser.add_argument('--path-to-model', action='store', help='path to DNN model', type=str, default=None, required=False)
parser.add_argument('--year', action='store', choices=['2016', '2017', '2018'], help='Year of data/MC samples', default='2017')
args = parser.parse_args()

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def partitionFileList(filelist, chunkSize=1):
    sampleFileList = np.loadtxt(filelist, dtype=str)
    return [sampleFileList[i:i+chunkSize] for i in range(0, len(sampleFileList), chunkSize)]   

import time
timestr = time.strftime("%Y%m%d-%H%M%S")
job_directory = os.path.join(os.getcwd(),f'logs_{timestr}')

# Make top level directories
mkdir_p(job_directory)

if len(args.samples)==1 and args.samples[0].endswith('.txt'):
  samples = [l.strip() for l in open(args.samples[0]).readlines()]
else:
  samples = args.samples

if not type(args.categories)==list:
    categories = [args.categories]
else:
    categories = args.categories

for s in samples:

    sample_directory = os.path.join(job_directory,s)
    mkdir_p(sample_directory)
    if 'Run' in s:
      subdir = 'Nano25Oct2019'
    else:
      if args.year=='2016':
        subdir = 'RunIISummer16NanoAODv6'
      elif args.year=='2017':
        subdir = 'RunIIFall17NanoAODv6'
      else:
        subdir = 'RunIIAutumn18NanoAODv6'
    par_files = partitionFileList(os.path.join(args.datasets,subdir,f'{s}.txt'), args.files_per_job * args.files_per_batch)

    for njob,f in enumerate(par_files): 

        job_file = os.path.join(sample_directory,f'{njob}.job')

        if args.batchSystem=="slurm_cpu":

            with open(job_file, "w") as fh:
                fh.write("#!/bin/bash\n")
                fh.write(f"#SBATCH --job-name={s}_{njob}.job\n")
                fh.write("#SBATCH -p wn\n")
                fh.write(f"#SBATCH -o {sample_directory}/slurm-{njob}.out\n")
                fh.write(f"#SBATCH -e {sample_directory}/slurm-{njob}.err\n")

                #fh.write("mkdir /scratch/c/\n")
                fh.write(f"PYTHONPATH=hepaccelerate:coffea:. python {os.getcwd()}/run_analysis.py ")
                fh.write("--categories ")
                fh.write(' '.join(map(str, categories)))
                fh.write(f" --boosted --sample {s} --files-per-batch {args.files_per_batch} --nthread {args.nthreads}  --outdir {os.path.join(args.outdir,args.year)} --year {args.year} --outtag _{njob} ")
                if args.DNN:
                    fh.write(f"--DNN {args.DNN} ")
                if args.from_cache:
                    fh.write("--from-cache ")
                fh.write(' '.join(map(str, f)))
                #fh.write('\n')
                #fh.write('rm -r /scratch/c/')

            os.system(f"sbatch {job_file}")

        elif args.batchSystem=="slurm_gpu":

            with open(job_file, "w") as fh:
                fh.write("#!/bin/bash\n")
                fh.write("#SBATCH --job-name={0}_{1}.job\n".format(s,njob))
                fh.write("#SBATCH --account=gpu_gres  # to access gpu resources\n")
                fh.write("#SBATCH --nodes=1       # request to run job on single node\n")
                fh.write("#SBATCH --ntasks=5     # request 10 CPU's (t3gpu01: balance between CPU and GPU : 5CPU/1GPU)\n")
                fh.write("#SBATCH --gres=gpu:1    # request 1 GPU's on machine\n")

                fh.write("PYTHONPATH=hepaccelerate:coffea:. python3 {0}/run_analysis.py ".format(os.getcwd()))
                fh.write("--categories ")
                fh.write(' '.join(map(str, categories)))
                fh.write(" --sample {0} --files-per-batch {1} --nthread {2} --cache-location {3} --outdir {4} --path-to-model {5} --year {6} ".format(s, args.files_per_batch, args.nthreads, args.cache_location, args.outdir, args.path_to_model, args.year))
                if args.DNN:
                    fh.write("--DNN {0} ".format(args.DNN))
                if args.from_cache:
                    fh.write("--from-cache ")
                if args.cache_only:
                    fh.write("--cache-only ")
                fh.write(' '.join(map(str, f)))

            os.system("sbatch -o {0}/slurm-{1}.out {2}".format(sample_directory,njob,job_file))

        else:
            print("Unknown batch system.")  

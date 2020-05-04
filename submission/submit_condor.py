#!/usr/bin/env python

import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
parser.add_argument('--datasets', help='directory with list of inputs', type=str, default='datasets')
parser.add_argument('--samples', nargs='+', help='List of samples to process', type=str, default=None)
parser.add_argument('--files-per-job', action='store', help='Number of files to process per job', type=int, default=5, required=False)
parser.add_argument('--rerun-failed-from', dest='rerun_failed', action='store_true', help='Search a submission directory for a list of failed .root files and reruns those files')

parser.add_argument('--from-cache', action='store_true', help='Load from cache (otherwise create it)')
parser.add_argument('--files-per-batch', action='store', help='Number of files to process per batch', type=int, default=1, required=False)
parser.add_argument('--nthreads', action='store', help='Number of CPU threads to use', type=int, default=4, required=False)
parser.add_argument('--cache-location', action='store', help='Path prefix for the cache, must be writable', type=str, default=os.path.join(os.getcwd(), 'cache'))
parser.add_argument('--outdir', action='store', help='directory to store outputs', type=str, default=os.path.join(os.getcwd(),'results'))
parser.add_argument('--DNN', action='store', choices=['save-arrays','cmb_binary', 'cmb_multiclass', 'ffwd_binary', 'ffwd_multiclass',False, 'mass_fit'], help='options for DNN evaluation / preparation', default=False)
parser.add_argument('--categories', nargs='+', help='categories to be processed (default: all)', default="all")
parser.add_argument('--path-to-model', action='store', help='path to DNN model', type=str, default=None, required=False)
parser.add_argument('--year', action='store', choices=['2016', '2017', '2018'], help='Year of data/MC samples', default='2017')
parser.add_argument('--version', action='store', help='tag added to the output directory', type=str, default='')
parser.add_argument('--dryrun', action='store_true', help='Option to only create job files, without submitting them')
args = parser.parse_args()

if not args.rerun_failed and args.samples==None:
  raise Exception('Must provide sample list')

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def partitionFileList(filelist, chunkSize=1):
    sampleFileList = np.loadtxt(filelist, dtype=str)
    return [sampleFileList[i:i+chunkSize] for i in range(0, len(sampleFileList), chunkSize)]   

import time
timestr = time.strftime("%Y%m%d-%H%M%S")
job_directory = os.path.join(os.getcwd(),f'logs_submission_{timestr}')

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
    submission_file = os.path.join(sample_directory, f'{s}_condor.sub')
    with open(submission_file, 'w') as sf:
      sf.write('universe = vanilla\n')
      sf.write(f'executable = {os.path.join(sample_directory,"$Fnx(script)")}\n')
      condor_outname = f'{s}_$Fn(script)'
      sf.write(f'log = {os.path.join(sample_directory, f"{condor_outname}.log")}\n')
      sf.write(f'error = {os.path.join(sample_directory, f"{condor_outname}.err")}\n')
      sf.write(f'output = {os.path.join(sample_directory, f"{condor_outname}.out")}\n')
      sf.write('getenv = True\n')
      #sf.write('+JobFlavour = "espresso"\n')
      sf.write('+JobFlavour = "tomorrow"\n')
      sf.write(f'queue script matching files ({os.path.join(sample_directory,"*.sh")})')

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

        job_file = os.path.join(sample_directory,f'{njob}.sh')

        with open(job_file, "w") as fh:
            fh.write("#!/bin/bash\n")
            fh.write('if [[ $USER == *"druini"* ]]; then\n')
            fh.write('    export PATH=/afs/cern.ch/cms/caf/scripts:/cvmfs/cms.cern.ch/common:/usr/sue/bin:/usr/lib64/qt-3.3/bin:/usr/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/afs/cern.ch/user/d/druini/bin\n')
            fh.write('    export X509_USER_PROXY=/afs/cern.ch/user/d/druini/.x509up_proxy\n')
            fh.write('    source /afs/cern.ch/work/d/druini/miniconda3/etc/profile.d/conda.sh\n')
            fh.write('    export PATH=/afs/cern.ch/work/d/druini/miniconda3/bin:$PATH\n')
            fh.write('    source activate hepaccelerate_cpu\n')
            fh.write('    cd /afs/cern.ch/work/d/druini/public/hepaccelerate\n')
            fh.write('else\n')
            fh.write('    export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148\n')
            fh.write('    source /afs/cern.ch/work/a/algomez/miniconda3/etc/profile.d/conda.sh\n')
            fh.write('    export PATH=/afs/cern.ch/work/a/algomez/miniconda3/bin:$PATH\n')
            fh.write('    source activate hepaccelerate_cpu\n')
            fh.write('    cd /afs/cern.ch/user/a/algomez/workingArea/ttH/hepaccelerate\n')
            fh.write('fi\n')

            fh.write('python --version\n')
            fh.write('echo PATH=$PATH\n')
            fh.write('which python\n')
            fh.write(f"time PYTHONPATH=hepaccelerate:coffea:. python {os.getcwd()}/run_analysis.py ")
            fh.write("--categories ")
            fh.write(' '.join(map(str, categories)))
            fh.write(f" --boosted --sample {s} --files-per-batch {args.files_per_batch} --nthread {args.nthreads}  --outdir {os.path.join(args.outdir,args.year)} --year {args.year} --outtag _{njob} --dir-for-fails {sample_directory} --version {args.version} ")
            if args.DNN:
                fh.write(f"--DNN {args.DNN} ")
            if args.from_cache:
                fh.write("--from-cache ")
            for fi in f:
                fh.write(f'{fi} ')
            #fh.write(f' '.join(f))
            #fh.write('\n')
            #fh.write('rm -r /scratch/c/')

    if not args.dryrun:
      os.system(f"condor_submit {submission_file}")

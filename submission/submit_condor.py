#!/usr/bin/env python

'''
USAGE
python submission/submit_condor.py --samples ttHTobb TTToSemiLeptonic --version v07 --year 2017
'''

import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
parser.add_argument('--datasets', help='directory with list of inputs', type=str, default='datasets')
parser.add_argument('--samples', nargs='+', help='List of samples to process', type=str, default=None)
parser.add_argument('--files-per-job', action='store', help='Number of files to process per job', type=int, default=5)
parser.add_argument('--condor-job-flavour', choices=['espresso','microcentury','longlunch','workday','tomorrow','testmatch','nextweek'], help='condor queue flavour as in https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor#Queue_Flavours', default='tomorrow')

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
    sampleFileList = [l.strip() for l in open(filelist).readlines()]
    return [sampleFileList[i:i+chunkSize] for i in range(0, len(sampleFileList), chunkSize)]   

import time
timestr = time.strftime("%Y%m%d-%H%M%S")
ver = '' if args.version=='' else f'{args.version}_'
job_directory = os.path.join(os.getcwd(),'logs_submission',f'{ver}{timestr}')

# Make top level directories
mkdir_p(job_directory)

if args.samples is None: 
  args.samples = f'submission/allSamples_{args.year}.txt'
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
    sf.write(f'+JobFlavour = \"{args.condor_job_flavour}\"\n')
    sf.write(f'queue script matching files ({os.path.join(sample_directory,"*.sh")})')

  if 'Run' in s:
    subdir = 'Nano02Apr2020'
  else:
    if args.year=='2016':
      subdir = 'RunIISummer16NanoAODv7'
    elif args.year=='2017':
      subdir = 'RunIIFall17NanoAODv7'
    else:
      subdir = 'RunIIAutumn18NanoAODv7'
  filelist = os.path.join(args.datasets,subdir,f'{s}_{args.year}.txt')
  par_files = partitionFileList(filelist, args.files_per_job * args.files_per_batch)

  for njob,f in enumerate(par_files): 
    job_file = os.path.join(sample_directory,f'{njob}.sh')

    with open(job_file, "w") as fh:
      fh.write("#!/bin/bash\n")
      fh.write('if [[ $USER == *"druini"* ]]; then\n')
      fh.write('    export PATH=/afs/cern.ch/cms/caf/scripts:/cvmfs/cms.cern.ch/common:/usr/sue/bin:/usr/lib64/qt-3.3/bin:/usr/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/afs/cern.ch/user/d/druini/bin\n')
      fh.write('    export X509_USER_PROXY=/afs/cern.ch/user/d/druini/.x509up_proxy\n')
      fh.write('    export XRD_NETWORKSTACK=IPv4\n')
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
      #fh.write('echo PATH=$PATH\n')
      #fh.write('which python\n')
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
      fh.write(f" --boosted --sample {s} --files-per-batch {args.files_per_batch} --nthread {args.nthreads}  --year {args.year} --outtag _{njob} --outdir {os.path.join(args.outdir,args.year)} ")
      for fi in f:
          fh.write(f'{fi} ')
      #fh.write(f' '.join(f))
      #fh.write('\n')
      #fh.write('rm -r /scratch/c/')

  if not args.dryrun:
    os.system(f"condor_submit {submission_file}")

import os, argparse
from glob import glob

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-j','--job-directory', type=str, required=True)
  parser.add_argument('-o','--output-directory', type=str, required=True)
  parser.add_argument('-s','--samples', default='submission/allSamples_2017.txt')
  args = parser.parse_args()

  samples = [l.strip() for l in open(args.samples).readlines()]
  for s in samples:
    jobs = glob(os.path.join(args.job_directory, s, '*job'))
#    jobs = [j for j in jobs if not 'rerunFail' in j]
    outs = glob(os.path.join(args.output_directory, f'*{s}*msd_nom*'))
#    outs = [o for o in outs if not 'rerunFail' in o]
    njobs = len(jobs)
    nouts = len(outs)
    if njobs!=nouts:
      print(f'{s}: {njobs} jobs and {nouts} outputs.')
      #nouts = [int(o.split('_')[-1].split('.')[0]) for o in outs]
      #failed = [j for j in jobs if not int(os.path.basename(j).split('.')[0]) in nouts]
      #print(failed)
      #submit_script = [l.strip() for l in open(os.path.join(args.job_directory, s, f'{s}_condor.sub')).readlines()]
      #new_submit = os.path.join(args.job_directory, s, f'{s}_condor_resubmit.sub')
      #with open(new_submit, 'w') as ns:
      #  for l in submit_script[:-2]:
      #    ns.write(l+'\n')
      #  ns.write('+JobFlavour = "tomorrow"\n')
      #  ns.write(f'queue script in ({" ".join(map(str, failed))})')
      #os.system(f"condor_submit {new_submit}")
      #for _ in range(2):
      #  print()

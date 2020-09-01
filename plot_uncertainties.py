import json,argparse,sys,os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from glob import glob

def load_mc(indir,histToLoad):
    json_file = glob( os.path.join(indir,'*ttHTobb*json') )
    if len(json_file)>1:
        json_file = input('Which file?\n'+'\n'.join(map(str,json_file))+'\n')
    elif len(json_file)==0:
        print('No file found...')
        sys.exit(0)
    else:
        json_file = json_file[0]
    with open(json_file) as f:
        return json.load(f)[histToLoad]

def rebin(bins, counts, yerr, rebin_factor):
    new_bins   = bins[::rebin_factor]
    new_counts = np.add.reduceat(counts, range(0, len(counts), rebin_factor))
    new_yerr   = np.add.reduceat(yerr, range(0, len(yerr), rebin_factor))
    return new_bins, new_counts, new_yerr

def plot_unc(indir, histName, uncName, outdir):
    hist = {
            'nominal' : load_mc(indir, histName),
            'up'      : load_mc(os.path.join(indir,f'{uncName}Up'), histName),
            'down'    : load_mc(os.path.join(indir,f'{uncName}Down'), histName)
            }

    color = {
            'nominal' : 'k',
            'up'      : 'r',
            'down'    : 'b'
            }
    ls    = {
            'nominal' : '-',
            'up'      : '--',
            'down'    : '-.'
            }

    rebin_factor = 5

    plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
    f, ax = plt.subplots()
    #ax = hep.cms.label(data=False, paper=False, year=args.year, ax=ax, loc=1)
    for hn,h in hist.items():
        h['edges'], h['contents'], h['contents_w2'] = rebin(h['edges'], h['contents'], h['contents_w2'], rebin_factor)
        hep.histplot(h['contents'], h['edges'], edges=True, label=(hn if hn=='nominal' else f'{uncName} {hn}'),ls=ls[hn],color=color[hn])
    plt.xlim(90,170)
    plt.legend()
    plt.xlabel('Leading AK8 jet softdrop mass [GeV]')
    plt.ylabel(f'Events / {rebin_factor} GeV')
    for ext in ['.png','.pdf']:
         plt.savefig(os.path.join(outdir,f'{args.variable}_{args.mask}_{uncName}{ext}'))
    #plt.show()

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  #parser.add_argument('-u', '--uncertainty', default='jer', help='which uncertainty to plot')
  parser.add_argument('-r', '--result-dir', default='results')
  parser.add_argument('-y', '--year', default='2017')
  parser.add_argument('-v', '--version')
  parser.add_argument('-s', '--selection', default='met20_btagDDBvL086')
  parser.add_argument('-p', '--path', default=None, help='path to target directory, overrides the r,y,v,s options')
  parser.add_argument('--variable', default='leadAK8JetMass')
  parser.add_argument('--mask', default='2J2WdeltaR_Pass')

  try: args = parser.parse_args()
  except:
    parser.print_help()
    sys.exit(0)
    
  for unc in ['jer','jesTotal','jmr','jms']:
      if args.path is not None:
          indir  = args.path
          outdir = os.path.join('plots', *args.path.split(os.path.sep)[1:])
      else:
          indir  = os.path.join(args.result_dir,args.year,args.version,args.selection)
          outdir = os.path.join('plots',args.year,args.version,args.selection)
      
      outdir = os.path.join(outdir,args.mask,'uncertainties')
      if not os.path.exists(outdir):
          os.makedirs(outdir)

      histName = f'hist_{args.variable}_{args.mask}'
      plot_unc(indir,histName, unc, outdir)
      #plot_unc(indir,histName, args.uncertainty, outdir)

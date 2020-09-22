import json,argparse,sys,os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from glob import glob
from itertools import cycle

def load_mc(indir,uncName,histToLoad):
    #json_file = glob( os.path.join(indir,f'*ttHTobb_{uncName}.json') )
    json_file = glob( os.path.join(indir,f'*signal_{uncName}*json') )
    if len(json_file)>1:
        json_file = input('Which file?\n'+'\n'.join(map(str,json_file))+'\n')
    elif len(json_file)==0:
        print(f'No file found in {indir}...')
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

def prepare_removeCorrections(indir, histName, correctionList):
  linestyle = cycle(["--","-.",":"])

  hist  = { name : load_mc(indir, f'{name}_merged', histName) for name in correctionList }
  color = { name : plt.cm.tab10(i) for i,name in enumerate(hist.keys()) }
  ls    = { name : next(linestyle) for name in hist }

  hist['msoftdrop_nom']  = load_mc(indir.replace(args.version, 'v10'), 'msd_nom_merged', histName)
  color['msoftdrop_nom'] = 'k'
  ls['msoftdrop_nom']    = '-'

  hist['msoftdrop_raw']  = load_mc(indir.replace(args.version, 'v09'), 'merged', histName)
  color['msoftdrop_raw'] = 'b'
  ls['msoftdrop_raw']    = '-' 

  return (hist, color, ls)

def plot_removeCorrections(hist, color, ls, outdir):
    rebin_factor = 5

    plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
    f, ax = plt.subplots()
    hep.cms.label(data=False, paper=False, year=args.year, ax=ax, loc=0)
    for hn,h in hist.items():
        h['edges'], h['contents'], h['contents_w2'] = rebin(h['edges'], h['contents'], h['contents_w2'], rebin_factor)

        hep.histplot(h['contents'], h['edges'], edges=True, label=hn,ls=ls[hn],color=color[hn])
        #print(uncName, hn, h['contents'][-3:])
    plt.xlim(90,170)
    #plt.xlim(0,1000)
    #ax.set_ylim(ymin=.5e-1)
    #plt.semilogy()
    plt.legend()
    #import pdb
    #pdb.set_trace()
    ax.set_xlabel('Leading AK8 jet softdrop mass [GeV]', ha='right', x=1)
    ax.set_ylabel(f'Events / {rebin_factor} GeV', ha='right', y=1)
    for ext in ['.png','.pdf']:
      plt.savefig(os.path.join(outdir,f'{args.variable}_{args.mask}{ext}'))

def prepare_unc(indir, histName, uncName, outdir):
    if uncName=='msd':
        hist = {
                'nominal' : load_mc(indir, 'msd_nanoAOD', histName),
                'up'      : load_mc(indir, 'msd_nom',     histName),
                'down'    : load_mc(indir, 'msd_raw',     histName)
                }
    else:
        hist = {
                'nominal' : load_mc(indir,  'msd_nom',   histName),
                'up'      : load_mc(indir, f'{uncName}Up',   histName),
                'down'    : load_mc(indir, f'{uncName}Down', histName)
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
    hep.cms.label(data=False, paper=False, year=args.year, ax=ax, loc=0)
    for hn,h in hist.items():
        h['edges'], h['contents'], h['contents_w2'] = rebin(h['edges'], h['contents'], h['contents_w2'], rebin_factor)
        if uncName=='msd':
            if   hn=='nominal': label = 'nanoAOD'
            elif hn=='up':      label = 'nom'
            elif hn=='down':    label = 'raw'
        else:
            if    hn=='nominal': label = hn
            else: label = f'{uncName} {hn}'

        hep.histplot(h['contents'], h['edges'], edges=True, label=label,ls=ls[hn],color=color[hn])
        #print(uncName, hn, h['contents'][-3:])
    plt.xlim(90,170)
    #plt.xlim(0,1000)
    #ax.set_ylim(ymin=.5e-1)
    #plt.semilogy()
    plt.legend()
    #import pdb
    #pdb.set_trace()
    ax.set_xlabel('Leading AK8 jet softdrop mass [GeV]', ha='right', x=1)
    ax.set_ylabel(f'Events / {rebin_factor} GeV', ha='right', y=1)
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

  corrList = ['no_PUPPI', 'no_JMS_JMR', 'no_JMS', 'no_PUPPI_JMS_JMR', 'no_JMR']
  hist, col, ls = prepare_removeCorrections(indir, histName, corrList)
  plot_removeCorrections(hist, col, ls, outdir) 

#  for unc in ['jer','jesTotal','jmr','jms', 'puWeight']:
#      plot_unc(indir,histName, unc, outdir)
#      #plot_unc(indir,histName, args.uncertainty, outdir)

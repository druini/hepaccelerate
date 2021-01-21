import numpy as np
import ROOT as r
import matplotlib.pyplot as plt
import mplhep as hep
import os,argparse
from pdb import set_trace

def getBias(rootFile, injr, rMin, rMax):
    print(rootFile)
    if not os.path.isfile(rootFile): return np.nan, np.nan

    lFile = r.TFile(rootFile)
    lTree = lFile.Get("tree_fit_sb")
    
    lH   = r.TH1D('h_bias','h_bias',50,-4,4)
    lH_1 = r.TH1D('h_bias_1','h_bias',50,-4,4)
    lH_2 = r.TH1D('h_bias_2','h_bias',50,-4,4)
    lTree.Project('h_bias_1',f'(r-{injr})/rLoErr', f'r>={injr}&&(rHiErr+r)<{rMax-1}&&(r-rLoErr)>{rMin+1}')
    lTree.Project('h_bias_2',f'(r-{injr})/rHiErr', f'r<{injr}&&(rHiErr+r)<{rMax-1}&&(r-rLoErr)>{rMin+1}')
    lH = lH_1
    lH.Add(lH_2)
    print(f'Tree Entries = {lTree.GetEntriesFast()}, pull entries = {lH.GetEntries()}')

    gaus_func = r.TF1("gaus_func","gaus(0)",-3,3)
    gaus_func.SetParameter(0,20)
    gaus_func.SetParameter(1,0)
    gaus_func.SetParameter(2,1)
    lH.Fit(gaus_func,"mlerqn")

    mean = round(gaus_func.GetParameter(1), 2)
    std  = round(gaus_func.GetParameter(2), 2)

    return mean, std

def plotMatrix(biasMatrix, funcs, funcDegs, outdir):
    hep.set_style(hep.style.CMS)
    fig, ax = plt.subplots()
    ax.set_xticklabels(funcs)
    ax.set_yticklabels(funcs)
    hep.hist2dplot(biasMatrix, labels=True, vmin=-1, vmax=1, cmap='RdBu_r')
    ax.set_xlabel('gen function', ha='right', x=1)
    ax.set_ylabel('fit function', ha='right', y=1)
    #ax.set_title('MC')
    if not os.path.isdir(outdir): os.makedirs(outdir)
    for ext in ['png','pdf']:
        fig.savefig(os.path.join(outdir,f'{args.dataOrMC}_'+'_'.join([f'{f}{d}' for f,d in funcDegs.items()])+f'_r{args.rInjected}_{args.signal}.{ext}'))
    #plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rootDir', default='output')
    parser.add_argument('-y', '--year', default='2018')
    parser.add_argument('-v', '--version', default='v14')
    parser.add_argument('-s', '--selection', default='met20_btagDDBvL086')
    parser.add_argument('-d', '--dataOrMC', default='data')
    parser.add_argument('--msd_start', default=100)
    parser.add_argument('--msd_stop', default=160)
    parser.add_argument('--msd_binSize', default=5)
    parser.add_argument('--signal', default='ttHTobb')
    parser.add_argument('--rMin', default=-20)
    parser.add_argument('--rMax', default=20)
    parser.add_argument('-r', '--rInjected', default=0)
    parser.add_argument('--degBern', default=3)
    parser.add_argument('--degCheb', default=3)
    parser.add_argument('--degpoly', default=4)

    args = parser.parse_args()

    funcNames = {
            'Bern' : 'Bern',
            'Cheb' : 'Cheb',
            'poly' : 'P',
            #'logN' : 'logN',
            }

    funcDegs = {
            'Bern' : args.degBern,
            'Cheb' : args.degCheb,
            'poly' : args.degpoly
            #'logN' : 'logN',
            }

    baseDir      = os.path.join( args.rootDir, args.year, args.version, args.selection)
    funcDirStart = f'{args.dataOrMC}_msd{args.msd_start}to{args.msd_stop}_msdbin{args.msd_binSize}'
    bkgDir       = f'bkgEstTests_r{args.rMin}to{args.rMax}'

    files = { f1 : { f2 : os.path.join( baseDir, f'{funcDirStart}_{f1}polyDegs{d1}_{args.signal}', bkgDir, f'biastoys_bias_{f1}1{d1}_vs_{f2}1{d2}_r{args.rInjected}_merged.root') for f2,d2 in funcDegs.items() } for f1,d1 in funcDegs.items() }
    biases = { f1+f2 : getBias(files[f1][f2],args.rInjected,args.rMin,args.rMax) for f1 in funcDegs for f2 in funcDegs }

    biasMatrix = np.array( [[biases[f2+f1][0] for f2 in funcDegs] for f1 in funcDegs] )

    ### now the plotting
    outdir = os.path.join( baseDir, 'plotsBias', funcDirStart)
    plotMatrix(biasMatrix, [n+str(d) for n,d in zip(funcNames.values(),funcDegs.values())], funcDegs, outdir)
    edges = range(4)

#biasMatrixMC = np.array([
#    [-.02,  .04, .69],
#    [2.13, -.09, .5],
#    [-.61, 4.27, -.03]
#    ])
#
#biasMatrixData = np.array([
#    [-.03, .08, .04],
#    [2.78, -.05, .01],
#    [3.02, 2.5, -.04]
#    ])
#
#fig, ax = plt.subplots()
#ax.set_xticklabels(funcs)
#ax.set_yticklabels(funcs)
#hep.hist2dplot(biasMatrixData, edges, edges, labels=True, vmin=-1, vmax=1, cmap='RdBu_r')
#ax.set_xlabel('gen function', ha='right', x=1)
#ax.set_ylabel('fit function', ha='right', y=1)
#ax.set_title('data')
#plt.show()

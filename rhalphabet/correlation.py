import ROOT
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import mplhep as hep
import pylab
import argparse
import numpy as np
import os
from makePlots import makeOutdir
from pdb import set_trace

#params = {'legend.fontsize': 'x-large', 'axes.labelsize': 'x-large', 'axes.titlesize':'x-large', 'xtick.labelsize':'x-large', 'ytick.labelsize':'x-large'}
#pylab.rcParams.update(params)

def plot(input, pars, outputname):

    # Load file
    file = ROOT.TFile(input)
    fit = file.Get("fit_s")

    # Read in correlations
    out = {"pars": pars}
    out["corrs"] = []
    
    for i in pars:
        for j in pars:
            out["corrs"] += [(i, j, fit.correlation(i, j))]
            print(i, j, fit.correlation(i, j))
            #if i == "r_ttH" and j == "r_ttbb":
            #    print fit.correlation(i,j)

    #print out["corrs"]

    # Build correlation matrix
    mat = np.zeros((len(out["pars"]), len(out["pars"])))
    print(len(out["pars"]))
    
    for i, ip in enumerate(out["pars"]):
        for j, jp in enumerate(out["pars"]):
            mat[i, j] = round(out["corrs"][i + len(out["pars"])*j][2], 2)

    # Plot correlation matrix
    #fig = plt.figure(figsize=(30,30))
    hep.set_style(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(30,30))
    fig.subplots_adjust(right=.85)
    ax = hep.cms.label(data=args.isData, paper=False, year=args.year, ax=ax, fontsize=60)
    pl = hep.hist2dplot(mat, labels=True, cmap = "seismic", vmin =-1, vmax =1)
    pl.cbar.ax.tick_params(labelsize=40)
    #m = plt.imshow(mat, cmap = "seismic", vmin =-1, vmax =1,interpolation='none')
    #for i in range(mat.shape[0]):
    #    for j in range(mat.shape[1]):
    #        plt.text(i, j, "{0: .2f}".format(mat[i,j]), color = "orange", fontsize = 11, ha = "center", va = "center")
    out["pars"][0] = r"$\mu_{ttH}$"
    plt.xticks(.5+np.array(range(len(out["pars"]))), out["pars"], rotation = 90, fontsize = 30)
    plt.yticks(.5+np.array(range(len(out["pars"]))), out["pars"], fontsize = 30)
    #plt.text(-0.5, -1,
    #    r"$\mathbf{CMS}$ Work In Progress",
    #    fontsize=38, ha="left", va="top", fontname="Helvetica"
    #)
    #plt.text(14.5, -1,
    #    "$41.5\ \mathrm{fb}^{-1}\ \mathrm{(13\ TeV)}$",
    #    fontsize=28, ha="right", va="top", fontname="Helvetica"
    #)
    #plt.gcf().subplots_adjust(left=0.3)
    #plt.gcf().subplots_adjust(bottom=0.3)
    #cbaxes = fig.add_axes([0.92, 0.3, 0.03, 0.6]) 
    #cbar = plt.colorbar(m, cax = cbaxes)   
    #cbar.set_label("post-fit correlation", size=18)
    plt.savefig(outputname+'.pdf', bbox_inches='tight')


if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(description = 'Plot correlation matrix')
    
    parser.add_argument('--input', '-i', help = '.root file containing correlation matrix', required = True)
    parser.add_argument('--outputname', '-o', default='covariance', help = 'name of output plot')
    parser.add_argument('--year', '-y')
    parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
    
    args = parser.parse_args()

    outdir = makeOutdir(args.input)

    #Normal

    pars = [
            'lumi_13TeV_2017',
            'CMS_ttHbb_AK4deepjetM',
            'CMS_ttHbb_AK8DDBvLM1',
            'CMS_ttHbb_jer',
            'CMS_ttHbb_jesAbsolute',
            f'CMS_ttHbb_jesAbsolute_{args.year}',
            'CMS_ttHbb_jesBBEC1',
            f'CMS_ttHbb_jesBBEC1_{args.year}',
            'CMS_ttHbb_jesEC2',
            f'CMS_ttHbb_jesEC2_{args.year}',
            'CMS_ttHbb_jesFlavorQCD',
            'CMS_ttHbb_jesHF',
            f'CMS_ttHbb_jesHF_{args.year}',
            'CMS_ttHbb_jesRelativeBal',
            f'CMS_ttHbb_jesRelativeSample_{args.year}',
            'CMS_ttHbb_jmr',
            'CMS_ttHbb_jms',
            'CMS_ttHbb_pdfWeight',
            'CMS_ttHbb_psWeight_FSR',
            'CMS_ttHbb_psWeight_ISR',
            'CMS_ttHbb_puWeight',
            f'boosted_bkg_paramX0_{args.year}',
            f'boosted_bkg_paramX1_{args.year}',
            f'boosted_bkg_paramX2_{args.year}',
            ]
    if args.year=='2017':
        pars += [ f'boosted_bkg_paramX3_{args.year}']
    plot(args.input, pars, os.path.join(outdir,args.outputname))
    print(len(pars))
    

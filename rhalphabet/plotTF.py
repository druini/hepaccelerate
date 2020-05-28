###########################################
##### Modified code taken from https://github.com/andrzejnovak/rhalphalib/tree/hxxdev
##### How to run: python plotTF.py -f ttH_met20_deepTagMD_bbvsLight08695/mc_msd90to160_msdbin5_pt2bin_polyDegs22/fitDiagnostics.root -o ttH_met20_deepTagMD_bbvsLight08695/mc_msd90to160_msdbin5_pt2bin_polyDegs22/ --isData --msd_start 90 --msd_stop 160
###########################################

import argparse
import os
from operator import methodcaller

import uproot
import ROOT as r

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import mplhep as hep

plt.switch_backend('agg')


# Benstein polynomial calculation
def bern_elem(x, v, n):
    # Bernstein element calculation
    normalization = 1. * math.factorial(n) / (math.factorial(v) * math.factorial(n - v))
    Bvn = normalization * (x**v) * (1-x)**(n-v)
    return float(Bvn)


def TF(pT, rho, n_pT=2, n_rho=2, par_map=np.ones((3, 3))):
    # Calculate TF Polynomial for (n_pT, n_rho) degree Bernstein poly
    val = 0
    for i_pT in range(0, n_pT+1):
        for i_rho in range(0, n_rho+1):
            val += (bern_elem(pT, i_pT, n_pT)
                    * bern_elem(rho, i_rho, n_rho)
                    * par_map[i_pT][i_rho])

    return val



###################################################################
# TF Plots
def plotTFsmooth(TF, msd, pt, mask=None, MC=False, raw=False, rhodeg=2, ptdeg=2, out=None,
           year="2017", label=None):
    """
    Parameters:
    TF: Transfer Factor array
    msd: Mass bins array (meshgrid-like)
    pt: pT bins array (meshgrid-like)
    """
    import matplotlib.pyplot as plt
    #import matplotlib.font_manager
    import mplhep as hep
    plt.style.use([hep.style.ROOT, {'font.size': 24}])
    plt.switch_backend('agg')
    fig, ax = plt.subplots()
    if mask is not None:
        TF = np.ma.array(TF, mask=~mask)

    zmin, zmax = np.floor(10*np.min(TF))/10, np.ceil(10*np.max(TF))/10
    zmin, zmax = zmin + 0.001, zmax - 0.001
    clim = np.round(np.min([abs(zmin - 1), abs(zmax - 1)]), 1)
    if clim < .3: clim = .3
    if clim > .5: clm = .5
    levels = np.linspace(1-clim, 1+clim, 500)

    if np.min(TF) < 1-clim and np.max(TF) > 1+clim: _extend = 'both '
    elif np.max(TF) > 1+clim: _extend = 'max'
    elif np.min(TF) < 1-clim: _extend = 'min'
    else: _extend = 'neither'

    if mask is not None:
        contf = ax.contourf(msd, pt, TF, levels=levels,
                            corner_mask=False, cmap='RdBu_r', extend=_extend)
    else:
        contf = ax.contourf(msd, pt, TF, levels=levels, cmap='RdBu_r', extend=_extend)
    cax = hep.make_square_add_cbar(ax, pad=0.2, size=0.5)
    if abs(1-zmin) > .3 and abs(1-zmax) > .3:
        c_extend = 'both'
    elif abs(1-zmin) > .3:
        c_extend = 'min'
    elif abs(1-zmax) > .3:
        c_extend = 'max'
    else:
        c_extend = 'neither'
    cbar = fig.colorbar(contf, cax=cax, extend=c_extend)
    cbar.set_ticks([np.arange(1-clim, 1+clim, 0.1)])

    def rho_bound(ms, rho):
        # rho = {minRho, maxRho}
        fpt = ms * np.e**(-rho/2)
        return fpt

    ### not needed for our mass range (this is the region rho<-6)
    #x = np.arange(40, 70)
    #ax.plot(x, rho_bound(x, minRho), 'black', lw=3)
    #ax.fill_between(x, rho_bound(x, minRho), 1200, facecolor="none", hatch="xx",
    #                edgecolor="black", linewidth=0.0)
    x = np.arange(50, 201)
    ax.plot(x, rho_bound(x, maxRho) + 5, 'black', lw=3)
    ax.fill_between(x, rho_bound(x, maxRho), facecolor="none", hatch="xx",
                    edgecolor="black", linewidth=0.0)

    _mbins, _pbins = np.array( range(minMass, maxMass, 5 )), np.array([minPt, 350, maxPt])
    #_mbins, _pbins = np.linspace(40, 201, 24), np.array([450, 500, 550, 600, 675, 800, 1200])
    sampling = np.meshgrid(_mbins[:-1] + 0.5 * np.diff(_mbins), _pbins[:-1] + 0.3 * np.diff(_pbins))
    valmask = (sampling[1] > rho_bound(sampling[0], maxRho)) & (sampling[1] < rho_bound(sampling[0], minRho))
    ax.scatter(sampling[0][valmask], sampling[1][valmask], marker='x', color='black', s=40, alpha=.4, )

    ax.set_xlim(minMass, maxMass)
    ax.set_ylim(minPt, maxPt)
    ###ax.invert_yaxis()

    tMC = "Tagger Response" if MC else "Data Residual"
    if raw:
        tMC = "Tagger Response (prefit)"
    if label is None:
        label = '{} TF({},{})'.format(tMC, rhodeg, ptdeg)
    ax.set_title(label,
                 pad=9,
                 fontsize=22,
                 loc='left')
    ax.set_title("({})".format(str(year)),
                 pad=9,
                 fontsize=22,
                 loc='right')
    ax.set_xlabel(r'Jet $\mathrm{m_{SD}}$', ha='right', x=1)
    ax.set_ylabel(r'Jet $\mathrm{p_{T}}$', ha='right', y=1)
    cbar.set_label(r'TF', ha='right', y=1)

    label = "MC" if MC else "Data"
    if raw:
        label = "MCRaw"
    import mplhep as hep
    hep.cms.cmslabel(loc=2, data=not raw, rlabel="", ax=ax)
    if out is not None:
      for ext in ['png','pdf']:
        fig.savefig(os.path.join(out,f'TF.{ext}'))#, bbox_inches="tight")
    else:
        return fig

######################################################################
def TF_params(xparlist, xparnames=None, nrho=None, npt=None):
    '''Returns a 2D array of parameters from the fit'''

    # TF map from param/name lists
    if xparnames is not None:
        from operator import methodcaller

        def _get(s):
            return s[-1][0]

        ptdeg = max(
            list(
                map(
                    int,
                    list(
                        map(_get, list(map(methodcaller("split", 'pt_par'),
                                           xparnames)))))))
        rhodeg = max(
            list(
                map(
                    int,
                    list(
                        map(_get, list(map(methodcaller("split", 'rho_par'),
                                           xparnames)))))))
    else:
        rhodeg, ptdeg = nrho, npt

    TF_cf_map = np.array(xparlist).reshape(rhodeg + 1, ptdeg + 1)

    return TF_cf_map, rhodeg, ptdeg


######################################################################
def TF_smooth_plot(_tfmap, _rhodeg, _ptdeg, minPt, maxPt, minMsd, maxMsd):
    # Define fine bins for smooth TF plots
    fptbins = np.arange(minPt, maxPt, 2)
    fmsdbins = np.arange(minMsd, maxMsd, .5)
    fptpts, fmsdpts = np.meshgrid(fptbins[:-1] + 0.3 * np.diff(fptbins),
                                  fmsdbins[:-1] + 0.5 * np.diff(fmsdbins),
                                  indexing='ij')
    frhopts = 2*np.log(fmsdpts/fptpts)
    fptscaled = (fptpts - minPt) / (maxPt - minPt)
    frhoscaled = (frhopts - (minRho)) / ((maxRho) - (minRho))
    fvalidbins = (frhoscaled >= 0) & (frhoscaled <= 1)
    frhoscaled[~fvalidbins] = 1  # we will mask these out later

    def wrapTF(pT, rho):
        return TF(pT, rho, n_pT=_ptdeg, n_rho=_rhodeg, par_map=_tfmap)

    TFres = np.array(list(map(wrapTF, fptscaled.flatten(), frhoscaled.flatten()))).reshape(fptpts.shape)

    # return TF, msd bins, pt bins, mask
    return TFres, fmsdpts, fptpts, fvalidbins


############################################################
## not used now, maybe in future
def plotTF_ratio(in_ratio, mask, region):
    fig, ax = plt.subplots()

    H = np.ma.masked_where(in_ratio * mask <= 0.01, in_ratio * mask)
    zmin, zmax = np.floor(10*np.min(TFres))/10, np.ceil(10*np.max(TFres))/10
    zmin, zmax = zmin + 0.001, zmax - 0.001
    clim = np.max([.3, np.min([abs(zmin - 1), abs(zmax - 1)])])
    ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    msdbins = np.linspace(40, 201, 24)
    hep.hist2dplot(H.T, msdbins, ptbins, vmin=1-clim, vmax=1+clim,
                   cmap='RdBu_r', cbar=False)
    cax = hep.make_square_add_cbar(ax, pad=0.2, size=0.5)
    if abs(1-zmin) > .3 and abs(1-zmax) > .3:
        c_extend = 'both'
    elif abs(1-zmin) > .3:
        c_extend = 'min'
    elif abs(1-zmax) > .3:
        c_extend = 'max'
    else:
        c_extend = 'neither'
    cbar = fig.colorbar(ax.get_children()[0], cax=cax, extend=c_extend)

    ax.set_xticks(np.arange(40, 220, 20))
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    ax.invert_yaxis()

    ax.set_title('{} QCD Ratio'.format(region), pad=15, fontsize=26)
    ax.set_xlabel(r'Jet $\mathrm{m_{SD}}$', ha='right', x=1)
    ax.set_ylabel(r'Jet $\mathrm{p_{T}}$', ha='right', y=1)
    cbar.set_label(r'(Pass QCD) / (Fail QCD * eff)', ha='right', y=1)

    fig.savefig('{}/{}{}.png'.format(args.output_folder, "TF_ratio_", region),
                bbox_inches="tight")




#############################################################################
if __name__ == '__main__':

    plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
    plt.switch_backend('agg')
    np.seterr(divide='ignore', invalid='ignore')

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fit",
                        default='fitDiagnostics.root',
                        dest='fit',
                        help="fitDiagnostics file, example: RhalphabetResults/fitDiagnostics.root")

    parser.add_argument("-o", "--output-folder",
                        default=None,
                        dest='output_folder',
                        help="Folder to store plots - will be created ? doesn't exist.")

    parser.add_argument("-y", "--year",
                        default="2017",
                        type=str,
                        help="year label")

    parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
    parser.add_argument('--msd_start', default=90, type=int, help='start of the mass range')
    parser.add_argument('--msd_stop', default=170, type=int, help='stop of the mass range')

    args = parser.parse_args()

    if not args.fit.endswith('.root'):
      raise Exception(f'Must supply ROOT filename, but got {args.fit}')
    if args.output_folder is None:
      if not (os.path.sep in args.fit): #if the file is in the current directory
        outdir = os.getcwd()
      else:
        outdir = os.path.join(*args.fit.split(os.path.sep)[:-1])
      outdir = os.path.join(outdir,'plots')
    else:
      outdir = args.output_folder
    if not os.path.exists(outdir):
      os.makedirs(outdir)

    minPt = 250
    maxPt = 1500
    minMass = args.msd_start
    maxMass = args.msd_stop
    minRho = -6.
    maxRho = -1.2

    # Get fitDiagnostics File
    rf = r.TFile.Open(args.fit)

    # Get TF parameters
    hmp = []
    par_names = rf.Get('fit_s').floatParsFinal().contentsString().split(',')
    par_names = [p for p in par_names if 'tf' in p]
    MCTF = []
    for pn in par_names:
        if "deco" not in pn:
            hmp.append(round(rf.Get('fit_s').floatParsFinal().find(pn).getVal(), 4))
        elif "deco" in pn:
            MCTF.append(round(rf.Get('fit_s').floatParsFinal().find(pn).getVal(), 4))

    def _get(s):
        # Helper
        return s[-1][0]

    par_names = [n for n in par_names if "deco" not in n]
    ptdeg = max(
        list(
            map(int, list(map(_get, list(map(methodcaller("split", 'pt_par'),
                                             par_names)))))))
    rhodeg = max(
        list(
            map(int, list(map(_get, list(map(methodcaller("split", 'rho_par'),
                                             par_names)))))))

    parmap = np.array(hmp).reshape(rhodeg+1, ptdeg+1)
    if len(MCTF) > 0:
        MCTF_map = np.array(MCTF).reshape(rhodeg+1, ptdeg+1)\

    ##### Smooth plots
    _values = hmp
    # TF Data
    plotTFsmooth(*TF_smooth_plot(*TF_params(_values, nrho=rhodeg, npt=ptdeg), minPt, maxPt, minMass, maxMass), MC=False,
                 out=outdir,
                 year=args.year,
                 label=('Data' if args.isData else 'MC')+' TF({},{})'.format(rhodeg, ptdeg))
                 #out=args.output_folder+'/TF_'+('Data' if args.isData else 'MC')+'_ptDeg'+str(ptdeg)+'rhoDeg'+str(rhodeg)+'_'+args.year, year=args.year,

#    # Effective TF (combination)  #### ALE: this is maybe needed later
#    _tf1, _, _, _ = TF_smooth_plot(*TF_params(hmp, nrho=2, npt=2), minPt, maxPt, minMass, maxMass)
#    _tf2, bit1, bit2, bit3 = TF_smooth_plot(*TF_params(_values, nrho=2, npt=2), minPt, maxPt, minMass, maxMass)
#    plotTFsmooth(_tf1*_tf2, bit1, bit2, bit3, MC=True, raw=False,
#                 out='{}/TF_eff'.format(args.output_folder), year=args.year,
#                 label='Effective TF({},{})'.format(2, 2))

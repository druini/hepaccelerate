from __future__ import print_function
import sys, os, argparse
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import ROOT
import json
rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def exec_me(command, dryRun=False, folder=False):
    print(command)
    if not dryRun:
        if folder: os.chdir(folder)
        os.system(command)

#ver          = 'v04'
#msd_start    = 95
#msd_stop     = 150
#polyDeg      = 2
#rebin_factor = 5
#ptbins = np.array([250,300,5000])
#ptbins = np.array([250,300,450,5000])

def load_from_json(indir, sample, ptStart, ptStop, msd_start_idx, msd_stop_idx, region, rebin_factor, obs):
  #filepath = '/afs/cern.ch/work/d/druini/public/hepaccelerate/results/2018/v05/'+ver+'/out_'+sample+'_merged.json'
  filepath = os.path.join(indir, 'out_'+sample+'_merged.json')
  with open(filepath) as json_file:
    data = json.load(json_file)
    data = data['hist_leadAK8JetMass_2J2WdeltaR_'+region+'_pt%sto%s' % (ptStart, ptStop)]
    rebin(data,rebin_factor)
  assert( np.all(np.array(obs.binning)==np.array(data['edges'])[msd_start_idx:msd_stop_idx+1]) )
  return(np.array(data['contents'])[msd_start_idx:msd_stop_idx], obs.binning, obs.name)

#def rebin(bins, counts, yerr, rebin_factor):
def rebin(hist, rebin_factor):
    bins   = np.array(hist['edges'])#[:-1])
    counts = np.array(hist['contents'])

    new_bins   = bins[::rebin_factor]
    new_counts = np.add.reduceat(counts, range(0, len(counts), rebin_factor))
#    new_yerr   = np.add.reduceat(yerr, range(0, len(yerr), rebin_factor))
#    return new_bins, new_counts#, new_yerr
    hist['edges']    = new_bins
    hist['contents'] = new_counts


def test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,isData=True,runBias=False):
    dataOrBkg = 'data' if isData else 'background'

    throwPoisson = False

    jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    massScale = rl.NuisanceParameter('CMS_msdScale', 'shape')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    tqqeffSF = rl.IndependentParameter('tqqeffSF', 1., 0, 10)
    tqqnormSF = rl.IndependentParameter('tqqnormSF', 1., 0, 10)

    #ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    #ptbins = np.append( np.arange(250,600,50), [600, 675, 800, 1200] )
    #ptbins = np.array([0,1500])
    npt = len(ptbins) - 1
    msdbins = np.linspace(0,300,301)[::rebin_factor]
    msd_start_idx = np.where(msdbins==msd_start)[0][0]
    msd_stop_idx  = np.where(msdbins==msd_stop)[0][0]
    msdbins = msdbins[msd_start_idx:msd_stop_idx+1]
    msd = rl.Observable('msd', msdbins)

    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(ptbins[:-1] + 0.3 * np.diff(ptbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    rhopts = 2*np.log(msdpts/ptpts)
    ptscaled = (ptpts - ptbins[0]) / (ptbins[-1] - ptbins[0])
    #import pdb
    #pdb.set_trace()
    rho_start = -6
    rho_stop  = -1.2
    rhoscaled = (rhopts - rho_start) / (rho_stop - rho_start)
    validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
    rhoscaled[~validbins] = 1  # we will mask these out later

    # Build bkg MC pass+fail model and fit to polynomial
    ### This model is for the prefit but helps to compute the ratio pass/fail
    bkgmodel = rl.Model("bkgmodel")
    bkgpass, bkgfail = 0., 0.
    for ptbin in range(npt):
        failCh = rl.Channel("ptbin%d%s" % (ptbin, 'fail'))
        passCh = rl.Channel("ptbin%d%s" % (ptbin, 'pass'))
        bkgmodel.addChannel(failCh)
        bkgmodel.addChannel(passCh)
        # mock template
        ptnorm = 1
#        failTempl = expo_sample(norm=ptnorm*1e5, scale=40, obs=msd)
#        passTempl = expo_sample(norm=ptnorm*1e3, scale=40, obs=msd)
        failTempl = load_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Fail', rebin_factor, msd)
        passTempl = load_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Pass', rebin_factor, msd)
        failCh.setObservation(failTempl)
        passCh.setObservation(passTempl)
        bkgfail += failCh.getObservation().sum()
        bkgpass += passCh.getObservation().sum()

    bkgeff = bkgpass / bkgfail
    if args.runPrefit:
        tf_MCtempl = rl.BernsteinPoly("tf_MCtempl", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-10, 10))
        tf_MCtempl_params = bkgeff * tf_MCtempl(ptscaled, rhoscaled)
        for ptbin in range(npt):
            failCh = bkgmodel['ptbin%dfail' % ptbin]
            passCh = bkgmodel['ptbin%dpass' % ptbin]
            failObs = failCh.getObservation()
            bkgparams = np.array([rl.IndependentParameter('bkgparam_ptbin%d_msdbin%d' % (ptbin, i), 0) for i in range(msd.nbins)])
            sigmascale = 10.
            scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**bkgparams
            fail_bkg = rl.ParametericSample('ptbin%dfail_bkg' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_bkg)
            pass_bkg = rl.TransferFactorSample('ptbin%dpass_bkg' % ptbin, rl.Sample.BACKGROUND, tf_MCtempl_params[ptbin, :], fail_bkg)
            passCh.addSample(pass_bkg)

            failCh.mask = validbins[ptbin]
            passCh.mask = validbins[ptbin]

        bkgfit_ws = ROOT.RooWorkspace('bkgfit_ws')
        simpdf, obs = bkgmodel.renderRoofit(bkgfit_ws)
        bkgfit = simpdf.fitTo(obs,
                              ROOT.RooFit.Extended(True),
                              ROOT.RooFit.SumW2Error(True),
                              ROOT.RooFit.Strategy(2),
                              ROOT.RooFit.Save(),
                              ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                              ROOT.RooFit.PrintLevel(-1),
                              )
        bkgfit_ws.add(bkgfit)
        if "pytest" not in sys.modules:
             bkgfit_ws.writeToFile(os.path.join(str(outdir), 'ttHbb_bkgfit.root'))
        if bkgfit.status() != 0:
            raise RuntimeError('Could not fit bkg')

        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', bkgfit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

    #### Actual transfer function for combine
    tf_dataResidual = rl.BernsteinPoly("tf_dataResidual", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-10, 10), coefficient_transform=(np.exp if runBias else None))
    tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
    #tf_params = bkgeff * (tf_MCtempl_params_final * tf_dataResidual_params if args.runPrefit else tf_dataResidual_params)
    tf_params = bkgeff * tf_dataResidual_params
    if args.runPrefit:
        tf_params *= tf_MCtempl_params_final

    # build actual fit model now
    model = rl.Model("ttHbb")

    for ptbin in range(npt):
        for region in ['Pass', 'Fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            templates = {
                'signal'     : load_from_json(indir, 'signal', ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor, msd),
                'background' : load_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor,  msd),
            }
            # some mock expectations
            templ = templates['signal']
            stype = rl.Sample.SIGNAL
            sample = rl.TemplateSample(ch.name + '_signal', stype, templ)

#            # mock systematics
#            jecup_ratio = np.random.normal(loc=1, scale=0.05, size=msd.nbins)
#            msdUp = np.linspace(0.9, 1.1, msd.nbins)
#            msdDn = np.linspace(1.2, 0.8, msd.nbins)
#
#            # for jec we set lnN prior, shape will automatically be converted to norm systematic
#            sample.setParamEffect(jec, jecup_ratio)
#            sample.setParamEffect(massScale, msdUp, msdDn)
            sample.setParamEffect(lumi, 1.027)      ### for bias/ftest/GOF signal needs at least one unc.

            ch.addSample(sample)

            #import pdb
            #pdb.set_trace()
            #yields = sum(tpl[0] for tpl in templates.values())
            #data_obs = (yields, msd.binning, msd.name)
            data_obs = templates['background']
            ch.setObservation(data_obs)

            # drop bins outside rho validity
            mask = validbins[ptbin]
            # blind bins 11, 12, 13
            # mask[11:14] = False
            if isData:
              msdWindow_start_idx = np.where(msd.binning==110)[0][0]
              msdWindow_stop_idx  = np.where(msd.binning==140)[0][0]
              mask[msdWindow_start_idx:msdWindow_stop_idx] = False
            ch.mask = mask

    for ptbin in range(npt):
        failCh = model['ptbin%dFail' % ptbin]
        passCh = model['ptbin%dPass' % ptbin]

        bkgparams = np.array([rl.IndependentParameter('bkgparam_ptbin%d_msdbin%d' % (ptbin, i), 0) for i in range(msd.nbins)])
        initial_bkg = failCh.getObservation().astype(float)  # was integer, and numpy complained about subtracting float from it
        for sample in failCh:
            initial_bkg -= sample.getExpectation(nominal=True)
        if np.any(initial_bkg < 0.):
            raise ValueError("initial_bkg negative for some bins..", initial_bkg)
        sigmascale = 10  # to scale the deviation from initial
        scaledparams = initial_bkg * (1 + sigmascale/np.maximum(1., np.sqrt(initial_bkg)))**bkgparams
        fail_bkg = rl.ParametericSample('ptbin%dFail_bkg' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams)
        failCh.addSample(fail_bkg)
        pass_bkg = rl.TransferFactorSample('ptbin%dPass_bkg' % ptbin, rl.Sample.BACKGROUND, tf_params[ptbin, :], fail_bkg)
        passCh.addSample(pass_bkg)

        #tqqpass = passCh['tqq']
        #tqqfail = failCh['tqq']
        #tqqPF = tqqpass.getExpectation(nominal=True).sum() / tqqfail.getExpectation(nominal=True).sum()
        #tqqpass.setParamEffect(tqqeffSF, 1*tqqeffSF)
        #tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
        #tqqpass.setParamEffect(tqqnormSF, 1*tqqnormSF)
        #tqqfail.setParamEffect(tqqnormSF, 1*tqqnormSF)

    #with open(os.path.join(str(outdir), 'ttHbb.pkl'), "wb") as fout:       ### ALE: still dont understand why we need this
    #    pickle.dump(model, fout)

    pref = ('biasTest_' if runBias else '') + ('data' if isData else 'mc')
    combineFolder = os.path.join(str(outdir), pref+'_msd%dto%d_msdbin%d_pt%dbin_polyDegs%d%d'%(msd_start,msd_stop,rebin_factor,len(ptbins)-1, polyDegPt,polyDegRho))
    model.renderCombine(combineFolder)
    exec_me('bash build.sh | combine -M FitDiagnostics ttHbb_combined.txt  --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace --setParameterRanges r=-1,1', folder=combineFolder)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
  parser.add_argument('-f', '--runPrefit', action='store_true', help='Run prefit on MC.' )
  parser.add_argument('-b', '--runBias', action='store_true', help='Run prefit on MC.' )
  parser.add_argument('-smr', '--scanMassRange', action='store_true', default=False, help='option to scan mass range')
  parser.add_argument('-spd', '--scanPolyDeg', action='store_true', default=False, help='option to scan degree of polynomial')
  parser.add_argument('--msd_start', default=90, type=int, help='start of the mass range')
  parser.add_argument('--msd_stop', default=170, type=int, help='stop of the mass range')
  parser.add_argument('--polyDegPt', default=2, type=int, help='degree of polynomial to fit pt')
  parser.add_argument('--polyDegRho', default=2, type=int, help='degree of polynomial to fit rho')
  parser.add_argument('-r', '--rebin_factor', default=5, type=int, help='rebin factor for json histograms, default mass bin size is 1GeV')
  parser.add_argument('--nptbins', default=2, type=int, help='number of pt bins')
  parser.add_argument('-y', '--year', default='2017', type=str, help='year to process, in file paths')
  parser.add_argument('-v', '--version', default='v05', help='version, in file paths')
  parser.add_argument('-s', '--selection', nargs='+', default=['met20_btagDDBvL_noMD07','met20_deepTagMD_bbvsLight05845','met20_deepTagMD_bbvsLight08695'], help='event selection, in file paths')
  parser.add_argument('-j', '--jsonpath', default='/afs/cern.ch/work/d/druini/public/hepaccelerate/results', help='path to json files')

  try: args = parser.parse_args()
  except:
    parser.print_help()
    sys.exit(0)

  if args.nptbins==1:
    ptbins = np.array([250,5000])
  elif args.nptbins==2:
    ptbins = np.array([250,300,5000])
  elif args.nptbins==3:
    ptbins = np.array([250,300,450,5000])
  else:
    raise Exception('invalid number of ptbins')

  msd_start    = args.msd_start
  msd_stop     = args.msd_stop
  polyDegPt    = args.polyDegPt
  polyDegRho   = args.polyDegRho
  rebin_factor = args.rebin_factor

  for sel in args.selection:
    indir = os.path.join(args.jsonpath, args.year, args.version, sel)
    if not os.path.exists(indir):
      raise Exception('invalid input path: %s'%indir)

    outdir = os.path.join('output', args.year, args.version, sel)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not (args.scanMassRange or args.scanPolyDeg):
      test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,args.isData,runBias=args.runBias)
    else:
      pref = 'data' if args.isData else 'mc'
      if args.scanMassRange:
        for msd_start in range(70,105,rebin_factor):
          for msd_stop in range(150,251,rebin_factor):
            with open(os.path.join(outdir,'scan_'+pref+'_fitRange_%dptbin_msdbin%dgev_polyDegs%d%d.txt'%(len(ptbins)-1, rebin_factor, polyDegPt,polyDegRho)),'a') as f:
              try:
                print('trying %d to %d'%(msd_start,msd_stop))
                test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,args.isData)
                f.write('msd range %d to %d worked!!\n' %(msd_start, msd_stop))
              except:
                f.write('msd range %d to %d failed\n' %(msd_start, msd_stop))

      if args.scanPolyDeg:
        for polyDegPt in range(0,5):
          for polyDegRho in range(0,10):
            with open(os.path.join(outdir,'scan_'+pref+'_polyDeg_%dptbin_msdbin%dgev_fitRange%dto%d.txt'%(len(ptbins)-1, rebin_factor, msd_start, msd_stop)),'a') as f:
              try:
                print('trying ptdeg %d and rhodeg %d'%(polyDegPt,polyDegRho))
                test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,args.isData)
                f.write('polyDegPt %d and polyDegRho %d worked!!\n'%(polyDegPt,polyDegRho))
              except:
                f.write('polyDegPt %d and polyDegRho %d failed\n'%(polyDegPt,polyDegRho))

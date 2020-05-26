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

#ver          = 'v04'
#msd_start    = 95
#msd_stop     = 150
#polyDeg      = 2
#rebin_factor = 5
#ptbins = np.array([250,300,5000])
#ptbins = np.array([250,300,450,5000])

def load_from_json(sample, ptStart, ptStop, msd_start_idx, msd_stop_idx, region, rebin_factor, obs, ver):
  filepath = '/afs/cern.ch/work/d/druini/public/hepaccelerate/results/2017/v05/'+ver+'/out_'+sample+'_merged.json'
  #filepath = 'json_histos/out_'+sample+'_'+ver+'.json'
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

def test_rhalphabet(tmpdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,ver,isData=True):
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
    rho_stop  = -1.6
    rhoscaled = (rhopts - rho_start) / (rho_stop - rho_start)
    validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
    rhoscaled[~validbins] = 1  # we will mask these out later

    # Build qcd MC pass+fail model and fit to polynomial
    ### This model is for the prefit but helps to compute the ratio pass/fail
    qcdmodel = rl.Model("qcdmodel")
    qcdpass, qcdfail = 0., 0.
    for ptbin in range(npt):
        failCh = rl.Channel("ptbin%d%s" % (ptbin, 'fail'))
        passCh = rl.Channel("ptbin%d%s" % (ptbin, 'pass'))
        qcdmodel.addChannel(failCh)
        qcdmodel.addChannel(passCh)
        # mock template
        ptnorm = 1
#        failTempl = expo_sample(norm=ptnorm*1e5, scale=40, obs=msd)
#        passTempl = expo_sample(norm=ptnorm*1e3, scale=40, obs=msd)
        failTempl = load_from_json(dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Fail', rebin_factor, msd, ver)
        passTempl = load_from_json(dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Pass', rebin_factor, msd, ver)
        failCh.setObservation(failTempl)
        passCh.setObservation(passTempl)
        qcdfail += failCh.getObservation().sum()
        qcdpass += passCh.getObservation().sum()

    qcdeff = qcdpass / qcdfail
    if args.runPrefit:
        tf_MCtempl = rl.BernsteinPoly("tf_MCtempl", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-50, 50))
        tf_MCtempl_params = qcdeff * tf_MCtempl(ptscaled, rhoscaled)
        for ptbin in range(npt):
            failCh = qcdmodel['ptbin%dfail' % ptbin]
            passCh = qcdmodel['ptbin%dpass' % ptbin]
            failObs = failCh.getObservation()
            qcdparams = np.array([rl.IndependentParameter('qcdparam_ptbin%d_msdbin%d' % (ptbin, i), 0) for i in range(msd.nbins)])
            sigmascale = 10.
            scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
            fail_qcd = rl.ParametericSample('ptbin%dfail_qcd' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('ptbin%dpass_qcd' % ptbin, rl.Sample.BACKGROUND, tf_MCtempl_params[ptbin, :], fail_qcd)
            passCh.addSample(pass_qcd)

            failCh.mask = validbins[ptbin]
            passCh.mask = validbins[ptbin]

        qcdfit_ws = ROOT.RooWorkspace('qcdfit_ws')
        simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
        qcdfit = simpdf.fitTo(obs,
                              ROOT.RooFit.Extended(True),
                              ROOT.RooFit.SumW2Error(True),
                              ROOT.RooFit.Strategy(2),
                              ROOT.RooFit.Save(),
                              ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                              ROOT.RooFit.PrintLevel(-1),
                              )
        qcdfit_ws.add(qcdfit)
        if "pytest" not in sys.modules:
             qcdfit_ws.writeToFile(os.path.join(str(tmpdir), 'ttHbb_qcdfit.root'))
        if qcdfit.status() != 0:
            raise RuntimeError('Could not fit qcd')

        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', qcdfit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

    #### Actual transfer function for combine
    tf_dataResidual = rl.BernsteinPoly("tf_dataResidual", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-50, 50))
    tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
    #tf_params = qcdeff * (tf_MCtempl_params_final * tf_dataResidual_params if args.runPrefit else tf_dataResidual_params)
    tf_params = qcdeff * tf_dataResidual_params
    if args.runPrefit:
        tf_params *= tf_MCtempl_params_final

    # build actual fit model now
    model = rl.Model("ttHbb")

    for ptbin in range(npt):
        for region in ['Pass', 'Fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            templates = {
                'signal'     : load_from_json('signal', ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor, msd, ver),
                'background' : load_from_json(dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor,  msd, ver),
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
#            sample.setParamEffect(lumi, 1.027)

            ch.addSample(sample)

            #import pdb
            #pdb.set_trace()
            yields = sum(tpl[0] for tpl in templates.values())
            data_obs = (yields, msd.binning, msd.name)
            #data_obs = templates[dataOrBkg]
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

        qcdparams = np.array([rl.IndependentParameter('qcdparam_ptbin%d_msdbin%d' % (ptbin, i), 0) for i in range(msd.nbins)])
        initial_qcd = failCh.getObservation().astype(float)  # was integer, and numpy complained about subtracting float from it
        for sample in failCh:
            initial_qcd -= sample.getExpectation(nominal=True)
        if np.any(initial_qcd < 0.):
            raise ValueError("initial_qcd negative for some bins..", initial_qcd)
        sigmascale = 10  # to scale the deviation from initial
        scaledparams = initial_qcd * (1 + sigmascale/np.maximum(1., np.sqrt(initial_qcd)))**qcdparams
        fail_qcd = rl.ParametericSample('ptbin%dFail_qcd' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams)
        failCh.addSample(fail_qcd)
        pass_qcd = rl.TransferFactorSample('ptbin%dPass_qcd' % ptbin, rl.Sample.BACKGROUND, tf_params[ptbin, :], fail_qcd)
        passCh.addSample(pass_qcd)

        #tqqpass = passCh['tqq']
        #tqqfail = failCh['tqq']
        #tqqPF = tqqpass.getExpectation(nominal=True).sum() / tqqfail.getExpectation(nominal=True).sum()
        #tqqpass.setParamEffect(tqqeffSF, 1*tqqeffSF)
        #tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
        #tqqpass.setParamEffect(tqqnormSF, 1*tqqnormSF)
        #tqqfail.setParamEffect(tqqnormSF, 1*tqqnormSF)

    with open(os.path.join(str(tmpdir), 'ttHbb.pkl'), "wb") as fout:
        pickle.dump(model, fout)

    pref = 'data' if isData else 'mc'
    model.renderCombine(os.path.join(str(tmpdir), pref+'_msd%dto%d_msdbin%d_pt%dbin_polyDegs%d%d'%(msd_start,msd_stop,rebin_factor,len(ptbins)-1, polyDegPt,polyDegRho)))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
  parser.add_argument('-f', '--runPrefit', action='store_true', help='Run prefit on MC.' )
  parser.add_argument('-smr', '--scanMassRange', action='store_true', default=False, help='option to scan mass range')
  parser.add_argument('-spd', '--scanPolyDeg', action='store_true', default=False, help='option to scan degree of polynomial')
  parser.add_argument('--msd_start', default=90, type=int, help='start of the mass range')
  parser.add_argument('--msd_stop', default=170, type=int, help='stop of the mass range')
  parser.add_argument('--polyDegPt', default=2, type=int, help='degree of polynomial to fit pt')
  parser.add_argument('--polyDegRho', default=2, type=int, help='degree of polynomial to fit rho')
  parser.add_argument('-r', '--rebin_factor', default=5, type=int, help='rebin factor for json histograms, default mass bin size is 1GeV')
  parser.add_argument('--nptbins', default=2, type=int, help='number of pt bins')
  parser.add_argument('-v', '--version', help='version, in file names')

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
    print('invalid number of ptbins')
    sys.exit(0)

  msd_start    = args.msd_start
  msd_stop     = args.msd_stop
  polyDegPt    = args.polyDegPt
  polyDegRho   = args.polyDegRho
  rebin_factor = args.rebin_factor
  ver          = args.version

  folder = 'ttH_'+ver
  if not os.path.exists(folder):
      os.mkdir(folder)

  if not (args.scanMassRange or args.scanPolyDeg):
    test_rhalphabet(folder,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,ver,args.isData)
  else:
    pref = 'data' if args.isData else 'mc'
    if args.scanMassRange:
      for msd_start in range(70,105,rebin_factor):
        for msd_stop in range(150,251,rebin_factor):
          with open(os.path.join(folder,'scan_'+pref+'_fitRange_%dptbin_msdbin%dgev_polyDegs%d%d.txt'%(len(ptbins)-1, rebin_factor, polyDegPt,polyDegRho)),'a') as f:
            try:
              print('trying %d to %d'%(msd_start,msd_stop))
              test_rhalphabet(folder,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,ver,args.isData)
              f.write('msd range %d to %d worked!!\n' %(msd_start, msd_stop))
            except:
              f.write('msd range %d to %d failed\n' %(msd_start, msd_stop))

    if args.scanPolyDeg:
      for polyDegPt in range(0,5):
        for polyDegRho in range(0,10):
          with open(os.path.join(folder,'scan_'+pref+'_polyDeg_%dptbin_msdbin%dgev_fitRange%dto%d.txt'%(len(ptbins)-1, rebin_factor, msd_start, msd_stop)),'a') as f:
            try:
              print('trying ptdeg %d and rhodeg %d'%(polyDegPt,polyDegRho))
              test_rhalphabet(folder,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,ver,args.isData)
              f.write('polyDegPt %d and polyDegRho %d worked!!\n'%(polyDegPt,polyDegRho))
            except:
              f.write('polyDegPt %d and polyDegRho %d failed\n'%(polyDegPt,polyDegRho))

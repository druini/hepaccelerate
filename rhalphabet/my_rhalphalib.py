from __future__ import print_function
import sys, os, argparse
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import ROOT
import json
ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()
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
  if ptStop==2000: ptStop = 5000
  with open(filepath) as json_file:
    data = json.load(json_file)
    data = data['hist_leadAK8JetMass_2J2WdeltaR_'+region+'_pt%sto%s' % (ptStart, ptStop)]
    rebin(data,rebin_factor)
  assert( np.all(np.array(obs.binning)==np.array(data['edges'])[msd_start_idx:msd_stop_idx+1]) )
  return(np.array(data['contents'])[msd_start_idx:msd_stop_idx], obs.binning, obs.name)

def loadTH1_from_json(indir, sample, ptStart, ptStop, msd_start_idx, msd_stop_idx, region, rebin_factor, obs):
  #filepath = '/afs/cern.ch/work/d/druini/public/hepaccelerate/results/2018/v05/'+ver+'/out_'+sample+'_merged.json'
  if sample.startswith('back'):
      if args.year.startswith('2016'): tmpName = '_noQCD_noDY'
      elif args.year.startswith('2017'): tmpName = ''
      elif args.year.startswith('2018'): tmpName = '_noDY'
  else: tmpName = ''
  filepath = os.path.join(indir, 'out_'+sample+tmpName+'_merged.json')
  if ptStop==2000: ptStop = 5000
  with open(filepath) as json_file:
    data = json.load(json_file)
    data = data['hist_leadAK8JetMass_2J2WdeltaR_'+region+'_pt%sto%s' % (ptStart, ptStop)]
    rebin(data,rebin_factor)
  assert( np.all(np.array(obs.binning)==np.array(data['edges'])[msd_start_idx:msd_stop_idx+1]) )
  tmpHisto = ROOT.TH1F( obs.name, 'hist_leadAK8JetMass_2J2WdeltaR_'+region+'_pt%sto%s' % (ptStart, ptStop), len(obs.binning)-1, data['edges'][msd_start_idx], data['edges'][msd_stop_idx])
  for nBin in range( len( data['contents'][msd_start_idx:msd_stop_idx] ) ):
      #print( data['contents'][msd_start_idx+nBin], ROOT.TMath.Sqrt(data['contents'][msd_start_idx+nBin]), ROOT.TMath.Sqrt(data['contents_w2'][msd_start_idx+nBin] ) )
      tmpHisto.SetBinContent( nBin+1, data['contents'][msd_start_idx+nBin] )
      tmpHisto.SetBinError( nBin+1, ROOT.TMath.Sqrt(data['contents_w2'][msd_start_idx+nBin] ) )

  return(tmpHisto)

#def rebin(bins, counts, yerr, rebin_factor):
def rebin(hist, rebin_factor):
    bins   = np.array(hist['edges'])#[:-1])
    counts = np.array(hist['contents'])
    countsErr = np.array(hist['contents_w2'])

    new_bins   = bins[::rebin_factor]
    new_counts = np.add.reduceat(counts, range(0, len(counts), rebin_factor))
    new_countsErr = np.add.reduceat(countsErr, range(0, len(countsErr), rebin_factor))
#    new_yerr   = np.add.reduceat(yerr, range(0, len(yerr), rebin_factor))
#    return new_bins, new_counts#, new_yerr
    hist['edges']    = new_bins
    hist['contents'] = new_counts
    hist['contents_w2'] = new_countsErr


def test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,isData=True,runExp=False):
    #dataOrBkg = 'data' if isData else ('all' if args.sig_and_bkg else 'background')
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
        failTempl = loadTH1_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Fail', rebin_factor, msd)
        passTempl = loadTH1_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, 'Pass', rebin_factor, msd)
        failCh.setObservation(failTempl)
        passCh.setObservation(passTempl)
        bkgfail += failCh.getObservation().sum()
        bkgpass += passCh.getObservation().sum()

    bkgeff = bkgpass / bkgfail
    if args.runPrefit:
        tf_MCtempl = rl.BernsteinPoly("tf_MCtempl", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-args.poly_limit, args.poly_limit), coefficient_transform=(np.exp if runExp else None))
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
        bkgfit.Print()

        bkgfit_ws.add(bkgfit)
        if "pytest" not in sys.modules:
             bkgfit_ws.writeToFile(os.path.join(str(outdir), 'ttHbb_bkgfit.root'))
        if bkgfit.status() != 0:
            raise RuntimeError('Could not fit bkg')

        bkgpar_names = [p for p in bkgfit.floatParsFinal().contentsString().split(',') if 'tf' in p]
        #for pn in bkgpar_names: print( pn, round(bkgfit.floatParsFinal().find(pn).getVal(), 4))
        prefit_bkgpar = np.array([ round(bkgfit.floatParsFinal().find(pn).getVal(), 4) for pn in bkgpar_names ])
        prefit_bkgparerror = [ round(bkgfit.floatParsFinal().find(pn).getError(), 4) for pn in bkgpar_names ]
        #decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', bkgfit, param_names)
        #tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        #tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

    #### Actual transfer function for combine
    polyLimit = max(prefit_bkgparerror)*5 if args.runPrefit else args.poly_limit
    #polyLimit = args.poly_limit
    tf_dataResidual = rl.BernsteinPoly("tf_dataResidual", (polyDegPt, polyDegRho), ['pt', 'rho'], limits=(-polyLimit, polyLimit), coefficient_transform=(np.exp if runExp else None), init_params=(prefit_bkgpar.reshape(polyDegPt+1, polyDegRho+1) if args.runPrefit else None))
    tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
    tf_params = bkgeff * tf_dataResidual_params
    #if args.runPrefit:
    #    tf_params *= tf_MCtempl_params_final

    # build actual fit model now
    str_polylims = ('data_' if args.isData else ('sig_' if args.sig_and_bkg else ''))+"polylims%ito%i"%(-args.poly_limit,args.poly_limit)
    model = rl.Model('ttHbb_'+str_polylims)

    for ptbin in range(npt):
        for region in ['Pass', 'Fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            templates = {
                'signal'     : loadTH1_from_json(indir, 'signal', ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor, msd),
                'background' : loadTH1_from_json(indir, dataOrBkg, ptbins[ptbin], ptbins[ptbin+1], msd_start_idx, msd_stop_idx, region, rebin_factor,  msd),
            }
            # some mock expectations
            templ = templates['signal']
            stype = rl.Sample.SIGNAL
            sample = rl.TemplateSample(ch.name + '_signal', stype, templ)

#            # mock systematics
            jecup_ratio = np.random.normal(loc=1, scale=0.05, size=msd.nbins)
            msdUp = np.linspace(0.9, 1.1, msd.nbins)
            msdDn = np.linspace(1.2, 0.8, msd.nbins)
#
#            # for jec we set lnN prior, shape will automatically be converted to norm systematic
            #sample.setParamEffect(jec, jecup_ratio)
            #sample.setParamEffect(massScale, msdUp, msdDn)
            sample.setParamEffect(lumi, 1.027)      ### for bias/ftest/GOF signal needs at least one unc.

            ch.addSample(sample)

            #import pdb
            #pdb.set_trace()
            #yields = sum(tpl[0] for tpl in templates.values())
            #data_obs = (yields, msd.binning, msd.name)
            if not isData and args.sig_and_bkg: templates['background'].Add( templates['signal'] )
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

    pref = ('data' if isData else 'mc'+( 'SB' if args.sig_and_bkg else '' ) )
    combineFolder = os.path.join(str(outdir), pref+'_msd%dto%d_msdbin%d_pt%dbin_%spolyDegs%d%d'%(msd_start,msd_stop,rebin_factor,len(ptbins)-1,('exp' if args.runExp else ''), polyDegPt,polyDegRho))
    model.renderCombine(combineFolder)
    #exec_me('bash build.sh | combine -M FitDiagnostics ttHbb_%s_combined.txt -n _r%ito%i_%si -t -1 --expectSignal 0 '%(str_polylims,args.rMin,args.rMax), folder=combineFolder)
    #exec_me('bash build.sh | combine -M FitDiagnostics ttHbb_%s_combined.txt -n _r%ito%i_%s --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace --expectSignal 0 -t -1 --toysFrequentist'%(str_polylims,args.rMin,args.rMax,str_polylims), folder=combineFolder)
    exec_me('bash build.sh | combine -M FitDiagnostics ttHbb_%s_combined.txt -n _r%ito%i_%s --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace --setParameterRanges r=%i,%i'%(str_polylims,args.rMin,args.rMax,str_polylims,args.rMin,args.rMax), folder=combineFolder)
    if args.runImpacts:
        exec_me( 'combineTool.py -M Impacts -d ttHbb_'+str_polylims+'_combined.root -m 125 --doInitialFit --robustFit 1' )
        exec_me( 'combineTool.py -M Impacts -d ttHbb_'+str_polylims+'_combined.root -m 125 --doFits --robustFit 1' )
        exec_me( 'combineTool.py -M Impacts -d ttHbb_'+str_polylims+'_combined.root -m 125 -o impacts.json' )
        exec_me( 'plotImpacts.py -i impacts.json -o impacts' )

    ##### Priting parameters
    rootFile = ROOT.TFile.Open(os.getcwd()+'/fitDiagnostics_r'+str(int(args.rMin))+'to'+str(int(args.rMax))+'_'+str_polylims+'.root')
    par_names = rootFile.Get('fit_s').floatParsFinal().contentsString().split(',')
    par_names = [p for p in par_names if 'tf' in p]
    for pn in par_names:
        print( 'fit_s', pn, round(rootFile.Get('fit_s').floatParsFinal().find(pn).getVal(), 4), '+/-', round(rootFile.Get('fit_s').floatParsFinal().find(pn).getError(), 4))
    par_names = rootFile.Get('fit_b').floatParsFinal().contentsString().split(',')
    par_names = [p for p in par_names if 'tf' in p]
    for pn in par_names:
        print( 'fit_b', pn, round(rootFile.Get('fit_b').floatParsFinal().find(pn).getVal(), 4), '+/-', round(rootFile.Get('fit_b').floatParsFinal().find(pn).getError(), 4))
    if args.runPrefit:
        print(prefit_bkgpar)
        print(prefit_bkgparerror)

def simpleFit(indir,outdir,msd_start,msd_stop,polyDegPt,rebin_factor,ptbins,isData=True,runExp=False):
    dataOrBkg = 'data' if isData else 'background'

    msdbins = np.linspace(0,300,301)[::rebin_factor]
    msd_start_idx = np.where(msdbins==msd_start)[0][0]
    msd_stop_idx  = np.where(msdbins==msd_stop)[0][0]
    msdbins = msdbins[msd_start_idx:msd_stop_idx+1]
    msd = rl.Observable('msd', msdbins)

    templates = {
        'signal'     : loadTH1_from_json(indir, 'signal', ptbins[0], ptbins[-1], msd_start_idx, msd_stop_idx, 'Pass', rebin_factor, msd),
        'background' : loadTH1_from_json(indir, dataOrBkg, ptbins[0], ptbins[-1], msd_start_idx, msd_stop_idx, 'Pass', rebin_factor,  msd),
    }
    if not isData and args.sig_and_bkg: templates['background'].Add( templates['signal'] )

    msd = ROOT.RooRealVar('msd', 'msd', msd_start, msd_stop)
    data_obs = ROOT.RooDataHist('data_obs', 'data_obs', ROOT.RooArgList(msd), templates['background'] )
    signal = ROOT.RooDataHist('signal', 'signal', ROOT.RooArgList(msd), templates['signal'] )

    polyArgList = ROOT.RooArgList( )
    rooDict = {}
    for i in range( int(polyDegPt) ):
        rooDict[ 'a'+str(i) ] = ROOT.RooRealVar('a'+str(i), 'a'+str(i), 1, -1000, 1000 )
        polyArgList.add( rooDict[ 'a'+str(i) ] )
    rooDict['bkgFunc'] = ROOT.RooBernstein("bkg", "bkg", msd, polyArgList ) if not args.runExp else ROOT.RooChebychev('bkg', 'bkg', msd, polyArgList)
    rooDict['bkg_norm'] = ROOT.RooRealVar( 'bkg_norm', 'bkg_norm', int(templates['background'].Integral()), 0, 100000 )

    #simpdf, obs = bkgmodel.renderRoofit(bkgfit_ws)
    bkgfit = rooDict['bkgFunc'].fitTo(data_obs,
                          #ROOT.RooFit.Extended(True),
                          ROOT.RooFit.SumW2Error(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Save(),
                          #ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                          ROOT.RooFit.PrintLevel(-1),
                          )
    bkgfit.Print()
    prefit_bkgpar = [ bkgfit.floatParsFinal().find('a'+str(i)).getVal() for i in range(polyDegPt) ]
    prefit_bkgparerror = [ bkgfit.floatParsFinal().find('a'+str(i)).getError() for i in range(polyDegPt) ]


    pref = ('data' if isData else 'mc'+( 'SB' if args.sig_and_bkg else '' ) )
    combineFolder = os.path.join(str(outdir), pref+'_msd%dto%d_msdbin%d_%spolyDegs%d'%(msd_start,msd_stop,rebin_factor,('exp' if args.runExp else ''), polyDegPt))
    try: os.makedirs(combineFolder)
    except OSError: print('|===>'+combineFolder+' Folder already exist')

    ##### simple plot of fit
    c1 = ROOT.TCanvas('c1', 'c1',  10, 10, 750, 500 )
    #c1.SetLogy()
    xframe = msd.frame()
    data_obs.plotOn( xframe )
    xframe.Draw()
    xframe.GetXaxis().SetTitle("Leading Jet Mass [GeV]")
    rooDict['bkgFunc'].plotOn( xframe )
    rooDict['bkgFunc'].paramOn( xframe, ROOT.RooFit.Layout(0.6,0.9,0.94))
    xframe.Draw()
    #xframe.SetMaximum(100000)
    #xframe.SetMinimum(0.00001)
    c1.SaveAs( combineFolder+'/test_simpleFit_rootfit.png')
    del c1

    ws = ROOT.RooWorkspace('ws')
    getattr(ws, 'import')( data_obs )
    getattr(ws, 'import')( signal )
    getattr(ws, 'import')( rooDict['bkgFunc'] )
    getattr(ws, 'import')( rooDict['bkg_norm'] )
    ws.writeToFile( combineFolder+'/ws_ttHbb.root' )
    ws.Print()

    combineLabel='_r'+str(args.rMin)+'to'+str(args.rMax)
    datacardLabel='ttHbb'+combineLabel+'.txt'
    datacard = open( combineFolder+'/'+datacardLabel, 'w')
    datacard.write("imax 1 number of bins \n")
    datacard.write("jmax * number of processes minus 1 \n")
    datacard.write("kmax * number of nuisance parameters \n")
    datacard.write("-------------------------------\n")
    datacard.write("shapes * * ws_ttHbb.root ws:$PROCESS\n")
    datacard.write("-------------------------------\n")
    datacard.write("bin           pass\n")
    datacard.write("observation   -1\n")
    datacard.write("-------------------------------\n")
    datacard.write("bin           pass      pass\n")
    datacard.write("process       bkg       signal\n")
    datacard.write("process       1         0\n")
    datacard.write('rate          1         -1\n')
    datacard.write("-------------------------------\n")
    datacard.write("lumi    lnN     1.025         -     \n")
    for q in range( int(polyDegPt) ):
        #datacard.write('a'+str(q)+"    flatParam\n")
        datacard.write('a'+str(q)+"    param    "+str(prefit_bkgpar[q])+"     "+str(prefit_bkgparerror[q])+"\n")
    datacard.close()

    exec_me('combine -M FitDiagnostics %s -n %s --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace --setParameterRanges r=%i,%i'%(datacardLabel,combineLabel,args.rMin,args.rMax), folder=combineFolder)

    if args.runImpacts:
        exec_me( 'combineTool.py -M Impacts -d '+datacardLabel+' -m 125 --doInitialFit --robustFit 1' )
        exec_me( 'combineTool.py -M Impacts -d '+datacardLabel+' -m 125 --doFits --robustFit 1' )
        exec_me( 'combineTool.py -M Impacts -d '+datacardLabel+' -m 125 -o impacts.json' )
        exec_me( 'plotImpacts.py -i impacts.json -o impacts' )

    ##### Priting parameters
    rootFile = ROOT.TFile.Open(os.getcwd()+'/fitDiagnostics_r'+str(int(args.rMin))+'to'+str(int(args.rMax))+'.root')
    par_names = rootFile.Get('fit_s').floatParsFinal().contentsString().split(',')
    par_names = [p for p in par_names if 'a' in p]
    for pn in par_names:
        print( 'fit_s', pn, round(rootFile.Get('fit_s').floatParsFinal().find(pn).getVal(), 4), '+/-', round(rootFile.Get('fit_s').floatParsFinal().find(pn).getError(), 4))
    par_names = rootFile.Get('fit_b').floatParsFinal().contentsString().split(',')
    par_names = [p for p in par_names if 'a' in p]
    for pn in par_names:
        print( 'fit_b', pn, round(rootFile.Get('fit_b').floatParsFinal().find(pn).getVal(), 4), '+/-', round(rootFile.Get('fit_b').floatParsFinal().find(pn).getError(), 4))
    print(prefit_bkgpar)
    print(prefit_bkgparerror)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
  parser.add_argument('--sig-and-bkg', action='store_true', default=False, help='sum signal and background samples when running on MC')
  parser.add_argument('--simpleFit', action='store_true', default=False, help='sum signal and background samples when running on MC')
  parser.add_argument('-f', '--runPrefit', action='store_true', help='Run prefit on MC.' )
  parser.add_argument('-i', '--runImpacts', action='store_true', help='Run impacts.' )
  #parser.add_argument('-e', '--runExp', action='store_true', help='Run with exponential Bernstein.' )
  parser.add_argument('--pdf', default='poly', choices=['poly','exp'])
  parser.add_argument('-smr', '--scanMassRange', action='store_true', default=False, help='option to scan mass range')
  parser.add_argument('-spd', '--scanPolyDeg', action='store_true', default=False, help='option to scan degree of polynomial')
  parser.add_argument('--msd_start', default=90, type=int, help='start of the mass range')
  parser.add_argument('--msd_stop', default=170, type=int, help='stop of the mass range')
  parser.add_argument('--polyDegPt', default=2, type=int, help='degree of polynomial to fit pt')
  parser.add_argument('--polyDegRho', default=2, type=int, help='degree of polynomial to fit rho')
  parser.add_argument('-r', '--rebin_factor', default=5, type=int, help='rebin factor for json histograms, default mass bin size is 1GeV')
  parser.add_argument('--nptbins', default=2, type=int, help='number of pt bins')
  parser.add_argument('--rMin', default=-20, type=float, help='minimum of r (signal strength) in combine fit')
  parser.add_argument('--rMax', default=20, type=float, help='maximum of r (signal strength) in combine fit')
  parser.add_argument('-l','--poly-limit', default=2, type=int, help='sets the limit for the parameters of the Bernsetin polynomial')
  parser.add_argument('-y', '--year', default='2017', type=str, help='year to process, in file paths')
  parser.add_argument('-v', '--version', default='v05', help='version, in file paths')
  parser.add_argument('-s', '--selection', nargs='+', default=['met20_btagDDBvL_noMD07','met20_deepTagMD_bbvsLight05845','met20_deepTagMD_bbvsLight08695'], help='event selection, in file paths')
  parser.add_argument('-j', '--jsonpath', default='/afs/cern.ch/work/d/druini/public/hepaccelerate/results', help='path to json files')
  parser.add_argument('-o','--outdir', default=None, help='specifiy a custom output directory')

  try: args = parser.parse_args()
  except:
    parser.print_help()
    sys.exit(0)

  if args.pdf=='poly':
    args.runExp = False
  else:
    args.runExp = True

  pt_start = 200 if args.version=='v06' else 250
  if args.nptbins==1:
    ptbins = np.array([pt_start,2000])
  elif args.nptbins==2:
    ptbins = np.array([pt_start,300,2000])
  elif args.nptbins==3:
    ptbins = np.array([pt_start,300,450,2000])
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

    if args.outdir is None: outdir = os.path.join('output', args.year, args.version, sel)
    else: outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.simpleFit:
        ptbins = np.array([pt_start,2000])
        simpleFit(indir,outdir,msd_start,msd_stop,polyDegPt,rebin_factor,ptbins,args.isData,runExp=args.runExp)
    else:
        if not (args.scanMassRange or args.scanPolyDeg):
          test_rhalphabet(indir,outdir,msd_start,msd_stop,polyDegPt,polyDegRho,rebin_factor,ptbins,args.isData,runExp=args.runExp)
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

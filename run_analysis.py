import os, glob
#os.environ["NUMBAPRO_NVVM"] = "/usr/local/cuda/nvvm/lib64/libnvvm.so"
#os.environ["NUMBAPRO_LIBDEVICE"] = "/usr/local/cuda/nvvm/libdevice/"
os.environ['KERAS_BACKEND'] = "tensorflow"
import argparse
import json
import numpy as np

import uproot
from uproot_methods import TLorentzVectorArray
#import hepaccelerate
from hepaccelerate.utils import Results, NanoAODDataset, Histogram, choose_backend

#import tensorflow as tf
#from keras.models import load_model
import itertools
#from lib_analysis import mse0,mae0,r2_score0

from definitions_analysis import histogram_settings

import lib_analysis
from lib_analysis import vertex_selection, lepton_selection, jet_selection, load_puhist_target, compute_pu_weights, compute_lepton_weights, compute_btag_weights, chunks, evaluate_DNN, calculate_variable_features, select_lepton_p4, hadronic_W, get_histogram

#config = tf.ConfigProto()
#config.gpu_options.allow_growth=True
#sess = tf.Session(config=config)

#This function will be called for every file in the dataset
def analyze_data(data, sample, NUMPY_LIB=None, parameters={}, samples_info={}, is_mc=True, lumimask=None, cat=False, boosted=False, DNN=False, DNN_model=None):
    #Output structure that will be returned and added up among the files.
    #Should be relatively small.
    ret = Results()

    muons = data["Muon"]
    electrons = data["Electron"]
    scalars = data["eventvars"]
    jets = data["Jet"]
    jets.p4 = TLorentzVectorArray.from_ptetaphim(jets.pt, jets.eta, jets.phi, jets.mass)
    fatjets = data["FatJet"]

    METp4 = TLorentzVectorArray.from_ptetaphim(scalars["MET_pt"], 0, scalars["MET_phi"], 0)
    nEvents = muons.numevents()

    indices = {
        "leading"    : NUMPY_LIB.zeros(nEvents, dtype=NUMPY_LIB.int32),
        "subleading" : NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.int32)
        }

    mask_events = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.bool)

    # apply event cleaning, PV selection and trigger selection
    flags = [
        "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"]#, "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter"]
    if not is_mc:
        flags.append("Flag_eeBadScFilter")
    for flag in flags:
        mask_events = mask_events & scalars[flag]
    if args.year.startswith('2016'):
        trigger = (scalars["HLT_Ele27_WPTight_Gsf"] | scalars["HLT_IsoMu24"]  | scalars["HLT_IsoTkMu24"])
    elif args.year.startswith('2017'):
        trigger = (scalars["HLT_Ele35_WPTight_Gsf"] | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"] | scalars["HLT_IsoMu27"] | scalars["HLT_IsoMu24_eta2p1"]) #FIXME for different runs
    elif args.year.startswith('2018'):
        trigger = (scalars["HLT_Ele32_WPTight_Gsf"] | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"] | scalars["HLT_IsoMu24"] )
    mask_events = mask_events & trigger
    mask_events = mask_events & (scalars["PV_npvsGood"]>0)

    # apply object selection for muons, electrons, jets
    good_muons, veto_muons = lepton_selection(muons, parameters["muons"])
    good_electrons, veto_electrons = lepton_selection(electrons, parameters["electrons"])
    good_jets = jet_selection(jets, muons, good_muons, parameters["jets"]) & jet_selection(jets, electrons, good_electrons, parameters["jets"])
#    good_jets = jet_selection(jets, muons, (veto_muons | good_muons), parameters["jets"]) & jet_selection(jets, electrons, (veto_electrons | good_electrons) , parameters["jets"])
#    bjets = good_jets & (getattr(jets, parameters["btagging algorithm"]) > parameters["btagging WP"])
    good_fatjets = jet_selection(fatjets, muons, good_muons, parameters["fatjets"]) & jet_selection(fatjets, electrons, good_electrons, parameters["fatjets"])
#    good_fatjets = jet_selection(fatjets, muons, (veto_muons | good_muons), parameters["fatjets"]) & jet_selection(fatjets, electrons, (veto_electrons | good_electrons), parameters["fatjets"]) #FIXME remove vet_leptons

#    higgs_candidates = good_fatjets & (fatjets.pt > 250)
#    nhiggs = ha.sum_in_offsets(fatjets, higgs_candidates, mask_events, fatjets.masks["all"], NUMPY_LIB.int8)
#    indices["best_higgs_candidate"] = ha.index_in_offsets(fatjets.pt, fatjets.offsets, 1, mask_events, higgs_candidates)
#    best_higgs_candidate = NUMPY_LIB.zeros_like(higgs_candidates)
#    best_higgs_candidate[ (fatjets.offsets[:-1] + indices["best_higgs_candidate"])[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] ] = True
#    best_higgs_candidate[ (fatjets.offsets[:-1] + indices["best_higgs_candidate"])[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] ] &= nhiggs.astype(NUMPY_LIB.bool)[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] # to avoid removing the leading fatjet in events with no higgs candidate

    good_jets_nohiggs = good_jets & ha.mask_deltar_first(jets, good_jets, fatjets, good_fatjets, 1.2, indices['leading'])
    bjets = good_jets_nohiggs & (getattr(jets, parameters["btagging algorithm"]) > parameters["btagging WP"])
    nonbjets = good_jets_nohiggs & (getattr(jets, parameters["btagging algorithm"]) < parameters["btagging WP"])

    # apply basic event selection -> individual categories cut later
    nmuons     = ha.sum_in_offsets(muons, good_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8)
    nelectrons = ha.sum_in_offsets(electrons, good_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8)
    nleps      = NUMPY_LIB.add(nmuons, nelectrons)
#    lepton_veto = NUMPY_LIB.add(ha.sum_in_offsets(muons, veto_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, veto_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
    njets = ha.sum_in_offsets(jets, nonbjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags = ha.sum_in_offsets(jets, bjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    nfatjets = ha.sum_in_offsets(fatjets, good_fatjets, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)
    #nhiggs = ha.sum_in_offsets(fatjets, higgs_candidates, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)

    # apply basic event selection
    #mask_events_higgs = mask_events & (nleps == 1) & (scalars["MET_pt"] > 20) & (nhiggs > 0) & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4) & (btags >=2) ) & (lepton_veto == 0)
    mask_events = mask_events & (scalars["MET_pt"] > 20) & (nfatjets > 0) #& (btags >=1)# & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4)  ) & (lepton_veto == 0)
    # for reference, this is the selection for the resolved analysis
    # mask_events = mask_events & (nleps == 1) & (lepton_veto == 0) & (njets >= 4) & (btags >=2) & met

############# calculate weights for MC samples
    weights = {}
    weights['ones'] = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.float32)
    weights["nominal"] = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.float32)

    if is_mc:
        weights["nominal"] = weights["nominal"] * scalars["genWeight"] * parameters["lumi"] * samples_info[sample]["XS"] / samples_info[sample]["ngen_weight"][args.year]

        # pu corrections
#        pu_weights = compute_pu_weights(parameters["pu_corrections_target"], weights["nominal"], scalars["Pileup_nTrueInt"], scalars["PV_npvsGood"])
        pu_weights = compute_pu_weights(parameters["pu_corrections_target"], weights["nominal"], scalars["Pileup_nTrueInt"], scalars["Pileup_nTrueInt"])
        weights["nominal"] = weights["nominal"] * pu_weights
        weights['no_lep']  = weights['nominal']

        # lepton SF corrections
        electron_weights = compute_lepton_weights(electrons, electrons.pt, (electrons.deltaEtaSC + electrons.eta), mask_events, good_electrons, evaluator, ["el_triggerSF", "el_recoSF", "el_idSF"])
        muon_weights = compute_lepton_weights(muons, muons.pt, NUMPY_LIB.abs(muons.eta), mask_events, good_muons, evaluator, ["mu_triggerSF", "mu_isoSF", "mu_idSF"])
        weights["nominal"] = weights["nominal"] * muon_weights * electron_weights

        # btag SF corrections
#        btag_weights = compute_btag_weights(jets, mask_events, good_jets, evaluator)
#        weights["nominal"] = weights["nominal"] * btag_weights

############# calculate basic variables
    leading_jet_pt        = ha.get_in_offsets(jets.pt, jets.offsets, indices['leading'], mask_events, nonbjets)
    leading_jet_eta       = ha.get_in_offsets(jets.eta, jets.offsets, indices['leading'], mask_events, nonbjets)
    leading_fatjet_SDmass = ha.get_in_offsets(fatjets.msoftdrop, fatjets.offsets, indices['leading'], mask_events, good_fatjets)
    leading_fatjet_pt     = ha.get_in_offsets(fatjets.pt, fatjets.offsets, indices['leading'], mask_events, good_fatjets)
    leading_fatjet_eta    = ha.get_in_offsets(fatjets.eta, fatjets.offsets, indices['leading'], mask_events, good_fatjets)
    leading_lepton_pt     = NUMPY_LIB.maximum(ha.get_in_offsets(muons.pt, muons.offsets, indices["leading"], mask_events, good_muons), ha.get_in_offsets(electrons.pt, electrons.offsets, indices["leading"], mask_events, good_electrons))
    leading_lepton_eta    = NUMPY_LIB.maximum(ha.get_in_offsets(muons.eta, muons.offsets, indices["leading"], mask_events, good_muons), ha.get_in_offsets(electrons.eta, electrons.offsets, indices["leading"], mask_events, good_electrons))

    lead_lep_p4        = select_lepton_p4(muons, good_muons, electrons, good_electrons, indices["leading"], mask_events)
    leading_fatjet_phi = ha.get_in_offsets(fatjets.phi, fatjets.offsets, indices['leading'], mask_events, good_fatjets)
    deltaRHiggsLepton  = ha.calc_dr(lead_lep_p4.phi, lead_lep_p4.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events)

############# masks for different selections
    mask_events_initial = mask_events
    mask_events = {}
    for lep, nlep in zip(['elec','muon'],[nelectrons,nmuons]):
      mask_events[f'{lep}_basic'] = mask_events_initial & (nleps == 1)

      mask_events[f'{lep}_2J']   = mask_events[f'{lep}_basic'] & (njets>1)

      #Ws reconstruction
      pznu = ha.METzCalculator(lead_lep_p4, METp4, mask_events[f'{lep}_2J'])
      neutrinop4 = TLorentzVectorArray.from_cartesian(METp4.x, METp4.y, pznu, NUMPY_LIB.sqrt( METp4.x**2 + METp4.y**2 + pznu**2 ))
      lepW = lead_lep_p4 + neutrinop4

      hadW = hadronic_W(jets, nonbjets, lepW, mask_events[f'{lep}_2J'])

      mask_events[f'{lep}_2J2W'] = mask_events[f'{lep}_2J'] & (hadW.mass>parameters['W']['min_mass']) & (hadW.mass<parameters['W']['max_mass']) & (lepW.mass>parameters['W']['min_mass']) & (lepW.mass<parameters['W']['max_mass'])

      #deltaR between objects
      deltaRlepWHiggs = ha.calc_dr(lepW.phi, lepW.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events[f'{lep}_2J2W'])
      deltaRhadWHiggs = ha.calc_dr(hadW.phi, hadW.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events[f'{lep}_2J2W'])

  #    mask_events['2J2WdeltaR'] = mask_events['2J2W'] & (deltaRlepWHiggs>1.5) & (deltaRhadWHiggs>1.5) & (deltaRlepWHiggs<4) & (deltaRhadWHiggs<4)
      mask_events[f'{lep}_2J2WdeltaR'] = mask_events[f'{lep}_2J2W'] & (deltaRlepWHiggs>1) & (deltaRhadWHiggs>1)# & (deltaRlepWHiggs<4) & (deltaRhadWHiggs<4)

      #boosted Higgs
      leading_fatjet_tau1 = ha.get_in_offsets(fatjets.tau1, fatjets.offsets, indices['leading'], mask_events[f'{lep}_2J2WdeltaR'], good_fatjets)
      leading_fatjet_tau2 = ha.get_in_offsets(fatjets.tau2, fatjets.offsets, indices['leading'], mask_events[f'{lep}_2J2WdeltaR'], good_fatjets)
      leading_fatjet_tau21 = NUMPY_LIB.divide(leading_fatjet_tau2, leading_fatjet_tau1)
      mask_events[f'{lep}_2J2WdeltaRTau21'] = mask_events[f'{lep}_2J2WdeltaR'] & (leading_fatjet_tau21<parameters["fatjets"]["tau21cut"][args.year])

      leading_fatjet_Hbb = ha.get_in_offsets(getattr(fatjets, parameters["bbtagging algorithm"]), fatjets.offsets, indices['leading'], mask_events[f'{lep}_2J2WdeltaRTau21'], good_fatjets)
      mask_events[f'{lep}_2J2WdeltaRTau21_Pass'] = mask_events[f'{lep}_2J2WdeltaRTau21'] & (leading_fatjet_Hbb>parameters['bbtagging WP'])
      mask_events[f'{lep}_2J2WdeltaRTau21_Fail'] = mask_events[f'{lep}_2J2WdeltaRTau21'] & (leading_fatjet_Hbb<=parameters['bbtagging WP'])

############# histograms
    vars_to_plot = {
      'nleps'             : nleps,
      'nelectrons'        : nelectrons,
      'nmuons'            : nmuons,
      'njets'             : njets,
      'btags'             : btags,
      'nfatjets'          : nfatjets,
      'met'               : scalars['MET_pt'],
      'leading_jet_pt'    : leading_jet_pt,
      'leading_jet_eta'   : leading_jet_eta,
      'leadAK8JetMass'    : leading_fatjet_SDmass,
      'leadAK8JetPt'      : leading_fatjet_pt,
      'leadAK8JetEta'     : leading_fatjet_eta,
      'leadAK8JetHbb'     : leading_fatjet_Hbb,
      'leadAK8JetTau21'   : leading_fatjet_tau21,
      'lepton_pt'         : leading_lepton_pt,
      'lepton_eta'        : leading_lepton_eta,
      'hadWPt'            : hadW.pt,
      'hadWEta'           : hadW.eta,
      'hadWMass'          : hadW.mass,
      'lepWPt'            : lepW.pt,
      'lepWEta'           : lepW.eta,
      'lepWMass'          : lepW.mass,
      'deltaRlepWHiggs'   : deltaRlepWHiggs,
      'deltaRhadWHiggs'   : deltaRhadWHiggs,
      'deltaRHiggsLepton' : deltaRHiggsLepton,
      'PV_npvsGood'       : scalars['PV_npvsGood'],
    }
    if is_mc:
      vars_to_plot['pu_weights'] = pu_weights

    var_name, var = 'leadAK8JetMass', leading_fatjet_SDmass
    ptbins = NUMPY_LIB.append( NUMPY_LIB.arange(250,600,50), [600, 675, 800, 1200] )
    for ipt in range( len(ptbins)-1 ):
      for lep in ['elec','muon']:
        for mask_name in [f'{lep}_2J2WdeltaRTau21_Pass', f'{lep}_2J2WdeltaRTau21_Fail']:
          mask = mask_events[mask_name] & (leading_fatjet_pt>ptbins[ipt]) & (leading_fatjet_pt<ptbins[ipt+1])
          ret[f'hist_leadAK8JetMass_{mask_name}_pt{ptbins[ipt]}to{ptbins[ipt+1]}'] = get_histogram( leading_fatjet_SDmass[mask], weights['nominal'][mask], NUMPY_LIB.linspace( *histogram_settings[var_name] ) )
          ret[f'hist_{var_name}_{mask_name}_pt{ptbins[ipt]}to{ptbins[ipt+1]}'] = get_histogram( var[mask], weights['nominal'][mask], NUMPY_LIB.linspace( *histogram_settings[var_name] ) )

    weight_names = {'' : 'nominal', '_NoWeights' : 'ones', '_NoLepWeights' : 'no_lep'}
    for weight_name, w in weight_names.items():
      for mask_name, mask in mask_events.items():
#        with open(f'/afs/cern.ch/work/d/druini/public/hepaccelerate/tests/events_pass_selection_{sample}_{mask_name}.txt','a+') as f:
#          for nevt, run, lumiBlock in zip(scalars['event'][mask], scalars['run'][mask], scalars['luminosityBlock']):
#            f.write(f'{nevt}, {run}, {lumiBlock}\n')
        ### tentative histogram of all jets in an event #FIXME
#        jet_feats = {'pt'  : NUMPY_LIB.array([]), 'eta' : NUMPY_LIB.array([])}
#        for evt_idx in NUMPY_LIB.where(mask)[0]:
#          start = jets.offsets[evt_idx]
#          stop  = jets.offsets[evt_idx+1]
#          for feat in ['pt','eta']:
#            jet_feats[feat] = NUMPY_LIB.append(jet_feats[feats], getattr(jets, feat)[start:stop][nonbjets[start:stop]])
#        for feat in ['pt','eta']:
#          ret[f'hist_jets_{feat}_{mask_name+weight_name}'] = get_histogram( jet_feats[feat],
        for var_name, var in vars_to_plot.items():
          if (not is_mc) and ('2J2WdeltaRTau21_Pass' in mask_name) and ('leadAK8JetMass' in var_name) : continue
          try:
            ret[f'hist_{var_name}_{mask_name+weight_name}'] = get_histogram( var[mask], weights[w][mask], NUMPY_LIB.linspace( *histogram_settings[var_name] ) )
          except KeyError:
            print(f'!!!!!!!!!!!!!!!!!!!!!!!! Please add variable {var_name} to the histogram settings')

    #synch
#    evts = [1550213, 1550290, 1550342, 1550361, 1550369, 1550387, 1550396, 1550467, 1550502, 1550566]
##    evts = [1550251, 1556872, 1557197, 1558222, 1558568, 1600001, 1602391, 3928629, 3930963, 3931311, 4086276]
#    mask = NUMPY_LIB.zeros_like(mask_events)
#    for iev in evts:
#      mask |= (scalars["event"] == iev)
#    print('nevt', scalars["event"][mask])
#    print('pass sel', mask_events[mask])
#    print('nleps', nleps[mask])
#    print('njets', njets[mask])
#    print('nfatjets', nfatjets[mask])
#    print('fatjet mass', leading_fatjet_SDmass[mask])
#    print('fatjet pt', leading_fatjet_pt[mask])
#    print('fatjet eta', leading_fatjet_eta[mask])
#    print('met', scalars['MET_pt'][mask])
#    print('lep_pt', leading_lepton_pt[mask])
#    print('lep_eta', leading_lepton_eta[mask])
#    print('pu_weight', pu_weights[mask])
#    print('lep_weight', muon_weights[mask] * electron_weights[mask])
#    print('nevents', np.count_nonzero(mask_events))

#    np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
#    evts = [1550290, 1550342, 1550361, 1550387, 1550467, 1550502, 1550607, 1550660]
#    for evt in evts:
#      evt_idx = NUMPY_LIB.where( scalars["event"] == evt )[0][0]
#      start = jets.offsets[evt_idx]
#      stop  = jets.offsets[evt_idx+1]
#      print(f'!!! EVENT {evt} !!!')
#      print(f'njets good {njets[evt_idx]}, total {stop-start}')
#      print('jets mask', nonbjets[start:stop])
#      print('jets pt', jets.pt[start:stop])
#      print('jets eta', jets.eta[start:stop])
#      print('jets btag', getattr(jets, parameters["btagging algorithm"])[start:stop])
#      print('jet Id', jets.jetId[start:stop]),
#      print('jet puId', jets.puId[start:stop])
#    with open('events_pass_selection.txt','w+') as f:
#      for nevt in scalars['event'][mask_events['2J']]:
#        f.write(str(nevt)+'\n')
#    import pdb
#    pdb.set_trace()

##
#    variables = [
##        ("jet", jets, good_jets, "leading", ["pt", "eta"]),
##        ("bjet", jets, bjets, "leading", ["pt", "eta"]),
#    ]
##
#    if boosted:
#        variables += [
##            ("fatjet", fatjets, good_fatjets, "leading",["pt", "eta", "mass", "msoftdrop", "tau32", "tau21"]),
##            ("fatjet", fatjets, good_fatjets, "subleading",["pt", "eta", "mass", "msoftdrop", "tau32", "tau21"]),
##            ("top_candidate", fatjets, top_candidates, "leading", ["pt", "eta", "mass", "msoftdrop", "tau32", "tau21"]),
##            ("WH_candidate", fatjets, WH_candidates, "inds_WHcandidates", ["pt", "eta", "mass", "msoftdrop", "tau32", "tau21"]),
##            ("higgs", genparts, higgs, "leading", ["pt", "eta"]),
##            ("tops", genparts, tops, "leading", ["pt", "eta"])
#            ("", fatjets, higgs_candidates, "best_higgs_candidate", ["pt", "msoftdrop", "tau21", parameters["bbtagging algorithm"]]),
##            ("boosted", fatjets, W_candidates, "best_W_candidate", ["pt", "msoftdrop", "tau21"]),
##            ("boosted", fatjets, top_candidates, "best_top_candidate", ["pt", "msoftdrop", "tau21"]),
#
#    ]
##
##    if boosted:
##      higgs = (genparts.pdgId == 25) & (genparts.status==62)
##      tops  = ( (genparts.pdgId == 6) | (genparts.pdgId == -6) ) & (genparts.status==62)
##      var["nfatjets"] = ha.sum_in_offsets(fatjets, good_fatjets, mask_events, fatjets.masks["all"], NUMPY_LIB.int8)
##      var["ntop_candidates"] = ha.sum_in_offsets(fatjets, tops, mask_events, fatjets.masks["all"], NUMPY_LIB.int8)
##
##    # special role of lepton
##    var["leading_lepton_pt"] = NUMPY_LIB.maximum(ha.get_in_offsets(muons.pt, muons.offsets, indices["leading"], mask_events, good_muons), ha.get_in_offsets(electrons.pt, electrons.offsets, indices["leading"], mask_events, good_electrons))
##    var["leading_lepton_eta"] = NUMPY_LIB.maximum(ha.get_in_offsets(muons.eta, muons.offsets, indices["leading"], mask_events, good_muons), ha.get_in_offsets(electrons.eta, electrons.offsets, indices["leading"], mask_events, good_electrons))
##
#    # all other variables
#    for v in variables:
#      calculate_variable_features(v, mask_events, indices, var)
#    for f in ['pt', 'mass']:
#      var['best_W_candidate_resolved_'+f] = getattr(hadW, f)
#      var['leptonic_W_candidate_'+f] = getattr(lepW, f)
#    var['deltaR_Wlep_WhadResolved'] = ha.calc_dr(lepW.phi, lepW.eta, hadW.phi, hadW.eta)
##    Whad_phi = ha.get_in_offsets(fatjets.phi, fatjets.offsets, indices["best_W_candidate"], mask_events, W_candidates)
##    Whad_eta = ha.get_in_offsets(fatjets.eta, fatjets.offsets, indices["best_W_candidate"], mask_events, W_candidates)
#    H_phi = ha.get_in_offsets(fatjets.phi, fatjets.offsets, indices["best_higgs_candidate"], mask_events, higgs_candidates)
#    H_eta = ha.get_in_offsets(fatjets.eta, fatjets.offsets, indices["best_higgs_candidate"], mask_events, higgs_candidates)
##    var['deltaR_Wlep_WhadBoosted'] = ha.calc_dr(lepW.phi, lepW.eta, Whad_phi, Whad_eta)
#    var['deltaR_Wlep_H'] = ha.calc_dr(lepW.phi, lepW.eta, H_phi, H_eta)
##    var['deltaR_H_WhadBoosted'] = ha.calc_dr(Whad_phi, Whad_eta, H_phi, H_eta)
#    var['deltaR_H_WhadResolved'] = ha.calc_dr(H_phi, H_eta, hadW.phi, hadW.eta)
#
#    #in case of data: check if event is in golden lumi file
#    if not is_mc and not (lumimask is None):
#        mask_lumi = lumimask(scalars["run"], scalars["luminosityBlock"])
#        mask_events = mask_events & mask_lumi
##
##    #evaluate DNN
##    if DNN:
##        DNN_pred = evaluate_DNN(jets, good_jets, electrons, good_electrons, muons, good_muons, scalars, mask_events, DNN, DNN_model)
##
#    # in case of tt+jets -> split in ttbb, tt2b, ttb, ttcc, ttlf
#    processes = {}
#    if sample.startswith("TT"):
#        ttCls = scalars["genTtbarId"]%100
#        processes["ttbb"] = mask_events & (ttCls >=53) & (ttCls <=56)
#        processes["tt2b"] = mask_events & (ttCls ==52)
#        processes["ttb"] = mask_events & (ttCls ==51)
#        processes["ttcc"] = mask_events & (ttCls >=41) & (ttCls <=45)
#        ttHF =  ((ttCls >=53) & (ttCls <=56)) | (ttCls ==52) | (ttCls ==51) | ((ttCls >=41) & (ttCls <=45))
#        processes["ttlf"] = mask_events & NUMPY_LIB.invert(ttHF)
#    else:
#        processes["unsplit"] = mask_events
#
#    for p in processes.keys():
#
#        mask_events_split = processes[p]
#
#        # Categories
#        categories = {}
#        if not boosted:
#          categories["sl_jge4_tge2"] = mask_events_split
#          categories["sl_jge4_tge3"] = mask_events_split & (btags >=3)
#
#          categories["sl_j4_tge3"] = mask_events_split & (njets ==4) & (btags >=3)
#          categories["sl_j5_tge3"] = mask_events_split & (njets ==5) & (btags >=3)
#          categories["sl_jge6_tge3"] = mask_events_split & (njets >=6) & (btags >=3)
#
#          categories["sl_j4_t3"] = mask_events_split & (njets ==4) & (btags ==3)
#          categories["sl_j4_tge4"] = mask_events_split & (njets ==4) & (btags >=4)
#          categories["sl_j5_t3"] = mask_events_split & (njets ==5) & (btags ==3)
#          categories["sl_j5_tge4"] = mask_events_split & (njets ==5) & (btags >=4)
#          categories["sl_jge6_t3"] = mask_events_split & (njets >=6) & (btags ==3)
#          categories["sl_jge6_tge4"] = mask_events_split & (njets >=6) & (btags >=4)
#        else:
##          categories['boosted_higgs_only'] = mask_events_split & (nhiggs>0)
#          categories['boosted_HandW'] = mask_events_split & (nhiggs>0) & ( ( (hadW.mass>65) & (hadW.mass<105) ) ) # (nW>0) |
#          for objlist in [['Wlep','H','WhadResolved']]:#, ['Wlep','H','WhadBoosted']]:
#            for obj in itertools.combinations(objlist,2):
#              categories['boosted_HandW_deltaRcut_{}_{}'.format(obj[0],obj[1])] = categories['boosted_HandW'] & (var['deltaR_{}_{}'.format(obj[0],obj[1])] > 1.5)
##          categories['boosted_higgs_and_W_dRcut'] = categories['boosted_higgs_and_W'] \
##              & (var['deltaR_Wlep_Whad_resolved'] < 3.5) \
##              & (var['deltaR_Wlep_Whad_boosted'] < 3.5) \
##              & (var['deltaR_Wlep_H'] < 3.5) \
##              & ( (var['deltaR_H_Whad_boosted'] > 2.) & (var['deltaR_H_Whad_boosted'] < 3.5) ) \
##              & ( (var['deltaR_H_Whad_resolved'] > 2.) & (var['deltaR_H_Whad_resolved'] < 3.5) )
#
#        #the following 2d histos are only needed before the cuts
#        cut = categories['boosted_HandW']
#        for objlist in [['Wlep','H','WhadResolved']]: #, ['Wlep','H','WhadBoosted']]:
#          for obj in itertools.combinations(itertools.combinations(objlist,2),2):
#            var1 = 'deltaR_{}_{}'.format(obj[0][0],obj[0][1])
#            var2 = 'deltaR_{}_{}'.format(obj[1][0],obj[1][1])
#            if (var1 in var) and (var2 in var):
#              hist, binsx, binsy = NUMPY_LIB.histogram2d(var[var1][cut], var[var2][cut],\
#                                                         bins=(\
#                                                               NUMPY_LIB.linspace(histogram_settings[var1][0], histogram_settings[var1][1], histogram_settings[var1][2]),\
#                                                               NUMPY_LIB.linspace(histogram_settings[var2][0], histogram_settings[var2][1], histogram_settings[var2][2]),\
#                                                              ),\
#                                                         weights=weights["nominal"][cut]\
#                                                        )
#              ret['hist2d_deltaR_{}_{}'.format(obj[0][0]+obj[0][1], obj[1][0]+obj[1][1])] =\
#                Histogram( hist, hist, (binsx[0],binsx[-1], binsy[0],binsy[-1]) )
#
#        if 'all' in cat:
#          cat=categories.keys()
##        elif not isinstance(cat, list):
##            cat = [cat]
#        for c in cat:
#            cut = categories[c]
#            cut_name = c
#
#            if p=="unsplit":
#                if "Run" in sample:
#                    name = "data" + "_" + cut_name
#                elif "process" in samples_info[sample].keys():
#                    name = samples_info[sample]["process"] + "_" + cut_name
#                else:
#                    name = sample.split('_')[0] + "_" + cut_name
#            else:
#                name = p + "_" + cut_name
#
#
#
#
##            if DNN:
##                if DNN.endswith("multiclass"):
##                    class_pred = NUMPY_LIB.argmax(DNN_pred, axis=1)
##                    for n, n_name in zip([0,1,2,3,4,5], ["ttH", "ttbb", "tt2b", "ttb", "ttcc", "ttlf"]):
##                        node = (class_pred == n)
##                        DNN_node = DNN_pred[:,n]
##                        hist_DNN = Histogram(*ha.histogram_from_vector(DNN_node[(cut & node)], weights["nominal"][(cut & node)], NUMPY_LIB.linspace(0.,1.,16)))
##                        ret["hist_{0}_DNN_{1}".format(name, n_name)] = hist_DNN
##                        hist_DNN_ROC = Histogram(*ha.histogram_from_vector(DNN_node[(cut & node)], weights["nominal"][(cut & node)], NUMPY_LIB.linspace(0.,1.,1000)))
##                        ret["hist_{0}_DNN_ROC_{1}".format(name, n_name)] = hist_DNN_ROC
##
##                else:
##                    hist_DNN = Histogram(*ha.histogram_from_vector(DNN_pred[cut], weights["nominal"][cut], NUMPY_LIB.linspace(0.,1.,16)))
##                    ret["hist_{0}_DNN".format(name)] = hist_DNN
##                    hist_DNN_ROC = Histogram(*ha.histogram_from_vector(DNN_pred[cut], weights["nominal"][cut], NUMPY_LIB.linspace(0.,1.,1000)))
##                    ret["hist_{0}_DNN_ROC".format(name)] = hist_DNN_ROC
##
##
##    #TODO: implement JECs
##
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
    parser.add_argument('--use-cuda', action='store_true', help='Use the CUDA backend')
    parser.add_argument('--from-cache', action='store_true', help='Load from cache (otherwise create it)')
    parser.add_argument('--nthreads', action='store', help='Number of CPU threads to use', type=int, default=4, required=False)
    parser.add_argument('--files-per-batch', action='store', help='Number of files to process per batch', type=int, default=1, required=False)
    parser.add_argument('--cache-location', action='store', help='Path prefix for the cache, must be writable', type=str, default=os.path.join(os.getcwd(), 'cache'))
    parser.add_argument('--outdir', action='store', help='directory to store outputs', type=str, default=os.getcwd())
    parser.add_argument('--filelist', action='store', help='List of files to load', type=str, default=None, required=False)
    parser.add_argument('--sample', action='store', help='sample name', type=str, default=None, required=True)
    parser.add_argument('--DNN', action='store', choices=['save-arrays','cmb_binary', 'cmb_multiclass', 'ffwd_binary', 'ffwd_multiclass',False], help='options for DNN evaluation / preparation', default=False)
    parser.add_argument('--categories', nargs='+', help='categories to be processed (default: sl_jge4_tge2)', default="sl_jge4_tge2")
    parser.add_argument('--path-to-model', action='store', help='path to DNN model', type=str, default=None, required=False)
    parser.add_argument('--boosted', action='store_true', help='Flag to include boosted objects', default=False)
    parser.add_argument('--year', action='store', choices=['2016', '2017', '2018'], help='Year of data/MC samples', default='2017')
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    # set CPU or GPU backend
    NUMPY_LIB, ha = choose_backend(args.use_cuda)
    lib_analysis.NUMPY_LIB, lib_analysis.ha = NUMPY_LIB, ha
    NanoAODDataset.numpy_lib = NUMPY_LIB

    if args.use_cuda:
        os.environ["HEPACCELERATE_CUDA"] = "1"
    else:
        os.environ["HEPACCELERATE_CUDA"] = "0"

    from coffea.util import USE_CUPY
    from coffea.lumi_tools import LumiMask, LumiData
    from coffea.lookup_tools import extractor

    # load definitions
    from definitions_analysis import parameters, eraDependentParameters, samples_info
    parameters.update(eraDependentParameters[args.year])

    outdir = args.outdir
    if not os.path.exists(outdir):
        print(os.getcwd())
        try:
          os.makedirs(outdir)
        except Exception as ex:
          print("An exception of type {0} occurred. Arguments:\n{1!r}".format(type(ex).__name__, ex.args))

    if "Run" in args.sample:
        is_mc = False
        lumimask = LumiMask(parameters["lumimask"])
    else:
        is_mc = True
        lumimask = None


    #define arrays to load: these are objects that will be kept together
    arrays_objects = [
        "Jet_pt", "Jet_eta", "Jet_phi", "Jet_btagDeepB", "Jet_btagDeepFlavB", "Jet_jetId", "Jet_puId", "Jet_mass",
        #"selectedPatJetsAK4PFPuppi_pt", "selectedPatJetsAK4PFPuppi_eta", "selectedPatJetsAK4PFPuppi_phi", "selectedPatJetsAK4PFPuppi_pfDeepCSVJetTags_probb", "selectedPatJetsAK4PFPuppi_pfDeepCSVJetTags_probbb", "selectedPatJetsAK4PFPuppi_jetId", "selectedPatJetsAK4PFPuppi_AK4PFPuppipileupJetIdEvaluator_fullId", "selectedPatJetsAK4PFPuppi_mass",
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_pfRelIso04_all", "Muon_tightId", "Muon_charge",
        "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass", "Electron_charge", "Electron_deltaEtaSC", "Electron_cutBased", "Electron_dz", "Electron_dxy",
    ]
    if args.boosted:
      arrays_objects += [
        "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_deepTagMD_bbvsLight", "FatJet_btagHbb", "FatJet_deepTagMD_HbbvsQCD", "FatJet_deepTagMD_ZHbbvsQCD", "FatJet_deepTagMD_TvsQCD", "FatJet_deepTag_H", "FatJet_btagDDBvL", "FatJet_deepTag_TvsQCD", "FatJet_jetId", "FatJet_mass", "FatJet_msoftdrop", "FatJet_tau1", "FatJet_tau2", "FatJet_tau3", "FatJet_tau4", "FatJet_n2b1"]

    #these are variables per event
    arrays_event = [
        "PV_npvsGood", "PV_ndof", "PV_npvs", "PV_score", "PV_x", "PV_y", "PV_z", "PV_chi2",
        "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter", "Flag_ecalBadCalibFilter",
        "MET_pt", "MET_phi", "MET_sumEt",
        "run", "luminosityBlock", "event"
    ]

    if args.year.startswith('2016'): arrays_event += [ "HLT_Ele27_WPTight_Gsf", "HLT_IsoMu24", "HLT_IsoTkMu24" ]
    elif args.year.startswith('2017'):
      arrays_event += [ "HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_IsoMu27", "HLT_IsoMu24_eta2p1" ]
      if not (args.sample.endswith('2017B') or args.sample.endswith('2017C')):
        arrays_event += ["HLT_Ele32_WPTight_Gsf"] #FIXME
    elif args.year.startswith('2018'): arrays_event += [ "HLT_Ele32_WPTight_Gsf", 'HLT_Ele28_eta2p1_WPTight_Gsf_HT150', "HLT_IsoMu24" ]

    if args.sample.startswith("TT"):
        arrays_event.append("genTtbarId")

    if is_mc:
        arrays_event += ["PV_npvsGood", "Pileup_nTrueInt", "genWeight", "nGenPart"]
        arrays_objects += [ "Jet_hadronFlavour", #"selectedPatJetsAK4PFPuppi_hadronFlavor",
                          #"GenPart_eta","GenPart_genPartIdxMother","GenPart_mass","GenPart_pdgId","GenPart_phi","GenPart_pt","GenPart_status","GenPart_statusFlags"
                         ]

    filenames = None
    if not args.filelist is None:
        filenames = [l.strip() for l in open(args.filelist).readlines()]
    else:
        filenames = args.filenames

    print("Number of files:", len(filenames))

    for fn in filenames:
        if not fn.endswith(".root"):
            print(fn)
            raise Exception("Must supply ROOT filename, but got {0}".format(fn))

    results = Results()


    for ibatch, files_in_batch in enumerate(chunks(filenames, args.files_per_batch)):
#      try:
        print(f'!!!!!!!!!!!!! loading {ibatch}: {files_in_batch}')
        #define our dataset
        structs = ["Jet", "Muon", "Electron"]#, "selectedPatJetsAK4PFPuppi"]
        if args.boosted:
          structs += ["FatJet"]#, "GenPart"]#, "MET"]
#          if is_mc:
#            structs += ['GenPart']
        dataset = NanoAODDataset(files_in_batch, arrays_objects + arrays_event, "Events", structs, arrays_event)
        dataset.get_cache_dir = lambda fn,loc=args.cache_location: os.path.join(loc, fn)

        if not args.from_cache:
            #Load data from ROOT files
            dataset.preload(nthreads=args.nthreads, verbose=True)

            #prepare the object arrays on the host or device
            dataset.make_objects()

            #print("preparing dataset cache")
            #save arrays for future use in cache
            #dataset.to_cache(verbose=True, nthreads=args.nthreads)  ###ALE: comment to run without cache


        #Optionally, load the dataset from an uncompressed format
        else:
          print("loading dataset from cache")
          dataset.from_cache(verbose=True, nthreads=args.nthreads)

        if is_mc:

            # add information needed for MC corrections
            parameters["pu_corrections_target"] = load_puhist_target(parameters["pu_corrections_file"])

            ext = extractor()
            print(parameters["corrections"])
            for corr in parameters["corrections"]:
                ext.add_weight_sets([corr])
            ext.finalize()
            evaluator = ext.make_evaluator()

        if ibatch == 0:
            print(dataset.printout())

        # in case of DNN evaluation: load model
        model = None
        if args.DNN:
            model = load_model(args.path_to_model, custom_objects=dict(itertools=itertools, mse0=mse0, mae0=mae0, r2_score0=r2_score0))

        print(args.categories)
        #### this is where the magic happens: run the main analysis
        results += dataset.analyze(analyze_data, NUMPY_LIB=NUMPY_LIB, parameters=parameters, is_mc = is_mc, lumimask=lumimask, cat=args.categories, sample=args.sample, samples_info=samples_info, boosted=args.boosted, DNN=args.DNN, DNN_model=model)
#      except Exception as ex:
#        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
#        message = template.format(type(ex).__name__, ex.args)
#        print(message)
#        print(f'!!!!!!!!!!!!!!! failed on {files_in_batch}')
#        #if is_mc:
#        #  folder = 'RunIIFall17NanoAODv5'
#        #else:
#        #  folder = 'Nano25Oct2019'
#        #with open(os.getcwd()+'/datasets/{folder}/{args.sample}.txt', 'a+') as f:
#        #with open(f'/afs/cern.ch/work/d/druini/public/hepaccelerate/datasets/{folder}/{args.sample}_fail.txt', 'a+') as f:
#        with open(args.filelist, 'a+') as f:
#          f.write(files_in_batch[0]+'\n')
#          continue

    print(results)

    #Save the results
    if not os.path.isdir(args.outdir):
      os.makedirs(args.outdir)
    results.save_json(os.path.join(outdir,"out_{}.json".format(args.sample)))

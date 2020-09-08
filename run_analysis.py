import os, glob
import argparse
import json
import numpy as np

import uproot
from uproot_methods import TLorentzVectorArray
#import hepaccelerate
from hepaccelerate.utils import Results, NanoAODDataset, Histogram, choose_backend

import itertools
#from lib_analysis import mse0,mae0,r2_score0

from definitions_analysis import histogram_settings

import lib_analysis
from lib_analysis import vertex_selection, lepton_selection, jet_selection, load_puhist_target, compute_pu_weights, compute_lepton_weights, compute_btag_weights, chunks, calculate_variable_features, select_lepton_p4, hadronic_W, get_histogram


#This function will be called for every file in the dataset
def analyze_data(data, sample, NUMPY_LIB=None, parameters={}, samples_info={}, is_mc=True, lumimask=None, cat=False, boosted=False, uncertainty=None):
    #Output structure that will be returned and added up among the files.
    #Should be relatively small.
    ret = Results()

    muons     = data["Muon"]
    electrons = data["Electron"]
    scalars   = data["eventvars"]
    jets      = data["Jet"]
    fatjets   = data["FatJet"]
    if is_mc:
      genparts  = data['GenPart']

    if uncertainty is not None:
        if len(uncertainty)!=2: raise Exception(f'Invalid uncertainty {uncertainty}')
        evUnc, objUnc = uncertainty
        if len(evUnc)%2!=0: raise Exception(f'Invalid uncertainty for events {evUnc}')
        for oldVar,newVar in zip(evUnc[::2],evUnc[1::2]):
            scalars[oldVar] = scalars[newVar]

        if len(objUnc)%3!=0: raise Exception(f'Invalid uncertainty for objects {objUnc}')
        for struct,oldBranch,newBranch in zip(objUnc[::3],objUnc[1::3],objUnc[2::3]):
            if struct=='FatJet':
                setattr(fatjets, oldBranch, getattr(fatjets,newBranch))
            elif struct=='Jet':
                setattr(jets, oldBranch, getattr(jets,newBranch))
            else:
                raise Exception(f'Problem with uncertainty on {struct}, {oldBranch}, {newBranch}')

    jets.p4 = TLorentzVectorArray.from_ptetaphim(jets.pt, jets.eta, jets.phi, jets.mass)

    if args.year=='2017':
      metstruct = 'METFixEE2017'
    else:
      metstruct = 'MET'

    METp4 = TLorentzVectorArray.from_ptetaphim(scalars[metstruct+"_pt"], 0, scalars[metstruct+"_phi"], 0)
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
    good_jets = jet_selection(jets, muons, (good_muons|veto_muons), parameters["jets"]) & jet_selection(jets, electrons, (good_electrons|veto_electrons), parameters["jets"])
#    good_jets = jet_selection(jets, muons, (veto_muons | good_muons), parameters["jets"]) & jet_selection(jets, electrons, (veto_electrons | good_electrons) , parameters["jets"])
    bjets_resolved = good_jets & (getattr(jets, parameters["btagging_algorithm"]) > parameters["btagging_WP"])
    good_fatjets = jet_selection(fatjets, muons, good_muons, parameters["fatjets"]) & jet_selection(fatjets, electrons, good_electrons, parameters["fatjets"])
#    good_fatjets = jet_selection(fatjets, muons, (veto_muons | good_muons), parameters["fatjets"]) & jet_selection(fatjets, electrons, (veto_electrons | good_electrons), parameters["fatjets"]) #FIXME remove vet_leptons

#    higgs_candidates = good_fatjets & (fatjets.pt > 250)
#    nhiggs = ha.sum_in_offsets(fatjets, higgs_candidates, mask_events, fatjets.masks["all"], NUMPY_LIB.int8)
#    indices["best_higgs_candidate"] = ha.index_in_offsets(fatjets.pt, fatjets.offsets, 1, mask_events, higgs_candidates)
#    best_higgs_candidate = NUMPY_LIB.zeros_like(higgs_candidates)
#    best_higgs_candidate[ (fatjets.offsets[:-1] + indices["best_higgs_candidate"])[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] ] = True
#    best_higgs_candidate[ (fatjets.offsets[:-1] + indices["best_higgs_candidate"])[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] ] &= nhiggs.astype(NUMPY_LIB.bool)[NUMPY_LIB.where( fatjets.offsets<len(best_higgs_candidate) )] # to avoid removing the leading fatjet in events with no higgs candidate

    good_jets_nohiggs = good_jets & ha.mask_deltar_first(jets, good_jets, fatjets, good_fatjets, 1.2, indices['leading'])
    bjets = good_jets_nohiggs & (getattr(jets, parameters["btagging_algorithm"]) > parameters["btagging_WP"])
    nonbjets = good_jets_nohiggs & (getattr(jets, parameters["btagging_algorithm"]) < parameters["btagging_WP"])

    # apply basic event selection -> individual categories cut later
    nleps          = NUMPY_LIB.add(ha.sum_in_offsets(muons, good_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, good_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
    lepton_veto    = NUMPY_LIB.add(ha.sum_in_offsets(muons, veto_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, veto_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
    njets          = ha.sum_in_offsets(jets, nonbjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    ngoodjets      = ha.sum_in_offsets(jets, good_jets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags          = ha.sum_in_offsets(jets, bjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags_resolved = ha.sum_in_offsets(jets, bjets_resolved, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    nfatjets       = ha.sum_in_offsets(fatjets, good_fatjets, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)
    #nhiggs = ha.sum_in_offsets(fatjets, higgs_candidates, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)

    # for reference, this is the selection for the resolved analysis
    mask_events_res = mask_events & (nleps == 1) & (lepton_veto == 0) & (ngoodjets >= 4) & (btags_resolved > 2) & (scalars[metstruct+"_pt"] > 20)
    # apply basic event selection
    #mask_events_higgs = mask_events & (nleps == 1) & (scalars[metstruct+"_pt"] > 20) & (nhiggs > 0) & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4) & (btags >=2) ) & (lepton_veto == 0)
    mask_events = mask_events & (nleps == 1) & (lepton_veto == 0) & (scalars[metstruct+"_pt"] > parameters['met']) & (nfatjets > 0) & (btags >= parameters['btags']) # & (btags_resolved < 3)# & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4)  )

############# calculate weights for MC samples
    weights = {}
    weights['ones'] = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.float32)
    weights["nominal"] = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.float32)

    if is_mc:
        weights["nominal"] = weights["nominal"] * scalars["genWeight"] * parameters["lumi"] * samples_info[sample]["XS"] / samples_info[sample]["ngen_weight"][args.year]

        # pu corrections
        if 'puWeight' in scalars:
            pu_weights = scalars['puWeight']
        else:
    #        pu_weights = compute_pu_weights(parameters["pu_corrections_target"], weights["nominal"], scalars["Pileup_nTrueInt"], scalars["PV_npvsGood"])
            pu_weights = compute_pu_weights(parameters["pu_corrections_target"], weights["nominal"], scalars["Pileup_nTrueInt"], scalars["Pileup_nTrueInt"])
        weights["nominal"] = weights["nominal"] * pu_weights

        # lepton SF corrections
        #electron_weights = compute_lepton_weights(electrons, electrons.pt, (electrons.deltaEtaSC + electrons.eta), mask_events, good_electrons, evaluator, ["el_triggerSF", "el_recoSF", "el_idSF"])
        #muon_weights = compute_lepton_weights(muons, muons.pt, NUMPY_LIB.abs(muons.eta), mask_events, good_muons, evaluator, ["mu_triggerSF", "mu_isoSF", "mu_idSF"])
        #weights["nominal"] = weights["nominal"] * muon_weights * electron_weights

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

    leading_fatjet_rho    = NUMPY_LIB.zeros_like(leading_lepton_pt)
    leading_fatjet_rho[mask_events] = NUMPY_LIB.log( leading_fatjet_SDmass[mask_events]**2 / leading_fatjet_pt[mask_events]**2 )

    lead_lep_p4        = select_lepton_p4(muons, good_muons, electrons, good_electrons, indices["leading"], mask_events)
    leading_fatjet_phi = ha.get_in_offsets(fatjets.phi, fatjets.offsets, indices['leading'], mask_events, good_fatjets)
    deltaRHiggsLepton  = ha.calc_dr(lead_lep_p4.phi, lead_lep_p4.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events)

############# masks for different selections
    mask_events = {
      'resolved' : mask_events_res,
      'basic'    : mask_events
    }
    mask_events['2J']   = mask_events['basic'] & (njets>1)

    #Ws reconstruction
    pznu = ha.METzCalculator(lead_lep_p4, METp4, mask_events['2J'])
    neutrinop4 = TLorentzVectorArray.from_cartesian(METp4.x, METp4.y, pznu, NUMPY_LIB.sqrt( METp4.x**2 + METp4.y**2 + pznu**2 ))
    lepW = lead_lep_p4 + neutrinop4

    hadW = hadronic_W(jets, nonbjets, lepW, mask_events['2J'])

    mask_events['2J2W'] = mask_events['2J'] & (hadW.mass>parameters['W']['min_mass']) & (hadW.mass<parameters['W']['max_mass']) & (lepW.mass>parameters['W']['min_mass']) & (lepW.mass<parameters['W']['max_mass'])

    #deltaR between objects
    deltaRlepWHiggs = ha.calc_dr(lepW.phi, lepW.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events['2J2W'])
    deltaRhadWHiggs = ha.calc_dr(hadW.phi, hadW.eta, leading_fatjet_phi, leading_fatjet_eta, mask_events['2J2W'])

#    mask_events['2J2WdeltaR'] = mask_events['2J2W'] & (deltaRlepWHiggs>1.5) & (deltaRhadWHiggs>1.5) & (deltaRlepWHiggs<4) & (deltaRhadWHiggs<4)
    mask_events['2J2WdeltaR'] = mask_events['2J2W'] & (deltaRlepWHiggs>1) & (deltaRhadWHiggs>1)# & (deltaRlepWHiggs<4) & (deltaRhadWHiggs<4)

    #boosted Higgs
    leading_fatjet_tau1     = ha.get_in_offsets(fatjets.tau1, fatjets.offsets, indices['leading'], mask_events['2J2WdeltaR'], good_fatjets)
    leading_fatjet_tau2     = ha.get_in_offsets(fatjets.tau2, fatjets.offsets, indices['leading'], mask_events['2J2WdeltaR'], good_fatjets)
    leading_fatjet_tau21    = NUMPY_LIB.divide(leading_fatjet_tau2, leading_fatjet_tau1)
    ### tau21DDT defined as in https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging#tau21DDT_0_43_HP_0_43_tau21DDT_0
#    leading_fatjet_tau21DDT = NUMPY_LIB.zeros_like(leading_fatjet_tau21)
#    if args.year=='2016':
#        leading_fatjet_tau21DDT[mask_events['2J2WdeltaR']] = leading_fatjet_tau21[mask_events['2J2WdeltaR']] * 0.063 * NUMPY_LIB.log(leading_fatjet_SDmass[mask_events['2J2WdeltaR']]**2 / leading_fatjet_pt[mask_events['2J2WdeltaR']])
#    elif args.year=='2017':
#        leading_fatjet_tau21DDT[mask_events['2J2WdeltaR']] = leading_fatjet_tau21[mask_events['2J2WdeltaR']] * 0.080 * NUMPY_LIB.log(leading_fatjet_SDmass[mask_events['2J2WdeltaR']]**2 / leading_fatjet_pt[mask_events['2J2WdeltaR']])
#    else:
#        leading_fatjet_tau21DDT[mask_events['2J2WdeltaR']] = leading_fatjet_tau21[mask_events['2J2WdeltaR']] * 0.082 * NUMPY_LIB.log(leading_fatjet_SDmass[mask_events['2J2WdeltaR']]**2 / leading_fatjet_pt[mask_events['2J2WdeltaR']])

#    mask_events['2J2WdeltaRTau21']    = mask_events['2J2WdeltaR'] & (leading_fatjet_tau21<parameters["fatjets"]["tau21cut"][args.year])
#    mask_events['2J2WdeltaRTau21DDT'] = mask_events['2J2WdeltaR'] & (leading_fatjet_tau21<parameters["fatjets"]["tau21DDTcut"][args.year])

    leading_fatjet_Hbb = ha.get_in_offsets(getattr(fatjets, parameters["bbtagging_algorithm"]), fatjets.offsets, indices['leading'], mask_events['2J2WdeltaR'], good_fatjets)
    #mask_events['2J2WdeltaRTau21_Pass'] = mask_events['2J2WdeltaRTau21'] & (leading_fatjet_Hbb>parameters['bbtagging_WP'])
    #mask_events['2J2WdeltaRTau21_Fail'] = mask_events['2J2WdeltaRTau21'] & (leading_fatjet_Hbb<=parameters['bbtagging_WP'])
    for m in ['2J2WdeltaR']:#, '2J2WdeltaRTau21']:#, '2J2WdeltaRTau21DDT']:
        mask_events[f'{m}_Pass'] = mask_events[m] & (leading_fatjet_Hbb>parameters['bbtagging_WP'])
        mask_events[f'{m}_Fail'] = mask_events[m] & (leading_fatjet_Hbb<=parameters['bbtagging_WP'])

    #mask_events['overlap'] = mask_events['2J2WdeltaR'] & mask_events['resolved']
    #mask_events['overlap'] = mask_events['2J2WdeltaR_Pass'] & mask_events['resolved']

############# overlap study
    for m in mask_events.copy():
        if m=='resolved': continue
        mask_events[m+'_orthogonal'] = mask_events[m] & (btags_resolved < 3)
        mask_events[m+'_overlap']    = mask_events[m] & mask_events['resolved']
    for mn,m in mask_events.items():
        ret['nevts_'+mn] = Histogram([sum(weights['nominal'][m])], 0,0)

    vars2d = {
            'ngoodjets' : ngoodjets,
            'njets'     : njets
            }
    for mn,m in mask_events.items():
        if 'overlap' in mn:
            #hist, binsx, binsy = NUMPY_LIB.histogram2d(njets[m], btags_resolved[m],\
            #        bins=(\
            #            NUMPY_LIB.linspace(*histogram_settings['njets']),\
            #            NUMPY_LIB.linspace(*histogram_settings['btags_resolved']),\
            #        ),\
            #        weights=weights["nominal"][m]\
            #        )
            #ret[f'hist2d_njetsVSbtags_{mn}'] = Histogram( hist, hist, (binsx[0],binsx[-1], binsy[0],binsy[-1]) )
            for vn,v in vars2d.items():
                hist, binsx, binsy = NUMPY_LIB.histogram2d(v[m], btags_resolved[m],\
                        bins=(\
                        NUMPY_LIB.linspace(*histogram_settings[vn]),\
                        NUMPY_LIB.linspace(*histogram_settings['btags_resolved']),\
                        ),\
                        weights=weights["nominal"][m]\
                        )
                ret[f'hist2d_{vn}VSbtags_{mn}'] = Histogram( hist, hist, (*histogram_settings[vn],*histogram_settings['btags_resolved']) )

############# histograms
    vars_to_plot = {
      'nleps'             : nleps,
      'njets'             : njets,
      'ngoodjets'         : ngoodjets,
      'btags'             : btags,
      'btags_resolved'    : btags_resolved,
      'nfatjets'          : nfatjets,
      'met'               : scalars[metstruct+'_pt'],
      'leading_jet_pt'    : leading_jet_pt,
      'leading_jet_eta'   : leading_jet_eta,
      'leadAK8JetMass'    : leading_fatjet_SDmass,
      'leadAK8JetPt'      : leading_fatjet_pt,
      'leadAK8JetEta'     : leading_fatjet_eta,
      'leadAK8JetHbb'     : leading_fatjet_Hbb,
      'leadAK8JetTau21'   : leading_fatjet_tau21,
      'leadAK8JetRho'     : leading_fatjet_rho,
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

    #var_name, var = 'leadAK8JetMass', leading_fatjet_SDmass
    vars_split = ['leadAK8JetMass', 'leadAK8JetRho']
    ptbins = NUMPY_LIB.append( NUMPY_LIB.arange(250,600,50), [600, 1000, 5000] )
    for var_name in vars_split:
      var = vars_to_plot[var_name]
      for ipt in range( len(ptbins)-1 ):
        for m in ['2J2WdeltaR']:#, '2J2WdeltaRTau21']:#, '2J2WdeltaRTau21DDT']:
          for r in ['Pass','Fail']:
            for o in ['','_orthogonal']:
              mask_name = f'{m}_{r}{o}'
              if not mask_name in mask_events: continue
              mask = mask_events[mask_name] & (leading_fatjet_pt>ptbins[ipt]) & (leading_fatjet_pt<ptbins[ipt+1])
              ret[f'hist_{var_name}_{mask_name}_pt{ptbins[ipt]}to{ptbins[ipt+1]}'] = get_histogram( var[mask], weights['nominal'][mask], NUMPY_LIB.linspace( *histogram_settings[var_name] ) )

    weight_names = {'' : 'nominal', '_NoWeights' : 'ones'}
    for weight_name, w in weight_names.items():
      if w=='ones': continue
      #ret[f'nevts_overlap{weight_name}'] = Histogram( [sum(weights[w]), sum(weights[w][mask_events['2J2WdeltaR']]), sum(weights[w][mask_events['resolved']]), sum(weights[w][mask_events['overlap']])], 0,0 )
      for mask_name, mask in mask_events.items():
#        with open(f'/afs/cern.ch/work/d/druini/public/hepaccelerate/tests/events_pass_selection_{sample}_{mask_name}.txt','a+') as f:
#          for nevt, run, lumiBlock in zip(scalars['event'][mask], scalars['run'][mask], scalars['luminosityBlock']):
#            f.write(f'{nevt}, {run}, {lumiBlock}\n')
        for var_name, var in vars_to_plot.items():
          if (not is_mc) and ('Pass' in mask_name) and (var_name=='leadAK8JetMass') : continue
          try:
            ret[f'hist_{var_name}_{mask_name+weight_name}'] = get_histogram( var[mask], weights[w][mask], NUMPY_LIB.linspace( *histogram_settings[var_name] ) )
          except KeyError:
            print(f'!!!!!!!!!!!!!!!!!!!!!!!! Please add variable {var_name} to the histogram settings')

############# genPart study: where are the b quarks?
    if sample=='ttHTobb':
        genH = (abs(genparts.pdgId)==25) & (genparts.status==22)
        nH   = ha.sum_in_offsets(genparts,genH,NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.bool),genparts.masks["all"], NUMPY_LIB.int8)
        if not NUMPY_LIB.all(nH==1):
          import pdb
          pdb.set_trace()
        for mn,m in mask_events.items():
            genH_pt = ha.get_in_offsets(genparts.pt, genparts.offsets, indices['leading'], m, genH)
            ret[f'hist_genH_pt_{mn}'] = get_histogram(genH_pt[m], weights['nominal'][m], NUMPY_LIB.linspace(*histogram_settings['leading_jet_pt']))
            genb = {}
            genb['top'] = ha.genPart_from_mother(genparts, 5, 6, m)
            genb['H']   = ha.genPart_from_mother(genparts, 5, 25, m)
            for mom,genmask in genb.items():
                genb_vars = {}
                genb_vars['pt']  = genparts.pt[genmask]
                genb_vars['phi'] = genparts.phi[genmask]
                genb_vars['eta'] = genparts.eta[genmask]
                nevs = NUMPY_LIB.sum(m)
                dr_b1fatjet = ha.calc_dr(genb_vars['phi'][::2], genb_vars['eta'][::2],leading_fatjet_phi[m],leading_fatjet_eta[m],NUMPY_LIB.ones(nevs))
                dr_b2fatjet = ha.calc_dr(genb_vars['phi'][1::2], genb_vars['eta'][1::2],leading_fatjet_phi[m],leading_fatjet_eta[m],NUMPY_LIB.ones(nevs))
                for weight_name, w in weight_names.items():
                    if w=='ones': continue
                    #ret[f'hist_dr_b1fatjet_{mn+weight_name}'] = get_histogram( dr_b1fatjet, weights[w][m], NUMPY_LIB.linspace(0,10,101) )
                    #ret[f'hist_dr_b2fatjet_{mn+weight_name}'] = get_histogram( dr_b2fatjet, weights[w][m], NUMPY_LIB.linspace(0,10,101) )
                    ret[f'hist_dr_genbfrom{mom}_fatjet_{mn+weight_name}'] = get_histogram( dr_b1fatjet, weights[w][m], NUMPY_LIB.linspace(0,10,101) ) + get_histogram( dr_b2fatjet, weights[w][m], NUMPY_LIB.linspace(0,10,101) )
                    for var in ['pt','eta']:
                        ret[f'hist_genbfrom{mom}_{var}_{mn+weight_name}'] = get_histogram( genb_vars[var][::2], weights[w][m], NUMPY_LIB.linspace(*histogram_settings[f'leading_jet_{var}']) ) + get_histogram( genb_vars[var][1::2], weights[w][m], NUMPY_LIB.linspace(*histogram_settings[f'leading_jet_{var}']) )

    ### next lines are to write event numbers of very high pt events
    #mask = mask_events['2J2WdeltaR'] & (leading_fatjet_pt>1500)
    #if 'Single' in sample:
    #  with open('/afs/cern.ch/work/d/druini/public/hepaccelerate/highptEvents.txt','a+') as f:
    #    for nevt, run, lumiBlock in zip(scalars['event'][mask], scalars['run'][mask], scalars['luminosityBlock']):
    #      f.write(f'{sample}, {nevt}, {run}, {lumiBlock}\n')

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
#    print('met', scalars[metstruct+'_pt'][mask])
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
#      print('jets btag', getattr(jets, parameters["btagging_algorithm"])[start:stop])
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
#            ("", fatjets, higgs_candidates, "best_higgs_candidate", ["pt", "msoftdrop", "tau21", parameters["bbtagging_algorithm"]]),
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
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs a simple array-based analysis')
    parser.add_argument('--use-cuda', action='store_true', help='Use the CUDA backend')
    parser.add_argument('--from-cache', action='store_true', help='Load from cache (otherwise create it)')
    parser.add_argument('--nthreads', action='store', help='Number of CPU threads to use', type=int, default=4, required=False)
    parser.add_argument('--files-per-batch', action='store', help='Number of files to process per batch', type=int, default=1, required=False)
    parser.add_argument('--cache-location', action='store', help='Path prefix for the cache, must be writable', type=str, default=os.path.join(os.getcwd(), 'cache'))
    parser.add_argument('--outdir', action='store', help='directory to store outputs', type=str, default=os.getcwd())
    parser.add_argument('--outtag', action='store', help='outtag added to output file', type=str, default="")
    parser.add_argument('--version', action='store', help='tag added to the output directory', type=str, default='')
    parser.add_argument('--filelist', action='store', help='List of files to load', type=str, default=None, required=False)
    parser.add_argument('--sample', action='store', help='sample name', type=str, default=None, required=True)
    parser.add_argument('--categories', nargs='+', help='categories to be processed (default: sl_jge4_tge2)', default="sl_jge4_tge2")
    parser.add_argument('--boosted', action='store_true', help='Flag to include boosted objects', default=False)
    parser.add_argument('--year', action='store', choices=['2016', '2017', '2018'], help='Year of data/MC samples', default='2017')
    parser.add_argument('--parameters', nargs='+', help='change default parameters, syntax: name value, eg --parameters met 40 bbtagging_algorithm btagDDBvL', default=None)
    parser.add_argument('--uncertainty', nargs='+', help='replace a branch with another, syntax: struct oldBranch newBranch, e.g. FatJet pt pt_jerUp', default=None)
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
    if args.parameters is not None:
      if len(args.parameters)%2 is not 0:
        raise Exception('incomplete parameters specified, quitting.')
      for p,v in zip(args.parameters[::2], args.parameters[1::2]):
        try: parameters[p] = type(parameters[p])(v) #convert the string v to the type of the parameter already in the dictionary
        except: print(f'invalid parameter specified: {p} {v}')

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
        "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_deepTagMD_bbvsLight", "FatJet_btagHbb", "FatJet_deepTagMD_HbbvsQCD", "FatJet_deepTagMD_ZHbbvsQCD", "FatJet_deepTagMD_TvsQCD", "FatJet_deepTag_H", "FatJet_btagDDBvL", "FatJet_btagDDBvL_noMD", "FatJet_deepTag_TvsQCD", "FatJet_jetId", "FatJet_mass", "FatJet_msoftdrop", "FatJet_tau1", "FatJet_tau2", "FatJet_tau3", "FatJet_tau4", "FatJet_n2b1"]

    #these are variables per event
    arrays_event = [
        "PV_npvsGood", "PV_ndof", "PV_npvs", "PV_score", "PV_x", "PV_y", "PV_z", "PV_chi2",
        "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_eeBadScFilter", "Flag_ecalBadCalibFilter",
        "run", "luminosityBlock", "event"
    ]

    metstruct = 'METFixEE2017' if args.year=='2017' else 'MET'
    arrays_event += [f'{metstruct}_{var}' for var in ['pt','phi','sumEt']]
    #if args.year.startswith('2017'): arrays_event += ["METFixEE2017_pt", "METFixEE2017_phi", "METFixEE2017_sumEt"]
    #else: arrays_event += ["MET_pt", "MET_phi", "MET_sumEt"]

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
                           "GenPart_eta","GenPart_genPartIdxMother","GenPart_mass","GenPart_pdgId","GenPart_phi","GenPart_pt","GenPart_status","GenPart_statusFlags"
                         ]

#    if len(args.uncertainty)%3 is not 0:
#        raise Exception('invalid uncertainty argument, quitting.')
#    arrays_objects += [f'{struct}_{newBranch}' for struct,newBranch in zip(args.uncertainty[::3], args.uncertainty[2::3])]
    arrays_objects += [f'FatJet_{var}_{unc}{ud}' for var in ['pt','mass'] for unc in ['jer','jesTotal'] for ud in ['Up','Down']]
    arrays_objects += [f'FatJet_msoftdrop_{unc}{ud}' for unc in ['jmr','jms'] for ud in ['Up','Down']]
    arrays_objects += [f'FatJet_{var}_{corr}' for var in ['msoftdrop','pt','mass'] for corr in ['raw','nom']]
    arrays_objects += [f'Jet_pt_nom']
    arrays_event   += [f'{metstruct}_{var}_nom' for var in ['pt','phi']]
    arrays_event   += [f'puWeight{var}' for var in ['','Up','Down']]

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

    #results = Results()
#    WPs_DAK8 = [0.8695]#, 0.9795]0.5845, 
#    WPs_DDB  = [0.86, 0.89, 0.91]#, 0.92]0.7, 
#    bbtags   = {'deepTagMD_bbvsLight': WPs_DAK8, 'btagDDBvL_noMD': WPs_DDB, 'btagDDBvL': WPs_DDB} #'deepTag_H': WPs_DAK8, 'btagDDBvL': WPs_DDB, 'btagDDBvL_noMD': WPs_DDB}
#    pars     = {f'met{met}_{bbAlg}0{str(bbWP).split(".")[-1]}' : (met,bbAlg,bbWP) for met in [20] for bbAlg,bbWPlist in bbtags.items() for bbWP in bbWPlist}
##    for p in pars.copy():
##        pars[f'{p}_1btag'] = pars[p] + (1,)
##        pars[p]            = pars[p] + (0,)
#
#    results  = {p : Results() for p in pars}

    uncertainties = {
            'jerUp'        : [[],['FatJet','pt','pt_jerUp','FatJet','mass','mass_jerUp']],
            'jerDown'      : [[],['FatJet','pt','pt_jerDown','FatJet','mass','mass_jerDown']],
            'jesTotalUp'   : [[],['FatJet','pt','pt_jesTotalUp','FatJet','mass','mass_jesTotalUp']],
            'jesTotalDown' : [[],['FatJet','pt','pt_jesTotalDown','FatJet','mass','mass_jesTotalDown']],
            'jmrUp'        : [[],['FatJet','msoftdrop','msoftdrop_jmrUp']],
            'jmrDown'      : [[],['FatJet','msoftdrop','msoftdrop_jmrDown']],
            'jmsUp'        : [[],['FatJet','msoftdrop','msoftdrop_jmsUp']],
            'jmsDown'      : [[],['FatJet','msoftdrop','msoftdrop_jmsDown']],
            'msd_nom'      : [[],['FatJet','msoftdrop','msoftdrop_nom']],
            #'msd_raw'      : ['FatJet','msoftdrop','msoftdrop_raw'],
            #'msd_nanoAOD'  : ['FatJet','msoftdrop','msoftdrop']
            'puWeightUp'   : [['puWeight','puWeightUp'],[]],
            'puWeightDown' : [['puWeight','puWeightDown'],[]],
            }
    for u in uncertainties.values():
        u[1] += ['Jet','pt','pt_nom']
        if 'msoftdrop' not in u[1]:
            u[1] += ['FatJet','msoftdrop','msoftdrop_nom']
        if not (('FatJet' in u[1]) and ('pt' in u[1])):
            u[1] += ['FatJet','pt','pt_nom',
                'FatJet','mass','mass_nom'
                ]

    results = {u : Results() for u in uncertainties}

    for ibatch, files_in_batch in enumerate(chunks(filenames, args.files_per_batch)):
        print(f'!!!!!!!!!!!!! loading {ibatch}: {files_in_batch}')
        #define our dataset
        structs = ["Jet", "Muon", "Electron"]#, "selectedPatJetsAK4PFPuppi"]
        if is_mc:
          structs += ['GenPart']
        if args.boosted:
          structs += ["FatJet"]#, "MET"]
#          if is_mc:
#            structs += ['GenPart']
        dataset = NanoAODDataset(files_in_batch, arrays_objects + arrays_event, "Events", structs, arrays_event)
        dataset.get_cache_dir = lambda fn,loc=args.cache_location: os.path.join(loc, fn)

        if not args.from_cache:
            #Load data from ROOT files
            dataset.preload(nthreads=args.nthreads, verbose=True)

            #prepare the object arrays on the host or device
            dataset.make_objects()

            #save arrays for future use in cache
#            print("preparing dataset cache")
#            dataset.to_cache(verbose=True, nthreads=args.nthreads)  ###ALE: comment to run without cache


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

#        for p in pars:
#          parameters['met'], parameters['bbtagging_algorithm'], parameters['bbtagging_WP'] = pars[p] #, parameters['btags']
        for un,u in uncertainties.items():
        #### this is where the magic happens: run the main analysis
          results[un] += dataset.analyze(analyze_data, NUMPY_LIB=NUMPY_LIB, parameters=parameters, is_mc = is_mc, lumimask=lumimask, cat=args.categories, sample=args.sample, samples_info=samples_info, boosted=args.boosted, uncertainty=u)

    #print(results)

    #Save the results
    outdir = args.outdir
    if args.version!='':
      outdir = os.path.join(outdir,args.version)
    auxdir = f"met{parameters['met']}_{parameters['bbtagging_algorithm']}0{str(parameters['bbtagging_WP']).split('.')[-1]}"
    outdir = os.path.join(outdir,auxdir)
    if not os.path.exists(outdir):
      os.makedirs(outdir)

    for r in results:
      results[r].save_json(os.path.join(outdir,f"out_{args.sample}_{r}{args.outtag}.json"))

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
def analyze_data(data, sample, NUMPY_LIB=None, parameters={}, samples_info={}, is_mc=True, lumimask=None, cat=False, boosted=False):

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
#    bjets = good_jets & (getattr(jets, parameters["btagging_algorithm"]) > parameters["btagging_WP"])
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
    nleps =  NUMPY_LIB.add(ha.sum_in_offsets(muons, good_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, good_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
#    lepton_veto = NUMPY_LIB.add(ha.sum_in_offsets(muons, veto_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, veto_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
    njets = ha.sum_in_offsets(jets, nonbjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags = ha.sum_in_offsets(jets, bjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    nfatjets = ha.sum_in_offsets(fatjets, good_fatjets, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)
    #nhiggs = ha.sum_in_offsets(fatjets, higgs_candidates, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)

    # for reference, this is the selection for the resolved analysis
    mask_events_res = mask_events & (nleps == 1) & (lepton_veto == 0) & (njets >= 4) & (btags >=2) & met
    # apply basic event selection
    #mask_events_higgs = mask_events & (nleps == 1) & (scalars["MET_pt"] > 20) & (nhiggs > 0) & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4) & (btags >=2) ) & (lepton_veto == 0)
    mask_events = mask_events & (nleps == 1) & (scalars["MET_pt"] > parameters['met']) & (nfatjets > 0) & (btags >= parameters['btags'])# & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4)  ) & (lepton_veto == 0)

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
      'basic' : mask_events
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

    #leading_fatjet_Hbb = ha.get_in_offsets(getattr(fatjets, parameters["bbtagging_algorithm"]), fatjets.offsets, indices['leading'], mask_events['2J2WdeltaRTau21'], good_fatjets)

    assert( len(mask_events['2J2WdeltaR'])==len(mask_events_res) )
    #return (total number of events, number of boosted events, number of resolved events, number of events in overlap)
    return (len(mask_events['2J2WdeltaR']), sum(mask_events['2J2WdeltaR']), sum(mask_events_res), sum(mask_events['2J2WdeltaR']&mask_events_res))

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

    #results = Results()
    WPs_DAK8 = [0.8695]#, 0.9795]0.5845, 
    WPs_DDB  = [0.86, 0.89, 0.91]#, 0.92]0.7, 
    bbtags   = {'btagDDBvL_noMD': WPs_DDB, 'btagDDBvL': WPs_DDB} #'deepTagMD_bbvsLight': WPs_DAK8, 'deepTag_H': WPs_DAK8, 'btagDDBvL': WPs_DDB, 'btagDDBvL_noMD': WPs_DDB}
    pars     = {f'met{met}_{bbAlg}0{str(bbWP).split(".")[-1]}' : (met,bbAlg,bbWP) for met in [20] for bbAlg,bbWPlist in bbtags.items() for bbWP in bbWPlist}
    #for p in pars.copy():
    #    pars[f'{p}_1btag'] = pars[p] + (1,)
    #    pars[p]            = pars[p] + (0,)

    results  = {}

    for ibatch, files_in_batch in enumerate(chunks(filenames, args.files_per_batch)):
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

        for p in pars:
          parameters['met'], parameters['bbtagging_algorithm'], parameters['bbtagging_WP'] = pars[p]
        #### this is where the magic happens: run the main analysis
          results[p] = dataset.analyze(analyze_data, NUMPY_LIB=NUMPY_LIB, parameters=parameters, is_mc = is_mc, lumimask=lumimask, cat=args.categories, sample=args.sample, samples_info=samples_info, boosted=args.boosted)

    #print(results)

    #Save the results
    outdir = args.outdir
    if args.version!='':
      outdir = os.path.join(outdir,args.version)
    if not os.path.exists(outdir):
      os.makedirs(outdir)

    with open(os.path.join(outdir,f'out_{args.sample}{args.outtag}.json'), 'w+') as f:
        json.dump(results,f)

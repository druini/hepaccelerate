import os, glob
import argparse
import json
import numpy as np

import uproot
from uproot_methods import TLorentzVectorArray
#import hepaccelerate
from hepaccelerate.utils import Results, NanoAODDataset, Histogram, choose_backend

#import itertools
#from lib_analysis import mse0,mae0,r2_score0

from definitions_analysis import histogram_settings

import lib_analysis
from lib_analysis import vertex_selection, lepton_selection, jet_selection, load_puhist_target, compute_pu_weights, compute_lepton_weights, compute_btag_weights, chunks, calculate_variable_features, select_lepton_p4, hadronic_W, get_histogram

from pdb import set_trace
import sys

#This function will be called for every file in the dataset
def analyze_data(data, sample, NUMPY_LIB=None, parameters={}, samples_info={}, is_mc=True, lumimask=None, cat=False, boosted=False, uncertainty=None, uncertaintyName=None, parametersName=None, extraCorrection=None):
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

    if args.year=='2017':
      metstruct = 'METFixEE2017'
    else:
      metstruct = 'MET'

    if is_mc and uncertaintyName.startswith('jes') and not sample.startswith('ttH'):
        jesMerged = {
                'Absolute'                    : ['AbsoluteMPFBias','AbsoluteScale','Fragmentation','PileUpDataMC','PileUpPtRef','RelativeFSR','SinglePionECAL','SinglePionHCAL'],
                f'Absolute_{args.year}'       : ['AbsoluteStat','RelativeStatFSR','TimePtEta'],
                'BBEC1'                       : ['PileUpPtBB','PileUpPtEC1','RelativePtBB'],
                'EC2'                         : ['PileUpPtEC2'],
                'HF'                          : ['PileUpPtHF','RelativeJERHF','RelativePtHF'],
                f'BBEC1_{args.year}'          : ['RelativeJEREC1','RelativePtEC1','RelativeStatEC'],
                f'EC2_{args.year}'            : ['RelativePtEC2'],
                f'RelativeSample_{args.year}' : ['RelativeSample'],
                f'HF_{args.year}'             : ['RelativeStatHF']
                }
        variation = 'Up' if uncertaintyName.endswith('Up') else 'Down'
        jesName = uncertaintyName.split('jes')[-1].split(variation)[0]
        if not hasattr(jets,f'pt_jes{jesName}{variation}'):
            def sumInQuadrature(arr):
                if len(arr)==1:
                    return arr
                else:
                    return NUMPY_LIB.sqrt(NUMPY_LIB.sum(NUMPY_LIB.square(NUMPY_LIB.array(arr,dtype=NUMPY_LIB.float32)), axis=0))

            for struct in [fatjets,jets]:
                for var in ['pt','mass']:
                    c = sumInQuadrature([getattr(struct,f'{var}_nom')-getattr(struct,f'{var}_jes{jesNameSplit}{variation}') for jesNameSplit in jesMerged[jesName]])
                    setattr(struct, f'{var}_jes{jesName}{variation}', getattr(struct, f'{var}_nom')+(c if variation=='Up' else -c) )

            c = sumInQuadrature([getattr(fatjets,'msoftdrop_nom')-getattr(fatjets,f'msoftdrop_jes{jesNameSplit}{variation}') for jesNameSplit in jesMerged[jesName]])
            setattr(fatjets, f'msoftdrop_jes{jesName}{variation}', fatjets.msoftdrop_nom+(c if variation=='Up' else -c) )
            for var in ['pt','phi']:
                c = sumInQuadrature([scalars[f'{metstruct}_{var}_jer']-scalars[f'{metstruct}_{var}_jes{jesNameSplit}{variation}'] for jesNameSplit in jesMerged[jesName]])
                scalars[f'{metstruct}_{var}_jes{jesName}{variation}'] = scalars[f'{metstruct}_{var}_nom']+(c if variation=='Up' else -c)
            set_trace()
    if uncertainty is not None:
        if len(uncertainty)!=2: raise Exception(f'Invalid uncertainty {uncertainty}')
        evUnc, objUnc = uncertainty
        if len(evUnc)%2!=0: raise Exception(f'Invalid uncertainty for events {evUnc}')
        for oldVar,newVar in zip(evUnc[::2],evUnc[1::2]):
            scalars[oldVar] = scalars[newVar].copy()

        if len(objUnc)%3!=0: raise Exception(f'Invalid uncertainty for objects {objUnc}')
        for struct,oldBranch,newBranch in zip(objUnc[::3],objUnc[1::3],objUnc[2::3]):
            if struct=='FatJet':
              if not oldBranch=='bbtagSF_DDBvL_M1':
                setattr(fatjets, oldBranch, getattr(fatjets,newBranch).copy())
            elif struct=='Jet':
                setattr(jets, oldBranch, getattr(jets,newBranch).copy())
            else:
                raise Exception(f'Problem with uncertainty on {struct}, {oldBranch}, {newBranch}')
    
    if extraCorrection is not None:
      for e in extraCorrection:
          fatjets.msoftdrop /= getattr(fatjets, f'msoftdrop_corr_{e}')
    jets.p4 = TLorentzVectorArray.from_ptetaphim(jets.pt, jets.eta, jets.phi, jets.mass)

    METp4 = TLorentzVectorArray.from_ptetaphim(scalars[metstruct+"_pt"], 0, scalars[metstruct+"_phi"], 0)
    nEvents = muons.numevents()

    indices = {
        "leading"    : NUMPY_LIB.zeros(nEvents, dtype=NUMPY_LIB.int32),
        "subleading" : NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.int32)
        }

    mask_events = NUMPY_LIB.ones(nEvents, dtype=NUMPY_LIB.bool)

    # apply event cleaning and  PV selection
    flags = [
        "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"]#, "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter"]
    if not is_mc:
        flags.append("Flag_eeBadScFilter")
    for flag in flags:
        mask_events = mask_events & scalars[flag]
    mask_events = mask_events & (scalars["PV_npvsGood"]>0)

    #in case of data: check if event is in golden lumi file
    if not is_mc and not (lumimask is None):
        mask_lumi = lumimask(scalars["run"], scalars["luminosityBlock"])
        mask_events = mask_events & mask_lumi

    # apply object selection for muons, electrons, jets
    good_muons, veto_muons = lepton_selection(muons, parameters["muons"], args.year)
    good_electrons, veto_electrons = lepton_selection(electrons, parameters["electrons"], args.year)
    good_jets = jet_selection(jets, muons, (good_muons|veto_muons), parameters["jets"]) & jet_selection(jets, electrons, (good_electrons|veto_electrons), parameters["jets"])
#    good_jets = jet_selection(jets, muons, (veto_muons | good_muons), parameters["jets"]) & jet_selection(jets, electrons, (veto_electrons | good_electrons) , parameters["jets"])
    bjets_resolved = good_jets & (getattr(jets, parameters["btagging_algorithm"]) > parameters["btagging_WP"])
    good_fatjets = jet_selection(fatjets, muons, good_muons, parameters["fatjets"]) & jet_selection(fatjets, electrons, good_electrons, parameters["fatjets"])
#    good_fatjets = jet_selection(fatjets, muons, (veto_muons | good_muons), parameters["fatjets"]) & jet_selection(fatjets, electrons, (veto_electrons | good_electrons), parameters["fatjets"]) #FIXME remove vet_leptons

    good_jets_nohiggs = good_jets & ha.mask_deltar_first(jets, good_jets, fatjets, good_fatjets, 1.2, indices['leading'])
    bjets = good_jets_nohiggs & (getattr(jets, parameters["btagging_algorithm"]) > parameters["btagging_WP"])
    nonbjets = good_jets_nohiggs & (getattr(jets, parameters["btagging_algorithm"]) < parameters["btagging_WP"])

    # apply basic event selection -> individual categories cut later
    nmuons         = ha.sum_in_offsets(muons, good_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8) 
    nelectrons     = ha.sum_in_offsets(electrons, good_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8)
    nleps          = NUMPY_LIB.add(nmuons, nelectrons)
    lepton_veto    = NUMPY_LIB.add(ha.sum_in_offsets(muons, veto_muons, mask_events, muons.masks["all"], NUMPY_LIB.int8), ha.sum_in_offsets(electrons, veto_electrons, mask_events, electrons.masks["all"], NUMPY_LIB.int8))
    njets          = ha.sum_in_offsets(jets, nonbjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    ngoodjets      = ha.sum_in_offsets(jets, good_jets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags          = ha.sum_in_offsets(jets, bjets, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    btags_resolved = ha.sum_in_offsets(jets, bjets_resolved, mask_events, jets.masks["all"], NUMPY_LIB.int8)
    nfatjets       = ha.sum_in_offsets(fatjets, good_fatjets, mask_events, fatjets.masks['all'], NUMPY_LIB.int8)

    # trigger logic
    trigger_el = (nleps==1) & (nelectrons==1)
    trigger_mu = (nleps==1) & (nmuons==1)
    if args.year.startswith('2016'):
        trigger_el &= scalars["HLT_Ele27_WPTight_Gsf"]
        trigger_mu &= (scalars["HLT_IsoMu24"] | scalars["HLT_IsoTkMu24"])
    elif args.year.startswith('2017'):
        #trigger = (scalars["HLT_Ele35_WPTight_Gsf"] | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"] | scalars["HLT_IsoMu27"] | scalars["HLT_IsoMu24_eta2p1"]) #FIXME for different runs
        if sample.endswith(('2017B','2017C')):
            trigger_tmp = scalars["HLT_Ele32_WPTight_Gsf_L1DoubleEG"] & any([scalars[f'L1_SingleEG{n}er2p5'] for n in (10,15,26,34,36,38,40,42,45,8)])
        else:
            trigger_tmp = scalars["HLT_Ele32_WPTight_Gsf"]
        trigger_el &= (trigger_tmp | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"])
        trigger_mu &= scalars["HLT_IsoMu27"]
    elif args.year.startswith('2018'):
        trigger = (scalars["HLT_Ele32_WPTight_Gsf"] | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"] | scalars["HLT_IsoMu24"] )
        trigger_el &= (scalars["HLT_Ele32_WPTight_Gsf"] | scalars["HLT_Ele28_eta2p1_WPTight_Gsf_HT150"])
        trigger_mu &= scalars["HLT_IsoMu24"]
    if "SingleMuon" in sample: trigger_el = NUMPY_LIB.zeros(nEvents, dtype=NUMPY_LIB.bool)
    if "SingleElectron" in sample: trigger_mu = NUMPY_LIB.zeros(nEvents, dtype=NUMPY_LIB.bool)
    mask_events = mask_events & (trigger_el | trigger_mu)

    # for reference, this is the selection for the resolved analysis
    mask_events_res = mask_events & (nleps == 1) & (lepton_veto == 0) & (ngoodjets >= 4) & (btags_resolved > 2) & (scalars[metstruct+"_pt"] > 20)
    # apply basic event selection
    #mask_events_higgs = mask_events & (nleps == 1) & (scalars[metstruct+"_pt"] > 20) & (nhiggs > 0) & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4) & (btags >=2) ) & (lepton_veto == 0)
    mask_events_boost = mask_events & (nleps == 1) & (lepton_veto == 0) & (scalars[metstruct+"_pt"] > parameters['met']) & (nfatjets > 0) & (btags >= parameters['btags']) # & (btags_resolved < 3)# & (njets > 1)  # & NUMPY_LIB.invert( (njets >= 4)  )

############# calculate basic variables
    mask_events = mask_events_res | mask_events_boost
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
      'basic'    : mask_events_boost
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
    for m in ['2J2WdeltaR']:#, '2J2WdeltaRTau21']:#, '2J2WdeltaRTau21DDT']:
        mask_events[f'{m}_Pass'] = mask_events[m] & (leading_fatjet_Hbb>parameters['bbtagging_WP'])
        mask_events[f'{m}_Fail'] = mask_events[m] & (leading_fatjet_Hbb<=parameters['bbtagging_WP'])

    #mask_events['overlap'] = mask_events['2J2WdeltaR'] & mask_events['resolved']
    #mask_events['overlap'] = mask_events['2J2WdeltaR_Pass'] & mask_events['resolved']

    jetsForEff = {
            'total_flavb'  : good_jets_nohiggs & (jets.hadronFlavour==5),
            'btags_flavb'  : bjets & (jets.hadronFlavour==5),
            'total_flavlc' : good_jets_nohiggs & (jets.hadronFlavour<5),
            'btags_flavlc' : bjets & (jets.hadronFlavour<5),
            'total_flavl'  : good_jets_nohiggs & (jets.hadronFlavour==0),
            'btags_flavl'  : bjets & (jets.hadronFlavour==0)
            }

    binsPt_comb_central = [20,1000]
    binsPt_comb_updown  = [20,30,50,70,100,140,200,300,600,1000]
    binsPt_incl         = [20,1000]

    binsForEff = {
            'flavb_central' : binsPt_comb_central,
            'flavb_updown'  : binsPt_comb_updown,
            'incl'          : binsPt_incl
            }

    for i in range(len(mask_events['2J2WdeltaR_Pass'])):
        if not mask_events['2J2WdeltaR_Pass'][i]:
            start = jets.offsets[i]
            end   = jets.offsets[i+1]
            for j in jetsForEff.values():
                j[start:end] = 0

    ret['total_flavb_central']  = get_histogram(jets.pt[jetsForEff['total_flavb']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavb'])), binsForEff['flavb_central'])
    ret['total_flavb_updown']   = get_histogram(jets.pt[jetsForEff['total_flavb']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavb'])), binsForEff['flavb_updown'])
    ret['btags_flavb_central']  = get_histogram(jets.pt[jetsForEff['btags_flavb']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavb'])), binsForEff['flavb_central'])
    ret['btags_flavb_updown']   = get_histogram(jets.pt[jetsForEff['btags_flavb']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavb'])), binsForEff['flavb_updown'])
    ret['total_flavlc_central'] = get_histogram(jets.pt[jetsForEff['total_flavlc']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavlc'])), binsForEff['incl'])
    ret['total_flavlc_updown']  = get_histogram(jets.pt[jetsForEff['total_flavlc']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavlc'])), binsForEff['incl'])
    ret['btags_flavlc_central'] = get_histogram(jets.pt[jetsForEff['btags_flavlc']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavlc'])), binsForEff['incl'])
    ret['btags_flavlc_updown']  = get_histogram(jets.pt[jetsForEff['btags_flavlc']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavlc'])), binsForEff['incl'])
    ret['total_flavl_central']  = get_histogram(jets.pt[jetsForEff['total_flavl']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavl'])), binsForEff['incl'])
    ret['total_flavl_updown']   = get_histogram(jets.pt[jetsForEff['total_flavl']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['total_flavl'])), binsForEff['incl'])
    ret['btags_flavl_central']  = get_histogram(jets.pt[jetsForEff['btags_flavl']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavl'])), binsForEff['incl'])
    ret['btags_flavl_updown']   = get_histogram(jets.pt[jetsForEff['btags_flavl']], NUMPY_LIB.ones(NUMPY_LIB.sum(jetsForEff['btags_flavl'])), binsForEff['incl'])

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
    parser.add_argument('--corrections', action='store_true', help='Flag to include corrections')
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
    from coffea.btag_tools import BTagScaleFactor

    # load definitions
    from definitions_analysis import parameters, eraDependentParameters, samples_info
    parameters.update(eraDependentParameters[args.year])
    if args.parameters is not None:
      if len(args.parameters)%2 is not 0:
        raise Exception('incomplete parameters specified, quitting.')
      for p,v in zip(args.parameters[::2], args.parameters[1::2]):
        try: parameters[p] = type(parameters[p])(v) #convert the string v to the type of the parameter already in the dictionary
        except: print(f'invalid parameter specified: {p} {v}')

    if "Single" in args.sample:
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
      if args.sample.endswith(('2017B','2017C')):
          arrays_event += ["HLT_Ele32_WPTight_Gsf_L1DoubleEG"]
          arrays_event += [f'L1_SingleEG{n}er2p5' for n in (10,15,26,34,36,38,40,42,45,8)]
      else:
        arrays_event += ["HLT_Ele32_WPTight_Gsf"] #FIXME
    elif args.year.startswith('2018'): arrays_event += [ "HLT_Ele32_WPTight_Gsf", 'HLT_Ele28_eta2p1_WPTight_Gsf_HT150', "HLT_IsoMu24" ]

    if args.sample.startswith("TT"):
        arrays_event.append("genTtbarId")

    if is_mc:
        arrays_event += ["PV_npvsGood", "Pileup_nTrueInt", "genWeight", "nGenPart", 'LHEPdfWeight', 'PSWeight']
        arrays_objects += [ "Jet_hadronFlavour", #"selectedPatJetsAK4PFPuppi_hadronFlavor",
                           "GenPart_eta","GenPart_genPartIdxMother","GenPart_mass","GenPart_pdgId","GenPart_phi","GenPart_pt","GenPart_status","GenPart_statusFlags"
                         ]

    if args.corrections: 
      arrays_objects += [f'FatJet_{var}_{corr}' for var in ['msoftdrop','pt','mass'] for corr in ['raw','nom']]
      arrays_objects += [f'Jet_{var}_nom' for var in ['pt','mass']]
      arrays_event   += [f'{metstruct}_{var}_nom' for var in ['pt','phi']]
      if is_mc:#args.sample.startswith('ttH'):
        arrays_event   += [f'{metstruct}_{var}_jer{ud}' for var in ['pt','phi'] for ud in ['','Up','Down']]
        arrays_objects += [f'Jet_btagSF_deepjet_M{var}' for var in ['','_up','_down']]
        arrays_objects += [f'FatJet_msoftdrop_{unc}{ud}' for unc in ['jmr','jms'] for ud in ['Up','Down']]
        arrays_event   += [f'puWeight{var}' for var in ['','Up','Down']]
        jesSources = ['Total','Absolute',f'Absolute_{args.year}','FlavorQCD','BBEC1',f'BBEC1_{args.year}','EC2',f'EC2_{args.year}','HF','RelativeBal',f'RelativeSample_{args.year}']
        jesSourcesSplit = ['Total','AbsoluteMPFBias','AbsoluteScale','AbsoluteStat','FlavorQCD','Fragmentation','PileUpDataMC','PileUpPtBB','PileUpPtEC1','PileUpPtEC2','PileUpPtHF','PileUpPtRef','RelativeFSR','RelativeJEREC1','RelativeJEREC2','RelativeJERHF','RelativePtBB','RelativePtEC1','RelativePtEC2','RelativePtHF','RelativeBal','RelativeSample','RelativeStatEC','RelativeStatFSR','RelativeStatHF','SinglePionECAL','SinglePionHCAL','TimePtEta']
 
        arrays_event   += [f'{metstruct}_{var}_{unc}{ud}' for var in ['pt','phi'] for unc in [f'jes{s}' for s in (jesSources if args.sample.startswith('ttH') else jesSourcesSplit)] for ud in ['Up','Down']]
        arrays_objects += [f'FatJet_{var}_{unc}{ud}' for var in ['pt','mass','msoftdrop'] for unc in ['jer']+[f'jes{s}' for s in (jesSources if args.sample.startswith('ttH') else jesSourcesSplit)] for ud in ['Up','Down']]
        arrays_objects += [f'Jet_{var}_{unc}{ud}' for var in ['pt','mass'] for unc in ['jer']+[f'jes{s}' for s in (jesSources if args.sample.startswith('ttH') else jesSourcesSplit)] for ud in ['Up','Down']]

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

    WPs_DAK8 = [0.8695, 0.9795]#0.5845, 
    WPs_DDB  = [0.86]#, 0.89, 0.91]#, 0.92]0.7, 
    bbtags   = {'btagDDBvL': WPs_DDB} #'deepTagMD_bbvsLight': WPs_DAK8, 'btagDDBvL_noMD': WPs_DDB, 'deepTag_H': WPs_DAK8, 'btagDDBvL': WPs_DDB, 'btagDDBvL_noMD': WPs_DDB}
    pars     = {f'met{met}_{bbAlg}0{str(bbWP).split(".")[-1]}' : (met,bbAlg,bbWP) for met in [20] for bbAlg,bbWPlist in bbtags.items() for bbWP in bbWPlist}
    #pars['met30_btagDDBvL086'] = (30, 'btagDDBvL', 0.86)
    for p in pars.copy():
        #pars[f'{p}_1btag'] = pars[p] + (1,)
        pars[p]            = pars[p] + (0,)

    results  = {p : {} for p in pars}

    if args.corrections: 
      uncertainties = {
            'nominal'      : [[],['FatJet','msoftdrop','msoftdrop_nom']],
            #'msd_raw'      : [[],['FatJet','msoftdrop','msoftdrop_raw']],
            #'msd_nanoAOD'  : ['FatJet','msoftdrop','msoftdrop']
            }
      if is_mc:#args.sample.startswith('ttH'):
        signalUncertainties = {
            #'jerUp'        : [[],['FatJet','pt','pt_jerUp','FatJet','mass','mass_jerUp']],
            #'jerDown'      : [[],['FatJet','pt','pt_jerDown','FatJet','mass','mass_jerDown']],
            #'jesTotalUp'   : [[],['FatJet','pt','pt_jesTotalUp','FatJet','mass','mass_jesTotalUp']],
            #'jesTotalDown' : [[],['FatJet','pt','pt_jesTotalDown','FatJet','mass','mass_jesTotalDown']],
            #'jmrUp'        : [[],['FatJet','msoftdrop','msoftdrop_jmrUp']],
            #'jmrDown'      : [[],['FatJet','msoftdrop','msoftdrop_jmrDown']],
            #'jmsUp'        : [[],['FatJet','msoftdrop','msoftdrop_jmsUp']],
            #'jmsDown'      : [[],['FatJet','msoftdrop','msoftdrop_jmsDown']],
            #'puWeightUp'   : [['puWeight','puWeightUp'],[]],
            #'puWeightDown' : [['puWeight','puWeightDown'],[]],
            }
        for variation in ['Up','Down']:
            signalUncertainties[f'psWeight_ISR{variation}']       = [[],[]]
            signalUncertainties[f'psWeight_FSR{variation}']       = [[],[]]
            signalUncertainties[f'pdfWeight{variation}']          = [[],[]]
            signalUncertainties[f'puWeight{variation}']           = [['puWeight',f'puWeight{variation}'],[]]
            #signalUncertainties[f'AK8jer{variation}']             = [[],['FatJet','pt','pt_jer{variation}','FatJet','mass','mass_jer{variation}','FatJet','msoftdrop','msoftdrop_jer{variation}']]
            #signalUncertainties[f'AK4jer{variation}']             = [[],['Jet','pt','pt_jer{variation}','Jet','mass','mass_jer{variation}']]
            signalUncertainties[f'jer{variation}']                = [
                    [f'{metstruct}_pt',f'{metstruct}_pt_jer{variation}',f'{metstruct}_phi',f'{metstruct}_phi_jer{variation}'],
                    ['FatJet','pt',f'pt_jer{variation}','FatJet','mass',f'mass_jer{variation}','FatJet','msoftdrop',f'msoftdrop_jer{variation}','Jet','pt',f'pt_jer{variation}','Jet','mass',f'mass_jer{variation}'] ]
            signalUncertainties[f'AK8DDBvLM1{variation}']         = [[],['FatJet','bbtagSF_DDBvL_M1',f'bbtagSF_DDBvL_M1_{variation.lower()}']]
            signalUncertainties[f'AK4deepjetM{variation}']        = [[],['Jet','btagSF_deepjet_M',f'btagSF_deepjet_M_{variation.lower()}']]
            for source in ['jmr','jms']:
                signalUncertainties[f'{source}{variation}']       = [[],['FatJet','msoftdrop',f'msoftdrop_{source}{variation}']]
            for source in jesSources:
                #signalUncertainties[f'AK8jes{source}{variation}'] = [[],['FatJet','pt','pt_jes{source}{variation}','FatJet','mass','mass_jes{source}{variation}','FatJet','msoftdrop','msoftdrop_jes{source}{variation}']]
                #signalUncertainties[f'AK4jes{source}{variation}'] = [[],['Jet','pt','pt_jes{source}{variation}','Jet','mass','mass_jes{source}{variation}']]
                signalUncertainties[f'jes{source}{variation}'] = [
                        [f'{metstruct}_pt',f'{metstruct}_pt_jes{source}{variation}',f'{metstruct}_phi',f'{metstruct}_phi_jes{source}{variation}'],
                        ['FatJet','pt',f'pt_jes{source}{variation}','FatJet','mass',f'mass_jes{source}{variation}','FatJet','msoftdrop',f'msoftdrop_jes{source}{variation}','Jet','pt',f'pt_jes{source}{variation}','Jet','mass',f'mass_jes{source}{variation}']]
        uncertainties.update(signalUncertainties)
      for u in uncertainties.values():
        if not f'{metstruct}_pt' in u[0]:
            metBranch = 'jer' if is_mc else 'nom'
            u[0] += [f'{metstruct}_pt',f'{metstruct}_pt_{metBranch}', f'{metstruct}_phi',f'{metstruct}_phi_{metBranch}']
        if not 'Jet' in u[1]:
            u[1] += ['Jet','pt','pt_nom', 'Jet','mass','mass_nom']
        if 'msoftdrop' not in u[1]:
            u[1] += ['FatJet','msoftdrop','msoftdrop_nom']
        if not ('FatJet','pt') in zip(u[1][::3],u[1][1::3]):
            u[1] += ['FatJet','pt','pt_nom', 'FatJet','mass','mass_nom']
      extraCorrections = {
          'no_PUPPI'         : ['PUPPI'],
          }
      arrays_objects += [f'FatJet_msoftdrop_corr_{e}' for e in sum(extraCorrections.values(),[])]
      #results = {u : Results() for u in uncertainties} 
      print(uncertainties)
      #with open('unc_dump.json','w') as f: json.dump(uncertainties, f, indent=4)
      #sys.exit()
      for p in results:
        results[p] = {u : Results() for u in uncertainties}
      ######### this block was for testing which corrections we should remove
#      if is_mc:
#        extraCorrections = {
#            'no_PUPPI_JMS_JMR' : ['PUPPI','JMS','JMR'],
#            'no_JMS_JMR'       : ['JMS','JMR'],
#            'no_PUPPI'         : ['PUPPI'],
#            #'no_JMR'           : ['JMR'],
#            #'no_JMS'           : ['JMS'],
#            }
#      else:
#        extraCorrections = {
#            'no_PUPPI'         : ['PUPPI'],
#            }
#
#      extraCorrections[''] = None
#      combinations = [('msd_nom',e) for e in extraCorrections]
#      combinations += [('msd_raw','')]
#      results = {(u[0] if u[1]=='' else u[1]) : Results() for u in combinations}
    else:
      uncertainties = {'' : None}
      results       = {'' : Results()}

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
            #parameters["btag_SF_target"] = BTagScaleFactor(parameters["btag_SF_{}".format(parameters["btagging_algorithm"])], BTagScaleFactor.RESHAPE, 'iterativefit,iterativefit,iterativefit', keep_df=True)
            ### this computes the lepton weights
            ext = extractor()
            print(parameters["corrections"])
            for corr in parameters["corrections"]:
                ext.add_weight_sets([corr])
            ext.finalize()
            evaluator = ext.make_evaluator()

        if ibatch == 0:
            print(dataset.printout())

        for p in pars:
          parameters['met'], parameters['bbtagging_algorithm'], parameters['bbtagging_WP'], parameters['btags'] = pars[p] #
          for un,u in uncertainties.items():
          #### this is where the magic happens: run the main analysis
            if un!='nominal': continue
            results[p][un] += dataset.analyze(analyze_data, NUMPY_LIB=NUMPY_LIB, parameters=parameters, is_mc = is_mc, lumimask=lumimask, cat=args.categories, sample=args.sample, samples_info=samples_info, boosted=args.boosted, uncertainty=u, uncertaintyName=un, parametersName=p, extraCorrection=extraCorrections['no_PUPPI'])

    #print(results)

    #Save the results
#    outdir = args.outdir
#    if args.version!='':
#      outdir = os.path.join(outdir,args.version)
#    auxdir = f"met{parameters['met']}_{parameters['bbtagging_algorithm']}0{str(parameters['bbtagging_WP']).split('.')[-1]}"
#    outdir = os.path.join(outdir,auxdir)
#    if not os.path.exists(outdir):
#      os.makedirs(outdir)
#
#    for r in results:
#      results[r].save_json(os.path.join(outdir,f"out_{args.sample}_{r}{args.outtag}.json"))
    for pn,res in results.items():
      for rn,r in res.items():
        outdir = args.outdir
        if args.version!='':
          outdir = os.path.join(outdir,args.version)
        outdir = os.path.join(outdir,pn,rn,'btagEfficiencyMaps')
        if not os.path.exists(outdir):
          os.makedirs(outdir)

        r.save_json(os.path.join(outdir,f"out_btagEfficiencyMaps_{args.sample}_{rn}{args.outtag}.json"))

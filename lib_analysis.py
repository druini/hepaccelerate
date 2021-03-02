import os, glob
import argparse
import json
import numpy as np

import uproot
import hepaccelerate

from hepaccelerate.utils import Results, NanoAODDataset, Histogram, choose_backend

from uproot_methods import TLorentzVectorArray

NUMPY_LIB = None
ha = None

############################################## OBJECT SELECTION ################################################

### Primary vertex selection
def vertex_selection(scalars, mask_events):

    PV_isfake = (scalars["PV_score"] == 0) & (scalars["PV_chi2"] == 0)
    PV_rho = NUMPY_LIB.sqrt(scalars["PV_x"]**2 + scalars["PV_y"]**2)
    mask_events = mask_events & (~PV_isfake) & (scalars["PV_ndof"] > 4) & (scalars["PV_z"]<24) & (PV_rho < 2)

    return mask_events


### Lepton selection
def lepton_selection(leps, cuts, year):

    passes_eta = (NUMPY_LIB.abs(leps.eta) < cuts["eta"])
    passes_subleading_pt = (leps.pt > cuts["subleading_pt"])
    passes_leading_pt = (leps.pt > cuts["leading_pt"][year])

    if cuts["type"] == "el":
        sca = NUMPY_LIB.abs(leps.deltaEtaSC + leps.eta)
        passes_id = (leps.cutBased >= 4)
        passes_SC = NUMPY_LIB.invert((sca >= 1.4442) & (sca <= 1.5660))
        # cuts taken from: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_92X_and_later
        passes_impact = ((leps.dz < 0.10) & (sca <= 1.479)) | ((leps.dz < 0.20) & (sca > 1.479)) | ((leps.dxy < 0.05) & (sca <= 1.479)) | ((leps.dxy < 0.1) & (sca > 1.479))

        #select electrons
        good_leps = passes_eta & passes_leading_pt & passes_id & passes_SC & passes_impact
        veto_leps = passes_eta & passes_subleading_pt & NUMPY_LIB.invert(good_leps) & passes_id & passes_SC & passes_impact

    elif cuts["type"] == "mu":
        passes_leading_iso = (leps.pfRelIso04_all < cuts["leading_iso"])
        passes_subleading_iso = (leps.pfRelIso04_all < cuts["subleading_iso"])
        passes_id = (leps.tightId == 1)

        #select muons
        good_leps = passes_eta & passes_leading_pt & passes_leading_iso & passes_id
        veto_leps = passes_eta & passes_subleading_pt & passes_subleading_iso & passes_id & NUMPY_LIB.invert(good_leps)

    return good_leps, veto_leps

### Jet selection
def jet_selection(jets, leps, mask_leps, cuts):

    jets_pass_dr = ha.mask_deltar_first(jets, jets.masks["all"], leps, mask_leps, cuts["dr"])
    jets.masks["pass_dr"] = jets_pass_dr
    good_jets = (jets.pt > cuts["pt"]) & (NUMPY_LIB.abs(jets.eta) < cuts["eta"]) & (jets.jetId >= cuts["jetId"]) & jets_pass_dr
    if cuts["type"] == "jet":
      good_jets &= ((jets.pt<50) & (jets.puId>=cuts["puId"]) ) | (jets.pt>=50) 

    return good_jets


###################################################### WEIGHT / SF CALCULATION ##########################################################

### PileUp weight
def compute_pu_weights(pu_corrections_target, weights, mc_nvtx, reco_nvtx):
    pu_edges, (values_nom, values_up, values_down) = pu_corrections_target

    src_pu_hist = get_histogram(mc_nvtx, weights, pu_edges)
    norm = sum(src_pu_hist.contents)
    src_pu_hist.contents = src_pu_hist.contents/norm
    src_pu_hist.contents_w2 = src_pu_hist.contents_w2/norm

#    fi = uproot.open('/afs/cern.ch/user/a/algomez/public/forDaniele/mcPileup2017.root')
#    h = fi['pu_mc']
#    mc_edges = np.array(h.edges)
#    mc_values = np.array(h.values)
#    mc_values /= np.sum(mc_values)
#    mc_values = np.append(mc_values, 1)

    ratio = values_nom / src_pu_hist.contents
#    ratio = values_nom / mc_values
    remove_inf_nan(ratio)
    pu_weights = NUMPY_LIB.zeros_like(weights)
    ha.get_bin_contents(reco_nvtx, NUMPY_LIB.array(pu_edges), NUMPY_LIB.array(ratio), pu_weights)
    #fix_large_weights(pu_weights)

    return pu_weights


def load_puhist_target(filename):
    fi = uproot.open(filename)

    h = fi["pileup"]
    edges = np.array(h.edges)
    values_nominal = np.array(h.values)
    values_nominal = values_nominal / np.sum(values_nominal)

    h = fi["pileup_plus"]
    values_up = np.array(h.values)
    values_up = values_up / np.sum(values_up)

    h = fi["pileup_minus"]
    values_down = np.array(h.values)
    values_down = values_down / np.sum(values_down)
    return edges, (values_nominal, values_up, values_down)


# lepton scale factors
def compute_lepton_weights(leps, lepton_pt, lepton_eta, mask_rows, mask_content, evaluator, SF_list, year=None):

    weights = NUMPY_LIB.ones(len(lepton_pt))

    for SF in SF_list:
        if SF.startswith('mu'):
            if year=='2016':
                if 'trigger' in SF:
                    x = lepton_pt
                    y = NUMPY_LIB.abs(lepton_eta)
                else:
                    x = lepton_eta
                    y = lepton_pt
            else:
                x = lepton_pt
                y = NUMPY_LIB.abs(lepton_eta)
        elif SF.startswith('el'):
            if 'trigger' in SF:
                x = lepton_pt
                y = lepton_eta
            else:
                x = lepton_eta
                y = lepton_pt
        else:
            raise Exception(f'unknown SF name {SF}')
        weights *= evaluator[SF](x, y)
    
    per_event_weights = ha.multiply_in_offsets(leps, weights, mask_rows, mask_content)
    return per_event_weights

def my_SF_extractor(weightdesc, variation):
    (local_name, name, thefile) = tuple(weightdesc.strip().split(" "))
    with uproot.open(thefile) as f:
        weights = f[name].values + (1. if variation=='Up' else -1.)*np.sqrt(f[name].variances)
        edges   = f[name].edges
    return (weights, edges)

# btagging scale factor 
def compute_btag_weights(jets, mask_rows, mask_content, systematic, parameters):
    btagalgorithm = parameters['btagging_algorithm']
    btagWP        = parameters['btagging_WP']

    tagged    = mask_content & (getattr(jets, btagalgorithm)>btagWP)
    nontagged = mask_content & (getattr(jets, btagalgorithm)<btagWP)
    tag_weight    = NUMPY_LIB.ones_like(jets.pt)
    nontag_weight = NUMPY_LIB.ones_like(jets.pt)
    bjets = mask_content & (jets.hadronFlavour==5)
    cjets = mask_content & (jets.hadronFlavour==4)
    ljets = mask_content & (jets.hadronFlavour==0)

    from coffea.btag_tools import BTagScaleFactor
    hfsf = BTagScaleFactor(parameters[f'btag_SF_{btagalgorithm}_YearCorrelation'], BTagScaleFactor.MEDIUM)
    lfsf = BTagScaleFactor(parameters[f'btag_SF_{btagalgorithm}'], BTagScaleFactor.MEDIUM)

    tag_weight[bjets] = hfsf.eval(systematic, 5, abs(jets.eta[bjets]), jets.pt[bjets], ignore_missing=True)
    #tag_weight[cjets] = sf.eval(systematic, 4, abs(jets.eta[cjets]), jets.pt[cjets], ignore_missing=True)
    if systematic.endswith('uncorrelated'):
        lfsystematic = systematic.split('_')[0]
    elif systematic.endswith('correlated'):
        lfsystematic = 'central'
    else:
        lfsystematic = systematic
    tag_weight[ljets] = lfsf.eval(lfsystematic, 0, abs(jets.eta[ljets]), jets.pt[ljets], ignore_missing=True)

    for flav in [bjets,ljets]:
        nontag_weight[flav] = (1 - tag_weight[flav]*jets.btag_MCeff[flav]) / (1 - jets.btag_MCeff[flav])

    evWeights_tagged    = ha.multiply_in_offsets(jets, tag_weight, mask_rows, tagged)
    evWeights_nontagged = ha.multiply_in_offsets(jets, nontag_weight, mask_rows, nontagged)

    per_event_weights = evWeights_tagged * evWeights_nontagged
    return per_event_weights

############################################# HIGH LEVEL VARIABLES (DNN evaluation, ...) ############################################

# calculate simple object variables
def calculate_variable_features(z, mask_events, indices, var):

    name, coll, mask_content, inds, feats = z
    idx = indices[inds]

    for f in feats:
        var[inds+"_"+name+"_"+f] = ha.get_in_offsets(getattr(coll, f), getattr(coll, "offsets"), idx, mask_events, mask_content)
    
####################################################### Simple helpers  #############################################################

def get_histogram(data, weights, bins):
    return Histogram(*ha.histogram_from_vector(data, weights, bins))

def remove_inf_nan(arr):
    arr[np.isinf(arr)] = 0
    arr[np.isnan(arr)] = 0
    arr[arr < 0] = 0

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

#import keras.backend as K
#import keras.losses
#import keras.utils.generic_utils
#
#def mse0(y_true,y_pred):
#    return K.mean( K.square(y_true[:,0] - y_pred[:,0]) )
#
#def mae0(y_true,y_pred):
#    return K.mean( K.abs(y_true[:,0] - y_pred[:,0]) )
#
#def r2_score0(y_true,y_pred):
#    return 1. - K.sum( K.square(y_true[:,0] - y_pred[:,0]) ) / K.sum( K.square(y_true[:,0] - K.mean(y_true[:,0]) ) )

def select_lepton_p4(objs1, mask1, objs2, mask2, indices, mask_rows):
  selected_obj1 = {}
  selected_obj2 = {}
  feats = ['pt','eta','phi','mass']
  for feat in feats:
    selected_obj1[feat] = ha.get_in_offsets(getattr(objs1,feat), objs1.offsets, indices, mask_rows, mask1)
    selected_obj2[feat] = ha.get_in_offsets(getattr(objs2,feat), objs2.offsets, indices, mask_rows, mask2)
  select_1_or_2 = (selected_obj1['pt'] > selected_obj2['pt'])
  selected_feats = {}
  for feat in feats:
    selected_feats[feat] = NUMPY_LIB.where(select_1_or_2, selected_obj1[feat], selected_obj2[feat])
  selected_p4 = TLorentzVectorArray.from_ptetaphim(selected_feats['pt'], selected_feats['eta'], selected_feats['phi'], selected_feats['mass'])
  return selected_p4

def hadronic_W(jets, jets_mask, lepWp4, mask_rows):
  from itertools import combinations
  init = -999.*np.zeros(len(jets.offsets) - 1, dtype=np.float32) 
  hadW = TLorentzVectorArray.from_ptetaphim(init.copy(), init.copy(), init.copy(), init.copy())
  for iev in range(jets.offsets.shape[0]-1):
    if not mask_rows[iev]: continue
    start = jets.offsets[iev]
    end = jets.offsets[iev + 1]
    smallestDiffW = 9999.
    for jpair in combinations(jets.p4[start:end][jets_mask[start:end]], 2):
      tmphadW = jpair[0] + jpair[1]
      tmpDiff = abs(lepWp4[iev].mass - tmphadW.mass)
      if tmpDiff<smallestDiffW:
        smallestDiffW = tmpDiff
        for feat in ['pt','eta','phi','mass']:
          getattr(hadW, feat)[iev] = getattr(tmphadW, feat)
  return hadW

def loadMCeff(jets, btag_MCeff_json, sysType):
    with open(btag_MCeff_json) as f:
        btag_MCeff = json.load(f)
    for e in btag_MCeff:
        btag_MCeff[e] = Histogram( *btag_MCeff[e].values() )
        #btag_MCeff[e].contents[ NUMPY_LIB.isnan(btag_MCeff[e].contents) ] = 1
        remove_inf_nan(btag_MCeff[e].contents)

    jets_btag_MCeff = NUMPY_LIB.ones_like(jets.pt)
    bjets = jets.hadronFlavour==5
    #cjets = jets.hadronFlavour==4
    ljets = jets.hadronFlavour==0
    ptbins_flavb = btag_MCeff[f'eff_flavb_{sysType}'].edges
    def idx(bins, pt):
        ### based on https://github.com/CoffeaTeam/coffea/blob/master/coffea/lookup_tools/dense_mapped_lookup.py#L34
        if len(bins)==2:
            return NUMPY_LIB.zeros_like(pt, dtype=NUMPY_LIB.uint)
        return NUMPY_LIB.clip(NUMPY_LIB.searchsorted(bins, pt, side='right')-1,0,len(bins)-2)
    eff_flavb = btag_MCeff[f'eff_flavb_{sysType}']
    eff_flavl = btag_MCeff[f'eff_flavl_{sysType}']
    jets_btag_MCeff[bjets] = eff_flavb.contents[ idx(eff_flavb.edges, jets.pt[bjets]) ]
    jets_btag_MCeff[ljets] = eff_flavl.contents[ idx(eff_flavl.edges, jets.pt[ljets]) ]
    return jets_btag_MCeff

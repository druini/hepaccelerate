

######################################## Selection criteria / Inputs for corrections ############################################

parameters = {
  "muons": {
    "type": "mu",
    "leading_pt": 29,
    "subleading_pt": 15,
    "eta": 2.4,
    "leading_iso": 0.15,
    "subleading_iso": 0.25,
  },

  "electrons" : {
    "type": "el",
    "leading_pt": 30,
    "subleading_pt": 15,
    "eta": 2.4
  },
  "jets": {
    "type": "jet",
    "dr": 0.4,
    "pt": 30,
    "eta": 2.4,
    "jetId": 2,
    "puId": 4
  },
  "fatjets": {
    "type": "fatjet",
    "dr": 0.8,
    "pt": 250,
    "eta": 2.4,
    "jetId": 2,
    "tau32cut": 0.4,
    "tau21cut": 0.35,
  },
  "W": {
    'min_mass': 60,
    'max_mass': 105
  }
}

eraDependentParameters = {
    "2016" : {  ###### NOT UPDATED!!!!
        "lumi":  36773.0,
        "lumimask": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt",
        "pu_corrections_file" : "data/pileup_profile_Summer16.root",
        "corrections" : [
            "el_triggerSF Ele27_WPTight_Gsf data/TriggerSF_Run2016All_v1.root",
            "el_recoSF EGamma_SF2D data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root",
            "el_idSF EGamma_SF2D data/2016LegacyReReco_ElectronTight_Fall17V2.root",
            "mu_triggerSF IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio data/EfficienciesAndSF_RunBtoF.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt data/MuonID_2016_RunBCDEF_SF_ISO.root",
            "mu_idSF NUM_TightID_DEN_genTracks_eta_pt data/MuonID_2016_RunBCDEF_SF_ID.root",
            "BTagSF * data/DeepCSV_Moriond17_B_H.csv"
        ],
        "btagging algorithm" : "btagDeepFlavB",#"btagDeepB",
        "btagging WP" : 0.3093,
        "bbtagging algorithm" :  "deepTagMD_bbvsLight",#"deepTagMD_HbbvsQCD", #"btagDDBvL",
        "bbtagging WP" : 0.8945,
    },
    "2017" : {
        "lumi":  41529.0,
        "lumimask": "./data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt",
        "pu_corrections_file" : "./data/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root",
        "corrections" : [
            "el_triggerSF SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF ./data/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5_0_histo.root",
            "el_recoSF EGamma_SF2D ./data/egammaEffi_EGM2D_runBCDEF_passingRECO_histo.root",
            "el_idSF EGamma_SF2D ./data/egammaEffi_EGM2D_runBCDEF_passingTight94X_histo.root",
            "mu_triggerSF IsoMu27_PtEtaBins/pt_abseta_ratio ./data/EfficienciesAndSF_RunBtoF_Nov17Nov2017_histo.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta ./data/RunBCDEF_SF_ISO_histo.root",
            "mu_idSF NUM_TightID_DEN_genTracks_pt_abseta ./data/RunBCDEF_SF_ID_histo.root",
            "BTagSF * ./data/deepCSV_sfs_v2_btag.csv"
        ],
        "btagging algorithm" : "btagDeepFlavB",#"btagDeepB",
        "btagging WP" : 0.3033, # 0.4941, # medium working point for btagDeepB
        "bbtagging algorithm" :  "deepTagMD_bbvsLight",#"deepTagMD_HbbvsQCD", #"btagDDBvL",
        "bbtagging WP" : 0.8695, ### https://indico.cern.ch/event/853828/contributions/3723593/attachments/1977626/3292045/lg-btv-deepak8v2-sf-20200127.pdf
    },
    "2018" : {
        "lumi":  58830.0,
        "lumimask": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt",
        "pu_corrections_file" : "./data/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root",
        "corrections" : [
            "el_triggerSF SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF ./data/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5_0.root",
            "el_recoSF EGamma_SF2D ./data/2018_ElectronTight.root",
            "el_idSF EGamma_SF2D ./data/egammaEffi_txt_EGM2D_updatedAll.root",
            "mu_triggerSF IsoMu24_PtEtaBins/pt_abseta_ratio ./data/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta ./data/MuonID_2018_RunABCD_SF_ISO.root",
            "mu_idSF NUM_TightID_DEN_TrackerMuons_pt_abseta ./data/MuonID_2018_RunABCD_SF_ID.root",
            "BTagSF * ./data/DeepCSV_Moriond17_B_H.csv"
        ],
        "btagging algorithm" : "btagDeepFlavB",#"btagDeepB",
        "btagging WP" : 0.2770, # medium working point for btagDeepB
        "bbtagging algorithm" :  "deepTagMD_bbvsLight",#"deepTagMD_HbbvsQCD", #"btagDDBvL",
        "bbtagging WP" : 0.8365,
    },

}


#################################### Samples info ##############################################################################

samples_info = {
    "Run_dummy": {
            "XS": 1,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 1.,
                '2018' : 1.,
                },
            },
    "mc_dummy": {
            "XS": 1,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 1.,
                '2018' : 1.,
                },
            },
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8": {
            "process": "ttHTobb",
            "XS": 0.2953,
            #"ngen_weight": 2410620.0644499995 #reduced file list
            "ngen_weight": {
                '2016' : 5253482.85,
                '2017' : 4216319.315883999,
                '2018' : 5046714.41,
                },
            },
    "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 365.4574,
            "ngen_weight": {
                '2016' : 32276798335.029,
                '2017' : 720253370.0403845,
                '2018' : 1.,
                },
            },
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8": {
            "process": "ttHToNonbb",
            "XS": 0.2118,
            #"ngen_weight": 913045.7391360003 #reduced file list
            "ngen_weight": {
                '2016' : 5248991.57,
                '2017' : 3095197.8117420007,
                '2018' : 3963935.78,
                },
            },
    "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 88.3419,
            "ngen_weight": {
                '2016' : 4746866164.44,
                '2017' : 283000430.5968169,
                '2018' : 1.,
                },
            },
    "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 377.9607,
            "ngen_weight": {
                '2016' : 21500086465.24,
                '2017' : 1647945788.3386502,
                '2018' : 1.,
                },
            },

    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 88.3419,
            "ngen_weight": {
                '2016' : 1,
                '2017' : 643135033.87298, #283000430.5968169
                '2018' : 4635769336.53,
                },
            },
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 377.9607,
            "ngen_weight": {
                '2016' : 1,
                '2017' : 41084368.0, #1647945788.338650
                '2018' : 41941959327.80,
                },
            },
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 365.4574,
            "ngen_weight": {
                '2016' : 1,
                '2017' : 10758424549.765732, #720253370.0403845
                '2018' : 30545238400.06,
                },
            },
    "TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8": {
      "XS": 0.6012,# 0.5297,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 383062.0686438,
                '2018' : 381537.24,
                },
      },
    "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8": {
      "XS": 0.6012,# 0.5297,
            "ngen_weight": {
                '2016' : 396340.93,
                '2017' : 1,
                '2018' : 1.,
                },
      },
    "ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8": {
      "process": "ST_s-channel",
      "XS": 3.36,#3.702224,
      #"ngen_weight": 24856809.513425056 #reduced file list
            "ngen_weight": {
                '2016' : 36768937.25,
                '2017' : 36781553.92694208,
                '2018' : 1.,
                },
      },
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8": {
      "process": "ST_s-channel",
      "XS": 3.36,#3.702224,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 1.,
                '2018' : 74634736.73,
                },
      },
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": {
      "process": "ST_tW_antitop",
      "XS": 35.85,
      #"ngen_weight": 182291193.36093727 #reduced file list
            "ngen_weight": {
                '2016' : 174109580.67,
                '2017' : 270762750.1725247,
                '2018' : 1.,
                },
      },
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8": {
      "process": "ST_tW_antitop",
      "XS": 35.85,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 1.,
                '2018' : 266470421.96,
                },
      },
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": {
      "process": "ST_tW_top",
      "XS": 35.85,
      #"ngen_weight": 241590614.9098064 #reduced file list
            "ngen_weight": {
                '2016' : 173908712.95,
                '2017' : 277241050.84022224,
                '2018' : 1.,
                },
      },
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8": {
      "process": "ST_tW_top",
      "XS": 35.85,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 1.,
                '2018' : 334874722.208,
                },
      },
    "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8": {
      "process": "ST_t-channel_antitop",
      "XS": 80.95,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 3675910.0,
                '2018' : 5125996535.38,
                },
      },
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": {
      "process": "ST_t-channel_antitop",
      "XS": 80.95,
            "ngen_weight": {
                '2016' : 17771478.65,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8": {
      "process": "ST_t-channel_top",
      "XS": 136.02,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 5863722.0,
                '2018' : 16603455266.97,
                },
      },
    "ST_t-channel_top_4f_inclusiveDecays_13TeV_PSweights-powhegV2-madspin": {
      "process": "ST_t-channel_top",
      "XS": 136.02,
            "ngen_weight": {
                '2016' : 67917907.24,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8": {
      "XS":  0.01517,#0.002879,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 4714331.0,
                '2018' : 14971606.149,
                },
      },
    "THW_ctcvcp_HIncl_M125_TuneCP5_13TeV-madgraph-pythia8": {
      "XS":  0.01517,#0.002879,
            "ngen_weight": {
                '2016' : 4989133.867,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 3.697,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 52309926.168262,
                '2018' : 5046714.41,
                },
      },
    "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 3.697,
            "ngen_weight": {
                '2016' : 33378088.32,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 0.3708, #0.4062,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 560315.1334,
                '2018' : 580758.845,
                },
      },
    "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 0.3708, #0.4062,
            "ngen_weight": {
                '2016' : 569424.145,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8": {
      "XS": 61526.7,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 33043732.0,
                '2018' : 70962105.741,
                },
      },
    "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": {
      "XS": 61526.7,
            "ngen_weight": {
                '2016' : 3481160029752.16,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "WW_TuneCP5_13TeV-pythia8": {
      "XS": 118.7,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 7791560.886900296,
                '2018' : 7811153.439,
                },
      },
    "WW_TuneCUETP8M1_13TeV-pythia8": {
      "XS": 118.7,
            "ngen_weight": {
                '2016' : 994032.181,
                '2017' : 1,
                '2018' : 1.,
                },
      },
    "WZ_TuneCP5_13TeV-pythia8": {
      "XS": 65.5443,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 3928630.0,
                '2018' : 3884167.004,
                },
      },
    "WZ_TuneCUETP8M1_13TeV-pythia8": {
      "XS": 65.5443,
            "ngen_weight": {
                '2016' : 1000000.0,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
    "ZZ_TuneCP5_13TeV-pythia8": {
      "XS": 15.8274,
            "ngen_weight": {
                '2016' : 0,
                '2017' : 1949768.0,
                '2018' : 1978776.751,
                },
      },
    "ZZ_TuneCUETP8M1_13TeV-pythia8": {
      "XS": 15.8274,
            "ngen_weight": {
                '2016' : 990064.0,
                '2017' : 1.,
                '2018' : 1.,
                },
      },
}


############################################################### Histograms ########################################################

histogram_settings = {
  'nleps'           : (0,5,6),
  'njets'           : (0,25,26),
  'nfatjets'        : (0,10,11),
  'met'             : (0,2000,201),
  'leading_jet_pt'  : (0,1500,151),
  'leading_jet_eta' : (-4,4,81),
  'leadAK8JetMass'  : (0,300,61),
  'leadAK8JetPt'    : (0,1500,151),
  'leadAK8JetEta'   : (-4,4,81),
  'leadAK8JetHbb'   : (0,1,41),
  'leadAK8JetTau21' : (0,1,41),
  'lepton_pt'       : (0,2000,201),
  'lepton_eta'      : (-4,4,81),
  'hadWPt'          : (0,300,61),
  'hadWEta'         : (-4,4,81),
  'hadWMass'        : (0,300,61),
  'lepWPt'          : (0,300,61),
  'lepWEta'         : (-4,4,81),
  'lepWMass'        : (0,300,61),
  'deltaRlepWHiggs' : (0,5,31),
  'deltaRhadWHiggs' : (0,5,31),
  'PV_npvsGood'     : (0,101,102),
  'pu_weights'      : (0,2,21),
#    "nleps" : (0,10,11),
#    "btags" : (0,8,9),
#    "leading_jet_pt" : (0,500,31),
#    "leading_jet_eta" : (-2.4,2.4,31),
#    "leading_lepton_pt" : (0,500,31),
#    "leading_lepton_eta" : (-2.4,2.4,31),
#    "leading_bjet_pt" : (0,500,31),
#    "leading_bjet_eta" : (-2.4,2.4,31),
#    "subleading_bjet_pt" : (0,500,31),
#    "subleading_bjet_eta" : (-2.4,2.4,31),
#
#    "higgs_pt": (0,500,31),
#    "higgs_eta": (-2.4,2.4,31),
#    "top_pt" : (0,500,31),
#    "top_eta": (-2.4,2.4,31),
#    "nbbtags": (0,4,5),
#    "ntop_candidates": (0,5,6),
#    "nWH_candidates": (0,5,6),
#    "leading_fatjet_pt": (200,500,31),
#    "leading_fatjet_eta": (-2.4,2.4,31),
#    "leading_fatjet_mass": (0,300,31),
#    "leading_fatjet_SDmass": (0,300,31),
#    "subleading_fatjet_pt": (200,500,31),
#    "subleading_fatjet_mass": (0,300,31),
#    "subleading_fatjet_SDmass": (0,300,31),
#    "leading_WHcandidate_SDmass": (0,300,31),
#    "leading_topcandidate_SDmass": (0,300,31),
#    "tau32_fatjets": (0,1,31),
#    "tau32_topcandidates": (0,1,31),
#    "tau32_WHcandidates": (0,1,31),
#    "tau21_fatjets": (0,1,31),
#    "tau21_topcandidates": (0,1,31),
#    "tau21_WHcandidates": (0,1,31),
#
#    "best_higgs_candidate__pt": (200, 500,31),
#    "best_higgs_candidate__msoftdrop": (0, 250, 26),
#    "best_higgs_candidate__tau21": (-1, 1, 31),
#    "best_higgs_candidate__"+eraDependentParameters["2017"]["bbtagging algorithm"]: (-1, 1, 31),
#
#    "best_W_candidate_boosted_pt": (200, 500,31),
#    "best_W_candidate_boosted_msoftdrop": (0, 250, 26),
#    "best_W_candidate_boosted_tau21": (-1, 1, 31),
#
#    "best_top_candidate_boosted_pt": (200, 500,31),
#    "best_top_candidate_boosted_msoftdrop": (0, 250, 26),
#    "best_top_candidate_boosted_tau21": (-1, 1, 31),
#
#    "best_W_candidate_resolved_pt": (0, 250, 26),
#    "best_W_candidate_resolved_mass": (0, 250, 26),
#
#    "leptonic_W_candidate_pt": (0, 250, 26),
#    "leptonic_W_candidate_mass": (0, 250, 26),
#
#    "deltaR_Wlep_WhadResolved": (0, 5, 31),
#    "deltaR_Wlep_WhadBoosted": (0, 5, 31),
#    "deltaR_Wlep_H": (0, 5, 31),
#    "deltaR_H_WhadBoosted": (0,5, 31),
#    "deltaR_H_WhadResolved": (0,5, 31),
}

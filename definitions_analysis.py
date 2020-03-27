

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
    "tau21cut": {
        '2016' : 0.35,
        '2017' : 0.45,
        '2018' : 0.45,
        },
  },
  "W": {
    'min_mass': 65,
    'max_mass': 105
  }
}

eraDependentParameters = {
    "2016" : {
        "lumi":  36773.0,
        "lumimask": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt",
        "pu_corrections_file" : "data/PileupData_GoldenJSON_Full2016.root",
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
        "lumimask": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt",
        "pu_corrections_file" : "data/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root",
        "corrections" : [
            "el_triggerSF SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF data/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5_0_histo.root",
            "el_recoSF EGamma_SF2D data/egammaEffi_EGM2D_runBCDEF_passingRECO_histo.root",
            "el_idSF EGamma_SF2D data/egammaEffi_EGM2D_runBCDEF_passingTight94X_histo.root",
            "mu_triggerSF IsoMu27_PtEtaBins/pt_abseta_ratio data/EfficienciesAndSF_RunBtoF_Nov17Nov2017_histo.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/RunBCDEF_SF_ISO_histo.root",
            "mu_idSF NUM_TightID_DEN_genTracks_pt_abseta data/RunBCDEF_SF_ID_histo.root",
            "BTagSF * data/deepCSV_sfs_v2_btag.csv"
        ],
        "btagging algorithm" : "btagDeepFlavB",#"btagDeepB",
        "btagging WP" : 0.3033, # 0.4941, # medium working point for btagDeepB
        "bbtagging algorithm" :  "deepTagMD_bbvsLight",#"deepTagMD_HbbvsQCD", #"btagDDBvL",
        "bbtagging WP" : 0.8695, ### https://indico.cern.ch/event/853828/contributions/3723593/attachments/1977626/3292045/lg-btv-deepak8v2-sf-20200127.pdf
    },
    "2018" : {
        "lumi":  58830.0,
        "lumimask": "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt",
        "pu_corrections_file" : "data/PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root",
        "corrections" : [
            "el_triggerSF SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF data/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5_0.root",
            "el_recoSF EGamma_SF2D data/2018_ElectronTight.root",
            "el_idSF EGamma_SF2D data/egammaEffi_txt_EGM2D_updatedAll.root",
            "mu_triggerSF IsoMu24_PtEtaBins/pt_abseta_ratio data/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/MuonID_2018_RunABCD_SF_ISO.root",
            "mu_idSF NUM_TightID_DEN_TrackerMuons_pt_abseta data/MuonID_2018_RunABCD_SF_ID.root",
            "BTagSF * data/DeepCSV_Moriond17_B_H.csv"
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
    "ttHTobb": {
            "process": "ttHTobb",
            "XS": 0.2953,
            #"ngen_weight": 2410620.0644499995 #reduced file list
            "ngen_weight": {
                '2016' : 5253482.85,
                '2017' : 4216319.31,
                '2018' : 5046714.41,
                },
            },
    "THW": {
      "XS":  0.01517,#0.002879,
            "ngen_weight": {
                '2016' : 4989133.86,
                '2017' : 4714331.00,
                '2018' : 14971606.14,
                },
      },
    "ttHToNonbb": {
            "process": "ttHToNonbb",
            "XS": 0.2118,
            #"ngen_weight": 913045.7391360003 #reduced file list
            "ngen_weight": {
                '2016' : 5248991.57,
                '2017' : 3095197.81,
                '2018' : 3963935.78,
                },
            },

    "TTToSemiLeptonic": {
            "XS": 365.4574,
            "ngen_weight": {
                '2016' : 32366940321.33,
                '2017' : 59457024911.66,
                '2018' : 60050384119.23,
                },
            },
    "TTTo2L2Nu": {
            "XS": 88.3419,
            "ngen_weight": {
                '2016' : 4784620999.11,
                '2017' : 648729877.29,
                '2018' : 4622080044.95,
                },
            },
    "TTToHadronic": {
            "XS": 377.9607,
            "ngen_weight": {
                '2016' : 21500086465.24,
                '2017' : 61932449366.28,
                '2018' : 62639466237.,
                },
            },
    "TTZToQQ": {
      "XS": 0.6012,# 0.5297,
            "ngen_weight": {
                '2016' : 396340.93,
                '2017' : 4564905.23,
                '2018' : 4534202.12,
                },
      },
    "ST_s-channel_4f_leptonDecays": {
      "process": "ST_s-channel",
      "XS": 3.36,#3.702224,
      #"ngen_weight": 24856809.513425056 #reduced file list
            "ngen_weight": {
                '2016' : 36768937.25,
                '2017' : 37052021.59,
                '2018' : 74634736.73,
                },
      },
    "ST_tW_antitop": {
      "process": "ST_tW_antitop",
      "XS": 35.85,
      #"ngen_weight": 182291193.36093727 #reduced file list
            "ngen_weight": {
                '2016' : 174109580.67,
                '2017' : 279005351.85,
                '2018' : 266470421.96,
                },
      },
    "ST_tW_top": {
      "process": "ST_tW_top",
      "XS": 35.85,
      #"ngen_weight": 241590614.9098064 #reduced file list
            "ngen_weight": {
                '2016' : 173908712.95,
                '2017' : 272081073.53,
                '2018' : 334874722.20,
                },
      },
    "ST_t-channel_antitop": {
      "process": "ST_t-channel_antitop",
      "XS": 80.95,
            "ngen_weight": {
                '2016' : 17771478.65,
                '2017' : 64689262.55,
                '2018' : 5125996535.38,
                },
      },
    "ST_t-channel_top": {
      "process": "ST_t-channel_top",
      "XS": 136.02,
            "ngen_weight": {
                '2016' : 67975483.38,
                '2017' : 5982064.0,
                '2018' : 16603455266.97,
                },
      },
    "TTGJets": {
      "XS": 3.697,
            "ngen_weight": {
                '2016' : 67622406.44,
                '2017' : 62364926.69,
                '2018' : 33778755.65,
                },
      },
    "TTWJetsToQQ": {
      "XS": 0.3708, #0.4062,
            "ngen_weight": {
                '2016' : 569424.14,
                '2017' : 560315.13,
                '2018' : 580758.84,
                },
      },
    "WJetsToLNu": {
      "XS": 61526.7,
            "ngen_weight": {
                '2016' : 57402435.0,
                '2017' : 107612500.0,
                '2018' : 70389866.80,
                },
      },
    "WJetsToLNu_HT-200To400": {
      "XS": 409.3,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 21192211.0,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_HT-400To600": {
      "XS": 57.91,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 14250114.0,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_HT-600To800": {
      "XS": 12.93,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 21582309.0,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_HT-800To1200": {
      "XS": 5.395,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 20272990.0,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_HT-1200To2500": {
      "XS": 1.081,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 19991892.0,
                '2018' : 1.,
                },
      },
    "WJetsToLNu_HT-2500ToInf": {
      "XS": 0.008060,
            "ngen_weight": {
                '2016' : 1.,
                '2017' : 20629585.0,
                '2018' : 1.,
                },
      },
    "WW": {
      "XS": 118.7,
            "ngen_weight": {
                '2016' : 6988278.14,
                '2017' : 7765891.02,
                '2018' : 7846135.95,
                },
      },
    "WZ": {
      "XS": 65.5443,
            "ngen_weight": {
                '2016' : 2997571.0,
                '2017' : 3928630.0,
                '2018' : 3884167.00,
                },
      },
    "ZZ": {
      "XS": 15.8274,
            "ngen_weight": {
                '2016' : 998034.0,
                '2017' : 1949768.0,
                '2018' : 1978776.75,
                },
      },
#    "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
#            "XS": 365.4574,
#            "ngen_weight": {
#                '2016' : 32366940321.33,
#                '2017' : 720253370.0403845,
#                '2018' : 1.,
#                },
#            "nanoAOD" : {
#                '2016' : '',
#                '2017' : '',
#                '2018' : '',
#                },
#            },
#    "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8": {
#            "XS": 88.3419,
#            "ngen_weight": {
#                '2016' : 4746866164.44,
#                '2017' : 283000430.5968169,
#                '2018' : 1.,
#                },
#            "nanoAOD" : {
#                '2016' : '',
#                '2017' : '',
#                '2018' : '',
#                },
#            },
#    "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
#            "XS": 377.9607,
#            "ngen_weight": {
#                '2016' : 21500086465.24,
#                '2017' : 1647945788.3386502,
#                '2018' : 1.,
#                },
#            "nanoAOD" : {
#                '2016' : '',
#                '2017' : '',
#                '2018' : '',
#                },
#            },
}


############################################################### Histograms ########################################################

histogram_settings = {
  'nleps'             : (0,5,6),
  'nelectrons'        : (0,5,6),
  'nmuons'            : (0,5,6),
  'njets'             : (0,25,26),
  'btags'             : (0,20,21),
  'nfatjets'          : (0,10,11),
  'met'               : (0,2000,201),
  'leading_jet_pt'    : (0,1500,151),
  'leading_jet_eta'   : (-4,4,81),
  'leadAK8JetMass'    : (0,300,61),
  'leadAK8JetPt'      : (0,1500,151),
  'leadAK8JetEta'     : (-4,4,81),
  'leadAK8JetHbb'     : (0,1,41),
  'leadAK8JetTau21'   : (0,1,41),
  'lepton_pt'         : (0,2000,201),
  'lepton_eta'        : (-4,4,81),
  'hadWPt'            : (0,300,61),
  'hadWEta'           : (-4,4,81),
  'hadWMass'          : (0,300,61),
  'lepWPt'            : (0,300,61),
  'lepWEta'           : (-4,4,81),
  'lepWMass'          : (0,300,61),
  'deltaRlepWHiggs'   : (0,5,31),
  'deltaRhadWHiggs'   : (0,5,31),
  'deltaRHiggsLepton' : (0,5,31),
  'PV_npvsGood'       : (0,101,102),
  'pu_weights'        : (0,2,21),
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

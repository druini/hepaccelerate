

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
}

eraDependentParameters = {
    "2016" : {
        "lumi":  35922.0,
        "lumimask": "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt",
        "pu_corrections_file" : "data/puData2017_withVar.root",
        "corrections" : [
            "el_triggerSF Ele27_WPTight_Gsf data/TriggerSF_Run2016All_v1.root",
            "el_recoSF EGamma_SF2D data/egammaEffi.txt_EGM2D.root",
            "el_idSF EGamma_SF2D data/egammaEffi.txt_EGM2D.root",
            "mu_triggerSF IsoMu27_PtEtaBins/pt_abseta_ratio data/EfficienciesAndSF_RunBtoF_Nov17Nov2017.histo.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/RunBCDEF_SF_ISO.histo.root",
            "mu_idSF NUM_TightID_DEN_genTracks_pt_abseta data/RunBCDEF_SF_ID.histo.root",
            "BTagSF * data/DeepCSV_Moriond17_B_H.csv"
        ]
    },
    "2017" : {
        "lumi":  41529.0,
        "lumimask": "data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt",
        "pu_corrections_file" : "data/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root",
        "corrections" : [
            "el_triggerSF SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF ./data/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5.0.histo.root",
            "el_recoSF EGamma_SF2D ./data/egammaEffi_EGM2D_runBCDEF_passingRECO.histo.root",
            "el_idSF EGamma_SF2D ./data/egammaEffi_EGM2D_runBCDEF_passingTight94X.histo.root",
            "mu_triggerSF IsoMu27_PtEtaBins/pt_abseta_ratio ./data/EfficienciesAndSF_RunBtoF_Nov17Nov2017.histo.root",
            "mu_isoSF NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta ./data/RunBCDEF_SF_ISO.histo.root",
            "mu_idSF NUM_TightID_DEN_genTracks_pt_abseta ./data/RunBCDEF_SF_ID.histo.root",
            "BTagSF * ./data/deepCSV_sfs_v2.btag.csv"
        ],
        "btagging algorithm" : "btagDeepFlavB",#"btagDeepB",
        "btagging WP" : 0.3033, # 0.4941, # medium working point for btagDeepB
        "bbtagging algorithm" : "deepTagMD_ZHbbvsQCD", # "deepTagMD_HbbvsQCD",#"btagDDBvL",
        "bbtagging WP" : 0.8945 #0.6795 # medium 2 working point for DeepDoubleB tagger
    }

}


#################################### Samples info ##############################################################################

samples_info = {
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_test": {
            "process": "ttHTobb",
            "XS": 0.29533504, #0.2934045,
            #"ngen_weight": 2410620.0644499995 #reduced file list
            "ngen_weight": 4216319.315883999
            },
    "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8": {
            "process": "ttHTobb",
            "XS": 0.29533504, #0.2934045,
            #"ngen_weight": 2410620.0644499995 #reduced file list
            "ngen_weight": 4216319.315883999
            },
    "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 365.45736135,
            "ngen_weight": 720253370.0403845
            },
    "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8": {
            "process": "ttHToNonbb",
            "XS": 0.21176496, #0.2150955,
            #"ngen_weight": 913045.7391360003 #reduced file list 
            "ngen_weight": 3095197.8117420007 
            },
    "TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 88.341903326,
            "ngen_weight": 283000430.5968169
            },
    "TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8": {
            "XS": 377.9607353256,
            "ngen_weight": 1647945788.3386502
            },

    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 88.34190333, #88.341903326,
            #"ngen_weight": 411032366.5312147 #reduced file list
            "ngen_weight": 643135033.87298 #283000430.5968169
            },
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 377.9607353, #377.9607353256,
            #"ngen_weight": 6141680.0 #reduced file list
            "ngen_weight": 41084368.0 #1647945788.3386502
            },
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8": {
            "XS": 365.4573613, #365.45736135,
            #"ngen_weight": 1394703463.7353685 #reduced file list
            "ngen_weight": 10758424549.765732 #720253370.0403845
            },
    "TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8": {
      "XS": 0.5297,
      "ngen_weight": 383062.0686438
      },
    "ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8": {
      "process": "ST_s-channel",
      "XS": 3.702224,
      #"ngen_weight": 24856809.513425056 #reduced file list
      "ngen_weight": 36781553.92694208
      },
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": {
      "process": "ST_tW_antitop",
      "XS": 35.85,
      #"ngen_weight": 182291193.36093727 #reduced file list
      "ngen_weight": 270762750.1725247
      },
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": {
      "process": "ST_tW_top",
      "XS": 35.85,
      #"ngen_weight": 241590614.9098064 #reduced file list
      "ngen_weight": 277241050.84022224
      },
    "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8": {
      "process": "ST_t-channel",
      "XS": 136.02,
      "ngen_weight": 5863722.0
      },
    "THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8": {
      "XS": 0.002879,
      #"ngen_weight": 2423222.0 #reduced file list
      "ngen_weight": 4714331.0
      },
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 3.697,
      #"ngen_weight": 29351538.669739988 #reduced file list
      "ngen_weight": 52309926.168262
      },
    "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8": {
      "XS": 0.4062,
      "ngen_weight": 560315.13342
      },
    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8": {
      "XS": 61526.7,
      #"ngen_weight": 6652665.0 #reduced file list
      "ngen_weight": 33043732.0
      },
    "WW_TuneCP5_13TeV-pythia8": {
      "XS": 118.7,
      "ngen_weight": 7791560.886900296
      },
}


############################################################### Histograms ########################################################

histogram_settings = {
  'nleps'          : (0,5,6),
  'njets'          : (0,25,26),
  'nfatjets'       : (0,10,11),
  'met'            : (0,2000,201),
  'leadAK8JetMass' : (0,300,61),
  'leadAK8JetPt'   : (0,1500,151),
  'leadAK8JetEta'  : (-4,4,81),
  'leadAK8JetHbb'  : (0,1,41),
  'leadAK8JetTau21': (0,1,41),
  'lepton_pt'      : (0,2000,201),
  'lepton_eta'     : (-4,4,81),
  'hadWPt'         : (0,300,61),
  'hadWEta'        : (-4,4,81),
  'hadWMass'       : (0,300,61),
  'lepWPt'         : (0,300,61),
  'lepWEta'        : (-4,4,81),
  'lepWMass'       : (0,300,61),
  'deltaRlepWHiggs': (0,5,31),
  'deltaRhadWHiggs': (0,5,31),
  'PV_npvsGood'    : (0,101,102),
  'pu_weights'     : (0,2,21),
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

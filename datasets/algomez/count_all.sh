echo "start running python script"

samples=(
/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/algomez-NANOAOD_v02-3ea2ff745e1084ea23260bd2ac726434/USER
/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/algomez-NANOAOD_v02-3ea2ff745e1084ea23260bd2ac726434/USER
/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/algomez-NANOAOD_v02-3ea2ff745e1084ea23260bd2ac726434/USER
#/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/algomez-NANOAOD_v02-3ea2ff745e1084ea23260bd2ac726434/USER
/ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
#/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
#/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
#/WW_TuneCP5_13TeV-pythia8/algomez-NANOAOD_v02-5157e087a222b5255c63dabe0cebaee6/USER
)
for sample in ${samples[@]}; do
  name=`echo $sample | awk -F'/' '{print $2}'`
  echo "counting sample $name"
  python ../../counts.py --filelist $name.txt >> weights.txt
  #PYTHONPATH=hepaccelerate:coffea:. python3 run_analysis.py --filelist datasets/RunIIFall17NanoAODv4/$sample.txt --sample $sample  --outdir results_studyBoostedStatistics/mc --from-cache --boosted --use-cuda #--DNN save-arrays --cache-location /scratch/druini/cache/ 
  #python run_analysis.py --filelist datasets/RunIIFall17NanoAODv4/$sample.txt --sample $sample  --cache-location /scratch/druini/cache/ --outdir /scratch/druini/results_studyBoostedStatistics/mc --from-cache --boosted #--DNN save-arrays
done

echo "done"

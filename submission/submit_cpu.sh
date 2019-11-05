echo "start running python script"

mc_samples=(
ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8
#ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8
#TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
#TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
#TTToHadronic_TuneCP5_13TeV-powheg-pythia8
#TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8
#ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8
#ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8
#ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8
#ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8
#THW_ctcvcp_5f_Hincl_13TeV_madgraph_pythia8
#TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8
#TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
#WW_TuneCP5_13TeV-pythia8
)
for sample in ${mc_samples[@]}; do
  echo "submitting sample $sample"
  #source submission/cpu_job.sh $sample
  sbatch submission/cpu_job.sh $sample
done

echo "done"

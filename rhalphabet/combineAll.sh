#todo=(
#mc_msd95to150_msdbin5_pt2bin_polyDegs12
#mc_msd95to155_msdbin5_pt2bin_polyDegs12
#mc_msd95to160_msdbin5_pt2bin_polyDegs12
#mc_msd95to165_msdbin5_pt2bin_polyDegs12
#mc_msd95to170_msdbin5_pt2bin_polyDegs12
#)
#for dir in ${todo[@]}; do
for dir in ./*/     # list directories in the form "/tmp/dirname/"
do
  dir=${dir%*/}      # remove the trailing "/"
  dir=${dir##*/}    # print everything after the final "/"
  echo ${dir}    # print everything after the final "/"
  cd ${dir}
  ls *png | wc -l
#  source build.sh
#  combine -M FitDiagnostics --rMin=-50 --rMax=50 --saveNormalizations --plot --saveShapes --saveWithUncertainties  -v 4 testModel_combined.txt
#  PostFitShapesFromWorkspace -w testModel_combined.root -o shapes_b.root --print --postfit --sampling -f fitDiagnostics.root:fit_b
#  PostFitShapesFromWorkspace -w testModel_combined.root -o shapes_s.root --print --postfit --sampling -f fitDiagnostics.root:fit_s
  cd ..
done

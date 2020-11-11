#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/ttH/CMSSW_10_2_13/src/
eval `scramv1 runtime -sh`
#cmsenv
cd /afs/cern.ch/work/a/algomez/ttH/CMSSW_10_6_5/src/TTH/Analyzer/hepaccelerate/rhalphabet/

#### Run all the datacards and GOF first
for j in Bern Cheb; do
    for k in 1 5; do
        for i in 2 3 4 5 6; do
            python my_rhalphalib.py --msd_start 90 --msd_stop 160 -s met20_btagDDBvL086 --polyDegPt $i -y ${1} --sig-and-bkg -v v13 -r $k --simpleFit --pdf $j -j /eos/home-a/algomez/tmpFiles/hepacc/results/
            if [[ $j == *"Bern"* ]]; then
                datacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_BernpolyDegs${i}/ttHbb_r-20to20.txt
            else
                datacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_ChebpolyDegs${i}/ttHbb_r-20to20.txt
            fi
            python bkgEstTests.py -M GoodnessOfFit -t 200 -a saturated  --selection met20_btagDDBvL086 -y ${1} --msd_start 90 --msd_stop 160  -d ${datacard} -s $RANDOM --pdf1 $j #--rMin 0 --rMax 20
            python bkgEstTests.py -M GoodnessOfFit -t 200 -a KS --selection met20_btagDDBvL086 -y ${1} --msd_start 90 --msd_stop 160  -d ${datacard} -s $RANDOM --pdf1 $j #--rMin 0 --rMax 20

        done
    done
done

for j in Bern Cheb; do
    for k in 1 5; do
        #for i in 3 4; do
        for i in 2 3 4 5; do
            if [[ $j == *"Bern"* ]]; then
                datacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_BernpolyDegs${i}/ttHbb_r-20to20.txt
                altdatacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_BernpolyDegs$(($i+1))/ttHbb_r-20to20.txt
            else
                datacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_ChebpolyDegs${i}/ttHbb_r-20to20.txt
                altdatacard=output/${1}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_ChebpolyDegs$(($i+1))/ttHbb_r-20to20.txt
            fi
            python bkgEstTests.py -M FTest --seed $RANDOM -t 500 --pt1 ${i} --pdf1 ${j} --pt2 $(($i+1)) --n $((70/$k)) --pdf2 ${j} --selection met20_btagDDBvL086 -y ${1} --msd_start 90 --msd_stop 160 -d $datacard --datacard-alt $altdatacard #--rMin -20 --rMax 20
        done
    done
done


#for y in 2016 2017 2018; do
#    for q in 1 2 3 4 5 6 7 8 9 10; do
#        for k in 1 5; do
#            for j in 0 1 3; do
#                datacard=output/${y}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_polyDegs4/ttHbb_r-20to20.txt
#                altdatacard=output/${y}/v13/met20_btagDDBvL086/mcSB_msd90to160_msdbin${k}_exppolyDegs3/ttHbb_r-20to20.txt
#                python bkgEstTests.py -M Bias --datacard-alt ${datacard} -t 100 -d ${datacard} --seed $RANDOM --rMin -20 --rMax 20 --toysFrequentist --selection met20_btagDDBvL086 -y ${y} --msd_start 90 --msd_stop 160 --p1 3 --pdf1 Bern --p2 3 --pdf2 Bern -r ${j}
#                python bkgEstTests.py -M Bias --datacard-alt ${altdatacard} -t 100 -d ${altdatacard} --seed $RANDOM --rMin -20 --rMax 20 --toysFrequentist --selection met20_btagDDBvL086 -y ${y} --msd_start 90 --msd_stop 160 --p1 2 --pdf1 Cheb --p2 2 --pdf2 Cheb -r ${j}
#                python bkgEstTests.py -M Bias --datacard-alt ${datacard} -t 100 -d ${altdatacard} --seed $RANDOM --rMin -20 --rMax 20 --toysFrequentist --selection met20_btagDDBvL086 -y ${y} --msd_start 90 --msd_stop 160 --p1 2 --pdf1 Cheb --p2 3 --pdf2 Bern -r ${j}
#            done
#        done
#    done
#done

#### Run all the datacards and GOF first
###for i in 0,1 0,2 1,0 1,1 1,2 1,3; do
#for i in 1,2 1,3 2,2; do
##for i in 1,2; do
#    IFS=","
#    set -- $i
#    echo $1 $2
#    #python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_btagDDBvL086 --polyDegPt $1 --polyDegRho $2 -y 2017 --sig-and-bkg --poly-limit 1000 -v v13 -r 5 --nptbins 2 --rMin 0 --rMax 20
#    python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_btagDDBvL086 --polyDegPt $1 --polyDegRho $2 -y 2017 --sig-and-bkg --poly-limit 1000 -v v13 -r 5 --nptbins 2 --rMin 0 --pdf exp
#done
#
##for i in 1,2 1,3 2,2; do
#for i in 1,2; do
#    IFS=","
#    set -- $i
#    echo $1 $2
#    #out=output/2017/v13/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}${2}/ttHbb_sig_polylims-1000to1000_combined.root
#    out=output/2017/v13/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs${1}${2}/ttHbb_sig_polylims-1000to1000_combined.root
#    if [ -f "$out" ]; then
#        python bkgEstTests.py -M GoodnessOfFit -t 100 -a saturated  --selection met20_btagDDBvL086 --pt1 $1 --rho1 $2 -y 2017 --msd_start 100 --msd_stop 150 --nmsdbins 10 --rMin 0 --rMax 20 -d $out  -s $RANDOM
#        python bkgEstTests.py -M FTest --seed $RANDOM -t 100 --pt1 $1 --rho1 $2 --pdf1 poly --rho2 $(($2 + 1)) --pt2 $1 --pdf2 poly --selection met20_btagDDBvL086 -y 2017 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d ${out} --datacard-alt output/2017/v13/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}$(($2 + 1))/ttHbb_sig_polylims-1000to1000_combined.root --rMin -100 --rMax 100
#        python bkgEstTests.py -M FTest --seed $RANDOM -t 100 --pt1 $1 --rho1 $2 --pdf1 poly --pt2 $(($1 + 1)) --rho2 $2 --pdf2 poly --selection met20_btagDDBvL086 -y 2017 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d ${out} --datacard-alt output/2017/v13/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs$(($1 + 1))${2}/ttHbb_sig_polylims-1000to1000_combined.root --rMin -100 --rMax 100
#        #python bkgEstTests.py -M GoodnessOfFit -t 100 -a saturated  --selection met20_btagDDBvL086 --pt1 $1 --rho1 $2 -y 2017 --msd_start 90 --msd_stop 150 --nmsdbins 10 --rMin 0 --rMax 20 -d $out  -s $RANDOM --pdf1 exp
#    fi
#done






#python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_btagDDBvL086 --polyDegPt $1 --polyDegRho $2 -y $3 --sig-and-bkg --rMin 0 --rMax 20 --poly-limit $5 -v $6 -r 5 #--runPrefit
#python bkgEstTests.py -M GoodnessOfFit -t $4 -a saturated  --selection met20_btagDDBvL086 --pt1 $1 --rho1 $2 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 --rMin -1 --rMax 1 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root -s 216741
#python bkgEstTests.py -M GoodnessOfFit -t $4 -a KS --selection met20_btagDDBvL086 --pt1 $1 --rho1 $2 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 --rMin -1 --rMax 1 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root -s 216741
#python bkgEstTests.py -M FTest --seed $RANDOM -t $4 --pt1 $1 --rho1 $2 --pdf1 poly --pt2 $(($1 + 1)) --rho2 $2  --pdf2 poly --selection met20_btagDDBvL086 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root --datacard-alt output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs$(($1 + 1))${2}/ttHbb_sig_polylims-$5to$5_combined.root --rMin -1 --rMax 1
#python bkgEstTests.py -M FTest --seed $RANDOM -t $4 --pt1 $1 --rho1 $2 --pdf1 poly --pt2 $1 --rho2 $(($2 + 1)) --pdf2 poly --selection met20_btagDDBvL086 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root --datacard-alt output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_polyDegs${1}$(($2 + 1))/ttHbb_sig_polylims-$5to$5_combined.root --rMin -1 --rMax 1


#python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_btagDDBvL086 --polyDegPt $1 --polyDegRho $2 -y $3 --sig-and-bkg --rMin 0 --rMax 1 --poly-limit $5 -v $6 --pdf exp #--runPrefit #--runImpacts
#python bkgEstTests.py -M GoodnessOfFit -t $4 -a saturated  --selection met20_btagDDBvL086 --pt1 $1 --rho1 $2 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root -s 216741 --rMin -1 --rMax 1
#python bkgEstTests.py -M FTest --seed 216741 -t $4 --pt1 $1 --rho1 $2 --pdf1 exp --pt2 $(($1 + 1)) --rho2 $2  --pdf2 exp --selection met20_btagDDBvL086 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root --datacard-alt output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs$(($1 + 1))${2}/ttHbb_sig_polylims-$5to$5_combined.root --rMin -1 --rMax 1
#python bkgEstTests.py -M FTest --seed 216741 -t $4 --pt1 $1 --rho1 $2 --pdf1 exp --pt2 $1 --rho2 $(($2 + 1)) --pdf2 exp --selection met20_btagDDBvL086 -y $3 --msd_start 100 --msd_stop 150 --nmsdbins 10 -d output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs${1}${2}/ttHbb_sig_polylims-$5to$5_combined.root --datacard-alt output/${3}/${6}/met20_btagDDBvL086/mcSB_msd100to150_msdbin5_pt2bin_exppolyDegs${1}$(($2 + 1))/ttHbb_sig_polylims-$5to$5_combined.root --rMin -1 --rMax 1

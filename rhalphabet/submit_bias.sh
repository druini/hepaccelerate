#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/d/druini/public/combine/CMSSW_10_2_13/src
eval `scramv1 runtime -sh`
#cmsenv
cd /afs/cern.ch/work/d/druini/public/hepaccelerate/rhalphabet

python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt 1 --polyDegRho 2 --rMin -3 --rMax 3 -v $2 -y $3
python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt 1 --polyDegRho 1 --pdf exp --rMin -3 --rMax 3 -v $2 -y $3
python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt 1 --polyDegRho 1 --pdf exp --rMin -3 --rMax 3 -l 3 -v $2 -y $3

python bkgEstTests.py -M Bias --seed 216741 -t 1000 --selection met20_deepTagMD_bbvsLight08695 --rMin -40 --rMax 40 --toysFrequentist --pt1 1 --rho1 2 --pdf2 exp --pt2 1 --rho2 1 -r $1  -v $2 -y $3

python bkgEstTests.py -M Bias --seed 216741 -t 1000 --selection met20_deepTagMD_bbvsLight08695 --rMin -50 --rMax 50 --toysFrequentist --pt2 1 --rho2 2 --pdf1 exp --pt1 1 --rho1 1 -r $1  -v $2 -y $3

python bkgEstTests.py -M Bias --seed 216741 -t 1000 --selection met20_deepTagMD_bbvsLight08695 --rMin -40 --rMax 40 --toysFrequentist --pt1 1 --rho1 2 --pt2 1 --rho2 2 -r $1  -v $2 -y $3

python bkgEstTests.py -M Bias --seed 216741 -t 1000 --selection met20_deepTagMD_bbvsLight08695 --rMin -50 --rMax 50 --toysFrequentist --pdf2 exp --pt2 1 --rho2 1 --pdf1 exp --pt1 1 --rho1 1 -r $1 --poly-limit 3 --polyAlt-limit 3 -v $2 -y $3

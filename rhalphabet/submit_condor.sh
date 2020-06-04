#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/d/druini/public/combine/CMSSW_10_2_13/src
eval `scramv1 runtime -sh`
#cmsenv
cd /afs/cern.ch/work/d/druini/public/hepaccelerate/rhalphabet

python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt $1 --polyDegRho $2
python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt $(($1 + 1)) --polyDegRho $2
python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt $1 --polyDegRho $(($2 + 1))
python my_rhalphalib.py --msd_start 100 --msd_stop 150 -s met20_deepTagMD_bbvsLight08695 --polyDegPt $(($1 + 1)) --polyDegRho $(($2 + 1))

time python bkgEstTests.py -M FTest --seed 216741 -t 1000 --pt1 $1 --rho1 $2 --pt2 $(($1 + 1)) --rho2 $2
time python bkgEstTests.py -M FTest --seed 216741 -t 1000 --pt1 $1 --rho1 $2 --pt2 $1 --rho2 $(($2 + 1))
time python bkgEstTests.py -M FTest --seed 216741 -t 1000 --pt1 $1 --rho1 $2 --pt2 $(($1 + 1)) --rho2 $(($2 + 1))

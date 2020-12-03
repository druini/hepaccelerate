#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/ttH/CMSSW_10_2_13/src/
eval `scramv1 runtime -sh`
#cmsenv
cd /afs/cern.ch/work/a/algomez/ttH/CMSSW_10_6_5/src/TTH/Analyzer/hepaccelerate/rhalphabet/
##### {1}: r value, {2} year, {3} Cheb degree, {4} Bern degree, {5} bin size, {6} ${6} or data

for i in 0 1 2 3 4 5 6 7 8; do
    python bkgEstTests.py -M Bias --datacard-alt output/${2}/v14/met20_btagDDBvL086/${6}_msd90to160_msdbin${5}_BernpolyDegs${4}/ttHbb_r-20to20.txt  -d output/${2}/v14/met20_btagDDBvL086/${6}_msd90to160_msdbin${5}_ChebpolyDegs${3}/ttHbb_r-20to20.txt --seed $RANDOM --rMin -20 --rMax 20 --toysFrequentist --selection met20_btagDDBvL086 -y ${2} --msd_start 90 --msd_stop 160 --pt1 ${3} --pdf1 Cheb --pt2 ${4} --pdf2 Bern -r $1 -t 300
done

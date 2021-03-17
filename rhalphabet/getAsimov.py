from __future__ import print_function
from my_rhalphalib import exec_me
import os, argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--year', default='2017', type=str, help='year to process, in file paths')
    parser.add_argument('-v', '--version', default='v14', help='version, in file paths')
    parser.add_argument('-s', '--selection',  default='met20_btagDDBvL086', help='event selection, in file paths')
    parser.add_argument('--pdf', default='Cheb', choices=['poly','exp', 'Cheb', 'Bern'])
    parser.add_argument('--polyDegPt', default=2, type=int, help='degree of polynomial to fit pt')
    parser.add_argument('-r', '--rebin_factor', default=5, type=int, help='rebin factor for json histograms, default mass bin size is 1GeV')
    parser.add_argument('--statOnly', action='store_true', default=False)
    parser.add_argument('--dryRun', action='store_true')
    args = parser.parse_args()

    folder_base = '/eos/home-d/druini/hepaccelerate/rhalphabet/output'
    folder = os.path.join(folder_base, args.year, args.version, args.selection, 'data_msd100to160_msdbin'+str(args.rebin_factor)+'_'+args.pdf+'polyDegs'+str(args.polyDegPt)+'_ttHTobb'+('_statOnly' if args.statOnly else ''))

    print(folder)

    if not os.path.isfile( os.path.join(folder,'ttHbb_r-20to20.txt') ):
        print('running my_rhalphalib.py')
        cmd = 'python my_rhalphalib.py --msd_start 100 --msd_stop 160 -s '+args.selection+' --pdf '+args.pdf+' --polyDegPt '+str(args.polyDegPt)+' -y '+str(args.year)+' -d -v '+args.version+' -r '+str(args.rebin_factor)+' --simpleFit -j /eos/home-d/druini/hepaccelerate/results_merged --signal ttHTobb -f '+('--statOnly' if args.statOnly else '')
        exec_me(cmd, dryRun=args.dryRun)

    cmd = "text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=.*/(TTH_PTH_GT300):TTH_PTH_GT300[1,-10,10]' -o 5D_Boosted_highpt.root ttHbb_r-20to20.txt --X-allow-no-signal"
    exec_me(cmd, dryRun=args.dryRun, folder=folder)

    cmd = "combine -M FitDiagnostics -m 125 --cminDefaultMinimizerTolerance 1e-2 --cminDefaultMinimizerStrategy 0 --setParameterRanges TTH_PTH_GT300=-10,10 --setParameters TTH_PTH_GT300=1  -t -1 -n _asimov_sig1_r_ttH_0_5D --redefineSignalPOIs TTH_PTH_GT300 5D_Boosted_highpt.root -v 4"
    exec_me(cmd, dryRun=args.dryRun, folder=folder)

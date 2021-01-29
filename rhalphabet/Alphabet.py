#!/usr/bin/env python
'''
./Rhalphabet.py -v v09 && combine -M FitDiagnostics datacard.txt  --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace
'''

from ROOT import *
import time, os, math, sys, copy
from array import array
import argparse
from collections import OrderedDict
import numpy as np
import json
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
from DrawHistogram import jsonToTH1

####gReset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)
gSystem.Load(os.getenv('CMSSW_BASE')+'/lib/'+os.getenv('SCRAM_ARCH')+'/libHiggsAnalysisCombinedLimit.so')

xline = array('d', [0,2000])
yline = array('d', [1,1])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

rhalPtList = [ '250', '300', '350', '400', '500', '2000' ]
#rhalPtList = [ '250','2000' ]

canvas = {}

def setSelection( listSel, xMin=0.65, yMax=0.65, align='right' ):

    for i in range( len( listSel ) ):
        textBox=TLatex()
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        if 'right' in align: textBox.SetTextAlign(31)
        textBox.SetTextFont(62) ### 62 is bold, 42 is normal
        textBox.DrawLatex(xMin, yMax, listSel[i])
        yMax = yMax -0.05

def exec_me(command, dryRun=False, folder=False):
    print(command)
    if not dryRun:
        if folder: os.chdir(folder)
        os.system(command)


######################################################
### mass decorrelation functions
### This was just a test, but maybe needed later
######################################################
def computeDDTMap(working_point,deepboosted_rho_pt):
    """docstring for computeDDTMap based on https://gitlab.cern.ch/abenecke/scripts/blob/master/diboson/taggerstudies/create_2D_decorrelation.py"""
    ddtMap= deepboosted_rho_pt.Project3D("yx");
#    ddtMap.SetName(ddtMapName);
    ddtMap.Reset();
    ddtMap.SetStats(0);
    ddtMap.SetDirectory(0);
    nbins_x = ddtMap.GetNbinsX();
    nbins_y = ddtMap.GetNbinsY();

    for xbin in range(1, nbins_x):
        for ybin in range(1,nbins_y):
            DeepBoostedproj = deepboosted_rho_pt.ProjectionZ("Rho2D", xbin, xbin, ybin, ybin);
            # Protection against default values
            # for bin in range(1,DeepBoostedproj.GetNbinsX()):
                # if DeepBoostedproj.GetBinCenter(bin)<0:
                #     DeepBoostedproj.SetBinContent(bin,0);

            if DeepBoostedproj.Integral() == 0:
                xval = deepboosted_rho_pt.GetXaxis().GetBinCenter(xbin)
                yval = deepboosted_rho_pt.GetYaxis().GetBinCenter(ybin)
#                print "Caution Integral 0!"
                ddtMap.SetBinContent(xbin,ybin,0)
                continue

            wp = array('d',[working_point])
            quantieles = array('d', [0.])
            DeepBoostedproj.GetQuantiles(1,quantieles , wp)
            ddtMap.SetBinContent(xbin,ybin,quantieles[0])

    return ddtMap;


######################################################
### This was just a test, but maybe needed later
def massDecorrelation( name  ):
    """docstring for massDecorrelation based on https://gitlab.cern.ch/abenecke/scripts/blob/master/diboson/taggerstudies/create_2D_decorrelation.py"""

    bkgHistos = OrderedDict()
    for bkgSamples in bkgFiles:
        if not bkgSamples.startswith("TTToSemi"): continue
        newName = name# .replace('___', ivar+'_').replace('__', '_'+side+'_')
        bkgHistos[ bkgSamples+'3D' ] = bkgFiles[ bkgSamples ][0].Get( 'tthbb13/'+newName )
        #bkgHistos[ bkgSamples+'3D' ].SetTitle(bkgSamples+ivar+side)
        bkgHistos[ bkgSamples+'3D' ].Scale( bkgFiles[ bkgSamples ][1] )
        #try: bkgHistos[ ivar+side ].Add( bkgHistos[ bkgSamples+ivar+side ].Clone() )
        #except (KeyError, AttributeError) as e: bkgHistos[ ivar+side ] = bkgHistos[ bkgSamples+ivar+side ].Clone()

        #create 2D projections
        DeepBoosted_v_rho = TH2D("DeepBoosted_v_rho","DeepBoosted_v_rho",14,-7.0,0,40,0.,1.0);
        DeepBoosted_v_pt = TH2D("DeepBoosted_v_pt","DeepBoosted_v_pt",150,0,1500,40,0.,1.0);
        NBins_rho = bkgHistos[ bkgSamples+'3D' ].GetNbinsX();
        NBins_pt = bkgHistos[ bkgSamples+'3D' ].GetNbinsY();
        NBins_DeepBoosted = bkgHistos[ bkgSamples+'3D' ].GetNbinsZ();

        print "Bins(x,y,z): " + str(NBins_rho) + "  ,  " + str(NBins_pt) + "  ,  " + str(NBins_DeepBoosted)

        for rhoBin in range(1,NBins_rho):
            for ptBin in range(1,NBins_pt):
                for DeepBoostedBin in range(1,NBins_DeepBoosted):
                    rho = bkgHistos[ bkgSamples+'3D' ].GetXaxis().GetBinCenter(rhoBin)
                    pt = bkgHistos[ bkgSamples+'3D' ].GetYaxis().GetBinCenter(ptBin)
                    DeepBoosted = bkgHistos[ bkgSamples+'3D' ].GetZaxis().GetBinCenter(DeepBoostedBin)
                    DeepBoosted_v_rho.Fill(rho,DeepBoosted,bkgHistos[ bkgSamples+'3D' ].GetBinContent(rhoBin,ptBin,DeepBoostedBin));
                    DeepBoosted_v_pt.Fill(pt,DeepBoosted,bkgHistos[ bkgSamples+'3D' ].GetBinContent(rhoBin,ptBin,DeepBoostedBin));

        print "Done!"

        canvas['deepBoostedrho'] = TCanvas('deepBoostedrho', 'deepBoostedrho',  10, 10, 750, 750 )
        DeepBoosted_v_rho.Draw("colz")
        canvas['deepBoostedrho'].SaveAs('Plots/deepBoostedrho.png')

        canvas['deepBoostedpt'] = TCanvas('deepBoostedpt', 'deepBoostedpt',  10, 10, 750, 750 )
        DeepBoosted_v_pt.Draw("colz")
        canvas['deepBoostedpt'].SaveAs('Plots/deepBoostedpt.png')

        for xbin in range(1,DeepBoosted_v_rho.GetNbinsX()+1):
            proj = DeepBoosted_v_rho.ProjectionY("_y",xbin,xbin)
            if not proj.Integral() == 0:
                proj.Scale(1/proj.Integral())
            if not xbin%10:
                print "rho bin  " + str(DeepBoosted_v_rho.GetXaxis().GetBinCenter(xbin))
                c1 = TCanvas("c"+str(xbin), "c"+str(xbin), 600, 600);
                gStyle.SetOptStat(kFALSE);
                gStyle.SetPadTickY(1);
                gStyle.SetPadTickX(1);
                gStyle.SetLegendBorderSize(0);
                gPad.SetBottomMargin(.2);
                gPad.SetRightMargin(.2);

                leg=TLegend(0.2,0.7,0.4,0.9,"","brNDC")
                leg.SetHeader("#splitline{rho = " + str(DeepBoosted_v_rho.GetXaxis().GetBinCenter(xbin))+"}{averaged in pT}")
                leg.SetBorderSize(0);
                leg.SetTextSize(0.035);
                leg.SetFillColor(0);
                leg.SetLineColor(1);
                leg.SetTextFont(42);


                proj.GetXaxis().SetRangeUser(0,1)
                proj.GetXaxis().SetTitle("DeepBoosted")
                proj.GetYaxis().SetTitle("#Delta N /N")
                #proj.GetXaxis().SetLabelSize(ldsize)

                proj.Draw()

                txbin_0p02 = -99
                txbin_0p05 = -99
                first = True
                for bin in range(1,proj.GetNbinsX()+1):
                    inte = proj.Integral(0,bin)

                    if inte >= 0.98:
                        txbin_0p02 = proj.GetBinCenter(bin)
                        break
                    if inte >= 0.95 and first:
                        txbin_0p05 = proj.GetBinCenter(bin)
                        first = False


                line_0p02 = TLine(txbin_0p02,0,txbin_0p02,1)
                line_0p02.SetLineColor(kBlue)
                line_0p02.Draw("same")

                leg.AddEntry(line_0p02,"2% mistag rate","l")

                line_0p05 = TLine(txbin_0p05,0,txbin_0p05,1)
                line_0p05.SetLineColor(kRed)
                line_0p05.Draw("same")

                leg.AddEntry(line_0p05,"5% mistag rate","l")
                leg.Draw()
                c1.SetLogy()
                c1.SaveAs("Plots/Proj"+str(xbin)+".png");

        #calculated Map
        ddt_0p05= computeDDTMap(0.95,bkgHistos[ bkgSamples+'3D' ]);

        ddt_0p05.SetTitle("Simple DeepBoosted-DDT map");
        ddt_0p05.SetTitle("Rho2D");
        c1 = TCanvas("c1", "c1", 600, 600);
        gPad.SetRightMargin(0.2);
        ddt_0p05.GetXaxis().SetRangeUser(-6.0,-1);
        ddt_0p05.GetYaxis().SetRangeUser(0,2000);
        ddt_0p05.GetZaxis().SetRangeUser(0,1);
        ddt_0p05.Draw("colz");
        c1.SaveAs("Plots/DDT_0p05.png");
        ddt_0p05.GetZaxis().SetRangeUser(0.6,1);
        c1.SaveAs("Plots/DDT_0p05_scale.png");
######################################################
######################################################

######################################################
### Functions for rhalphabet
######################################################
def buildPolynomialArray( label0, label1, iNVar0, iNVar1, pXMin, pXMax ):
    """docstring for buildPolynomialArray"""
    polyArray = []
    for i0 in range(iNVar0+1):
        for i1 in range(iNVar1+1):
            pVar = label0+str(i0)+label1+str(i1)
            pRooVar = RooRealVar(pVar, pVar, 0.0, pXMin, pXMax)
            polyArray.append(pRooVar)
    return polyArray


def buildRooPolyArray( iPt, maxPolyPt, iMass, maxPolyMass, lUnity, iVars ):
    """docstring for buildRooPolyArray"""

    print( '-'*20, "MYINFO - buildRooPolyArray")
    # print len(iVars);
    striPt = str(int(iPt))
    striMass = str(int(iMass))

    lPt = RooConstVar( "Var_Pt_"+striPt+"_"+striMass, "Var_Pt_"+striPt+"_"+striMass, iPt )
    lMass = RooConstVar( "Var_Mass_"+striPt+"_"+striMass, "Var_Mass_"+striPt+"_"+striMass, iMass )
    lMassArray = RooArgList()
    lNCount = 0
    for pRVar in range(0, maxPolyMass+1):
        lTmpArray = RooArgList()
        for pVar in range(0, maxPolyPt+1):
            if lNCount == 0: lTmpArray.add(lUnity)  # for the very first constant (e.g. p0r0), just set that to 1
            else: lTmpArray.add(iVars[lNCount])
            lNCount = lNCount + 1
        pLabel = "Var_Pol_Bin_"+striPt+"_"+striMass+"_"+str(pRVar)#+suffix
        pPol = RooPolyVar(pLabel , pLabel , lPt, lTmpArray)
        #pPol.Print()
        lMassArray.add(pPol)

    lLabel = "Var_MassPol_Bin_"+striPt+"_"+striMass#+suffix
    lMassPol = RooPolyVar(lLabel , lLabel , lMass, lMassArray)
    #lMassPol.Print()
    return lMassPol


def generate_bernstein_string(n):
    # x = @(n+1)
    monomials = []
    for v in xrange(0, n+1):
            normalization = 1. * math.factorial(n) / (math.factorial(v) * math.factorial(n - v))
            monomials.append("({} * @{} * (@{}**{}) * ((1.-@{})**{}))".format(normalization, v, n+1, v, n+1, n-v))
    return " + ".join(monomials)

def buildRooPolyRhoArrayBernsteinExp( iPt, maxPolyPt, iRho, maxPolyRho, iQCD, iZero, iVars ):

    print( '-'*20, "MYINFO - buildRooPolyArrayBernsteinExp")

    striPt = str(int(iPt))
    striRho = str(int(iRho))
    lPt = RooConstVar("Var_Pt_"+striPt+"_"+striRho, "Var_Pt_"+striPt+"_"+striRho, (iPt))
    lPt_rescaled = RooConstVar("Var_Pt_rescaled_"+striPt+"_"+striRho, "Var_Pt_rescaled_"+striPt+"_"+striRho, ((iPt - float(rhalPtList[0])) / (float(rhalPtList[-1]) - float(rhalPtList[0]))))
    lRho = RooConstVar("Var_Rho_"+striPt+"_"+striRho, "Var_Rho_"+striPt+"_"+striRho, iRho )
    lRho_rescaled = RooConstVar("Var_Rho_rescaled_"+striPt+"_"+striRho,
"Var_Rho_rescaled_"+striPt+"_"+striRho, ((iRho - (-6)) / ((-2) - (-6))))

    ptPolyString = generate_bernstein_string(maxPolyPt)
    rhoPolyString = generate_bernstein_string(maxPolyRho)

    lRhoArray = RooArgList()
    lNCount = 0
    for pRVar in range(0, maxPolyRho + 1):
        lTmpArray = RooArgList()
        for pVar in range(0, maxPolyPt + 1):
            #if lNCount == 0:
            #    lTmpArray.add(iQCD)  # for the very first constant (e.g. p0r0), just set that to 1
            #else:
            print "lNCount = " + str(lNCount)
            lTmpArray.add(iVars[lNCount])
            print "iVars[lNCount]: ", iVars[lNCount]
            print "iVars[lNCount]"
            iVars[lNCount].Print()
            lNCount = lNCount + 1
        pLabel = "Var_Pol_Bin_exp_"+striPt+"_"+striRho+"_"+str(pRVar)
        lTmpArray.add(lPt_rescaled)
        print "lTmpArray: "
        lTmpArray.Print()
        pPol = RooFormulaVar(pLabel, pLabel, ptPolyString, lTmpArray)
        print "pPol:"
        pPol.Print("V")
        lRhoArray.add(pPol)

    lLabel = "Var_RhoPol_Bin_exp_"+striPt+"_"+striRho
    lRhoArray.add(lRho_rescaled)
    print "lRhoArray: "
    #lRhoArray.Print()
    lRhoPol = RooFormulaVar(lLabel, lLabel, 'exp('+rhoPolyString+')', lRhoArray)
    print('exp('+rhoPolyString+')')
    return lRhoPol
######################################################


######################################################
######################################################
def manualAlphabet( passHisto, passLabel, failHisto, failLabel, passHistoSignal, failHistoSignal, name, outfolder, xmin='', xmax=''):
    """"Run manually alphabet"""

    outName = outfolder+'/'+name+'_'+passLabel+failLabel+'_test.png'

    binWidth = passHisto.GetBinWidth(1)

    legend=TLegend(0.60,0.75,0.90,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    #histo.SetFillColor(48)
    passHisto.SetFillStyle(1001)

    tdrStyle.SetPadRightMargin(0.05)
    canvas[name] = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    legend.AddEntry( passHisto, passLabel, 'l' )
    legend.AddEntry( failHisto, failLabel, 'l' )
    if xmax and xmin: passHisto.GetXaxis().SetRangeUser( xmin, xmax )
    passHisto.GetYaxis().SetTitleOffset(0.90)
    passHisto.SetMaximum( 1.3*max( passHisto.GetMaximum(), failHisto.GetMaximum() )  )
    passHisto.SetLineColor(kBlue)
    failHisto.SetLineColor(kRed)
    passHisto.Draw('hist')
    failHisto.Draw('hist same')
    passHisto.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    legend.Draw()
    canvas[name].SaveAs( outName.replace('_test', '_testPassFailRegions') )


    #### Ratio
    hRatio = passHisto.Clone()
    hRatio.Reset()
    hRatio.Divide( passHisto, failHisto )
    canvas[name+'ratio'] = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    hRatio.Draw('P')
    fitLine = TF1( 'fitLine', 'pol1', xmin, xmax ) #800, 5000)
    hRatio.Fit( 'fitLine', 'MIRS')
    fitLine.Draw("same")
    hRatio.GetYaxis().SetTitle( 'Ratio pass/fail' )
    fitResults = hRatio.GetFunction('fitLine')
    parameters = [ 'p0 = '+str(round(fitResults.GetParameter(0),2))+' #pm '+str(round(fitResults.GetParError(0),2)), 'p1 = '+str(round(fitResults.GetParameter(1),4))+' #pm '+str(round(fitResults.GetParError(1),4)) ]
    setSelection( parameters, 0.7, 0.8, 'left' )
    canvas[name+'ratio'].SaveAs( outName.replace( '_test', '_testRatio' ) )

    ### alphabet
    alphabet = failHisto.Clone()
    alphabet.Multiply( fitLine )
    alphabet.SetLineColor(kViolet)
    alphabet.SetLineWidth(2)

    canvas[name+'alpha'] = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    legend2=TLegend(0.60,0.75,0.90,0.90)
    legend2.SetFillStyle(0)
    legend2.SetTextSize(0.03)
    legend2.AddEntry( alphabet, 'Alphabet estimation', 'l' )
    legend2.AddEntry( passHisto, passLabel, 'l' )
    if xmax and xmin: passHisto.GetXaxis().SetRangeUser( xmin, xmax )
    passHisto.GetYaxis().SetTitleOffset(0.90)
    passHisto.SetMaximum( 1.3*max( alphabet.GetMaximum(), passHisto.GetMaximum() )  )
    passHisto.SetLineColor(kBlue)
    passHisto.SetLineWidth(2)
    passHisto.Draw('histe')
    alphabet.Draw('histe same')
    passHisto.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    legend2.Draw()
    canvas[name+'alpha'].SaveAs( outName.replace( '_test', '_testAlphabet' ) )


    canvas[name] = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    legend3=TLegend(0.60,0.75,0.90,0.90)
    legend3.SetFillStyle(0)
    legend3.SetTextSize(0.03)
    legend3.AddEntry( passHistoSignal, passLabel, 'l' )
    legend3.AddEntry( failHistoSignal, failLabel, 'l' )
    if xmax and xmin: passHistoSignal.GetXaxis().SetRangeUser( xmin, xmax )
    passHistoSignal.GetYaxis().SetTitleOffset(0.90)
    passHistoSignal.SetMaximum( 1.3*max( passHistoSignal.GetMaximum(), failHistoSignal.GetMaximum() )  )
    passHistoSignal.SetLineColor(kRed)
    failHistoSignal.SetLineColor(kBlue)
    passHistoSignal.Draw('hist')
    failHistoSignal.Draw('hist same')
    passHistoSignal.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    #labelAxis( name, passHistoSignal, '' )
    legend.Draw()
    canvas[name].SaveAs( outName.replace( '_test', '_testSignal' ) )
####################################################


######################################################
### Main bkg estimation
######################################################
def BkgEstimation( ): #name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50 ):
    """
    Bkg Estimation with rhalphabet method based on https://github.com/DAZSLE/ZPrimePlusJet/blob/35ca072541e8bf9ebbddca523658aad81fea97bc/fitting/rhalphabet_builder.py#L28
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/SWGuideNonStandardCombineUses#GammaN_with_shapes_using_RooPara
    """

    #### Creating output directory
    pref = ('data' if args.isData else 'mc'+( 'SB' if args.sig_and_bkg else '' ) )
    combineFolder = os.path.join(str(outdir), pref+'_alpha_msd%dto%d_msdbin%d_%spolyDegs%d'%(msd_start,msd_stop,rebin_factor,args.pdf, args.polyDegPt))
    try: os.makedirs(combineFolder)
    except OSError: print('|===>'+combineFolder+' Folder already exist')

    ### Initializing RooFit variables
    msd = RooRealVar( 'msd', 'msd', args.msd_start, args.msd_stop)
    msd.Print()
    rooDict = OrderedDict()     ## dict of roofit objects

    ### Opening plots
    dataOrBkg = 'data_merged' if args.isData else 'background'+( '' if args.year.startswith('2017') else ( '_noQCD_noDY' if args.year.startswith('2016') else '_noDY'  ) )
    templates = {}
    for X in [ 'Pass', 'Fail' ]:
        if X.startswith('Fail'): tmpdir = indir.replace('1btag', 'ex0btags')
        else: tmpdir = indir
        templates[ 'ttH'+X ] = jsonToTH1( tmpdir+'/nominal/out_signal_nominal_merged.json', [ 'leadAK8JetMass_2J2WdeltaR_Pass_pt'+str(ptbins[0])+'to'+str(ptbins[-1]) ] )
        templates[ 'background'+X ] = jsonToTH1( tmpdir+'/nominal/out_'+dataOrBkg+'_nominal_merged.json', [ 'leadAK8JetMass_2J2WdeltaR_Pass_pt'+str(ptbins[0])+'to'+str(ptbins[-1]) ] )

        for unc in uncList:
            for UpDown in [ 'Up', 'Down' ]:
                if args.year.startswith('allyears') and unc.endswith(('2016','2017','2018')): tmpdir = indir.replace('allyears', unc.split('_')[1])
                else: tmpdir = indir
                templates['CMS_ttHbb_'+unc+UpDown+X] = jsonToTH1( indir+'/'+unc+UpDown+'/out_signal_'+unc+UpDown+'_merged.json', [ 'leadAK8JetMass_2J2WdeltaR_Pass_pt'+str(ptbins[0])+'to'+str(ptbins[-1]) ] )
    print(templates)
    if args.rebin_factor > 0:
        for ih in templates: templates[ih].Rebin(args.rebin_factor)

    ##### test plot pass and fail
    manualAlphabet( templates['backgroundPass'], '1btag', templates['backgroundFail'], '0btag', templates['ttHPass'], templates['ttHFail'], 'leadAK8JetMass', combineFolder,  xmin=msd_start , xmax=msd_stop  )
    sys.exit(0)



    ##############################
    bkgEFF = RooRealVar( "bkgeff", "bkgeff", 0.01, 0., 10.)
    bkgEFF.setVal( templates['backgroundPass'].Integral()/templates['backgroundFail'].Integral() )
    bkgEFF.setConstant(True)
    bkgEFF.Print()

    ### polynomial
    polyArgList = RooArgList( )
    for i in range( int(args.polyDegPt)+1 ):
        if args.pdf.startswith('Bern'): rooDict[ 'boosted_bkg_paramX'+str(i) ] = RooRealVar('boosted_bkg_paramX'+str(i), 'boosted_bkg_paramX'+str(i), 1./TMath.Power(1,i), -100., 100. )
        else: rooDict[ 'boosted_bkg_paramX'+str(i) ] = RooRealVar('boosted_bkg_paramX'+str(i), 'boosted_bkg_paramX'+str(i), 1., -1000., 1000. )
        #rooDict[ 'boosted_bkg_paramX'+str(i) ] = RooRealVar('boosted_bkg_paramX'+str(i), 'boosted_bkg_paramX'+str(i), 0.001, -1000000, 1000000 )
        polyArgList.add( rooDict[ 'boosted_bkg_paramX'+str(i) ] )

    rooDict[ 'boosted_bkg_Fail_bins' ] = RooArgList( )
    rooDict[ 'boosted_bkg_Pass_bins' ] = RooArgList( )
    for ibin in range(1, templates['backgroundFail'].GetNbinsX()+1):

        ### Fail workspace
        iCont = templates['backgroundFail'].GetBinContent(ibin)
        iContErr = templates['backgroundFail'].GetBinError(ibin)

        rooDict[ 'boosted_bkg_Fail_bin'+str(ibin) ] = RooRealVar( 'boosted_bkg_Fail_bin'+str(ibin), 'boosted_bkg_Fail_bin'+str(ibin), iCont ) #, max(iCont-(5*iContErr),0), max(iCont+(5*iContErr),0) )   ###### CHECK if 5 is needed
        rooDict[ 'boosted_bkg_Fail_bin'+str(ibin) ].Print()
        rooDict[ 'boosted_bkg_Fail_bins' ].add( rooDict[ 'boosted_bkg_Fail_bin'+str(ibin) ] )

        ### Pass workspace
        iCenter = templates['backgroundFail'].GetXaxis().GetBinCenter(ibin)
        #print iCenter
        #rooDict[ 'Var_Mass_bin'+str(ibin) ] = RooConstVar('Var_Mass_bin'+str(ibin), 'Var_Mass_bin'+str(ibin), iCenter)
        rooDict[ 'Var_Mass_bin'+str(ibin) ] = RooRealVar('Var_Mass_bin'+str(ibin), 'Var_Mass_bin'+str(ibin), iCenter)
        if args.pdf.startswith('poly') : rooDict[ 'poly_bin'+str(ibin) ] = RooPolyVar("poly_bin"+str(ibin),"poly_bin"+str(ibin), rooDict[ 'Var_Mass_bin'+str(ibin) ], polyArgList ) ## simple polynomial
        else:
            ptPolyString = generate_bernstein_string(args.polyDegPt)
            polyArgList.add( rooDict['Var_Mass_bin'+str(ibin)]  )
            if args.pdf.startswith('exp'): rooDict[ 'poly_bin'+str(ibin) ] = RooFormulaVar("poly_bin"+str(ibin),"poly_bin"+str(ibin), 'exp('+ptPolyString+')', polyArgList )
            else: rooDict[ 'poly_bin'+str(ibin) ] = RooFormulaVar("poly_bin"+str(ibin),"poly_bin"+str(ibin), ptPolyString, polyArgList )
            #rooDict[ 'poly_bin'+str(ibin) ] = RooChebychev("poly_bin"+str(ibin),"poly_bin"+str(ibin), rooDict[ 'Var_Mass_bin'+str(ibin) ], polyArgList ) ##does not work
            #rooDict[ 'poly_bin'+str(ibin) ] = RooBernstein("poly_bin"+str(ibin),"poly_bin"+str(ibin), rooDict[ 'Var_Mass_bin'+str(ibin) ], polyArgList ) ##does not work
            rooDict[ 'poly_bin'+str(ibin) ].Print()

        rooDict[ 'boosted_bkg_Pass_bin'+str(ibin) ] = RooFormulaVar( 'boosted_bkg_Pass_bin'+str(ibin), 'boosted_bkg_Pass_bin'+str(ibin), "@0*max(@1,0)", RooArgList( rooDict[ 'boosted_bkg_Fail_bin'+str(ibin) ], rooDict[ 'poly_bin'+str(ibin) ] ) )
        #rooDict[ 'boosted_bkg_Pass_bin'+str(ibin) ] = RooFormulaVar( 'boosted_bkg_Pass_bin'+str(ibin), 'boosted_bkg_Pass_bin'+str(ibin), "@0*max(@1,0)*@2", RooArgList( rooDict[ 'boosted_bkg_Fail_bin'+str(ibin) ], rooDict[ 'poly_bin'+str(ibin) ], bkgEFF ) )  #### NOT SURE IF bkgEFF does something
        rooDict[ 'boosted_bkg_Pass_bin'+str(ibin) ].Print()
        rooDict[ 'boosted_bkg_Pass_bins' ].add( rooDict[ 'boosted_bkg_Pass_bin'+str(ibin) ] )

    rooDict[ 'boosted_bkg_Pass_bins' ].Print()
    rooDict[ 'boosted_bkg_Pass' ] = RooParametricHist( 'boosted_bkg_Pass', 'boosted_bkg_Pass', msd, rooDict[ 'boosted_bkg_Pass_bins' ], templates['backgroundPass'] )   ### last option is only need to initialize bins
    rooDict[ 'boosted_bkg_Pass' ].Print()
    rooDict[ 'boosted_bkg_Pass_norm' ] = RooAddition( 'boosted_bkg_Pass_norm', 'boosted_bkg_Pass_norm', rooDict[ 'boosted_bkg_Pass_bins' ] )
    rooDict[ 'boosted_bkg_Pass_bins' ].Print()
    rooDict[ 'boosted_bkg_Fail' ] = RooParametricHist( 'boosted_bkg_Fail', 'boosted_bkg_Fail', msd, rooDict[ 'boosted_bkg_Fail_bins' ], templates['backgroundFail'])
    rooDict[ 'boosted_bkg_Fail' ].Print()
    rooDict[ 'boosted_bkg_Fail_norm' ] = RooAddition( 'boosted_bkg_Fail_norm', 'boosted_bkg_Fail_norm', rooDict[ 'boosted_bkg_Fail_bins' ] )

    WS_Fail = RooWorkspace("WS_Fail")
    getattr(WS_Fail, 'import')(rooDict['boosted_bkg_Fail'] ) #, RooFit.RecycleConflictNodes() ) #, RooFit.RenameAllVariablesExcept('_', 'msd'))
    getattr(WS_Fail, 'import')(rooDict['boosted_bkg_Fail_norm'], RooFit.RecycleConflictNodes() ) #, RooFit.RenameAllVariablesExcept('_', 'msd'))
    WS_Fail.Print('V')
    WS_Fail.writeToFile( combineFolder+'/alphabetWS.root' )
    WS_Pass = RooWorkspace("WS_Pass")
    getattr(WS_Pass, 'import')(rooDict['boosted_bkg_Pass'] ) #, RooFit.RecycleConflictNodes() ) #, RooFit.RenameAllVariablesExcept('_', 'msd'))
    getattr(WS_Pass, 'import')(rooDict['boosted_bkg_Pass_norm'], RooFit.RecycleConflictNodes() ) #, RooFit.RenameAllVariablesExcept('_', 'msd'))
    WS_Pass.Print('V')
    WS_Pass.writeToFile( combineFolder+'/alphabetWS.root', False )

    ### data signal workspace
    WS_data = RooWorkspace("WS_data")
    for X in [ 'Pass', 'Fail' ]:
        rooDict[ 'data_obs_'+X ] = RooDataHist("data_obs_"+X,"data_obs_"+X, RooArgList( msd ), templates['background'+X] )
        getattr(WS_data,'import')(rooDict[ 'data_obs_'+X ] )
        rooDict[ 'TTH_PTH_GT300_'+X ] = RooDataHist("TTH_PTH_GT300_"+X,"TTH_PTH_GT300_"+X, RooArgList( msd ), templates['ttH'+X] )
        getattr(WS_data,'import')(rooDict[ 'TTH_PTH_GT300_'+X ] )
        for iunc in uncList:
            for UpDown in  ['Up', 'Down']:
                rooDict[ 'TTH_PTH_GT300_'+iunc+UpDown+'_'+X ] = RooDataHist('TTH_PTH_GT300_'+X+'_CMS_ttHbb_'+iunc+UpDown,'TTH_PTH_GT300_'+X+'_CMS_ttHbb_'+iunc+UpDown, RooArgList( msd ), templates['CMS_ttHbb_'+iunc+UpDown+X] )
                getattr(WS_data,'import')(rooDict[ 'TTH_PTH_GT300_'+iunc+UpDown+'_'+X ] )
    WS_data.writeToFile( combineFolder+'/alphabetWS_base.root', False )
    WS_data.Print()


    #### Combine data card
    combineLabel='_r'+str(args.rMin)+'to'+str(args.rMax)
    datacardLabel='ttHbb'+combineLabel+'.txt'
    datacard = open( combineFolder+'/'+datacardLabel, 'w')
    datacard.write("imax 2 number of bins \n")
    datacard.write("jmax * number of processes minus 1 \n")
    datacard.write("kmax * number of nuisance parameters \n")
    datacard.write("-------------------------------\n")
    datacard.write("shapes data_obs           *     alphabetWS_base.root WS_data:$PROCESS_$CHANNEL\n")
    datacard.write("shapes TTH_PTH_GT300      *     alphabetWS_base.root WS_data:$PROCESS_$CHANNEL WS_data:$PROCESS_$CHANNEL_$SYSTEMATIC\n")
    datacard.write("shapes boosted_bkg        Fail     alphabetWS.root WS_Fail:$PROCESS_Fail\n")
    datacard.write("shapes boosted_bkg        Pass     alphabetWS.root WS_Pass:$PROCESS_Pass\n")
    datacard.write("-------------------------------\n")
    datacard.write("bin           Pass  Fail\n")
    datacard.write("observation   -1    -1\n")
    datacard.write("-------------------------------\n")
    datacard.write("bin           Fail              Fail                Pass                Pass\n")
    datacard.write("process       boosted_bkg       TTH_PTH_GT300       boosted_bkg         TTH_PTH_GT300\n")
    datacard.write("process       1                 0                   1                   0\n")
    datacard.write('rate          1                 -1                  1                   -1\n')
    #datacard.write('rate          1           '+str(round(templates['ttHFail'].Integral(),2))+'        1         '+str(round(templates['ttHPass'].Integral(),2))+'\n')   #### similar as previous line. Check if needed
    datacard.write("-------------------------------\n")
    datacard.write("lumi_13TeV_2017    lnN     -        1.023     -        1.023\n")
    datacard.write("bkgeff flatParam\n")
    for iunc in uncList: datacard.write("CMS_ttHbb_"+iunc+"     shape   -        1       -     1\n")
    for q, k in rooDict.iteritems():
        if q.startswith('boosted_bkg_paramX'): datacard.write(q+"    flatParam\n")
        if q.startswith('boosted_bkg_Fail_bin') and not q.endswith(('func', 'In', 'unc', 's')):
            datacard.write(q+"    flatParam\n")
    datacard.close()
    ##############################

    combineCmd = 'combine -M FitDiagnostics %s -n %s --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --saveShapes --saveWorkspace --setParameterRanges r=%i,%i'%(datacardLabel,combineLabel,args.rMin,args.rMax)
    #if not isData: combineCmd += '--plot '
    combineCmd += ' --plot'
    exec_me(combineCmd, folder=combineFolder)

    sys.exit(0)

#####    Manual rhalphabet DOES NOT WORK YET
#####    ### Loop over pt bins
#####    print( '-'*20, "MYINFO - number of pt bins: ", len(rhalPtList))
#####    for k, ptbin in enumerate(rhalPtList):
#####        if ptbin.endswith(rhalPtList[-1]): continue  ### last
#####        print( '-'*20, "MYINFO - pt bin number: ", k)
#####
#####        ### Converting TH1 into list
#####        #hbins_inPtBin = []
#####        #hPass_inPtBin = []
#####        #hFail_inPtBin = []
#####        #for ibin in range(bkgHistos['MassPtPass'+ptbin].GetNbinsX()):
#####        #    hbins_inPtBin.append( bkgHistos['MassPtPass'+ptbin].GetBinLowEdge(ibin) )
#####        #    hPass_inPtBin.append( bkgHistos['MassPtPass'+ptbin].GetBinContent(ibin) )
#####        #    hFail_inPtBin.append( bkgHistos['MassPtFail'+ptbin].GetBinContent(ibin) )
#####
#####        #################################################################
#####        ### Make RooDataset, RooPdfs, and histograms
#####        ### In principle this can be extended for diff bkgs
#####        cats = RooCategory( "sample", "sample" )
#####        cats.defineType("Pass",1)
#####        cats.defineType("Fail",0)
#####        rooDict[ 'PassBkg'+ptbin ] = RooDataHist("bkg_Pass_"+ptbin,"bkg_Pass_"+ptbin, RooArgList( MSD ), bkgHistos['MassPtPass'+ptbin] )
#####        rooDict[ 'FailBkg'+ptbin ] = RooDataHist("bkg_Fail_"+ptbin,"bkg_Fail_"+ptbin, RooArgList( MSD ), bkgHistos['MassPtFail'+ptbin] )
#####        rooDict[ 'Bkg'+ptbin ] = RooDataHist("comb_bkg_"+ptbin,"comb_bkg_"+ptbin, RooArgList( MSD ), RooFit.Index(cats), RooFit.Import( "Pass", rooDict['PassBkg'+ptbin] ), RooFit.Import( "Fail", rooDict['FailBkg'+ptbin] ) )
#####
#####        ### Normalization: RooExtendPdfs are coupled via their normalizations, N*eff or N*(1-eff).
#####        totalN = bkgHistos['MassPtPass'+ptbin].Integral()+bkgHistos['MassPtFail'+ptbin].Integral()
#####        rooDict[ 'bkgMC_total_norm_'+ptbin ] = RooRealVar( "bkgMC_norm"+ptbin, "bkgMC_norm"+ptbin, totalN, 0., 5*totalN )
#####        rooDict[ 'bkgMC_Pass_norm_'+ptbin ] = RooFormulaVar( "bkgMC_fPass"+ptbin, "bkgMC_norm"+ptbin+"*(veff)", RooArgList( rooDict['bkgMC_total_norm_'+ptbin], EFF ) )
#####        rooDict[ 'bkgMC_Fail_norm_'+ptbin ] = RooFormulaVar( "bkgMC_fqail"+ptbin, "bkgMC_norm"+ptbin+"*(1-veff)", RooArgList( rooDict['bkgMC_total_norm_'+ptbin], EFF ) )
#####
#####        ### Shapes
#####        rooDict[ 'bkgMC_Pass_'+ptbin ] = RooDataHist("bkgMC_Pass_"+ptbin,"bkgMC_Pass_"+ptbin, RooArgList( MSD ), bkgHistos['MassPtPass'+ptbin] )
#####        rooDict[ 'bkgMC_Fail_'+ptbin ] = RooDataHist("bkgMC_Fail_"+ptbin,"bkgMC_Fail_"+ptbin, RooArgList( MSD ), bkgHistos['MassPtFail'+ptbin] )
#####        rooDict[ 'bkgMC_Passh_'+ptbin ] = RooHistPdf( 'bkgMC_Passh_'+ptbin, 'bkgMC_Passh_'+ptbin, RooArgList(SHIFT), RooArgList(MSD), rooDict['bkgMC_Pass_'+ptbin], 0 )
#####        rooDict[ 'bkgMC_Failh_'+ptbin ] = RooHistPdf( 'bkgMC_Failh_'+ptbin, 'bkgMC_Failh_'+ptbin, RooArgList(SHIFT), RooArgList(MSD), rooDict['bkgMC_Fail_'+ptbin], 0 )
#####
#####        ### extended likelihood from normalization and shape above
#####        rooDict[ 'bkgMC_Passe_'+ptbin ] = RooExtendPdf( 'bkgMC_Passe_'+ptbin, 'bkgMC_Passe_'+ptbin, rooDict['bkgMC_Passh_'+ptbin], rooDict['bkgMC_Pass_norm_'+ptbin] )
#####        rooDict[ 'bkgMC_Passe_'+ptbin ].Print()
#####        rooDict[ 'bkgMC_Faile_'+ptbin ] = RooExtendPdf( 'bkgMC_Faile_'+ptbin, 'bkgMC_Faile_'+ptbin, rooDict['bkgMC_Failh_'+ptbin], rooDict['bkgMC_Fail_norm_'+ptbin] )
#####        rooDict[ 'bkgMC_Faile_'+ptbin ].Print()
#####
#####
#####        ### Add all bkg in RooAddpdf
#####        ### needed if differnt bkg components
#####        rooDict[ 'totalMC_pdf_Pass'+ptbin ] = RooAddPdf( 'totalMC_Pass'+ptbin, 'totalMC_Pass'+ptbin, RooArgList( rooDict[ 'bkgMC_Passe_'+ptbin ] ) )
#####        rooDict[ 'totalMC_pdf_Fail'+ptbin ] = RooAddPdf( 'totalMC_Fail'+ptbin, 'totalMC_Fail'+ptbin, RooArgList( rooDict[ 'bkgMC_Faile_'+ptbin ] ) )
#####
#####        ### Make RooSimultaneous
#####        rooDict[ 'total_simulpdf'+ptbin ] = RooSimultaneous( 'tot'+ptbin, 'tot'+ptbin, cats )
#####        rooDict[ 'total_simulpdf'+ptbin ].addPdf( rooDict[ 'totalMC_pdf_Pass'+ptbin ], 'Pass' )
#####        rooDict[ 'total_simulpdf'+ptbin ].addPdf( rooDict[ 'totalMC_pdf_Fail'+ptbin ], 'Fail' )
#####        rooDict[ 'total_simulpdf'+ptbin ].Print()
#####
#####
#####        #################################################################
#####        ### Make Rhalphabet
#####        meanPtBin = (float(rhalPtList[k+1]) + float(ptbin))/2
#####        PT.setVal(meanPtBin)
#####        print( '-'*20, "MYINFO - this pt bin value: ", meanPtBin)
#####
#####        ### Build polynomial
#####        maxPolyPt = 0
#####        maxPolyRho = 0
#####        rooDict[ 'polyArray'+ptbin ] = buildPolynomialArray( 'p', 'r', maxPolyPt, maxPolyRho, -30, 30 )
#####        print( '-'*20, "MYINFO - polynomial_variables: ", rooDict[ 'polyArray'+ptbin ])
#####
#####        ### Now build the function
#####        lUnity = RooConstVar("unity","unity",1.)
#####        lZero  = RooConstVar("lZero","lZero",0.)
#####
#####        for massbin in range(bkgHistos['MassPtPass'+ptbin].GetNbinsX()):
#####            if bkgHistos['MassPtPass'+ptbin].GetXaxis().GetBinLowEdge(massbin) < 50: continue  ### skipping low mass
#####            MSD.setVal( bkgHistos['MassPtPass'+ptbin].GetXaxis().GetBinCenter(massbin) )
#####
#####            #rooDict[ 'lPass'+ptbin+str(massbin) ] = buildRooPolyArray( PT.getVal(), maxPolyPt, MSD.getVal(), maxPolyRho, lUnity, rooDict[ 'polyArray'+ptbin ])
#####            rooDict[ 'roopolyarray'+ptbin+str(massbin) ] = buildRooPolyRhoArrayBernsteinExp( PT.getVal(), maxPolyPt, MSD.getVal(), maxPolyRho, lUnity, lZero, rooDict[ 'polyArray'+ptbin ])
#####            #rooDict[ 'roopolyarray'+ptbin+str(massbin) ].Print()
#####            #sys.exit(0)
#####
#####            Fail_bin_content = bkgHistos[ 'MassPtFail'+ptbin ].GetBinContent(massbin)
#####            #print Fail_bin_content, bkgHistos[ 'MassPtFail'+ptbin ].GetBinCenter(massbin)
#####
#####        Pass_workspace = RooWorkspace('w_Pass_cat'+ptbin)
#####        Fail_workspace = RooWorkspace('w_Fail_cat'+ptbin)
#####        #getattr(Pass_workspace, 'import')(Pass_rparh, RooFit.RecycleConflictNodes(), r.RooFit.RenameAllVariablesExcept(self._suffix.replace('_',''),'x'))
#####        #getattr(Pass_workspace, 'import')(Pass_norm, r.RooFit.RecycleConflictNodes(), r.RooFit.RenameAllVariablesExcept(self._suffix.replace('_',''),'x'))
#####
#####    print '-'*10
#####    print 'To make it run: (dont forget to set the combine environment)'
#####    print "combine -M FitDiagnostics datacard.txt  --robustFit 1 --setRobustFitAlgo Minuit2,Migrad --saveNormalizations --plot --saveShapes --saveWorkspace"


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='1D', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
    parser.add_argument('-v', '--version', action='store', default='v14', help='Version: v01, v02.' )
    #parser.add_argument('-d', '--degreePoly', action='store', default='4', help='Degree of polynominal.' )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-s', '--selection', nargs='+', default=['met20_btagDDBvL_noMD07','met20_deepTagMD_bbvsLight05845','met20_deepTagMD_bbvsLight08695'], help='event selection, in file paths')
    parser.add_argument('-j', '--jsonpath', default='/eos/home-a/algomez/tmpFiles/hepacc/results/', help='path to json files')
    parser.add_argument('-o','--outdir', default=None, help='specifiy a custom output directory')
    parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')
    parser.add_argument('-y', '--year', default='2017', type=str, help='year to process, in file paths')
    parser.add_argument('--msd_start', default=90, type=int, help='start of the mass range')
    parser.add_argument('--msd_stop', default=170, type=int, help='stop of the mass range')
    parser.add_argument('--polyDegPt', default=2, type=int, help='degree of polynomial to fit pt')
    parser.add_argument('--polyDegRho', default=2, type=int, help='degree of polynomial to fit rho')
    parser.add_argument('-r', '--rebin_factor', default=5, type=int, help='rebin factor for json histograms, default mass bin size is 1GeV')
    parser.add_argument('--nptbins', default=1, type=int, help='number of pt bins')
    parser.add_argument('--sig-and-bkg', action='store_true', default=False, help='sum signal and background samples when running on MC')
    parser.add_argument('--pdf', default='Cheb', choices=['poly','exp', 'Cheb', 'Bern'])
    parser.add_argument('--rMin', default=-20, type=float, help='minimum of r (signal strength) in combine fit')
    parser.add_argument('--rMax', default=20, type=float, help='maximum of r (signal strength) in combine fit')

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if not os.path.exists('Plots/'): os.makedirs('Plots/')
    CMS_lumi.extraText = "Preliminary Simulation"
    #CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, 2017"
    CMS_lumi.lumi_13TeV = "13 TeV, 2017"

    if args.nptbins==1:
        ptbins = np.array([300,5000])
    elif args.nptbins==2:
        ptbins = np.array([250,300,5000])
    elif args.nptbins==3:
        ptbins = np.array([250,300,450,5000])
    else:
        raise Exception('invalid number of ptbins')

    msd_start    = args.msd_start
    msd_stop     = args.msd_stop
    polyDegRho   = args.polyDegRho
    rebin_factor = args.rebin_factor
    uncList = [ 'puWeight' ] # 'AK4deepjetM', 'AK8DDBvLM1', 'jer', 'jesAbsolute', 'jesAbsolute_'+args.year, 'jesBBEC1', 'jesBBEC1_'+args.year, 'jesEC2', 'jesEC2_'+args.year, 'jesFlavorQCD', 'jesHF', 'jesRelativeBal', 'jesRelativeSample_'+args.year, 'jmr', 'jms', 'pdfWeight', 'psWeight_FSR', 'psWeight_ISR' ]
    for iunc in uncList:
      if iunc.endswith('allyears'):
          j = uncList.index(iunc)
          del uncList[j]
          tmp = iunc.split("_")[0]
          uncList+=[ tmp+'_2016', tmp+'_2017', tmp+'_2018' ]
    print(uncList)

    for sel in args.selection:
        indir = os.path.join(args.jsonpath, args.year, args.version, sel)
        if not os.path.exists(indir):
            raise Exception('invalid input path: %s'%indir)
        if args.outdir is None: outdir = os.path.join('output', args.year, args.version, sel)
        else: outdir = args.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    #massDecorrelation( 'leadAK8JetRhoPtHbb_2JdeltaR2WTau21DDT_boostedHiggs' )  ### this was just a test, but maybe needed later
    #for i in plotList:
    BkgEstimation( ) #i[0], i[2], i[3], i[4], axisX=i[1])

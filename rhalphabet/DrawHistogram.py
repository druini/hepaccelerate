#!/usr/bin/env python
'''
Example: for a histogram named: leadAK8JetMass_2J2WdeltaR_overlap_weights_nominaldata_nominal
./DrawHistogram.py -p signalBkgData -v met20_btagDDBvL086 -y 2018 -s leadAK8JetMass -c 2J2WdeltaR_overlap_weights_nominal -j -L
'''

from ROOT import *
import json, glob
import time, os, math, sys, copy
from array import array
import argparse
from collections import OrderedDict
import subprocess
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle

####gReset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)

xline = array('d', [0,2000])
yline = array('d', [1,1])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

canvas = {}

def rootHistograms( version, lumi, year):
    """Information about the samples included in the plots.
    If running from root files, it will open the root histogram described below.
    Contains also the color of the histrogram and the label.
    If running from json files, it only loops over the processes.
    """

    dataFiles = OrderedDict()
    bkgFiles = OrderedDict()
    signalFiles = OrderedDict()
    extra='_boosted_'+year

    bkgFiles["THW"] = [ TFile('Rootfiles/'+version+'/histograms_THW'+extra+'.root'), 46, 'tHW' ]
    bkgFiles["TTGJets"] = [ TFile('Rootfiles/'+version+'/histograms_TTGJets'+extra+'.root'), kGreen+3, 'tt+V/Gluon' ]
    bkgFiles["TTWJetsToQQ"] = [ TFile('Rootfiles/'+version+'/histograms_TTWJetsToQQ'+extra+'.root'), kGreen+3, 'tt+V/Gluon' ]
    bkgFiles["TTWJetsToLNu"] = [ TFile('Rootfiles/'+version+'/histograms_TTWJetsToLNu'+extra+'.root'), kGreen+3, 'tt+V/Gluon' ]
    bkgFiles["TTZToQQ"] = [ TFile('Rootfiles/'+version+'/histograms_TTZToQQ'+extra+'.root'),  kGreen+3, 'tt+V/Gluon' ]
    bkgFiles["TTZToLLNuNu"] = [ TFile('Rootfiles/'+version+'/histograms_TTZToLLNuNu'+extra+'.root'),  kGreen+3, 'tt+V/Gluon' ]
    bkgFiles["WW"] = [ TFile('Rootfiles/'+version+'/histograms_WW'+extra+'.root'), kYellow+2, 'Dibosons' ]
    bkgFiles["WZ"] = [ TFile('Rootfiles/'+version+'/histograms_WZ'+extra+'.root'), kYellow+2, 'Dibosons' ]
    bkgFiles["ZZ"] = [ TFile('Rootfiles/'+version+'/histograms_ZZ'+extra+'.root'), kYellow+2, 'Dibosons' ]
    bkgFiles["QCD_HT500to700"] = [ TFile('Rootfiles/'+version+'/histograms_QCD_HT500to700'+extra+'.root'), kMagenta+3 , 'QCD' ]
    bkgFiles["QCD_HT700to1000"] = [ TFile('Rootfiles/'+version+'/histograms_QCD_HT700to1000'+extra+'.root'), kMagenta+3 , 'QCD' ]
    bkgFiles["QCD_HT1000to1500"] = [ TFile('Rootfiles/'+version+'/histograms_QCD_HT1000to1500'+extra+'.root'), kMagenta+3 , 'QCD' ]
    bkgFiles["QCD_HT1500to2000"] = [ TFile('Rootfiles/'+version+'/histograms_QCD_HT1500to2000'+extra+'.root'), kMagenta+3 , 'QCD' ]
    bkgFiles["QCD_HT2000toInf"] = [ TFile('Rootfiles/'+version+'/histograms_QCD_HT2000toInf'+extra+'.root'), kMagenta+3 , 'QCD' ]
    bkgFiles["TTTo2L2Nu"] = [ TFile('Rootfiles/'+version+'/histograms_TTTo2L2Nu'+extra+'.root'), kBlack, 'Non-SL t#bar{t}' ]
    bkgFiles["TTToHadronic"] = [ TFile('Rootfiles/'+version+'/histograms_TTToHadronic'+extra+'.root'), kBlack, 'Non-SL t#bar{t}' ]
    bkgFiles["ST_s-channel_4f_leptonDecays"] = [ TFile('Rootfiles/'+version+'/histograms_ST_s-channel_4f_leptonDecays'+extra+'.root'),  40, 'Single top' ]
    bkgFiles["ST_t-channel_top"] = [ TFile('Rootfiles/'+version+'/histograms_ST_t-channel_top'+extra+'.root'),  kCyan+1, 'Single top' ]
    bkgFiles["ST_t-channel_antitop"] = [ TFile('Rootfiles/'+version+'/histograms_ST_t-channel_antitop'+extra+'.root'),  kCyan+1, 'Single antitop' ]
    bkgFiles["ST_tW_antitop"] = [ TFile('Rootfiles/'+version+'/histograms_ST_tW_antitop'+extra+'.root'), kCyan+1, 'Single top' ]
    bkgFiles["ST_tW_top"] = [ TFile('Rootfiles/'+version+'/histograms_ST_tW_top'+extra+'.root'), kCyan+1, 'Single top' ]
    bkgFiles["DYJetsToLL"] = [ TFile('Rootfiles/'+version+'/histograms_DYJetsToLL'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["ZJetsToQQ_HT400to600"] = [ TFile('Rootfiles/'+version+'/histograms_ZJetsToQQ_HT400to600'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["ZJetsToQQ_HT600to800"] = [ TFile('Rootfiles/'+version+'/histograms_ZJetsToQQ_HT600to800'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["ZJetsToQQ_HT-800toInf"] = [ TFile('Rootfiles/'+version+'/histograms_ZJetsToQQ_HT-800toInf'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-200To400"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-200To400'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-400To600"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-400To600'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-600To800"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-600To800'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-800To1200"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-800To1200'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-1200To2500"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-1200To2500'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToLNu_HT-2500ToInf"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToLNu_HT-2500ToInf'+extra+'.root'), kAzure+2, '(W/Z/DY)+Jets' ]
    ##bkgFiles["WJetsToQQ_HT400to600"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToQQ_HT400to600'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToQQ_HT600to800"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToQQ_HT600to800'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["WJetsToQQ_HT-800toInf"] = [ TFile('Rootfiles/'+version+'/histograms_WJetsToQQ_HT-800toInf'+extra+'.root'),  kAzure+2, '(W/Z/DY)+Jets' ]
    bkgFiles["TTToSemiLeptonic"] = [ TFile('Rootfiles/'+version+'/histograms_TTToSemiLeptonic'+extra+'.root'), kWhite, 'SL t#bar{t}' ]
    #bkgFiles[""] = [ TFile('Rootfiles/'+version+'/'), 1 ]

    signalFiles["ttHTobb"] = [ TFile('Rootfiles/'+version+'/histograms_ttHTobb'+extra+'.root'), kRed, 'ttH' ]
    signalFiles["ttHToNonbb"] = [ TFile('Rootfiles/'+version+'/histograms_ttHToNonbb_M125'+extra+'.root'), kRed, 'ttH' ]
    #signalFiles[""] = [ TFile('Rootfiles/'+version+'/'), 1 ]

    dataFiles['SingleElectron'] = TFile.Open('Rootfiles/'+version+'/histograms_SingleElectron_Run'+year+'ALL'+extra+'.root')
    dataFiles['SingleMuon'] = TFile.Open('Rootfiles/'+version+'/histograms_SingleMuon_Run'+year+'ALL'+extra+'.root')

    return bkgFiles, signalFiles, dataFiles

##########################################################
def jsonToTH1( jsonFile, variables, debug=False ):
    """Open json files and convert those histograms into TH1 root"""

    ## opening the json file
    with open(jsonFile) as json_file:
        data = json.load(json_file)

    ## priting the list of histograms in json
    if debug:
        print("In jsonFile: ", jsonFile, "the histograms found are: ")
        for i in data.keys(): print(i)

    histoDict = ''
    if 'Single' in jsonFile: isData = True
    else: isData = False

    ## creating histograms with the information in json
    for xvar in data:
        for jvar in variables:
            if xvar.endswith(jvar):
                histoDict = TH1F( xvar+jsonFile.split('out_')[1], xvar, len(data[xvar]["edges"])-1, data[xvar]["edges"][0], data[xvar]["edges"][-1] )
                histoDict.Sumw2()
                for icont in range(len(data[xvar]["contents"])):
                    if data[xvar]["contents"][icont]<0: continue
                    histoDict.SetBinContent( icont+1, data[xvar]["contents"][icont] )
                    histoDict.SetBinError( icont+1, TMath.Sqrt(data[xvar]["contents"][icont]) )

    return histoDict


##########################################################
def setSelection( listSel, xMin=0.65, yMax=0.65, align='right' ):
    """Not used yet. It will add a label with the selection"""

    for i in range( len( listSel ) ):
        textBox=TLatex()
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        if 'right' in align: textBox.SetTextAlign(31)
        textBox.SetTextFont(62) ### 62 is bold, 42 is normal
        textBox.DrawLatex(xMin, yMax, listSel[i])
        yMax = yMax -0.05


##########################################################
def plotQuality( nameInRoot, label, xmin, xmax, rebinX, labX, labY ):
    """It only checks shapes, i.e. all data vs all bkgs"""

    outputFileName = nameInRoot+'_dataQualityPlots_'+args.year+'_'+args.version+'.'+args.ext
    print 'Processing.......', outputFileName

    histos = {}
    for idataLabel in dataFiles:
        if args.json:
            #for iSamData in glob.glob(folder+'/*'+idataLabel+'*'):
            #    histos[ iSamData.split('out_')[1].split('.json')[0] ] = jsonToTH1( iSamData, [nameInRoot] )
            histos[ 'tmp' ] = jsonToTH1( folder+'out_data_nominal_merged.json', [nameInRoot] )
        else: histos[ idataLabel ] = dataFiles[idataLabel].Get( 'tthbb13/'+nameInRoot )
    for ihdata in histos.keys():
        try: histos[ 'AllData' ].Add( histos[ ihdata ].Clone() )
        except (KeyError, AttributeError) as e:
            histos[ 'AllData' ] = histos[ ihdata ].Clone()
    print histos['AllData'].Integral()

    histos[ 'Bkg' ] = histos[ 'AllData' ].Clone()
    histos[ 'Bkg' ].Reset()
    for isamLabel in bkgFiles:
        if args.json:
            try: histos[ isamLabel ] = jsonToTH1( folder+'/out_'+isamLabel+'_nominal_merged.json', [nameInRoot] )
            except IOError:
                print 'Sample missing: ', isamLabel
                bkgFiles.pop( isamLabel )
                continue
        else:
            histos[ isamLabel ] = bkgFiles[ isamLabel ][0].Get( 'tthbb13/'+nameInRoot )
            try: tmpScale = args.lumi*checkDict( isamLabel, dictSamples )['XS']/checkDict( isamLabel, dictSamples )[args.year][1]
            except KeyError:
                tmpScale = 1
                print 'Sample missing :', isamLabel
            histos[ isamLabel ].Scale( tmpScale )
        if isinstance(histos[ isamLabel ], TH1): histos[ 'Bkg' ].Add( histos[ isamLabel ] )
        else: print 'Sample missing: ', isamLabel

    if rebinX != 1:
        histos[ 'AllData' ].Rebin( rebinX )
        histos[ 'Bkg' ].Rebin( rebinX )
    hData = histos[ 'AllData' ].Clone()
    hBkg = histos[ 'Bkg' ].Clone()

    hRatio = TGraphAsymmErrors()
    hRatio.Divide( hData, hBkg, 'pois' )
    hRatioStatErr = hBkg.Clone()
    hRatioStatErr.Divide( hBkg )
    hRatioStatErr.SetFillColor(kBlack)
    hRatioStatErr.SetFillStyle(3004)

    binWidth = histos['AllData'].GetBinWidth(1)

    if (labY < 0.5) and ( labX < 0.5 ): legend=TLegend(0.20,0.50,0.50,0.62)
    elif (labX < 0.5): legend=TLegend(0.20,0.75,0.50,0.87)
    else: legend=TLegend(0.70,0.75,0.90,0.87)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.AddEntry( hData, 'DATA' , 'ep' )
    legend.AddEntry( hBkg, 'All Bkg', 'lp' )

    hBkg.SetLineColor(kRed-4)
    hBkg.SetLineWidth(2)
    #hBkg.SetFillColor(kBlack)
    hBkg.SetFillStyle(3004)
    hData.SetMarkerStyle(8)

    tdrStyle.SetPadRightMargin(0.05)
    tdrStyle.SetPadLeftMargin(0.15)
    can = TCanvas('c1', 'c1',  10, 10, 750, 750 )
    pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
    pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    if args.log: pad1.SetLogy()
    hData.Draw("E")
    hBkg.Draw('hist same E1')
    hData.Draw("same E")
    hData.SetMaximum( 1.5* max( hData.GetMaximum(), hBkg.GetMaximum() )  )
    hData.SetMinimum( 1 )
    #hData.GetYaxis().SetTitleOffset(1.2)
    if xmax: hData.GetXaxis().SetRangeUser( xmin, xmax )
    #hData.GetYaxis().SetTitle( 'Normalized' )
    #hData.GetYaxis().SetTitle( 'Normalized / '+str(int(binWidth))+' GeV' )
    hData.GetYaxis().SetTitle(  'Events / '+str(int(binWidth))+' GeV' )

    #CMS_lumi.relPosX = 0.13
    if args.final:
        CMS_lumi.cmsTextOffset = 0.1
        CMS_lumi.relPosX = 0.15
    else:
        CMS_lumi.cmsTextOffset = 0.0
        CMS_lumi.relPosX = 0.13
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    #labelAxis( name, hData, '' )
    legend.Draw()

    pad2.cd()
    gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame(xmin, ( 0.5 if args.addFit else 0.5), xmax,1.5)
    #labelAxis( name.replace( args.cut, ''), tmpPad2, ( 'softDrop' if 'Puppi' in args.grooming else Groom ) )
    tmpPad2.GetXaxis().SetTitle( label )
    tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    hRatio.SetMarkerStyle(8)
    hRatio.Draw('P')
    hRatioStatErr.Draw('same e2')
    if args.addFit:
        fitLine = TF1( 'fitLine', 'pol1', 0, 2 ) #800, 5000)
        hRatio.Fit( 'fitLine', 'MIR')
        fitLine.Draw("same")
        pad2.Update()
        st1 = hRatio.GetListOfFunctions().FindObject("stats")
        st1.SetX1NDC(.65)
        st1.SetX2NDC(.95)
        st1.SetY1NDC(.75)
        st1.SetY2NDC(.95)
        #st1.SetTextColor(kRed)
        pad2.Modified()

    can.SaveAs( 'Plots/'+ outputFileName.replace('Plots', ( 'Fit' if args.addFit else '') ) )
    del can

########################################################################
def plotSimpleComparison( inFile1, sample, inFile2, sample2, name, rebinX=1, xmin='', xmax='', labX=0.92, labY=0.50, axisX='', axisY='', log=False, ext='png', Norm=False ):
    """"Take two root files, make simple comparison plot"""

    outputFileName = name+'_'+sample+sample2+'_simpleComparisonPlot'+args.version+'.'+ext
    print('Processing.......', outputFileName)

    histo =  inFile1 if isinstance( inFile1, TH1F ) else inFile1.Get( 'tthbb13/'+name )
    print histo
    if rebinX!=1: histo.Rebin( rebinX )
    histo2 =  inFile2 if isinstance( inFile2, TH1F ) else inFile2.Get( 'tthbb13/'+name )
    if rebinX!=1: histo2.Rebin( rebinX )

    binWidth = histo.GetBinWidth(1)

    legend=TLegend(0.60,0.75,0.90,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    #histo.SetFillColor(48)
    histo.SetFillStyle(1001)

    tdrStyle.SetPadRightMargin(0.05)
    canvas[name] = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    if log:
        canvas[name].SetLogy()
        outName = outputFileName.replace('_simplePlot','_Log_simplePlot')
    else: outName = outputFileName

    legend.AddEntry( histo, sample, 'f' )
    legend.AddEntry( histo2, sample2, 'f' )
    if xmax and xmin: histo.GetXaxis().SetRangeUser( xmin, xmax )
    histo.GetYaxis().SetTitleOffset(0.90)
    histo.SetMaximum( 1.3*max( histo.GetMaximum(), histo2.GetMaximum() )  )
    histo.SetLineColor(kRed)
    histo2.SetLineColor(kBlue)
    histo.Draw('hist')
    histo2.Draw('hist same')
    if not axisY: histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: histo.GetXaxis().SetTitle( axisX )

    #labelAxis( name, histo, '' )
    legend.Draw()

    canvas[name].SaveAs( 'Plots/'+outName )
    #del can

########################################################################
def plotSignalBkg( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False, addRatioFit=False, Norm=False, ext='png' ):
    """function to plot s and b histos"""

    outputFileName = name+'_PlusBkg_AnalysisPlots_'+args.year+'_'+args.version+'.'+ext
    if args.process.endswith('Data'): outputFileName = outputFileName.replace('PlusBkg','BkgData')
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    if Norm: outputFileName = outputFileName.replace('Plots','Plots_Normalized')
    print('Processing.......', outputFileName)

    if args.process.endswith('Data'): legend=TLegend(0.69,0.48,0.90,0.88)
    else: legend=TLegend(0.60,0.60,0.90,0.90)
    if name.startswith("leadAK8JetMass"):
        legend=TLegend(0.45,0.65,0.90,0.88)
        legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    ######################## Data histos
    dataHistos = {}
    if args.process.endswith('Data'):
        for idata in dataFiles:
            if args.json:
                dataHistos[ 'tmp' ] = jsonToTH1( folder+'out_data_nominal_merged.json', [name] )
            else: dataHistos[ idata ] = dataFiles[idata].Get( 'tthbb13/'+name )
        for ihdata in dataHistos.keys():
            try: dataHistos[ 'AllData' ].Add( dataHistos[ ihdata ].Clone() )
            except (KeyError, AttributeError) as e:
                dataHistos[ 'AllData' ] = dataHistos[ ihdata ].Clone()
        if rebinX > 1: dataHistos[ "AllData" ] = dataHistos[ "AllData" ].Rebin( rebinX )
        if Norm: dataHistos[ "AllData" ].Scale( 1 /dataHistos["AllData"].Integral() )
        legend.AddEntry( dataHistos[ 'AllData' ], 'Data', 'lep' )

    ######################## Bkg histos
    bkgHistos = OrderedDict()
    binWidth = 0
    maxList = []
    bkgInMassWindow = 0
    bkgInMassWindowErr = 0
    if len(bkgFiles) > 0:
        for bkgSamples in bkgFiles:
            if args.json:
                try:
                    bkgHistos[ bkgSamples ] = jsonToTH1( folder+'/out_'+bkgSamples+'_nominal_merged.json', [name])
                except IOError:
                    print 'Sample missing: ', bkgSamples
                    bkgFiles.pop( bkgSamples )
                    continue
            else:
                bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( 'tthbb13/'+name )
                try: tmpScale = args.lumi*checkDict( bkgSamples, dictSamples )['XS']/checkDict( bkgSamples, dictSamples )[args.year][1]
                except KeyError:
                    tmpScale = 1
                    print 'Sample missing :', bkgSamples
                bkgHistos[ bkgSamples ].Scale( tmpScale )
            bkgHistos[ bkgSamples ].SetTitle(bkgSamples)
            #print(bkgSamples, round(bkgHistos[ bkgSamples ].Integral(), 2) )
            if rebinX > 1: bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( rebinX )

            if Norm:
                bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][1] )
                bkgHistos[ bkgSamples ].SetLineWidth( 2 )
                try: bkgHistos[ bkgSamples ].Scale( 1 / bkgHistos[ bkgSamples ].Integral() )
                except ZeroDivisionError: pass
                maxList.append( bkgHistos[ bkgSamples ].GetMaximum() )
            else:
                bkgHistos[ bkgSamples ].SetFillStyle( 1001 )
                bkgHistos[ bkgSamples ].SetFillColor( int(bkgFiles[ bkgSamples ][1]) )

    #### Merging samples
    for bkg in bkgFiles:
        if bkg.endswith(('WZ','ZZ')):
            bkgHistos['WW'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        elif bkg.startswith('TTToHad'):
            bkgHistos['TTTo2L2Nu'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        elif bkg.startswith('ST_') and not bkg.endswith('tW_antitop'):
            bkgHistos['ST_tW_antitop'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        elif bkg.startswith(('WJets','ZJets', 'DY')) and not bkg.endswith('WJetsToLNu_HT-200To400'):
            bkgHistos['WJetsToLNu_HT-200To400'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        elif bkg.startswith('QCD_HT') and not bkg.endswith('500to700'):
            bkgHistos['QCD_HT500to700'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        elif bkg.startswith( ('TTZ', 'TTWJets', 'TTGJets') ) and not bkg.endswith('TTZToQQ'):
            bkgHistos['TTZToQQ'].Add( bkgHistos[bkg] )
            bkgHistos.pop(bkg, None)
        else:
            legend.AddEntry( bkgHistos[ bkg ], bkgFiles[bkg][2], 'l' if Norm else 'f' )

    hBkg = bkgHistos[next(iter(bkgHistos))].Clone()
    hBkg.Reset()
    for bkgSamples in bkgHistos:
        print bkgSamples, round(bkgHistos[ bkgSamples ].Integral(), 2)

    ######################## Signal histos
    signalHistos = OrderedDict()
    if len(signalFiles) > 0:
        dummySig=0
        for sigSamples in signalFiles:
            if args.json:
                try: signalHistos[ sigSamples ] = jsonToTH1( folder+'/out_'+sigSamples+'_nominal_merged.json', [name] )
                except IOError:
                    print 'Sample missing: ', sigSamples
                    signalFiles.pop( sigSamples )
                    continue
            else:
                signalHistos[ sigSamples ] = signalFiles[ sigSamples ][0].Get( 'tthbb13/'+name )
                try: tmpScale = args.lumi*checkDict( sigSamples, dictSamples )['XS']/checkDict( sigSamples, dictSamples )[args.year][1]
                except KeyError:
                    tmpScale = 1
                    print 'Sample missing :', sigSamples
                signalHistos[ sigSamples ].Scale( tmpScale )
            print(sigSamples, round(signalHistos[ sigSamples ].Integral(), 2) )
            if rebinX > 1: signalHistos[ sigSamples ] = signalHistos[ sigSamples ].Rebin( rebinX )
            if Norm:
                signalHistos[ sigSamples ].SetLineColor( signalFiles[ sigSamples ][1] )
                signalHistos[ sigSamples ].SetLineWidth( 3 )
                signalHistos[ sigSamples ].SetLineStyle( 10-dummySig )
                signalHistos[ sigSamples ].Scale( 1 / signalHistos[ sigSamples ].Integral() )
                maxList.append( signalHistos[ sigSamples ].GetMaximum() )
            else:
                signalHistos[ sigSamples ].SetLineColorAlpha( signalFiles[ sigSamples ][1], 1 )
                signalHistos[ sigSamples ].SetLineWidth(3)
                signalHistos[ sigSamples ].SetLineStyle(2+dummySig)
            binWidth = int(signalHistos[ sigSamples ].GetBinWidth( 1 ))
            dummySig+=8

    #### Merging samples
    for signal in signalFiles:
        if signal.endswith('Nonbb'):
            signalHistos['ttHTobb'].Add( signalHistos[signal] )
            signalHistos.pop(signal, None)
        else:
            legend.AddEntry( signalHistos[ signal ], signalFiles[signal][2], 'l' if Norm else 'f' )

    ################### Start ploting
    if not Norm:

        stackHisto = THStack('stackHisto'+name, 'stack'+name)
        for samples in bkgHistos:
            stackHisto.Add( bkgHistos[ samples ].Clone() )
            hBkg.Add( bkgHistos[ samples ].Clone() )
        stackSigHisto = THStack('stackSigHisto'+name, 'stackSigHisto'+name)
        for samples in signalHistos:
            stackSigHisto.Add( signalHistos[ samples ].Clone() )

        canvas[outputFileName] = TCanvas('c1'+name, 'c1'+name,  10, 10, 750, (750 if args.process.endswith('Data') else 500 ) )
        if args.process.endswith('Data'):
            tdrStyle.SetPadRightMargin(0.05)
            tdrStyle.SetPadLeftMargin(0.15)
            pad1 = TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
            pad2 = TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
            pad1.Draw()
            pad2.Draw()

            pad1.cd()
            if log: pad1.SetLogy()
        elif log: canvas[outputFileName].SetLogy()
        stackHisto.Draw('hist')
        stackSigHisto.Draw('hist same')

        if xmax: stackHisto.GetXaxis().SetRangeUser( xmin, xmax )
        stackHisto.SetMinimum( 1. )

        hBkg.SetFillStyle(0)
        hBkg.SetLineColor(kBlack)
        hBkg.SetLineStyle(1)
        hBkg.SetLineWidth(1)
        #hBkg.SetFillStyle(3004)
        hBkg.SetFillColor( kRed )
        hBkg.Draw("same")

        stackHisto.GetYaxis().SetTitle( 'Events / '+str(binWidth)+' GeV' )
        stackHisto.GetXaxis().SetTitle( axisX )

        tmpHisto = {}
        for sample in signalHistos:
            tmpHisto[ sample ] = signalHistos[ sample ].Clone()
            tmpHisto[ sample ].SetFillColor(0)
            tmpHisto[ sample ].SetLineStyle(2)
            tmpHisto[ sample ].SetLineWidth(3)
            #tmpHisto[ sample ].Draw("hist same")

        legend.Draw()
        if args.process.endswith('Data'):
            stackHisto.SetMaximum( max(hBkg.GetMaximum(), dataHistos['AllData'].GetMaximum() )*1.5 )
            dataHistos['AllData'].SetMarkerStyle(8)
            dataHistos['AllData'].Draw('E same')
            CMS_lumi.extraText = "Preliminary"
            CMS_lumi.relPosX = 0.14
            CMS_lumi.CMS_lumi( pad1, 4, 0)
        else:
            stackHisto.SetMaximum( hBkg.GetMaximum()*1.5 )
            stackHisto.GetYaxis().SetTitleOffset( 0.8 )
            CMS_lumi.CMS_lumi( canvas[outputFileName], 4, 0)


        if args.process.endswith('Data'):
           pad2.cd()
           pad2.SetGrid()
           pad2.SetTopMargin(0)
           pad2.SetBottomMargin(0.3)

           tmpPad2= pad2.DrawFrame(xmin,0.5,xmax,1.5)
           tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
           tmpPad2.GetXaxis().SetTitle( axisX )
           tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
           tmpPad2.GetYaxis().CenterTitle()
           tmpPad2.SetLabelSize(0.12, 'x')
           tmpPad2.SetTitleSize(0.12, 'x')
           tmpPad2.SetLabelSize(0.12, 'y')
           tmpPad2.SetTitleSize(0.12, 'y')
           tmpPad2.SetNdivisions(505, 'x')
           tmpPad2.SetNdivisions(505, 'y')
           pad2.Modified()
           hRatio = TGraphAsymmErrors()
           hRatio.Divide( dataHistos[ 'AllData' ], hBkg, 'pois' )
           hRatio.SetMarkerStyle(8)
           hRatio.Draw('P')
           hRatioStatErr = hBkg.Clone()
           hRatioStatErr.Divide( hBkg )
           hRatioStatErr.SetFillColor(kBlack)
           hRatioStatErr.SetFillStyle(3004)
           hRatioStatErr.Draw("same e2")

        canvas[outputFileName].SaveAs( 'Plots/'+outputFileName )

    else:

        tdrStyle.SetPadRightMargin(0.05)
        canvas[outputFileName]= TCanvas('c1', 'c1', 750, 500 )
        if log: canvas[outputFileName].SetLogy()
        signalHistos[next(iter(signalHistos))].GetYaxis().SetTitleOffset(1.0)
        signalHistos[next(iter(signalHistos))].GetYaxis().SetTitle( ( 'Normalized / '+str(int(binWidth))+' GeV' ) )
        if xmax: signalHistos[next(iter(signalHistos))].GetXaxis().SetRangeUser( xmin, xmax )
        signalHistos[next(iter(signalHistos))].Draw('hist')
        for signalSamples in signalHistos: signalHistos[ signalSamples ].Draw('hist same')
        for bkgSamples in bkgHistos: bkgHistos[ bkgSamples ].Draw('hist same')
        if 'DATA' in args.process:
                dataHistos[ 'DATA' ].SetMarkerStyle(8)
                dataHistos[ 'DATA' ].Draw('same')
                CMS_lumi.extraText = ""#"Preliminary"
        signalHistos[next(iter(signalHistos))].SetMaximum( 1.1 * max( maxList ) )

        if not 'DATA' in args.process: CMS_lumi.lumi_13TeV = ''
        CMS_lumi.relPosX = 0.11
        CMS_lumi.CMS_lumi(canvas[outputFileName], 4, 0)
        legend.Draw()

        canvas[outputFileName].SaveAs( 'Plots/'+outputFileName )
    del canvas[outputFileName]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='1D', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
    parser.add_argument('-v', '--version', action='store', default='v0', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument('-c', '--cut', action='store', nargs='+', default='2J2WdeltaR', help='cut, example: "2J 2J2W"' )
    parser.add_argument('-s', '--single', action='store', default='all', help='single histogram, example: massAve_cutDijet.' )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-j', '--json', action='store_true', default=False, dest='json',  help='Plot from json (true) or not (false)' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument('-f', '--final', action='store_true', default=False, dest='final',  help='If plot is final' )
    parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if not os.path.exists('Plots/'): os.makedirs('Plots/')
    if args.year.endswith('2016'): args.lumi = 35920.
    elif args.year.endswith('2017'): args.lumi = 41530.
    elif args.year.endswith('2018'): args.lumi = 59740.
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+args.year

    taulabX = 0.90
    taulabY = 0.85
    massMinX = 0
    massMaxX = 400

    plotList = [

            #### 0: process, 1: name of histogram in file, 2: X label, 3: minX, 4: maxX, 5: rebining
            [ 'signalBkg', 'PV_npvsGood', 'Number of PV', 0, 100, 2 ],
            [ 'signalBkg', 'nleps', 'Number of leptons', 0, 5, 1 ],
            [ 'signalBkg', 'lepton_pt', 'Lepton pT [GeV]', 0, 500, 2 ],
            [ 'signalBkg', 'lepton_eta', 'Lepton #eta', -3, 3, 2 ],
            [ 'signalBkg', 'njets', 'Number of AK4 jets', 0, 10, 1 ],
            [ 'signalBkg', 'btags', 'Number of AK4 b-tagged jets', 0, 10, 1 ],
            [ 'signalBkg', 'leading_jet_pt', 'Leading AK4 jet pT [GeV]', 0, 500, 1 ],
            [ 'signalBkg', 'leading_jet_eta', 'Leading AK4 jets #eta', -3, 3, 2 ],
            [ 'signalBkg', 'nfatjets', 'Number of AK8 jets', 0, 10, 1 ],
            [ 'signalBkg', 'met', 'MET [GeV]', 0, 800, 2 ],
            [ 'signalBkg', 'lepWMass', 'Leptonic W mass [GeV]', 50, 150, 1 ],
            [ 'signalBkg', 'lepWPt', 'Leptonic W pT [GeV]', 0, 300, 2 ],
            [ 'signalBkg', 'hadWMass', 'Hadronic W mass [GeV]', 50, 150, 1 ],
            [ 'signalBkg', 'hadWPt', 'Hadronic W pT [GeV]', 0, 300, 2 ],
            [ 'signalBkg', 'leadAK8JetPt', 'Leading AK8 jet pT [GeV]', 100, 1500, 5 ],
            [ 'signalBkg', 'leadAK8JetEta', 'Leading AK8 jet #eta', -3, 3, 2 ],
            #[ 'signalBkg', 'leadAK8JetMass', 'Leading AK8 jet mass [GeV]', 0, 500, 5 ],
            [ 'signalBkg', 'leadAK8JetMass', 'Leading AK8 jet mass [GeV]', 90, 160, 5 ],
            [ 'signalBkg', 'leadAK8JetTau21', 'Leading AK8 jet #tau_{21}', 0, 1, 2 ],
            [ 'signalBkg', 'leadAK8JetHbb', 'Leading AK8 jet Hbb', 0, 1, 2 ],
            [ 'signalBkg', 'deltaRlepWHiggs', '#Delta R(W_{lep}, AK8 jet)', 0, 5, 1 ],
            [ 'signalBkg', 'deltaRhadWHiggs', '#Delta R(W_{had}, AK8 jet)', 0, 5, 1 ],

            [ 'qual', 'nPVs', 'Number of PV', 0, 100, 1,  0.85, 0.70 ],
            [ 'qual', 'PV_npvsGood', 'Number of PV', 0, 100, 2,  0.85, 0.70 ],
            [ 'qual', 'nleps', 'Number of leptons', 0, 10, 1,  0.85, 0.70 ],
            [ 'qual', 'lepton_pt', 'Lepton pT [GeV]', 0, 500, 2,  0.85, 0.70 ],
            [ 'qual', 'lepton_eta', 'Lepton #eta', -3, 3, 2,  0.85, 0.70 ],
            [ 'qual', 'lepton_phi', 'Lepton #phi', -3, 3, 4,  0.85, 0.70 ],
            [ 'qual', 'njets', 'Number of AK4 jets', 0, 10, 1,  0.85, 0.70 ],
            [ 'qual', 'jets_pt', 'AK4 jets pT [GeV]', 0, 500, 1,  0.85, 0.70 ],
            [ 'qual', 'jets_eta', 'AK4 jets #eta', -3, 3, 2,  0.85, 0.70 ],
            [ 'qual', 'jets_phi', 'AK4 jets #phi', -3, 3, 2,  0.85, 0.70 ],
            [ 'qual', 'nBjets', 'Number of AK4 bjets', 0, 10, 1,  0.85, 0.70 ],
            [ 'qual', 'nAK8jets', 'Number of AK8 jets', 0, 10, 1,  0.85, 0.70 ],
            [ 'qual', 'METPt', 'MET [GeV]', 0, 800, 2,  0.85, 0.70 ],
            [ 'qual', 'lepWMass', 'Leptonic W mass [GeV]', 50, 250, 1,  0.85, 0.70 ],
            [ 'qual', 'lepWPt', 'Leptonic W pT [GeV]', 0, 300, 2,  0.85, 0.70 ],
            [ 'qual', 'resolvedWCandMass', 'Hadronic W mass [GeV]', 0, 200, 1,  0.85, 0.70 ],
            [ 'qual', 'resolvedWCandPt', 'Hadronic W pT [GeV]', 0, 300, 2,  0.85, 0.70 ],
            [ 'qual', 'leadAK8JetPt', 'Leading AK8 jet pT [GeV]', 100, 1500, 5, 0.85, 0.70 ],
            [ 'qual', 'leadAK8JetMass', 'Leading AK8 jet mass [GeV]', 30, 250, 2, 0.85, 0.70 ],
            [ 'qual', 'leadAK8JetTau21', 'Leading AK8 jet #tau_{21}', 0, 1, 2, 0.85, 0.70 ],
            [ 'qual', 'leadAK8JetHbb', 'Leading AK8 jet Hbb', 0, 1, 2, 0.85, 0.70 ],

            [ 'simple', 'PV_npvsGood', 'Number of PV', 0, 100, 2,  0.85, 0.70, False ],
    ]

    if 'all' in args.single: Plots = [ x[1:] for x in plotList if ( ( args.process.startswith( x[0] ) ) ) ]
    else: Plots = [ y[1:] for y in plotList if ( args.process.startswith(y[0]) and y[1].startswith(args.single) )  ]

    VER = args.version.split('_')[1] if '_' in args.version else args.version
    bkgFiles, signalFiles, dataFiles = rootHistograms( VER, args.lumi, args.year )
    folder = '/eos/home-a/algomez/tmpFiles/hepacc/results/'+args.year+'/v14/'+args.version+"/nominal/"
    #folder = '/afs/cern.ch/work/d/druini/public/hepaccelerate/results/'+args.year+'/v14/'+args.version+"/"
    if args.json:
        print '|-----> Ignore errors above this.'
        args.version='hepacc_'+args.version

    if args.norm:
        bkgFiles.pop('TTTo2L2Nu', None)
        #bkgFiles.pop('ST_s-channel', None)
        #bkgFiles.pop('ST_t-channel', None)
        #bkgFiles.pop('ST_tW_top', None)
        bkgFiles.pop('WW', None)
        bkgFiles.pop('WZ', None)
        bkgFiles.pop('ZZ', None)
        bkgFiles.pop('TTGJets', None)

    for i in Plots:

        ######### qual makes all bkg vs data plots
        if ( 'qual' in args.process ):
            for icut in args.cut:
                plotQuality(
                    i[0]+'_'+icut, i[1], i[2], i[3], i[4], i[5], i[6] )

        ######## signalBkg makes stuck plots of bkg vs signal. If adding signalBkgData, it includes data
        elif args.process.startswith( 'signalBkg'):
            for icut in args.cut:
                plotSignalBkg( i[0]+'_'+icut, i[2], i[3], i[4], log=args.log, axisX=i[1], Norm=args.norm, ext=args.ext)

        ######### simple comparison plots for testing
        elif ( 'simple' in args.process ):
            for icut in args.cut:
                plotSimpleComparison(
                    jsonToTH1( '/afs/cern.ch/user/a/algomez/cernbox/tmpFiles/hepacc/results/v02/2017/out_TTToSemiLeptonic.json', [i[0]+'_'+icut]  ), 'allWeights',
                    jsonToTH1( '/afs/cern.ch/user/a/algomez/cernbox/tmpFiles/hepacc/results/v02/2017_noLepWeight/2017/out_TTToSemiLeptonic.json', [i[0]+'_'+icut]  ), 'noLepWeights',
                    i[0]+'_'+icut, xmin=i[2], xmax=i[3], rebinX=i[4], log=i[7], axisX=i[1] )


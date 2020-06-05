import ROOT as rt
import tdrstyle
import sys, os, argparse
#from ctypes import c_double as double # possible replacement for ROOT.Double
from math import sqrt
from glob import glob
from pdb import set_trace

rt.gROOT.SetBatch()
rt.gROOT.ForceStyle()
tdrstyle.setTDRStyle()
rt.gStyle.SetOptStat(0)

rt.tdrStyle.SetPadLeftMargin(0.15)
rt.tdrStyle.SetPadRightMargin(0.05)

def TH1fromRooObj(RooObj, errorBand=None):
  if type(RooObj) == rt.RooHist:
    bins = range(RooObj.GetN())
  elif type(RooObj) == rt.RooCurve:
    bins = range(RooObj.GetN())[3:-4:2]
  else:
    raise Exception(f'unknown object type: {type(RooObj)}')

  hist = [[rt.Double(0),rt.Double(0)] for b in bins]
  for n,b in enumerate(bins):
    RooObj.GetPoint(b, *hist[n])
  binWidth = round(hist[1][0] - hist[0][0])
  nbins    = len(bins)

  if errorBand is None:
    err = [sqrt(v[1]) for v in hist]
  else:
    idxUp   = range(int(errorBand.GetN()/2))[3:-4:2]
    idxDown = range(int(errorBand.GetN()/2),errorBand.GetN())[-4:3:-2]
    assert( len(idxUp)==nbins )
    assert( len(idxDown)==nbins )
    errUp   = [[rt.Double(0),rt.Double(0)] for b in range(nbins)]
    errDown = [[rt.Double(0),rt.Double(0)] for b in range(nbins)]
    for n,i in enumerate(idxUp):
      errorBand.GetPoint(i, *errUp[n])
    for n,i in enumerate(idxDown):
      errorBand.GetPoint(i, *errDown[n])
      assert( errUp[n][0]-errDown[n][0]<.5 ) # check that up and down errors correspond to the same bin
      assert( hist[n][1]-(errUp[n][1]+errDown[n][1])/2 < 1e-3 ) # check that histogram is in the middle of the error band
    err = [ (errUp[i][1] - errDown[i][1])/2 for i in range(nbins) ]

  if type(RooObj) == rt.RooHist:
    h_start = round(hist[0][0],2)-binWidth/2
    h_stop  = round(hist[-1][0],2)+binWidth/2
  else:
    h_start = round(hist[0][0],2)
    h_stop  = round(hist[-1][0],2)+binWidth
  TH1 = rt.TH1D('','', nbins, h_start, h_stop)
  for i in range(nbins):
    TH1.Fill( *hist[i] )
    TH1.SetBinError( i+1, err[i] ) # stupid root starts indexing at 1 instead of 0
  return TH1.Clone()

def plotRhalphaShapes(rootFilePath, nptbins):
  if not rootFilePath.endswith(".root"):
    raise Exception(f'Must supply ROOT filename, but got {rootFilePath}')
  if not (os.path.sep in rootFilePath): #if the file is in the current directory
    outdir = os.getcwd()
  else:
    outdir = os.path.join(*rootFilePath.split(os.path.sep)[:-1])
  outdir = os.path.join(outdir,'plots')
  if not os.path.exists(outdir):
    os.mkdir(outdir)

  rootfile = rt.TFile(rootFilePath)

  # draw options
  drawOpts = {
    'data'       : 'same p',
    'errs'       : 'af',
    'signal'     : 'same',
    'bkg'        : 'same',
    'sigPlusBkg' : 'same'
  }
  #legend labels
  leglab = {
    'data'       : 'Data',
    'signal'     : 'Signal',
    'bkg'        : 'Background',
    'sigPlusBkg' : 'Signal+Bkg'
  }
  #legend options
  legopts = {
    'data'       : 'ep',
    'signal'     : 'lp',
    'bkg'        : 'lp',
    'sigPlusBkg' : 'lp',
  }

  for ibin in range(nptbins):
    for region in ['Pass','Fail']:
      for SorB in ['s','b']:
        p = f'ptbin{ibin}{region}_msd_fit_{SorB}'
        can = rt.TCanvas('c2', 'c1',  10, 10, 750, 750 )
        pad1 = rt.TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
        pad2 = rt.TPad("pad2", "Pull",0,0.00,1.00,0.28,-1)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()

        leg = rt.TLegend(0.70,0.7,0.90,0.87)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        plot   = rootfile.Get(p)
        b_suffix = '_model_bonly__' if SorB=='b' else ''
        RooObj = {
          'data'       : plot.getHist(f'h_ptbin{ibin}{region}'),
          'errs'       : plot.getCurve(f'pdf_binptbin{ibin}{region}_{b_suffix}Norm[msd]_errorband'),
          'signal'     : plot.getCurve(f'pdf_binptbin{ibin}{region}_{b_suffix}Norm[msd]_Comp[shapeSig*]'),
          'bkg'        : plot.getCurve(f'pdf_binptbin{ibin}{region}_{b_suffix}Norm[msd]_Comp[shapeBkg*]'),
          'sigPlusBkg' : plot.getCurve(f'pdf_binptbin{ibin}{region}_{b_suffix}Norm[msd]')
        }

        rmin = RooObj['data'].GetXaxis().GetXmin()
        rmax = RooObj['data'].GetXaxis().GetXmax()
        xmin = (rmax+rmin)/2 - (rmax-rmin)/2.4
        xmax = (rmax+rmin)/2 + (rmax-rmin)/2.4
        RooObj['errs'].GetYaxis().SetTitle( f'Events / {RooObj["data"].getNominalBinWidth()} GeV' )
        #RooObj['errs'].GetYaxis().SetTitleOffset( 0.5 )
        #RooObj['errs'].GetYaxis().SetLabelSize(0.12)
        #RooObj['errs'].GetYaxis().SetTitleSize(0.12)
        for s in ['errs', 'signal', 'bkg', 'sigPlusBkg', 'data']:
          RooObj[s].GetXaxis().SetLabelOffset(999)
          RooObj[s].GetXaxis().SetLabelSize(0)
          RooObj[s].GetXaxis().SetLimits(xmin,xmax)
          RooObj[s].Draw(drawOpts[s])
          if s in leglab:
            leg.AddEntry( RooObj[s], leglab[s], legopts[s] )

        leg.Draw()
    #    #ymax = h['data_obs'].GetMaximum() + 1.5*sqrt(h['data_obs'].GetMaximum())
    #    #h[namesInRoot[0]].GetYaxis().SetRangeUser(0, ymax)


        hRatio = rt.TGraphAsymmErrors()
        data_h = TH1fromRooObj(RooObj['data'])
        sigPlusBkg_h = TH1fromRooObj(RooObj['sigPlusBkg'],RooObj['errs'])
        hRatio.Divide( data_h,sigPlusBkg_h, 'pois' )
        hRatioStatErr = sigPlusBkg_h.Clone()
        hRatioStatErr.Divide( sigPlusBkg_h )
        hRatioStatErr.SetFillColor(rt.kBlack)
        hRatioStatErr.SetFillStyle(3004)

        pad2.cd()
        rt.gStyle.SetOptFit(1)
        pad2.SetGrid()
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.3)
        tmpPad2= pad2.DrawFrame(xmin, 0.5, xmax,1.5)
        tmpPad2.GetXaxis().SetTitle( 'softdrop mass [GeV]' )
        tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
        tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
        tmpPad2.GetYaxis().CenterTitle()
        tmpPad2.SetLabelSize(0.12, 'x')
        tmpPad2.SetTitleSize(0.12, 'x')
        tmpPad2.SetLabelSize(0.12, 'y')
        tmpPad2.SetTitleSize(0.12, 'y')
        tmpPad2.SetNdivisions( data_h.GetNdivisions('x'), 'x')
        tmpPad2.SetNdivisions(505, 'y')
        pad2.Modified()
        hRatio.SetMarkerStyle(8)
        hRatio.Draw('P')
        hRatioStatErr.Draw('same e2')

        for ext in ['pdf','png']:
          can.SaveAs( os.path.join(outdir,f'{p}.{ext}') )
        del can

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--rootfile', action='store', dest='rootfile', help='root file containing the output from CombineHarvester')
  parser.add_argument('-n', '--nptbins', action='store', default=2, type=int, dest='nptbins', help='number of pt bins')

  try: args = parser.parse_args()
  except:
    parser.print_help()
    sys.exit(0)

  plotRhalphaShapes(args.rootfile, args.nptbins)
  #for f in glob('output/201*/v05/*/mc_msd100to150*/fitDiagnostics.root'):
  #  try:
  #    plotRhalphaShapes(f, args.nptbins)
  #  except:
  #    print(f'{f} failed')

import ROOT as rt
import tdrstyle
import CMS_lumi as CMS_lumi
import sys, os, argparse
#from ctypes import c_double as double # possible replacement for ROOT.Double
from math import sqrt
import numpy as np
from glob import glob
from pdb import set_trace

rt.gROOT.SetBatch()
rt.gROOT.ForceStyle()
tdrstyle.setTDRStyle()
rt.gStyle.SetOptStat(0)

rt.tdrStyle.SetPadLeftMargin(0.15)
rt.tdrStyle.SetPadRightMargin(0.05)

def TH1fromRooObj(RooObj, errorBand=None, returnHistArray=False):
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

  if returnHistArray: return np.array(hist)

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
    'data'       : 'Data' if args.isData else 'MC as Data',
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

        CMS_lumi.cmsTextOffset = 0.0
        CMS_lumi.relPosX = 0.13
        CMS_lumi.CMS_lumi(pad1, 4, 0)
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

def plotTF(rootFilePath, nptbins):
  if not rootFilePath.endswith(".root"):
    raise Exception(f'Must supply ROOT filename, but got {rootFilePath}')
  if not (os.path.sep in rootFilePath): #if the file is in the current directory
    outdir = os.getcwd()
  else:
    outdir = os.path.join(*rootFilePath.split(os.path.sep)[:-1])
  outdir = os.path.join(outdir,'plots')
  if not os.path.exists(outdir):
    os.mkdir(outdir)

  import matplotlib.pyplot as plt
  import mplhep as hep
  plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
  plt.switch_backend('agg')
  rootfile = rt.TFile(rootFilePath)
  ptdeg  = rootFilePath.split('polyDegs')[1][0]
  rhodeg = rootFilePath.split('polyDegs')[1][1]
  tf = []
  for ibin in range(2):#[1,0]:
    msd_hist = {}
    for region in ['Pass','Fail']:
      SorB = 's'
      p = f'ptbin{ibin}{region}_msd_fit_{SorB}'
      b_suffix = '_model_bonly__' if SorB=='b' else ''
      hist_name = f'pdf_binptbin{ibin}{region}_{b_suffix}Norm[msd]'
      msd_hist[region] = TH1fromRooObj( rootfile.Get(p).getCurve(hist_name), returnHistArray=True )
    tf.append(msd_hist['Pass'][:,1]/msd_hist['Fail'][:,1])
  tf = np.array(tf)

  msdbins  = np.round(msd_hist['Pass'][:,0])
  msdbins  = np.append(msdbins, msdbins[-1]+np.diff(msdbins)[0])
  msdsampl = msdbins[:-1] + 0.5*np.diff(msdbins)
  ptbins   = np.array([250,300,2000])
  ptsampl  = ptbins[:-1] + 0.3*np.diff(ptbins)
  sampling = np.meshgrid(msdsampl, ptsampl)
  rho = 2*np.log(sampling[0]/sampling[1])
  mask = (rho>-1.2) | (rho<-6)
  fig, ax = plt.subplots()
  ax = hep.cms.cmslabel(data=args.isData, paper=False, year=args.year, ax=ax, loc=1)
  ptbins[-1] = 400
  vmin = np.floor(100*min(tf[~mask]))/100
  vmax = np.ceil(100*max(tf[~mask]))/100
  hist2d = hep.hist2dplot(tf.T, msdbins, ptbins, vmin=vmin, vmax=vmax, cmap='inferno', ax=ax)
  for ibin in range(len(ptbins)-1):
    pt = ptbins[ibin]
    ptnext = ptbins[ibin+1]
    ax.fill_between(msdbins, np.full_like(msdbins,pt), np.full_like(msdbins,ptnext), where=np.append(mask[ibin],True), facecolor='w', edgecolor='k', linewidth=0, hatch='xx')
  ax.set_title(fr'MC TF, $\deg( p_\mathrm{{T}}, \rho ) = ({ptdeg},{rhodeg})$', pad=9, fontsize=22, loc='left')
  #ax.set_title('2017',pad=9, fontsize=22, loc='right')
  ax.set_xlabel(r'Jet $\mathrm{m_{SD}}$ [GeV]', ha='right', x=1)
  ax.set_ylabel(r'Jet $\mathrm{p_{T}}$ [GeV]', ha='right', y=1)
  ax.yaxis.set_ticks([250,300])
  cax = fig.axes[1]
  cax.set_ylabel('TF', ha='right', y=1)
  for ext in ['pdf','png']:
    fig.savefig(os.path.join(outdir,f'TF.{ext}'))

  ######### trying the pt-rho plot
#  ptbins[-1] = 2000
#  msdptgrid = np.meshgrid(msdbins,ptbins)
#  rho_vertices = 2*np.log(msdptgrid[0]/msdptgrid[1])# these are the vertices of the bins. Idea: make rectangular bins centered in the centre of these skewed bins and with the correct bin width
#  #plt.scatter(rho,msdptgrid[1])
#  rho_ptbin0 = 2*np.log(msdsampl/ptsampl[0])-np.diff(rho_vertices)[0]/2
#  rho_ptbin0 = np.append(rho_ptbin0, rho_ptbin0[-1]+np.diff(rho_vertices)[0][-1]/2)
#  rho_ptbin1 = 2*np.log(msdsampl/ptsampl[1])-np.diff(rho_vertices)[0]/2
#  rho_ptbin1 = np.append(rho_ptbin1, rho_ptbin1[-1]+np.diff(rho_vertices)[0][-1]/2)
#  fig, ax = plt.subplots()
#  rhobins = np.append(rho_ptbin1, rho_ptbin0)
#  tf = np.array([np.append(np.zeros(11),tf[0]),np.append(tf[1],np.zeros(11))])
#  ptbins[-1] = 400
#  hist2d = hep.hist2dplot(tf.T, rhobins, ptbins, vmin=vmin, vmax=vmax, cmap='inferno', ax=ax)
#  ax.set_xlim(right=-1.2)
#  for ibin in range(len(ptbins)-1):
#    pt = ptbins[ibin]
#    ptnext = ptbins[ibin+1]
#    wh = np.append(tf[ibin]==0,True) if ibin==1 else np.append(True,tf[ibin]==0)
#    ax.fill_between(rhobins, np.full_like(rhobins,pt), np.full_like(rhobins,ptnext), where=wh, facecolor='w', edgecolor='k', linewidth=0, hatch='xx')
#  ax.set_title(fr'MC TF, $\deg( p_\mathrm{{T}}, \rho ) = ({ptdeg},{rhodeg})$', pad=9, fontsize=22, loc='left')
#  ax.set_title('2017',pad=9, fontsize=22, loc='right')
#  ax.set_xlabel(r'Jet $\rho$', ha='right', x=1)
#  ax.set_ylabel(r'Jet $\mathrm{p_{T}}$ [GeV]', ha='right', y=1)
#  ax.yaxis.set_ticks([250,300])
#  cax = fig.axes[1]
#  cax.set_ylabel('TF', ha='right', y=1)
#  #set_trace()
#  #hep.hist2dplot(tf[0], rho_ptbin0, ptbins[0:2], ax=ax)
#  #hep.hist2dplot(tf[1], rho_ptbin0, ptbins[1:3], ax=ax)
#  for ext in ['pdf','png']:
#    fig.savefig(f'test/TF_rho.{ext}')


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--rootfile', action='store', dest='rootfile', help='root file containing the output from CombineHarvester')
  parser.add_argument('-n', '--nptbins', action='store', default=2, type=int, dest='nptbins', help='number of pt bins')
  parser.add_argument('-y', '--year', default='2017', type=str, help='year to process, in file paths')
  parser.add_argument('-d', '--isData', action='store_true', default=False, help='flag to run on data or mc')

  try: args = parser.parse_args()
  except:
    parser.print_help()
    sys.exit(0)

  if args.year.endswith('2016'): args.lumi = 35920.
  elif args.year.endswith('2017'): args.lumi = 41530.
  elif args.year.endswith('2018'): args.lumi = 59740.
  CMS_lumi.extraText = "Preliminary" if args.isData else "Simulation Preliminary"
  CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, "+args.year+" (13 TeV)"

  plotTF(args.rootfile, args.nptbins)
  plotRhalphaShapes(args.rootfile, args.nptbins)
  #for f in glob('output/201*/v05/*/mc_msd100to150*/fitDiagnostics.root'):
  #  try:
  #    plotRhalphaShapes(f, args.nptbins)
  #  except:
  #    print(f'{f} failed')

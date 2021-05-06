from __future__ import print_function
import ROOT as rt
import sys
from pdb import set_trace

#print(f'loading file {sys.argv[1]}...')
inFile = rt.TFile(sys.argv[1])
ws     = inFile.Get('ws')
ws.Print()

#msd = ws.var('msd')

#canvas = rt.TCanvas()
#xframe = msd.frame()
#ws.data('TTH_PTH_GT300').plotOn(xframe)
#ws.data('TTH_PTH_GT300_CMS_ttHbb_jms_2016Down').plotOn(xframe,rt.RooFit.LineColor(rt.kBlue),rt.RooFit.MarkerColor(rt.kBlue))
#ws.data('TTH_PTH_GT300_CMS_ttHbb_jms_2016Up').plotOn(xframe,rt.RooFit.LineColor(rt.kRed),rt.RooFit.MarkerColor(rt.kRed))
#xframe.Draw()
#set_trace()

print(ws.data('TTH_PTH_GT300_CMS_ttHbb_jms_2016Down').numEntries())
print(ws.data('TTH_PTH_GT300_CMS_ttHbb_jms_2016Down').sumEntries())

h = inFile.Get('TTH_PTH_GT300__CMS_ttHbb_jms_2016Down')
print(h.GetEntries())
print(h.Integral())

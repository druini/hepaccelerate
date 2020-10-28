import numpy as np
from hepaccelerate.utils import Histogram, Results
from glob import glob
import json,os,argparse
from pdb import set_trace

flist = glob('results/201*/v12/met20_btagDDBvL086/nominal/btagEfficiencyMaps/out_btagEfficiencyMaps_*json')

def divide(h1,h2):
  contents    = h1.contents/h2.contents
  contents_w2 = h1.contents_w2/h2.contents_w2
  edges       = h1.edges
  return Histogram(contents, contents_w2, edges)

for fn in flist:
  with open(fn) as f:
    data = json.load(f)
  for h in data:
    data[h] = Histogram( *data[h].values() )

  for flav in ['b','l','lc']:
    for var in ['central','updown']:
      data[f'eff_flav{flav}_{var}'] = divide( data[f'btags_flav{flav}_{var}'], data[f'total_flav{flav}_{var}'] )

  ret = Results(data)
  ret.save_json(fn)


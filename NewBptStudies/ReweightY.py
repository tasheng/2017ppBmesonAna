import ROOT as r
import numpy as np
import pandas as pd
import array
from os import listdir
from os.path import join


# if 'can' not in locals():
#     can = r.TCanvas('can', 'can')
r.gROOT.SetBatch(True)

foutstr = 'FONLLFine/FONLL_rapidity.root'
fout = r.TFile(foutstr, 'recreate')

rapdir = 'FONLLFine/rapidity/'
ylim = ['0', '0p5', '1', '1p5', '2', '2p4']

glist = []
# for fin in listdir(rapdir):
#     fin = join(rapdir, fin)
for i in range(len(ylim) - 1):
    fstr = f'fonll_pp_B_5p03TeV_y-{ylim[i]}-{ylim[i+1]}.dat'
    fin = rapdir + fstr
    colnames = ['pt', 'cross']
    df = pd.read_csv(fin, names=colnames, delimiter='\s+')
    # select the pT range we are using
    df = df[df.pt >= 5]
    x = df['pt'].to_numpy(dtype=np.float64)
    y = df['cross'].to_numpy(dtype=np.float64)
    g = r.TGraph(len(x), x, y)
    g.SetName(f'gaeSigmaBplus{i}')
    g.Write()
    glist.append(g)
fout.Close()

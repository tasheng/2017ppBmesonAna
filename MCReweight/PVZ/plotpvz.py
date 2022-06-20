import ROOT as r
can = r.TCanvas('can', 'can')

fstr = "pvzmc0.log"
with open(fstr, "r") as f:
    x = f.readlines()
print(x[-1])
fnstr = x[-1].replace('PVz', 'x')
print(fnstr)
fnratio = r.TF1("ratio", fnstr, -30, 30)
fnratio.Draw()
can.Show()

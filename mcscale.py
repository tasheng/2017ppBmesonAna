# Obtain the scale factors of low pt bins from MC
import ROOT as r
# import uproot

bpfilestr = "CutSkim/BPMC.root"
bptree = 'ntKp'

bsfilestr = "CutSkim/BsMC.root"
bstree = 'ntphi'

ptbinsbp = [5, 7, 10]
ptbinsbs = [7, 10]
ynarrow = 'By < 2.4 && By > 1.5'
yall = 'By < 2.4'

sample = ['BP', 'Bs']
dbp = r.RDataFrame(bptree, bpfilestr)
dbs = r.RDataFrame(bstree, bsfilestr)
dflist = [dbp, dbs]
ptblist = [ptbinsbp, ptbinsbs]
# get results for all bins
s, ptini, ptfin, narrow, scaled, factor = 'sample', 'ptini', 'ptfin', 'narrow', 'scaled', 'factor'
# print(f"{s:>8} {ptini:>8} {ptfin:>8} {narrow:>8} {scaled:>8} {factor:>8}")
print("{stype:>8} {ptlow:>8} {pthigh:>8} {ent_narrow:>8} {ent_all:>8} {factor:>8}".format(
    stype=s, ptlow=ptini, pthigh=ptfin, ent_narrow=narrow, ent_all=scaled, factor=factor))
for stype, df, ptbins in (zip(sample, dflist, ptblist)):
    for ipt in range(len(ptbins) - 1):
        ptlow = ptbins[ipt]
        pthigh = ptbins[ipt + 1]
        ptcut = 'Bpt > {ptlow} && Bpt < {pthigh}'.format(ptlow=ptlow, pthigh=pthigh)
        dpt = df.Filter(ptcut)
        entries = df.Count().GetValue()
        # print(f"{entries} events in {ptlow} < pt < {pthigh}")

        ent_narrow = dpt.Filter(ynarrow).Count().GetValue()
        ent_all = dpt.Filter(yall).Count().GetValue()
        factor = ent_all / float(ent_narrow)
        # print(f"narrow: {ent_narrow}, scaled: {ent_all}, factor: {ent_all / ent_narrow :.2f}")
        # print(f"{stype:>8} {ptlow:>8} {pthigh:>8} {ent_narrow:>8} {ent_all:>8} {factor:8.4g}")
        print("{stype:>8} {ptlow:>8} {pthigh:>8} {ent_narrow:>8} {ent_all:>8} {factor:8.4g}".format(
            stype=stype, ptlow=ptlow, pthigh=pthigh, ent_narrow=ent_narrow, ent_all=ent_all, factor=factor))

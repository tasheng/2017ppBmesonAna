import sys
import ROOT as r
import numpy as np

# if 'can' not in locals():
#     can = r.TCanvas('can', 'can')
#     r.gROOT.ProcessLine('.L ~/ecal/color.hpp')

r.gStyle.SetOptStat(0)


btypes = ['BP', 'Bs']
btypes_low = ['Bp', 'Bs']
# all: all cuts, nominal
# nommc: no pre-selections




cut = "c3m"
# cutlist = ['c3m', 'cchi2', 'cnjmass', 'c3mdls', 'c3mjmass']
# for cut in cutlist:
for bind, btype in enumerate(btypes):
    # if bind==0: continue
    # if bind==1: continue

    fio = r.TFile('2DMapSyst/OutFiles/{btype}Error2D.root'.format(btype=btype), "update")
    hBDT = fio.Get('BDTSyst')
    hMCData = hBDT.Clone('MCDataSyst')
    hPre = hBDT.Clone('PreSelectionSyst')
    fcorr_yield_nom = r.TFile('{btype}pre/EffAna/FinalFiles/{btype}_all.root'.format(btype=btype))
    fcorr_yield_pre = r.TFile('{btype}pre/EffAna/FinalFiles/{btype}_{cut}.root'.format(
        btype=btype, cut=cut))

    fraw_yield_nom = r.TFile('henri2022_pre/ROOTfiles/yields_{blow}_binned_pt_all.root'.format(
        blow=btypes_low[bind],))
    fraw_yield_pre = r.TFile('henri2022_pre/ROOTfiles/yields_{blow}_binned_pt_{cut}.root'.format(
        blow=btypes_low[bind], cut=cut))
    # fraw_yield_nommc = r.TFile('henri2022_pre/ROOTfiles/yields_Bp_binned_pt_nommc.root')

    feff_nom = r.TFile('{btype}pre/EffAna/NewEff2DMaps/EffFineBDTall.root'.format(
        btype=btype, cut=cut))
    feff_pre = r.TFile('{btype}pre/EffAna/NewEff2DMaps/EffFineBDT{cut}.root'.format(
        btype=btype, cut=cut))

    fcorrs = [fcorr_yield_nom, fcorr_yield_pre]
    feffs = [feff_nom, feff_pre]
    fraws = [fraw_yield_nom, fraw_yield_pre]
    # frawsmc = [fraw_yield_nom, fraw_yield_nommc]

    yds = [f.Get('hPtSigma') for f in fcorrs]
    fineeffs = [f.Get('invEff2D') for f in feffs]
    raws = [f.Get('hPt') for f in fraws]
    # rawsmc = [f.Get('hPt') for f in frawsmc]
    effs = [f.Get('hInvEff') for f in fcorrs]

    nbins = 7

    def comp(histlist, inverse=False, update=False):
        for ib in range(1, raws[0].GetNbinsX() + 1):
            nom = histlist[0].GetBinContent(ib)
            pre = histlist[1].GetBinContent(ib)
            nomerr = histlist[0].GetBinError(ib)
            preerr = histlist[1].GetBinError(ib)
            if inverse:
                nom, pre = 1/nom, 1/pre
            ptlow = histlist[0].GetBinLowEdge(ib)
            pthi = histlist[0].GetBinLowEdge(ib + 1)

            diff = abs(pre / nom * 100 - 100)
            bdterr = hMCData.GetBinContent(ib)
            if (update):
                hPre.SetBinContent(ib, diff)
            if (diff > bdterr and update):
                hMCData.SetBinContent(ib, diff)

            if sys.version_info[0] >= 3:
                eval("print(f'{ptlow:>2.0f} < pT < {pthi:>2.0f} loose: {pre:>8.3g} , nominal: {nom:>8.3g}, \
    difference: {pre / nom * 100 - 100 :>8.1f}%, nominal stat: {nomerr / nom * 100:>5.2g}%, \
    loose stat: {preerr / pre * 100:>5.2g}%')")

    print(cut)
    print('raw yield')
    comp(raws)
    # # print('raw yield (nominal MC)')
    # # comp(rawsmc)
    # print('efficiency')
    # comp(effs, True)
    print('inverse efficiency')
    comp(effs)
    print('corrected yield')
    # comp(yds, update=True)
    comp(yds, update=False)
    fio.cd()
    hMCData.Write()
    hPre.Write()
    fio.Write()
    fio.Close()

    def table(ylist, elist, clist):
        binlist = range(1, raws[0].GetNbinsX() + 1)
        ynom = np.array([ ylist[0].GetBinContent(ib) for ib in binlist])
        ypre = np.array([ ylist[1].GetBinContent(ib) for ib in binlist])
        enom = np.array([ elist[0].GetBinContent(ib) for ib in binlist])
        epre = np.array([ elist[1].GetBinContent(ib) for ib in binlist])
        cnom = np.array([ clist[0].GetBinContent(ib) for ib in binlist])
        cpre = np.array([ clist[1].GetBinContent(ib) for ib in binlist])
        ynom_err = np.array([ ylist[0].GetBinError(ib) for ib in binlist])
        ypre_err = np.array([ ylist[0].GetBinError(ib) for ib in binlist])
        enom_err = np.array([ elist[0].GetBinError(ib) for ib in binlist])
        epre_err = np.array([ elist[0].GetBinError(ib) for ib in binlist])
        # cnom_err = np.array([ clist[0].GetBinError(ib) for ib in binlist])
        # cpre_err = np.array([ clist[0].GetBinError(ib) for ib in binlist])
        ynom_err /= ynom
        ypre_err /= ypre
        enom_err /= enom
        epre_err /= epre
        # cnom_err /= cnom
        # cpre_err /= cpre

        title = ['loose raw yield', 'error',
                 'nominal raw yield', 'error',
                 'loose inverse eff.', 'error',
                 'nominal inverse eff.', 'error',
                 'loose corrected yield',
                 'nominal corrected yield',
                 'difference'
                 ]

        numlist = [ypre, ypre_err,
                   ynom, ynom_err,
                   epre, epre_err,
                   enom, enom_err,
                   cpre, cnom, (cpre - cnom) / cnom
                   ]
        percent = [0, 1,
                   0, 1,
                   0, 1,
                   0, 1,
                   0, 0, 1]
        for ti, nu, per in zip(title, numlist, percent):
            linestr = ['{s:>24}'.format(s=ti)] + ['{n:>7.2f}'.format(n = x) for x in nu]
            if ti.rfind('corrected') > 0:
                linestr[1:] = ['{n:>7.0f}'.format(n = x) for x in nu]
            if per:
                linestr = ['{s:>24}'.format(s=ti)] + \
                    ['{n:>5.1f}\\%'.format(n = x * 100) for x in nu]
            line = ' & '.join(linestr) +  ' \\\\'
            print(line)
    table(raws, effs, yds)





    # fineeffs[0].Draw('color')
    # fineeffs[0].SetTitle('1/#epsilon_{nominal}; pT; y')
    # can.SaveAs('efftable_nom.png')
    # fineeffs[1].Draw('color')
    # fineeffs[1].SetTitle('1/#epsilon_{loose}; pT; y')
    # can.SaveAs('efftable_pre.png')
    # finediff = fineeffs[0].Clone()
    # finediff.Add(fineeffs[1], -1)
    # finediff.Draw('colorz')
    # finediff.SetTitle('1/#epsilon_{nominal} - 1/#epsilon_{loose}; pT; y')
    # can.Update()
    # can.SaveAs('efftablediff.png')

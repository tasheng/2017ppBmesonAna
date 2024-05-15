import ROOT as r
import numpy as np

trkErrorFile = "../../../syst_track_sel.root"

BPerrorFile = "../../../2DMapSyst/OutFiles/BPError2D.root";
BPpdfErrorFile = "../../../bp_pdf.root";
BPoutFile = "BPsyst.org";

BserrorFile = "../../../2DMapSyst/OutFiles/BsError2D.root";
BspdfErrorFile = "../../../bs_pdf.root";
BsoutFile = "Bssyst.org";

sep = ' | '
# suf = r' \\'
suf = sep
pre = sep

pp_tracking_error = 2.4

def get_line(name, array):
    return pre + sep.join([name] + ['{:.2g}'.format(e) for e in array]) + suf + '\n'
def make_table(meson, errf, pdff, tracking_error, outFile):
    fError = r.TFile(errf);
    fPdfError = r.TFile(pdff);
    fTrkError = r.TFile(trkErrorFile);

    TnPSyst = fError.Get("TnPSyst");
    BptSyst = fError.Get("BptSyst");
    MCDataSyst = fError.Get("MCDataSyst");
    BDTSyst = fError.Get("BDTSyst");
    PreSyst = fError.Get("PreSelectionSyst");
    nBins = TnPSyst.GetNbinsX()

    pdfSyst = fPdfError.Get(meson + "_error");
    trkSyst = fTrkError.Get(meson + "_track_sel_error");
    bkgeffSyst = fTrkError.Get(meson + "_bkgeff_error");

    tnp = np.array([TnPSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    bpt = np.array([BptSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    mcdata = np.array([MCDataSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    tracking = np.full(nBins, tracking_error)
    pdf = np.array(pdfSyst.GetY())
    trk_sel = np.array(trkSyst.GetY())
    bkgeff = np.array(bkgeffSyst.GetY())

    total = np.sqrt(np.sum([e**2 for e in [bpt, tnp, mcdata, tracking, trk_sel, bkgeff, pdf]], axis=0))

    bdt = np.array([BDTSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    pre = np.array([PreSyst.GetBinContent(i) for i in range(1, nBins + 1)])

    with open(outFile, "w") as f:
        f.write(get_line("Hadron tracking efficiency", tracking))
        f.write(get_line("Track selection", trk_sel))
        f.write(get_line("Data-MC discrepancy", mcdata))
        f.write(get_line(r"\pt shape", bpt))
        f.write(get_line("PDF variation", pdf))
        f.write(get_line("Muon efficiency", tnp))
        f.write(get_line("Data averaging of efficiency", bkgeff))
        f.write(get_line("Sum", total))
        f.write('\n\n')
        f.write(get_line("BDT syst.", bdt))
        f.write(get_line("Pre-selection syst.", pre))
        f.write(get_line("Data--MC syst.", mcdata))

make_table("bp", BPerrorFile, BPpdfErrorFile, pp_tracking_error, BPoutFile)
make_table("bs", BserrorFile, BspdfErrorFile, 2 * pp_tracking_error, BsoutFile)

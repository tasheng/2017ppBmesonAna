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
    return pre + sep.join([name] + ['{:.3f}'.format(e) for e in array]) + suf + '\n'
def make_table(meson, errf, pdff, tracking_error, outFile):
    fError = r.TFile(errf);
    fPdfError = r.TFile(pdff);
    fTrkError = r.TFile(trkErrorFile);

    TnPSyst = fError.Get("TnPSyst");
    BptSyst = fError.Get("BptSyst");
    MCDataSyst = fError.Get("MCDataSyst");
    nBins = TnPSyst.GetNbinsX()

    pdfSyst = fPdfError.Get(meson + "_error");
    trkSyst = fTrkError.Get(meson + "_track_sel_error");

    tnp = np.array([TnPSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    bpt = np.array([BptSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    bdt = np.array([MCDataSyst.GetBinContent(i) for i in range(1, nBins + 1)])
    tracking = np.full(nBins, tracking_error)
    pdf = np.array(pdfSyst.GetY())
    trk_sel = np.array(trkSyst.GetY())

    total = np.sqrt(np.sum([e**2 for e in [bpt, tnp, bdt, tracking, trk_sel, pdf]], axis=0))

    with open(outFile, "w") as f:
        f.write(get_line("Hadron tracking efficiency", tracking))
        f.write(get_line("Track selection", trk_sel))
        f.write(get_line("Data-MC discrepancy", bdt))
        f.write(get_line(r"\pt shape", bpt))
        f.write(get_line("PDF variation", pdf))
        f.write(get_line("TnP systematics", tnp))
        f.write(get_line("Sum", total))

make_table("bp", BPerrorFile, BPpdfErrorFile, pp_tracking_error, BPoutFile)
make_table("bs", BserrorFile, BspdfErrorFile, 2 * pp_tracking_error, BsoutFile)

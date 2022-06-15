# list of produced files

# from subprocess import run, PIPE, Popen
import os
from time import ctime
import re
import ROOT as r
import numpy as np

# examples
# p1 = subprocess.Popen(['echo', link], stdout=PIPE)
# p2 = run(['xclip', '-i', '-selection', 'clipboard'], input=link)

bp_list = [
    # skim output
    "CutSkim/BPMC.root",
    "CutSkim/BPData.root",
    # fit output
    "BP/RawYieldFits/ROOTfiles/yields_Bp_binned_pt.root",
    # cross section
]

bp_pdf_list = [
    "BP/RawYieldFits/signal_systematics_table_Bpt_ntKp.tex",
    "BP/RawYieldFits/background_systematics_table_Bpt_ntKp.tex",
    # "BP/RawYieldFits/general_systematics_table_Bpt_ntKp.tex",
    ]

bp_list.extend(bp_pdf_list)

bs_pdf_list = [
    "Bs/RawYieldFits/signal_systematics_table_Bpt_ntphi.tex",
    "Bs/RawYieldFits/background_systematics_table_Bpt_ntphi.tex",
    # "BP/RawYieldFits/general_systematics_table_Bpt_ntKp.tex",
    ]

def print_time():
    """print the mtime of all the files"""
    for file in bp_list:
        # p1 = Popen(['ls', '-lh', file])
        time = os.path.getmtime(file)
        print(ctime(time), file)


def get_pdf_syst(inFileList, outfile, hname, nbins):
    """
    get the biggest error from signal and background pdf variation
    return the root mean squared error
    """
    error_sig_bkg = np.zeros([2, nbins])
    nshape = [2, 2]
    for itype, file in enumerate(inFileList):
        with open(file) as fin:
            lines = fin.readlines()[2:]
        xgroups = re.findall('\d+\D+\d+', lines[0])
        xs = [np.mean([float(num) for num in re.findall('\d+', g)]) for g in xgroups]
        err = np.zeros([nshape[itype], nbins])
        for ishape in range(nshape[itype]):
            err[ishape] = np.array(re.findall('\d+\.\d*', lines[ishape + 1]), 'd')
        bigger_error = np.amax(err, axis=0)
        error_sig_bkg[itype] = bigger_error
    error_sum = np.sqrt(np.sum(error_sig_bkg ** 2, axis=0))
    # write output
    fout = r.TFile(outfile, "recreate")
    h = r.TGraph(len(xs), np.array(xs), error_sum)
    h.SetName(hname)
    h.Write()
    print(np.array(h.GetY()))
    fout.Write()
    fout.Close()
    return

get_pdf_syst(bp_pdf_list, "bp_pdf.root", "bp_error", 7)
get_pdf_syst(bs_pdf_list, "bs_pdf.root", "bs_error", 4)

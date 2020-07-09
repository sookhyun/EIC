from __future__ import absolute_import

import os
import sys
import ROOT

__all__ = [
    'eicsmear'
]



def eicsmear(infile, outfile):

    curdir = os.getcwd()
    ROOT.gSystem.Load("libeicsmear")
    ROOT.gROOT.ProcessLine(".L " + curdir + "/smear/SmearHandBook.cxx")
    ROOT.gROOT.ProcessLine("SmearTree(BuildHandBookDetector(), \""+infile+"\", \""+outfile+"\")")
    return 1;

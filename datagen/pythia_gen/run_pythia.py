from __future__ import absolute_import

import os
import sys

__all__ = [
    'run_pythiaeRHIC'
]

def run_pythiaeRHIC(ofile, e_nucleus, e_lepton, min_Q2, max_Q2, min_x,max_x, min_y,max_y,nevents):

    # defualt settings
    curdir =os.getcwd()
    dfile = open(curdir+"/pythia_gen/pythia_eRHIC_Hera.cfg", "r")
    default_cfg = dfile.read()
    
    # header
    output = ofile + " ! output file name\n"
    pid_lepton = "11 ! lepton beam type\n"
    energies = e_nucleus + "," + e_lepton + " ! proton and electron beam energy\n"
    events = nevents + ",1 ! number of events, number of events to print to stdout\n"

    # write a config file
    cfile = open("pythia.cfg","w")
    cfile.write(output)
    cfile.write(pid_lepton)
    cfile.write(energies)
    cfile.write(events)
    
    # add the x, q, and y, which we might want to add as arguments later
    if min_x == '0':
        print("Using default x_min = 1e-09")
        min_x='1e-09'
    if max_x == '0':
	print("Using default x_max = 0.99")
        max_x='0.99'
    cfile.write(min_x + "," + max_x + "      ! xmin and xmax\n")
    if min_y == '0':
        min_y='1e-04'
    if max_y == '0':
        max_y='1.00'
    cfile.write(min_y + "," + max_y + "      ! ymin and ymax\n")

    if max_Q2 == '0':
        max_Q2 = '20000'
    cfile.write(min_Q2 +","+max_Q2+"       ! Q2min and Q2max\n")
    
    # add all the other parameters
    cfile.write(default_cfg)
    cfile.close()

    # run Pythia
    executable = "$EICDIRECTORY/bin/pythiaeRHIC < pythia.cfg"
    os.system(executable)

    #convert an ASCII to a ROOT file
    import ROOT
    ROOT.gSystem.Load("libeicsmear")
    ROOT.BuildTree(ofile, ".", -1)

    # clean up
    #os.system("rm pythia.cfg")
    #os.system("rm " + ofile)
    return 1

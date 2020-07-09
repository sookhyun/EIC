import os
import sys
import subprocess

if len(sys.argv)<6:
    print("Enter arguments w/o ',' : 'truename' 'freconame' E_proton E_lepton min_Q2 max_Q2 x_min x_max y_min y_max nevents ")
    print("EIC : python run_sim.py 'ftrue.root' 'freco.root' 275 18 2500 2500 0 0 0.2518 0.2522 500  ")
    print("H1 : python run_sim.py 'ftrue.root' 'freco.root' 820 27.5 6399 6401 0 0 0.35 0.36 1000  ")
#								      Q2        x     y	   
# above is to compare with 1303.6952, 
# where sqrt(s)=300 GeV, E_p= 820 GeV, E_e = 27.5, x=0.2, Q=80GeV, y=6400/0.2/90000,
# other option is x=0.05, Q=50GeV, y= 2500/0.1/90000=25/90=0.2777

# if you fix y and Q2 then x is determined by xys=Q2, x is determined by y, not the other way around!
# Q2>25 to remove photo production
# for EIC: 
# sqrt(s)= 140.716 GeV, E_p=275, E_e=18, x=0.5, Q=50GeV, y=2500/0.5/19801=0.252
# python run_sim.py 'ftrue.root' 'freco.root' 275 18 2500 2500 0 0 0.2518 0.2522 1000 

 
# index starts from 1. 
ftrue = sys.argv[1]
freco = sys.argv[2]
e_nucl = sys.argv[3]
e_lep = sys.argv[4]
min_Q2 = sys.argv[5]
max_Q2 = sys.argv[6]
min_x = sys.argv[7]
max_x = sys.argv[8]
min_y = sys.argv[9]
max_y = sys.argv[10]
nevents = sys.argv[11]

from pythia_gen import *
from smear import *

run_pythiaeRHIC(ftrue,
      		e_nucl,
      		e_lep,
      		min_Q2,
		max_Q2,
		min_x,
		max_x,
		min_y,
		max_y,
      		nevents)
eicsmear(ftrue,freco)

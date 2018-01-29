#! /usr/bin/env python
################################################################################
#TEMPLATE: 2D-mesh of the DM particles spectrum from events in LHE-format
#
#input:  the user must provide
#        - pdg code of the DM particle
#        - lhe input file (default: unweighted_events.lhe.gz)
#        - efficiency function (modelling the particular experimental setup)
#          (default = 1) 
#         
#output: required input files for the maddump plugin
#        ehist.dat: it contains the table of the energy 1D-histogram
#                   (after integration over angular variable)
#                   it is normalized to 1
#                   format: 4 columns 
#                           E_i (bin startingpoint), E_i+1 (bin endpoint),
#                           phi(E) (bin height), error on phi(E)
#                                                (as given by Poisson statics)  
#        cell_fortran.dat: it contains the fitted 2D-mesh as a list of cells
#                   format 
################################################################################

if '__main__' == __name__:
    import sys
    sys.path.append('../../')
    sys.path.append('../../madgraph/various')

import os
import logging
import math
import lhe_parser
import madgraph.various.banner as banner_mod
import numpy as np
from meshfitter2D import *

pjoin = os.path.join
pwd_path=os.getcwd()
logging.basicConfig()

lpp2 = {'electron' : 0, 'proton' : 1}
ebeam2 = {'electron' : 0.000511, 'proton' : 0.938}

#INPUT
#DM_pdgcode = [5000521,-5000521]
DM_pdgcode = [22006,-22006]
events_lhefile = 'unweighted_events.lhe.gz'
detector_particle = 'proton' 

def eff_function(E,Theta):
    return 1

# Construction of the data sample suited for the 2D histogramming
# by parsing the LHE events file.
# Data is a list of WeightedPoint (defined in the 2DMeshFitter class)
# (E, theta, w(E,theta)), where w is the weight given by the eff_function 
lhe = lhe_parser.EventFile(events_lhefile)
data = []
E_min = 1e20
E_max = 0.
theta_min = np.pi/2.
theta_max = 0.

for event in lhe:
    for particle in event:
        if particle.status == 1: # stable final state 
            if particle.pid in DM_pdgcode:
                p = lhe_parser.FourMomentum(particle)
                theta = np.arccos(p.pz/p.norm)
                data.append(WeightedPoint(p.E,theta,eff_function(p.E,theta)))
                if p.E < E_min:
                    E_min = p.E
                if p.E > E_max:
                    E_max = p.E
                if theta < theta_min:
                    theta_min = theta
                if theta > theta_max:
                    theta_max = theta

print(E_min,theta_min,E_max,theta_max,len(data))
# Generation of the 2DMesh and of the output file cell_fortran.dat
# the mesh is also plotted in a 2D graph 
hist2D = CellHistogram(Point(E_min,theta_min),E_max-E_min,theta_max-theta_min,20)
hist2D.add_pts(data)
hist2D.fit()
hist2D.export_histogram('cell_fortran.dat')
hist2D.plot()

# Generation of the 1D distribution integrated in angles
# and of the output ehist.dat
hist2D.fit('1D_x')
hist2D.export_histogram('ehist.dat','1D_x')

# Update parameters in the run card
run_card = banner_mod.RunCard(pjoin(pwd_path,'run_card_default.dat'))
run_card['lpp2'] = lpp2[detector_particle]
run_card['ebeam1'] = E_max
run_card['ebeam2'] = ebeam2[detector_particle]

run_card.write(pjoin(pwd_path, 'run_card.dat'))

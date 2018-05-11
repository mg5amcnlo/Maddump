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
    sys.path.append('../../../')
    sys.path.append('../../../madgraph/various')

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


def get_proc_characteristics(path):
    proc_file = open(path,'r')
    proc_characteristics = {} 
    for line in proc_file:
        if '=' not in line:
            continue
        else:
            args = line.split()
            proc_characteristics[args[0]] = args[2]
    return proc_characteristics


proc_characteristics = get_proc_characteristics('../SubProcesses/proc_characteristics')

print(proc_characteristics)
DM = int(proc_characteristics['DM'])
DM_pdgcode = [DM,-DM]
events_lhefile = 'unweighted_events.lhe.gz'
detector_particle = 'proton' 


#distance target - detector in cm
d_target_detector = 3804.75

#circular shape
th_min=0.
th_max=.0118

#rectangular shape
#coordinate of the center wrt to the beam axis in cm 
xc = 0.
yc = 0.
#side in cm
x_side=37.45*2
y_side=45.15*2
depth = 321.0
xmin = xc - x_side/2.
xmax = xc + x_side/2.
ymin = yc - y_side/2.
ymax = yc + y_side/2.

ncores=4

def heaviside(x):
    if x>0:
        return 1
    else:
        return 0

def max_travel_distance(theta,cphi,sphi):

    sth = np.sin(theta)
    # cphi = np.cos(phi)
    # sphi = np.sin(phi)
    z1 = d_target_detector
    z2 = z1 + depth

    x1 = z1*sth*cphi
    y1 = z1*sth*sphi
    
    in_z1 = heaviside(x1-xmin)*heaviside(xmax-x1)*heaviside(y1-ymin)*heaviside(ymax-y1)
    if in_z1 == 0.:
        return 0.

    x2 = z2*sth*cphi
    y2 = z2*sth*sphi

    in_z2 = heaviside(x2-xmin)*heaviside(xmax-x2)*heaviside(y2-ymin)*heaviside(ymax-y2)

    if in_z2 != 0.:
        return np.sqrt((x2-x1)**2+(y2-y1)**2+depth**2)

    if x2 > xmax:
        x3 =  xmax
        y3 =  xmax*sphi/cphi
        z3 =  xmax/sth/cphi 
        return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    if x2 < xmin:
        x3 =  xmin          
        y3 =  xmin*sphi/cphi
        z3 =  xmin/sth/cphi 
        return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    if y2 > ymax:
        x3 =  ymax*cphi/sphi
        y3 =  ymax
        z3 =  ymax/sth/sphi
        return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    if y2 < ymin:
        x3 =  ymin*cphi/sphi
        y3 =  ymin
        z3 =  ymin/sth/sphi
        return np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)

def eff_function(E,theta,cphi,sphi):
    return max_travel_distance(theta,cphi,sphi)*np.cos(theta)/depth


# def eff_function(E,theta,cphi,sphi):
#     x = d_target_detector*np.sin(theta)*cphi
#     y = d_target_detector*np.sin(theta)*sphi
#     return heaviside(x-xmin)*heaviside(xmax-x)*heaviside(y-ymin)*heaviside(ymax-y)

# def eff_function(E,theta):
#     return heaviside(theta-th_min)*heaviside(th_max-theta)

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

xsec =1.#lhe.cross
nevt=len(lhe)
npass = 0
for event in lhe:
    for particle in event:
        if particle.status == 1: # stable final state 
            if particle.pid in DM_pdgcode:
                p = lhe_parser.FourMomentum(particle)
                theta = np.arccos(p.pz/p.norm)
                cphi = p.px/np.sqrt(p.px**2+p.py**2)
                sphi = p.py/np.sqrt(p.px**2+p.py**2)
                data.append(WeightedPoint(p.E,theta,xsec/nevt*eff_function(p.E,theta,cphi,sphi)))
                
                if eff_function(p.E,theta,cphi,sphi) != 0:                    
                    npass+=1
                    if p.E < E_min:
                        E_min = p.E
                    if p.E > E_max:
                        E_max = p.E
                    if theta < theta_min:
                        theta_min = theta
                    if theta > theta_max:
                        theta_max = theta
#                 if p.E < E_min:
#                     E_min = p.E
#                 if p.E > E_max:
#                     E_max = p.E
#                 if theta < theta_min:
#                     theta_min = theta
#                 if theta > theta_max:
#                     theta_max = theta

# if(theta_max > th_max): theta_max= th_max
# if(theta_min < th_min): theta_max= th_min

print(E_min,theta_min,E_max,theta_max,len(data))
    
# Generation of the 2DMesh and of the output file cell_fortran.dat
# the mesh is also plotted in a 2D graph 
hist2D = CellHistogram(Point(E_min,theta_min),E_max-E_min,theta_max-theta_min,50)
hist2D.add_pts(data)
hist2D.fit(ncores=ncores)
#hist2D.export_histogram('cell_fortran.dat')
plot('cell_fortran.dat')

# Generation of the 1D distribution integrated in angles
# and of the output ehist.dat
hist2D.fit('1D_x',ncores=ncores)

# Update parameters in the run card
run_card = banner_mod.RunCard(pjoin(pwd_path,'run_card_default.dat'))
run_card['lpp2'] = lpp2[detector_particle]
run_card['ebeam1'] = E_max
run_card['ebeam2'] = ebeam2[detector_particle]

run_card.write(pjoin(pwd_path, 'run_card.dat'))

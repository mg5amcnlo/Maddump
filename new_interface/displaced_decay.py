# if '__main__' == __name__:
#     import sys
#     sys.path.append('../../../../')
#     sys.path.append('../../../../madgraph/various')
#     sys.path.append('../../../../models')

import madgraph.various.lhe_parser as lhe_parser
import models.check_param_card as param_card_mod
import numpy as np
import random as random
import copy as copy

class displaced_decay():

    c = 29979245800 # speed of light in cm/s
    hbar_times_c = 0.19732697e-13 # in GeV*cm
    
    def __init__(self, lhe_path, param_card_path, mother, daughter):
        self.lhe_input = lhe_parser.EventFile(lhe_path)
        self.total_events = 0.
        self.pdg_code = mother
        param_card = param_card_mod.ParamCard(param_card_path)
        self.mass = param_card['mass'].get(self.pdg_code).value
        self.width = param_card['decay'].get(self.pdg_code).value
        for BR in param_card['decay'].decay_table[self.pdg_code]:
            if abs(BR.lhacode[1]) == abs(daughter):
                self.BR = BR.value
        if not self.BR:
            raise Exception, 'No decay in %d, check the particle id and the couplings!' % daughter
        self.norm = 2e20

    def finalize_output(self, path):
        #out_file = open(path,'w')
        out = lhe_parser.EventFile(path,mode='w')
        nLLP = 0 
        for event in self.lhe_input:
            for particle in event:
                if particle.pid == self.pdg_code:
                    nLLP += 1
                    p = lhe_parser.FourMomentum(particle)
                    weight,displ = self.compute_weight_displacement(p)
                    if weight == 0.:
                        particle.vtim = -1.
                    else:
                        particle.vtim = displ
                        self.total_events += weight
            out.write_events(event)

        # divided by  the total number ntot = len(out) of events 
        # to take into account meson multiplicity!
        self.total_events = self.total_events/len(out) * self.norm
        
                    
    def compute_weight_displacement(self,p):
        ctheta = p.pz/p.norm
        if ctheta < 0 :
            return [0.,0.]
        if ctheta != 1.:
            cphi = p.px/(p.norm**2-p.pz**2)
            sphi = p.py/(p.norm**2-p.pz**2)
        else:
            cphi = 0.
            sphi = 0.
        theta = np.arccos(ctheta)
        l1,l2 = self.travel_distance(theta,cphi,sphi)
        if l2 < 0.:
            return [0.,0.]
        beta= p.norm/p.E
        gamma = p.E/self.mass
        Lambda =  gamma*beta*self.hbar_times_c/self.width
        #print(Lambda)
        if Lambda == 0.:
            print(gamma,beta)
        #print(gamma,Lambda)
        weight = np.exp(-l1/Lambda) - np.exp(-l2/Lambda) *self.BR
        displacement = -Lambda*np.log( np.exp(-l2/Lambda) + (np.exp(-l1/Lambda)-np.exp(-l2/Lambda))*random.random() )
        return weight,displacement
    
    def travel_distance(self,theta,cphi,sphi):
        d = 7000 # in cm
        delta = 5500 # in cm
        lmin = d/np.cos(theta)
        lmax = (d+delta)/np.cos(theta)
        return lmin,lmax


    

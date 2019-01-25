# if '__main__' == __name__:
#     import sys
#     sys.path.append('../../../../')
#     sys.path.append('../../../../madgraph/various')
#     sys.path.append('../../../../models')


import models.check_param_card as param_card_mod
import numpy as np
import random
import copy

import madgraph.various.lhe_parser as lhe_parser
import madgraph.various.misc as misc

class displaced_decay():

    c = 29979245800 # speed of light in cm/s
    hbar_times_c = 0.19732697e-13 # in GeV*cm
    def __init__(self, lhe_path, param_card_path, mother):
        self.lhe_input = lhe_parser.EventFile(lhe_path)
        self.lhe_input.allow_empty_event = True
        self.nevts = len(self.lhe_input)
        self.total_events = 0.
        self.pdg_code = mother
        self.param_card = param_card_mod.ParamCard(param_card_path)
        self.mass = self.param_card['mass'].get(self.pdg_code).value
        self.width = self.param_card['decay'].get(self.pdg_code).value
        self.decay_list=[]
        for BR in self.param_card['decay'].decay_table[self.pdg_code]:
            pdg_BR=[BR.lhacode[i] for i in range(1,BR.lhacode[0]+1)]            
            self.decay_list.append([pdg_BR,False,BR.value])

        #print(self.decay_list)
        
        #jet = [1,2,3,4,5,21]
        # for BR in param_card['decay'].decay_table[self.pdg_code]:
        #     if abs(BR.lhacode[1]) == abs(daughter):
        #         if abs(BR.lhacode[2]) in jet:
        #             if abs(BR.lhacode[3]) in jet:
        #                 self.BR += BR.value
        # if not self.BR:
        #     raise Exception, 'No decay in %d, check the particle id and the couplings!' % daughter
        self.norm = 2e20

    def finalize_output(self, path):
        #out_file = open(path,'w')
        out = lhe_parser.EventFile(path,mode='w')
        for event in self.lhe_input:
            self.update_BR(event,self.pdg_code)
            for particle in event:
                if particle.pid == self.pdg_code:
                    p = lhe_parser.FourMomentum(particle)
                    weight,displ = self.compute_weight_displacement(p)
                    weight= weight
                    if weight == 0.:
                        particle.vtim = -1.
                    else:
                        particle.vtim = displ
                        self.total_events += weight
            out.write_events(event)
        BRfac = 0
        for k in range(len(self.decay_list)):
            if self.decay_list[k][1] == True:
                BRfac += self.decay_list[k][2] 
        # divided by  the total number ntot = self.nevts of events 
        # to take into account meson multiplicity!
        #print(BRfac)
        self.total_events = self.total_events/self.nevts * self.norm * BRfac

        
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
        if Lambda == 0.:
            print(gamma,beta)
        weight = np.exp(-l1/Lambda) - np.exp(-l2/Lambda)
        displacement = -Lambda*np.log( np.exp(-l2/Lambda) + (np.exp(-l1/Lambda)-np.exp(-l2/Lambda))*random.random() )
        return weight,displacement

    
    # def get_BR(self,event,pdgcode):

    #     decay_evt=event.get_decay(pdg_code=pdgcode)
    #     print(decay_evt)
    #     ndaughters=len(decay_evt)-1
    #     pdg_daughters=[]
    #     for particle in decay_evt:
    #         if particle.status == 1:
    #             pdg_daughters.append(particle.pid)

    #     found = False


    #     for BR in self.param_card['decay'].decay_table[pdgcode]:
    #         if BR.lhacode[0] == ndaughters:
    #             pdg_BR=[BR.lhacode[i] for i in range(1,ndaughters+1)]
    #             if sorted(pdg_BR) == sorted(pdg_daughters):
    #                 x=BR.value
    #                 found=True
    #                 break
    
    #     if not found:
    #         for BR in self.param_card['decay'].decay_table[pdgcode]:
    #             pdg_daughters_opposite = [-pdg_daughters[i] for i in range(ndaughters)]
    #             if BR.lhacode[0] == ndaughters:            
    #                 pdg_BR=[BR.lhacode[i] for i in range(1,ndaughters+1)]
    #                 if sorted(pdg_BR) == sorted(pdg_daughters_opposite):
    #                     x=BR.value
    #                     found=True
    #                     break
                    
    #     if not found:
    #         raise Exception, 'Error in getting the BR for the current event, check the BRs have been computed correclty' 

    #     return x

    def update_BR(self,event,pdgcode):

        pdg_daughters=[]

        for i in range(len(event)):
            if (event[i].pid==pdgcode) and (event[i].status==2):
                j=1 
                while(event[i+j].status==1):
                    pdg_daughters.append(event[i+j].pid)
                    j+=1
                    if(i+j==len(event)):
                        break
                        
                pdg_daughters_opposite = [-pdg_daughters[i] for i in range(len(pdg_daughters))]
                for k in range(len(self.decay_list)):
                    if (sorted(pdg_daughters) == sorted(self.decay_list[k][0]) or sorted(pdg_daughters_opposite) == sorted(self.decay_list[k][0])):
                        self.decay_list[k][1] = True

                        
    def travel_distance(self,theta,cphi,sphi):
        d = 7000 # in cm
        delta = 5500 # in cm
        lmin = d/np.cos(theta)
        lmax = (d+delta)/np.cos(theta)
        return lmin,lmax

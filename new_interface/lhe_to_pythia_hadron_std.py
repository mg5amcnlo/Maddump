# if '__main__' == __name__:
#     import sys
#     sys.path.append('../../../../')
#     sys.path.append('../../../../madgraph/various')

import madgraph.various.lhe_parser as lhe_parser
import numpy as np
import random as random
import copy as copy

class LHEtoPYTHIAHadronSTD():

    quark_conv_table = { 1: 1,
                         2: 2,
                         3: 1,
                         4: 2,
                         -1: 1,
                         -2: 2,
                         -3: 1,
                         -4: 2,}

    mass = { 'proton':  0.938,
             'neutron': 0.9396}
    
    quark_list = [1,2,3,4,-1,-2,-3,-4]

    
    def __init__(self, path, target='proton'):
        self.lhe_input = lhe_parser.EventFile(path)
        self.target = target
        #self.target = lhe_input[1].pdg        

        
    def write_PYTHIA_input(self, path):
        #out_file = open(path,'w')
        out = lhe_parser.EventFile(path,mode='w')
        out.write('<LesHouchesEvents version="3.0">\n')
        out.write('<init></init>\n')
        for event in self.lhe_input:
            hadronic_event = self.parton_level_reconstruction(event)  
            # for particle in hadronic_event:
            #     print(particle.pid, particle.E,particle.helicity)
            # exit(1)
            out.write_events(hadronic_event)
            
    def parton_level_reconstruction(self,event):
        event_to_hadronize = lhe_parser.Event()
        event_to_hadronize.tag = event.tag        
        diquark = lhe_parser.Particle(event=event_to_hadronize)
        quark = lhe_parser.Particle(event=event_to_hadronize)
        for particle in event:
            if particle.pid in self.quark_list:
                if particle.status ==-1:
                    diquark.pid = self.assign_diquark(self.quark_conv_table[particle.pid])
                    diquark.status = 1
                    diquark.mother1= 0
                    diquark.mother2= 0                    
                    diquark.px = -particle.px
                    diquark.py = -particle.py
                    diquark.pz = -particle.pz
                    diquark.E = self.mass[self.target]-particle.px
                    p = lhe_parser.FourMomentum(diquark)
                    diquark.mass = p.mass
                    diquark.vtim = 0
                    diquark.helicity = 9
                elif particle.status == 1:
                    quark = copy.deepcopy(particle)
                    quark.pid = self.quark_conv_table[particle.pid]
                    quark.mother1= 0
                    quark.mother2= 0
                    if quark.color1 == 0:
                        quark.color1 = quark.color2
                        quark.color2 = 0
                        diquark.color1 = 0
                        diquark.color2 = quark.color1
                    else:
                        diquark.color1 = 0
                        diquark.color2 = quark.color1                        
        event_to_hadronize.nexternal = 2
        event_to_hadronize.append(diquark)
        event_to_hadronize.append(quark)
        return event_to_hadronize
        
        
                    
    def assign_diquark(self,quark):
        if self.target == 'proton':
            if quark == 1:   # d quark
                return 2203  # -> (uu)_1
            elif quark == 2: # u quark
                if random.random() < 0.75: 
                    return 2101 # -> (ud)_0 75%
                else: 
                    return 2103 # -> (ud)_1 25%
            else:
                print('Error, quark do not allowed!')
                exit(-1)
        else:
            if quark == 2:   # d quark
                return 1103  # -> (uu)_1
            elif quark == 1: # u quark
                if random.random() < 0.75: 
                    return 2101 # -> (ud)_0 75%
                else: 
                    return 2103 # -> (ud)_1 25%
            else:
                print('Error, quark do not allowed!')
                exit(-1)

# if '__main__' == __name__:
#     import sys
#     sys.path.append('../../../../')
#     sys.path.append('../../../../madgraph/various')
#     sys.path.append('../../../../models')


from __future__ import absolute_import
from __future__ import print_function
import models.check_param_card as param_card_mod
import numpy as np
import random
import copy

import madgraph.various.lhe_parser as lhe_parser
import madgraph.various.misc as misc
from six.moves import range

class displaced_decay():

    c = 29979245800 # speed of light in cm/s
    hbar_times_c = 0.19732697e-13 # in GeV*cm
    def __init__(self, lhe_path, param_card_path, mothers):
        self.lhe_input = lhe_parser.EventFile(lhe_path)
        self.lhe_input.allow_empty_event = True
        self.nevts = len(self.lhe_input)
        self.total_events = 0.
        self.pdg_code=[]
        self.param_card = param_card_mod.ParamCard(param_card_path)
        self.mass = []
        self.width = []
        self.decay_list = []
        self.vertices = []
        self.index_particle = 0
        
        for mother in mothers:
            self.pdg_code.append(mother)
            self.mass.append(self.param_card['mass'].get(mother).value)
            self.width.append(self.param_card['decay'].get(mother).value)
            self.decay_list.append([])
            for BR in self.param_card['decay'].decay_table[mother]:
                pdg_BR=[BR.lhacode[i] for i in range(1,BR.lhacode[0]+1)]            
                self.decay_list[-1].append([pdg_BR,False,BR.value])

            #print(self.decay_list)
        
        #jet = [1,2,3,4,5,21]
        # for BR in param_card['decay'].decay_table[self.pdg_code]:
        #     if abs(BR.lhacode[1]) == abs(daughter):
        #         if abs(BR.lhacode[2]) in jet:
        #             if abs(BR.lhacode[3]) in jet:
        #                 self.BR += BR.value
        # if not self.BR:
        #      raise Exception, 'No decay in %d, check the particle id and the couplings!' % daughter
        self.norm = 3.e16

    def finalize_output(self, path):
        #out_file = open(path,'w')
        out = lhe_parser.EventFile(path,mode='w')
        # print BRs on top of the evts file
        out.write("<BRs>\n")
        for imoth in range(len(self.pdg_code)):
            # write all the BRS on top of the final file
            out.write("#  pdg_code: "+ str(self.pdg_code[imoth]))
            out.write("#  BR \t \t NDA \t ID1 \t ID2 \t ... \n")
            for decay in self.decay_list[imoth]:
                out.write('   ' + str(decay[2])+'\t'+str(len(decay[0]))+'\t')
                for id in decay[0]:
                    out.write(str(id)+'\t')
                out.write('\n')
        out.write("</BRs>\n")
        ia=0
        for event in self.lhe_input:
            weight_list = []
            for imoth in range(len(self.pdg_code)):
                self.index_part = ia
                weight_list.append([])
                # update the BRs used
                self.update_BR(event,imoth)
                for particle in event:
                    if particle.pid == self.pdg_code[imoth]:
                        p = lhe_parser.FourMomentum(particle)
                        weight,displ = self.compute_weight_displacement(p,imoth)
                        weight_list[imoth].append(weight)
                        if weight > 0.:
                            particle.vtim = displ
                            # self.total_events += weight
                        else:
                            particle.vtim = -1.
                        self.index_part +=1
            ia=self.index_part
            out.write_events(event)
            out.write('<wgt>\n')
            for j in range(len(weight_list[0])):
                weight = 1.
                for imoth in range(len(weight_list)):
                    weight *= weight_list[imoth][j]
                self.total_events += weight
                out.write(str(weight/self.nevts)+'\n')
            out.write('<\wgt>\n')

            
        BRfac = 1
        for imoth in range(len(self.pdg_code)):
            BRfac_imoth = 0
            for k in range(len(self.decay_list[imoth])):
                #print(self.decay_list[imoth][k][1],self.decay_list[imoth][k][2])
                if self.decay_list[imoth][k][1] == True:
                    BRfac_imoth += self.decay_list[imoth][k][2]
            BRfac *= BRfac_imoth
                    
        # divided by  the total number ntot = self.nevts of events 
        # to take into account meson multiplicity!
        # print('++++++++++++++++++++')
        # print(self.vertices[0])
        # print(self.vertices[1])        
        # print(BRfac)
        self.total_events = self.total_events/self.nevts * self.norm * BRfac

        
    def compute_weight_displacement(self,p,imoth):
        if(self.index_part==0):
            self.vertices.append([])            
        ctheta = p.pz/p.norm
        # print p.px,p.py,p.pz,p.E,ctheta
        # print ctheta
        if imoth==0 and ctheta < 0 :
            self.vertices[imoth].append([-1.,-1.,-1.])
            return [0.,0.]
        if ctheta != 1.:
            cphi = p.px/np.sqrt((p.norm**2-p.pz**2))
            sphi = p.py/np.sqrt((p.norm**2-p.pz**2))
        else:
            cphi = 0.
            sphi = 0.
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Theta $$$$$$$$$$$$$$$$$$$$$$ = 0 ") 
        theta = np.arccos(ctheta)
        l1,l2 = self.travel_distance(p,theta,cphi,sphi,imoth)

        
        if l2 < 0.:
            self.vertices[imoth].append([-1.,-1.,-1.])
            return [0.,0.]
        beta= p.norm/p.E
        gamma = p.E/self.mass[imoth]
        Lambda =  gamma*beta*self.hbar_times_c/self.width[imoth]
        if Lambda == 0.:
            print((gamma,beta))
        if (l1/Lambda>300.):
            self.vertices[imoth].append([-1.,-1.,-1.])            
            return 0.,-1.
        else:
            weight = np.exp(-l1/Lambda) - np.exp(-l2/Lambda)
            # displacement = -Lambda*np.log( np.exp(-l2/Lambda) + (np.exp(-l1/Lambda)-np.exp(-l2/Lambda))*random.random() )
            displacement = -Lambda*np.log( np.exp(-l1/Lambda) - (np.exp(-l1/Lambda)-np.exp(-l2/Lambda))*random.random() )

        if (imoth==0):
            e1 = np.sin(theta)*cphi
            e2 = np.sin(theta)*sphi
            e3 = np.cos(theta)
            self.vertices[imoth].append([e1*displacement,e2*displacement,e3*displacement])
        else:
            e1 = np.sin(theta)*cphi
            e2 = np.sin(theta)*sphi
            e3 = np.cos(theta)
            x0 = self.vertices[imoth-1][self.index_part][0]
            y0 = self.vertices[imoth-1][self.index_part][1]
            z0 = self.vertices[imoth-1][self.index_part][2]

            x=x0+e1*displacement
            y=y0+e2*displacement
            z=z0+e3*displacement
            
            self.vertices[imoth].append([x,y,z])
            
        # print "+++++++++++++++++++++++++++++", theta, l1, l2, displacement, beta, gamma, Lambda, p.norm, p.E, self.mass, self.width
        
        
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

    def update_BR(self,event,imoth):
        
        pdgcode=self.pdg_code[imoth]
        pdg_daughters=[]

        for i in range(len(event)):
            if (event[i].pid==pdgcode) and (event[i].status==2):
                start = False
                end = False
                j=1
                while(event[i+j].status==1 or event[i+j].status==2):
                    #print i,event[i+j].pid,event[i+j].mother1.pid
                    if(event[i+j].mother1.pid == pdgcode):
                        if start == False:
                            start = True
                        pdg_daughters.append(event[i+j].pid)
                    elif (start==True):
                        end = True
                    if(end == True or i+j+1==len(event)):
                        break
                    if(event[i+j].status==2):
                        break
                    j+=1
                #print(pdg_daughters)
                pdg_daughters_opposite = [-pdg_daughters[i] for i in range(len(pdg_daughters))]
                for k in range(len(self.decay_list[imoth])):
                    if (sorted(pdg_daughters) == sorted(self.decay_list[imoth][k][0]) or sorted(pdg_daughters_opposite) == sorted(self.decay_list[imoth][k][0])):
                        self.decay_list[imoth][k][1] = True

                        
    def travel_distance(self,p,theta,cphi,sphi,imoth):
        d = 10000 # in cm
        delta = 12000 # in cm
        theta_max = 0.01

        # print('*******')
        # print(imoth,self.index_part)

        if(imoth==0):
            if theta > theta_max:
                return 0.,-1 
            lmin = d/np.cos(theta)
            lmax = (d+delta)/np.cos(theta)
            return lmin, lmax
        
        lmin = 0.
        e1 = np.sin(theta)*cphi
        e2 = np.sin(theta)*sphi
        e3 = np.cos(theta)
        x0 = self.vertices[imoth-1][self.index_part][0]
        y0 = self.vertices[imoth-1][self.index_part][1]
        z0 = self.vertices[imoth-1][self.index_part][2]

        # print(x0,y0,z0)

        if z0<0:
            return 0.,-1.

        k = np.tan(theta_max)**2
        
        if e3<0:
            print('Should not be here: child particle is not in the forward direction! Case not implemented. Skip decay.')
            return 0.,-1.

        z=d+delta
        t = (z-z0)/e3
        x=x0+e1*t
        y=y0+e2*t


        if (x**2+y**2 < k*z**2):
            lmax = np.sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2)
            return lmin,lmax

        #t= -((e1*x0 + e2*y0 - k*e3*z0 + np.sqrt(4*(e1*x0 + e2*y0 - k*e3*z0)**2 - 4*(e1**2 + e2**2 - k*e3**2)*(x0**2 + y0**2 - k*z0**2))/2.)/(e1**2 + e2**2 - k*e3**2))
        t=(-(e1*x0) - e2*y0 + k*e3*z0 + np.sqrt(4*(e1*x0 + e2*y0 - k*e3*z0)**2 - 4*(e1**2 + e2**2 - k*e3**2)*(x0**2 + y0**2 - k*z0**2))/2.)/(e1**2 + e2**2 - k*e3**2)

        x=x0+e1*t
        y=y0+e2*t
        z=z0+e3*t                

        # print(e1,e2,e3,k)
        # print(t,z)

        lmax = np.sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2)
            
        return lmin,lmax

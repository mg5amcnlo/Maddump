if '__main__' == __name__:
    import sys
    sys.path.append('../../')
    sys.path.append('../../madgraph/various')

import numpy as np
import cmath
import os
import matplotlib.pyplot as plt
from meshfitter2D import *
import multiprocessing
import fit2D_card as fit2D
import ebeampdffit_card as ebeamfit
import models.check_param_card as param_card_mod

class fit2D_ebeampdf(CellHistogram):

    _ebeamdict = {'electron_pdf': 11, 'positron_pdf': -11, 'gamma_pdf': 22} 
    
    def __init__(self,input_lhe_evts,label):
        
        self.ebeamfit_card = ebeamfit.EbeamPdfFitCard(pjoin('../Cards/ebeampdf_fit_card.dat'))
        self.fit2D_card = fit2D.Fit2DCard(pjoin('../Cards/fit2D_card.dat'))
        self.param_card = param_card_mod.ParamCard(pjoin('../Cards/param_card.dat'))

        self.label = label

        self.E_min,self.E_max,self.theta_min,self.theta_max,self.data = \
                    self.get_data(input_lhe_evts)
        
        super(fit2D_ebeampdf,self).__init__(Point(self.E_min,self.theta_min), \
                                                self.E_max-self.E_min,self.theta_max-self.theta_min,self.ebeamfit_card['ebeam_nexit'])

        self.add_pts(self.data)
        self.resol_fac = 3.

    def get_data(self,input_lhe_evts):

        E_min = 1e20
        E_max = 0.
        theta_min = np.pi/2.
        theta_max = 0.

        data = []
        lhe_evts = lhe_parser.EventFile(input_lhe_evts) 
        nevts = len(lhe_evts)
        nfail = 0

        if self._ebeamdict[self.label] == 11:
            norm = self.ebeamfit_card['flux_norm_electron']
        elif self._ebeamdict[self.label] == -11:
            norm = self.ebeamfit_card['flux_norm_positron']
        elif self._ebeamdict[self.label] == 22:
            norm = self.ebeamfit_card['flux_norm_gamma']
        else:
            exit(-1)
            
        for event in lhe_evts:
            for particle in event:
                if particle.status == -1: # initial state 
                    if particle.pid ==  self._ebeamdict[self.label]:
                        p = lhe_parser.FourMomentum(particle)
                        ctheta = (p.pz/p.norm)
                        stheta = np.sqrt(1.-ctheta**2)
                        # cphi = p.px/np.sqrt(p.px**2+p.py**2)
                        # sphi = p.py/np.sqrt(p.px**2+p.py**2)
                        
                        if ctheta < 0.:
                            weight = 0.
                            nfail +=1
                        else:
                            weight = 1.
                                                
                        theta = np.arccos(ctheta)

                        data.append(WeightedPoint(p.E,theta, \
                                                      norm/nevts*weight))
                        if p.E < E_min:
                            E_min = p.E
                        if p.E > E_max:
                            E_max = p.E
                        if theta < theta_min:
                            theta_min = theta
                        if theta > theta_max:
                            theta_max = theta
                            
        # print E_min,E_max,theta_min,theta_max
        return E_min,E_max,theta_min,theta_max,data

                
    def do_fit(self):
        self.ncores = self.ebeamfit_card['ebeam_ncores']
        # Generate the 2DMesh, output file: cell_fortran.dat
        self.fit(ncores=self.ncores)
        outname = 'cell_fortran_'+self.label+'.dat'
        self.fit('1D_x',ncores=self.ncores)
        
        os.system('mv cell_fortran.dat '+outname)
        if self.ebeamfit_card['ebeam_testplot']:
            meshname = 'mesh2D_'+self.label+'.png'
            plot(outname,meshname)

        outname = 'ehist_'+self.label+'.dat'
        os.system('mv ehist.dat '+outname)
        
        self.ebeamfit_card["ebeampdf"] = True
        self.ebeamfit_card.write('../production/Cards/ebeampdf_fit_card.dat',template='../Cards/ebeampdf_fit_card.dat')



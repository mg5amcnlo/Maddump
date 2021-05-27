from __future__ import absolute_import
import numpy as np
import madgraph.various.banner as bannermod

import os
pjoin = os.path.join

class BremsstrahlungCard(bannermod.RunCard):
    """an object to handle in a nice way the proton-bremsstrahlung information"""

    filename = 'bremsstrahlung_card'
    default_include_file = 'bremsstrahlung_card.inc'
    
    def __new__(cls, finput=None, **opt):
        """Bypass the standard RunCard one"""
        return super(bannermod.RunCard, cls).__new__(cls, finput, **opt)

    def default_setup(self):
        """default value for the bremstrahlung_card.dat"""

        #self.add_param()
        self.add_param("pbeam",100.,comment = "momentum of the primary beam")
        self.add_param("npot",1e18,comment = "number of proton on target")
        self.add_param("z_min",0.1,comment = "minimum value for the energy fraction z")
        self.add_param("z_max",0.9,comment = "maximum value for the energy fraction z")
        self.add_param("pt2_min",1e-14,comment = "minimum value for the square of the transverse momentum pt2")
        self.add_param("pt2_max",1.,comment = "maximum value for the square of the transverse momentum pt2")


        self.add_param("ngen",10000,comment = "number of events to be generated")        
        self.add_param("nfit",100000,comment = "number of points to be generate for the fit ")
        self.add_param("nexit",1000,comment = "exit condition in the fit procedure")
        
    def check_validity(self):
        """ """        
        super(BremsstrahlungCard, self).check_validity()
        
               
    def write(self, output_file, template=None, python_template=False,
              **opt):
        """Write the bremsstrahlung_card in output_file according to template 
           (a path to a valid bremsstrahlung_card)"""

        if not template:
            files = os.listdir('.')
            if 'PLUGIN' in files:
                template = pjoin('PLUGIN/maddump/Templates', 'Cards', 
                                 'bremsstrahlung_card.dat')
            else:
                template = pjoin('../PLUGIN/maddump/Templates', 'Cards', 
                                 'bremsstrahlung_card.dat')
            python_template = True

        super(BremsstrahlungCard, self).write(output_file, template=template,
                                    python_template=python_template, **opt)

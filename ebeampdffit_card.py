from __future__ import absolute_import
import numpy as np
import madgraph.various.banner as bannermod

import os
pjoin = os.path.join

class EbeamPdfFitCard(bannermod.RunCard):
    """an object to handle in a nice way the ebeam pdf fitx information"""

    filename = 'ebeampdf_fit_card'
    default_include_file = 'ebeampdf_fit_card.inc'
    
    def __new__(cls, finput=None, **opt):
        """Bypass the standard RunCard one"""
        return super(bannermod.RunCard, cls).__new__(cls, finput, **opt)

    def default_setup(self):
        """default value for the ebeampdf_fit_card.dat"""

        
        self.add_param("ebeam_dofit", True, comment= 'Turn on/off ebeampdf fit')
        self.add_param("flux_norm_electron", 1., comment= 'Overall normalization of the incoming flux of electrons')
        self.add_param("flux_norm_positron", 1., comment= 'Overall normalization of the incoming flux of positrons')
        self.add_param("flux_norm_gamma", 1., comment= 'Overall normalization of the incoming flux of gammas')
        self.add_param("ebeampdf_interp_method", 'hist', comment= 'method of the interpolation: "hist" (piecewise), "lagrangian"(2nd order polynomial), "spline"(smooth)')                
        self.add_param("ebeam_ncores",4,comment = "number of cores for the fit")
        self.add_param("ebeam_testplot",True,comment = "turn on/off 2D mesh plotting")        
        self.add_param("ebeam_nexit",500,comment = "exit condition in the fit procedure")

        self.add_param("ebeampdf",False,comment = "turn on/off ebeampdf, technical parameter" )
        
    def check_validity(self):
        """ """        
        super(EbeamPdfFitCard, self).check_validity()
        
               
    def write(self, output_file, template=None, python_template=False,
              **opt):
        """Write the ebeampdf_fit_card in output_file according to template 
           (a path to a valid ebeampdf_fit_card)"""

        if not template:
            files = os.listdir('.')
            if 'PLUGIN' in files:
                template = pjoin('PLUGIN/maddump/Templates', 'Cards', 
                                 'ebeampdf_fit_card.dat')
            else:
                template = pjoin('../PLUGIN/maddump/Templates', 'Cards', 
                                 'ebeampdf_fit_card.dat')
            python_template = True

        super(EbeamPdfFitCard, self).write(output_file, template=template,
                                    python_template=python_template, **opt)

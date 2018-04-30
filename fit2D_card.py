import numpy as np
import madgraph.various.banner as bannermod

import os
pjoin = os.path.join

class Fit2DCard(bannermod.RunCard):
    """an object to handle in a nice way the 2d fit and detector information"""

    filename = 'fit2D_card'
    default_include_file = 'fit2D_card.inc'
    
    def __new__(cls, finput=None, **opt):
        """Bypass the standard RunCard one"""
        return super(bannermod.RunCard, cls).__new__(cls, finput, **opt)

    def default_setup(self):
        """default value for the fit2D_card.dat"""
        
        #flux normalization
        self.add_param("flux_norm", 1., comment= 'Overall normalization of the incoming flux without applying cuts')
    
        #detector parameters
        self.add_param("d_target_detector", 1.)
        self.add_param("total_mass", 10.)
        
        #cut
        self.add_param("xc", 0., comment= 'x coordinate central point wrt the beam axis')
        self.add_param("yc", 0., comment= 'y coordinate central point wrt the beam axis')
        
        self.add_param("cylinder", False)
        self.add_param("theta_min", 0.)
        self.add_param("theta_max", np.pi/2.)

        self.add_param("parallelepiped", True)
        self.add_param("x_side", 1.)
        self.add_param("y_side", 1.)

        self.add_param("depth", 1.)

        #technical parameters
        self.add_param("ncores", 4)
        
        
    def check_validity(self):
        """ """        
        super(Fit2DCard, self).check_validity()
        
               
    def write(self, output_file, template=None, python_template=False,
              **opt):
        """Write the fit2D_card in output_file according to template 
           (a path to a valid fit2D_card)"""

        if not template:
            files = os.listdir('.')
            if 'PLUGIN' in files:
                template = pjoin('PLUGIN/maddump/Templates', 'Cards', 
                                 'fit2D_card.dat')
            else:
                template = pjoin('../PLUGIN/maddump/Templates', 'Cards', 
                                 'fit2D_card.dat')
            python_template = True
                   
        super(Fit2DCard, self).write(output_file, template=template,
                                    python_template=python_template, **opt)

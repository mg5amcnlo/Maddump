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

        self.add_param("nevts_interaction", -1, comment= 'nevts_interaction ! Number of events to be generated within the interaction process. A negative value stands for the default behavior, which corresponds to the choice 3* *effective number of input events')

        self.add_param("npoints_cell", 50, comment= 'Minimum numebr of point per cell: internal parameter which specifies the exit condition of the fitting procedure')
        
        self.add_param("rescale_fac", 0, comment= 'to study the fit uncertainties on DM energy spectrum: 0 for central value, 1 for the upper limit, -1 for the lower limit')

        self.add_param("fit_syst", False, comment= 'enable the computation and the storing of fit systematics analyses. They are time comsuming, the default behavior is to take them off')

        self.add_param("e_arr", str(1.) , comment= 'Energy values for which angular distribution will be computed. Format E1,E2,... ')

        self.add_param("nbins", 40, comment= 'number of bins for the 1D angular histograms')
            
        #flux normalization
        self.add_param("flux_norm", 1., comment= 'Overall normalization of the incoming flux without applying cuts')
        self.add_param("prod_xsec_in_norm", False, comment= 'If True the weight is multiplied by the production cross section')
    
        #detector parameters
        self.add_param("d_target_detector", 1.)
        self.add_param("detector_density", 1.)
        
        #cut
        self.add_param("off_axis", False,comment= 'select off-axis mode')
        #self.add_param("xc", 0., comment= 'x coordinate of the central point wrt the beam axis')
        #self.add_param("yc", 0., comment= 'y coordinate central point wrt the beam axis in cm')
        self.add_param("thetac", 0., comment= 'angular coordinate  of the cone wrt the beam axis in rad')
        self.add_param("theta_aperture", 0., comment= 'angular aperture of the cone in rad')
        self.add_param("off_axis_depth", 0., comment= 'cone depth in cm ')
        
        self.add_param("cylinder", False)
        #self.add_param("theta_min", 0.)
        self.add_param("theta_max", np.pi/2.)

        self.add_param("parallelepiped", True)
        self.add_param("x_side", 1.)
        self.add_param("y_side", 1.)

        self.add_param("depth", 1.)

        #technical parameters
        self.add_param("ncores", 4)
        self.add_param("testplot", False)
        
        
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

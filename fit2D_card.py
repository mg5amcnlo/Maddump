import numpy as np
import madgraph.various.banner as bannermod

class Fit2D_Detector(bannermod.RunCard):
    """an object to handle in a nice way the 2d fit and detector information"""
    
    def default_setup(self):
        """default value for the run_card.dat"""
        
        #flux normalization
        self.add_param("flux_norm", 1., comment= 'Overall normalization of the incoming flux without applying cuts')
    
        #detector parameters
        self.add_param("d_target_detector", 1.)
        self.add_param("total_mass", 10d7)
        
        #cut
        self.add_param("xc", 0., comment= 'x coordinate central point wrt the beam axis')
        self.add_param("yc", 0., comment= 'y coordinate central point wrt the beam axis')
        
        self.add_param("circular", False)
        self.add_param("theta_min", 0.)
        self.add_param("theta_max", np.pi/2.)

        self.add_param("rectangular", True)
        self.add_param("x_side", 1.)
        self.add_param("y_side", 1.)

        
    def check_validity(self):
        """ """
        
        super(Fit2D_Detector, self).check_validity()
        
               
    def write(self, output_file, template=None, python_template=False,
              **opt):
        """Write the fit2D_card in output_file according to template 
           (a path to a valid fit2D_card)"""

        if not template:
            template = pjoin('Templates', 'Cards', 
                             'fit2D_card.dat')
            python_template = True
                   
        super(Fit2D_Detector, self).write(output_file, template=template,
                                    python_template=python_template, **opt)

      double precision flux_norm, d_target_detector,detector_density,
     &           thetac,theta_aperture,off_axis_depth,
     &           theta_max,rescale_fac,
     &		 x_side,y_side,depth
      integer nevts_interaction,npoints_cell,ncores,nbins
      logical prod_xsec_in_norm,off_axis,cylinder,parallelepiped,testplot,fit_syst
      character *200 e_arr
      common/fit2D_detector/flux_norm,d_target_detector,detector_density,
     &           thetac,theta_aperture,off_axis_depth,
     &           theta_max,x_side,y_side,depth,rescale_fac,nevts_interaction,
     &		 npoints_cell,ncores,nbins,prod_xsec_in_norm,off_axis,cylinder,
     &		 parallelepiped,testplot,fit_syst,e_arr
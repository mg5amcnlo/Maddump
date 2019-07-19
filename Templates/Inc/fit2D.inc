      double precision flux_norm, d_target_detector,detector_density,depth,
     &           yc,radius,
     &           theta_max,rescale_fac,z_average,a_average,
     &		 x_side,y_side
      integer nevts_interaction,nevts_norm,npoints_cell,ncores,nbins
      logical prod_xsec_in_norm,off_axis,cylinder,parallelepiped,testplot,fit_syst
      character *200 e_arr, interpolation_method
      common/fit2D_detector/flux_norm,d_target_detector,detector_density,depth,
     &           yc,radius,z_average,a_average,
     &           theta_max,x_side,y_side,rescale_fac,nevts_interaction,
     &		 nevts_norm,npoints_cell,ncores,nbins,
     &		 prod_xsec_in_norm,off_axis,cylinder,parallelepiped,testplot,fit_syst,e_arr,
     &		 interpolation_method
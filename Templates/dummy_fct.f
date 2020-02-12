      subroutine get_dummy_x1(sjac, X1, R, pbeam1, pbeam2, stot, shat)
      implicit none
      include 'maxparticles.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'ebeampdf_fit.inc'
c     include 'genps.inc'
      double precision sjac ! jacobian. should be updated not reinit
      double precision X1   ! bjorken X. output
      double precision R    ! random value after grid transfrormation. between 0 and 1
      double precision pbeam1(0:3) ! momentum of the first beam (output)
      double precision pbeam2(0:3) ! momentum of the second beam (output)
      double precision stot        ! total energy  (input and /or output)
      double precision shat        ! output

      double precision pmass(nexternal)
      common/to_mass/  pmass

      double precision xmin,xmax,m1,m2,p0,p3
      double precision Emin,DMpdf         ! Dark Matter minimum energy and pdf

c     global variable to set (or not)
      double precision cm_rap
      logical set_cm_rap
      common/to_cm_rap/set_cm_rap,cm_rap
      
      logical firstcall
      data firstcall/.true./
      save firstcall,Emin

      include '../../Source/ebeampdf_fit_card.inc'
      
c     default behaviour
      set_cm_rap=.false.   ! then cm_rap will be set as 
                           ! .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
      
c     Initialization of the dark matter pdf

      Emin = 0d0
      if(firstcall.and.(.not.ebeampdf)) then
         call init_DMpdf(Emin)
         firstcall = .false.
      endif      
      
c     note: ebeam(1),pmass(1) and ebeam(2),pmass(2) are taken from 'run.inc'
      m1=pmass(1)
      m2=ebeam(2)

c     x fraction for DM pdf
      xmin=Emin/ebeam(1)
      xmax= 1d0
      X1 = (xmax-xmin)*R+xmin 
c     we add in the jacobian the DMpdf

      if (.not.ebeampdf) then
         sjac = sjac*(xmax-xmin)*ebeam(1)*DMpdf(X1*ebeam(1))
      else
         sjac = sjac*(xmax-xmin)*ebeam(1)
      endif
      
c     Set CM rapidity for use in the rap() function
      p0=X1*ebeam(1)+ebeam(2)
      p3=sqrt( (X1*ebeam(1))**2-m1**2 )!-sqrt( ebeam(2)**2-m2**2 )
      cm_rap=.5d0*dlog((p0+p3)/(p0-p3))
      set_cm_rap=.true.

c     Set shat
      shat = m1**2 + m2**2 + 2d0*X1*ebeam(1)*ebeam(2) 
c     $       + sqrt((X1*ebeam(1))**2-m1**2)*sqrt(ebeam(2)**2-m2**2) )

c     Boost in the 'partonic' centre-of-mass frame
c     pbeam(1) and pbeam(2) are returned to the caller
      call mom2cx(sqrt(shat),m1,m2,1d0,0d0,pbeam1(0),pbeam2(0))
      
      return 
      end


      subroutine get_dummy_x1_x2(sjac, X, R, pbeam1, pbeam2, stot,shat)
      implicit none
      include 'maxparticles.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'ebeampdf_fit.inc' 
c      include 'genps.inc'
      double precision sjac     ! jacobian. should be updated not reinit
      double precision X(2)     ! bjorken X. output
      double precision R(2)     ! random value after grid transfrormation. between 0 and 1
      double precision pbeam1(0:3) ! momentum of the first beam (output)
      double precision pbeam2(0:3) ! momentum of the second beam (output)
      double precision stot        ! total energy
      double precision shat        ! output

      double precision pmass(nexternal)
      common/to_mass/  pmass

c     global variable to set (or not)
      double precision cm_rap
      logical set_cm_rap
      common/to_cm_rap/set_cm_rap,cm_rap

      double precision xmin(2),xmax(2),m1,p0,p3
      double precision Emin,DMpdf     ! Dark Matter minimum energy and pdf

      logical firstcall
      data firstcall/.true./
      save firstcall,Emin
      
      include '../../Source/ebeampdf_fit_card.inc'
      
c     default behaviour
      set_cm_rap=.false.   ! then cm_rap will be set as 
                           ! .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
      
c     Initialization of the dark matter pdf

c      if(firstcall.and.(.not.ebeampdf)) then
      if(firstcall) then
         call init_DMpdf(Emin)
         firstcall = .false.
      endif      
      
c     note: ebeam(1),pmass(1) and ebeam(2),pmass(2) are taken from 'run.inc'
c     it is assumed m2=0
      m1=pmass(1)

c     x fraction for DM pdf
      xmin(1)=Emin/ebeam(1)
      xmax(1)=1d0
      X(1) = (xmax(1)-xmin(1))*R(1)+xmin(1) 
c     we add in the jacobian the DMpdf

      if (.not.ebeampdf) then
         sjac = sjac*(xmax(1)-xmin(1))*ebeam(1)*DMpdf(X(1)*ebeam(1))
      else
         sjac = sjac*(xmax(1)-xmin(1))*ebeam(1)
      endif

c     x fraction for parton pdf     
      xmin(2)=0d0
      xmax(2)=1d0
      X(2) = (xmax(2)-xmin(2))*R(2)+xmin(2) 
      sjac = sjac*(xmax(2)-xmin(2))

c     Set CM rapidity for use in the rap() function
      p0=X(1)*ebeam(1)+X(2)*ebeam(2)
      p3=sqrt( (X(1)*ebeam(1))**2 - m1**2)-X(2)*ebeam(2)
      cm_rap=.5d0*dlog((p0+p3)/(p0-p3))
      set_cm_rap=.true.

c     Set shat
      shat = m1**2 + 2d0*(X(1)*ebeam(1)+sqrt((X(1)*ebeam(1))**2-m1**2))
     $     *X(2)*ebeam(2)

c     Boost in the partonic centre-of-mass frame
c     pbeam(1) and pbeam(2) are returned to the caller
      call mom2cx(sqrt(shat),m1,0d0,1d0,0d0,pbeam1(0),pbeam2(0))

      return 
      end

      subroutine init_DMpdf(min)
      implicit none
      include 'fit2D.inc'
      double precision min
      include '../../Source/fit2D_card.inc'
      if (interpolation_method .eq. 'hist') then
c      if (.true.) then
         call init_DMpdf_hist(min)
      else
         call init_DMpdf_spline(min)
      endif
      return
      end

      subroutine init_DMpdf_hist(min)
***   Read the data table for phitilde(E) and load the 1D hist ***  
      implicit none
      include 'fit2D.inc'
      integer nbin,nmax,i,k,ios
      parameter (nmax=100000)
      double precision edges(nmax),hist(nmax),errhist(nmax),a(4),min
      common/fit_hist/nbin,edges,hist,errhist
      
      open(unit=200,file='../ehist.dat',status='old',
     $        err=999)

      k=1
      do 
         read(200,*,iostat=ios) (a(i), i=1,4)
         if (ios.gt.0) then
            write(*,*) 'Something wrong in reading phitilde
     $                  table! Exit!'
            write(*,*) (a(i), i=1,4)
            call exit(-1)
         else if (ios.lt.0) then
            edges(k) = a(2)
            close(200)
            exit
         else
            edges(k) = a(1)
            hist(k) = a(3)+rescale_fac*a(4)
            errhist(k) = a(4)
            k=k+1
         endif
      enddo
      
      nbin = k-1
      min = edges(1)

      return
      
 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end
      

      subroutine init_DMpdf_spline(min)
***   Read the data table for phitilde(E) and perform the 1D fit ***  
      implicit none
      include 'fit2D.inc'
      integer ios
      integer nmax,n,i
      parameter (nmax=100000) 
      double precision x(nmax),y(nmax),w(nmax),a(4)
c      common/phitilde_table/x,y,w,n
      integer nxest,lwrk,iwrk(100000)
      double precision s,min,max,xb,xe,fp,wrk(1000000),resfac
      integer kx,nx,ier
      double precision tx(nmax),c(nmax)
      double precision E,fitvalue
      real ran1
      external ran1
      integer iseed
      data iseed /10/
      common/fit1dim/kx,nx,tx,c,xe,xb,wrk,resfac
      include '../../Source/fit2D_card.inc'
      
      open(unit=200,file='../ehist.dat',status='old',
     $        err=999)

c     store the data table energy, phitilde in the arrays x(n),y(n)
      x=0d0
      y=0d0
      w=0d0

c     loop over infile lines until EoF is reached
      n=3
      do 
         read(200,*,iostat=ios) (a(i), i=1,4)
         if (ios.gt.0) then
            write(*,*) 'Something wrong in reading phitilde
     $                  table! Exit!'
            write(*,*) (a(i), i=1,4)
            call exit(-1)
         else if (ios.lt.0) then
            n=n+1
            x(2) = 0.5d0*(min+x(3))
            y(2) = y(3)/2d0
            w(2) = w(3)*2d0
            x(1) = min
            y(1) = y(2)/2d0
            w(1) = w(2)*2d0
            x(n-1) = 0.5d0*(max+x(n-2)) 
            y(n-1) = y(n-2)/2d0
            w(n-1)=w(n-2)*2d0
            x(n)=max
            y(n)=y(n-1)/2d0
            w(n)=w(n-1)*2d0
            close(200)
            exit
         else
            if (n.eq.3) min=a(1)  
            x(n) = 0.5d0*(a(1)+a(2))
            y(n) = a(3)+rescale_fac*a(4)
            w(n) = 1d0/a(4)
            max=a(2)
            n=n+1
         endif
      enddo

      
c     rescaling for numerical stability
      resfac = maxval(y)
      y(:) = y(:)/resfac
      w(:) = w(:)*resfac
               
c     init parameters for bi-splines fitting  
      xb=min                       !xmin range
      xe=max                       !xmin range
      kx=3                         !x spline order
      s = n                     !smoothing parameter
      nxest= n+kx+1 !n/2        
      lwrk = n*(kx+1)+nxest*(7+3*kx)

c     curve fitting using the FITPACK by Dierckx
      call curfit(0,n,x,y,w,xb,xe,kx,s,nxest,nx,tx,c,fp,wrk,
     *     lwrk,iwrk,ier)

      if(testplot) then
         open(unit=230,file='../../../Cards/fit1D.dat',status='unknown')
         do i=1,100000
c     evaluate the spline interpolation
            E = min + ran1(iseed)*(max-min)
            call splev(tx,nx,c,kx,E,fitvalue,1,ier)
            fitvalue = fitvalue*resfac
            write(230,*) E, fitvalue
         enddo
         close(230)
      endif
      
      return

 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end

      double precision function DMpdf(E)
      implicit none
      include 'fit2D.inc'
      double precision E, DMpdf_hist, DMpdf_spline
      if (interpolation_method .eq. 'hist') then
c      if (.true.) then
         DMpdf = DMpdf_hist(E)
      else
         DMpdf = DMpdf_spline(E)
      endif
      
      return
      end
      
      double precision function DMpdf_hist(E)
      implicit none
      
      integer nbin,nmax,i,k,ia,ib
      parameter (nmax=100000)
      double precision E,edges(nmax),hist(nmax),errhist(nmax),a(4)
      common/fit_hist/nbin,edges,hist,errhist     

      ia = 1
      ib = nbin + 1
      DMpdf_hist = 0d0
      
      if (E < edges(1) .or. E > edges(nbin+1)) return

      do
         k = (ib+ia)/2
         if ( E > edges(k) ) then
            ia = k
            write(*,*) ia,ib,k,E,edges(k),edges(k+1)
            
            if ( E < edges(k+1) ) then
               DMpdf_hist = hist(k)
               write(*,*) k,E,edges(k+1),ia,hist(k),'found'
               return
            endif
            
         else
            ib = k
            write(*,*) ia,ib,k,E,edges(k),edges(k-1)
            if ( E > edges(k-1) ) then
               DMpdf_hist = hist(k-1)
               write(*,*) k,E,edges(k-1),ib,hist(k-1),'found'
               return
            endif
         endif
      enddo
      end

      double precision function DMpdf_spline(E)
      implicit none
      double precision E
      integer nmax,i
      parameter (nmax=100000) 
      integer kx,nx,ier
      double precision tx(nmax),c(nmax),xe,xb,wrk(1000000)
      double precision norm,splint,resfac
      common/fit1dim/kx,nx,tx,c,xe,xb,wrk,resfac
      external splint
      
c     evaluate normalization
c      norm = splint(tx,nx,c,kx,xb,xe,wrk)
c      write(*,*) '#' , norm

c     evaluate the spline interpolation
      call splev(tx,nx,c,kx,E,DMpdf_spline,1,ier)
      DMpdf_spline = DMpdf_spline*resfac

      if (DMpdf_spline.lt.0d0) DMpdf_spline = 0d0
      
      end
      
      double precision function get_ebeampdf(pdg,x)
      implicit none
      include 'maxparticles.inc'
      include 'nexternal.inc'
      include 'run.inc'
      include 'ebeampdf_fit.inc' 
      double precision x,E
      integer nmax,i,npart,pdg,ipart,j
      parameter (nmax=100000,npart=3) 
      integer kx(npart),nx(npart),ier,ipdg(npart)
      double precision tx(npart,nmax),c(npart,nmax),xe(npart),xb(npart),wrk(npart,1000000)
      double precision norm,splint,resfac(npart)
      common/ebeamfit1dim/kx,nx,tx,c,xe,xb,wrk,resfac
      external splint
      integer k,ia,ib,nbin(npart),narr(npart)
      double precision edges(npart,nmax),hist(npart,nmax),errhist(npart,nmax)
      common/ebeampdf_fit_hist/nbin,edges,hist,errhist
      double precision xarr(nmax,npart),yarr(nmax,npart),dy2arr(nmax,npart),dy
      common/ebeampdf_fit_common/xarr,yarr,narr
      common/ebeampdf_fit_spline/dy2arr
      
      logical firstcall
      data firstcall/.true./
      save firstcall,ipdg
      
c     Initialization of the ebeam pdf
      if(firstcall) then
         call init_ebeampdf(ipdg)
         firstcall = .false.
      endif
      
      E = x*ebeam(1)

      do i = 1,npart
         if (ipdg(i) == pdg) ipart=i  
      enddo
      
c     evaluate normalization
c     norm = splint(tx,nx,c,kx,xb,xe,wrk)
c     write(*,*) '#' , norm
      get_ebeampdf = 0d0

      if (ebeampdf_interp_method .eq. 'hist') then

         ia = 1
         ib = nbin(ipart) + 1
      
         if (E < edges(ipart,1) .or. E > edges(ipart,nbin(ipart)+1)) return

         do
            k = (ib+ia)/2
            if ( E > edges(ipart,k) ) then
               ia = k
               if ( E < edges(ipart,k+1) ) then
                 get_ebeampdf  = hist(ipart,k)
                 return
              endif
           else
              ib = k
              if ( E > edges(ipart,k-1) ) then
                 get_ebeampdf = hist(ipart,k-1)
                 return
              endif
           endif
        enddo
 
      elseif (ebeampdf_interp_method .eq. 'lagrangian') then

         if (dlog(E) < xarr(1,ipart) .or. dlog(E) > xarr(narr(ipart),ipart)) return
         call locate(xarr(1,ipart),narr(ipart),dlog(E),j)
         if (j+1>narr(ipart)) then
            j=narr(ipart)-1
         endif
         call polint(xarr(j,ipart),yarr(j,ipart),2,dlog(E),get_ebeampdf,dy)
           
      else
         if (dlog(E) < xarr(1,ipart) .or. dlog(E) > xarr(narr(ipart),ipart)) return

         call spline_int(xarr(1,ipart),yarr(1,ipart),dy2arr(1,ipart),narr(ipart),dlog(E),get_ebeampdf)
c--- old spline routine         
c$$$         if ( (E .gt.xb(ipart) ) .and. (E .lt. xe(ipart) ) ) then 
c$$$c     evaluate the spline interpolation
c$$$            call splev(tx(ipart,:),nx(ipart),c(ipart,:),kx(ipart),E,get_ebeampdf,1,ier)
c$$$            get_ebeampdf = get_ebeampdf*resfac(ipart)
c$$$         endif        
      endif
      
      if (get_ebeampdf.lt.0d0) get_ebeampdf = 0d0
      
      end

      subroutine init_ebeampdf(ipdg)
      implicit none
      include 'ebeampdf_fit.inc' 
      integer npart
      parameter(npart=3)
      integer ipdg(npart)
      if (ebeampdf_interp_method .eq. 'hist') then
         call init_ebeampdf_hist(ipdg)
      elseif (ebeampdf_interp_method .eq. 'lagrangian') then
         call init_ebeampdf_lagrangian(ipdg)
      elseif (ebeampdf_interp_method .eq. 'spline') then
         call init_ebeampdf_spline(ipdg)
      else
         write(*,*) "Wrong interpolation method for ebeam_pdf" 
         call exit(-1)
      endif
      
      end

      subroutine init_ebeampdf_hist(ipdg)
***   Read the data table for phitilde(E) and load the 1D hist ***  
      implicit none
      include 'fit2D.inc'
      integer npart,ipart
      parameter(npart=3)
      integer nbin(npart),nmax,i,k,ios,iun(npart),ipdg(npart),nfit
      parameter (nmax=100000)
      double precision edges(npart,nmax),hist(npart,nmax),errhist(npart,nmax),a(4)
      common/ebeampdf_fit_hist/nbin,edges,hist,errhist

      ipdg = -999999999

      nfit = 0
      open(unit=200,file='../../../../EbeamPdfFit/ehist_electron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 200
         ipdg(nfit) = 11 
      endif
      open(unit=210,file='../../../../EbeamPdfFit/ehist_positron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 210
         ipdg(nfit) = -11 
      endif
      open(unit=220,file='../../../../EbeamPdfFit/ehist_gamma_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 220
         ipdg(nfit) = 22 
      endif
      
      do ipart=1,nfit
         k=1
         do 
            read(iun(ipart),*,iostat=ios) (a(i), i=1,4)
            if (ios.gt.0) then
               write(*,*) 'Something wrong in reading ebeampdf
     $table! Exit!'
               write(*,*) (a(i), i=1,4)
               call exit(-1)
            else if (ios.lt.0) then
               edges(ipart,k) = a(2)
               close(iun(ipart))
               exit
            else
               edges(ipart,k) = a(1)
               hist(ipart,k) = a(3)+rescale_fac*a(4)
               errhist(ipart,k) = a(4)
               k=k+1
            endif
         enddo
         nbin(ipart) = k-1
      enddo
      

      return
      
 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end


      subroutine init_ebeampdf_lagrangian(ipdg)
***   Read the data table for phitilde(E) and load the 1D hist ***  
      implicit none
      include 'fit2D.inc'
      include 'ebeampdf_fit.inc' 
      integer npart,ipart
      parameter(npart=3)
      integer narr(npart),nmax,i,k,ios,iun(npart),ipdg(npart),nfit,j
      parameter (nmax=100000)
      double precision xarr(nmax,npart),yarr(nmax,npart),fitvalue,dy,a(4)
      double precision E,xx(nmax),x(3),y(3)
      common/ebeampdf_fit_common/xarr,yarr,narr
      real ran1
      external ran1
      integer iseed
      data iseed /10/
      character *50 filename

      ipdg = -999999999

      nfit = 0
      open(unit=200,file='../../../../EbeamPdfFit/ehist_electron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 200
         ipdg(nfit) = 11 
      endif
      open(unit=210,file='../../../../EbeamPdfFit/ehist_positron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 210
         ipdg(nfit) = -11 
      endif
      open(unit=220,file='../../../../EbeamPdfFit/ehist_gamma_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 220
         ipdg(nfit) = 22 
      endif
      
      do ipart=1,nfit
         k=1
         do 
            read(iun(ipart),*,iostat=ios) (a(i), i=1,4)
            if (ios.gt.0) then
               write(*,*) 'Something wrong in reading ebeampdf
     $table! Exit!'
               write(*,*) (a(i), i=1,4)
               call exit(-1)
            else if (ios.lt.0) then
               xarr(k,ipart) = dlog(a(2))
               yarr(k,ipart) = 0d0
               close(iun(ipart))
               exit
            else
               if(k==1) then
                  xarr(k,ipart) = dlog(a(1)+1d-10)
                  yarr(k,ipart) = a(3)+rescale_fac*a(4)
                  k=k+1
               endif
               xarr(k,ipart) = 0.5d0* ( dlog(a(1)+1d-10) + dlog(a(2)) )
               yarr(k,ipart) = a(3)+rescale_fac*a(4)
               k=k+1
            endif
         enddo
         narr(ipart)=k
         
         if(ebeam_testplot) then
            write (filename, "(A24,I1,A4)") "../../../Cards/ebeampdf-", ipart,".dat"
            open(unit=250,file=trim(filename),status='unknown')
            do i=1,100000
               E = xarr(1,ipart) + ran1(iseed)*(xarr(narr(ipart),ipart)-xarr(1,ipart))
               call locate(xarr(1,ipart),narr(ipart),E,j)
               if (j+1>narr(ipart)) then
                  j=narr(ipart)-1
               endif
               call polint(xarr(j,ipart),yarr(j,ipart),2,E,fitvalue,dy)
               write(250,*) exp(E), fitvalue
            enddo
            close(250)
         endif

      enddo

      
      return
      
 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end


      subroutine init_ebeampdf_spline(ipdg)
***   Read the data table for phitilde(E) and load the 1D hist ***  
      implicit none
      include 'fit2D.inc'
      include 'ebeampdf_fit.inc' 
      integer npart,ipart
      parameter(npart=3)
      integer narr(npart),nmax,i,k,ios,iun(npart),ipdg(npart),nfit,j
      parameter (nmax=100000)
      double precision xarr(nmax,npart),yarr(nmax,npart),dy2arr(nmax,npart),fitvalue,dy,a(4)
      double precision E,xx(nmax),x(3),y(3)
      common/ebeampdf_fit_common/xarr,yarr,narr
      common/ebeampdf_fit_spline/dy2arr
      real ran1
      external ran1
      integer iseed
      data iseed /10/
      character *50 filename

      ipdg = -999999999

      nfit = 0
      open(unit=200,file='../../../../EbeamPdfFit/ehist_electron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 200
         ipdg(nfit) = 11 
      endif
      open(unit=210,file='../../../../EbeamPdfFit/ehist_positron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 210
         ipdg(nfit) = -11 
      endif
      open(unit=220,file='../../../../EbeamPdfFit/ehist_gamma_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 220
         ipdg(nfit) = 22 
      endif
      
      do ipart=1,nfit
         k=1
         do 
            read(iun(ipart),*,iostat=ios) (a(i), i=1,4)
            if (ios.gt.0) then
               write(*,*) 'Something wrong in reading ebeampdf
     $table! Exit!'
               write(*,*) (a(i), i=1,4)
               call exit(-1)
            else if (ios.lt.0) then
               xarr(k,ipart) = dlog(a(2))
               yarr(k,ipart) = 0d0
               close(iun(ipart))
               exit
            else
               if(k==1) then
                  xarr(k,ipart) = dlog(a(1)+1d-10)
                  yarr(k,ipart) = a(3)+rescale_fac*a(4)
                  k=k+1
               endif
               xarr(k,ipart) = 0.5d0* ( dlog(a(1)+1d-10) + dlog(a(2)) )
               yarr(k,ipart) = a(3)+rescale_fac*a(4)
               k=k+1
            endif
         enddo
         narr(ipart)=k
         
         call spline(xarr(1,ipart),yarr(1,ipart),narr(ipart),1d99,1d99,dy2arr(1,ipart))
         if(ebeam_testplot) then
            write (filename, "(A24,I1,A4)") "../../../Cards/ebeampdf-", ipart,".dat"
            open(unit=250,file=trim(filename),status='unknown')
            do i=1,100000
               E = xarr(1,ipart) + ran1(iseed)*(xarr(narr(ipart),ipart)-xarr(1,ipart))
               call spline_int(xarr(1,ipart),yarr(1,ipart),dy2arr(1,ipart),narr(ipart),E,fitvalue)
               write(250,*) exp(E), fitvalue
            enddo
            close(250)
         endif

      enddo

      
      return
      
 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end
      

      

      subroutine init_ebeampdf_spline_old(ipdg)
***   Read the data table for the ebeam pdf and perform the 1D fit ***  
      implicit none
      include 'fit2D.inc'
      include 'ebeampdf_fit.inc' 
      include 'nexternal.inc'
      include 'maxamps.inc'
      integer ios,ipart,nfit
      integer nmax,n,i,npart
      parameter (nmax=100000,npart=3) 
      double precision x(npart,nmax),y(npart,nmax),w(npart,nmax),a(4)
c      common/phitilde_table/x,y,w,n
      integer nxest,lwrk,iwrk(100000),iun(npart),ipdg(npart)
      double precision s,min,max,xb(npart),xe(npart),fp,wrk(npart,1000000),resfac(npart)
      integer kx(npart),nx(npart),ier
      double precision tx(npart,nmax),c(npart,nmax)
      double precision E,fitvalue
      real ran1
      external ran1
      integer iseed
      data iseed /10/
      common/ebeamfit1dim/kx,nx,tx,c,xe,xb,wrk,resfac
      character *50 filename

      ipdg = -999999999

      nfit = 0
      open(unit=200,file='../../../../EbeamPdfFit/ehist_electron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 200
         ipdg(nfit) = 11 
      endif
      open(unit=210,file='../../../../EbeamPdfFit/ehist_positron_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 210
         ipdg(nfit) = -11 
      endif
      open(unit=220,file='../../../../EbeamPdfFit/ehist_gamma_pdf.dat',status='old',
     $     iostat=ios)
      if (ios.gt.0) then
         continue
      else
         nfit = nfit+1
         iun(nfit) = 220
         ipdg(nfit) = 22 
      endif
      

c     store the data table energy, phitilde in the arrays x(n),y(n)
      x=0d0
      y=0d0
      w=0d0

      do ipart=1,nfit
c     loop over infile lines until EoF is reached
      n=3
      do 
         read(iun(ipart),*,iostat=ios) (a(i), i=1,4)
         if (ios.gt.0) then
            write(*,*) 'Something wrong in reading phitilde
     $                  table! Exit!'
            write(*,*) (a(i), i=1,4)
            call exit(-1)
         else if (ios.lt.0) then
            n=n+1
            x(ipart,2) = 0.5d0*(min+x(ipart,3))
            y(ipart,2) = y(ipart,3)/2d0
            w(ipart,2) = w(ipart,3)*2d0
            x(ipart,1) = min
            y(ipart,1) = y(ipart,2)/2d0
            w(ipart,1) = w(ipart,2)*2d0
            x(ipart,n-1) = 0.5d0*(max+x(ipart,n-2)) 
            y(ipart,n-1) = y(ipart,n-2)/2d0
            w(ipart,n-1)=w(ipart,n-2)*2d0
            x(ipart,n)=max
            y(ipart,n)=y(ipart,n-1)/2d0
            w(ipart,n)=w(ipart,n-1)*2d0
            close(iun(ipart))
            exit
         else
            if (n.eq.3) min=a(1)  
            x(ipart,n) = 0.5d0*(a(1)+a(2))
            y(ipart,n) = a(3)+rescale_fac*a(4)
            w(ipart,n) = 1d0/a(4)
            max=a(2)
            n=n+1
         endif
      enddo

      
c     rescaling for numerical stability
      resfac(ipart) = maxval(y(ipart,:))
      y(ipart,:) = y(ipart,:)/resfac(ipart)
      w(ipart,:) = w(ipart,:)*resfac(ipart)
               
c     init parameters for bi-splines fitting  
      xb(ipart)=min                       !xmin range
      xe(ipart)=max                       !xman range
      kx(ipart)=3                         !x spline order
      s = n                     !smoothing parameter
      nxest = n+kx(ipart)+1 !n/2        
      lwrk = n*(kx(ipart)+1)+nxest*(7+3*kx(ipart))

c     curve fitting using the FITPACK by Dierckx
      call curfit(0,n,x(ipart,:),y(ipart,:),w(ipart,:),xb(ipart),xe(ipart),
     *     kx(ipart),s,nxest,nx(ipart),tx(ipart,:),c(ipart,:),fp,
     *     wrk(ipart,:),lwrk,iwrk,ier)
      
      if(testplot) then
         write (filename, "(A24,I1,A4)") "../../../Cards/ebeampdf-", ipart,".dat"
         open(unit=250,file=trim(filename),status='unknown')
         do i=1,100000
c     evaluate the spline interpolation
            E = min + ran1(iseed)*(max-min)
            call splev(tx(ipart,:),nx(ipart),c(ipart,:),kx(ipart),E,fitvalue,1,ier)
            fitvalue = fitvalue*resfac(ipart)
            write(250,*) E, fitvalue
         enddo
         close(250)
      endif
      enddo
      
      return

 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)
      end


cz  subroutine POLINT that is printed in Numerical Recipes.
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      double precision x,xx(n)
c      Given an array xx(1:n), and given a value x, returns a value j such that x is between
c      xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0 or j=n is returned to indicate that x is out of range.
      INTEGER jl,jm,ju
c      write(*,*) xx,'locate'
c      write(*,*) n,x,'locate'
      jl=0 
      ju=n+1 
 10   if(ju-jl.gt.1)then
         jm=(ju+jl)/2 
         if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
            jl=jm 
         else
            ju=jm 
         endif
         goto 10 
      endif 
      if(x.eq.xx(1)) then 
         j=1
      else if(x.eq.xx(n)) then
         j=n-1
      else
         j=jl
      endif
      return 
      END

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
c$$$      Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
c$$$      x1 < x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
c$$$      function at points 1 and n, respectively, this routine returns an array y2(1:n) of
c$$$      length n which contains the second derivatives of the interpolating function at the tabulated
c$$$      points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set
c$$$      the corresponding boundary condition for a natural spline, with zero second derivative on
c$$$      that boundary.
c$$$  Parameter: NMAX is the largest anticipated value of n.
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then 
         y2(1)=0d0
         u(1)=0d0
      else 
         y2(1)=-0.5d0
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1 
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=( 6d0*( ( y(i+1)-y(i) ) / ( x(i+1)-x(i) ) -( y(i)-y(i-1) )
     .        / ( x(i)-x(i-1) ) ) / ( x(i+1)-x(i-1) )- sig*u(i-1) )/p
      enddo
      if (ypn.gt..99e30) then 
         qn=0d0
         un=0d0
      else 
         qn=0.5d0
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1 
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo 
      return
      END      

      SUBROUTINE spline_int(xa,ya,y2a,n,x,y)
      INTEGER n
      double precision x,y,xa(n),y2a(n),ya(n)
c$$$      Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
c$$$      xai ’s in order), and given the array y2a(1:n), which is the output from spline above,
c$$$      and given a value of x, this routine returns a cubic-spline interpolated value y.
      integer k,khi,klo
      double precision a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif 
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splin'
      a=(xa(khi)-x)/h 
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     .     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END      
      

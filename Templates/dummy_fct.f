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
      parameter (nmax=1000)
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
      parameter (nmax=1000) 
      double precision x(nmax),y(nmax),w(nmax),a(4)
c      common/phitilde_table/x,y,w,n
      integer nxest,lwrk,iwrk(1000)
      double precision s,min,max,xb,xe,fp,wrk(100000),resfac
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
      nxest= n/2        
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
      parameter (nmax=1000)
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
      parameter (nmax=1000) 
      integer kx,nx,ier
      double precision tx(nmax),c(nmax),xe,xb,wrk(100000)
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
      integer nmax,i,npart,pdg,ipart
      parameter (nmax=1000,npart=3) 
      integer kx(npart),nx(npart),ier,ipdg(npart)
      double precision tx(npart,nmax),c(npart,nmax),xe(npart),xb(npart),wrk(npart,100000)
      double precision norm,splint,resfac(npart)
      common/ebeamfit1dim/kx,nx,tx,c,xe,xb,wrk,resfac
      external splint
      integer k,ia,ib,nbin(npart)
      double precision edges(npart,nmax),hist(npart,nmax),errhist(npart,nmax)
      common/ebeampdf_fit_hist/nbin,edges,hist,errhist
      
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
 
      else
         if ( (E .gt.xb(ipart) ) .and. (E .lt. xe(ipart) ) ) then 
c     evaluate the spline interpolation
            call splev(tx(ipart,:),nx(ipart),c(ipart,:),kx(ipart),E,get_ebeampdf,1,ier)
            get_ebeampdf = get_ebeampdf*resfac(ipart)
         endif
         
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
c      if (.true.) then
         call init_ebeampdf_hist(ipdg)
      else
         call init_ebeampdf_spline(ipdg)
      endif
      
      end

      subroutine init_ebeampdf_hist(ipdg)
***   Read the data table for phitilde(E) and load the 1D hist ***  
      implicit none
      include 'fit2D.inc'
      integer npart,ipart
      parameter(npart=3)
      integer nbin(npart),nmax,i,k,ios,iun(npart),ipdg(npart),nfit
      parameter (nmax=1000)
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


      subroutine init_ebeampdf_spline(ipdg)
***   Read the data table for the ebeam pdf and perform the 1D fit ***  
      implicit none
      include 'fit2D.inc'
      include 'ebeampdf_fit.inc' 
      include 'nexternal.inc'
      include 'maxamps.inc'
      integer ios,ipart,nfit
      integer nmax,n,i,npart
      parameter (nmax=1000,npart=3) 
      double precision x(npart,nmax),y(npart,nmax),w(npart,nmax),a(4)
c      common/phitilde_table/x,y,w,n
      integer nxest,lwrk,iwrk(1000),iun(npart),ipdg(npart)
      double precision s,min,max,xb(npart),xe(npart),fp,wrk(npart,100000),resfac(npart)
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
      nxest = n/2        
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

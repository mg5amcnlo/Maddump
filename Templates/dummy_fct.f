      subroutine get_dummy_x1(sjac, X1, R, pbeam1, pbeam2, stot, shat)
      implicit none
      include 'maxparticles.inc'
      include 'nexternal.inc'
      include 'run.inc'
c      include 'genps.inc'
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
      
c     default behaviour
      set_cm_rap=.false.   ! then cm_rap will be set as 
                           ! .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
      
c     Initialization of the dark matter pdf

      if(firstcall) then
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
      sjac = sjac*(xmax-xmin)*ebeam(1)*DMpdf(X1*ebeam(1))
      
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
      
      
c     default behaviour
      set_cm_rap=.false.   ! then cm_rap will be set as 
                           ! .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
      
c     Initialization of the dark matter pdf

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
      sjac = sjac*(xmax(1)-xmin(1))*ebeam(1)*DMpdf(X(1)*ebeam(1))

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
     $              err=999)

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
            y(n) = a(3)
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
      call splev(tx,nx,c,kx,E,DMpdf,1,ier)
      DMpdf = Dmpdf*resfac
      end



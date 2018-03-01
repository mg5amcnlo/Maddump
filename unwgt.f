      SUBROUTINE write_leshouche(p,wgt,numproc,do_write_events)
C**************************************************************************
C     Writes out information for event
C**************************************************************************
      IMPLICIT NONE
c
c     Constants
c
      double precision zero
      parameter       (ZERO = 0d0)
      include 'genps.inc'
      include 'nexternal.inc'
      include 'maxamps.inc'
      include 'message.inc'
      include 'cluster.inc'
      include 'run.inc'
      include 'run_config.inc'

c
c     Arguments
c
      double precision p(0:3,nexternal),wgt
      integer numproc
      logical do_write_events
c
c     Local
c
      integer i,j,k,iini,ifin
      double precision sum_wgt,sum_wgt2, xtarget,targetamp(maxflow)
      integer ip, np, ic, nc
      integer ida(2),ito(-nexternal+3:nexternal),ns,nres,ires,icloop
      integer iseed
      double precision pboost(0:3)
      double precision beta, get_betaz
      double precision ebi(0:3), ebo(0:3)
      double precision ptcltmp(nexternal), pdum(0:3)
     
      integer idup(nexternal,maxproc,maxsproc)
      integer mothup(2,nexternal)
      integer icolup(2,nexternal,maxflow,maxsproc)

      integer nsym

      integer ievent
      logical flip

      real ran1
      external ran1

      character*40 cfmt

      double precision E, theta, rv(2), phi
      double precision pi
      parameter (pi = 3.1415926d0)      

C     
C     GLOBAL
C
      double precision twgt, maxwgt,swgt(maxevents)
      integer                             lun, nw, itmin
      common/to_unwgt/twgt, maxwgt, swgt, lun, nw, itmin

      integer          IPSEL
      COMMON /SubProc/ IPSEL

      Double Precision amp2(maxamps), jamp2(0:maxflow)
      common/to_amps/  amp2,       jamp2

      character*101       hel_buf
      common/to_helicity/hel_buf

      integer           mincfig, maxcfig
      common/to_configs/mincfig, maxcfig

      double precision stot,m1,m2
      common/to_stot/stot,m1,m2
c
c     Data
c
      include 'leshouche.inc'
      data iseed /0/

      double precision pmass(nexternal), tmp
      common/to_mass/  pmass

c      integer ncols,ncolflow(maxamps),ncolalt(maxamps)
c      common/to_colstats/ncols,ncolflow,ncolalt,ic
c      data ncolflow/maxamps*0/
c      data ncolalt/maxamps*0/

      include 'coupl.inc'

      include 'lhe_event_infos.inc'
      data AlreadySetInBiasModule/.False./

      include 'symswap.inc'
C-----
C  BEGIN CODE
C-----
      
      if ((nw .ge. maxevents).and.do_write_events) return

C     if all the necessary inputs to write the events have already been
C     computed in the bias module, then directly jump to write_events
      if (AlreadySetInBiasModule) then
        goto 1123
      endif

c
c     In case of identical particles symmetry, choose assignment
c
      xtarget = ran1(iseed)*nsym
      jsym = 1
      do while (xtarget .gt. jsym .and. jsym .lt. nsym)
         jsym = jsym+1
      enddo
c
c     Fill jpart color and particle info
c
      do i=1,nexternal
         jpart(1,isym(i,jsym)) = idup(i,ipsel,numproc)
         jpart(2,isym(i,jsym)) = mothup(1,i)
         jpart(3,isym(i,jsym)) = mothup(2,i)
c        Color info is filled in mothup
         jpart(4,isym(i,jsym)) = 0
         jpart(5,isym(i,jsym)) = 0
         jpart(6,isym(i,jsym)) = 1
      enddo
      do i=1,nincoming
         jpart(6,isym(i,jsym))=-1
      enddo

c   Set helicities
c      write(*,*) 'Getting helicity',hel_buf(1:50)
      read(hel_buf,'(20i5)') (jpart(7,isym(i, jsym)),i=1,nexternal)
c      write(*,*) 'ihel',jpart(7,1),jpart(7,2)

c   Fix ordering of ptclus
      do i=1,nexternal
        ptcltmp(isym(i,jsym)) = ptclus(i)
      enddo
      do i=1,nexternal
        ptclus(i) = ptcltmp(i)
      enddo

c     Check if we have flipped particle 1 and 2, and flip back
      flip = .false.
      if (p(3,1).lt.0) then
         do j=0,3
            pdum(j)=p(j,1)
            p(j,1)=p(j,2)
            p(j,2)=pdum(j)
         enddo
         flip = .true.
      endif

c
c     Boost momentum to lab frame
c
      pboost(0)=1d0
      pboost(1)=0d0
      pboost(2)=0d0
      pboost(3)=0d0

      if (nincoming.eq.2) then
         
         if (xbk(1) .gt. 0d0 .and. xbk(1) .le. 1d0 .and.
     $        xbk(2) .gt. 0d0 .and. xbk(2) .le. 1d0) then            
            if(lpp(2).ne.0.and.
     $           (xbk(1).eq.1d0.or.pmass(1).eq.0d0)) then
!construct the beam momenta in each frame and compute the related (z)boost
               ebi(0) = p(0,1)/xbk(1) ! this assumes that particle 1 is massless or mass equal to beam
               ebi(1) = 0
               ebi(2) = 0
               ebi(3) = DSQRT(ebi(0)**2-m1**2)
               ebo(0) = ebeam(1)
               ebo(1) = 0
               ebo(2) = 0
               ebo(3) = DSQRT(ebo(0)**2-m1**2)
               beta = get_betaz(ebi, ebo)
            else
               ebi(0) = p(0,2)/xbk(2) ! this assumes that particle 2 is massless or mass equal to beam
               ebi(1) = 0
               ebi(2) = 0
               ebi(3) = -1d0*DSQRT(ebi(0)**2-m2**2)
               ebo(0) = ebeam(2)
               ebo(1) = 0
               ebo(2) = 0
               ebo(3) = -1d0*DSQRT(ebo(0)**2-m2**2)
               beta = get_betaz(ebi, ebo)
!     wrong boost if both parton are massive!
            endif
         else
            write(*,*) 'Warning bad x1 or x2 in write_leshouche',
     $           xbk(1),xbk(2)
         endif
         do j=1,nexternal
            call zboost_with_beta(p(0,j),beta,pb(0,isym(j,jsym)))
            pb(4,isym(j,jsym))=pmass(j)
         enddo
c**     new incoming beam alternative start
         if(lpp(1).eq.9) then

            E=xbk(1)*ebeam(1)   ! energy of beam 1 particle in the lab
            
c     Boost to the lab frame
            call momntx(E,pmass(1),1d0,0d0,pboost)
            if(lpp(2) .eq. 1) then ! DIS
               pboost(0) = pboost(0) + xbk(2)*ebeam(2)
               pboost(3) = pboost(3) - xbk(2)*ebeam(2)
            elseif(lpp(2) .eq. 0) then ! elastic scattering 
               pboost(0) = pboost(0) + ebeam(2)
               pboost(3) = pboost(3)
            endif
            
            do j=1,nexternal
               call boostx(p(0,j),pboost,pb(0,isym(j,jsym)))
               pb(4,isym(j,jsym))=pmass(j)
            enddo

c            write(*,*) 'Energy DM:', E
c            write(*,*) (p(:,j),NEW_LINE('A'), j=1,nexternal)
c            write(*,*) (pb(:,isym(j,jsym)),NEW_LINE('A'), j=1,nexternal)

c     Pick a theta value at fixed E distributed according to the 2D DM 
c     histograms flux. The azimuthal angle is set at random. 
            rv(1)=ran1(iseed+1)
            rv(2)=ran1(iseed+2)
            call picktheta(E,rv,theta)
            phi = 2d0*pi*ran1(iseed+3)

c     Rotate the event:
c     pboost is now the reference 4-vector for the rotation
c     it corresponds to the rotated dark matter 4-vector 
            call momntx(E,pmass(1),dcos(theta),phi,pboost)
           
            do j=1,nexternal
               p(:,j)=pb(0:3,isym(j,jsym))
               call rotxxx(p(0,j),pboost,pb(0,isym(j,jsym)))
            enddo

c            write(*,*) 'theta,phi: ', theta,phi 
c            write(*,*) (pb(:,isym(j,jsym)),NEW_LINE('A'), j=1,nexternal)

         endif
c**   new incoming beam alternative end

      else
         do j=1,nexternal
            call boostx(p(0,j),pboost,pb(0,isym(j,jsym)))
!     Add mass information in pb(4)
            pb(4,isym(j,jsym))=pmass(j)
         enddo
      endif
c     
c     Add info on resonant mothers
c     
      call addmothers(ipsel,jpart,pb,isym,jsym,sscale,aaqcd,aaqed,buff,
     $     npart,numproc,flip)
      
      if (nincoming.eq.1)then
         do i=-nexternal+3,2*nexternal-3
            if (jpart(2,i).eq.1)then
               jpart(3,i) = 0
            endif
         enddo
      endif
c     
c     Write events to lun
c     
      if(q2fact(1).gt.0.and.q2fact(2).gt.0)then
         sscale = sqrt(max(q2fact(1),q2fact(2)))
      else if(q2fact(1).gt.0)then
         sscale = sqrt(q2fact(1))
      else if(q2fact(2).gt.0)then
         sscale = sqrt(q2fact(2))
      else
         sscale = 0d0
      endif
      aaqcd = g*g/4d0/3.1415926d0
      aaqed = gal(1)*gal(1)/4d0/3.1415926d0
      
      if (btest(mlevel,3)) then
         write(*,*)' write_leshouche: SCALUP to: ',sscale
      endif
      
c     Write out buffer for systematics studies
      ifin=1
      if(use_syst)then
c         print *,'Systematics:'
c         print *,'s_scale: ',s_scale
c         print *,'n_qcd,n_alpsem: ',n_qcd,n_alpsem
c         print *,'s_qalps: ',(s_qalps(I),I=1,n_alpsem) 
c         print *,'n_pdfrw: ',n_pdfrw
c         print *,'i_pdgpdf: ',((i_pdgpdf(i,j),i=1,n_pdfrw(j)),j=1,2)
c         print *,'s_xpdf: ',((s_xpdf(i,j),i=1,n_pdfrw(j)),j=1,2)
c         print *,'s_qpdf: ',((s_qpdf(i,j),i=1,n_pdfrw(j)),j=1,2)
         s_buff(1) = '<mgrwt>'
         write(s_buff(2), '(a,I3,E15.8,a)') '<rscale>',n_qcd-n_alpsem,
     $        s_scale,'</rscale>'
         if(n_alpsem.gt.0) then
            write(cfmt,'(a,I1,a)') '(a,I3,',n_alpsem,'E15.8,a)'
            write(s_buff(3), cfmt) '<asrwt>',n_alpsem,
     $           (s_qalps(I),I=1,n_alpsem) ,'</asrwt>'
         else
            write(s_buff(3), '(a)') '<asrwt>0</asrwt>'
         endif
         if(n_pdfrw(1).gt.0)then
            if(2*n_pdfrw(1).lt.10) then
               write(cfmt,'(a,I1,a,I1,a)') '(a,I3,',
     $              n_pdfrw(1),'I9,',2*n_pdfrw(1),'E15.8,a)'
            else
               write(cfmt,'(a,I1,a,I2,a)') '(a,I3,',
     $              n_pdfrw(1),'I9,',2*n_pdfrw(1),'E15.8,a)'
            endif
            write(s_buff(4), cfmt) '<pdfrwt beam="1">',
     $           n_pdfrw(1),(i_pdgpdf(i,1),i=1,n_pdfrw(1)),
     $           (s_xpdf(i,1),i=1,n_pdfrw(1)),
     $           (s_qpdf(i,1),i=1,n_pdfrw(1)),
     $           '</pdfrwt>'
         else
            write(s_buff(4), '(a)') '<pdfrwt beam="1">0</pdfrwt>'
         endif
         if(n_pdfrw(2).gt.0)then
            if(2*n_pdfrw(2).lt.10) then
               write(cfmt,'(a,I1,a,I1,a)') '(a,I3,',
     $              n_pdfrw(2),'I9,',2*n_pdfrw(2),'E15.8,a)'
            else
               write(cfmt,'(a,I1,a,I2,a)') '(a,I3,',
     $              n_pdfrw(2),'I9,',2*n_pdfrw(2),'E15.8,a)'
            endif
            write(s_buff(5), cfmt) '<pdfrwt beam="2">',
     $           n_pdfrw(2),(i_pdgpdf(i,2),i=1,n_pdfrw(2)),
     $           (s_xpdf(i,2),i=1,n_pdfrw(2)),
     $           (s_qpdf(i,2),i=1,n_pdfrw(2)),
     $           '</pdfrwt>'
         else
            write(s_buff(5), '(a)') '<pdfrwt beam="2">0</pdfrwt>'
         endif
         write(s_buff(6), '(a,E15.8,a)') '<totfact>',s_rwfact,
     $        '</totfact>'
         s_buff(7) = '</mgrwt>'
      endif

c     Write out buffers for clustering info
      nclus=0
      if(icluster(1,1).ne.0 .and. ickkw.ne.0 .and. clusinfo)then
         nclus=nexternal
         write(buffclus(1),'(a)')'<clustering>'
         do i=1,nexternal-2
            write(buffclus(i+1),'(a13,f9.3,a2,4I3,a7)') '<clus scale="',
     $           dsqrt(pt2ijcl(i)),'">',(icluster(j,i),j=1,4),'</clus>'
         enddo
         write(buffclus(nexternal),'(a)')'</clustering>'
      endif
      
C     If the arguments of write_event have already been set in the
C     bias module, then the beginning of the routine will directly
C     jump here.

 1123 continue
      if (.not.do_write_events) then
         return
      endif

c     Store weight for event
      nw = nw+1
      swgt(nw)=wgt
      
      call write_event(lun,pb(0,1),wgt,npart,jpart(1,1),ngroup,
     &   sscale,aaqcd,aaqed,buff,use_syst,s_buff,nclus,buffclus)
      if(btest(mlevel,1))
     &   call write_event(6,pb(0,1),wgt,npart,jpart(1,1),ngroup,
     &   sscale,aaqcd,aaqed,buff,use_syst,s_buff,nclus,buffclus)

      end
            
      subroutine picktheta(E,rv,theta)
      implicit none
      integer maxcell
      parameter (maxcell=1000000)
      integer ncells,n,i,j,icell(maxcell),ios
      double precision E,rv(2),theta
      double precision cells(maxcell,4),a(4),w(maxcell),wtot,Emin,Emax
      double precision area,eps,Em,Ep,s
      common/celltable/cells,ncells
      logical fstcall
      data fstcall/.true./
      save fstcall,eps

c     At the first call, read and store the cell parameters from the 
c     cell_fortran.dat file
      if (fstcall) then
         open(unit=210,file='../cell_fortran.dat',status='old',
     *        err=999)
c     loop over infile lines until EoF is reached
         eps = 1d9
         ncells=1
         do 
            read(210,*,iostat=ios) a(1),a(2),a(3),a(4)
            if (ios.gt.0) then
               write(*,*) 'Something wrong in reading cell
     $                     table! Exit!'
               write(*,*) a(1),a(2),a(3),a(4)
               call exit(-1)
            else if (ios.lt.0) then
               ncells=ncells-1
               close(210)
               fstcall = .false.
               do i= 1,ncells
                  if(cells(i,3).lt.eps)  eps=cells(i,3)
               enddo
               eps= eps/10d0
               exit
            else
               cells(ncells,:)= a(:)
               ncells=ncells+1

               if(ncells.gt.maxcell) then
                  write(*,*) 'Error: the number of cells of the 2D mesh exceeds
     $                        the allowed maxcell value. Exit!'
                  call exit(-1)
               endif
            endif
         enddo
      endif
      
c     look for the cells in which the energy interval [E-eps,E+eps]
c     is contained. The case in which there is crossing with 2 cells is 
c     considered.
c     The weights and their normalization are computed 
      n=1
      w=0d0
      wtot=0d0
      icell = 0
      do i=1,ncells
         Emin = cells(i,1)
         Emax = cells(i,1)+cells(i,3)
         Ep = E+eps
         Em = E-eps
         if (Em.ge.Emin.and.Em.lt.Emax) then
            icell(n) = i
            if(Ep.le.Emax) then
               w(n) = (2d0*eps)/cells(i,3)
               wtot = wtot + w(n)
            else 
               w(n) = (Emax-Em)/cells(i,3)
               wtot = wtot + w(n)
            endif
            n = n+1
         elseif (Em.le.Emin.and.Ep.gt.Emin) then
            icell(n) = i
            if(Ep.le.Emax) then
               w(n) = (Ep-Emin)/cells(i,3)
               wtot = wtot + w(n)
            else 
               w(n) = (Emax-Emin)/cells(i,3)
               wtot = wtot + w(n)
            endif
            n = n+1
         endif
      enddo
      n=n-1

c     check if something went wrong
      if (n.eq.0) then 
         write(*,*) 'Error has occurred with the energy value: ', E
         stop
      endif

c     pick a theta value according to the hit cells and their weights;
c     once a cell is selected, a value of theta is taken uniformly inside it 
      s=0d0
      do j = 1,n
         s = s + w(j)/wtot
         if (rv(1).lt.s) then
            theta = cells(icell(j),2) + rv(2)*cells(icell(j),4)
            return
         endif
      enddo

 999  write(*,*) 'Cannot open input data file, exit!'
      call exit(-1)

      end

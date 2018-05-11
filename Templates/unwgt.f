      subroutine store_events(force_max_wgt, scale_to_xsec)
C**************************************************************************
C     Takes events from scratch file (lun) and writes them to a permanent
c     file  events.dat
c     if force_max_weight =-1, then get it automatically (for a given truncation)
c     if xscale=0 then the sum of the weight will be reweighted to the cross-section.
c     computed from the last 3 iteration. otherwise the weight of each event
c     will be multiply by that value.
C**************************************************************************
      IMPLICIT NONE
c
c     Constants
c
      include 'genps.inc'
      include 'nexternal.inc'
      include 'run_config.inc'
c
c     Arguments
c
      double precision force_max_wgt 
      logical scale_to_xsec
c
c     Local
c
      integer i, lunw, ic(7,2*nexternal-3), n, j
      logical done
      double precision wgt,p(0:4,2*nexternal-3)
      double precision xsec,xsecabs,xerr,xtot
      double precision xsum, xover, target_wgt
      double precision orig_Wgt(maxevents)
      double precision xscale
      logical store_event(maxevents)
      integer iseed, nover, nstore
      double precision scale,aqcd,aqed
      double precision random
      integer ievent
      character*1000 buff
      logical u_syst
      character*(s_bufflen) s_buff(7)
      integer nclus
      character*(clus_bufflen) buffclus(nexternal)

      logical use_vert_displacement
      character*(30) buff_vert_displacement(3)      

C     
C     GLOBAL
C
      double precision twgt, maxwgt,swgt(maxevents)
      integer                             lun, nw, itmin
      common/to_unwgt/twgt, maxwgt, swgt, lun, nw, itmin

      integer                   neventswritten
      common /to_eventswritten/ neventswritten
      
      integer th_nunwgt
      double precision th_maxwgt
      common/theoretical_unwgt_max/th_maxwgt, th_nunwgt

c      save neventswritten

      integer ngroup
      common/to_group/ngroup

c
c     external
c
      real xran1

      data iseed/0/
      data neventswritten/0/
C-----
C  BEGIN CODE
C-----
c
c     First scale all of the events to the total cross section
c

      if (nw .le. 0) return
      if (scale_to_xsec) then
         call sample_result(xsecabs,xsec,xerr,itmin)
         if (xsecabs .le. 0) return !Fix by TS 12/3/2010
      else
         xscale = nw*twgt
      endif
      xtot=0
      call dsort(nw, swgt)
      do i=1,nw
         xtot=xtot+dabs(swgt(i))
      enddo
c
c     Determine minimum target weight given truncation parameter
c
      xsum = 0d0
      i = nw
      do while (xsum-dabs(swgt(i))*(nw-i) .lt. xtot*trunc_max
     $ .and. i .gt. 2)
         xsum = xsum + dabs(swgt(i))
         i = i-1
      enddo
      if (i .lt. nw) i=i+1
      th_maxwgt = dabs(swgt(i))
      if ( force_max_wgt.lt.0)then
         target_wgt = dabs(swgt(i))
      elseif (.not.scale_to_xsec) then
         target_wgt = force_max_wgt / xscale
      else
         stop 1
      endif
c
c     Select events which will be written
c
      xsum = 0d0
      nstore = 0
      th_nunwgt = 0
      rewind(lun)
      done = .false. 
      do i=1,nw
         if (.not. done) then
            call read_event(lun,P,wgt,n,ic,ievent,scale,aqcd,aqed,buff,
     $           u_syst,s_buff,nclus,buffclus,use_vert_displacement,buff_vert_displacement,done)
         else
            wgt = 0d0
         endif
         random = xran1(iseed)
         if (dabs(wgt) .gt. target_wgt*random) then
            xsum=xsum+max(dabs(wgt),target_Wgt)
            store_event(i)=.true.
            nstore=nstore+1
         else
            store_event(i) = .false.
         endif
c           we use the same seed for the two evaluation of the unweighting efficiency
         if (dabs(wgt) .gt. th_maxwgt*random) then
            th_nunwgt = th_nunwgt +1
         endif
      enddo
      if (scale_to_xsec)then
         xscale = xsecabs/xsum
      endif
      target_wgt = target_wgt*xscale
      th_maxwgt = th_maxwgt*xscale

      rewind(lun)
c     JA 8/17/2011 Don't check for previously stored events
c      if (nstore .le. neventswritten) then
c         write(*,*) 'No improvement in events',nstore, neventswritten
c         return
c      endif
      lunw = 25
      open(unit = lunw, file='events.lhe', status='unknown')
      done = .false.
      i=0      
      xtot = 0
      xover = 0
      nover = 0
      do j=1,nw
         if (.not. done) then
            call read_event(lun,P,wgt,n,ic,ievent,scale,aqcd,aqed,buff,
     $           u_syst,s_buff,nclus,buffclus,use_vert_displacement,buff_vert_displacement,done)
         else
            write(*,*) 'Error done early',j,nw
         endif
         if (store_event(j) .and. .not. done) then
            wgt=wgt*xscale
            wgt = dsign(max(dabs(wgt), target_wgt),wgt)
            if (dabs(wgt) .gt. target_wgt) then
               xover = xover + dabs(wgt) - target_wgt
               nover = nover+1
            endif
            xtot = xtot + dabs(wgt)
            i=i+1
            call write_event(lunw,p,wgt,n,ic,ngroup,scale,aqcd,aqed,
     $           buff,u_syst,s_buff,nclus,buffclus,use_vert_displacement,buff_vert_displacement)
         endif
      enddo
      write(*,*) 'Found ',nw,' events.'
      write(*,*) 'Wrote ',i ,' events.'
      if (scale_to_xsec)then
         write(*,*) 'Actual xsec ',xsec
         write(*,*) 'Correct abs xsec ',xsecabs
         write(*,*) 'Event xsec ', xtot
      endif
      write(*,*) 'Events wgts > 1: ', nover
      write(*,*) '% Cross section > 1: ',xover, xover/xtot*100.
      neventswritten = i
      maxwgt = target_wgt
      if (force_max_wgt.lt.0)then
         th_maxwgt = target_wgt
         th_nunwgt = neventswritten
      endif

 99   close(lunw)
      
c      close(lun)
      end

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
      double precision travel_dist, max_travel_distance
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

      logical use_vert_displacement
      data use_vert_displacement/.true./
      character*30 buff_vert_displacement(3)      
      
c      integer ncols,ncolflow(maxamps),ncolalt(maxamps)
c      common/to_colstats/ncols,ncolflow,ncolalt,ic
c      data ncolflow/maxamps*0/
c      data ncolalt/maxamps*0/

      include 'coupl.inc'

      include 'lhe_event_infos.inc'
      data AlreadySetInBiasModule/.False./

      include 'symswap.inc'

      include 'fit2D.inc'
      include '../../Source/fit2D_card.inc'
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
            rv(2)=ran1(iseed+1)
            call picktheta(E,rv,theta)

            rv(1)=ran1(iseed+2)
            rv(2)=ran1(iseed+2)
            if(cylinder) then
               phi = 2d0*pi*rv(1)
            elseif (parallelepiped) then
               call pickphi(theta,rv,phi)
            endif

            travel_dist = max_travel_distance(theta,dcos(phi),dsin(phi))
     &           *ran1(iseed+3)
            
c start debug ---            
c            write(230,*) theta,dcos(phi),dsin(phi)
c            write(*,*) max_travel_distance(theta,dcos(phi),dsin(phi))
c     &                  ,travel_dist
c end debug ---
            
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

c     write out buffer for interaction vertex placement 
      if(use_vert_displacement) then
         buff_vert_displacement(1) = '<vert_displacement>'
         write(buff_vert_displacement(2),*) travel_dist
         buff_vert_displacement(3) = '</vert_displacement>'
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
     &     sscale,aaqcd,aaqed,buff,use_syst,s_buff,nclus,buffclus,
     &     use_vert_displacement,buff_vert_displacement)
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
      double precision eps,Em,Ep,s
      common/celltable/cells,ncells
      double precision theta_min,theta_max
      logical fstcall
      data fstcall/.true./
      save fstcall,eps
      double precision pi
      parameter (pi=3.1415926d0)

      double precision fac_eps
      fac_eps= 1.5d0
      
c     At the first call, read and store the cell parameters from the 
c     cell_fortran.dat file
      if (fstcall) then
         open(unit=210,file='../cell_fortran.dat',status='old',
     *        err=999)
         open(unit=215,file='../gen_theta.stat',status='unknown')
         
c     The parameter eps play the role of the energy resolution.
c     Here, it is treated as a dimensional variable. Its value
c     should be a multiple of the smallest energy interval in the 2D
c     cell fit (this can be reasonably thought as the precision of the fit). 
c     This parameter is related to possible error in generation of the theta
c     angle when no cell corresponds to the given energy value.
c     The statistics about how many failures and the corresponding
c     problematic energy values are printed in the file "gen_theta.stat"
c     The user can set the multiplicative factor: raising it values
c     should result in a smaller number of failures.
         eps = 1d9
         ncells=1               !total number of cells
         theta_min = pi/2d0
         theta_max = 0d0
c     loop over infile lines until EoF is reached
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
                  if(cells(i,3).lt.eps) eps=cells(i,3)
               enddo
               eps= eps*fac_eps
               exit
            else
               cells(ncells,:)= a(:)
               ncells=ncells+1

c     update theta_min, theta_max
               if(a(2).lt.theta_min) theta_min = a(2)
               if(a(2)+a(4).gt.theta_max) theta_max = a(2)+a(4)

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
         write(215,*) 'Error has occurred with the energy value: ', E
c     in this situation, we generate flat in [theta_min,theta_max]
         theta = theta_min + rv(1)*(theta_max-theta_min)
         return         
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

      write(*,*) 'theta do not generate! Check the consistency of cell_fortran.dat grid!'
      call exit(-1)
      return
      
 999  write(*,*) 'Cannot open input data file: cell_fortran, exit!'
      call exit(-1)

      end

      subroutine pickphi(theta,rv,phi)
*     extracts an azimuthal angle (phi) keeping the correlations theta-phi
*     for the case of a rectangular hit surface. 
      implicit none
      include 'fit2D.inc'
      real * 8 rv(2),theta,phi,x_side_on2,y_side_on2
      real * 8 radius,phi0,phi1,phi_min,phi_max
      double precision pi
      parameter (pi=3.1415926d0)
      
      x_side_on2 = x_side/2d0
      y_side_on2 = y_side/2d0
      
      radius = d_target_detector*dtan(theta)
      
      if (y_side_on2 .gt. x_side_on2) then
         
         if(radius.le.x_side_on2) then 
            phi = 2d0*pi*rv(1)
         else if (radius.le.y_side_on2) then
            phi0 = dacos(x_side_on2/radius)
            if (rv(1).lt.0.5d0) then
               phi_min = phi0
               phi_max = pi-phi0
            else
               phi_min = pi+phi0
               phi_max = 2d0*pi-phi0
            endif
            phi = phi_min + rv(2)*(phi_max-phi_min)
         else
            phi0 = dacos(x_side_on2/radius)
            phi1 = dasin(y_side_on2/radius)
            if (rv(1) .lt.0.25d0) then
               phi_min = phi0 
               phi_max = phi1
            else if (rv(1) .lt. 0.5d0) then
               phi_min = pi-phi1
               phi_max = pi-phi0
            else if (rv(1) .lt. 0.75d0) then
               phi_min = pi+phi0
               phi_max = pi+phi1
            else 
               phi_min = 2d0*pi-phi1
               phi_max = 2d0*pi-phi0
            endif
            phi = phi_min + rv(2)*(phi_max-phi_min)
         endif
         
      else
         
         if(radius.le.y_side_on2) then 
            phi = 2d0*pi*rv(1)
         else if (radius.lt.x_side_on2) then
            phi0 = dasin(y_side_on2/radius)
            if (rv(1).lt.0.5d0) then
               phi_min = pi-phi0
               phi_max = pi+phi0
            else
               phi_min = -phi0
               phi_max = +phi0
            endif
            phi = phi_min + rv(2)*(phi_max-phi_min)
            if (phi.lt.0d0) phi = 2d0*pi+phi
         else
            phi0 = dacos(x_side_on2/radius)
            phi1 = dasin(y_side_on2/radius)
            if (rv(1) .lt.0.25d0) then
               phi_min = phi0 
               phi_max = phi1
            else if (rv(1) .lt. 0.5d0) then
               phi_min = pi-phi1
               phi_max = pi-phi0
            else if (rv(1) .lt. 0.75d0) then
               phi_min = pi+phi0
               phi_max = pi+phi1
            else 
               phi_min = 2d0*pi-phi1
               phi_max = 2d0*pi-phi0
            endif
            phi = phi_min + rv(2)*(phi_max-phi_min)
         endif
      endif
      
      return
      end


      real * 8 function heaviside(x)
      implicit none
      real * 8 x
      if (x.gt.0d0) then
         heaviside = 1d0
      else
         heaviside = 0d0
      endif
      return
      end

      
      real *8 function max_travel_distance(theta,cphi,sphi)
      implicit none
      include 'fit2D.inc'
      real * 8 theta,cphi,sphi
      real * 8 radius,tgth,z1,z2,x1,y1,x2,y2,in_z1,in_z2,x3,y3,z3
      real * 8 heaviside,theta_star,xmax,xmin,ymax,ymin 

      z1 = d_target_detector
      z2 = z1 + depth

c--- cylinder detector
      if(cylinder) then
         if (theta.gt.theta_max) then
            max_travel_distance =0d0
            return
         endif
         radius = z1*dtan(theta_max)
         theta_star = datan(radius/z2)
         if (theta.lt.theta_star) then
            max_travel_distance = depth/dcos(theta)
            return
         else
            max_travel_distance = z1*(dtan(theta_max)-dtan(theta))
     &                               /dsin(theta)
            return
         endif
      endif
      
c--- parallelepiped detector
      if(parallelepiped) then
         
         xmax = 0.5d0 * x_side 
         xmin = -0.5d0 * x_side 
      
         ymax = 0.5d0 * y_side 
         ymin = -0.5d0 * y_side 
      
         tgth = dtan(theta)
         x1 = z1*tgth*cphi
         y1 = z1*tgth*sphi      
      
         in_z1 = heaviside(x1-xmin)*heaviside(xmax-x1)
     &        *heaviside(y1-ymin)*heaviside(ymax-y1)
         
         if (in_z1 .eq. 0d0) then
            max_travel_distance =0d0
            return
         endif

         x2 = z2*tgth*cphi
         y2 = z2*tgth*sphi

         in_z2 = heaviside(x2-xmin)*heaviside(xmax-x2)
     &        *heaviside(y2-ymin)*heaviside(ymax-y2)

         if (in_z2 .ne. 0d0) then
            max_travel_distance = dsqrt((x2-x1)**2+(y2-y1)**2+depth**2)
            return 
         endif
      
         if (x2 .gt. xmax) then
            x3 =  xmax
            y3 =  xmax*sphi/cphi
            z3 =  xmax/tgth/cphi
            max_travel_distance=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            return 
         endif
         if (x2 .lt. xmin) then 
            x3 =  xmin          
            y3 =  xmin*sphi/cphi
            z3 =  xmin/tgth/cphi
            max_travel_distance=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            return 
         endif
         if (y2 .gt. ymax) then
            x3 =  ymax*cphi/sphi
            y3 =  ymax
            z3 =  ymax/tgth/sphi
            max_travel_distance=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            return 
         endif
         if (y2 .lt. ymin) then
            x3 =  ymin*cphi/sphi
            y3 =  ymin
            z3 =  ymin/tgth/sphi
            max_travel_distance=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
            return 
         endif
      endif

      end
      


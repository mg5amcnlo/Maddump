      logical function setclscales(p, keepq2bck)
c**************************************************
c     Calculate dynamic scales based on clustering
c     Also perform xqcut and xmtc cuts
c     keepq2bck allow to not reset the parameter q2bck
c**************************************************
      implicit none

      logical keepq2bck
      include 'message.inc'
      include 'genps.inc'
      include 'maxconfigs.inc'
      include 'nexternal.inc'
      include 'maxamps.inc'
      include 'cluster.inc'
      include 'run.inc'
      include 'coupl.inc'
      include 'run_config.inc'
C   
C   ARGUMENTS 
C   
      DOUBLE PRECISION P(0:3,NEXTERNAL)
C   global variables
C     Present process number
      INTEGER IMIRROR,IPROC
      COMMON/TO_MIRROR/IMIRROR, IPROC
C     ICONFIG has this config number
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
c     Common block for reweighting info
c     q2bck holds the central q2fact scales
      integer jlast(2)
      integer njetstore(lmaxconfigs),iqjetstore(nexternal-2,lmaxconfigs)
      real*8 q2bck(2)
      integer njets,iqjets(nexternal)
      common /to_rw/jlast,njetstore,iqjetstore,njets,iqjets,q2bck
      data njetstore/lmaxconfigs*-1/
      real*8 xptj,xptb,xpta,xptl,xmtc
      real*8 xetamin,xqcut,deltaeta
      common /to_specxpt/xptj,xptb,xpta,xptl,xmtc,xetamin,xqcut,deltaeta
      double precision stot,m1,m2
      common/to_stot/stot,m1,m2

C   local variables
      integer i, j, idi, idj, k
      real*8 PI
      parameter( PI = 3.14159265358979323846d0 )
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      double precision asref, pt2prev(n_max_cl),pt2min
      integer n, ibeam(2), iqcd(0:2)
      integer idfl, idmap(-nexternal:nexternal)
      integer ipart(2,n_max_cl)
      double precision xnow(2),etot
      integer jfirst(2),jcentral(2),nwarning
      logical qcdline(2),partonline(2)
      logical failed,first
      data first/.true./
      data nwarning/0/
      integer nqcd(lmaxconfigs)
      include 'config_nqcd.inc'

c     Variables for keeping track of jets
      logical goodjet(n_max_cl)
      integer fsnum(2),ida(2),imo,jcode
      logical chclusold,fail,increasecode
      save chclusold
      integer tmpindex

      logical isqcd,isjet,isparton,cluster,isjetvx
      integer ifsno
      double precision alphas
      external isqcd, isjet, isparton, cluster, isjetvx, alphas, ifsno
      setclscales=.true.

      if(ickkw.le.0.and.xqcut.le.0d0.and.q2fact(1).gt.0.and.scale.gt.0) then
         if(use_syst)then
            s_scale=scale
            n_qcd=nqcd(iconfig)
            n_alpsem=0
            do i=1,2
               n_pdfrw(i)=0
            enddo
            s_rwfact=1d0
         endif
      return
      endif
c   
c   Cluster the configuration
c   
      
c     First time, cluster according to this config and store jets
c     (following times, only accept configurations if the same partons
c      are flagged as jets)
      chclusold=chcluster
      if(njetstore(iconfig).eq.-1)then
         chcluster=.true.
      endif
 100  clustered = cluster(p(0,1))
      if(.not.clustered) then
         open(unit=26,file='../../../error',status='unknown',err=999)
         write(26,*) 'Error: Clustering failed in cluster.f.'
         write(*,*) 'Error: Clustering failed in cluster.f.'
         stop
 999     write(*,*) 'error'
         setclscales=.false.
         clustered = .false.
         return
      endif
c     Reset chcluster to run_card value
      chcluster=chclusold

      if (btest(mlevel,1)) then
        write(*,*)'setclscales: identified tree {'
        do i=1,nexternal-2
          write(*,*)'  ',i,': ',idacl(i,1),'(',ipdgcl(idacl(i,1),igraphs(1),iproc),')',
     $       '&',idacl(i,2),'(',ipdgcl(idacl(i,2),igraphs(1),iproc),')',
     $       ' -> ',imocl(i),'(',ipdgcl(imocl(i),igraphs(1),iproc),')',
     $       ', ptij = ',dsqrt(pt2ijcl(i))
          write(*,*)'   icluster(',i,')=',(icluster(j,i),j=1,4)
        enddo
        write(*,*)'  process: ',iproc
        write(*,*)'  graphs (',igraphs(0),'):',(igraphs(i),i=1,igraphs(0))
        write(*,*)'}'
        write(*,*)'iconfig is ',iconfig
      endif

C   If we have fixed factorization scale, for ickkw>0 means central
C   scale, i.e. last two scales (ren. scale for these vertices are
C   anyway already set by "scale" above)
      if(ickkw.gt.0) then
         if(fixed_fac_scale.and.first)then
            q2bck(1)=q2fact(1)
            q2bck(2)=q2fact(2)
            first=.false.
         else if(fixed_fac_scale) then
            q2fact(1)=q2bck(1)
            q2fact(2)=q2bck(2)
         endif
      endif

c   Preparing graph particle information (ipart, needed to keep track of
c   external particle clustering scales)

c   ipart gives the external particle number corresponding to the present
c   quark or gluon line. 
c   For t-channel lines, ipart(1) contains the connected beam. 
c   For s-channel lines, it depends if it is quark or gluon line:
c   For quark lines, ipart(2) is 0 and ipart(1) connects to the corresponding
c   final-state quark. For gluons, if it splits into two gluons, 
c   it connects to the hardest gluon. If it splits into qqbar, it ipart(1) is
c   the hardest and ipart(2) is the softest.

      do i=1,nexternal
         ipart(1,ishft(1,i-1))=i
         ipart(2,ishft(1,i-1))=0
      enddo
      do n=1,nexternal-3
        call ipartupdate(p,imocl(n),idacl(n,1),idacl(n,2),
     $       ipdgcl(1,igraphs(1),iproc),ipart)
      enddo

c     Prepare beam related variables for scale and jet determination
      do i=1,2
         ibeam(i)=ishft(1,i-1)
c        jfirst is first parton splitting on this side
         jfirst(i)=0
c        jlast is last parton on this side This means
c        the last cluster which is still QCD.
         jlast(i)=0
c        jcentral is the central scale vertex on this side. i.e it stops
c        when the T channel particles is not colored anymore.
         jcentral(i)=0
c        qcdline gives whether this IS line is QCD
         qcdline(i)=isqcd(ipdgcl(ibeam(i),igraphs(1),iproc))
c        partonline gives whether this IS line is parton (start out true for any QCD)
         partonline(i)=qcdline(i)
c        goodjet gives whether this cluster line is considered a jet
c        i.e. if all related/previous clustering are jet
         goodjet(ibeam(i))=partonline(i)
      enddo

      do i=3,nexternal
         j=ishft(1,i-1)
         goodjet(j)=isjet(ipdgcl(j,igraphs(1),iproc))
      enddo

c     Go through clusterings and set factorization scale points for use in dsig
c     as well as which FS particles count as jets (from jet vertices)
      do i=1,nexternal
         iqjets(i)=0
      enddo
      if (nexternal.eq.3) goto 10
c     jcode helps keep track of how many QCD/non-QCD flips we have gone through
      jcode=1
c     increasecode gives whether we should increase jcode at next vertex
      increasecode=.false.
      do n=1,nexternal-2
        do i=1,2 ! index of the child in the interaction
          do j=1,2 ! j index of the beam
            if(idacl(n,i).eq.ibeam(j))then
c             IS clustering
              ibeam(j)=imocl(n)
c             Determine which are beam particles based on n
              if(n.lt.nexternal-2) then
                 ida(i)=idacl(n,i)
                 ida(3-i)=idacl(n,3-i)
                 imo=imocl(n)
              else
                 ida(i)=idacl(n,i)
                 ida(3-i)=imocl(n)
                 imo=idacl(n,3-i)
              endif
c             
              if(partonline(j))then
c             If jfirst not set, set it
                 if(jfirst(j).eq.0) jfirst(j)=n
c             Stop fact scale where parton line stops
                 jlast(j)=n
                 partonline(j)=goodjet(ida(3-i)).and.
     $                isjet(ipdgcl(imo,igraphs(1),iproc))
              else if (jfirst(j).eq.0) then
                 jfirst(j) = n
                 goodjet(imo)=.false.
              else
                 goodjet(imo)=.false.
              endif
c             If not jet vertex, increase jcode. This is needed
c             e.g. in VBF if we pass over to the other side and hit
c             parton vertices again.
              if(.not.goodjet(ida(3-i)).or.
     $             .not.isjet(ipdgcl(ida(i),igraphs(1),iproc)).or.
     $             .not.isjet(ipdgcl(imo,igraphs(1),iproc))) then
                  jcode=jcode+1
                  increasecode=.true.
               else if(increasecode) then
                  jcode=jcode+1
                  increasecode=.false.
               endif
c             Consider t-channel jet radiations as jets only if
c             FS line is a jet line
              if(goodjet(ida(3-i))) then
                 if(partonline(j).or.
     $ ipdgcl(ida(3-i),igraphs(1),iproc).eq.21)then
c                   Need to include gluon to avoid soft singularity
                    iqjets(ipart(1,ida(3-i)))=1 ! 1 means for sure jet
                 else
                    iqjets(ipart(1,ida(3-i)))=jcode ! jcode means possible jet
                 endif
              endif
c             Trace QCD line through event
              if(qcdline(j))then
                 jcentral(j)=n
                 qcdline(j)=isqcd(ipdgcl(imo,igraphs(1),iproc))
              endif
            endif
          enddo
        enddo
        if (imocl(n).ne.ibeam(1).and.imocl(n).ne.ibeam(2)) then
c          FS clustering
c          Check QCD jet, take care so not a decay
           if(.not.isjetvx(imocl(n),idacl(n,1),idacl(n,2),
     $        ipdgcl(1,igraphs(1),iproc),ipart,n.eq.nexternal-2)) then
c          Remove non-gluon jets that lead up to non-jet vertices
           if(ipart(1,imocl(n)).gt.2)then ! ipart(1) set and not IS line
c          The ishft gives the FS particle corresponding to imocl
              if(ipdgcl(ishft(1,ipart(1,imocl(n))-1),igraphs(1),iproc).ne.21)then
                 iqjets(ipart(1,imocl(n)))=0
              else if (ipdgcl(imocl(n),igraphs(1),iproc).eq.21)then
c                special case for g > g h remove also the hardest gluon
                 iqjets(ipart(1,imocl(n)))=0
              endif
           endif
           if(ipart(2,imocl(n)).gt.2)then ! ipart(1) set and not IS line
c             The ishft gives the FS particle corresponding to imocl
              if(ipdgcl(ishft(1,ipart(2,imocl(n))-1),igraphs(1),iproc).ne.21.or.
     $                                   ipdgcl(imocl(n),igraphs(1),iproc).ne.21) then
c                 The second condition is to prevent the case of ggh where the gluon split in quark later.
c                 The first quark is already remove so we shouldn't remove this one.      
              iqjets(ipart(2,imocl(n)))=0
              endif
           endif
c          Set goodjet to false for mother
              goodjet(imocl(n))=.false.
              cycle
           endif

c          This is a jet vertex, so set jet flag for final-state jets
c          ifsno gives leg number if daughter is FS particle, otherwise 0
           fsnum(1)=ifsno(idacl(n,1),ipart)
           if(isjet(ipdgcl(idacl(n,1),igraphs(1),iproc)).and.
     $          fsnum(1).gt.0) then
              iqjets(fsnum(1))=1
           endif
           fsnum(1)=ifsno(idacl(n,2),ipart)
           if(isjet(ipdgcl(idacl(n,2),igraphs(1),iproc)).and.
     $          fsnum(1).gt.0) then
              iqjets(fsnum(1))=1
           endif
c          Flag mother as good jet if PDG is jet and both daughters are jets
           goodjet(imocl(n))=
     $          (isjet(ipdgcl(imocl(n),igraphs(1),iproc)).and.
     $          goodjet(idacl(n,1)).and.goodjet(idacl(n,2)))
        endif
      enddo

      if (btest(mlevel,4))then
         write(*,*) 'QCD jet status (before): ',(iqjets(i),i=3,nexternal)
      endif
c     Emissions with code 1 are always jets
c     Now take care of possible jets (i.e., with code > 1)
      if(.not. partonline(1).or..not.partonline(2))then
c       First reduce jcode by one if one remaining partonline
c       (in that case accept all jets with final jcode)
        if(partonline(1).or.partonline(2)) jcode=jcode-1
c       There parton emissions with code <= jcode are not jets
         do i=3,nexternal
            if(iqjets(i).gt.1.and.iqjets(i).le.jcode)then
               iqjets(i)=0
            endif
         enddo
      endif

 10   if(jfirst(1).le.0) jfirst(1)=jlast(1)
      if(jfirst(2).le.0) jfirst(2)=jlast(2)

      if (btest(mlevel,3))
     $     write(*,*) 'jfirst is ',jfirst(1),jfirst(2),
     $     ' jlast is ',jlast(1),jlast(2),
     $     ' and jcentral is ',jcentral(1),jcentral(2)

      if (btest(mlevel,3)) then
         write(*,'(a$)') 'QCD jets (final): '
         do i=3,nexternal
            if(iqjets(i).gt.0) write(*,'(i3$)') i
         enddo
         write(*,*)
      endif
      if(njetstore(iconfig).eq.-1) then
c     Store external jet numbers if first time
         njets=0
         do i=3,nexternal
            if(iqjets(i).gt.0)then
               njets=njets+1
               iqjetstore(njets,iconfig)=i
            endif
         enddo
         njetstore(iconfig)=njets
         if (btest(mlevel,4))
     $        write(*,*) 'Storing jets: ',(iqjetstore(i,iconfig),i=1,njets)
c     Recluster without requiring chcluster
         goto 100
      else
c     Otherwise, check that we have the right jets
c     if not, recluster according to iconfig
         fail=.false.
         njets=0
         do i=1,nexternal
            if(iqjets(i).gt.0)then
               njets=njets+1
c               if (iqjetstore(njets,iconfig).ne.i) fail=.true.
            endif
         enddo
         if(njets.ne.njetstore(iconfig)) fail=.true.
         if (fail) then
            if (igraphs(1).eq.iconfig) then
               open(unit=26,file='../../../error',status='unknown',err=999)
               write(*,*) 'Error: Failed despite same graph: ',iconfig
               write(*,*) 'Have jets (>0)',(iqjets(i),i=1,nexternal)
               write(*,*) 'Should be ',
     $              (iqjetstore(i,iconfig),i=1,njetstore(iconfig))
               write(26,*) 'Error: Failed despite same graph: ',iconfig,
     $              '. Have jets (>0)',(iqjets(i),i=1,nexternal),
     $              ', should be ',
     $              (iqjetstore(i,iconfig),i=1,njetstore(iconfig))
               stop
            endif
            if (btest(mlevel,3))
     $           write(*,*) 'Bad clustering, jets fail. Reclustering ',
     $           iconfig
            chcluster=.true.
            goto 100
         endif
      endif
      
c     If last clustering is s-channel QCD (e.g. ttbar) use mt2last instead
c     (i.e. geom. average of transverse mass of t and t~)
        if(mt2last.gt.4d0 .and. nexternal.gt.3) then
           if(jlast(1).eq.nexternal-2.and.jlast(2).eq.nexternal-2.and.
     $        isqcd(ipdgcl(idacl(nexternal-3,1),igraphs(1),iproc)).and.
     $        isqcd(ipdgcl(idacl(nexternal-3,2),igraphs(1),iproc)).and.
     $        isqcd(ipdgcl(imocl(nexternal-3),igraphs(1),iproc)))then
              mt2ij(nexternal-2)=mt2last
              mt2ij(nexternal-3)=mt2last
              if (btest(mlevel,3)) then
                 write(*,*)' setclscales: set last vertices to mtlast: ',sqrt(mt2last)
              endif
           endif
        endif

c     Set central scale to mT2
      if(jcentral(1).gt.0) then
         if(mt2ij(jcentral(1)).gt.0d0)
     $        pt2ijcl(jcentral(1))=mt2ij(jcentral(1))
      endif
      if(jcentral(2).gt.0)then
         if(mt2ij(jcentral(2)).gt.0d0)
     $     pt2ijcl(jcentral(2))=mt2ij(jcentral(2))
      endif
      if(btest(mlevel,4))then
         write(*,*) 'jlast, jcentral: ',(jlast(i),i=1,2),(jcentral(i),i=1,2)
         if(jlast(1).gt.0) write(*,*)'pt(jlast 1): ', sqrt(pt2ijcl(jlast(1)))
         if(jlast(2).gt.0) write(*,*)'pt(jlast 2): ', sqrt(pt2ijcl(jlast(2)))
         if(jcentral(1).gt.0) write(*,*)'pt(jcentral 1): ', sqrt(pt2ijcl(jcentral(1)))
         if(jcentral(2).gt.0) write(*,*)'pt(jcentral 2): ', sqrt(pt2ijcl(jcentral(2)))
      endif
c     Check xqcut for vertices with jet daughters only
      ibeam(1)=ishft(1,0)
      ibeam(2)=ishft(1,1)
      if(xqcut.gt.0) then
         do n=1,nexternal-3
c        Check if any of vertex daughters among jets
            do i=1,2
c              ifsno gives leg number if daughter is FS particle, otherwise 0
               fsnum(1)=ifsno(idacl(n,i),ipart)
               if(fsnum(1).gt.0)then
                  if(iqjets(fsnum(1)).gt.0)then
c                    Daughter among jets - check xqcut
                     if(sqrt(pt2ijcl(n)).lt.xqcut)then
                        if (btest(mlevel,3))
     $                       write(*,*) 'Failed xqcut: ',n,
     $                       ipdgcl(idacl(n,1),igraphs(1),iproc),
     $                       ipdgcl(idacl(n,2),igraphs(1),iproc),
     $                       sqrt(pt2ijcl(n))
                        setclscales=.false.
                        clustered = .false.
                        return
                     endif
                  endif
               endif
            enddo
         enddo
      endif
c     JA: Check xmtc cut for central process
      if(xmtc**2.gt.0) then
         if(jcentral(1).gt.0.and.pt2ijcl(jcentral(1)).lt.xmtc**2
     $      .or.jcentral(2).gt.0.and.pt2ijcl(jcentral(2)).lt.xmtc**2)then
            setclscales=.false.
            clustered = .false.
            if(btest(mlevel,3)) write(*,*)'Failed xmtc cut ',
     $           sqrt(pt2ijcl(jcentral(1))),sqrt(pt2ijcl(jcentral(1))),
     $           ' < ',xmtc
            return
         endif
      endif
      
      if(ickkw.eq.0.and.(fixed_fac_scale.or.q2fact(1).gt.0).and.
     $     (fixed_ren_scale.or.scale.gt.0)) return

c     Ensure that last scales are at least as big as first scales
      if(jlast(1).gt.0)
     $     pt2ijcl(jlast(1))=max(pt2ijcl(jlast(1)),pt2ijcl(jfirst(1)))
      if(jlast(2).gt.0)
     $     pt2ijcl(jlast(2))=max(pt2ijcl(jlast(2)),pt2ijcl(jfirst(2)))

      if(ickkw.gt.0.and.q2fact(1).gt.0) then
c     Use the fixed or previously set scale for central scale
         if(jcentral(1).gt.0) pt2ijcl(jcentral(1))=q2fact(1)
         if(jcentral(2).gt.0.and.jcentral(2).ne.jcentral(1))
     $        pt2ijcl(jcentral(2))=q2fact(2)
      endif

      if(nexternal.eq.3.and.nincoming.eq.2.and.q2fact(1).eq.0) then
         q2fact(1)=pt2ijcl(nexternal-2)
         q2fact(2)=pt2ijcl(nexternal-2)
      endif

      if(q2fact(1).eq.0d0) then
c     Use the geom. average of central scale and first non-radiation vertex
         if(jlast(1).gt.0) q2fact(1)=sqrt(pt2ijcl(jlast(1))*pt2ijcl(jcentral(1)))
         if(jlast(2).gt.0) q2fact(2)=sqrt(pt2ijcl(jlast(2))*pt2ijcl(jcentral(2)))
         if(jcentral(1).gt.0.and.jcentral(1).eq.jcentral(2))then
c     We have a qcd line going through the whole event, use single scale
            q2fact(1)=max(q2fact(1),q2fact(2))
            q2fact(2)=q2fact(1)
         endif
      endif
      if(.not. fixed_fac_scale) then
         q2fact(1)=scalefact**2*q2fact(1)
         q2fact(2)=scalefact**2*q2fact(2)
         if (.not.keepq2bck)then
            q2bck(1)=q2fact(1)
            q2bck(2)=q2fact(2)
         endif
         if (btest(mlevel,3))
     $      write(*,*) 'Set central fact scales to ',sqrt(q2bck(1)),sqrt(q2bck(2))
      endif
         
c     Set renormalization scale to geom. aver. of relevant scales
      if(scale.eq.0d0) then
         if(jlast(1).gt.0.and.jlast(2).gt.0)then
c           Use geom. average of last and central scales
            scale=(pt2ijcl(jlast(1))*pt2ijcl(jcentral(1))*
     $             pt2ijcl(jlast(2))*pt2ijcl(jcentral(2)))**0.125
         elseif(jlast(1).gt.0)then
c           Use geom. average of last and central scale
            scale=(pt2ijcl(jlast(1))*pt2ijcl(jcentral(1)))**0.25
         elseif(jlast(2).gt.0)then
c           Use geom. average of last and central scale
            scale=(pt2ijcl(jlast(2))*pt2ijcl(jcentral(2)))**0.25
         elseif(jcentral(1).gt.0.and.jcentral(2).gt.0) then
c           Use geom. average of central scales
            scale=(pt2ijcl(jcentral(1))*pt2ijcl(jcentral(2)))**0.25d0
         elseif(jcentral(1).gt.0) then
            scale=sqrt(pt2ijcl(jcentral(1)))
         elseif(jcentral(2).gt.0) then
            scale=sqrt(pt2ijcl(jcentral(2)))
         else
            scale=sqrt(pt2ijcl(nexternal-2))
         endif
         scale=scalefact*scale
         if(scale.gt.0)
     $        G = SQRT(4d0*PI*ALPHAS(scale))
      endif
      if (btest(mlevel,3))
     $     write(*,*) 'Set ren scale to ',scale


c     Take care of case when jcentral are zero
      if(jcentral(1).eq.0.and.jcentral(2).eq.0)then
         if(q2fact(1).gt.0)then
            pt2ijcl(nexternal-2)=q2fact(1)
            if(nexternal.gt.3) pt2ijcl(nexternal-3)=q2fact(1)
         else
            q2fact(1)=scalefact**2*pt2ijcl(nexternal-2)
            q2fact(2)=scalefact**2*q2fact(1)
         endif
      elseif(jcentral(1).eq.0)then
            q2fact(1) = scalefact**2*pt2ijcl(jfirst(1))
      elseif(jcentral(2).eq.0)then
            q2fact(2) = scalefact**2*pt2ijcl(jfirst(2))
      elseif(ickkw.eq.2.or.(pdfwgt.and.ickkw.gt.0))then
c     Total pdf weight is f1(x1,pt2E)*fj(x1*z,Q)/fj(x1*z,pt2E)
c     f1(x1,pt2E) is given by DSIG, just need to set scale.
c     Use the minimum scale found for fact scale in ME
         if(jlast(1).gt.0.and.jfirst(1).le.jlast(1))
     $        q2fact(1)=scalefact**2*min(pt2ijcl(jfirst(1)),q2fact(1))
         if(jlast(2).gt.0.and.jfirst(2).le.jlast(2))
     $        q2fact(2)=scalefact**2*min(pt2ijcl(jfirst(2)),q2fact(2))
      endif

      if (lpp(1) .eq. 9 .and. lpp(2) .eq. 0) then
         continue
c     Check that factorization scale is >= 2 GeV
      elseif(lpp(1).ne.0.and.q2fact(1).lt.4d0.or.
     $        lpp(2).ne.0.and.q2fact(2).lt.4d0)then
         if(nwarning.le.10) then
            nwarning=nwarning+1
            write(*,*) 'Warning: Too low fact scales: ',
     $           sqrt(q2fact(1)), sqrt(q2fact(2))
         endif
         if(nwarning.eq.11) then
            nwarning=nwarning+1
            write(*,*) 'No more warnings written out this run.'
         endif
         setclscales=.false.
         clustered = .false.
         return
      endif

      if (btest(mlevel,3))
     $     write(*,*) 'Set fact scales to ',sqrt(q2fact(1)),sqrt(q2fact(2))

c
c     Store jet info for matching
c
      etot=sqrt(stot)
      do i=1,nexternal
         ptclus(i)=0d0
      enddo

      do n=1,nexternal-2
         if(n.lt.nexternal-2) then
            ida(1)=idacl(n,1)
            ida(2)=idacl(n,2)
            imo=imocl(n)
         else
            ida(1)=idacl(n,1)
            ida(2)=imocl(n)
            imo=idacl(n,2)
         endif
         do i=1,2
            do j=1,2
c              First adjust goodjet based on iqjets
               if(goodjet(ida(i)).and.ipart(j,ida(i)).gt.2)then
                  if(iqjets(ipart(j,ida(i))).eq.0) goodjet(ida(i))=.false.
               endif
c              Now reset ptclus if jet vertex
               if(ipart(j,ida(i)).gt.2) then
                  if(isjetvx(imocl(n),idacl(n,1),idacl(n,2),
     $               ipdgcl(1,igraphs(1),iproc),ipart,n.eq.nexternal-2)
     $              .and.goodjet(ida(i))) then
                     ptclus(ipart(j,ida(i)))=
     $                    max(ptclus(ipart(j,ida(i))),dsqrt(pt2ijcl(n)))
                  else if(ptclus(ipart(j,ida(i))).eq.0d0) then
                     ptclus(ipart(j,ida(i)))=etot
                  endif
                  if (btest(mlevel,3))
     $                 write(*,*) 'Set ptclus for ',ipart(j,ida(i)),
     $                 ' to ', ptclus(ipart(j,ida(i))),ida(i),goodjet(ida(i))
               endif
            enddo
         enddo
      enddo
c
c     Store information for systematics studies
c

      if(use_syst)then
         s_scale=scale
         n_qcd=nqcd(igraphs(1))
         n_alpsem=0
         do i=1,2
            n_pdfrw(i)=0
         enddo
         s_rwfact=1d0
      endif
      return
      end

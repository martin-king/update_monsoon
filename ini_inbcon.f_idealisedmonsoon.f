 
      SUBROUTINE INBCON (grav0,radlat)
C--
C--   SUBROUTINE INBCON (grav0,radlat)
C--
C--   Purpose : Read topography and climatological boundary conditions 
C--   Input :   grav0  = gravity accel.
C--             radlat = grid latitudes in radiants
 									
      include "atparam.h"

      PARAMETER ( NLON=IX, NLAT=IL )

      include "com_tsteps.h" 
      include "com_cpl_flags.h"
 
      include "com_surfcon.h"    

      include "com_cli_land.h" 
      include "com_cli_sea.h" 

      real   radlat(il)

      real*4 r4inp(ix,il), dummy4
      real*4 veg(ix,il), swl1(ix,il), swl2(ix,il)

c mpkaqua for calculating gaussian thermal forcing
         real ymlon(ix), ymlat(il)

      iitest=1

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      thrsh = 0.1

C--   1. Read topographical fields (orography, land-sea mask)

      if (iitest.ge.1) print*,' read orography' 

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
c mpkaqua          phi0(i,j) = grav0*r4inp(i,j)
c mpkaqua just an initial perturbation; set flat again in ini_invars.f
c          phi0(i,j)=2.*grav0*abs(sin(dble(i)/dble(ix)*6.28))
           phi0(i,j)=2.*grav0
        enddo
      enddo
c mpk BEGIN putting in Gaussian mountain here; also used in specifying thermal forcing,
c mpk see below.
c mpk this is the Himalaya mountains
      tpi = 2.0*3.1415926536
      dtr = tpi/360.
      ymolon = 78.*dtr
      ymolat = 31.*dtr
      ymwlon = 15.*dtr
      ymwlat = 2.*dtr
      ymrlon = 10.*dtr
      ymrlat = 2.*dtr
      ymheight = 6000.d0

      do i = 1, ix
           ymlon(i) = float(i-1)/float(ix)*360.*dtr
      enddo
      do j = 1, il/2
           ymlat(j) = (float(j-1)/float(il)*180.-90.)*dtr
           ymlat(il-j+1) = -1.*ymlat(j)
      enddo
 
      do j = 1, il
            ymdy = abs(ymlat(j) - ymolat)         ! dist from y origin
            ymyy = max(0., ymdy-ymrlat)/ymwlat
            do i = 1, ix
              ymdx = abs(ymlon(i) - ymolon)       ! dist from x origin
              ymdx = min(ymdx, abs(ymdx-tpi))     ! To ensure that: -pi <= dx <= pi
              ymxx = max(0., ymdx-ymrlon)/ymwlon
              phi0(i,j) = phi0(i,j) + 
     &                       grav0*ymheight*exp(-ymxx**2-ymyy**2)
            enddo
       enddo

c mpk this is east Africa highlands
      tpi = 2.0*3.1415926536
      dtr = tpi/360.
      ymolon = 35.*dtr
      ymolat = 5.*dtr
      ymwlon = 2.*dtr
      ymwlat = 4.*dtr
      ymrlon = 2.*dtr
      ymrlat = 4.*dtr
      ymheight = 2000.d0

      do i = 1, ix
           ymlon(i) = float(i-1)/float(ix)*360.*dtr
      enddo
      do j = 1, il/2
           ymlat(j) = (float(j-1)/float(il)*180.-90.)*dtr
           ymlat(il-j+1) = -1.*ymlat(j)
      enddo
 
      do j = 1, il
            ymdy = abs(ymlat(j) - ymolat)         ! dist from y origin
            ymyy = max(0., ymdy-ymrlat)/ymwlat
            do i = 1, ix
              ymdx = abs(ymlon(i) - ymolon)       ! dist from x origin
              ymdx = min(ymdx, abs(ymdx-tpi))     ! To ensure that: -pi <= dx <= pi
              ymxx = max(0., ymdx-ymrlon)/ymwlon
              phi0(i,j) = phi0(i,j) + 
     &                       grav0*ymheight*exp(-ymxx**2-ymyy**2)
            enddo
       enddo

c mpk END

      call truncg (ntrun,phi0,phis0)
 
      if (iitest.ge.1) print*,' read fractional land-sea mask'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
c mpkaqua          fmask(i,j) = r4inp(i,j)
c mpk put a big polar continent in the NH
         if (j.ge.29 .and. j.lt.44) then       !asia and europe
           if(i.le.32)then
             fmask(i,j) = 1.d0
           elseif(i.gt.65.and.i.lt.78)then     !north america
             fmask(i,j) = 1.d0
           else
             fmask(i,j) = 0.d0
           endif
         elseif(j.ge.15 .and. j.le.29 .and. i.lt.13)then !africa
           fmask(i,j) = 1.d0
         elseif(j.ge.12 .and. j.le.29 .and. i.gt.75      !south america
     &          .and. i.lt.84)then
           fmask(i,j) = 1.d0
         elseif(j.ge.15 .and. j.le.21 .and. i.ge.32    !australia
     &           .and. i.lt.42) then
           fmask(i,j) = 1.d0
         else
           fmask(i,j) = 0.d0
         endif
         if(i.le.4 .and. j.ge.15 .and. j.le.25)then !add back gulf of guinea
           fmask(i,j) = 0.d0
         endif
c all aqua
c          fmask(i,j)=0.d0
c
        enddo
      enddo

C--   2. Initialize land-sfc boundary conditions

C--   2.1 Fractional and binary land masks

      do j=1,il
       do i=1,ix

         fmask_l(i,j) = fmask(i,j)

         if (fmask_l(i,j).ge.thrsh) then
           bmask_l(i,j) = 1.
           if (fmask(i,j).gt.(1.-thrsh)) fmask_l(i,j) = 1.
         else
           bmask_l(i,j) = 0.
           fmask_l(i,j) = 0.
         endif

         fmask1(i,j) = fmask_l(i,j)

       enddo
      enddo

C--   2.2 Annual-mean surface albedo

      if (iitest.ge.1) print*,' read surface albedo' 
 
      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
c mpkaqua          alb0(i,j) = 0.01*r4inp(i,j)
          alb0(i,j) = (1-fmask1(i,j))*0.07d0 + fmask1(i,j)*0.3d0
         if(j.ge.21 .and. j.le.27 .and. i.lt.13)then ! 10S-10N Africa
            alb0(i,j) = 0.1d0
         endif
        enddo
      enddo

C--   2.3 Land-surface temp.

      if (iitest.ge.1) print*,' reading land-surface temp.'
  
      do it = 1,12
        read (23) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            stl12(i,j,it) = r4inp(i,j)
c mpkaqua this will be corrected below
          enddo
        enddo
      enddo

      if (iitest.eq.1) print*,' checking land-surface temp.'

      CALL FORCHK (bmask_l,stl12,ix*il,12,0.,400.,273.)

C     Correction for model-to-actual topography
      do it = 1,12

        call ftland (stl12(1,1,it),phi0,phis0,bmask_l)

        if (iitest.gt.1) then
          do j = 1,il
            do i = 1,ix
              r4inp(i,j) = stl12(i,j,it)
            enddo
          enddo
          write (18) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        endif

      enddo

C--   2.4 Snow depth

      if (iitest.ge.1) print*,' reading snow depth'  

      do it = 1,12
        read (24) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
c mpkaqua            snowd12(i,j,it) = r4inp(i,j)
            snowd12(i,j,it) = 0.d0             
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking snow depth'

      CALL FORCHK (bmask_l,snowd12,ix*il,12,0.,20000.,0.)

C--  2.5 Read soil moisture and compute soil water availability 
C--      using vegetation fraction

      if (iitest.ge.1) print*,' reading soil moisture'  

c     read vegetation fraction (in %)
      read (25) ((veg(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
c mpkaqua          veg(i,j)=max(0.,0.01*veg(i,j))
            veg(i,j) = 0.d0
            if(j.ge.21 .and. j.le.27 .and. i.lt.13)then ! 10S-10N Africa
             veg(i,j) = 0.9d0 
            endif
        enddo
      enddo

      sdep1 = 70.
      idep2 = 3
      sdep2 = idep2*sdep1

      swwil2= sdep2*swwil
      rsw   = 1./(sdep1*swcap+sdep2*(swcap-swwil))

      do it = 1,12

        read (26) ((swl1(i,j),i=1,ix),j=il,1,-1)
        read (26) ((swl2(i,j),i=1,ix),j=il,1,-1)
        read (26) dummy4

        do j = 1,il
          do i = 1,ix
c mpkaqua            swroot = idep2*swl2(i,j)
c mpkaqua            soilw12(i,j,it) = min(1.,rsw*(swl1(i,j)+
c mpkaqua     &                        veg(i,j)*max(0.,swroot-swwil2)))		
            soilw12(i,j,it) = 0.d0
            if(j.ge.21 .and. j.le.27 .and. i.lt.13)then ! 10S-10N Africa
              swroot = idep2*18.d0
              soilw12(i,j,it) = min(1.,rsw*(15.d0+
     &                        veg(i,j)*max(0.,swroot-swwil2)))
            endif
          enddo
        enddo

      enddo

      if (iitest.ge.1) print*,' checking soil moisture'

      CALL FORCHK (bmask_l,soilw12,ix*il,12,0.,10.,0.)


C--   3. Initialize sea-sfc boundary conditions

C--   3.1 Fractional and binary sea masks

      do j=1,il
       do i=1,ix

         fmask_s(i,j) = 1.-fmask(i,j)

         if (fmask_s(i,j).ge.thrsh) then
           bmask_s(i,j) = 1.
           if (fmask_s(i,j).gt.(1.-thrsh)) fmask_s(i,j) = 1.
         else
           bmask_s(i,j) = 0.
           fmask_s(i,j) = 0.
         endif

       enddo
      enddo

C     Grid latitudes for sea-sfc. variables
      rad2deg = 90./asin(1.)
      do j=1,il
         deglat_s(j) = rad2deg*radlat(j)
      enddo

C--   3.2 SST 

      if (iitest.ge.1) print*,' reading sst' 

      do it = 1,12
        read (21) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            sst12(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo
c mpkaqua BEGIN: **********************************
      do imm = 1,6
        immp6 = imm+6
        do j = 1,il/2
          temp = 0.d0
          nlonsea = 0
          do i = 1,ix
c mpkaqua   for the SH
c mpkaqua   the if condition is to leave out undefined values designated land points 
            if (sst12(i,j,imm) .gt. 150.d0) then  
             nlonsea = nlonsea + 1 
             temp = temp + sst12(i,j,imm) 
            endif
            if (sst12(i,il+1-j,immp6) .gt. 150.d0) then 
             nlonsea = nlonsea + 1 
             temp = temp + sst12(i,il+1-j,immp6) 
            endif
            if (stl12(i,j,imm) .gt. 150.d0) then  
             nlonsea = nlonsea + 1 
             temp = temp + stl12(i,j,imm) 
            endif
            if (stl12(i,il+1-j,immp6) .gt. 150.d0) then 
             nlonsea = nlonsea + 1 
             temp = temp + stl12(i,il+1-j,immp6) 
            endif         
          enddo
          do i = 1, ix
            sst12(i,j,imm) = temp/float(nlonsea)
          enddo

          temp = 0.d0
          nlonsea = 0 
          do i = 1, ix
c mpkaqua   for the NH
            if (sst12(i,j,immp6) .gt. 150.d0) then 
              nlonsea = nlonsea + 1 
              temp = temp + sst12(i,j,immp6) 
            endif 
            if (sst12(i,il+1-j,imm) .gt. 150.d0) then
              nlonsea = nlonsea + 1
              temp = temp + sst12(i,il+1-j,imm) 
            endif 
            if (stl12(i,j,immp6) .gt. 150.d0) then 
              nlonsea = nlonsea + 1 
              temp = temp + stl12(i,j,immp6) 
            endif 
            if (stl12(i,il+1-j,imm) .gt. 150.d0) then
              nlonsea = nlonsea + 1
              temp = temp + stl12(i,il+1-j,imm) 
            endif 
          enddo 
          do i = 1, ix   
           sst12(i,il+1-j,imm) = temp/float(nlonsea)
          enddo  

        enddo
      enddo
      
      do imm=7,12
       do j=1, il/2
        do i=1, ix
         sst12(i,j,imm)=sst12(i,il+1-j,imm-6)
         sst12(i,il+1-j,imm)=sst12(i,j,imm-6)
        enddo
       enddo
      enddo

c        height -->   ___________________________
c                   /                           \
c                  /              |              \
c    gaussian     /               |               \
c      sides --> /                |                \
c               /               olon                \
c         _____/                olat                 \______
c
c              |    |             |
c              |<-->|<----------->|
c              |wlon|    rlon     |
c               wlat     rlat
   
         tpi = 2.0*3.1415926536
         dtr = tpi/360.
         ymolon = 120.*dtr
         ymolat = 10.*dtr
         ymrlon = 15.*dtr
         ymrlat = 5.*dtr
         ymwlon = 20.*dtr
         ymwlat = 5.*dtr
         ymheight = 0.d0
         do i = 1, ix
           ymlon(i) = float(i-1)/float(ix)*360.*dtr
         enddo
         do j = 1, il/2
           ymlat(j) = (float(j-1)/float(il)*180.-90.)*dtr
           ymlat(il-j+1) = -1.*ymlat(j)
         enddo
         do j = 1, il
            ymdy = abs(ymlat(j) - ymolat)         ! dist from y origin
            ymyy = max(0., ymdy-ymrlat)/ymwlat
            do i = 1, ix
              ymdx = abs(ymlon(i) - ymolon)       ! dist from x origin
              ymdx = min(ymdx, abs(ymdx-tpi))     ! To ensure that: -pi <= dx <= pi
              ymxx = max(0., ymdx-ymrlon)/ymwlon
              sst12(i,j,6) = sst12(i,j,6) + 
     &                       ymheight*exp(-ymxx**2-ymyy**2)
              sst12(i,j,7) = sst12(i,j,7) + 
     &                       ymheight*exp(-ymxx**2-ymyy**2)
              sst12(i,j,8) = sst12(i,j,8) + 
     &                       ymheight*exp(-ymxx**2-ymyy**2)
            enddo
         enddo


c mpkaqua the following is to make sure it is always running the surface 
c mpkaqua climate calculated above 
         
         do it = 1, 12 
          do j = 1, il/2
           do i = 1, ix
c mpkaqua this should not be necessary since fmask is set to sea everywhere. 
c mpkaqua but just in case...
            stl12(i,j,it) = sst12(i,j,it)
            stl12(i,il+1-j,it) = sst12(i,il+1-j,it)
           enddo
          enddo
         enddo

c mpk temperature corrected with topography
        do it = 1,12
         do j = 1, il
          do i = 1, ix
           stl12(i,j,it) = stl12(i,j,it) - 0.0065d0 * phis0(i,j) / grav0
          enddo
         enddo
        enddo

c mpkaqua END: ********************************

      if (iitest.ge.1) print*,' checking sst'

      CALL FORCHK (bmask_s,sst12,ix*il,12,100.,400.,273.)

c     3.2 Sea ice fraction

      if (iitest.ge.1) print*,' reading sea ice'  

      do it = 1,12
        read (22) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
c mpkaqua             
c             sice12(i,j,it) = r4inp(i,j)
             if (sst12(i,j,it) .LT. 270.
     &        .and. bmask_s(i,j).eq.1) then 
c mpkaqua too severe?             oice12(i,j,it) = 1.d0
              sice12(i,j,it) = 0.5d0
             else
              sice12(i,j,it) = 0.d0
             endif    
             if (stl12(i,j,it).LT.270. 
     &        .and. bmask_s(i,j).eq.0)then 
c mpkaqua too severe?             oice12(i,j,it) = 1.d0
              snowd12(i,j,it) = 100.
             else
              snowd12(i,j,it) = 0.d0
             endif       
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sea ice'

      CALL FORCHK (bmask_s,sice12,ix*il,12,0.,1.,0.)

C--   3.3 SST anomalies for initial and prec./following months

      if (isstan.gt.0) then

        if (iitest.ge.1) print*,' reading sst anomalies' 

        do jrec=1,isst0-2
          read (30) dummy4
        enddo

        do it=1,3

c         NB If isst0 = 1, SST_an(it=1) is set = SST_an(it=2)  
          if (it.ne.2.or.isst0.gt.1)
     &       read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

          do j = 1,il
            do i = 1,ix
              sstan3(i,j,it) = r4inp(i,j)  
            enddo
          enddo

        enddo

        if (iitest.ge.1) print*,' checking sst anomalies'

        CALL FORCHK (bmask_s,sstan3,ix*il,3,-50.,50.,0.)

      endif

C--   3.4. Annual-mean heat flux into sea-surface

      do j = 1,il
        do i = 1,ix
          hfseacl(i,j) = 0.
        enddo
      enddo

      if (icsea.ge.1) then

        if (iitest.ge.1) print*,' reading sfc heat fluxes' 

        irecl = 4*ix*il
        irec = 0

        open ( unit=31, file='fort.31', status='old', 
     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,12

          irec=irec+2
          read (31,rec=irec) r4inp

          do j = 1,il
            do i = 1,ix
              hfseacl(i,j) = hfseacl(i,j)+r4inp(i,j)
            enddo
          enddo  

        enddo

        do j = 1,il
          do i = 1,ix
            if (bmask_s(i,j).gt.0.) then
                hfseacl(i,j) = hfseacl(i,j)/(12.*fmask_s(i,j))
            else
                hfseacl(i,j) = 0.
            endif
          enddo
        enddo  

        if (iitest.ge.1) print*,' checking sfc heat fluxes'

        CALL FORCHK (bmask_s,hfseacl,ix*il,1,-1000.,1000.,0.)

      endif

C--   3.5. Ocean model SST climatology:
C--        defined by adding SST model bias to obs. climatology
C--        (bias may be defined in a different period from climatology)

      if (icsea.ge.3) then

        if (iitest.ge.1) print*,' reading ocean model SST bias' 

c        irecl = 4*ix*il
c        irec = 0

c        open ( unit=32, file='fort.32', status='old', 
c     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,12

c          irec=irec+1
c          read (32,rec=irec) r4inp
           read (32) r4inp

          do j = 1,il
            do i = 1,ix
              sstom12(i,j,it) = sst12(i,j,it)+r4inp(i,j)
            enddo
          enddo  

        enddo

        if (iitest.ge.1) print*,' checking ocean model SST'

        CALL FORCHK (bmask_s,sstom12,ix*il,12,100.,400.,273.)

      endif
C--
      RETURN
      END

      SUBROUTINE FORCHK (FMASK,FIELD,NGP,NF,FMIN,FMAX,FSET)
										
C--   Aux. routine FORCHK: Check consistency of sfc fields with land-sea mask 
C--   and set undefined values to a constant (to avoid over/underflow)

      real fmask(ngp), field(ngp,nf)

      do jf = 1,nf

        nfault=0

        do jgp = 1,ngp
          if (fmask(jgp).gt.0.0) then
            if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax)
     *          nfault = nfault+1
          else
            field(jgp,jf) = fset
          endif
        enddo

        print *, ' field: ', jf, '   no. of faulty points:', nfault

      enddo

      print *, ' undefined values set to', fset

      RETURN
      END

      SUBROUTINE FTLAND (STL,PHI0,PHIS0,FMASKL)

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL)

      include "com_dyncon0.h" 
      include "com_dyncon1.h" 

      REAL STL(NLON,NLAT), PHI0(NLON,NLAT), PHIS0(NLON,NLAT),
     &     FMASKL(NLON,NLAT)

      REAL STL2(NLON,NLAT)

      NL8 = NLAT/8
      GAM = 0.001*GAMMA/GRAV

      NLAT1 = 1
      NLAT2 = NL8

      DO JBAND=1,8

        SUMT=0.
        SUMW=0.

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           STL(I,J)=STL(I,J)+GAM*PHI0(I,J)
           SUMT=SUMT+GCOS(J)*FMASKL(I,J)*STL(I,J)
           SUMW=SUMW+GCOS(J)*FMASKL(I,J)
         ENDDO
        ENDDO

        SUMT=SUMT/SUMW

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=SUMT
         ENDDO
        ENDDO
  
        NLAT1=NLAT1+NL8
        NLAT2=NLAT2+NL8

      ENDDO

      ITR=7
      IDTR=(NTRUN-6)/3

      DO JFIL=1,4

        CALL TRUNCG (ITR,STL,STL2)

        DO J=1,NLAT
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=STL2(I,J)
         ENDDO
        ENDDO

        ITR=MIN(ITR+IDTR,NTRUN)

      ENDDO

      CALL TRUNCG (ITR,STL,STL2)

      DO J=1,NLAT
       DO I=1,NLON
         STL(I,J)=STL2(I,J)-GAM*PHIS0(I,J)
       ENDDO
      ENDDO       

      RETURN
      END

      SUBROUTINE TRUNCG (ITR,FG1,FG2)

C--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
C--   Purpose : compute a spectrally-filtered grid-point field
C--   Input   : ITR : spectral truncation (triangular)
C--           : FG1 : original grid-point field
C--   Output  : FG2 : filtered grid-point field

      include "atparam.h"

      REAL FG1 (IX,IL), FG2(IX,IL)
      COMPLEX FSP(MX,NX), ZERO 

      ZERO = (0.,0.)

      CALL SPEC (FG1,FSP)

      DO N=1,NX
        DO M=1,MX
          ITWN=ISC*(M-1)+N-1
          IF (ITWN.GT.ITR) FSP(M,N)=ZERO
        ENDDO
      ENDDO

      CALL GRID (FSP,FG2,1)

      RETURN
      END

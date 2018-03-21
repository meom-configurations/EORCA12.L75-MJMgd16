MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of mbathy and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.6  ! 2015-05  (P. Mathiot) ISF: add wmask,wumask and wvmask
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk        : compute land/ocean mask
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         !
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
!{ DRAKKAR
   USE iom             ! For shlat2d
   USE fldread         ! for sn_shlat2d
!}

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk    ! routine called by inidom.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
   REAL(wp)        :: rn_shlat   ! type of lateral boundary condition on velocity
   LOGICAL, PUBLIC :: ln_vorlat  !  consistency of vorticity boundary condition 
   !                                            with analytical eqs.

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: dommsk.F90 4624 2014-04-28 12:09:03Z acc $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   INTEGER FUNCTION dom_msk_alloc()
      !!---------------------------------------------------------------------
      !!                 ***  FUNCTION dom_msk_alloc  ***
      !!---------------------------------------------------------------------
      dom_msk_alloc = 0
#if defined key_noslip_accurate
      ALLOCATE(icoord(jpi*jpj*jpk,3), STAT=dom_msk_alloc)
#endif
      IF( dom_msk_alloc /= 0 )   CALL ctl_warn('dom_msk_alloc: failed to allocate icoord array')
      !
   END FUNCTION dom_msk_alloc


   SUBROUTINE dom_msk
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zonal velocity points (u & v), vorticity points (f) points.
      !!
      !! ** Method  :   The ocean/land mask is computed from the basin bathy-
      !!      metry in level (mbathy) which is defined or read in dommba.
      !!      mbathy equals 0 over continental T-point 
      !!      and the number of ocean level over the ocean.
      !!
      !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
      !!      t-point : 0. IF mbathy( ji ,jj) =< 0
      !!                1. IF mbathy( ji ,jj) >= jk
      !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
      !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
      !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
      !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
      !!                  and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
      !!      tmask_i : interior ocean mask at t-point, i.e. excluding duplicated
      !!                rows/lines due to cyclic or North Fold boundaries as well
      !!                as MPP halos.
      !!
      !!        The lateral friction is set through the value of fmask along
      !!      the coast and topography. This value is defined by rn_shlat, a
      !!      namelist parameter:
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
      !!      are defined with the proper value at lateral domain boundaries.
      !!
      !!      In case of open boundaries (lk_bdy=T):
      !!        - tmask is set to 1 on the points to be computed bay the open
      !!          boundaries routines.
      !!
      !! ** Action :   tmask    : land/ocean mask at t-point (=0. or 1.)
      !!               umask    : land/ocean mask at u-point (=0. or 1.)
      !!               vmask    : land/ocean mask at v-point (=0. or 1.)
      !!               fmask    : land/ocean mask at f-point (=0. or 1.)
      !!                          =rn_shlat along lateral boundaries
      !!               tmask_i  : interior ocean mask
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      INTEGER  ::   iif, iil, ii0, ii1, ii   ! local integers
      INTEGER  ::   ijf, ijl, ij0, ij1       !   -       -
      INTEGER  ::   ios
      INTEGER  ::   isrow                    ! index for ORCA1 starting row
      INTEGER , POINTER, DIMENSION(:,:) ::  imsk
      REAL(wp), POINTER, DIMENSION(:,:) ::  zwf

!{ DRAKKAR
      INTEGER  :: inum             !  logical unit for shlat2d
      REAL(wp) :: zshlat           !: locally modified shlat for some strait
      REAL(wp), POINTER, DIMENSION(:,:) :: zshlat2d
      LOGICAL  :: ln_shlat2d
      TYPE(FLD_N) :: sn_shlat2d
      !!
      NAMELIST/namlbc/ rn_shlat, ln_vorlat, ln_shlat2d, sn_shlat2d
!}
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_msk')
      !
      CALL wrk_alloc( jpi, jpj, imsk )
      CALL wrk_alloc( jpi, jpj, zwf  )
      !
      REWIND( numnam_ref )              ! Namelist namlbc in reference namelist : Lateral momentum boundary condition
      READ  ( numnam_ref, namlbc, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namlbc in configuration namelist : Lateral momentum boundary condition
      READ  ( numnam_cfg, namlbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namlbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namlbc )
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat 
      ENDIF

      IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  free-slip '
      ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  no-slip '
      ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  partial-slip '
      ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ocean lateral  strong-slip '
      ELSE
         WRITE(ctmp1,*) ' rn_shlat is negative = ', rn_shlat
         CALL ctl_stop( ctmp1 )
      ENDIF

!{ DRAKKAR
      IF ( ln_shlat2d ) THEN
         IF(lwp) WRITE(numout,*) '         READ shlat as a 2D coefficient in a file '
         CALL wrk_alloc( jpi, jpj, zshlat2d  )
         CALL iom_open(sn_shlat2d%clname, inum)
         CALL iom_get (inum, jpdom_data, sn_shlat2d%clvar, zshlat2d, 1) !
         CALL iom_close(inum)
      ENDIF
!}

      ! 1. Ocean/land mask at t-point (computed from mbathy)
      ! -----------------------------
      ! N.B. tmask has already the right boundary conditions since mbathy is ok
      !
      tmask(:,:,:) = 0._wp
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( REAL( mbathy(ji,jj) - jk, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
            END DO  
         END DO  
      END DO  
      
      ! (ISF) define barotropic mask and mask the ice shelf point
      ssmask(:,:)=tmask(:,:,1) ! at this stage ice shelf is not masked
      
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( REAL( misfdep(ji,jj) - jk, wp ) - 0.1_wp >= 0._wp )   THEN
                  tmask(ji,jj,jk) = 0._wp
               END IF
            END DO  
         END DO  
      END DO  

      ! Interior domain mask (used for global sum)
      ! --------------------
      tmask_i(:,:) = ssmask(:,:)            ! (ISH) tmask_i = 1 even on the ice shelf

      tmask_h(:,:) = 1._wp                 ! 0 on the halo and 1 elsewhere
      iif = jpreci                         ! ???
      iil = nlci - jpreci + 1
      ijf = jprecj                         ! ???
      ijl = nlcj - jprecj + 1

      tmask_h( 1 :iif,   :   ) = 0._wp      ! first columns
      tmask_h(iil:jpi,   :   ) = 0._wp      ! last  columns (including mpp extra columns)
      tmask_h(   :   , 1 :ijf) = 0._wp      ! first rows
      tmask_h(   :   ,ijl:jpj) = 0._wp      ! last  rows (including mpp extra rows)

      ! north fold mask
      ! ---------------
      tpol(1:jpiglo) = 1._wp 
      fpol(1:jpiglo) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0._wp
         fpol(     1    :jpiglo) = 0._wp
         IF( mjg(nlej) == jpjglo ) THEN                  ! only half of the nlcj-1 row
            DO ji = iif+1, iil-1
               tmask_h(ji,nlej-1) = tmask_h(ji,nlej-1) * tpol(mig(ji))
            END DO
         ENDIF
      ENDIF
     
      tmask_i(:,:) = tmask_i(:,:) * tmask_h(:,:)

      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0._wp
         fpol(jpiglo/2+1:jpiglo) = 0._wp
      ENDIF

      ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)
      ! -------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
            END DO
            DO ji = 1, jpim1      ! NO vector opt.
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO
      ! (ISF) MIN(1,SUM(umask)) is here to check if you have effectively at least 1 wet cell at u point
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector loop
            ssumask(ji,jj)  = ssmask(ji,jj) * ssmask(ji+1,jj  )  * MIN(1._wp,SUM(umask(ji,jj,:)))
            ssvmask(ji,jj)  = ssmask(ji,jj) * ssmask(ji  ,jj+1)  * MIN(1._wp,SUM(vmask(ji,jj,:)))
         END DO
         DO ji = 1, jpim1      ! NO vector opt.
            ssfmask(ji,jj) =  ssmask(ji,jj  ) * ssmask(ji+1,jj  )   &
               &            * ssmask(ji,jj+1) * ssmask(ji+1,jj+1) * MIN(1._wp,SUM(fmask(ji,jj,:)))
         END DO
      END DO

!{ DRAKKAR
      IF( cp_cfg == "orca" .AND. jp_cfg == 25 ) THEN  ! ORCA R025 configuration
          !                                            !  =======================
          ii0 = 212    ;   ii1 = 212        ! East of Ombai strait
          ij0 = 464    ;   ij1 = 465   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the East Ombai Strait'
          ii0 = 210    ;   ii1 = 211        ! West of Ombai strait
          ij0 = 466    ;   ij1 = 466   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the West Ombai Strait '
          ii0 = 210    ;   ii1 = 210        ! exit of Ombai strait
          ij0 = 464    ;   ij1 = 465   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the exit of Ombai Strait '
          ii0 = 172    ;   ii1 = 175        ! Lombok strait
          ij0 = 463    ;   ij1 = 463   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the Lombok Strait'
          ii0 = 278    ;   ii1 = 279        ! Torres Strait
          ij0 = 456    ;   ij1 = 461   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  4.0
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 4 at the Torres Strait'
 
          !
       ENDIF
!}
      CALL lbc_lnk( umask  , 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( vmask  , 'V', 1._wp )
      CALL lbc_lnk( fmask  , 'F', 1._wp )
      CALL lbc_lnk( ssumask, 'U', 1._wp )      ! Lateral boundary conditions
      CALL lbc_lnk( ssvmask, 'V', 1._wp )
      CALL lbc_lnk( ssfmask, 'F', 1._wp )

      ! 3. Ocean/land mask at wu-, wv- and w points 
      !----------------------------------------------
      wmask (:,:,1) = tmask(:,:,1)     ! surface
      wumask(:,:,1) = umask(:,:,1)
      wvmask(:,:,1) = vmask(:,:,1)
      DO jk = 2, jpk                   ! interior values
         wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
         wumask(:,:,jk) = umask(:,:,jk) * umask(:,:,jk-1)   
         wvmask(:,:,jk) = vmask(:,:,jk) * vmask(:,:,jk-1)
      END DO

      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------     
      DO jk = 1, jpk
         zwf(:,:) = fmask(:,:,jk)         
!{ DRAKKAR
         IF (  ln_shlat2d ) THEN
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fmask(ji,jj,jk) == 0. ) THEN
                     fmask(ji,jj,jk) = zshlat2d(ji,jj) * MIN( 1._wp, MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                        &                                              zwf(ji-1,jj), zwf(ji,jj-1)  )  )
                  ENDIF
               END DO
            END DO

         ELSE
!}
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fmask(ji,jj,jk) == 0._wp ) THEN
                     fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                       &                                           zwf(ji-1,jj), zwf(ji,jj-1)  )  )
                  ENDIF
               END DO
            END DO
         ENDIF
!}

!{ DRAKKAR
       ! Locally modify shlat :
       IF( cp_cfg == "orca" .AND. jp_cfg == 025 ) THEN
         !                                           ! =======================
         ! Increased lateral friction in             !  ORCA_R025 configuration
         ! the vicinity of some straits              ! =======================
         !
            !! Gibraltar strait and Gulf of Cadiz
            ij0 = 652    ;   ij1 = 654
            ii0 = 1125   ;   ii1 = 1127
            zshlat=3
         DO jj = mj0(ij0),mj1(ij1)
            DO ji = mi0(ii0), mi1(ii1)
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = zshlat * MIN( 1., MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                     &                                     zwf(ji-1,jj), zwf(ji,jj-1)  )  )
               ENDIF
            END DO
         END DO
         !
       ENDIF
!}
         DO jj = 2, jpjm1
            IF( fmask(1,jj,jk) == 0._wp ) THEN
               fmask(1  ,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
            ENDIF
            IF( fmask(jpi,jj,jk) == 0._wp ) THEN
               fmask(jpi,jj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
            ENDIF
         END DO         
         DO ji = 2, jpim1
            IF( fmask(ji,1,jk) == 0._wp ) THEN
               fmask(ji, 1 ,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
            ENDIF
            IF( fmask(ji,jpj,jk) == 0._wp ) THEN
               fmask(ji,jpj,jk) = rn_shlat * MIN( 1._wp , MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
            ENDIF
         END DO
      END DO
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN   ! ORCA_R2 configuration
         !                                                 ! Increased lateral friction near of some straits
         !                                ! Gibraltar strait  : partial slip (fmask=0.5)
         ij0 = 101   ;   ij1 = 101
         ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
         ij0 = 102   ;   ij1 = 102
         ii0 = 139   ;   ii1 = 140   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.5_wp
         !
         !                                ! Bab el Mandeb : partial slip (fmask=1)
         ij0 =  87   ;   ij1 =  88
         ii0 = 160   ;   ii1 = 160   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
         ij0 =  88   ;   ij1 =  88
         ii0 = 159   ;   ii1 = 159   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  1._wp
         !
         !                                ! Danish straits  : strong slip (fmask > 2)
! We keep this as an example but it is instable in this case 
!         ij0 = 115   ;   ij1 = 115
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
!         ij0 = 116   ;   ij1 = 116
!         ii0 = 145   ;   ii1 = 146   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 4._wp
         !
      ENDIF
      !
      IF( cp_cfg == "orca" .AND. jp_cfg == 1 ) THEN   ! ORCA R1 configuration
         !                                                 ! Increased lateral friction near of some straits
         ! This dirty section will be suppressed by simplification process:
         ! all this will come back in input files
         ! Currently these hard-wired indices relate to configuration with
         ! extend grid (jpjglo=332)
         !
         isrow = 332 - jpjglo
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   orca_r1: increase friction near the following straits : '
         IF(lwp) WRITE(numout,*) '      Gibraltar '
         ii0 = 282           ;   ii1 = 283        ! Gibraltar Strait 
         ij0 = 241 - isrow   ;   ij1 = 241 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Bhosporus '
         ii0 = 314           ;   ii1 = 315        ! Bhosporus Strait 
         ij0 = 248 - isrow   ;   ij1 = 248 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Makassar (Top) '
         ii0 =  48           ;   ii1 =  48        ! Makassar Strait (Top) 
         ij0 = 189 - isrow   ;   ij1 = 190 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  

         IF(lwp) WRITE(numout,*) '      Lombok '
         ii0 =  44           ;   ii1 =  44        ! Lombok Strait 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Ombai '
         ii0 =  53           ;   ii1 =  53        ! Ombai Strait 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      Timor Passage '
         ii0 =  56           ;   ii1 =  56        ! Timor Passage 
         ij0 = 164 - isrow   ;   ij1 = 165 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 2._wp  

         IF(lwp) WRITE(numout,*) '      West Halmahera '
         ii0 =  58           ;   ii1 =  58        ! West Halmahera Strait 
         ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  

         IF(lwp) WRITE(numout,*) '      East Halmahera '
         ii0 =  55           ;   ii1 =  55        ! East Halmahera Strait 
         ij0 = 181 - isrow   ;   ij1 = 182 - isrow   ;   fmask( mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1:jpk ) = 3._wp  
         !
      ENDIF
      !
      CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
      !
      ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 )
      !
      CALL wrk_dealloc( jpi, jpj, imsk )
      CALL wrk_dealloc( jpi, jpj, zwf  )
!{ DRAKKAR
      IF ( ln_shlat2d ) THEN
        CALL wrk_dealloc( jpi, jpj, zshlat2d  )
      ENDIF
!}
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_msk')
      !
   END SUBROUTINE dom_msk
   
   !!======================================================================
END MODULE dommsk

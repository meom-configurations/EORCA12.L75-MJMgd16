MODULE ldfdyn
   !!======================================================================
   !!                       ***  MODULE  ldfdyn  ***
   !! Ocean physics:  lateral viscosity coefficient 
   !!=====================================================================
   !! History :  OPA  ! 1997-07  (G. Madec)  multi dimensional coefficients
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_dyn_init  : initialization, namelist read, and parameters control
   !!   ldf_dyn       : update lateral eddy viscosity coefficients at each time step 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE ldfc1d_c2d      ! lateral diffusion: 1D and 2D cases
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module for ehanced bottom friction file
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_dyn_init   ! called by nemogcm.F90
   PUBLIC   ldf_dyn        ! called by step.F90

   !                                                !!* Namelist namdyn_ldf : lateral mixing on momentum *
   LOGICAL , PUBLIC ::   ln_dynldf_lap   !: laplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_blp   !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_lev   !: iso-level direction
   LOGICAL , PUBLIC ::   ln_dynldf_hor   !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_dynldf_iso   !: iso-neutral direction
   INTEGER , PUBLIC ::   nn_ahm_ijk_t    !: choice of time & space variations of the lateral eddy viscosity coef.
   REAL(wp), PUBLIC ::   rn_ahm_0        !: lateral laplacian eddy viscosity            [m2/s]
   REAL(wp), PUBLIC ::   rn_ahm_b        !: lateral laplacian background eddy viscosity [m2/s]
   REAL(wp), PUBLIC ::   rn_bhm_0        !: lateral bilaplacian eddy viscosity          [m4/s]

   LOGICAL , PUBLIC ::   l_ldfdyn_time   !: flag for time variation of the lateral eddy viscosity coef.

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahmt, ahmf   !: eddy diffusivity coef. at U- and V-points   [m2/s or m4/s]

   REAL(wp) ::   r1_12   = 1._wp / 12._wp    ! =1/12
   REAL(wp) ::   r1_4    = 0.25_wp           ! =1/4
   REAL(wp) ::   r1_288  = 1._wp / 288._wp   ! =1/( 12^2 * 2 )

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: ldfdyn.F90 6140 2015-12-21 11:35:23Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_init  ***
      !!                   
      !! ** Purpose :   set the horizontal ocean dynamics physics
      !!
      !! ** Method  :   the eddy viscosity coef. specification depends on:
      !!              - the operator:
      !!             ln_dynldf_lap = T     laplacian operator
      !!             ln_dynldf_blp = T   bilaplacian operator
      !!              - the parameter nn_ahm_ijk_t:
      !!    nn_ahm_ijk_t  =  0 => = constant
      !!                  = 10 => = F(z) :     = constant with a reduction of 1/4 with depth 
      !!                  =-20 => = F(i,j)     = shape read in 'eddy_viscosity.nc' file
      !!                  = 20    = F(i,j)     = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_viscosity.nc'  file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!                  = 31    = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                           or |u|e^3/12 bilaplacian operator )
      !!----------------------------------------------------------------------
      INTEGER  ::   jk                ! dummy loop indices
      INTEGER  ::   ierr, inum, ios   ! local integer
      REAL(wp) ::   zah0              ! local scalar
      !
      NAMELIST/namdyn_ldf/ ln_dynldf_lap, ln_dynldf_blp,                  &
         &                 ln_dynldf_lev, ln_dynldf_hor, ln_dynldf_iso,   &
         &                 nn_ahm_ijk_t , rn_ahm_0, rn_ahm_b, rn_bhm_0
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namdyn_ldf in reference namelist : Lateral physics
      READ  ( numnam_ref, namdyn_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_ldf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namdyn_ldf in configuration namelist : Lateral physics
      READ  ( numnam_cfg, namdyn_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdyn_ldf )

      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf : set lateral mixing parameters'
         !
         WRITE(numout,*) '      type :'
         WRITE(numout,*) '         laplacian operator                   ln_dynldf_lap = ', ln_dynldf_lap
         WRITE(numout,*) '         bilaplacian operator                 ln_dynldf_blp = ', ln_dynldf_blp
         !
         WRITE(numout,*) '      direction of action :'
         WRITE(numout,*) '         iso-level                            ln_dynldf_lev = ', ln_dynldf_lev
         WRITE(numout,*) '         horizontal (geopotential)            ln_dynldf_hor = ', ln_dynldf_hor
         WRITE(numout,*) '         iso-neutral                          ln_dynldf_iso = ', ln_dynldf_iso
         !
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         type of time-space variation         nn_ahm_ijk_t  = ', nn_ahm_ijk_t
         WRITE(numout,*) '         lateral laplacian eddy viscosity     rn_ahm_0_lap  = ', rn_ahm_0, ' m2/s'
         WRITE(numout,*) '         background viscosity (iso case)      rn_ahm_b      = ', rn_ahm_b, ' m2/s'
         WRITE(numout,*) '         lateral bilaplacian eddy viscosity   rn_ahm_0_blp  = ', rn_bhm_0, ' m4/s'
      ENDIF

      !                                ! Parameter control
      IF( .NOT.ln_dynldf_lap .AND. .NOT.ln_dynldf_blp ) THEN
         IF(lwp) WRITE(numout,*) '   No viscous operator selected. ahmt and ahmf are not allocated'
         l_ldfdyn_time = .FALSE.
         RETURN
      ENDIF
      !
      IF( ln_dynldf_blp .AND. ln_dynldf_iso ) THEN     ! iso-neutral bilaplacian not implemented
         CALL ctl_stop( 'dyn_ldf_init: iso-neutral bilaplacian not coded yet' ) 
      ENDIF

      ! ... Space/Time variation of eddy coefficients
      !                                               ! allocate the ahm arrays
      ALLOCATE( ahmt(jpi,jpj,jpk) , ahmf(jpi,jpj,jpk) , STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate arrays')
      !
      ahmt(:,:,jpk) = 0._wp                           ! last level always 0  
      ahmf(:,:,jpk) = 0._wp
      !
      !                                               ! value of eddy mixing coef.
      IF    ( ln_dynldf_lap ) THEN   ;   zah0 =      rn_ahm_0         !   laplacian operator
      ELSEIF( ln_dynldf_blp ) THEN   ;   zah0 = ABS( rn_bhm_0 )       ! bilaplacian operator
      ELSE                                                                  ! NO viscous  operator
         CALL ctl_warn( 'ldf_dyn_init: No lateral viscous operator used ' )
      ENDIF
      !
      l_ldfdyn_time = .FALSE.                         ! no time variation except in case defined below
      !
      IF( ln_dynldf_lap .OR. ln_dynldf_blp ) THEN     ! only if a lateral diffusion operator is used
         !
         SELECT CASE(  nn_ahm_ijk_t  )                ! Specification of space time variations of ahmt, ahmf
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = constant '
            ahmt(:,:,:) = zah0 * tmask(:,:,:)
            ahmf(:,:,:) = zah0 * fmask(:,:,:)
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( depth )'
            ahmt(:,:,1) = zah0 * tmask(:,:,1)                      ! constant surface value
            ahmf(:,:,1) = zah0 * fmask(:,:,1)
            CALL ldf_c1d( 'DYN', r1_4, ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )
            !
         CASE ( -20 )      !== fixed horizontal shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F(i,j) read in eddy_viscosity.nc file'
            CALL iom_open( 'eddy_viscosity_2D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_2d', ahmt(:,:,1) )
            CALL iom_get ( inum, jpdom_data, 'ahmf_2d', ahmf(:,:,1) )
            CALL iom_close( inum )
!!gm Question : info for LAP or BLP case  to take into account the SQRT in the bilaplacian case ???
!!              do we introduce a scaling by the max value of the array, and then multiply by zah0 ????
!!              better:  check that the max is <=1  i.e. it is a shape from 0 to 1, not a coef that has physical dimension
            DO jk = 2, jpkm1
               ahmt(:,:,jk) = ahmt(:,:,1) * tmask(:,:,jk)
               ahmf(:,:,jk) = ahmf(:,:,1) * fmask(:,:,jk)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( e1, e2 ) or F( e1^3, e2^3 ) (lap. or blp. case)'
            IF( ln_dynldf_lap )   CALL ldf_c2d( 'DYN', 'LAP', zah0, ahmt, ahmf )    ! surface value proportional to scale factor
            IF( ln_dynldf_blp )   CALL ldf_c2d( 'DYN', 'BLP', zah0, ahmt, ahmf )    ! surface value proportional to scale factor^3
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F(i,j,k) read in eddy_diffusivity_3D.nc file'
            CALL iom_open( 'eddy_viscosity_3D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_3d', ahmt )
            CALL iom_get ( inum, jpdom_data, 'ahmf_3d', ahmf )
            CALL iom_close( inum )
!!gm Question : info for LAP or BLP case  to take into account the SQRT in the bilaplacian case ????
!!              do we introduce a scaling by the max value of the array, and then multiply by zah0 ????
            DO jk = 1, jpkm1
               ahmt(:,:,jk) = ahmt(:,:,jk) * tmask(:,:,jk)
               ahmf(:,:,jk) = ahmf(:,:,jk) * fmask(:,:,jk)
            END DO
            !
         CASE(  30  )       !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth )'
            IF( ln_dynldf_lap )   CALL ldf_c2d( 'DYN', 'LAP', zah0, ahmt, ahmf )    ! surface value proportional to scale factor
            IF( ln_dynldf_blp )   CALL ldf_c2d( 'DYN', 'BLP', zah0, ahmt, ahmf )    ! surface value proportional to scale factor
            !                                                    ! reduction with depth
            CALL ldf_c1d( 'DYN', r1_4, ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )
            !
         CASE(  31  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '                                proportional to the velocity : |u|e/12 or |u|e^3/12'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_dyn_init: wrong choice for nn_ahm_ijk_t, the type of space-time variation of ahm')
         END SELECT
         !
! JMM TEST : no SQRT here case 20   nn_ahm_ijk_t
!        IF( ln_dynldf_blp .AND. .NOT. l_ldfdyn_time .AND. nn_ahm_ijk_t /= 20 ) THEN       ! bilapcian and no time variation:
         IF( ln_dynldf_blp .AND. .NOT. l_ldfdyn_time ) THEN       ! bilapcian and no time variation:
            ahmt(:,:,:) = SQRT( ahmt(:,:,:) )                     ! take the square root of the coefficient
            ahmf(:,:,:) = SQRT( ahmf(:,:,:) )
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_dyn_init


   SUBROUTINE ldf_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn  ***
      !! 
      !! ** Purpose :   update at kt the momentum lateral mixing coeff. (ahmt and ahmf)
      !!
      !! ** Method  :   time varying eddy viscosity coefficients:
      !!
      !!    nn_ahm_ijk_t = 31    ahmt, ahmf = F(i,j,k,t) = F(local velocity) 
      !!                         ( |u|e /12  or  |u|e^3/12 for laplacian or bilaplacian operator )
      !!                BLP case : sqrt of the eddy coef, since bilaplacian is en re-entrant laplacian
      !!
      !! ** action  :    ahmt, ahmf   update at each time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu2pv2_ij_p1, zu2pv2_ij, zu2pv2_ij_m1, zetmax, zefmax   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_dyn')
      !
      SELECT CASE(  nn_ahm_ijk_t  )       !== Eddy vicosity coefficients ==!
      !
      CASE(  31  )       !==  time varying 3D field  ==!   = F( local velocity )
         !
         IF( ln_dynldf_lap   ) THEN          !   laplacian operator : |u| e /12 = |u/144| e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zetmax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     zefmax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zetmax * tmask(ji,jj,jk)      ! 288= 12*12 * 2
                     ahmf(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zefmax * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ELSEIF( ln_dynldf_blp ) THEN      ! bilaplacian operator : sqrt( |u| e^3 /12 ) = sqrt( |u/144| e ) * e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zetmax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     zefmax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zetmax  ) * zetmax * tmask(ji,jj,jk)
                     ahmf(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zefmax  ) * zefmax * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !
         CALL lbc_lnk( ahmt, 'T', 1. )   ;   CALL lbc_lnk( ahmf, 'F', 1. )
         !
      END SELECT
      !
      CALL iom_put( "ahmt_2d", ahmt(:,:,1) )   ! surface u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_2d", ahmf(:,:,1) )   ! surface v-eddy diffusivity coeff.
      CALL iom_put( "ahmt_3d", ahmt(:,:,:) )   ! 3D      u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_3d", ahmf(:,:,:) )   ! 3D      v-eddy diffusivity coeff.
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_dyn')
      !
   END SUBROUTINE ldf_dyn

   !!======================================================================
END MODULE ldfdyn

MODULE domain
   !!==============================================================================
   !!                       ***  MODULE domain   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================
   !! History :  OPA  !  1990-10  (C. Levy - G. Madec)  Original code
   !!                 !  1992-01  (M. Imbard) insert time step initialization
   !!                 !  1996-06  (G. Madec) generalized vertical coordinate 
   !!                 !  1997-02  (G. Madec) creation of domwri.F
   !!                 !  2001-05  (E.Durand - G. Madec) insert closed sea
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.3  !  2010-11  (G. Madec)  initialisation in C1D configuration
   !!            3.6  !  2013     ( J. Simeon, C. Calone, G. Madec, C. Ethe ) Online coarsening of outputs
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dom_init       : initialize the space and time domain
   !!   dom_nam        : read and contral domain namelists
   !!   dom_ctl        : control print for the ocean domain
   !!   dom_stiff      : diagnose maximum grid stiffness/hydrostatic consistency (s-coordinate)
   !!----------------------------------------------------------------------
   USE oce             ! ocean variables
   USE dom_oce         ! domain: ocean
   USE sbc_oce         ! surface boundary condition: ocean
   USE phycst          ! physical constants
   USE closea          ! closed seas
   USE domhgr          ! domain: set the horizontal mesh
   USE domzgr          ! domain: set the vertical mesh
   USE domstp          ! domain: set the time-step
   USE dommsk          ! domain: set the mask system
   USE domwri          ! domain: write the meshmask file
   USE domvvl          ! variable volume
   USE c1d             ! 1D vertical configuration
   USE dyncor_c1d      ! Coriolis term (c1d case)         (cor_c1d routine)
   !
   USE in_out_manager  ! I/O manager
   USE wrk_nemo        ! Memory Allocation
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_init   ! called by opa.F90

   !!-------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domain.F90 4624 2014-04-28 12:09:03Z acc $
   !! Software governed by the CeCILL licence        (NEMOGCM/NEMO_CeCILL.txt)
   !!-------------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!                    
      !! ** Purpose :   Domain initialization. Call the routines that are 
      !!              required to create the arrays which define the space 
      !!              and time domain of the ocean model.
      !!
      !! ** Method  : - dom_msk: compute the masks from the bathymetry file
      !!              - dom_hgr: compute or read the horizontal grid-point position
      !!                         and scale factors, and the coriolis factor
      !!              - dom_zgr: define the vertical coordinate and the bathymetry
      !!              - dom_stp: defined the model time step
      !!              - dom_wri: create the meshmask file if nmsh=1
      !!              - 1D configuration, move Coriolis, u and v at T-point
      !!----------------------------------------------------------------------
      INTEGER ::   jk          ! dummy loop indices
      INTEGER ::   iconf = 0   ! local integers
      REAL(wp), POINTER, DIMENSION(:,:) ::   z1_hu_0, z1_hv_0
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dom_init')
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init : domain initialization'
         WRITE(numout,*) '~~~~~~~~'
      ENDIF
      !
      !                       !==  Reference coordinate system  ==!
      !
                     CALL dom_nam               ! read namelist ( namrun, namdom )
                     CALL dom_clo               ! Closed seas and lake
                     CALL dom_hgr               ! Horizontal mesh
                     CALL dom_zgr               ! Vertical mesh and bathymetry
                     CALL dom_msk               ! Masks
      IF( ln_sco )   CALL dom_stiff             ! Maximum stiffness ratio/hydrostatic consistency
      !
      ht_0(:,:) = e3t_0(:,:,1) * tmask(:,:,1)   ! Reference ocean thickness
      hu_0(:,:) = e3u_0(:,:,1) * umask(:,:,1)
      hv_0(:,:) = e3v_0(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpk
         ht_0(:,:) = ht_0(:,:) + e3t_0(:,:,jk) * tmask(:,:,jk)
         hu_0(:,:) = hu_0(:,:) + e3u_0(:,:,jk) * umask(:,:,jk)
         hv_0(:,:) = hv_0(:,:) + e3v_0(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      !              !==  time varying part of coordinate system  ==!
      !
      IF( ln_linssh ) THEN          ! Fix in time : set to the reference one for all
         !       before        !          now          !       after         !
            gdept_b = gdept_0  ;   gdept_n = gdept_0   !        ---          ! depth of grid-points
            gdepw_b = gdepw_0  ;   gdepw_n = gdepw_0   !        ---          !
                                   gde3w_n = gde3w_0   !        ---          !
         !                                                                  
              e3t_b =   e3t_0  ;     e3t_n =   e3t_0   ;   e3t_a =  e3t_0    ! scale factors
              e3u_b =   e3u_0  ;     e3u_n =   e3u_0   ;   e3u_a =  e3u_0    !
              e3v_b =   e3v_0  ;     e3v_n =   e3v_0   ;   e3v_a =  e3v_0    !
                                     e3f_n =   e3f_0   !        ---          !
              e3w_b =   e3w_0  ;     e3w_n =   e3w_0   !        ---          !
             e3uw_b =  e3uw_0  ;    e3uw_n =  e3uw_0   !        ---          !
             e3vw_b =  e3vw_0  ;    e3vw_n =  e3vw_0   !        ---          !
         !
         CALL wrk_alloc( jpi,jpj,   z1_hu_0, z1_hv_0 )
         !
         z1_hu_0(:,:) = ssumask(:,:) / ( hu_0(:,:) + 1._wp - ssumask(:,:) )     ! _i mask due to ISF
         z1_hv_0(:,:) = ssvmask(:,:) / ( hv_0(:,:) + 1._wp - ssvmask(:,:) )
         !
         !        before       !          now          !       after         !
                                      ht_n =    ht_0   !                     ! water column thickness
               hu_b =    hu_0  ;      hu_n =    hu_0   ;    hu_a =    hu_0   ! 
               hv_b =    hv_0  ;      hv_n =    hv_0   ;    hv_a =    hv_0   !
            r1_hu_b = z1_hu_0  ;   r1_hu_n = z1_hu_0   ; r1_hu_a = z1_hu_0   ! inverse of water column thickness
            r1_hv_b = z1_hv_0  ;   r1_hv_n = z1_hv_0   ; r1_hv_a = z1_hv_0   !
         !
         CALL wrk_dealloc( jpi,jpj,   z1_hu_0, z1_hv_0 )
         !
      ELSE                         ! time varying : initialize before/now/after variables
         !
         CALL dom_vvl_init 
         !
      ENDIF
      !
      IF( lk_c1d         )   CALL cor_c1d       ! 1D configuration: Coriolis set at T-point
      !
                             CALL dom_stp       ! time step
      IF( nmsh /= 0 .AND. .NOT. ln_iscpl )                         CALL dom_wri      ! Create a domain file
      IF( nmsh /= 0 .AND.       ln_iscpl .AND. .NOT. ln_rstart )   CALL dom_wri      ! Create a domain file
      IF( .NOT.ln_rstart )   CALL dom_ctl       ! Domain control
      !
      IF( nn_timing == 1 )   CALL timing_stop('dom_init')
      !
   END SUBROUTINE dom_init


   SUBROUTINE dom_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read domaine namelists and print the variables.
      !!
      !! ** input   : - namrun namelist
      !!              - namdom namelist
      !!              - namnc4 namelist   ! "key_netcdf4" only
      !!----------------------------------------------------------------------
      USE ioipsl
      NAMELIST/namrun/ cn_ocerst_indir, cn_ocerst_outdir, nn_stocklist, ln_rst_list, cn_dirout ,     &
         &             nn_no   , cn_exp   , cn_ocerst_in, cn_ocerst_out, ln_rstart , nn_rstctl ,     &
         &             nn_it000, nn_itend , nn_date0    , nn_time0     , nn_leapy  , nn_istate ,     &
         &             nn_stock, nn_write , ln_mskland  , ln_clobber   , nn_chunksz, nn_euler  ,     &
         &             ln_cfmeta, ln_iscpl
      NAMELIST/namdom/ nn_bathy, rn_bathy , rn_e3zps_min, rn_e3zps_rat, nn_msh, rn_hmin, rn_isfhmin, &
         &             rn_atfp , rn_rdt   , nn_closea   , ln_crs      , jphgr_msh ,                  &
         &             ppglam0, ppgphi0, ppe1_deg, ppe2_deg, ppe1_m, ppe2_m,                         &
         &             ppsur, ppa0, ppa1, ppkth, ppacr, ppdzmin, pphmax, ldbletanh,                  &
         &             ppa2, ppkth2, ppacr2
#if defined key_netcdf4
      NAMELIST/namnc4/ nn_nchunks_i, nn_nchunks_j, nn_nchunks_k, ln_nc4zip
#endif
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      CHARACTER(LEN=255) :: cl_no
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namrun in reference namelist : Parameters of the run
      READ  ( numnam_ref, namrun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namrun in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namrun, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namrun )
      ! Add extension (job number to the restart dir. Differ for restart input and restart output
!{ DRAKKAR modification : NEMO reads restart files :<CN_OCERST_INDIR>.<<nn_no-1>>/<CN_OCERST_IN>-<<nn_no -1 >>_<RANK>.nc
      WRITE(cl_no,*) nn_no-1 ; cl_no = TRIM(ADJUSTL(cl_no) )
      cn_ocerst_indir=TRIM(cn_ocerst_indir)//'.'//TRIM(cl_no)
      cn_ocerst_in= TRIM(cn_ocerst_in)//'-'//TRIM(cl_no)

      !  DRAKKAR modification : NEMO writes restart files :<CN_OCERST_INDIR>.<<nn_no>>/<CN_OCERST_IN>-<<nn_no  >>_<RANK>.nc
      WRITE(cl_no,*) nn_no   ; cl_no = TRIM(ADJUSTL(cl_no) )
      cn_ocerst_outdir=TRIM(cn_ocerst_outdir)//'.'//TRIM(cl_no)
      cn_ocerst_out= TRIM(cn_ocerst_out)//'-'//TRIM(cl_no)
!}

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_nam  : domain initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namrun'
         WRITE(numout,*) '      job number                      nn_no      = ', nn_no
         WRITE(numout,*) '      experiment name for output      cn_exp     = ', cn_exp
         WRITE(numout,*) '      file prefix restart input       cn_ocerst_in= ', cn_ocerst_in
         WRITE(numout,*) '      restart input directory         cn_ocerst_indir= ', cn_ocerst_indir
         WRITE(numout,*) '      file prefix restart output      cn_ocerst_out= ', cn_ocerst_out
         WRITE(numout,*) '      restart output directory        cn_ocerst_outdir= ', cn_ocerst_outdir
         WRITE(numout,*) '      restart logical                 ln_rstart  = ', ln_rstart
         WRITE(numout,*) '      start with forward time step    nn_euler   = ', nn_euler
         WRITE(numout,*) '      control of time step            nn_rstctl  = ', nn_rstctl
         WRITE(numout,*) '      number of the first time step   nn_it000   = ', nn_it000
         WRITE(numout,*) '      number of the last time step    nn_itend   = ', nn_itend
         WRITE(numout,*) '      initial calendar date aammjj    nn_date0   = ', nn_date0
         WRITE(numout,*) '      initial time of day in hhmm     nn_time0   = ', nn_time0
         WRITE(numout,*) '      leap year calendar (0/1)        nn_leapy   = ', nn_leapy
         WRITE(numout,*) '      initial state output            nn_istate  = ', nn_istate
         IF( ln_rst_list ) THEN
            WRITE(numout,*) '      list of restart dump times      nn_stocklist   =', nn_stocklist
         ELSE
            WRITE(numout,*) '      frequency of restart file       nn_stock   = ', nn_stock
         ENDIF
         WRITE(numout,*) '      frequency of output file        nn_write   = ', nn_write
         WRITE(numout,*) '      mask land points                ln_mskland = ', ln_mskland
         WRITE(numout,*) '      additional CF standard metadata ln_cfmeta  = ', ln_cfmeta
         WRITE(numout,*) '      overwrite an existing file      ln_clobber = ', ln_clobber
         WRITE(numout,*) '      NetCDF chunksize (bytes)        nn_chunksz = ', nn_chunksz
         WRITE(numout,*) '      Output directory                cn_dirout  = ', cn_dirout
         WRITE(numout,*) '      IS coupling at the restart step ln_iscpl   = ', ln_iscpl
      ENDIF

      no = nn_no                    ! conversion DOCTOR names into model names (this should disappear soon)
      cexper = cn_exp
      nrstdt = nn_rstctl
      nit000 = nn_it000
      nitend = nn_itend
      ndate0 = nn_date0
      nleapy = nn_leapy
      ninist = nn_istate
      nstock = nn_stock
      nstocklist = nn_stocklist
      nwrite = nn_write
      neuler = nn_euler
      IF ( neuler == 1 .AND. .NOT. ln_rstart ) THEN
         WRITE(ctmp1,*) 'ln_rstart =.FALSE., nn_euler is forced to 0 '
         CALL ctl_warn( ctmp1 )
         neuler = 0
      ENDIF

      !                             ! control of output frequency
      IF ( nstock == 0 .OR. nstock > nitend ) THEN
         WRITE(ctmp1,*) 'nstock = ', nstock, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nstock = nitend
      ENDIF
      IF ( nwrite == 0 ) THEN
         WRITE(ctmp1,*) 'nwrite = ', nwrite, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nwrite = nitend
      ENDIF

      IF( Agrif_Root() ) THEN
      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE (  1 ) 
         CALL ioconf_calendar('gregorian')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "gregorian", i.e. leap year'
      CASE (  0 )
         CALL ioconf_calendar('noleap')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "noleap", i.e. no leap year'
      CASE ( 30 )
         CALL ioconf_calendar('360d')
         IF(lwp) WRITE(numout,*) '   The IOIPSL calendar is "360d", i.e. 360 days in a year'
      END SELECT
      ENDIF

      REWIND( numnam_ref )              ! Namelist namdom in reference namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_ref, namdom, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in reference namelist', lwp )
  
      !
      REWIND( numnam_cfg )              ! Namelist namdom in configuration namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_cfg, namdom, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdom in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdom )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namdom : space & time domain'
         WRITE(numout,*) '      flag read/compute bathymetry      nn_bathy     = ', nn_bathy
         WRITE(numout,*) '      Depth (if =0 bathy=jpkm1)         rn_bathy     = ', rn_bathy
         WRITE(numout,*) '      min depth of the ocean    (>0) or    rn_hmin   = ', rn_hmin
         WRITE(numout,*) '      min number of ocean level (<0)       '
         WRITE(numout,*) '      treshold to open the isf cavity   rn_isfhmin   = ', rn_isfhmin, ' (m)'
         WRITE(numout,*) '      minimum thickness of partial      rn_e3zps_min = ', rn_e3zps_min, ' (m)'
         WRITE(numout,*) '         step level                     rn_e3zps_rat = ', rn_e3zps_rat
         WRITE(numout,*) '      create mesh/mask file(s)          nn_msh       = ', nn_msh
         WRITE(numout,*) '           = 0   no file created           '
         WRITE(numout,*) '           = 1   mesh_mask                 '
         WRITE(numout,*) '           = 2   mesh and mask             '
         WRITE(numout,*) '           = 3   mesh_hgr, msh_zgr and mask'
         WRITE(numout,*) '      ocean time step                       rn_rdt    = ', rn_rdt
         WRITE(numout,*) '      asselin time filter parameter         rn_atfp   = ', rn_atfp
         WRITE(numout,*) '      suppression of closed seas (=0)       nn_closea = ', nn_closea
         WRITE(numout,*) '      online coarsening of dynamical fields ln_crs    = ', ln_crs
         WRITE(numout,*) '      type of horizontal mesh jphgr_msh           = ', jphgr_msh
         WRITE(numout,*) '      longitude of first raw and column T-point ppglam0 = ', ppglam0
         WRITE(numout,*) '      latitude  of first raw and column T-point ppgphi0 = ', ppgphi0
         WRITE(numout,*) '      zonal      grid-spacing (degrees) ppe1_deg        = ', ppe1_deg
         WRITE(numout,*) '      meridional grid-spacing (degrees) ppe2_deg        = ', ppe2_deg
         WRITE(numout,*) '      zonal      grid-spacing (degrees) ppe1_m          = ', ppe1_m
         WRITE(numout,*) '      meridional grid-spacing (degrees) ppe2_m          = ', ppe2_m
         WRITE(numout,*) '      ORCA r4, r2 and r05 coefficients  ppsur           = ', ppsur
         WRITE(numout,*) '                                        ppa0            = ', ppa0
         WRITE(numout,*) '                                        ppa1            = ', ppa1
         WRITE(numout,*) '                                        ppkth           = ', ppkth
         WRITE(numout,*) '                                        ppacr           = ', ppacr
         WRITE(numout,*) '      Minimum vertical spacing ppdzmin                  = ', ppdzmin
         WRITE(numout,*) '      Maximum depth pphmax                              = ', pphmax
         WRITE(numout,*) '      Use double tanf function for vertical coordinates ldbletanh = ', ldbletanh
         WRITE(numout,*) '      Double tanh function parameters ppa2              = ', ppa2
         WRITE(numout,*) '                                      ppkth2            = ', ppkth2
         WRITE(numout,*) '                                      ppacr2            = ', ppacr2
      ENDIF
      !
      ntopo     = nn_bathy          ! conversion DOCTOR names into model names (this should disappear soon)
      e3zps_min = rn_e3zps_min
      e3zps_rat = rn_e3zps_rat
      nmsh      = nn_msh
      atfp      = rn_atfp
      rdt       = rn_rdt

#if defined key_netcdf4
      !                             ! NetCDF 4 case   ("key_netcdf4" defined)
      REWIND( numnam_ref )              ! Namelist namnc4 in reference namelist : NETCDF
      READ  ( numnam_ref, namnc4, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namnc4 in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namnc4 in configuration namelist : NETCDF
      READ  ( numnam_cfg, namnc4, IOSTAT = ios, ERR = 908 )
908   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namnc4 in configuration namelist', lwp )
      IF(lwm) WRITE( numond, namnc4 )

      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namnc4 - Netcdf4 chunking parameters'
         WRITE(numout,*) '      number of chunks in i-dimension      nn_nchunks_i   = ', nn_nchunks_i
         WRITE(numout,*) '      number of chunks in j-dimension      nn_nchunks_j   = ', nn_nchunks_j
         WRITE(numout,*) '      number of chunks in k-dimension      nn_nchunks_k   = ', nn_nchunks_k
         WRITE(numout,*) '      apply netcdf4/hdf5 chunking & compression ln_nc4zip = ', ln_nc4zip
      ENDIF

      ! Put the netcdf4 settings into a simple structure (snc4set, defined in in_out_manager module)
      ! Note the chunk size in the unlimited (time) dimension will be fixed at 1
      snc4set%ni   = nn_nchunks_i
      snc4set%nj   = nn_nchunks_j
      snc4set%nk   = nn_nchunks_k
      snc4set%luse = ln_nc4zip
#else
      snc4set%luse = .FALSE.        ! No NetCDF 4 case
#endif
      !
   END SUBROUTINE dom_nam


   SUBROUTINE dom_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_ctl  ***
      !!
      !! ** Purpose :   Domain control.
      !!
      !! ** Method  :   compute and print extrema of masked scale factors
      !!----------------------------------------------------------------------
      INTEGER ::   iimi1, ijmi1, iimi2, ijmi2, iima1, ijma1, iima2, ijma2
      INTEGER, DIMENSION(2) ::   iloc   ! 
      REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
      !!----------------------------------------------------------------------
      !
      IF(lk_mpp) THEN
         CALL mpp_minloc( e1t(:,:), tmask_i(:,:), ze1min, iimi1,ijmi1 )
         CALL mpp_minloc( e2t(:,:), tmask_i(:,:), ze2min, iimi2,ijmi2 )
         CALL mpp_maxloc( e1t(:,:), tmask_i(:,:), ze1max, iima1,ijma1 )
         CALL mpp_maxloc( e2t(:,:), tmask_i(:,:), ze2max, iima2,ijma2 )
      ELSE
         ze1min = MINVAL( e1t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze2min = MINVAL( e2t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze1max = MAXVAL( e1t(:,:), mask = tmask_i(:,:) == 1._wp )    
         ze2max = MAXVAL( e2t(:,:), mask = tmask_i(:,:) == 1._wp )    

         iloc  = MINLOC( e1t(:,:), mask = tmask_i(:,:) == 1._wp )
         iimi1 = iloc(1) + nimpp - 1
         ijmi1 = iloc(2) + njmpp - 1
         iloc  = MINLOC( e2t(:,:), mask = tmask_i(:,:) == 1._wp )
         iimi2 = iloc(1) + nimpp - 1
         ijmi2 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e1t(:,:), mask = tmask_i(:,:) == 1._wp )
         iima1 = iloc(1) + nimpp - 1
         ijma1 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e2t(:,:), mask = tmask_i(:,:) == 1._wp )
         iima2 = iloc(1) + nimpp - 1
         ijma2 = iloc(2) + njmpp - 1
      ENDIF
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, iima1, ijma1
         WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, iimi1, ijmi1
         WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, iima2, ijma2
         WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, iimi2, ijmi2
      ENDIF
      !
   END SUBROUTINE dom_ctl


   SUBROUTINE dom_stiff
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_stiff  ***
      !!                     
      !! ** Purpose :   Diagnose maximum grid stiffness/hydrostatic consistency
      !!
      !! ** Method  :   Compute Haney (1991) hydrostatic condition ratio
      !!                Save the maximum in the vertical direction
      !!                (this number is only relevant in s-coordinates)
      !!
      !!                Haney, R. L., 1991: On the pressure gradient force
      !!                over steep topography in sigma coordinate ocean models. 
      !!                J. Phys. Oceanogr., 21, 610???619.
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk 
      REAL(wp) ::   zrxmax
      REAL(wp), DIMENSION(4) ::   zr1
      !!----------------------------------------------------------------------
      rx1(:,:) = 0._wp
      zrxmax   = 0._wp
      zr1(:)   = 0._wp
      !
      DO ji = 2, jpim1
         DO jj = 2, jpjm1
            DO jk = 1, jpkm1
               zr1(1) = ABS(  ( gdepw_0(ji  ,jj,jk  )-gdepw_0(ji-1,jj,jk  )               & 
                    &          +gdepw_0(ji  ,jj,jk+1)-gdepw_0(ji-1,jj,jk+1) )             &
                    &       / ( gdepw_0(ji  ,jj,jk  )+gdepw_0(ji-1,jj,jk  )               &
                    &          -gdepw_0(ji  ,jj,jk+1)-gdepw_0(ji-1,jj,jk+1) + rsmall )  ) * umask(ji-1,jj,jk)
               zr1(2) = ABS(  ( gdepw_0(ji+1,jj,jk  )-gdepw_0(ji  ,jj,jk  )               &
                    &          +gdepw_0(ji+1,jj,jk+1)-gdepw_0(ji  ,jj,jk+1) )             &
                    &       / ( gdepw_0(ji+1,jj,jk  )+gdepw_0(ji  ,jj,jk  )               &
                    &          -gdepw_0(ji+1,jj,jk+1)-gdepw_0(ji  ,jj,jk+1) + rsmall )  ) * umask(ji  ,jj,jk)
               zr1(3) = ABS(  ( gdepw_0(ji,jj+1,jk  )-gdepw_0(ji,jj  ,jk  )               &
                    &          +gdepw_0(ji,jj+1,jk+1)-gdepw_0(ji,jj  ,jk+1) )             &
                    &       / ( gdepw_0(ji,jj+1,jk  )+gdepw_0(ji,jj  ,jk  )               &
                    &          -gdepw_0(ji,jj+1,jk+1)-gdepw_0(ji,jj  ,jk+1) + rsmall )  ) * vmask(ji,jj  ,jk)
               zr1(4) = ABS(  ( gdepw_0(ji,jj  ,jk  )-gdepw_0(ji,jj-1,jk  )               &
                    &          +gdepw_0(ji,jj  ,jk+1)-gdepw_0(ji,jj-1,jk+1) )             &
                    &       / ( gdepw_0(ji,jj  ,jk  )+gdepw_0(ji,jj-1,jk  )               &
                    &          -gdepw_0(ji,jj  ,jk+1)-gdepw_0(ji,jj-1,jk+1) + rsmall )  ) * vmask(ji,jj-1,jk)
               zrxmax = MAXVAL( zr1(1:4) )
               rx1(ji,jj) = MAX( rx1(ji,jj) , zrxmax )
            END DO
         END DO
      END DO
      CALL lbc_lnk( rx1, 'T', 1. )
      !
      zrxmax = MAXVAL( rx1 )
      !
      IF( lk_mpp )   CALL mpp_max( zrxmax ) ! max over the global domain
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_stiff : maximum grid stiffness ratio: ', zrxmax
         WRITE(numout,*) '~~~~~~~~~'
      ENDIF
      !
   END SUBROUTINE dom_stiff

   !!======================================================================
END MODULE domain

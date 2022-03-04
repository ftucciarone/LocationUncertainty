!!===========================================================================
!!                       ***  MODULE  tlusvd  ***
!! Transport under Location Uncertainty: stochastic noise SVD
!!                       noise generation based on pseudo-observations and SVD.
!!===========================================================================
!! History : 0.0  !  2017-..  (P. DERIAN)  Original code
!!
!! [TODO]    - read namelist from _cfg, not only _ref
!!           - write error messages to appropriate unit
!!           - use NEMO's work arrays?
!!---------------------------------------------------------------------------
MODULE tlusvd
   ! [mod_dep]
   USE par_kind           ! data types defined in par_kind module
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce            ! ocean dynamics parameters
   USE oce                ! ocean dynamics and active tracers
   USE dom_oce            ! ocean space and time domain
   USE lib_mpp            ! MPP library
   USE timing             ! Timing
   USE lbclnk             ! ocean lateral boundary conditions (or mpp link)
   USE tlu
   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_svd_init    ! called by tlu_init in tlu.f90
   PUBLIC tlu_noi         ! called by step.f90
   PUBLIC read_vel_mod
   PUBLIC read_wnd_mod
   !
   INCLUDE 'mpif.h'

   ! [tlu_switches] 
   LOGICAL,     PUBLIC                                      :: ln_tlu_svd = .FALSE.    !> @public   SVD-based noise
   LOGICAL,     PUBLIC                                      :: ln_tlu_bia = .FALSE.    !> @public   Mean-component bias
   ! [tlu_switches] 
   !
   ! [tlu_nmodes]
   INTEGER(i4), PUBLIC                                      :: nm                      !> @public   Number of modes
   ! [tlu_nmodes]
   !
   ! [tlu_spmodes] 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spx_mod                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spy_mod                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spz_mod                 !> @public   W-velocity modes 
   ! [tlu_spmodes]
   !
   ! [tlu_windmodes]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wsx_mod                 !> @public   Wind U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wsy_mod                 !> @public   Wind V-velocity modes
   ! [tlu_windmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:)     :: brwn_rv                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_ken_modes]
   INTEGER(i4), PUBLIC                                      :: nn_tlu_nmod             !> @public   Number of ocean POD modes (static) 
   INTEGER(i4), PUBLIC                                      :: nn_tlu_wmod             !> @public   Number of wind POD modes
   ! [tlu_ken_modes]
   !
   ! [tlu_kenspect]
   CHARACTER(lc), PUBLIC                                    :: nm_mod_nc               !> @public   Name of the NetCDF file containing p. comp
   REAL(wp), PUBLIC                                         :: tlu_lam_spc             !> @public   Percentage of the energy spectrum
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: eigens                  !> @public   Principal components
   INTEGER(i4), PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: nmz                     !> @public   Number of modes in each layer
   INTEGER(i4), PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: midz                    !> @public   Index of first mode in new layer
   ! [tlu_kenspect]

   !
   ! Things to be modified, understood or even taken away
   !
   REAL(wp),    PUBLIC                                      :: tkelr_tlu
   REAL(wp),    PUBLIC                                      :: tkehr_tlu
   REAL(wp),    PUBLIC                                      :: ls_tlu
   REAL(wp),    PUBLIC                                      :: cof_tke      ! Coefficients for noise/variance and cof_tke is coefficient for tke scaling according to TKE in tkelr_tlu or tkehr_tlu
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     :: std_dev      ! time for noise based on best time-match obtained by mode match
   REAL(wp),    PUBLIC                                      :: sv_spm_mat
   INTEGER(i4), PUBLIC                                      :: nn_tco_dat   !Number of the first time-coeff available from reference POD


CONTAINS

   SUBROUTINE tlu_svd_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_svd_init  ***
      !!
      !! ** Purpose :   Initialize module variables.
      !!
      !! ** Method  :   Read namelist, compute other parameters and allocate arrays.
      !!
      !! ** Action  :   Arrays allocated.
      !!
      !!------------------------------------------------------------------------
      INTEGER(i4) :: ios, ierr(4), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Read namelist: transport under location uncertainty parametrization
      !
      NAMELIST/namtlu_noi/   ln_tlu_svd,   ln_tlu_bia,  nm_mod_nc,   tlu_lam_spc,    & 
                         &   nn_tlu_wmod,  nn_tlu_nmod, tkelr_tlu,      &
                         &   tkehr_tlu,    ls_tlu,       nn_tco_dat, sv_spm_mat, biaSIGN
      !
      ! Read namelist
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu_noi, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu_noi, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_svd_init : TLU, data-based noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '      data-based noise model           ln_tlu_svd  = ', ln_tlu_svd
      END IF
      !
      ! Check multiple method error:
      !
      chk = 0
      IF ( ln_tlu_svd ) chk = chk+1

      IF ( chk .lt. 1) THEN
            WRITE(numout,*) ' E R R O R! No model chosen for noise formulation! ', ierr
            STOP
      ELSEIF ( chk .gt. 1) THEN
            WRITE(numout,*) ' E R R O R! Multiple models chosen for noise formulation! ', ierr
            STOP
      ENDIF
      !
      ! SVD-Based noise
      !
      IF ( ln_tlu_svd ) THEN
         !
         ! Allocate index vectors
         !
         ALLOCATE( nmz(jpk-1), midz(jpk) )
         !
         ! Read (if present) eigenvalues data layer by layer (nm_mod_nc in namelist_cfg)
         !
         INQUIRE(FILE=nm_mod_nc, EXIST=file_exists)

         IF (file_exists) THEN
            !
            ! Compute the number of modes layer by layer (nmz) and global (nm)
            !
            CALL tlu_nmz(  )
            nm = SUM(nmz)
            !
            ! Allocate Spatial modes
            !
            ierr = 0
            ALLOCATE( spx_mod(jpi,jpj,nm), &
                   &  spy_mod(jpi,jpj,nm), &
                   &  spz_mod(jpi,jpj,nm), stat=ierr(1) )   ! the singular values
            !
            ! Read Spatial modes
            !
            CALL read_vel_mod( )
            !
            !
         ELSE
            !
            !  Compute the number of modes with a fixed number for all layers
            !
            IF (lwp) THEN
               WRITE(numout,*) '      number of spatial modes to read   nn_tlu_nmod = ', nn_tlu_nmod
            END IF
            !
            ! Compute index vector in fixed number case
            !
            DO ji = 1, jpkm1
               midz(ji) = (ji-1) * nn_tlu_nmod +1
            END DO
            nmz = nn_tlu_nmod
            nm = SUM(nmz)
            !
            ! Allocate Spatial modes
            !
            ierr = 0
            ALLOCATE( spx_mod(jpi,jpj,nm), &
                   &  spy_mod(jpi,jpj,nm), &
                   &  spz_mod(jpi,jpj,nm), stat=ierr(1) )
            !
            ! Read Spatial modes
            !
            CALL read_vel_mod( )
            !
            !
         ENDIF

         IF (ln_wnd) THEN
            !
            ! Allocate Wind modes
            !
            ierr = 0
            ALLOCATE( wsx_mod(jpi,jpj,nn_tlu_wmod), &
                   &  wsy_mod(jpi,jpj,nn_tlu_wmod), stat=ierr(1) )
            !
            ! Read Wind modes
            ! 
            CALL read_wnd_mod( )
            !
            !
         END IF
         !
         ! Allocate Brownian motion
         !
         ALLOCATE( brwn_rv(nm), std_dev(nm),  stat=ierr(2) ) 

         IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
            WRITE(numout,*) ' tlu_svd_init(): allocation failed = ', ierr   ! [TODO] should be ctmp1? instead of numout
            !CALL ctl_stop( ctmp1 )   ! [TODO] enable
            STOP
         END IF
         !
         !
      END IF

   END SUBROUTINE tlu_svd_init


   SUBROUTINE read_vel_mod( )
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE read_spa_mod  ***
      !!
      !! ** Purpose :   Read from file the spatial bases for forming noise
      !!
      !! ** Method :    Read from pre-computed file the bases.
      !!
      !! ** Action :    spx/y/z_mod contains the three component spatially organised
      !!                modes with the total number specified from the namelist
      !!
      !! ** Note:       Reading the variable using fortran read may not work!
      !!                For mode matching, make sure to normalize spatial basis
      !!                with singular value
      !!------------------------------------------------------------------------
      !
      INTEGER ::   jk, jp, idx, ierr ! dummy loop arguments
      INTEGER ::   endloop

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwa1, zwa2, zwa3
      INTEGER                                 :: ucid, vcid, wcid
      CHARACTER (len=90)                      :: ubas, vbas, wbas

      REAL(wp)  :: val

      !
      ! Allocate dummy matrices
      !
      ALLOCATE( zwa1(jpi,jpj,jpk), &
             &  zwa2(jpi,jpj,jpk), &
             &  zwa3(jpi,jpj,jpk),  stat=ierr ) 
      IF ( ierr /= 0 ) THEN  
         WRITE(numout,*) ' read_spa_mod2(): allocation failed = ', ierr   
         STOP
      END IF  
      !
      ! Get netCDF file ID
      !
      CALL iom_open('spmu.nc', ucid )
      CALL iom_open('spmv.nc', vcid )
      CALL iom_open('spmw.nc', wcid )

      val = 1._wp
      endloop = 0 


      DO jp = 1, MAXVAL(nmz)
         !
         ! Initialize dummy matrices
         !
         zwa1 = 0._wp
         zwa2 = 0._wp
         zwa3 = 0._wp
         !
         ! Write variable name 
         !
         WRITE (ubas, '(a, i3.3)' )  "spat_basis_u_", jp
         WRITE (vbas, '(a, i3.3)' )  "spat_basis_v_", jp
         WRITE (wbas, '(a, i3.3)' )  "spat_basis_w_", jp
         !
         ! Read modes (in X, Y, Z)
         !
         CALL iom_get(ucid, jpdom_autoglo, ubas, zwa1, ldxios = lrxios)
         CALL iom_get(vcid, jpdom_autoglo, vbas, zwa2, ldxios = lrxios)
         CALL iom_get(wcid, jpdom_autoglo, wbas, zwa3, ldxios = lrxios)

         !                                ! ================
         DO jk = 1, jpkm1                 ! Horizontal slab
            !                             ! ================
            IF ( jp <= nmz(jk) ) THEN
               idx = midz(jk) + jp -1
               spx_mod(:,:,idx) = zwa1(:,:,jk) * val
               spy_mod(:,:,idx) = zwa2(:,:,jk) * val
               spz_mod(:,:,idx) = zwa3(:,:,jk) * val
            END IF
            !                             ! ================
         END DO                           !   End of slab
         !                                ! ================
!         IF ( jp <= nmz(1) ) THEN
!            endloop = endloop + 1
!            CALL cml_int2d_ene3c(  zwa1(:,:,1),  zwa2(:,:,1),  zwa3(:,:,1), 1._wp, val)
!         END IF
      ENDDO
!      print '(A20, I4, A12, E16.7)', ' Sum of first ',endloop, ' modes (int)', val
      !
      ! Read bias (if flagged)
      !
      IF (ln_tlu_bia) CALL iom_get(ucid,jpdom_autoglo,'u_mean',ubia, ldxios = lrxios)
      IF (ln_tlu_bia) CALL iom_get(vcid,jpdom_autoglo,'v_mean',vbia, ldxios = lrxios)
      IF (ln_tlu_bia) CALL iom_get(wcid,jpdom_autoglo,'w_mean',wbia, ldxios = lrxios)
      !
      ! Close netCDF files
      !
      CALL iom_close(ucid)
      CALL iom_close(vcid)
      CALL iom_close(wcid)
      !
      ! Deallocate dummy matrices
      !
      DEALLOCATE(zwa1, zwa2, zwa3)
   END SUBROUTINE read_vel_mod

   SUBROUTINE read_wnd_mod( )
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE read_spa_mod  ***
      !!
      !! ** Purpose :   Read from file the spatial bases for forming noise
      !!
      !! ** Method :    Read from pre-computed file the bases.
      !!
      !! ** Action :    spx/y/z_mod contains the three component spatially organised
      !!                modes with the total number specified from the namelist
      !!
      !! ** Note:       Reading the variable using fortran read may not work!
      !!                For mode matching, make sure to normalize spatial basis
      !!                with singular value
      !!------------------------------------------------------------------------
      !
      INTEGER ::   jp, ierr ! dummy loop arguments


      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zwa1, zwa2
      INTEGER                               :: ucid, vcid
      CHARACTER (len=90)                    :: ubas, vbas

      !
      ! Allocate dummy matrices
      !
      ALLOCATE( zwa1(jpi,jpj), &
             &  zwa2(jpi,jpj), stat=ierr ) 
      IF ( ierr /= 0 ) THEN  
         WRITE(numout,*) ' read_wnd_mod(): allocation failed = ', ierr   
         STOP
      END IF  
      !
      ! Get netCDF file ID
      !
      CALL iom_open('spmu.nc', ucid )
      CALL iom_open('spmv.nc', vcid )

      DO jp = 1, nn_tlu_wmod
         !
         ! Initialize dummy matrices
         !
         zwa1 = 0._wp
         zwa2 = 0._wp
         !
         ! Write variable name 
         !
         WRITE (ubas, '(a, i3.3)' )  "spat_basis_uwnd_", jp
         WRITE (vbas, '(a, i3.3)' )  "spat_basis_vwnd_", jp
         !
         ! Read modes (in X, Y, Z)
         !
         CALL iom_get(ucid, jpdom_autoglo, ubas, zwa1, ldxios = lrxios)
         CALL iom_get(vcid, jpdom_autoglo, vbas, zwa2, ldxios = lrxios)

         wsx_mod(:,:,jp) = zwa1(:,:)
         wsy_mod(:,:,jp) = zwa2(:,:)
      
      ENDDO
      !
      ! Read bias (if flagged)
      !
      ! IF (ln_tlu_bia) CALL iom_get(ucid,jpdom_autoglo,'uwnd_mean',ubia, ldxios = lrxios)
      ! IF (ln_tlu_bia) CALL iom_get(vcid,jpdom_autoglo,'vwnd_mean',vbia, ldxios = lrxios)
      !
      ! Close netCDF files
      !
      CALL iom_close(ucid)
      CALL iom_close(vcid)
      !
      ! Deallocate dummy matrices
      !
      DEALLOCATE(zwa1, zwa2)

   END SUBROUTINE read_wnd_mod

   SUBROUTINE tlu_noi( kt )
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_noi  ***
      !!
      !! ** Purpose :   Construct one realisation of noise from spatial basis
      !!
      !! ** Method :    For each time-step draw from the spatial bases (kmod in number)
      !!                a realisation of the noise (sigma dBt) taking in to account
      !!                energy of small scales with (kene) and the ratio of length
      !!                scales (ls_tlu) and the time-step (rdt)
      !!
      !! ** Action :    unoi, vnoi, wnoi      : the three component noise ;
      !!                spx/y/z_mod      : the orthogonal spatial modes;
      !!                var_ten : the six component variance tensor (a)
      !!                order of a : uu, vv, ww, uv, uw, vw
      !!
      !! ** Note:       Var_ten is calculated at the beginning only as it is considered
      !!                stationary. For non-stationary, to be modified.
      !! ** TODO :      Make sure var_ten is defined on the T grid!
      !!
      !!------------------------------------------------------------------------
      !
      INTEGER, INTENT(in   ) ::   kt         ! ocean time-step index
      !
      INTEGER                                 :: ierr
      !!------------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_draw_real') ![NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_svd : SVD-based stochastic noise'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      
      !
      ! Draw k random (standard) gaussian distribution
      !
      CALL r8vec_normal_01( nm, brwn_rv )  
      CALL MPI_BCAST( brwn_rv, nm, MPI_INTEGER,0,MPI_COMM_WORLD, ierr )
      ! 
      ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
      !
      CALL tlu_noi_pod( brwn_rv )
      ! 
      ! Compute the variance tensor a = sum_{k} \phi_{k}\phi_{k} * dt
      !
      IF( kt == nit000 )  THEN                      ! Stationary assumption here
         std_dev(:) = 1._wp                         ! Equal weighting for classic POD based noise
         CALL tlu_var_pod( rn_rdt, std_dev )   ! Computation of variance tensor weighted by standard deviation
         CALL tlu_isd
      ENDIF
      ! 
      ! Rescale the noise by some manipulation: All this part should go away at some point!
      !               dBt = C * (      sum_{k} \phi_{k} \xi_{k} )
      !                 a = C * ( dt * sum_{k} \phi_{k}\phi_{k} )
      !
      call tlu_norm_tke ( kt, tkelr_tlu, 4, .TRUE.)  !Normalise as a random noise with energy conservation
!     call tlu_norm_tke ( kt, tkelr_tlu, 2, .FALSE., tkehr_tlu)  !Normalise as a random noise with higher dissipation

      IF( kt == nit000 )  THEN !Print variance tensor (a) at time t0 i.e. stationary
         CALL iom_put( 'var_uu_once', var_ten(:,:,:,1) )  
         CALL iom_put( 'var_vv_once', var_ten(:,:,:,2) )  
         CALL iom_put( 'var_ww_once', var_ten(:,:,:,3) )  
         CALL iom_put( 'var_uv_once', var_ten(:,:,:,4) )  
         CALL iom_put( 'var_uw_once', var_ten(:,:,:,5) )  
         CALL iom_put( 'var_vw_once', var_ten(:,:,:,6) ) 
         !
         ! Print Ito-Stokes Drift
         !
         CALL iom_put( 'U_isd', uisd_n )  
         CALL iom_put( 'V_isd', visd_n )  
         CALL iom_put( 'W_isd', wisd_n )   
      END IF
   
      !Lateral boundary condition transfer across nodes
      CALL lbc_lnk_multi( 'tlu_noi_crt', unoi , 'U', -1., vnoi , 'V', -1., wnoi, 'W', -1.)

      CALL iom_put( "unoi", unoi )        ! Output noise fields
      CALL iom_put( "vnoi", vnoi )        ! Output noise fields
      CALL iom_put( "wnoi", wnoi )        ! Output noise fields
      !
      IF( ln_timing )  CALL timing_stop('tlu_draw_real')   ! [NEMO] check
      !
   END SUBROUTINE tlu_noi

   SUBROUTINE tlu_noi_pod( tcoeff )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE Create one realisation of noise  ***
      !!
      !! ** Purpose :   Formulate noise from spatial modes and time coeff/random numbers
      !!
      !! ** Method  :   Project spatial modes on corresponding temporal coefficient or
      !!                on gaussian random number as provided by input.
      !!                The noise is divided by the time step because later on in the 
      !!                NEMO deterministic kernel it will multiplied by rn_dt
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_tlu_nmod), INTENT( in ) :: tcoeff
      INTEGER                                        :: ji, jj, jk, jm
      INTEGER                                        :: stt, stp, idt

      !
      ! Initialize sigma dBt
      !
      unoi = 0._wp
      vnoi = 0._wp
      wnoi = 0._wp
      DO jk = 1, jpkm1
         !
         ! Define start and stop indexes in the compact storage form of spX_mod
         !
         stt = midz(jk  )
         stp = midz(jk  ) + nmz(jk)   ! It's not midz(jk+1) for consistency with variance
         DO jm= stt, stp
            !
            ! Define temporal index (i.e. rescale to mod-1)
            !
            idt = jm - stt + 1
            DO jj = 1, jpjm1
               DO ji= 1, jpim1
                  !
                  ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
                  !
                  unoi(ji,jj,jk) = unoi(ji,jj,jk) + spx_mod(ji,jj,jm) * tcoeff(idt)
                  vnoi(ji,jj,jk) = vnoi(ji,jj,jk) + spy_mod(ji,jj,jm) * tcoeff(idt)
                  wnoi(ji,jj,jk) = wnoi(ji,jj,jk) + spz_mod(ji,jj,jm) * tcoeff(idt)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
   END SUBROUTINE tlu_noi_pod  

   SUBROUTINE tlu_var_pod( dt, var_tco )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE Create one realisation of noise  ***
      !!
      !! ** Purpose :   Formulate variance tensor from spatial modes and variance of time_coeff
      !!
      !! ** Method  :   Average spatial modes onto corresponding mesh weighted by the variance
      !!                calculated from each corresponding time-coefficient value (MM method)
      !!                or a constant (POD method)
      !!------------------------------------------------------------------------
      REAL(wp),                         INTENT( in ) :: dt        ! Time interval
      REAL(wp), DIMENSION(nn_tlu_nmod), INTENT( in ) :: var_tco
      INTEGER :: ji, jj, jk, jm
      INTEGER :: stt, stp, idt, dkp1, dkm1
      !
      ! Initialize variance tensor
      !
      var_ten(:,:,:,:) = 0._wp
      !
      ! Boundary condiditions
      !
      !                                ! =====================
      ! jk = 1                         ! First Horizontal slab
         !                             ! =====================
         dkp1 = nmz(1)
         dkm1 = 0
         !
         stt = 1
         stp = nmz(1)
         !
         DO jm = stt, stp
            idt = jm - stt + 1
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  ! 
                  ! a_{xx}(x_{i},y,z) = 0.5 * \left[  a_{xx}(x_{i+1/2},y,z) +  a_{xx}(x_{i-1/2},y,z) \right]
                  !
                  !   with  a_{xx}(x,y,z) = \sum_{n=1)^{N} \phi^{(x)}_{n}(x,y,z) \phi^{(x)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, 1, 1) =  var_ten(ji  ,jj  , 1,1) +    0.5_wp * var_tco(idt)       *     &
                         &              ( (spx_mod(ji  ,jj  ,jm  ) * spx_mod(ji  ,jj  ,jm  ))       +     &
                         &                (spx_mod(ji-1,jj  ,jm  ) * spx_mod(ji-1,jj  ,jm  )) )
                  ! 
                  ! a_{yy}(x,y_{i},z) = 0.5 * \left[  a_{yy}(x,y_{j+1/2},z) +  a_{yy}(x,y_{j-1/2},z) \right]
                  !
                  !   with  a_{yy}(x,y,z) = \sum_{n=1)^{N} \phi^{(y)}_{n}(x,y,z) \phi^{(y)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, 1, 2) =  var_ten(ji  ,jj  , 1,2) +    0.5_wp * var_tco(idt)       *     &
                         &              ( (spy_mod(ji  ,jj  ,jm  ) * spy_mod(ji  ,jj  ,jm  ))       +     &
                         &                (spy_mod(ji  ,jj-1,jm  ) * spy_mod(ji  ,jj-1,jm  )) )
                  ! 
                  ! a_{zz}(x,y,z_{k-1/2}) = \sum_{n=1)^{N(k-1/2)} \phi^{(z)}_{n}(x,y,z) \phi^{(z)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, 1, 3) =  var_ten(ji  ,jj  , 1,3) +    0.5_wp * var_tco(idt)       *     &
                         &              ( (spz_mod(ji  ,jj  ,jm + dkp1) * spz_mod(ji  ,jj  ,jm + dkp1))   +     &
                         &                (spz_mod(ji  ,jj  ,jm       ) * spz_mod(ji  ,jj  ,jm       )) )
   
         
                  var_ten(ji, jj, 1, 4) =  var_ten(ji  ,jj  , 1,4) +    0.25_wp * var_tco(idt)       *     &
                         &              ( (spx_mod(ji  ,jj+1,jm  ) + spx_mod(ji  ,jj  ,jm  ))        *     &
                         &                (spy_mod(ji+1,jj  ,jm  ) + spy_mod(ji  ,jj  ,jm  )) )      ! uv on f grid        

                  ! Components in z that requires backwards interpolation: this is the standard part (set wind contribution
                  ! to zero in interpolation procedure)
                  var_ten(ji, jj, 1, 5) =  var_ten(ji  ,jj  , 1,5) +   0.25_wp * var_tco(idt)        *     &
                         &              ( (spx_mod(ji  ,jj  ,jm  )                            )      *     &
                         &                (spz_mod(ji+1,jj  ,jm  ) + spz_mod(ji  ,jj  ,jm  )) )      ! uw on uw grid
                  var_ten(ji, jj, 1, 6) =  var_ten(ji  ,jj  , 1,6) +   0.25_wp * var_tco(idt)        *     &
                         &              ( (spy_mod(ji  ,jj  ,jm  )                            )      *     &
                         &                (spz_mod(ji  ,jj+1,jm  ) + spz_mod(ji  ,jj  ,jm  )) )      ! vw on vw grid
               END DO
            END DO
         END DO
         IF (ln_wnd) THEN
            DO jm = 1, nn_tlu_wmod
               idt = jm 
               DO jj = 2, jpjm1 
                  DO ji = 2, jpim1
                     ! Components in z that requires backwards interpolation: wind contibution
                     var_ten(ji, jj, 1, 5) =  var_ten(ji  ,jj  , 1,5) +   0.25_wp * var_tco(idt)        *     &
                            &              ( (wsx_mod(ji  ,jj  ,jm  )                            )      *     &
                            &                (spz_mod(ji+1,jj  ,jm  ) + spz_mod(ji  ,jj  ,jm  )) )      ! uw on uw grid
                     var_ten(ji, jj, 1, 6) =  var_ten(ji  ,jj  , 1,6) +   0.25_wp * var_tco(idt)        *     &
                            &              ( (wsy_mod(ji  ,jj  ,jm  )                            )      *     &
                            &                (spz_mod(ji  ,jj+1,jm  ) + spz_mod(ji  ,jj  ,jm  )) )      ! vw on vw grid
                  END DO
               END DO
            END DO
         END IF
         !                             ! =====================
         !                             !   End of first slab
      !                                ! =====================

      !
      ! Inner layers
      !
      !                                ! ===============
      DO jk = 2, jpkm1                 ! Horizontal slab
         !                             ! ===============
         dkp1 = nmz(jk  )
         dkm1 = nmz(jk-1)
         !
         stt = midz(jk  )
         stp = midz(jk  ) + nmz(jk)
         !
         DO jm = stt, stp
            idt = jm - stt + 1
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  ! 
                  ! a_{xx}(x_{i},y,z) = 0.5 * \left[  a_{xx}(x_{i+1/2},y,z) +  a_{xx}(x_{i-1/2},y,z) \right]
                  !
                  !   with  a_{xx}(x,y,z) = \sum_{n=1)^{N} \phi^{(x)}_{n}(x,y,z) \phi^{(x)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, jk, 1) =  var_ten(ji  ,jj  ,jk  ,1) +    0.5_wp * var_tco(idt)       *     &
                         &               ( (spx_mod(ji  ,jj  ,jm  ) * spx_mod(ji  ,jj  ,jm  ))       +     &
                         &                 (spx_mod(ji-1,jj  ,jm  ) * spx_mod(ji-1,jj  ,jm  )) )
                  ! 
                  ! a_{yy}(x,y_{i},z) = 0.5 * \left[  a_{yy}(x,y_{j+1/2},z) +  a_{yy}(x,y_{j-1/2},z) \right]
                  !
                  !   with  a_{yy}(x,y,z) = \sum_{n=1)^{N} \phi^{(y)}_{n}(x,y,z) \phi^{(y)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, jk, 2) =  var_ten(ji  ,jj  ,jk  ,2) +    0.5_wp * var_tco(idt)       *     &
                         &               ( (spy_mod(ji  ,jj  ,jm  ) * spy_mod(ji  ,jj  ,jm  ))       +     &
                         &                 (spy_mod(ji  ,jj-1,jm  ) * spy_mod(ji  ,jj-1,jm  )) )
                  ! 
                  ! a_{zz}(x,y,z_{k-1/2}) = \sum_{n=1)^{N(k-1/2)} \phi^{(z)}_{n}(x,y,z) \phi^{(z)}_{n}(x,y,z)
                  ! 
                  var_ten(ji, jj, jk, 3) =  var_ten(ji  ,jj  ,jk  ,3) +   var_tco(idt)       *     &
                         &              ( (spz_mod(ji  ,jj  ,jm + dkp1) * spz_mod(ji  ,jj  ,jm + dkp1))   +     &
                         &                (spz_mod(ji  ,jj  ,jm       ) * spz_mod(ji  ,jj  ,jm       )) )
         
                  var_ten(ji, jj, jk, 4) =  var_ten(ji  ,jj  ,jk  ,4) +    0.25_wp * var_tco(idt)     *     &
                         &               ( (spx_mod(ji  ,jj+1,jm  )   + spx_mod(ji  ,jj  ,jm  ))      *     &
                         &                 (spy_mod(ji+1,jj  ,jm  )   + spy_mod(ji  ,jj  ,jm  )) )    ! uv on f grid        

                  ! Components in z that requires backwards interpolation: Here one could decide to do the same procedure 
                  !    used for the velocity and compute some temporal modes for the wind stress.

                  var_ten(ji, jj, jk, 5) =  var_ten(ji  ,jj  ,jk  ,5) +   0.25_wp * var_tco(idt)        *     &
                         &               ( (spx_mod(ji  ,jj  ,jm  )   + spx_mod(ji  ,jj  ,jm - dkm1))   *     &
                         &                 (spz_mod(ji+1,jj  ,jm  )   + spz_mod(ji  ,jj  ,jm       )) ) ! uw on uw grid
                  var_ten(ji, jj, jk, 6) =  var_ten(ji  ,jj  ,jk  ,6) +   0.25_wp * var_tco(idt)        *     &
                         &               ( (spy_mod(ji  ,jj  ,jm  )   + spy_mod(ji  ,jj  ,jm - dkm1))   *     &
                         &                 (spz_mod(ji  ,jj+1,jm  )   + spz_mod(ji  ,jj  ,jm       )) ) ! vw on vw grid
               END DO
            END DO
         END DO

         DO jj = 2, jpjm1 
            DO ji = 2, jpim1
               var_ten(ji, jj, jk-1, 3) = 0.5_wp * ( var_ten(ji, jj, jk-1, 3) + var_ten(ji, jj, jk, 3) )! w on T grid      
            END DO
         END DO
         !                             ! ===============
      END DO                           !   End of slab
      !                                ! ===============
      var_ten(:,:,jpkm1,3) = 0.5_wp * var_ten(:,:,jpkm1,3) 
      !
      ! Multiplication by dt
      !
      var_ten(:,:,:,:) = var_ten(:,:,:,:) * dt
      ! 
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,1:3) , 'T', -1.)
      CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,4)   , 'F', -1.)
      CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,5)   , 'U', -1.)
      CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,6)   , 'V', -1.)

   END SUBROUTINE tlu_var_pod



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_isd ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     kcmp: Define velocity type for calculation
   !> @param[in]     cdtype: Define velocity type
   !> @param[inout]  pcn, pca:  modified advection term for this component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_isd
   ! [tlu_isd]
   SUBROUTINE tlu_isd
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      uisd_n = 0._wp
      visd_n = 0._wp
      wisd_n = 0._wp
      !                                ! ================
      DO jk = 1, jpkm1                 ! Horizontal slab
         !                             ! ================
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               ! First component
               !
               ! dx(a11)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) + r1_e1e2u(ji,jj) *                           & 
                              & ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11)    &
                              & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) )  &
                              & * umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               ! dy(a21)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) + r1_e1e2u(ji,jj) *                           & 
                              & ( e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia12)    &
                              & - e1f(ji,jj-1) * e3f_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia12) )  &
                              & * umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               ! dz(a31)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) +                                             & 
                              & ( var_ten(ji,jj,jk  ,ia13) - var_ten(ji,jj,jk+1,ia13) ) *        &
                              &   umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               !
               ! Second component
               !
               ! dx(a12)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) + r1_e1e2v(ji,jj) *                           & 
                              & ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia12)    &
                              & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia12) )  &
                              & * vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               ! dy(a22)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) + r1_e1e2v(ji,jj) *                           & 
                              & ( e1f(ji,jj+1) * e3f_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia22)    &
                              & - e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia22) )  &
                              & * vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               ! dz(a32)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) +                                             & 
                              & ( var_ten(ji,jj,jk  ,ia23) - var_ten(ji,jj,jk+1,ia23) ) *        &
                              &   vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               !
               ! Third component
               !
               ! dx(a13)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) + r1_e1e2t(ji,jj) *                           & 
                              & ( e2t(ji  ,jj) * e3uw_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia13)   &
                              & - e2t(ji-1,jj) * e3uw_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia13) ) &
                              & * wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
               ! dy(a23)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) + r1_e1e2t(ji,jj) *                           & 
                              & ( e1t(ji,jj  ) * e3vw_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia23)   &
                              & - e1t(ji,jj-1) * e3vw_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia23) ) &
                              & * wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            !
            ! Third component Surface boundary condition
            !
            ! dz(a33)
            wisd_n(ji,jj, 1) =  wisd_n(ji,jj,1)                                                  & 
                           & - var_ten(ji,jj,1,ia33) * wmask(ji,jj,1) / e3w_n(ji,jj,1)
         END DO
      END DO

      !                                ! ================
      DO jk = 2, jpkm1                 ! Horizontal slab
         !                             ! ================
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               ! Third component
               !
               ! dz(a33)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) +                                            & 
                              & ( var_ten(ji,jj,jk  ,ia33) - var_ten(ji,jj,jk+1,ia33) ) *       &
                              &   wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      uisd_n = 0.5_wp * uisd_n 
      visd_n = 0.5_wp * visd_n
      wisd_n = 0.5_wp * wisd_n
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_isd', uisd_n , 'U', 1., visd_n , 'V', 1., wisd_n , 'T', 1.  )
      !
   END SUBROUTINE tlu_isd
   ! [tlu_isd]

   SUBROUTINE tlu_norm_tke ( kt, kene, norm_typ, prt_noi, kene_hr   )
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_svd_tke  ***
      !!
      !! ** Purpose :   Normalise noise and variance of LU on basis of TKE from energy spectrum
      !!
      !! ** Method :    Given the energy of small-scales from data, calculate 0.5*trace(var)
      !!                Calculate ratio of the two energy
      !!                Scale the noise and variance by calculated ratio
      !!
      !! ** Action :    var_ten : multiplied by cof_tke^2.*dt (now a diffusion tensor)
      !!                *noi : multiplied by cof_tke*dt (now a random velocity noise)
      !!                norm_typ: 1 - POD stationary noise - same coeff noise/var
      !!                          2 - POD stationary noise - HR var/ LR noise
      !!                          3 - SVD nonstationary noise (need to limit)
      !!                          4 - MM nonstationary noise with no scaling
      !!                          5 - MM nonstationary noise with tke scaling
      !!
      !!------------------------------------------------------------------------
      INTEGER                           , INTENT(in   ) ::   kt         ! ocean time-step index
      LOGICAL                           , INTENT(in   ) ::   prt_noi    ! ocean time-step index
      REAL(wp)                          , INTENT(in   ) ::   kene       ! the time-step
      REAL(wp) ,optional                , INTENT(in   ) ::   kene_hr
      REAL(wp) ::   zwt1,zwt2       ! local scalar
      REAL(wp), DIMENSION(9) :: zmaxl, zmaxg
      REAL(wp), DIMENSION(mppsize) :: ztau
      INTEGER  ::   ierr, norm_typ                   ! dummy loop argument
      !!------------------------------------------------------------------------
      !

      IF (norm_typ .ne. 4) THEN
        IF (norm_typ .eq. 1 .AND. kt .ne. nit000) THEN
                !empty to loop to avoid recalculation for stationary variance tensor
        ELSEIF (norm_typ .eq. 2 .AND. kt .ne. nit000) THEN
                !empty to loop to avoid recalculation for stationary variance tensor
        ELSE
          zwt1 = 0.
          zwt1 =  ( 0.5_wp * SUM (var_ten(:,:,:,1:3)) )/( jpi*jpj*jpk )   ! [TODO]

          call MPI_BARRIER(mpi_comm_oce, ierr)
          call MPI_ALLGATHER( zwt1, 1, mpi_double_precision, ztau, 1, mpi_double_precision, mpi_comm_oce, ierr)
          call MPI_BARRIER(mpi_comm_oce, ierr)
          zwt1 = SUM(ztau)/mppsize      !Tke as calculated from the variance tensor
          cof_tke = (kene/zwt1)**(0.5)   !Ratio of tke of small scales and tke of variance tensor
        ENDIF
      ENDIF
      
      !Scale the variables with the calculation ratio cof_tke
      IF (norm_typ .eq. 1) THEN !POD stationary noise
          unoi(:,:,:) = unoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          vnoi(:,:,:) = vnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          wnoi(:,:,:) = wnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          IF (kt == nit000) var_ten(:,:,:,:) = var_ten(:,:,:,:) * (cof_tke**2.) !*(ls_tlu**(7./3.))
      ELSEIF (norm_typ .eq. 2) THEN !POD stationary noise with higher dissipation
          unoi(:,:,:) = unoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          vnoi(:,:,:) = vnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          wnoi(:,:,:) = wnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          IF (kt == nit000) zwt2 = (kene_hr/zwt1)**(0.5)   !Ratio of tke of small scales and tke of variance tensor
          IF (kt == nit000) var_ten(:,:,:,:) = var_ten(:,:,:,:) * (zwt2**2.) !*(ls_tlu**(7./3.))
      ELSEIF (norm_typ .eq. 3) THEN !SVD Noise with limite
          IF (cof_tke .gt. 10) cof_tke = 10
!          unoi(:,:,:) = unoi(:,:,:)*cof_tke*pdt*(nn_tlu_psiz**(-1./3.))
!          vnoi(:,:,:) = vnoi(:,:,:)*cof_tke*pdt*(nn_tlu_psiz**(-1./3.))
!          wnoi(:,:,:) = wnoi(:,:,:)*cof_tke*pdt*(nn_tlu_psiz**(-1./3.))
!          var_ten(:,:,:,:) = var_ten(:,:,:,:)*(cof_tke**2.)*pdt*(nn_tlu_psiz**(-2./3.))
      ELSEIF (norm_typ .eq. 4) THEN !MM nonstationary noise
          unoi(:,:,:) = unoi(:,:,:)!*(ls_tlu**(7./3.))
          vnoi(:,:,:) = vnoi(:,:,:)!*(ls_tlu**(7./3.))
          wnoi(:,:,:) = wnoi(:,:,:)!*(ls_tlu**(7./3.))
          var_ten(:,:,:,:) = var_ten(:,:,:,:)!*(ls_tlu**(7./3.))
      ELSEIF (norm_typ .eq. 5) THEN !MM nonstationary noise
          unoi(:,:,:) = unoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          vnoi(:,:,:) = vnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          wnoi(:,:,:) = wnoi(:,:,:)*cof_tke!*(ls_tlu**(7./3.))
          var_ten(:,:,:,:) = var_ten(:,:,:,:) * (cof_tke**2.) !*(ls_tlu**(7./3.))
      ENDIF

      IF (prt_noi) THEN
        zmaxl(1) = maxval(unoi); zmaxl(2) = maxval(vnoi); zmaxl(3) = maxval(wnoi);
        zmaxl(4) = maxval(var_ten(:,:,:,1)); zmaxl(5) = maxval(var_ten(:,:,:,2)); zmaxl(6) = maxval(var_ten(:,:,:,3));
        zmaxl(7) = maxval(var_ten(:,:,:,4)); zmaxl(8) = maxval(var_ten(:,:,:,5)); zmaxl(9) = maxval(var_ten(:,:,:,6));
    
        call MPI_BARRIER(mpi_comm_oce, ierr)
        call MPI_REDUCE(zmaxl, zmaxg, 9, mpi_double_precision, MPI_MAX, 0, mpi_comm_oce, ierr)
        IF (lwp .AND. mod(kt,10) .eq. 0) THEN
            print *, '    For time (with coeff): ', kt, cof_tke
            print *, '    Max noise: ', zmaxg(1),zmaxg(2),zmaxg(3)
            print *, 'Max trace var: ', zmaxg(4),zmaxg(5),zmaxg(6)
            print *, 'Max cross var: ', zmaxg(7),zmaxg(8),zmaxg(9)
        ENDIF
      ENDIF

   END SUBROUTINE tlu_norm_tke

   SUBROUTINE  r8vec_normal_01 ( n, x )
   !!-------------------------------------------------------------------------------
   !!
   !! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
   !!
   !!  Discussion:
   !!
   !!    An R8VEC is an array of double precision real values.
   !!
   !!    The standard normal probability distribution function (PDF) has
   !!    mean 0 and standard deviation 1.
   !!
   !!  Licensing:
   !!
   !!    This code is distributed under the GNU LGPL license.
   !!
   !!  Modified:
   !!
   !!    18 May 2014
   !!
   !!  Author:
   !!
   !!    John Burkardt
   !!
   !!  Parameters:
   !!
   !!    Input, integer N, the number of values desired.
   !!
   !!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
   !!
   !!  Local:
   !!
   !!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
   !!    random values.  Its dimension is N+1, but really it is only needed
   !!    to be the smallest even number greater than or equal to N.
   !!
   !!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range
   !!    of entries of X that we need to compute
   !!
   !!-------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Dummy variables
      INTEGER,  INTENT( in  )                 :: n
      REAL(wp), INTENT( out ), DIMENSION( n ) :: x 
      ! Local variables
      INTEGER                                 :: m
      INTEGER                                 :: x_lo_index
      INTEGER                                 :: x_hi_index
      REAL(wp),                DIMENSION(n+1) :: r
      REAL(wp), PARAMETER                     :: r8_pi = 4._wp * atan(1._wp)

      !
      !  Record the range of X we need to fill in.
      !
      x_lo_index = 1
      x_hi_index = n
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
      IF ( x_hi_index - x_lo_index + 1 == 1 ) THEN

         CALL random_number ( harvest = r(1:2) )

         x(x_hi_index) = SQRT ( - 2._wp * LOG ( r(1) ) ) * COS ( 2._wp * r8_pi * r(2) )
      !
      !  If we require an even number of values, that's easy.
      !
      ELSE IF ( MOD ( x_hi_index - x_lo_index, 2 ) == 1 ) THEN

         m = ( x_hi_index - x_lo_index + 1 ) / 2

         call random_number ( harvest = r(1:2*m) )

         x(x_lo_index:x_hi_index-1:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-1:2) ) ) &
                                    & *  COS (   2._wp * r8_pi * r(2:2*m:2) )

         x(x_lo_index+1:x_hi_index:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-1:2) ) ) &
                                    & *  SIN (   2._wp * r8_pi * r(2:2*m:2) )
      !
      !  If we require an odd number of values, we generate an even number,
      !  and handle the last pair specially, storing one in X(N), and
      !  saving the other for later.
      !
      ELSE

         x_hi_index = x_hi_index - 1

         m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

         CALL Random_number ( harvest = r(1:2*m) )

         x(x_lo_index:x_hi_index-1:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-3:2) ) ) &
                                    & *  COS (   2._wp * r8_pi * r(2:2*m-2:2) )

         x(x_lo_index+1:x_hi_index:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-3:2) ) )  &
                                    & * SIN  (   2._wp * r8_pi * r(2:2*m-2:2) )

         x(n) = SQRT ( - 2._wp * LOG ( r(2*m-1) ) ) * COS ( 2._wp * r8_pi * r(2*m) )

      END IF

      RETURN
   END SUBROUTINE 


   SUBROUTINE tlu_nmz(  )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE   ***
      !!
      !! ** Purpose :  
      !!
      !! ** Method  :   
      !!                
      !!                
      !!------------------------------------------------------------------------
      INTEGER                             :: ji, jj
      integer, allocatable                ::   supp(:)
      real(wp) :: r
      INTEGER(i4) ::  ierr(2) 

 
      ALLOCATE( supp(nn_tco_dat),       &
             &  eigens(nn_tco_dat, jpk-1), stat=ierr(1) )
      OPEN(131, FILE=TRIM(nm_mod_nc),FORM='FORMATTED', STATUS = 'OLD')
      !
      ! Read Eigenvalues and cumulate:
      !      eigens(k) = \sum_{j=1}^{k} \sqrt(\lambda_{j})
      !
      READ(131,*) (eigens(1,jj) , jj = 1, jpk-1) 
      DO ji=2, nn_tco_dat
         READ(131,*)(eigens(ji,jj) , jj = 1, jpk-1) 
         eigens(ji,:) = eigens(ji,:) + eigens(ji-1,:)
      ENDDO
      !
      ! Column by column (i.e. layer by layer) normalize by \sum_{k}^{N}\lambda_{k}
      !
      DO ji = 1, jpk-1
         supp(:) = 0
         eigens(:,ji) = eigens(:,ji)/eigens(nn_tco_dat,ji)
         !
         ! Set the treshold "tlu_lam_spec" 
         !
         WHERE (eigens(:,ji) <= tlu_lam_spc)  supp = 1
         nmz(ji) = sum( supp ) 
      END DO
      !
      ! Shifts the array down one position, to ensure the right amount of energy in
      ! those components requiring a vertical interpolation between modes
      !
      nmz = [nmz(1), nmz(1:jpk-1)]
      !
      ! Compute index vector
      !
      midz(1) = 1
      DO ji = 2, jpk
         midz(ji) = midz(ji - 1) + nmz(ji)
      END DO
      !
      ! Verifies the input energy
      !
      r = 0._wp
      DO ji = 1, jpk-1
         r = r + eigens(nmz(ji),ji)/ (jpk-1)
      END DO

      IF (lwp) THEN
         WRITE(numout,*) '                                reading data from  : ', TRIM(nm_mod_nc)
         WRITE(numout,*) '              percentage of energy (least) intake  = ', tlu_lam_spc
         WRITE(numout,*) '                 real percentage of energy intake  = ', r
      END IF

      CLOSE(131)
      DEALLOCATE( eigens, supp, stat=ierr(2) )

   END SUBROUTINE tlu_nmz


END MODULE tlusvd




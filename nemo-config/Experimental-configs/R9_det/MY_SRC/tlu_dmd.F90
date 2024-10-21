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
MODULE tludmd
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
!   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_init_dmd
   PUBLIC tlu_tcoef_dmd
   PUBLIC tlu_bia_dmd
   !
   INCLUDE 'mpif.h'
   !
   ! [tlu_kenspect]
   CHARACTER(lc), PUBLIC                                      :: nm_mod_nc             !> @public   Name of the NetCDF file containing modes
   ! [tlu_kenspect]
   !
   ! [tlu_nmodes]
   INTEGER(i4),   PUBLIC                                      :: nn_tlu_nmod_r         !> @public   Number of Random modes
   INTEGER(i4),   PUBLIC                                      :: nn_tlu_nmod_c         !> @public   Number of Correlated modes
   ! [tlu_nmodes]
   !
   ! [tlu_frequencies]
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: omega_r               !> @public   Frequency of Random modes
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: omega_c               !> @public   Frequency of Correlated modes
   ! [tlu_frequencies]
   !
   ! [tlu_randgauss] 
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: rbrwn_rv              !> @public   Gaussian Random Variables (real part)
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: ibrwn_rv              !> @public   Gaussian Random Variables (imaginary part)
   ! [tlu_randgauss]
   !
   ! [tlu_randgauss] 
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: rtime_coef            !> @public   Time coefficient for modes real part
   REAL(wp),      PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: itime_coef            !> @public   Time coefficient for modes imaginary part
   ! [tlu_randgauss]
   !
   ! [tlu_noisemodes] 
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmd_rxmd_r, dmd_ixmd_r  !> @public   U-velocity modes, real and imaginary
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmd_rymd_r, dmd_iymd_r  !> @public   V-velocity modes, real and imaginary
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmd_rzmd_r, dmd_izmd_r  !> @public   W-velocity modes, real and imaginary
   ! [tlu_noiseodes]
   !
   ! [tlu_noisemodes] 
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmd_rxmd_c, dmd_ixmd_c  !> @public   U-velocity modes, real and imaginary
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmd_rymd_c, dmd_iymd_c  !> @public   V-velocity modes, real and imaginary
   ! [tlu_noiseodes]
   !
   ! [tlu_noisemodes] 
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: umean                 !> @public   U-velocity modes, real and imaginary
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: vmean                 !> @public   V-velocity modes, real and imaginary
   REAL(wp),      PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wmean                 !> @public   W-velocity modes, real and imaginary
   ! [tlu_noiseodes]
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_init_dmd ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_init_dmd
   ! [tlu_init_dmd]
   SUBROUTINE tlu_init_dmd
      !
      INTEGER(i4) :: ios, ierr(4), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Read namelist: transport under location uncertainty parametrization
      !
      NAMELIST/namtlu_dmd/   nm_mod_nc,       &
                         &   nn_tlu_nmod_r,   &
                         &   nn_tlu_nmod_c,   & 
                         &   biaSIGN
      !
      ! Read namelist
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu_dmd, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu_dmd, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_dmd_init : TLU, Dynamical Mode Decomposition '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '      data-based noise model             ln_tlu_dmd = ', ln_tlu_dmd
         WRITE(numout,*) '     number of uncorrelated modes     nn_tlu_nmod_r = ', nn_tlu_nmod_r
         WRITE(numout,*) '     number of   correlated modes     nn_tlu_nmod_c = ', nn_tlu_nmod_c
         WRITE(numout,*) '        file opened to read modes         nm_mod_nc = ', nm_mod_nc
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      ALLOCATE(   dmd_rxmd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &           dmd_rymd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &           dmd_rzmd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &           dmd_ixmd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &           dmd_iymd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &           dmd_izmd_r(jpi,jpj,nn_tlu_nmod_r * jpk), &
      &              omega_r(        nn_tlu_nmod_r      ), stat=ierr(1) )
      !
      ! Allocate Brownian motion
      !
      ALLOCATE(   rbrwn_rv(nn_tlu_nmod_r), &
      &           ibrwn_rv(nn_tlu_nmod_r), &
      &         rtime_coef(nn_tlu_nmod_r), &
      &         itime_coef(nn_tlu_nmod_r),   stat=ierr(2) )  
      !
      ! Read bias (if flagged)
      !
      IF (ln_tlu_bia) THEN
         !
         ! Allocate Spatial modes
         !
         ALLOCATE( dmd_rxmd_c(jpi,jpj,nn_tlu_nmod_c * jpk), &
         &         dmd_rymd_c(jpi,jpj,nn_tlu_nmod_c * jpk), &
         &         dmd_ixmd_c(jpi,jpj,nn_tlu_nmod_c * jpk), &
         &         dmd_iymd_c(jpi,jpj,nn_tlu_nmod_c * jpk), &
         &            omega_c(        nn_tlu_nmod_c      ),  stat=ierr(3) )
         !
         ! Allocate Spatial modes
         !
         ALLOCATE( umean(jpi,jpj,jpk), &
         &         vmean(jpi,jpj,jpk), &
         &         wmean(jpi,jpj,jpk), stat=ierr(4) )
         !
      END IF
      !
      ! Read Baroclinic modes
      !
      CALL read_vel_mod( )
      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_dmd_init(): allocation failed = ', ierr   ! [TODO] should be ctmp1? instead of numout
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !
   END SUBROUTINE tlu_init_dmd
   ! [tlu_init_dmd]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE read_vel_mod ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this read_vel_mod
   ! [read_vel_mod]
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
      INTEGER ::   jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at

      INTEGER                                 :: ncid

      REAL(wp)  :: val_r, val_c

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'read_vel_mod : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      ! Get netCDF file ID
      !
      val_r = 1.0_wp
      val_c = 0.125_wp
      !
      dmd_rxmd_r = 0._wp
      dmd_rymd_r = 0._wp
      dmd_rzmd_r = 0._wp
      dmd_ixmd_r = 0._wp
      dmd_iymd_r = 0._wp
      dmd_izmd_r = 0._wp
      !
      !
      CALL readFileID ("read_vel_md", '', nm_mod_nc, ncID)
      !  
      DO jm = 1,nn_tlu_nmod_r
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Read real parts of Random components
         !
         CALL read3Df4D_REALvar ("read_vel_md", "umode_real_r", ncID, dmd_rxmd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         CALL read3Df4D_REALvar ("read_vel_md", "vmode_real_r", ncID, dmd_rymd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         CALL read3Df4D_REALvar ("read_vel_md", "wmode_real_r", ncID, dmd_rzmd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         !
         ! Read imaginary parts of Random component
         !
         CALL read3Df4D_REALvar ("read_vel_md", "umode_imag_r", ncID, dmd_ixmd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         CALL read3Df4D_REALvar ("read_vel_md", "vmode_imag_r", ncID, dmd_iymd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         CALL read3Df4D_REALvar ("read_vel_md", "wmode_imag_r", ncID, dmd_izmd_r(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
         !
         dmd_rxmd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_rxmd_r(:,:,m_idx + 1 : m_idx + jpk ) * umask * val_r
         dmd_rymd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_rymd_r(:,:,m_idx + 1 : m_idx + jpk ) * vmask * val_r
         dmd_rzmd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_rzmd_r(:,:,m_idx + 1 : m_idx + jpk ) * wmask * val_r
         dmd_ixmd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_ixmd_r(:,:,m_idx + 1 : m_idx + jpk ) * umask * val_r
         dmd_iymd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_iymd_r(:,:,m_idx + 1 : m_idx + jpk ) * vmask * val_r
         dmd_izmd_r(:,:,m_idx + 1 : m_idx + jpk ) = dmd_izmd_r(:,:,m_idx + 1 : m_idx + jpk ) * wmask * val_r
         !
      END DO
      !
      ! Read frequency vector
      !
      CALL read1D_REALvar ("read_vel_md", "omega_r", ncID, omega_r)
      !
      ! Read bias (if flagged)
      !
      IF (ln_tlu_bia) THEN
         !
         dmd_rxmd_c = 0._wp
         dmd_rymd_c = 0._wp
         dmd_ixmd_c = 0._wp
         dmd_iymd_c = 0._wp
         !
         DO jm = 1,nn_tlu_nmod_c
            !
            ! Define zero-th indexed mode
            !
            m_idx = ( jm - 1 ) * jpk 
            !
            ! Read real parts of Random components
            !
            CALL read3Df4D_REALvar ("read_vel_md", "umode_real_c", ncID, dmd_rxmd_c(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
            CALL read3Df4D_REALvar ("read_vel_md", "vmode_real_c", ncID, dmd_rymd_c(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
            !
            ! Read imaginary parts of Random component
            !
            CALL read3Df4D_REALvar ("read_vel_md", "umode_imag_c", ncID, dmd_ixmd_c(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
            CALL read3Df4D_REALvar ("read_vel_md", "vmode_imag_c", ncID, dmd_iymd_c(:,:,m_idx + 1 : m_idx + jpk ), (jm-1)*2 + 1 )
            !
            dmd_rxmd_c(:,:,m_idx + 1 : m_idx + jpk ) = dmd_rxmd_c(:,:,m_idx + 1 : m_idx + jpk ) * umask * val_c
            dmd_rymd_c(:,:,m_idx + 1 : m_idx + jpk ) = dmd_rymd_c(:,:,m_idx + 1 : m_idx + jpk ) * vmask * val_c
            dmd_ixmd_c(:,:,m_idx + 1 : m_idx + jpk ) = dmd_ixmd_c(:,:,m_idx + 1 : m_idx + jpk ) * umask * val_c
            dmd_iymd_c(:,:,m_idx + 1 : m_idx + jpk ) = dmd_iymd_c(:,:,m_idx + 1 : m_idx + jpk ) * vmask * val_c
            !
         END DO
         !
         ! Read frequency vector
         !
         CALL read1D_REALvar ("read_vel_md", "omega_c", ncID, omega_c)
         !
         ! Read velocity average
         !
         CALL read3Df4D_REALvar ("read_vel_md", "uco", ncID, umean, 1 )
         CALL read3Df4D_REALvar ("read_vel_md", "vco", ncID, vmean, 1 )
         !
         umean = umean * umask * 0.01_wp
         vmean = vmean * vmask * 0.01_wp
         !
      END IF
      ! 
   END SUBROUTINE read_vel_mod
   ! [read_vel_mod]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tcoef_dmd ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_tcoef_dmd
   ! [tlu_tcoef_dmd]
   SUBROUTINE tlu_tcoef_dmd( kt, rBtcos_iBtsin, rBtsin_iBtcos )
      USE tlurnd
      !
      INTEGER,                            INTENT(in   ) :: kt
      REAL(wp), DIMENSION(nn_tlu_nmod_r), INTENT(  out) :: rBtcos_iBtsin
      REAL(wp), DIMENSION(nn_tlu_nmod_r), INTENT(  out) :: rBtsin_iBtcos
      !
      INTEGER :: ierr
      !
      ! Draw k random (standard) gaussian distribution
      !
      CALL Standard_Gaussian( nn_tlu_nmod_r, rbrwn_rv ) !, allones = .FALSE. )
      CALL Standard_Gaussian( nn_tlu_nmod_r, ibrwn_rv ) !, allones = .FALSE. )
      CALL MPI_BCAST( rbrwn_rv, nn_tlu_nmod_r, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      CALL MPI_BCAST( ibrwn_rv, nn_tlu_nmod_r, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      !
      ! Compute the time coefficients
      !
      rBtcos_iBtsin = ( ( rbrwn_rv * cos( omega_r * REAL( kt, wp ) * rn_rdt ) )   &
      &             -   ( ibrwn_rv * sin( omega_r * REAL( kt, wp ) * rn_rdt ) ) ) * (+2._wp)
      rBtsin_iBtcos = ( ( rbrwn_rv * sin( omega_r * REAL( kt, wp ) * rn_rdt ) )   &
      &             +   ( ibrwn_rv * cos( omega_r * REAL( kt, wp ) * rn_rdt ) ) ) * (-2._wp) 
      !
   END SUBROUTINE tlu_tcoef_dmd    
   ! [tlu_tcoef_dmd]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_bia_dmd ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_bia_dmd
   ! [tlu_bia_dmd]
   SUBROUTINE tlu_bia_dmd( kt, freq )
      INTEGER,                            INTENT(in   ) :: kt
      REAL(wp), DIMENSION(nn_tlu_nmod_c), INTENT(in   ) :: freq
      INTEGER                                           :: jm
      INTEGER                                           :: m_idx
      !
      ubia_n = umean
      vbia_n = vmean
      !     
      DO jm = 1, nn_tlu_nmod_c
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
         !
         ubia_n(:,:,:) = ubia_n(:,:,:) + 2._wp * dmd_rxmd_c(:,:,m_idx + 1 : m_idx + jpk ) * cos( freq(jm) * REAL( kt, wp ) * rn_rdt ) & 
         &                             - 2._wp * dmd_ixmd_c(:,:,m_idx + 1 : m_idx + jpk ) * sin( freq(jm) * REAL( kt, wp ) * rn_rdt )
         !
         vbia_n(:,:,:) = vbia_n(:,:,:) + 2._wp * dmd_rymd_c(:,:,m_idx + 1 : m_idx + jpk ) * cos( freq(jm) * REAL( kt, wp ) * rn_rdt ) &
         &                             - 2._wp * dmd_iymd_c(:,:,m_idx + 1 : m_idx + jpk ) * sin( freq(jm) * REAL( kt, wp ) * rn_rdt )
         !
      ENDDO
      !
      ! Mask the noise to enforce No-Flux boundary conditions 
      !
      ubia_n = ubia_n * umask
      vbia_n = vbia_n * vmask
      !
   END SUBROUTINE tlu_bia_dmd
   ! [tlu_bia_dmd]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE readFileID ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this readFileID
   ! [readFileID]
   SUBROUTINE readFileID (subnam, iodir, ncnam, ncID)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: iodir    ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: ncnam    ! Parent subroutine name (for EH)
      INTEGER,             INTENT(  OUT) :: ncID     ! NetCDF file ID 

      ! Internal  
      INTEGER :: ncstat

      ! Open netCDF file 
      ncstat = nf_open (TRIM(iodir)//TRIM(ncnam), nf_nowrite, ncID) 
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//subnam//'] '//ncnam)

   END SUBROUTINE readFileID


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 2D variable 
   !!            
   !!
   !!  @details
   !!  
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE read1D_REALvar (subnam, varnam, ncID, varOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )               :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )               :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )               :: ncID     ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:) :: varOUT   ! NetCDF file ID 
      REAL (KIND = wp),    ALLOCATABLE,   DIMENSION(:) :: zwVar    ! Temporary variable
      INTEGER                                          :: varID    ! Variable ID
      INTEGER,DIMENSION(1)                             :: starts
      INTEGER,DIMENSION(1)                             :: counts
      ! Internal  
      INTEGER :: ncstat

      ! Check existence
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      starts(:) = 1
      counts(:) = 2*size(varOUT)

      ALLOCATE( zwVar(2*size(varOUT)) )

      ! Read variable
      ncstat = nf_get_vara_double(ncID, varID, starts, counts, zwVar)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      varOUT =  zwVar(1:2*size(varOUT):2)
      DEALLOCATE( zwVar )

   END SUBROUTINE read1D_REALvar



   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable from a 4D array
   !!            
   !!
   !!  @details
   !!  
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE read3Df4D_REALvar (subnam, varnam, ncID, varOUT, idt)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                   :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                   :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                   :: ncID     ! NetCDF file ID 
      INTEGER,             INTENT(IN   )                   :: idt      ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                              :: varID    ! Variable ID    

      ! Internal  
      INTEGER                            :: ncstat  
      INTEGER                            :: starts(4), counts(4)
      CHARACTER (LEN = lc)               :: xID      ! x-dimension ID
      CHARACTER (LEN = lc)               :: yID      ! y-dimension ID
      CHARACTER (LEN = lc)               :: zID      ! z-dimension ID

      ! Set starting indexes
      starts(1) = mig(1)
      starts(2) = mjg(1)
      starts(3) = 1
      starts(4) = idt

      ! Set ending indexes
      counts(1) = jpi
      counts(2) = jpj
      counts(3) = jpk
      counts(4) =       1

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, 'x', xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'x')
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, 'y', yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'y')
      ! Check existence of dimension z
      ncstat = nf_inq_dimid (ncID, 'z', zID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'z')

!      ! Read dimension x
!      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
!      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'x')
!      ! Read dimension y
!      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
!      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'y')
!      ! Read dimension z
!      ncstat = nf_inq_dimlen (ncID, zID , counts(3) )
!      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// 'z')


      ! Read variable
      ncstat = nf_get_vara_double (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read3Df4D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    
   !!            
   !!
   !!  @details
   !!  
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE handle_err (ncstat, fromst)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Subroutine arguments
      INTEGER,             INTENT(IN   )           :: ncstat
      CHARACTER (LEN = *), INTENT(IN   ), OPTIONAL :: fromst

!     fromst is an optional string indicating where the CALL came
!     from that caused the netCDF problem (e.g. subroutine name).

!     Routine which interprets errors from netCDF output functions,
!     prints them to standard output and then kills the whole run.
      IF ( ncstat  .ne.  NF_NOERR ) THEN
         IF ( PRESENT(fromst) ) THEN
            PRINT *, TRIM(fromst)//'  '//TRIM( nf_strerror(ncstat) )
         ELSE
            PRINT *, TRIM( nf_strerror(ncstat) )
         ENDIF
         STOP 'netCDF:: STOPPED'
      ENDIF

   END SUBROUTINE handle_err


END MODULE tludmd





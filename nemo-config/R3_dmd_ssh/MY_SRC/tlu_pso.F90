MODULE tlupso
   !!===========================================================================
   !!                       ***  MODULE  tlupso  ***
   !! Transport under Location Uncertainty: stochastic noise SVD
   !!                       noise generation based on pseudo-observations and SVD.
   !!===========================================================================
   !! History : 0.0  !  2017-..  (P. DERIAN)  Original code
   !!
   !! [TODO]    - read namelist from _cfg, not only _ref
   !!           - write error messages to appropriate unit
   !!           - use NEMO's work arrays?
   !!---------------------------------------------------------------------------
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
   USE tluflt
   !   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_init_pso    ! called by tlu_init in tlu.f90
   PUBLIC gen_vel_mod     ! called by tlu_fields in tlu_noi.F90
   !
   INCLUDE 'mpif.h'
   !
   ! [tlu_ken_modes]
   INTEGER(i4), PUBLIC                                      :: nn_tlu_psiz =  3      !> @public size of patch
   INTEGER(i4), PUBLIC                                      :: nn_tlu_hsiz =  1      !> @public size of halo
   INTEGER(i4), PUBLIC                                      :: nn_tlu_pmes =  9      !> @public measure of patch
   INTEGER(i4), PUBLIC                                      :: nn_tlu_pobs           !> @public number of pseudo-observations
   LOGICAL,     PUBLIC                                      :: nn_tlu_exft = .TRUE.  !> @public extra filtering
   INTEGER(i4), PUBLIC                                      :: nn_tlu_fwdt =  3      !> @public width of filter ( odd number ) 
   ! [tlu_ken_modes]   
   !
   ! [tlu_spmodes] 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: pso_xmd           !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: pso_ymd           !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: pso_zmd           !> @public   W-velocity modes
   ! Bias
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psx_bia           !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psy_bia           !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psz_bia           !> @public   W-velocity modes
   ! [tlu_spmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),    PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: gaus_rv           !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_psob] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: u_psob            !> @private  U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: v_psob            !> @private  V-velocity pseudo-observations
   real(wp),            allocatable, save, dimension(:,:)   :: u_temp            !> @private  u-velocity pseudo-observations
   real(wp),            allocatable, save, dimension(:,:)   :: v_temp            !> @private  v-velocity pseudo-observations
   real(wp),            allocatable, save, dimension(:,:)   :: u_psot            !> @private  u-velocity pseudo-observations
   real(wp),            allocatable, save, dimension(:,:)   :: v_psot            !> @private  v-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: u_mean            !> @private  U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: v_mean            !> @private  V-velocity pseudo-observations 
   ! [tlu_psob]
   !
   ! [tlu_eigen] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Cglo              !> @private  U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Cloc              !> @private  U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Evec              !> @private  V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Tmod              !> @private  W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: Eval              !> @private  W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: work
   INTEGER,                          SAVE                   :: lwork
   ! [tlu_eigen]  
   !
   ! [tlu_area]
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     ::    dAu,    dAv    !> @private  U,V-point reshaped area factor
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     ::    dVu,    dVv    !> @private  U,V-point reshaped volume factor
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: rdA_Au, rdA_Av    !> @private  U,V-point dA/sum(dA) of the patch
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     ::  r1_Au,  r1_Av    !> @private  U,V-point  1/sum(dA) of the patch
   ! [tlu_area]
   !
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: real_unif_01      !> @private
   INTEGER(i4),         ALLOCATABLE,       DIMENSION(:,:)   :: pick_id           !> @private
   ! [tlu_randgauss] 
   !
   INTEGER(i4)                                              :: extra_halo        !> @private extra halo for u and v velocity
   INTEGER(i4)                                              :: jpi_ext, jpj_ext  !> @private extra halo max dimension 
   INTEGER(i4)                                              :: fih, fjh          !> @private first index after extra halo
   INTEGER(i4)                                              :: lih, ljh          !> @private last index begore extra halo
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: un_ext, vn_ext    !> @private working vectors
   !
   ! [tlu_filtering]
   INTEGER(i4)                                              :: flt_hsiz          !> @private extra halo for u and v velocity
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: un_flt, vn_flt    !> @private working vectors
   ! [tlu_filtering]

   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: e1e2u_ext, e1e2v_ext    !> @private working vectors
   !
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk1, wrk2        !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk3, wrk4        !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk5, wrk6        !> @private working vectors
   !
   !
   TYPE(filter)                                             :: boxfilter
   !
CONTAINS

   SUBROUTINE tlu_init_pso
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_init_pso  ***
      !!
      !! ** Purpose :   Initialize module variables.
      !!
      !! ** Method  :   Read namelist, compute other parameters and allocate arrays.
      !!
      !! ** Action  :   Arrays allocated.
      !!
      !!------------------------------------------------------------------------
      INTEGER(i4) :: ios, ierr(9), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Read namelist: transport under location uncertainty parametrization
      !
      NAMELIST/namtlu_svd/   nn_tlu_pobs,  &
      &                      biaSIGN
      !
      ! Read namelist
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu_svd, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu_svd, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_svd_init : TLU, Pseudo-Observation noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '           data-free noise model        ln_tlu_pso = ', ln_tlu_pso
         WRITE(numout,*) '   size of halo arond obervation       nn_tlu_hsiz = ', nn_tlu_hsiz
         WRITE(numout,*) '              width of the patch       nn_tlu_psiz = ', nn_tlu_psiz
         WRITE(numout,*) '  (integer) measure of the patch       nn_tlu_pmes = ', nn_tlu_pmes
         WRITE(numout,*) '   number of pseudo-observations       nn_tlu_pobs = ', nn_tlu_pobs
      END IF
      !
      ! Sanity check on number of pseudo observations used (always bigger than the patch members number) 
      ! 
      IF ( nn_tlu_pobs .lt. nn_tlu_pmes) THEN
         nn_tlu_pobs = nn_tlu_pmes
         IF (lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   WARNING: number of pseudo observations given was smaller than '
            WRITE(numout,*) '            the (integer) patch measure and thus increased.'
            WRITE(numout,*) '   number of pseudo-observations       nn_tlu_pobs = ', nn_tlu_pobs
         END IF
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      !
      ALLOCATE( pso_xmd( jpi, jpj, nn_tlu_pobs * jpk ), &
      &         pso_ymd( jpi, jpj, nn_tlu_pobs * jpk ), &
      &         pso_zmd( jpi, jpj, nn_tlu_pobs * jpk ), stat=ierr(1) )   ! the singular values
      !
      ALLOCATE( psx_bia( jpi, jpj, jpk ), &
      &         psy_bia( jpi, jpj, jpk ), &
      &         psz_bia( jpi, jpj, jpk ), stat=ierr(1) )   ! the singular values
      pso_xmd = 0._wp
      pso_ymd = 0._wp
      pso_zmd = 0._wp
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( gaus_rv( nn_tlu_pobs ), stat=ierr(2) )
      !
      gaus_rv = 0._wp 
      !
      ! Allocate Internal variables 
      !
      ALLOCATE(  u_psob( nn_tlu_pobs, jpi * jpj * (jpk-1)              ), &
      &          v_psob( nn_tlu_pobs, jpi * jpj * (jpk-1)              ), &
      &          u_temp( nn_tlu_pobs, jpi * jpj * (jpk-1)              ), &
      &          v_temp( nn_tlu_pobs, jpi * jpj * (jpk-1)              ), &
      &          u_psot(              jpi * jpj * (jpk-1), nn_tlu_pobs ), &
      &          v_psot(              jpi * jpj * (jpk-1), nn_tlu_pobs ), stat=ierr(3) )   ! the singular values
      !
      ! Allocate Internal variables 
      !
      ALLOCATE(  u_mean(              jpi * jpj * (jpk-1) ), &
      &          v_mean(              jpi * jpj * (jpk-1) ), stat=ierr(4) )     
      !
      ! Allocate scale factors for online POD
      !
      ALLOCATE(     dAu(              jpi * jpj * (jpk-1) ), &
      &             dAv(              jpi * jpj * (jpk-1) ), &
      &             dVu(              jpi * jpj * (jpk-1) ), &
      &             dVv(              jpi * jpj * (jpk-1) ), &
      &           r1_Au(              jpi * jpj * (jpk-1) ), &
      &           r1_Av(              jpi * jpj * (jpk-1) ), &
      &          rdA_Au( nn_tlu_pmes, jpi * jpj * (jpk-1) ), &
      &          rdA_Av( nn_tlu_pmes, jpi * jpj * (jpk-1) ), stat=ierr(5) )
      !
      ! Allocate the unifor random choices for the patches
      !
      ALLOCATE( real_unif_01( nn_tlu_pobs, jpi * jpj ), &
      &              pick_id( nn_tlu_pobs, jpi * jpj ), stat=ierr(6) )   ! the singular values
      !
      ! Allocate Matrix for eigen problem
      !
      ALLOCATE( Cglo( nn_tlu_pobs, nn_tlu_pobs ), &
      &         Cloc( nn_tlu_pobs, nn_tlu_pobs ), &
      &         Evec( nn_tlu_pobs, nn_tlu_pobs ), &
      &         Tmod( nn_tlu_pobs, nn_tlu_pobs ), &
      &         Eval( nn_tlu_pobs),               stat=ierr(7) )
      !
      ! Allocate workspaces
      !
      ALLOCATE( wrk3( jpi  ,jpj  ,jpk  ), &
      &         wrk4( jpi  ,jpj  ,jpk  ), stat=ierr(8) )
      !
      ! Allocate extra halo workspaces
      !
      extra_halo = nn_tlu_hsiz + 1
      fih = extra_halo + 1
      fjh = extra_halo + 1
      lih = fih + jpim1 
      ljh = fjh + jpjm1


      jpi_ext = jpi + 2*extra_halo
      jpj_ext = jpj + 2*extra_halo
      !
      ALLOCATE( e1e2u_ext( jpi_ext  , jpj_ext  )       , &
      &         e1e2v_ext( jpi_ext  , jpj_ext  )       , & 
      &            un_ext( jpi_ext  , jpj_ext  , jpk  ), &
      &            vn_ext( jpi_ext  , jpj_ext  , jpk  ), & 
      &              wrk1( jpi      , jpj      , jpk-1), &
      &              wrk2( jpi      , jpj      , jpk-1), stat=ierr(9) )
      !
      IF (nn_tlu_exft) THEN
         !
         flt_hsiz = ( nn_tlu_fwdt - 1 ) / 2
         !
         ! Initialize filter
         !
         CALL boxfilter%kern_init( "box", flt_hsiz )       
         !
         ALLOCATE(   wrk5( jpi + flt_hsiz , jpj + flt_hsiz , jpk  ), &
         &           wrk6( jpi + flt_hsiz , jpj + flt_hsiz , jpk  ), stat=ierr(8) )
         !
         ALLOCATE( un_flt( jpi , jpj , jpk  ), &
         &         vn_flt( jpi , jpj , jpk  ), stat=ierr(9) )
         !
      END IF
      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         IF(lwp) WRITE(numout,*) ' tlu_pod_init(): allocation failed = ', ierr  
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !
      ! Horizontal scale factors does not change in time, can be set now
      !
      CALL compute_areas( nit000 + dt_delay, 'hor' )
      CALL compute_areas( nit000 + dt_delay, 'ver' )  
      !     
   END SUBROUTINE tlu_init_pso


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE gen_vel_mod ***
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
   ! @snippet this gen_vel_mod
   ! [gen_vel_mod]
   SUBROUTINE gen_vel_mod( kt )
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
      INTEGER ::   kt
      INTEGER ::   ncid

      REAL(wp)  :: tic(4), toc(4)

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'gen_vel_mod : TLU, Dynamical Mode Decomposition reading modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      ! Update scale factors for online POD
      !
      IF( .NOT.ln_linssh ) CALL compute_areas( kt, 'ver' )
      !
      ! Generate Pseudo-Observations
      !
      !call cpu_time(tic(1))
      CALL gen_psobs( kt )
      !call cpu_time(toc(1))
      !
      ! Build Correlation matrix
      !
      !call cpu_time(tic(2))
      CALL build_Cmat( kt )
      !call cpu_time(toc(2))
      !
      ! Solve eigenvalue problem
      !
      !call cpu_time(tic(3))
      CALL solve_Cmat( kt )
      !call cpu_time(toc(3))
      !
      ! Build spatial modes
      !
      !call cpu_time(tic(4))
      CALL build_spmodes( kt )
      !call cpu_time(toc(4))
      !
      !print *, toc(1)-tic(1), toc(2)-tic(2), toc(3)-tic(3), toc(4)-tic(4), toc(4)-tic(1)
      !CALL check_spmodes( kt )
      ! 
   END SUBROUTINE gen_vel_mod
   ! [gen_vel_mod]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE compute_areas ***
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
   ! @snippet this compute_areas
   ! [compute_areas]
   SUBROUTINE compute_areas( kt, grid )
      USE tluMPI
      INTEGER,                                  INTENT(in   ) :: kt                     ! ocean time-step index
      CHARACTER(len=3),                         INTENT(in   ) :: grid                 ! = G (Model indicator)  
      INTEGER ::   jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   i, j, k, jk, jjii
      INTEGER ::   rel_idx(nn_tlu_psiz)

      INTEGER                                 :: patch_size

      REAL(wp)  :: val_r, val_c

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'compute_areas : TLU, generation of pseudo-observations '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      SELECT CASE ( grid )
      !
      CASE( 'hor' )
         !
         ! Collect extra information
         !
         CALL add_extra_halo_2D( kt, e1e2u, e1e2u_ext, extra_halo, 'e1e2u' )
         CALL add_extra_halo_2D( kt, e1e2v, e1e2v_ext, extra_halo, 'e1e2v' )
         !
         ! Global area factors for POD
         !
         dAu = reshape( spread( e1e2u, 3, jpk-1), (/ jpi * jpj * (jpk-1) /) )
         dAv = reshape( spread( e1e2v, 3, jpk-1), (/ jpi * jpj * (jpk-1) /) )
         !
         ! Local area factors for patch average
         !
         rel_idx = [( i, i = - nn_tlu_hsiz, nn_tlu_hsiz )]
         !
         rdA_Au = 0._wp
         rdA_Av = 0._wp
         !
         k = 1 
         !
         DO i = 1, nn_tlu_psiz
            !
            DO j = 1, nn_tlu_psiz
               !
               wrk1 = 0._wp
               wrk2 = 0._wp
               !
               wrk1( :, :, :) = spread( e1e2u_ext( fih+rel_idx(i):lih+rel_idx(i) ,           &
               &                                   fjh+rel_idx(j):ljh+rel_idx(j) ), 3, jpk-1 )
               wrk2( :, :, :) = spread( e1e2v_ext( fih+rel_idx(i):lih+rel_idx(i) ,           &
               &                                   fjh+rel_idx(j):ljh+rel_idx(j) ), 3, jpk-1 )
               !
               rdA_Au(k, :) = reshape( wrk1, (/ jpi * jpj * (jpk-1) /) )
               rdA_Av(k, :) = reshape( wrk2, (/ jpi * jpj * (jpk-1) /) )
               !
               k = k + 1
               !
            END DO
            !
         END DO
         !
         rdA_Au = rdA_Au / spread( SUM(rdA_Au, 1), 1, nn_tlu_pmes)
         rdA_Av = rdA_Av / spread( SUM(rdA_Av, 1), 1, nn_tlu_pmes)
         !
      CASE( 'ver' )
         !
         ! Recompute scale factors
         !
         dVu = dAu * reshape( e3u_n(:, :, 1:jpk-1), (/ jpi * jpj * (jpk-1) /) )
         dVv = dAv * reshape( e3v_n(:, :, 1:jpk-1), (/ jpi * jpj * (jpk-1) /) )
         !
      END SELECT
      !
   END SUBROUTINE compute_areas
   ! [compute_areas]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE gen_psobs ***
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
   ! @snippet this gen_psobs
   ! [gen_psobs]
   SUBROUTINE gen_psobs( kt )
      USE tluMPI
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


      INTEGER :: kt ! dummy loop arguments
      INTEGER ::   jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   i, j, k, ji, jj, jk, jjii, jjiikk
      INTEGER ::   rel_idx(nn_tlu_psiz)

      INTEGER                                 :: patch_size

      REAL(wp)  :: val_r, val_c

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'gen_psobs : TLU, generation of pseudo-observations '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      END IF
      !
      ! Filter the velocity if requested
      !
      IF (nn_tlu_exft) THEN
         !
         CALL add_extra_halo( kt, un, wrk5, flt_hsiz, 'un' )
         CALL add_extra_halo( kt, vn, wrk6, flt_hsiz, 'vn' )
         !
         ! Box filtering velocity
         !
         CALL boxfilter%kern_apply_2fld( wrk5, wrk6, un_flt, vn_flt )
         !
         ! Collect extra information
         !
         CALL add_extra_halo( kt, un_flt, un_ext, extra_halo, 'un' )
         CALL add_extra_halo( kt, vn_flt, vn_ext, extra_halo, 'vn' )
         ! 
      ELSE
         !
         ! Collect extra information
         !
         CALL add_extra_halo( kt, un, un_ext, extra_halo, 'un' )
         CALL add_extra_halo( kt, vn, vn_ext, extra_halo, 'vn' )
         !      
      END IF
      !
      ! Assign Indexes [ ..., -1, 0, +1, ...] and compute patch size
      !
      rel_idx = [( i, i = - nn_tlu_hsiz, nn_tlu_hsiz )]
      !
      !
      u_psob = 0._wp
      v_psob = 0._wp
      u_temp = 0._wp
      v_temp = 0._wp
      u_mean = 0._wp
      v_mean = 0._wp
      !
      call random_number( real_unif_01 )
      pick_id = INT( 1 + FLOOR( ( nn_tlu_pmes ) * real_unif_01 ) , i4)  
      !
      !
      jjiikk = 1
      !
      DO jk = 1, jpk -1
         !
         DO jj = fjh, ljh; DO ji = fih, lih
            ! 
            k = 1
            DO i = 1, nn_tlu_psiz; DO j = 1, nn_tlu_psiz
               !
               u_temp(k, jjiikk) = un_ext( ji + rel_idx(i), jj + rel_idx(j), jk )   
               v_temp(k, jjiikk) = vn_ext( ji + rel_idx(i), jj + rel_idx(j), jk )   
               !
               u_mean( jjiikk ) = u_mean( jjiikk ) + u_temp( k,  jjiikk ) * rdA_Au( k, jjiikk ) 
               v_mean( jjiikk ) = v_mean( jjiikk ) + v_temp( k,  jjiikk ) * rdA_Av( k, jjiikk )

               psx_bia(ji,jj,jk) = u_mean( jjiikk ) 
               psy_bia(ji,jj,jk) = v_mean( jjiikk ) 
               !    
               k = k + 1
               !
            END DO; END DO
            !
            k = 1
            DO i = 1, nn_tlu_psiz; DO j = 1, nn_tlu_psiz
               !
               u_temp(k, jjiikk) = u_temp(k, jjiikk) - u_mean( jjiikk ) 
               v_temp(k, jjiikk) = v_temp(k, jjiikk) - v_mean( jjiikk ) 
               !
               k = k + 1
               !
            END DO; END DO
            !
            jjiikk = jjiikk + 1
            !
         END DO; END DO; !END DO
         !
         ! Define zero-th indexed position
         !
         m_idx = ( jk - 1 ) * jpi * jpj 
         !
         DO jjii = 1, jpi * jpj
            !
            u_psob(:, m_idx + jjii) = u_temp( pick_id(:, jjii), m_idx + jjii)
            v_psob(:, m_idx + jjii) = v_temp( pick_id(:, jjii), m_idx + jjii)
            !
            u_psot(m_idx + jjii, :) = u_psob(:, m_idx + jjii) !u_psob * SPREAD( dVu, 1, nn_tlu_pobs) * dVu( k, jjiikk ) 
            v_psot(m_idx + jjii, :) = v_psob(:, m_idx + jjii) !v_psob * SPREAD( dVv, 1, nn_tlu_pobs) * dVu( k, jjiikk ) 

            !
         END DO
         !
      END DO
      !
   END SUBROUTINE gen_psobs
   ! [gen_psobs]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE build_Cmat ***
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
   ! @snippet this build_Cmat
   ! [build_Cmat]
   SUBROUTINE build_Cmat( kt )
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

      INTEGER :: kt ! dummy loop arguments
      INTEGER ::   jm, m_idx, ierr, rank ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at
      INTEGER ::   i, j
      INTEGER                                 :: ncid

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'build_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Cglo
      !
      Cloc = 0._wp
      !
!      u_psot = 0._wp! TRANSPOSE( u_psob * SPREAD( dVu, 1, nn_tlu_pobs) )
!      v_psot = 0._wp! TRANSPOSE( v_psob * SPREAD( dVv, 1, nn_tlu_pobs) )
      !
      Cloc = Cloc + MATMUL( u_psob, u_psot) / nn_tlu_pobs
      Cloc = Cloc + MATMUL( v_psob, v_psot) / nn_tlu_pobs
      !
      ! Communicate the values of Cglo across all processors (with sum)
      !
      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLREDUCE(Cloc, Cglo, nn_tlu_pobs**2, mpi_double_precision, MPI_SUM, mpi_comm_oce, ierr)
      !
   END SUBROUTINE build_Cmat
   ! [build_Cmat]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE solve_Cmat ***
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
   ! @snippet this solve_Cmat
   ! [solve_Cmat]
   SUBROUTINE solve_Cmat( kt )
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
      INTEGER :: kt ! dummy loop arguments

      real(wp) :: val(nn_tlu_pobs)
      real(wp) :: dummy(1)
      integer ::  nrA, ncA, nrB, ncB
      INTEGER              :: info
      INTEGER ::   i, j
      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'solve_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Eigenvalue matrix by constructing upper triangular of C
      !
      Evec = Cglo
      ! 
      ! Allocate the workspaces for LAPACK's dsyev
      !
      IF ( kt == nit000 + dt_delay ) THEN 
         !
         !    dsyev( JOBZ, UPLO,           N,    A,         LDA,    W, WORK, LWORK, INFO )
         CALL dsyev(  'V',  'L', nn_tlu_pobs, Evec, nn_tlu_pobs, Eval, dummy, -1, info )
         lwork = INT( dummy(1) )
         ALLOCATE( work( lwork ) )
         !
      END IF
      !
      ! Solve   C a=\lambda a
      !
      !    dsyev( JOBZ, UPLO,           N,    A,         LDA,    W, WORK, LWORK, INFO )
      CALL dsyev(  'V',  'L', nn_tlu_pobs, Evec, nn_tlu_pobs, Eval, work, lwork, info )
      !
      DO i = 1, nn_tlu_pobs
         !
         !
         Evec(:, i) = Evec(:, i) * SQRT( nn_tlu_pobs * Eval(i) )
         Tmod(i, :) = Evec(:, i) / ( nn_tlu_pobs * SQRT( Eval(i) ) ) 
         !
      END DO
      !
      ! Transpose components ( https://stackoverflow.com/questions/55222312/do-most-compilers-optimize-matmultransposea-b )
      !
      !nrA = SIZE( Evec, 1)
      !ncA = SIZE( Evec, 2)
      !nrB = SIZE( u_psob, 1)
      !ncB = SIZE( u_psob, 2)
      !     DGEMM( TRANSA, TRANSB,   M,   N,   K, ALPHA,    A, LDA,      B, LDB,  BETA,      C, LDC)
      !CALL DGEMM(    'T',    'N', nrA, ncB, ncA, 1._wp, Evec, nrA, u_psob, nrB, 0._wp, u_psob, nrB)
      !CALL DGEMM(    'T',    'N', nrA, ncB, ncA, 1._wp, Evec, nrA, v_psob, nrB, 0._wp, v_psob, nrB)
      !
      u_psob = matmul( Tmod, u_psob)
      v_psob = matmul( Tmod, v_psob)
      !
   END SUBROUTINE solve_Cmat
   ! [solve_Cmat]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE build_spmodes ***
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
   ! @snippet this build_spmodes
   ! [build_spmodes]
   SUBROUTINE build_spmodes( kt )
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
      
      INTEGER :: kt ! dummy loop arguments
      INTEGER ::   jo, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at

      INTEGER        :: ncid

      REAL(wp)  :: val

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'build_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      val = 0.25_wp! 2*SQRT(REAL(nn_tlu_pmes, wp))**(1._wp/3._wp)
      !
      DO jo = 1, nn_tlu_pobs
         !
         ! Initialize array
         !
         wrk3 = 0._wp
         wrk3 = 0._wp
         !
         ! Define zero-th indexed position (in target array, that has horizontal and vertical halo)
         !
         m_idx = ( jo - 1 ) * jpk 
         !
!         wrk3( :, :, 1:jpk-1) = 0._wp!reshape( u_psob( jo, :), (/ jpi, jpj, jpk-1 /) ) * umask * val
!         wrk4( :, :, 1:jpk-1) = 0._wp!reshape( v_psob( jo, :), (/ jpi, jpj, jpk-1 /) ) * vmask * val
         !
         ! Lateral boundary condition transfer across nodes
         !
!         CALL lbc_lnk_multi( 'build_vmodes', wrk3 , 'U', -1., wrk4 , 'V', -1.   )         
         !
         pso_xmd(:,:, m_idx + 1 : m_idx + jpk -1 ) = reshape( u_psob( jo, :), (/ jpi, jpj, jpk-1 /) ) * umask * val
         pso_ymd(:,:, m_idx + 1 : m_idx + jpk -1 ) = reshape( v_psob( jo, :), (/ jpi, jpj, jpk-1 /) ) * vmask * val
         !
         ! Build the vertical mode
         !
         CALL build_vmodes( pso_zmd(:,:,m_idx + 1 : m_idx + jpk ), &
         &                  pso_xmd(:,:,m_idx + 1 : m_idx + jpk ), &
         &                  pso_ymd(:,:,m_idx + 1 : m_idx + jpk )  )
         !
         ! Apply safety mask
         !
         pso_zmd(:,:,m_idx + 1 : m_idx + jpk ) = pso_zmd(:,:,m_idx + 1 : m_idx + jpk ) * wmask * val
         !
      END DO
      !      
   END SUBROUTINE build_spmodes
   ! [build_spmodes]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE build_vmodes ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Original code <br> 
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
   ! @snippet this build_vmodes
   ! [build_vmodes]
   SUBROUTINE build_vmodes(wout, uin, vin)
      REAL(wp), DIMENSION(:,:,:), INTENT(in   )  :: uin
      REAL(wp), DIMENSION(:,:,:), INTENT(in   )  :: vin
      REAL(wp), DIMENSION(:,:,:), INTENT(  out)  :: wout
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !
      !!----------------------------------------------------------------------
      !
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         wout( 2:jpi-1, 2:jpj-1, jk) =     wout( 2:jpi-1, 2:jpj-1, jk + 1) -   &
         !
         !  dx( dy*dz*u_in )
         &                           (      e2u( 2:jpi-1, 2:jpj-1 )            &
         &                           *    e3u_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                           *      uin( 2:jpi-1, 2:jpj-1, jk )        &
         &                           -      e2u( 1:jpi-2, 2:jpj-1 )            &
         &                           *    e3u_n( 1:jpi-2, 2:jpj-1, jk )        &
         &                           *      uin( 1:jpi-2, 2:jpj-1, jk )        &
         !  dy ( dx*dz*v_isd )
         &                           +      e1v( 2:jpi-1, 2:jpj-1 )            &
         &                           *    e3v_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                           *      vin( 2:jpi-1, 2:jpj-1, jk )        &
         &                           -      e1v( 2:jpi-1, 1:jpj-2 )            &
         &                           *    e3v_n( 2:jpi-1, 1:jpj-2, jk )        &
         &                           *      vin( 2:jpi-1, 1:jpj-2, jk ) )      &
         !  .../dV
         &                           *    tmask( 2:jpi-1, 2:jpj-1, jk)         & 
         &                           * r1_e1e2t( 2:jpi-1, 2:jpj-1 )
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      !
      ! Lateral boundary condition transfer across nodes
      !
     ! CALL lbc_lnk_multi( 'build_vmodes', wout , 'T', -1.  )
      !
   END SUBROUTINE build_vmodes   
   ! [tlu_wzvcmp]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE build_spmodes ***
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
   ! @snippet this check_spmodes
   ! [check_spmodes]
   SUBROUTINE check_spmodes( kt )
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
      
      INTEGER :: kt !
      INTEGER ::   jo, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   i, j, k

      INTEGER        :: ncid

      REAL(wp)  ::  eigenval(nn_tlu_pobs, nn_tlu_pobs)
      REAL(wp)  ::  delta(nn_tlu_pobs, nn_tlu_pobs)

      REAL(wp)  ::  u_psob_t_locl( jpi * jpj * jpk, nn_tlu_pobs)


      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'check_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      eigenval = 1._wp
      delta = 0._wp
      !
      eigenval = matmul(  TRANSPOSE( Evec ), Evec) / nn_tlu_pobs
      print *, " "
      print *, " "
      do i = 1, nn_tlu_pobs
         !
         print *, (eigenval(i, j) , j = 1, nn_tlu_pobs), "|", Eval(i)         !
      end do
      !
       !
      ! Initialize Cglo
      !
      Cloc = 0._wp
      !
      u_psot = TRANSPOSE( u_psob * SPREAD( dVu, 1, nn_tlu_pobs) )
      v_psot = TRANSPOSE( v_psob * SPREAD( dVv, 1, nn_tlu_pobs) )
      !
      Cloc = Cloc + MATMUL( u_psob, u_psot)
      Cloc = Cloc + MATMUL( v_psob, v_psot)
      !
      ! Communicate the values of Cglo across all processors (with sum)
      !
      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLREDUCE(Cloc, Cglo, nn_tlu_pobs**2, mpi_double_precision, MPI_SUM, mpi_comm_oce, ierr)
      !
      print *, " "
      print *, " "
      do i = 1, nn_tlu_pobs
         !
         print *, ( Cglo(i, j) / Eval(i)  , j = 1, nn_tlu_pobs)
         !
      end do

      stop
     !      
   END SUBROUTINE check_spmodes
   ! [check_spmodes]


END MODULE tlupso




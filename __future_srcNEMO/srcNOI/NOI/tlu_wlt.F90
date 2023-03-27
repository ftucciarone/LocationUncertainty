MODULE tluwlt
   !!===========================================================================
   !!                       ***  MODULE  tluwlt  ***
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
   !   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_init_wlt    ! called by tlu_init in tlu.f90
   PUBLIC gen_wlt_mod     ! called by tlu_fields in tlu_noi.F90
   !
   INCLUDE 'mpif.h'

   !
   ! [tlu_ken_modes]
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_nlev = 3    !: size of patch
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_hsiz = 1    !: size of patch
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_nlev        !: number of pseudo-observations
   ! [tlu_ken_modes]   
   !
   ! [tlu_spmodes] 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wlx_mod                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wly_mod                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wlz_mod                 !> @public   W-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wlx_noi                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wly_noi                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wlz_noi                 !> @public   W-velocity modes
   ! [tlu_spmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),    PUBLIC, ALLOCATABLE,       DIMENSION(:,:)   :: rnd_fld                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_wtsf] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: u_wtsf                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: v_wtsf                 !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: w_wtsf                 !> @public   W-velocity pseudo-observations
   ! [tlu_wtsf]
   !
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:)     :: wlt_coef               !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: msk_fld                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk1, wrk2

CONTAINS

   SUBROUTINE tlu_init_wlt
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_init_wlt  ***
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
      NAMELIST/namtlu_wlt/   nn_tlu_nlev,  &
      &                      nn_daub_ord,  &
      &                      biaSIGN
      !
      ! Read namelist
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu_wlt, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu_wlt, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_wlt_init : TLU, Wavelet decomposition noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '           wavelet decomposition        ln_tlu_wlt = ', ln_tlu_wlt
         WRITE(numout,*) '          number of detail level       nn_tlu_nlev = ', nn_tlu_nlev
         WRITE(numout,*) '        Daubechies wavelet order       nn_daub_ord = ', nn_daub_ord
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      !
      ALLOCATE( wlx_mod( jpi, jpj, jpk ), &
      &         wly_mod( jpi, jpj, jpk ), &
      &         wlz_mod( jpi, jpj, jpk ), &
      &         wlx_noi( jpi, jpj, jpk ), &
      &         wly_noi( jpi, jpj, jpk ), &
      &         wlz_noi( jpi, jpj, jpk ), stat=ierr(1) )   ! the singular values
      !
      wlx_mod = 0._wp
      wly_mod = 0._wp
      wlz_mod = 0._wp
      wlx_noi = 0._wp
      wly_noi = 0._wp
      wlz_noi = 0._wp
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( rnd_fld( jpi, jpj ), stat=ierr(2) )
      !
      rnd_fld = 0._wp 
      !
      ! Allocate mask
      !
      ALLOCATE( msk_fld( jpi, jpj ), stat=ierr(3) )
      !
      msk_fld = 0._wp       
      !
      ! Allocate and set Wavelet coefficients
      !
      CALL set_Daubechies ( nn_daub_ord )
      !
      ! Allocate Internal variables 
      !
      ALLOCATE( u_wtsf( jpi, jpj, jpk ), &
      &         v_wtsf( jpi, jpj, jpk ), &
      &         w_wtsf( jpi, jpj, jpk ), stat=ierr(4) )   ! the singular values
      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_wlt_init(): allocation failed = ', ierr  
         STOP
      END IF
      !
   END SUBROUTINE tlu_init_wlt


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE gen_wlt_mod ***
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
   ! @snippet this gen_wlt_mod
   ! [gen_wlt_mod]
   SUBROUTINE gen_wlt_mod( kt )
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

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'gen_wlt_mod : TLU, Wavelet-based generation of modes '
         WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      ! Direct Discrete Wavelet transform
      !
      DO jk = 1, jpk
         !
         ! One way 
         !
         DO jj = 1, jpj
            !
            CALL wlt_transform( un(:, jj, jk ), u_wtsf(:, jj, jk ), nn_tlu_nlev, wlt_coef, "dir")
            CALL wlt_transform( vn(:, jj, jk ), v_wtsf(:, jj, jk ), nn_tlu_nlev, wlt_coef, "dir")
            CALL wlt_transform( wn(:, jj, jk ), w_wtsf(:, jj, jk ), nn_tlu_nlev, wlt_coef, "dir")
            !
         END DO
         !
         ! Other way around         
         !
         DO ji = 1, jpi
            !
            CALL wlt_transform( u_wtsf(ji, :, jk ), u_wtsf(ji, :, jk ), nn_tlu_nlev, wlt_coef, "dir")
            CALL wlt_transform( v_wtsf(ji, :, jk ), v_wtsf(ji, :, jk ), nn_tlu_nlev, wlt_coef, "dir")
            CALL wlt_transform( w_wtsf(ji, :, jk ), w_wtsf(ji, :, jk ), nn_tlu_nlev, wlt_coef, "dir")
            !
         END DO
         !
      END DO    
      !
      ! Generation of the random gaussian field (it varies in space as well this time!)
      !
      rnd_fld = 0._wp
      msk_fld = 0._wp
      mx = jpi
      my = jpj
      !
      DO jl = 1, nn_tlu_nlev
         !
         rnd_fld(    1: mx/2,    1: my/2 ) = 0._wp
         rnd_fld( mx/2:   mx,    1: my/2 ) = 0._wp
         rnd_fld( mx/2:   mx, my/2:   my ) = 0._wp
         !
         msk_fld(    1: mx/2,    1: my/2 ) = 1._wp
         msk_fld( mx/2:   mx,    1: my/2 ) = 1._wp
         msk_fld( mx/2:   mx, my/2:   my ) = 1._wp
         !
         mx = mx / 2
         my = my / 2
         !
      END DO
      !
      ! Inverse Discrete Wavelet transform
      !
      DO jk = 1, jpk
         !
         ! Multiply by the mask and random field
         !
         u_wtsf(:,:,jk) = u_wtsf(:,:,jk) * msk_fld(:,:)
         u_wtsf(:,:,jk) = v_wtsf(:,:,jk) * msk_fld(:,:)
         u_wtsf(:,:,jk) = w_wtsf(:,:,jk) * msk_fld(:,:)
         !
         wlx_noi(:,:,jk) = u_wtsf(:,:,jk) * rnd_fld(:,:)
         wly_noi(:,:,jk) = v_wtsf(:,:,jk) * rnd_fld(:,:)
         wlz_noi(:,:,jk) = w_wtsf(:,:,jk) * rnd_fld(:,:)
         !
         ! One way 
         ! 
         DO jj = 1, jpj
            !
            ! Transform the variance
            !
            CALL wlt_transform(  u_wtsf(:, jj, jk ), wlx_mod(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform(  v_wtsf(:, jj, jk ), wly_mod(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform(  w_wtsf(:, jj, jk ), wlz_mod(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            !
            ! Transform the noise
            !
            CALL wlt_transform( wlx_noi(:, jj, jk ), wlx_noi(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wly_noi(:, jj, jk ), wly_noi(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wlz_noi(:, jj, jk ), wlz_noi(:, jj, jk ), nn_tlu_nlev, wlt_coef, "inv")
            !
         END DO
         !
         ! Other way around         
         !
         DO ji = 1, jpi
            !
            ! Transform the variance
            !
            CALL wlt_transform( wlx_mod(ji, :, jk ), wlx_mod(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wly_mod(ji, :, jk ), wly_mod(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wlz_mod(ji, :, jk ), wlz_mod(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            !
            ! Transform the noise
            !
            CALL wlt_transform( wlx_noi(ji, :, jk ), wlx_noi(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wly_noi(ji, :, jk ), wly_noi(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            CALL wlt_transform( wlz_noi(ji, :, jk ), wlz_noi(ji, :, jk ), nn_tlu_nlev, wlt_coef, "inv")
            !
         END DO
         !
      END DO
      !
   END SUBROUTINE gen_wlt_mod
   ! [gen_wlt_mod]




 

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE wlt_transform ***
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
   ! @snippet this wlt_transform
   ! [wlt_transform]
   SUBROUTINE wlt_transform(sig_in, sig_out, ndl, wlt_coef, direction)
      CHARACTER(3),               INTENT(in   ) :: direction
      INTEGER(4),                 INTENT(in   ) :: ndl
      REAL(wp),     DIMENSION(:), INTENT(in   ) :: wlt_coef
      REAL(wp),     DIMENSION(:), INTENT(in   ) :: sig_in
      REAL(wp),     DIMENSION(:), INTENT(  out) :: sig_out
      !
      REAL(wp),     DIMENSION(:)                :: wrk_sig
      INTEGER(4)                                :: ierr
      INTEGER(4)                                :: i, j, k, m, n, p, q, s
      !
      n = SIZE(sig_in)
      p = SIZE(wlt_coef) - 1
      !
      ! Dimension check
      !
      IF ( MOD( n, LSHIFT(1, ndl) ) .NE. 0 ) THEN
         WRITE(numerr, *) 'wlt_transform: number of levels (ndl) incompatible with given vector size (N)'
         STOP
      END IF
      !
      ! Allocate workspace
      !
      ALLOCATE( wrk_sig(n), stat=ierr)
      !
      sig_out = sig_in
      s = 0   ! the relative scale counter
      q = ( p - 1 ) / 2
      !
      !
      SELECT CASE( direction )
      !
      ! Direct transform
      !
      CASE( "dir" )
         !
         m = n
         !
         DO WHILE ( s < ndl )   ! for each decomposition level
            !
            i = 1
            wrk_sig(1:m) = 0._wp
            !
            DO j = 1, m - 1, 2
               !
               DO k = 0, p - 1, 2
                  j0 = i4_wrap ( j + k,     1, m )
                  j1 = i4_wrap ( j + k + 1, 1, m )
                  wrk_sig(i)     = wrk_sig(i)     + wlt_coef(  k) * sig_out(j0) + wlt_coef(  k+1) * sig_out(j1)
                  wrk_sig(i+m/2) = wrk_sig(i+m/2) + wlt_coef(p-k) * sig_out(j0) - wlt_coef(p-k-1) * sig_out(j1)
               END DO
               !
               i = i + 1
               !
            END DO
            !
            sig_out(1:m) = wrk_sig(1:m)
            m = m / 2
            s = s + 1
            !
         END DO
      !
      ! Inverse transform
      !
      CASE( "inv" )
         !
         m = 2 * rshift(n, ndl)
         !
         DO WHILE ( s < ndl )   ! for each decomposition level
            !
            j = 1
            wrk_sig(1:m) = 0._wp
            !
            do i = - q + 1, m / 2 - q
               !
               do k = 0, p - 1, 2
                  i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
                  i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
                  wrk_sig(j)   = wrk_sig(j)   + wlt_coef(p-k-1) * sig_out(i0) + wlt_coef(k+1) * sig_out(i1)
                  wrk_sig(j+1) = wrk_sig(j+1) + wlt_coef(p-k)   * sig_out(i0) - wlt_coef(k)   * sig_out(i1)
               end do
               !
               j = j + 2
               !
            end do
            !
            sig_out(1:m) = wrk_sig(1:m)
            m = m * 2
            s = s + 1
            !
         END DO

      CASE DEFAULT                                             ! error
         STOP 'wlt_transform: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE wlt_transform


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE set_Daubechies ***
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
   !> @details       https://www.di.ens.fr/~mallat/College/WaveletTourChap7.pdf
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
   ! @snippet this set_Daubechies
   ! [set_Daubechies]
   SUBROUTINE set_Daubechies ( order )
      INTEGER(4),  INTENT(in   ) :: order
      !
      INTEGER(4)                 :: ierr

      ALLOCATE( wlt_coef(0:order-1), STAT=ierr)



      SELECT CASE( order )
      !
      !
      !
      CASE( 4 )
         wlt_coef( 0) =  0.482962913145_wp
         wlt_coef( 1) =  0.836516303738_wp
         wlt_coef( 2) =  0.224143868042_wp
         wlt_coef( 3) = -0.129409522551_wp
      !
      !
      !
      CASE( 6 )
         wlt_coef( 0) =  0.332670552950_wp
         wlt_coef( 1) =  0.806891509311_wp
         wlt_coef( 2) =  0.459877502118_wp
         wlt_coef( 3) = -0.135011020010_wp
         wlt_coef( 4) = -0.085441273882_wp
         wlt_coef( 5) =  0.035226291882_wp      
      !
      !
      !
      CASE( 8 )
         wlt_coef( 0) =  0.230377813309_wp 
         wlt_coef( 1) =  0.714846570553_wp
         wlt_coef( 2) =  0.630880767930_wp
         wlt_coef( 3) = -0.027983769417_wp
         wlt_coef( 4) = -0.187034811719_wp
         wlt_coef( 5) =  0.030841381836_wp 
         wlt_coef( 6) =  0.032883011667_wp
         wlt_coef( 7) = -0.010597401785_wp
      !
      !
      !
      CASE( 10 )
         wlt_coef( 0) =  0.160102397974_wp
         wlt_coef( 1) =  0.603829269797_wp
         wlt_coef( 2) =  0.724308528438_wp
         wlt_coef( 3) =  0.138428145901_wp
         wlt_coef( 4) = -0.242294887066_wp
         wlt_coef( 5) = -0.032244869585_wp
         wlt_coef( 6) =  0.077571493840_wp
         wlt_coef( 7) = -0.006241490213_wp
         wlt_coef( 8) = -0.012580751999_wp
         wlt_coef( 9) =  0.003335725285_wp
      !
      !
      !
      CASE( 12 )
         wlt_coef( 0) =  0.111540743350_wp
         wlt_coef( 1) =  0.494623890398_wp
         wlt_coef( 2) =  0.751133908021_wp
         wlt_coef( 3) =  0.315250351709_wp
         wlt_coef( 4) = -0.226264693965_wp
         wlt_coef( 5) = -0.129766867567_wp
         wlt_coef( 6) =  0.097501605587_wp
         wlt_coef( 7) =  0.027522865530_wp
         wlt_coef( 8) = -0.031582039317_wp
         wlt_coef( 9) =  0.000553842201_wp
         wlt_coef(10) =  0.004777257511_wp
         wlt_coef(11) = -0.001077301085_wp
      !
      !
      !
      CASE( 14 )
         wlt_coef( 0) =  0.077852054085_wp
         wlt_coef( 1) =  0.396539319482_wp
         wlt_coef( 2) =  0.729132090846_wp
         wlt_coef( 3) =  0.469782287405_wp
         wlt_coef( 4) = -0.143906003929_wp
         wlt_coef( 5) = -0.224036184994_wp
         wlt_coef( 6) =  0.071309219267_wp
         wlt_coef( 7) =  0.080612609151_wp
         wlt_coef( 8) = -0.038029936935_wp
         wlt_coef( 9) = -0.016574541631_wp
         wlt_coef(10) =  0.012550998556_wp
         wlt_coef(11) =  0.000429577973_wp
         wlt_coef(12) = -0.001801640704_wp
         wlt_coef(13) =  0.000353713800_wp
      !
      !
      !
      CASE( 16 )
         wlt_coef( 0) =  0.054415842243_wp
         wlt_coef( 1) =  0.312871590914_wp
         wlt_coef( 2) =  0.675630736297_wp
         wlt_coef( 3) =  0.585354683654_wp
         wlt_coef( 4) = -0.015829105256_wp
         wlt_coef( 5) = -0.284015542962_wp
         wlt_coef( 6) =  0.000472484574_wp
         wlt_coef( 7) =  0.128747426620_wp
         wlt_coef( 8) = -0.017369301002_wp
         wlt_coef( 9) = -0.044088253930_wp
         wlt_coef(10) =  0.013981027917_wp
         wlt_coef(11) =  0.008746094047_wp
         wlt_coef(12) = -0.004870352993_wp
         wlt_coef(13) = -0.000391740373_wp
         wlt_coef(14) =  0.000675449406_wp
         wlt_coef(15) = -0.000117476784_wp
      !
      !
      !
      CASE( 18 )
         wlt_coef( 0) =  0.038077947364_wp
         wlt_coef( 1) =  0.243834674613_wp
         wlt_coef( 2) =  0.604823123690_wp
         wlt_coef( 3) =  0.657288078051_wp
         wlt_coef( 4) =  0.133197385825_wp
         wlt_coef( 5) = -0.293273783279_wp
         wlt_coef( 6) = -0.096840783223_wp
         wlt_coef( 7) =  0.148540749338_wp
         wlt_coef( 8) =  0.030725681479_wp
         wlt_coef( 9) = -0.067632829061_wp
         wlt_coef(10) =  0.000250947115_wp
         wlt_coef(11) =  0.022361662124_wp
         wlt_coef(12) = -0.004723204758_wp
         wlt_coef(13) = -0.004281503682_wp
         wlt_coef(14) =  0.001847646883_wp
         wlt_coef(15) =  0.000230385764_wp
         wlt_coef(16) = -0.000251963189_wp
         wlt_coef(17) =  0.000039347320_wp
      !
      !
      !
      CASE( 20 )
         wlt_coef( 0) =  0.026670057901_wp
         wlt_coef( 1) =  0.188176800078_wp
         wlt_coef( 2) =  0.527201188932_wp
         wlt_coef( 3) =  0.688459039454_wp
         wlt_coef( 4) =  0.281172343661_wp
         wlt_coef( 5) = -0.249846424327_wp
         wlt_coef( 6) = -0.195946274377_wp
         wlt_coef( 7) =  0.127369340336_wp
         wlt_coef( 8) =  0.093057364604_wp
         wlt_coef( 9) = -0.071394147166_wp
         wlt_coef(10) = -0.029457536822_wp
         wlt_coef(11) =  0.033212674059_wp
         wlt_coef(12) =  0.003606553567_wp
         wlt_coef(13) = -0.010733175483_wp
         wlt_coef(14) =  0.001395351747_wp
         wlt_coef(15) =  0.001992405295_wp
         wlt_coef(16) = -0.000685856695_wp
         wlt_coef(17) = -0.000116466855_wp
         wlt_coef(18) =  0.000093588670_wp
         wlt_coef(19) =  0.000013264203_wp
      !
      !
      !
      CASE DEFAULT                                             ! error
         STOP 'wlt_transform: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE set_Daubechies




END MODULE tluwlt


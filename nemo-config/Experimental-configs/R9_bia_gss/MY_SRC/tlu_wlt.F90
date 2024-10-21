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
   USE tluflt   
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
   INTEGER(i4)                                              ::   nn_tlu_psiz =  7   !> @public width of box filter
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_nlev        !> number of pseudo-observations
   INTEGER(i4), PUBLIC                                      ::   nn_daub_ord        !> Order of Daubechies Wavelet
   ! [tlu_ken_modes]   
   !
   ! [ ]
   INTEGER(i4), PUBLIC                                      :: nn_tlu_hsiz = 2 !> @public size of halo
   INTEGER(i4), PUBLIC                                      :: nn_tlu_hbox = 2 !> @public 
   ! [ ]   

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
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:)     :: wlt_coef               !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: msk_fld                !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   !
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk1, wrk2
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk3, wrk4
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: u_int, v_int
   !
   TYPE(filter)                                             :: gaussian
   TYPE(filter)                                             :: box_3by3
   !
CONTAINS

   SUBROUTINE tlu_init_wlt
      USE tluMPI
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
      INTEGER(i4) :: ios, ierr(9), chk   ! namelist output, allocation statistics
      INTEGER     :: i, ki, kj
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
      ! Allocate interior fields
      !
      ALLOCATE(   u_int( jpi-2, jpj-2, jpk ), &
      &           v_int( jpi-2, jpj-2, jpk ), &
      &            wrk1( jpi-2, jpj-2, jpk ), &
      &            wrk2( jpi-2, jpj-2, jpk ), &
      &            wrk3( jpi-2, jpj-2, jpk ), &
      &            wrk4( jpi-2, jpj-2, jpk ), stat=ierr(2) )
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( rnd_fld( jpi-2, jpj-2 ), stat=ierr(3) )
      !
      rnd_fld = 0._wp 
      !
      ! Allocate mask
      !
      ALLOCATE( msk_fld( jpim1, jpjm1 ), stat=ierr(4) )
      !
      msk_fld = 1._wp
      !
      ! Allocate and set Wavelet coefficients
      !
      CALL set_Daubechies ( nn_daub_ord )

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
      USE mpi
      USE tlurnd
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
      INTEGER ::   jm, m_idx, ierr, ierror ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at
      INTEGER ::   kt
      INTEGER ::   ncid 
      INTEGER(i4) :: i,j,k, ji, jj, jk, jl, ki, kj 
      INTEGER(i4) :: mx, my, ni_2, nj_2, mi, mj
 
      INTEGER(i4)                                              :: fim, fjm          !> @private first index after extra halo
      INTEGER(i4)                                              :: lim, ljm          !> @private last index begore extra halo
      INTEGER ::   rel_idx(nn_tlu_psiz)
      INTEGER :: my_rank
      INTEGER :: status(MPI_STATUS_SIZE)

      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'gen_wlt_mod : TLU, Wavelet-based generation of modes '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      u_int = 0._wp
      v_int = 0._wp
      !
      wrk1 = 0._wp
      wrk2 = 0._wp
      wrk3 = 0._wp
      wrk4 = 0._wp
      !
      wlx_noi = 0._wp
      wly_noi = 0._wp
      wlz_noi = 0._wp
      !
      ! Reshape field
      !
      u_int = un( 2:jpi-1, 2:jpj-1, : )
      v_int = vn( 2:jpi-1, 2:jpj-1, : )
      !
      ! Apply 2D direct wavelet transform
      !
      DO jk = 1, jpk
         !
         !
         !
         CALL wavelet_2D_ext( u_int(:, :, jk ), wrk1(:, :, jk ), nn_tlu_nlev, wlt_coef, nn_daub_ord, "dir" )
         CALL wavelet_2D_ext( v_int(:, :, jk ), wrk2(:, :, jk ), nn_tlu_nlev, wlt_coef, nn_daub_ord, "dir" )
         !
      END DO
      !
      ! Gaussian field generation
      !
      CALL Standard_Gaussian_2D( jpi-2, jpj-2, rnd_fld ) !, allones = .False. )
      !rnd_fld = 1._wp
      rnd_fld( 1 : (jpi-2)/(2**nn_tlu_nlev) + 1, 1 : (jpj-2)/(2**nn_tlu_nlev) + 1 ) = 0._wp
      rnd_fld = 0.75_wp * rnd_fld
      !
      ! Inverse Wavelet Transform
      !
      DO jk = 1, jpk
         !
         !
         wrk1(:,:,jk) = wrk1(:,:,jk) * rnd_fld
         wrk2(:,:,jk) = wrk2(:,:,jk) * rnd_fld
         !
         CALL wavelet_2D_ext( wrk1(:, :, jk ), wrk3(:, :, jk ), nn_tlu_nlev, wlt_coef, nn_daub_ord, "inv" )
         CALL wavelet_2D_ext( wrk2(:, :, jk ), wrk4(:, :, jk ), nn_tlu_nlev, wlt_coef, nn_daub_ord, "inv" )
         !
      END DO  
      !
      wlx_noi( 2:jpi-1, 2:jpj-1, : ) = wrk3
      wly_noi( 2:jpi-1, 2:jpj-1, : ) = wrk4
      !
      ! Transfer horizontal fields across nodes
      !
      CALL lbc_lnk_multi( 'gen_wlt_mod', wlx_noi , 'U', -1., wly_noi , 'V', -1. )
      !
      ! Compute vertical velocity
      !
      CALL build_vmodes(wlz_noi, wlx_noi, wly_noi)
      !
      ! Transfer vertical fields across nodes
      !
      CALL lbc_lnk_multi( 'gen_wlt_mod', wlz_noi , 'T', -1. )      
      !
   END SUBROUTINE gen_wlt_mod
   ! [gen_wlt_mod]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE wlt_2D_transform ***
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
   ! @snippet this wavelet_2D_ext
   ! [wavelet_2D_ext]
   SUBROUTINE wavelet_2D_ext(sig_in, sig_out, ndl, wlt_coef, d, direction)
      USE tluMPI
      CHARACTER(3),                            INTENT(in   ) :: direction
      INTEGER(4),                              INTENT(in   ) :: d
      INTEGER(4),                              INTENT(in   ) :: ndl
      REAL(wp),              DIMENSION(0:d-1), INTENT(in   ) :: wlt_coef
      REAL(wp),                DIMENSION(:,:), INTENT(in   ) :: sig_in
      REAL(wp),                DIMENSION(:,:), INTENT(  out) :: sig_out
      !
      INTEGER(4)                                             :: ierr
      INTEGER(4)                                             :: i0, i1, j0, j1
      INTEGER(4)                                             ::  i,  j, k, m, n, p, q, s, ji, jj
      INTEGER(i4)                                            :: ni, nj, mi, mj !, k, m, n, p, q, s
      INTEGER(i4)                                            :: is, js, ie, je     
      INTEGER(i4)                                            :: dm2 
      
      REAL(wp),                DIMENSION(d-2,jpj-2) :: sig_lr
      REAL(wp),                DIMENSION(jpi-2,d-2) :: sig_tb


      Real(wp), dimension(jpi-2,jpj-2) :: dummy
      REAL(wp), DIMENSION(  d-2,jpj-2) :: dummy_lr
      REAL(wp), DIMENSION(jpi-2,  d-2) :: dummy_tb
      INTEGER ::  ierror, my_rank
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

      !
      ni = SIZE(sig_in, 1)
      nj = SIZE(sig_in, 2)
      !
      s = 0   ! the relative scale counter
      !
      p = d - 1
      dm2 = d - 2
      q = ( p - 1 ) / 2
      !
      sig_out = sig_in
      !
      sig_lr = 0._wp
      sig_tb = 0._wp

      do jj = 1, jpj-2
         dummy(:,jj) = mig(2) + jj - 2
      end do

      !
      SELECT CASE( direction )
      !
      ! Direct transform
      !
      CASE( "dir" )
         !
         mi = ni ! not working on halo
         mj = nj ! not working on halo
         !
         DO WHILE ( s < ndl )   ! for each decomposition level
            sig_lr = 0._wp
            sig_tb = 0._wp        
            !
            ! Collect the data from top 2 bottom
            !
            is = 1 
            ie = jpi - 2
            js = 1
            je = js + dm2 - 1
            !
            CALL  MPI_collect_t2b( 1, sig_out, jpi-2, jpj-2, sig_tb, is, js, ie, je, .true.)
            CALL     MPI_wrap_b2t( 1, sig_out, jpi-2, jpj-2, sig_tb, is, js, ie, je, .true.)
            !
            DO ji = 1, mi
               !
              CALL wavelet_sl_ext(sig_out(ji, 1:mj), sig_tb(ji, :), sig_out(ji, 1:mj), mj, wlt_coef, d, "dir")
               !
            END DO           
            !
            ! Collect the data from right 2 left
            !
            is = 1 
            ie = is + dm2 - 1
            js = 1
            je = jpj - 2
            !
            CALL MPI_collect_r2l( 1, sig_out, jpi-2, jpj-2,   sig_lr, is, js, ie, je, .true.)
            CALL    MPI_wrap_l2r( 1, sig_out, jpi-2, jpj-2,   sig_lr, is, js, ie, je, .true.)
            !
            ! Perform the wavelet by rows
            !
            DO jj = 1, mj
               !
               CALL wavelet_sl_ext(sig_out(1:mi, jj ), sig_lr(:, jj ), sig_out(1:mi, jj ), mi, wlt_coef, d, "dir")
               !
            END DO 
            !
            mi = mi / 2
            mj = mj / 2
            s = s + 1
            !
         END DO
      !
      ! Inverse transform
      !
      CASE( "inv" )
         mi = 2 * rshift(ni, ndl) !ni
         mj = 2 * rshift(nj, ndl) !nj
         !
         DO WHILE ( s < ndl )   ! for each decomposition level
            sig_lr = 0._wp
            sig_tb = 0._wp         
            !
            ! Collect the data from bottom 2 top
            !
            is = 1
            ie = jpi - 2
            js = mj/2 - dm2/2 + 1
            je = mj/2
            CALL MPI_collect_b2t( 1, sig_out, jpi-2, jpj-2, sig_tb(:,1:dm2/2), is, js, ie, je, .true.)
            CALL    MPI_wrap_t2b( 1, sig_out, jpi-2, jpj-2, sig_tb(:,1:dm2/2), is, js, ie, je, .true.)
            !
            is = 1 
            ie = jpi - 2
            js = mj - dm2/2 + 1
            je = mj
            !
            CALL MPI_collect_b2t( 1, sig_out, jpi-2, jpj-2, sig_tb(:,1+dm2/2:dm2), is, js, ie, je, .true.)
            CALL    MPI_wrap_t2b( 1, sig_out, jpi-2, jpj-2, sig_tb(:,1+dm2/2:dm2), is, js, ie, je, .true.)
            !
            DO ji = 1, mi
               !
               CALL wavelet_sl_ext(sig_out(ji, 1:mj), sig_tb(ji, :), sig_out(ji, 1:mj), mj, wlt_coef, d, "inv")
               !
            END DO   
            !
            ! Collect the data from left 2 right
            !
            is = mi/2 - dm2/2 + 1 
            ie = mi/2
            js = 1
            je = jpj - 2
            !
            CALL MPI_collect_l2r( 11, sig_out, jpi-2, jpj-2, sig_lr( 1      :dm2/2,: ), is, js, ie, je, .true.)
            CALL    MPI_wrap_r2l( 11, sig_out, jpi-2, jpj-2, sig_lr( 1      :dm2/2,: ), is, js, ie, je, .true.)
            !
            is = mi - dm2/2 + 1 
            ie = mi
            js = 1
            je = jpj - 2
            !
            CALL MPI_collect_l2r( 11, sig_out, jpi-2, jpj-2, sig_lr( 1+dm2/2:dm2  ,: ), is, js, ie, je, .true.)
            CALL    MPI_wrap_r2l( 11, sig_out, jpi-2, jpj-2, sig_lr( 1+dm2/2:dm2  ,: ), is, js, ie, je, .true.)
            !
            ! Perform the wavelet by rows
            !
            DO jj = 1, mj
               ! 
               CALL wavelet_sl_ext(sig_out(1:mi, jj ), sig_lr(:, jj ), sig_out(1:mi, jj ), mi, wlt_coef, d, "inv")
               !
            END DO 
            !
            mi = mi * 2
            mj = mj * 2
            s = s + 1
            !
         END DO
         
      CASE DEFAULT                                             ! error
         STOP 'wavelet_2D_ext: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE wavelet_2D_ext 


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE sl_ext_wavelet ***
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
   ! @snippet this sl_ext_wavelet
   ! [sl_ext_wavelet]
   SUBROUTINE wavelet_sl_ext(sig_in, sig_ex, sig_out, n, wlt_coef, d, direction)
      CHARACTER(3),                   INTENT(in   ) :: direction
      INTEGER(4),                     INTENT(in   ) :: n
      INTEGER(4),                     INTENT(in   ) :: d
      REAL(wp),     DIMENSION(0:d-1), INTENT(in   ) :: wlt_coef
      REAL(wp),     DIMENSION( n   ), INTENT(in   ) :: sig_in
      REAL(wp),     DIMENSION( d-2 ), INTENT(in   ) :: sig_ex
      REAL(wp),     DIMENSION( n   ), INTENT(  out) :: sig_out
      !
      REAL(wp),     DIMENSION( n+d )                :: d_sig_in
      REAL(wp),     DIMENSION(1-(d-2)/2: n+(d-2)/2 ):: i_sig_in
      REAL(wp),     DIMENSION( n   )                :: wrk_sig
      INTEGER(4)                                    :: i0, i1
      INTEGER(4)                                    :: j0, j1
      INTEGER(4)                                    :: i, j, k, m, p, q
      INTEGER(i4)                                   :: dm2 
      !
      ! Wavelet parameters
      !
      m = n
      p = d - 1
      dm2 = d - 2
      q = ( p - 1 ) / 2
      !
      SELECT CASE( direction )
      !
      ! Direct transform
      !
      CASE( "dir" )
         !
         i = 1
         wrk_sig(1:m) = 0._wp
         ! Extending the work signal with extra data bypass wrap-around needs
         d_sig_in(    1:m          ) = sig_in(1 : m  )
         d_sig_in(m + 1:m + dm2 - 2) = sig_ex(1 : dm2)
         !
         DO j = 1, m - 1, 2
            !
            DO k = 0, p - 1, 2
               j0 = j + k        ! instead of i4_wrap ( j + k,     1, m )
               j1 = j + k + 1    ! instead of i4_wrap ( j + k +1,  1, m )
               wrk_sig(i)     = wrk_sig(i)     + wlt_coef(  k) * d_sig_in(j0) + wlt_coef(  k+1) * d_sig_in(j1)
               wrk_sig(i+m/2) = wrk_sig(i+m/2) + wlt_coef(p-k) * d_sig_in(j0) - wlt_coef(p-k-1) * d_sig_in(j1)
            END DO
            !
            i = i + 1
            !
         END DO
         !
         sig_out(1:m) = wrk_sig(1:m)
      !
      ! Inverse transform
      !
      CASE( "inv" )
         !
         j = 1
         wrk_sig(1:m) = 0._wp
         ! Extending the work signal with extra data bypass wrap-around needs
         i_sig_in(1         - q : 0         ) = sig_ex(1         : q     )
         i_sig_in(1             : m / 2     ) = sig_in(1         : m / 2 )
         i_sig_in(1 + m / 2     : m / 2 + q ) = sig_ex(1 + q     : dm2   )
         i_sig_in(1 + m / 2 + q : m + q     ) = sig_in(1 + m / 2 : m     )
         !
         do i = - q + 1, m / 2 - q
            !
            do k = 0, p - 1, 2
               ! here the q is necessary to shift the index from original (m/2:m) to (q+m/2:q+m)
               i0 =     i         + k / 2 ! instead of i4_wrap ( i         + k / 2,     1,         m / 2 )
               i1 = q + i + m / 2 + k / 2 ! instead of i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
               wrk_sig(j)   = wrk_sig(j)   + wlt_coef(p-k-1) * i_sig_in(i0) + wlt_coef(k+1) * i_sig_in(i1)
               wrk_sig(j+1) = wrk_sig(j+1) + wlt_coef(p-k)   * i_sig_in(i0) - wlt_coef(k)   * i_sig_in(i1)
            end do
            !
            j = j + 2
            !
         end do
         !
        ! stop
         sig_out(1:m) = wrk_sig(1:m)      
         !
      CASE DEFAULT                                             ! error
         STOP 'wavelet_sl_ext: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE wavelet_sl_ext


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE MPI_wrap_l2r ***
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
   ! @snippet this MPI_wrap_l2r
   ! [MPI_wrap_l2r]
   SUBROUTINE MPI_wrap_l2r(kt, field_in, ni, nj, field_out, is, js, ie, je, hold)
  !
   INCLUDE 'mpif.h'
      INTEGER(i4),                              INTENT(in   ) :: kt
      INTEGER(i4),                              INTENT(in   ) :: ni, nj
      REAL(wp),   DIMENSION( ni, nj),           INTENT(in   ) :: field_in
      LOGICAL,                                  INTENT(in   ) :: hold
      INTEGER(i4),                              INTENT(in   ) :: is
      INTEGER(i4),                              INTENT(in   ) :: js
      INTEGER(i4),                              INTENT(in   ) :: ie
      INTEGER(i4),                              INTENT(in   ) :: je
      REAL(wp),   DIMENSION( ie-is+1, je-js+1), INTENT(  out) :: field_out
      REAL(wp),   DIMENSION( ie-is+1, je-js+1)                :: lint, rext 
      !
      ! MPI variables
      INTEGER(i4) :: ierror, status, tag, my_rank      
      !
      lint = 0._wp
      rext = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      lint = field_in( is : ie, js : je)
      tag = 65464
      !
      !
      !
      IF (jpni .ne. 1) THEN
         !
         ! Get my rank and do the corresponding job
         CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
         !
         ! Send from West to East
         !
         IF ( MOD(my_rank, jpni) .eq.      0 ) THEN
            !print *, INT(my_rank/jpni),my_rank, "-->", my_rank+jpni-1 
            CALL MPI_send(lint, size(lint), MPI_double_precision, my_rank+jpni-1,         tag, MPI_COMM_WORLD, ierror)
         END IF

         IF ( MOD(my_rank, jpni) .eq. jpni-1 ) THEN
            !print *, INT(my_rank/jpni),my_rank, "<--", my_rank-jpni+1
            field_out = 0._wp
            CALL MPI_Recv(rext, size(rext), MPI_double_precision, my_rank-jpni+1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
            field_out = rext
         END IF
         !
         ! Safely wait
         !
         IF ( hold ) call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      ELSE
         field_out = lint
      END IF
      !
   END SUBROUTINE MPI_wrap_l2r

   SUBROUTINE MPI_wrap_r2l( kt, field_in, ni, nj, field_out, is, js, ie, je, hold)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE MPI_wrap_r2l  ***
      !!
      !! ** Purpose :  Collect information from LEFT to RIGHT (l2r) 
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!-----------------------------------------------------------------------
      INTEGER(i4),                              INTENT(in   ) :: kt
      INTEGER(i4),                              INTENT(in   ) :: ni, nj
      REAL(wp),   DIMENSION( ni, nj),           INTENT(in   ) :: field_in
      LOGICAL,                                  INTENT(in   ) :: hold
      INTEGER(i4),                              INTENT(in   ) :: is
      INTEGER(i4),                              INTENT(in   ) :: js
      INTEGER(i4),                              INTENT(in   ) :: ie
      INTEGER(i4),                              INTENT(in   ) :: je
      REAL(wp),   DIMENSION( ie-is+1, je-js+1), INTENT(  out) :: field_out
      REAL(wp),   DIMENSION( ie-is+1, je-js+1)                :: rint, lext 
      !
      ! MPI variables
      INTEGER(i4) :: ierror, status, tag, my_rank      
      !
      rint = 0._wp
      lext = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      rint = field_in( is : ie, js : je )
      tag = 65464
      !
      !
      !
      IF (jpni .ne. 1) THEN      
         !
         ! Get my rank and do the corresponding job
         CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
         !
         ! West to East
         !
         IF ( MOD(my_rank, jpni) .eq. jpni-1 ) THEN
            !print *, MOD(my_rank, jpni),my_rank, "-->", my_rank-jpni+1 
            CALL MPI_send(rint, size(rint), MPI_double_precision, my_rank-jpni+1,         tag, MPI_COMM_WORLD, ierror)
         END IF

         IF ( MOD(my_rank, jpni) .eq.      0 ) THEN
            !print *, MOD(my_rank, jpni),my_rank, "<--", my_rank+jpni-1
            field_out = 0._wp
            CALL MPI_Recv(lext, size(lext), MPI_double_precision, my_rank+jpni-1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
            field_out = lext
         END IF
         !
         ! Safely wait
         !
         IF ( hold ) call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      ELSE
         field_out = rint
      END IF
      !
      END SUBROUTINE MPI_wrap_r2l


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE MPI_wrap_t2b ***
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
   ! @snippet this MPI_wrap_b2t
   ! [MPI_wrap_b2t]
   SUBROUTINE MPI_wrap_b2t(kt, field_in, ni, nj, field_out, is, js, ie, je, hold)
  !
   INCLUDE 'mpif.h'
      INTEGER(i4),                              INTENT(in   ) :: kt
      INTEGER(i4),                              INTENT(in   ) :: ni, nj
      REAL(wp),   DIMENSION( ni, nj),           INTENT(in   ) :: field_in
      LOGICAL,                                  INTENT(in   ) :: hold
      INTEGER(i4),                              INTENT(in   ) :: is
      INTEGER(i4),                              INTENT(in   ) :: js
      INTEGER(i4),                              INTENT(in   ) :: ie
      INTEGER(i4),                              INTENT(in   ) :: je
      REAL(wp),   DIMENSION( ie-is+1, je-js+1), INTENT(  out) :: field_out
      REAL(wp),   DIMENSION( ie-is+1, je-js+1)                :: bint, text 
      !
      ! MPI variables
      INTEGER(i4) :: ierror, status, tag, my_rank      
      !
      bint = 0._wp
      text = 0._wp
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      bint = field_in( is : ie, js : je)
      tag = 65464
      !
      !
      !
      IF (jpnj .ne. 1) THEN
         !
         ! Get my rank and do the corresponding job
         CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
         !
         ! Send from West to East
         !
         IF ( INT(my_rank/jpni)     .eq.    0 ) THEN
            !print *, INT(my_rank/jpni),my_rank, "-->", my_rank+jpni*(jpnj-1)
            CALL MPI_send(bint, size(bint), MPI_double_precision, my_rank+jpni*(jpnj-1),         tag, MPI_COMM_WORLD, ierror)
         END IF

         IF ( INT(my_rank/jpni) +1  .eq. jpnj ) THEN
            !print *, INT(my_rank/jpni),my_rank, "<--", my_rank-jpni*(jpnj-1)
            field_out = 0._wp
            CALL MPI_Recv(text, size(text), MPI_double_precision, my_rank-jpni*(jpnj-1), MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
            field_out = text
         END IF
         !
         ! Safely wait
         !
         IF ( hold ) call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      ELSE
         field_out = bint
      END IF
      !
   END SUBROUTINE MPI_wrap_b2t 

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE MPI_wrap_t2b ***
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
   ! @snippet this MPI_wrap_t2b
   ! [MPI_wrap_t2b]
   SUBROUTINE MPI_wrap_t2b(kt, field_in, ni, nj, field_out, is, js, ie, je, hold)
  !
   INCLUDE 'mpif.h'
      INTEGER(i4),                              INTENT(in   ) :: kt
      INTEGER(i4),                              INTENT(in   ) :: ni, nj
      REAL(wp),   DIMENSION( ni, nj),           INTENT(in   ) :: field_in
      LOGICAL,                                  INTENT(in   ) :: hold
      INTEGER(i4),                              INTENT(in   ) :: is
      INTEGER(i4),                              INTENT(in   ) :: js
      INTEGER(i4),                              INTENT(in   ) :: ie
      INTEGER(i4),                              INTENT(in   ) :: je
      REAL(wp),   DIMENSION( ie-is+1, je-js+1), INTENT(  out) :: field_out
      REAL(wp),   DIMENSION( ie-is+1, je-js+1)                :: tint, bext 
      !
      ! MPI variables
      INTEGER :: ierror, status, tag, my_rank
      !
      tint = 0._wp
      bext = 0._wp
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      tint = field_in( is : ie, js : je)
      tag = 65464
      !
      !
      !
      IF (jpnj .ne. 1) THEN
         !
         ! Get my rank and do the corresponding job
         CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
         !
         ! Send from West to East
         !
         IF ( INT(my_rank/jpni) +1  .eq. jpnj ) THEN         
            !print *, INT(my_rank/jpni),my_rank, "-->", my_rank-jpni*(jpnj-1)
            CALL MPI_send(tint, size(tint), MPI_double_precision, my_rank-jpni*(jpnj-1),         tag, MPI_COMM_WORLD, ierror)
         END IF

         IF ( INT(my_rank/jpni)     .eq.    0 ) THEN
            !print *, INT(my_rank/jpni),my_rank, "<--", my_rank+jpni*(jpnj-1)
            field_out = 0._wp
            CALL MPI_Recv(bext, size(bext), MPI_double_precision, my_rank+jpni*(jpnj-1), MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
            field_out = bext
         END IF
         !
         ! Safely wait
         !
         IF ( hold ) call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      ELSE
         field_out = tint
      END IF
      !
   END SUBROUTINE MPI_wrap_t2b


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE wavelet_2D_wrp ***
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
   ! @snippet this wavelet_2D_wrp
   ! [wavelet_2D_wrp]
   SUBROUTINE wavelet_2D_wrp(sig_in, sig_out, ndl, wlt_coef, d, direction)
      USE tluMPI
      CHARACTER(3),                            INTENT(in   ) :: direction
      INTEGER(4),                              INTENT(in   ) :: d
      INTEGER(4),                              INTENT(in   ) :: ndl
      REAL(wp),              DIMENSION(0:d-1), INTENT(in   ) :: wlt_coef
      REAL(wp),                DIMENSION(:,:), INTENT(in   ) :: sig_in
      REAL(wp),                DIMENSION(:,:), INTENT(  out) :: sig_out
      !
      INTEGER(4)                                             :: ierr
      INTEGER(4)                                             :: i0, i1, j0, j1
      INTEGER(4)                                             ::  i,  j, k, m, n, p, q, s, ji, jj
      INTEGER(i4)                                            :: ni, nj, mi, mj !, k, m, n, p, q, s
      INTEGER(i4)                                            :: is, js, ie, je     
      INTEGER(i4)                                            :: dm2 
      
      REAL(wp),                DIMENSION(d-2,jpj-2) :: sig_lr
      REAL(wp),                DIMENSION(jpi-2,d-2) :: sig_tb
      !
      ni = SIZE(sig_in, 1)
      nj = SIZE(sig_in, 2)
      !
      s = 0   ! the relative scale counter
      !
      p = d - 1
      dm2 = d - 2
      q = ( p - 1 ) / 2
      !
      sig_out = sig_in
      !
      sig_lr = 0._wp
      sig_tb = 0._wp
      !
      SELECT CASE( direction )
      !
      ! Direct transform
      !
      CASE( "dir" )
         !
         mi = ni
         mj = nj
         !
         DO WHILE ( s < ndl )   ! for each decomposition level
            !
            ! Other way around         
            !
            DO ji = 1, mi
               !
               CALL wavelet_sl_wrp(  sig_out(ji, 1:mj ), sig_lr, sig_out(ji, 1:mj ), mj, wlt_coef, d, "dir")
               !
            END DO
            !
            ! One way 
            !
            DO jj = 1, mj
               !
               CALL wavelet_sl_wrp(  sig_out(1:mi, jj ), sig_lr, sig_out(1:mi, jj ), mi, wlt_coef, d, "dir")
               !
            END DO
            !
            mi = mi / 2
            mj = mj / 2
            s = s + 1
            !
         END DO
      !
      ! Inverse transform
      !
      CASE( "inv" )
         mi = 2 * rshift(ni, ndl) !ni
         mj = 2 * rshift(nj, ndl) !nj
         !
         sig_out = sig_in
         !
         DO WHILE ( s < ndl )   ! for each decomposition level

            DO ji = 1, mi
               !
               CALL wavelet_sl_wrp(  sig_out(ji, 1:mj ), sig_lr, sig_out(ji, 1:mj ), mj, wlt_coef, d, "inv")
               !
            END DO
            !
            ! One way 
            !
            DO jj = 1, mj
               !
               CALL wavelet_sl_wrp( sig_out(1:mi, jj ), sig_lr, sig_out(1:mi, jj ), mi, wlt_coef, d, "inv")
               !
            END DO
            !
            mi = mi * 2
            mj = mj * 2
            s = s + 1
            !
         END DO
         
      CASE DEFAULT                                             ! error
         STOP 'wavelet_2D_wrp: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE wavelet_2D_wrp 


   SUBROUTINE wavelet_sl_wrp(sig_in, sig_ex, sig_out, n, wlt_coef, d, direction)
      CHARACTER(3),                   INTENT(in   ) :: direction
      INTEGER(4),                     INTENT(in   ) :: n
      INTEGER(4),                     INTENT(in   ) :: d
      REAL(wp),     DIMENSION(0:d-1), INTENT(in   ) :: wlt_coef
      REAL(wp),     DIMENSION( n   ), INTENT(in   ) :: sig_in
      REAL(wp),     DIMENSION( d-2 ), INTENT(in   ) :: sig_ex
      REAL(wp),     DIMENSION( n   ), INTENT(  out) :: sig_out
      !
      REAL(wp),     DIMENSION( n+d )                :: d_wrk_sig
      REAL(wp),     DIMENSION(1-(d-2)/2: n+(d-2)/2 ):: i_wrk_sig
      REAL(wp),     DIMENSION( n   )                :: wrk_sig
      INTEGER(4)                                    :: i0, i1
      INTEGER(4)                                    :: j0, j1
      INTEGER(4)                                    :: i, j, k, m, p, q
      INTEGER(i4)                                   :: dm2 

      !
      ! Wavelet parameters
      !
      m = n
      p = d - 1
      dm2 = d - 2
      q = ( p - 1 ) / 2
      !
      SELECT CASE( direction )
      !
      ! Direct transform
      !
      CASE( "dir" )
         !
         i = 1
         wrk_sig(1:m) = 0._wp
         !
         DO j = 1, m - 1, 2
            !
            DO k = 0, p - 1, 2
               j0 = i4_wrap ( j + k,     1, m )
               j1 = i4_wrap ( j + k + 1, 1, m )
               wrk_sig(i)     = wrk_sig(i)     + wlt_coef(  k) * sig_in(j0) + wlt_coef(  k+1) * sig_in(j1)
               wrk_sig(i+m/2) = wrk_sig(i+m/2) + wlt_coef(p-k) * sig_in(j0) - wlt_coef(p-k-1) * sig_in(j1)
            END DO
            !
            i = i + 1
            !
            END DO
         !
         sig_out(1:m) = wrk_sig(1:m)
      !
      ! Inverse transform
      !
      CASE( "inv" )
         !
         j = 1
         wrk_sig(1:m) = 0._wp
         !
         do i = - q + 1, m / 2 - q
            !
            do k = 0, p - 1, 2
               i0 = i4_wrap ( i         + k / 2,     1,         m / 2 )
               i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2 + 1, m     )
               wrk_sig(j)   = wrk_sig(j)   + wlt_coef(p-k-1) * sig_in(i0) + wlt_coef(k+1) * sig_in(i1)
               wrk_sig(j+1) = wrk_sig(j+1) + wlt_coef(p-k)   * sig_in(i0) - wlt_coef(k)   * sig_in(i1)
               !IF (lwp) print *,j, "i0=", i0, "i1=", i1, "s(i0)", sig_in(i0), "s(i1)", sig_in(i1)
            end do
            !
            j = j + 2
            !
         end do
         !
         sig_out(1:m) = wrk_sig(1:m)
         !
      CASE DEFAULT                                             ! error
         STOP 'wavelet_sl_wrp: wrong direction for transform, must be dir or inv'
      END SELECT
      !
   END SUBROUTINE wavelet_sl_wrp


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
      CALL lbc_lnk_multi( 'build_vmodes', wout , 'T', -1.  )
      !
   END SUBROUTINE build_vmodes   
   ! [tlu_wzvcmp]

  !!===========================================================================
  !!                LOW-LEVEL
  !!===========================================================================

  function i4_modp ( i, j )

  !*****************************************************************************80
  !
  !! I4_MODP returns the nonnegative remainder of I4 division.
  !
  !  Discussion:
  !
  !    If
  !      NREM = I4_MODP ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
  !
  !    An I4 is an integer ( kind = 4 ) value.
  !
  !  Example:
  !
  !        I     J     MOD I4_MODP    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) I, the number to be divided.
  !
  !    Input, integer ( kind = 4 ) J, the number that divides I.
  !
  !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
  !    divided by J.
  !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) j
    integer ( kind = 4 ) value

    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_MODP - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
      stop
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
      value = value + abs ( j )
    end if

    i4_modp = value

    return
  end

  function i4_wrap ( ival, ilo, ihi )

  !*****************************************************************************80
  !
  !! I4_WRAP forces an I4 to lie between given limits by wrapping.
  !
  !  Discussion:
  !
  !    An I4 is an integer ( kind = 4 ) value.
  !
  !    There appears to be a bug in the GFORTRAN compiler which can lead to
  !    erroneous results when the first argument of I4_WRAP is an expression.
  !    In particular:
  !
  !    do i = 1, 3
  !      if ( test ) then
  !        i4 = i4_wrap ( i + 1, 1, 3 )
  !      end if
  !    end do
  !
  !    was, when I = 3, returning I4 = 3.  So I had to replace this with
  !
  !    do i = 1, 3
  !      if ( test ) then
  !        i4 = i + 1
  !        i4 = i4_wrap ( i4, 1, 3 )
  !      end if
  !    end do
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  Value
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 September 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) IVAL, a value.
  !
  !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
  !
  !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
  !
    implicit none

    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) ival
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    integer ( kind = 4 ) value
    integer ( kind = 4 ) wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
      value = jlo
    else
      value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value

    return
  end

  
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
      CASE( 2 )
         wlt_coef( 0) =  0.7071067811865475_wp
         wlt_coef( 1) =  0.7071067811865475_wp
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
         PRINT *, "received value: ", order
         STOP 'set_Daubechies: wrong dimension for the wavelet, choose upon 4, 6, 8, 10, 12, 14, 16, 18, 20'
      END SELECT
      !
   END SUBROUTINE set_Daubechies

END MODULE tluwlt


MODULE tlusvd
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
   PUBLIC tlu_init_svd    ! called by tlu_init in tlu.f90
   PUBLIC gen_vel_mod     ! called by tlu_fields in tlu_noi.F90
   !
   INCLUDE 'mpif.h'
   !
   ! [tlu_ken_modes]
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_psiz = 3    !: size of patch
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_hsiz = 1    !: size of patch
   INTEGER(i4), PUBLIC                                      ::   nn_tlu_pobs        !: number of pseudo-observations
   ! [tlu_ken_modes]   
   !
   ! [tlu_spmodes] 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spx_mod                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spy_mod                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spz_mod                 !> @public   W-velocity modes
   ! [tlu_spmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),    PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: brwn_rv                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_psob] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: u_psob                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: v_psob                 !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: w_psob                 !> @public   W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: u_mean                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: v_mean                 !> @public   V-velocity pseudo-observations
   ! [tlu_psob]
   !
   
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dAu, dAv, dAt           !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dzu, dzv, dzw           !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dVu, dVv, dVw           !> @public   V-velocity pseudo-observations
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: real_unif_01 
   INTEGER(i4),         ALLOCATABLE,       DIMENSION(:,:)   :: pick_id 
   ! [tlu_randgauss] 
   !
   !
CONTAINS

   SUBROUTINE tlu_init_svd
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_init_svd  ***
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
      NAMELIST/namtlu_svd/   nn_tlu_psobs,  &
                         &   biaSIGN
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
         WRITE(numout,*) '           data-free noise model        ln_tlu_svd = ', ln_tlu_svd
         WRITE(numout,*) '   number of pseudo-observations       nn_tlu_nmod = ', nn_tlu_pobs
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      !
      ALLOCATE( spx_mod(jpi, jpj, nn_tlu_pobs * jpk), &
      &         spy_mod(jpi, jpj, nn_tlu_pobs * jpk), &
      &         spz_mod(jpi, jpj, nn_tlu_pobs * jpk), stat=ierr(1) )   ! the singular values
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( brwn_rv(nn_tlu_pobs), stat=ierr(2) ) 
      !
      ! Allocate Internal variables 
      !
      ALLOCATE( u_psob(jpi * jpj * jpk, nn_tlu_pobs), &
      &         v_psob(jpi * jpj * jpk, nn_tlu_pobs), &
      &         w_psob(jpi * jpj * jpk, nn_tlu_pobs), &
      &         u_mean(jpi * jpj * jpk),              &
      &         v_mean(jpi * jpj * jpk),              stat=ierr(3) )   ! the singular values
      !
      ! Allocate scale factors for online POD
      !
      ALLOCATE(    dAu(jpi * jpj * jpk),              &
      &            dAv(jpi * jpj * jpk),              &
      &            dAt(jpi * jpj * jpk),              &
      &            dzu(jpi * jpj * jpk),              &
      &            dzv(jpi * jpj * jpk),              &
      &            dzw(jpi * jpj * jpk),              &
      &            dVu(jpi * jpj * jpk),              &
      &            dVv(jpi * jpj * jpk),              &
      &            dVw(jpi * jpj * jpk),              stat=ierr(4) )   ! the singular values
      !
      ! Area factors does not change in time, can be set now
      !
      dAu = reshape( spread( e1e2u,3,jpk), jpi * jpj * jpk )
      dAv = reshape( spread( e1e2v,3,jpk), jpi * jpj * jpk )
      dAt = reshape( spread( e1e2t,3,jpk), jpi * jpj * jpk )
      !
      ALLOCATE( real_unif_01(jpi * jpj, nn_tlu_pobs), &
      &              pick_id(jpi * jpj, nn_tlu_pobs), stat=ierr(5) )   ! the singular values
      !
      ALLOCATE( Cmat(nn_tlu_pobs, nn_tlu_pobs),       &
      &         Evec(nn_tlu_pobs, nn_tlu_pobs),       &
      &         Eval(nn_tlu_pobs),                    stat=ierr(6) )
      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_pod_init(): allocation failed = ', ierr  
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !
      !
   END SUBROUTINE tlu_init_svd


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
   SUBROUTINE gen_vel_mod( )
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
         WRITE(numout,*) 'gen_vel_mod : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      ! Update scale factors for online POD
      !
      dAdzu = dAu * reshape( e3u_n, jpi * jpj * jpk )
      dAdzv = dAv * reshape( e3v_n, jpi * jpj * jpk )
      dAdzw = dAt * reshape( e3w_n, jpi * jpj * jpk )

      rdAdz_dVu = dAdzu / SUM( dAdzu )
      rdAdz_dVv = dAdzv / SUM( dAdzv )
      !
      ! Generate Pseudo-Observations
      !
      CALL gen_psobs( )
      !
      ! Build Correlation matrix
      !
      CALL build_Cmat( )
      !
      ! Solve eigenvalue problem
      !
      CALL solve_Cmat( )
      !
      ! Build spatial modes
      !
      CALL build_spmodes( )
      IF ( .FALSE. ) CALL check_spmodes( ) 
      ! 
   END SUBROUTINE gen_vel_mod
   ! [gen_vel_mod]


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
   SUBROUTINE gen_psobs( )
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
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw1, zw2


      INTEGER ::   jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at

      INTEGER                                 :: ncid

      REAL(wp)  :: val_r, val_c

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'gen_psobs : TLU, generation of pseudo-observations '
         WRITE(numout,*) '~~~~~~~~~'
      END IF
      !
      ALLOCATE( zwa1(jpi,jpj,jpk), zwa2(jpi,jpj,jpk), stat=ierr )
      !
      ! Assign Indexes [ ..., -1, 0, +1, ...] and compute patch size
      !
      rel_idx = [( i, i = - nn_tlu_hsiz, nn_tlu_hsiz )]
      patch_size = nn_tlu_psiz**2
      !
      rdAdz_dVu = dAdzu / SUM( dAdzu )
      rdAdz_dVv = dAdzv / SUM( dAdzv )
      !
      u_psob = 0._wp
      v_psob = 0._wp
      w_psob = 0._wp      
      u_mean = 0._wp
      v_mean = 0._wp
      !
      k = 1 
      !
      DO i = 1, SIZE( rel_idx )
         !
         DO j = 1, SIZE( rel_idx )
            !
            zw1 = 0._wp
            zw2 = 0._wp
            !
            zw1( 2:jpim1, 2:jpjm1, :) = un( 2+rel_idx(i):jpim1+rel_idx(i), 2+rel_idx(j):jpim1+rel_idx(j), :)
            zw2( 2:jpim1, 2:jpjm1, :) = vn( 2+rel_idx(i):jpim1+rel_idx(i), 2+rel_idx(j):jpim1+rel_idx(j), :)
            !
	    ! This communication is fundamental to compute w_psobs near the boundaries
	    !
            CALL lbc_lnk_multi( 'gen_psobs', zw1 , 'U', -1., zw2 , 'V', -1.)
            !
            u_psob(:, k) = reshape( zw1, jpi * jpj * jpk )
            v_psob(:, k) = reshape( zw2, jpi * jpj * jpk )
            !
            ! Compute patch mean 
	    !
            u_mean = u_mean + u_psob(:,k) * rdAdz_dVu
            v_mean = v_mean + v_psob(:,k) * rdAdz_dVv
            !
            k = k + 1
	    !
         END DO
         !
      END DO
      !
      ! Compute fluctuations by Removing average
      !
      u_psob(:,1:patch_size) = u_psob(:,1:patch_size) - spread( u_mean, 2, patch_size)
      v_psob(:,1:patch_size) = v_psob(:,1:patch_size) - spread( v_mean, 2, patch_size)
      !
      call random_number( real_unif_01 )
      pick_id = INT( 1 + FLOOR( ( patch_size ) * real_unif_01 ) , i4)  
      !
      DO jk = 1, jpk
         !
         ! Define zero-th indexed position
         !
         m_idx = ( jk - 1 ) * jpi * jpj 
         !
         DO jjii = 1, jpi * jpj
            !
            u_psob( m_idx + jjii, : ) = u_psob( m_idx + jjii, pick_id(jjii, :) )
            v_psob( m_idx + jjii, : ) = v_psob( m_idx + jjii, pick_id(jjii, :) )
            !
         END DO
         !
      END DO
      !
      CALL uv2w_reshaped(u_psob, v_psob, w_psob)
      !
      DEALLOCATE( zwa1, zwa2 )
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
   SUBROUTINE build_Cmat( )
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
         WRITE(numout,*) 'build_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Cmat
      !
      Cmat = 0._wp
      !
      ! Transpose components ( https://stackoverflow.com/questions/55222312/do-most-compilers-optimize-matmultransposea-b )
      !
      u_psob_t = transpose( u_psobs * dAdzu )
      v_psob_t = transpose( v_psobs * dAdzv )
      w_psob_t = transpose( w_psobs * dAdzt )
      !
      ! Add the inner product: this is euclidean, unless changed in gen_psobs
      !
      Cmat = Cmat + matmul( u_psob_t, u_psob)
      Cmat = Cmat + matmul( v_psob_t, v_psob)
      Cmat = Cmat + matmul( w_psob_t, w_psob)
      Cmat = Cmat / nn_tlu_pobs
      !
      !
      ! COMMUNICATION NEEDED HERE!!!!!!!!!!!
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
   SUBROUTINE solve_Cmat( )
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
         WRITE(numout,*) 'solve_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Eigenvalue matrix by constructing upper triangular of C
      !
      Evec = 0._wp      
      DO i = 1, nn_tlu_pobs
         Evec(i,:) = ( Cmat(i,j), j=i:nn_tlu_pobs )
      END DO
      !
      ! Solve   C a=\lambda a
      !
      !    dsyev( JOBZ, UPLO,           N,    A,         LDA,    W, WORK, LWORK, INFO )
      CALL dsyev(  'V',  'U', nn_tlu_pobs, Evec, nn_tlu_pobs, Eval, WORK, LWORK, INFO )
      !
      ! Transpose components ( https://stackoverflow.com/questions/55222312/do-most-compilers-optimize-matmultransposea-b )
      !
      u_psobs_t = matmul( Evec, u_psobs_t)
      v_psobs_t = matmul( Evec, v_psobs_t)
      w_psobs_t = matmul( Evec, w_psobs_t)
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
   SUBROUTINE build_spmodes( )
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

      REAL(wp)  :: val

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'build_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      DO jo = 1, nn_tlu_pobs
         !
         ! Define zero-th indexed position
         !
         m_idx = ( jo - 1 ) * jpk 
         !
         spx_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( u_psobs_t( jo, :), (jpi, jpj, jpk) ) * umask * val
         spy_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( v_psobs_t( jo, :), (jpi, jpj, jpk) ) * vmask * val
         spz_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( w_psobs_t( jo, :), (jpi, jpj, jpk) ) * wmask * val
         !
      END DO
      !      
   END SUBROUTINE build_spmodes
   ! [build_spmodes]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE check_spmodes ***
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
   SUBROUTINE check_spmodes( )
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

      REAL(wp)  :: val

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'check_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      DO jo = 1, nn_tlu_pobs
         !
         ! Define zero-th indexed position
         !
         m_idx = ( jo - 1 ) * jpk 
         !
         spx_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( u_psobs_t( jo, :), (jpi, jpj, jpk) ) * umask * val
         spy_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( v_psobs_t( jo, :), (jpi, jpj, jpk) ) * vmask * val
         spz_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( w_psobs_t( jo, :), (jpi, jpj, jpk) ) * wmask * val
         !
      END DO
      !      
   END SUBROUTINE check_spmodes
   ! [check_spmodes]

END MODULE tlusvd


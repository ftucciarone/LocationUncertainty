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
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psx_mod                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psy_mod                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: psz_mod                 !> @public   W-velocity modes
   ! [tlu_spmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),    PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: gaus_rv                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_psob] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: u_psob                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: v_psob                 !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: w_psob                 !> @public   W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: u_psot                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: v_psot                 !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: w_psot                 !> @public   W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: u_mean                 !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: v_mean                 !> @public   V-velocity pseudo-observations 
   ! [tlu_psob]
   !
   ! [tlu_eigen] 
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Cglo                   !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Cloc                   !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Evec                   !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: Eval                   !> @public   W-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: work
   INTEGER,                          SAVE                   :: lwork
   ! [tlu_eigen]   
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dAu, dAv, dAt          !> @public   U-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dzu, dzv, dzw          !> @public   V-velocity pseudo-observations
   REAL(wp),            ALLOCATABLE, SAVE, DIMENSION(:)     :: dVu, dVv, dVw          !> @public   V-velocity pseudo-observations   
   REAL(wp),            ALLOCATABLE,       DIMENSION(:)     :: rdAdz_dVu, rdAdz_dVv
   ! [tlu_randgauss] 
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: real_unif_01 
   INTEGER(i4),         ALLOCATABLE,       DIMENSION(:,:)   :: pick_id 
   ! [tlu_randgauss] 
   !
   !
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk1, wrk2


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
      INTEGER(i4) :: ios, ierr(8), chk   ! namelist output, allocation statistics
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
         WRITE(numout,*) '           data-free noise model        ln_tlu_svd = ', ln_tlu_svd
         WRITE(numout,*) '   number of pseudo-observations       nn_tlu_nmod = ', nn_tlu_pobs
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      !
      ALLOCATE( psx_mod( jpi, jpj, nn_tlu_pobs * jpk ), &
      &         psy_mod( jpi, jpj, nn_tlu_pobs * jpk ), &
      &         psz_mod( jpi, jpj, nn_tlu_pobs * jpk ), stat=ierr(1) )   ! the singular values
      !
      psx_mod = 0._wp
      psy_mod = 0._wp
      psz_mod = 0._wp
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( gaus_rv( nn_tlu_pobs ), stat=ierr(2) )
      !
      gaus_rv = 0._wp 
      !
      ! Allocate Internal variables 
      !
      ALLOCATE(    u_psob( nn_tlu_pobs, jpi * jpj * jpk ), &
      &            v_psob( nn_tlu_pobs, jpi * jpj * jpk ), &
      &            w_psob( nn_tlu_pobs, jpi * jpj * jpk ), &
      &            u_mean(              jpi * jpj * jpk ), &
      &            v_mean(              jpi * jpj * jpk ), stat=ierr(3) )   ! the singular values
      !
      ! Allocate Internal variables 
      !
      ALLOCATE(  u_psot( jpi * jpj * jpk, nn_tlu_pobs ), &
      &          v_psot( jpi * jpj * jpk, nn_tlu_pobs ), &
      &          w_psot( jpi * jpj * jpk, nn_tlu_pobs ), stat=ierr(4) )     
      !
      ! Allocate scale factors for online POD
      !
      ALLOCATE(       dAu( jpi * jpj * jpk ), &
      &               dAv( jpi * jpj * jpk ), &
      &               dAt( jpi * jpj * jpk ), &
      &               dzu( jpi * jpj * jpk ), &
      &               dzv( jpi * jpj * jpk ), &
      &               dzw( jpi * jpj * jpk ), &
      &               dVu( jpi * jpj * jpk ), &
      &               dVv( jpi * jpj * jpk ), &
      &               dVw( jpi * jpj * jpk ), &
      &         rdAdz_dVu( jpi * jpj * jpk ), &
      &         rdAdz_dVv( jpi * jpj * jpk ), stat=ierr(5) )
      !
      ! Allocate the unifor random choices for the patches
      !
      ALLOCATE( real_unif_01( nn_tlu_pobs, jpi * jpj ), &
      &              pick_id( nn_tlu_pobs, jpi * jpj ), stat=ierr(6) )   ! the singular values
      !
      ! Allocate Matrix for eigen problem
      !
      ALLOCATE( Cglo( nn_tlu_pobs, nn_tlu_pobs ),       &
      &         Cloc( nn_tlu_pobs, nn_tlu_pobs ),       &
      &         Evec( nn_tlu_pobs, nn_tlu_pobs ),       &
      &         Eval( nn_tlu_pobs),                     stat=ierr(7) )
      !
      ! Allocate workspaces
      !
      ALLOCATE( wrk1(jpi,jpj,jpk), wrk2(jpi,jpj,jpk), stat=ierr(8) )


      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_pod_init(): allocation failed = ', ierr  
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !
      ! Area factors does not change in time, can be set now
      !
      dAu = reshape( spread( e1e2u,3,jpk), (/ jpi * jpj * jpk /) )
      dAv = reshape( spread( e1e2v,3,jpk), (/ jpi * jpj * jpk /) )
      dAt = reshape( spread( e1e2t,3,jpk), (/ jpi * jpj * jpk /) )
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

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'gen_vel_mod : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      ! Update scale factors for online POD
      !
      dVu = dAu * reshape( e3u_n, (/ jpi * jpj * jpk /) )
      dVv = dAv * reshape( e3v_n, (/ jpi * jpj * jpk /) )
      dVw = dAt * reshape( e3w_n, (/ jpi * jpj * jpk /) )
      !
      ! Generate Pseudo-Observations
      !
      call cpu_time(tic(1))
      CALL gen_psobs( )
      call cpu_time(toc(1))
      !
      ! Build Correlation matrix
      !
      call cpu_time(tic(2))
      CALL build_Cmat( )
      call cpu_time(toc(2))
      !
      ! Solve eigenvalue problem
      !
      call cpu_time(tic(3))
      CALL solve_Cmat( kt )
      call cpu_time(toc(3))
      !
      ! Build spatial modes
      !
      call cpu_time(tic(4))
      CALL build_spmodes( )
      call cpu_time(toc(4))
      !
      !print *, toc(1)-tic(1), toc(2)-tic(2), toc(3)-tic(3), toc(4)-tic(4), toc(4)-tic(1)
      !CALL check_spmodes( )
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


      INTEGER ::   jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   i, j, k, jk, jjii
      INTEGER ::   rel_idx(nn_tlu_psiz)

      INTEGER                                 :: patch_size

      REAL(wp)  :: val_r, val_c

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'gen_psobs : TLU, generation of pseudo-observations '
         WRITE(numout,*) '~~~~~~~~~'
      END IF
      !
      !
      ! Assign Indexes [ ..., -1, 0, +1, ...] and compute patch size
      !
      rel_idx = [( i, i = - nn_tlu_hsiz, nn_tlu_hsiz )]
      patch_size = nn_tlu_psiz**2
      !
      rdAdz_dVu = dVu / SUM( dVu )
      rdAdz_dVv = dVv / SUM( dVv )
      !
      u_psob = 0._wp
      v_psob = 0._wp
      w_psob = 0._wp      
      u_mean = 0._wp
      v_mean = 0._wp
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
            wrk1( 2:jpim1, 2:jpjm1, :) = un( 2+rel_idx(i):jpim1+rel_idx(i),  &
            &                                2+rel_idx(j):jpim1+rel_idx(j), :)
            wrk2( 2:jpim1, 2:jpjm1, :) = vn( 2+rel_idx(i):jpim1+rel_idx(i),  &
            &                                2+rel_idx(j):jpim1+rel_idx(j), :)
            !
            ! This communication is fundamental to compute w_psobs near the boundaries
            !
            CALL lbc_lnk_multi( 'gen_psobs', wrk1 , 'U', -1., wrk2 , 'V', -1.)
            !
            u_psob(k, :) = reshape( wrk1, (/ jpi * jpj * jpk /) )
            v_psob(k, :) = reshape( wrk2, (/ jpi * jpj * jpk /) )
            !
            k = k + 1
            !
         END DO
         !
      END DO
      !
      ! Compute patch mean 
      !
      u_mean = SUM( u_psob, 1) * rdAdz_dVu
      v_mean = SUM( v_psob, 1) * rdAdz_dVv      
      !
      ! Compute fluctuations by Removing average
      !
      u_psob(1:patch_size, :) = u_psob(1:patch_size, :) - spread( u_mean, 1, patch_size)
      v_psob(1:patch_size, :) = v_psob(1:patch_size, :) - spread( v_mean, 1, patch_size)
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
            u_psob(:, m_idx + jjii) = u_psob(pick_id(:, jjii), m_idx + jjii)
            v_psob(:, m_idx + jjii) = v_psob(pick_id(:, jjii), m_idx + jjii)
            !
         END DO
         !
      END DO
      !
!      CALL uv2w_reshaped(u_psob, v_psob, w_psob)
      !
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

      INTEGER ::   jm, m_idx, ierr, rank ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at
      INTEGER ::   i, j
      INTEGER                                 :: ncid

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'build_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Cglo
      !
      Cloc = 0._wp
      !
      u_psot = TRANSPOSE( u_psob * SPREAD( dVu, 1, nn_tlu_pobs) )
      v_psot = TRANSPOSE( v_psob * SPREAD( dVv, 1, nn_tlu_pobs) )
!     w_psot = TRANSPOSE( w_psob * SPREAD( dVw, 1, nn_tlu_pobs) )
      !
      Cloc = Cloc + MATMUL( u_psob, u_psot) / nn_tlu_pobs
      Cloc = Cloc + MATMUL( v_psob, v_psot) / nn_tlu_pobs
!     Cloc = Cloc + MATMUL( w_psob, w_psot) / nn_tlu_pobs
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
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'solve_Cmat : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~'
      END IF
      !
      ! Initialize Eigenvalue matrix by constructing upper triangular of C
      !
      Evec = Cglo
      ! 
      ! Allocate the workspaces for LAPACK's dsyev
      !
      IF ( (kt == nit000) ) THEN 
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
         Cglo(i, :) = Evec(:, i) / ( nn_tlu_pobs * SQRT( Eval(i) ) ) 
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
      
      u_psob = matmul( Cglo, u_psob)
      v_psob = matmul( Cglo, v_psob)
      w_psob = 0._wp !matmul( Evec, w_psobs)

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
      INTEGER ::   jo, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at

      INTEGER        :: ncid

      REAL(wp)  :: val

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'build_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~~'
      END IF
      !
      val = 1._wp
      !
      DO jo = 1, nn_tlu_pobs
         !
         ! Define zero-th indexed position
         !
         m_idx = ( jo - 1 ) * jpk 
         !
         psx_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( u_psob( jo, :), (/ jpi, jpj, jpk /) ) * umask * val
         psy_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( v_psob( jo, :), (/ jpi, jpj, jpk /) ) * vmask * val
         psz_mod(:,:,m_idx + 1 : m_idx + jpk ) = reshape( w_psob( jo, :), (/ jpi, jpj, jpk /) ) * wmask * val
         !
      END DO
      !      
   END SUBROUTINE build_spmodes
   ! [build_spmodes]


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
      INTEGER ::   jo, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   i, j, k

      INTEGER        :: ncid

      REAL(wp)  ::  eigenval(nn_tlu_pobs, nn_tlu_pobs)
      REAL(wp)  ::  delta(nn_tlu_pobs, nn_tlu_pobs)

      REAL(wp)  ::  u_psob_t_locl( jpi * jpj * jpk, nn_tlu_pobs)


      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'check_spmodes : TLU, Dynamical Mode Decomposition reading modes '
         WRITE(numout,*) '~~~~~~~~~~~~~'
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
         print *, (eigenval(i, j) /  Eval(i), j = 1, nn_tlu_pobs), Eval(i)         !
      end do
      !
      delta = delta + matmul(u_psob, TRANSPOSE( u_psob * SPREAD( dVu, 1, nn_tlu_pobs)  ) )
      delta = delta + matmul(v_psob, TRANSPOSE( v_psob * SPREAD( dVv, 1, nn_tlu_pobs)  ) )
!     delta = delta + matmul(w_psob, TRANSPOSE( w_psob * SPREAD( dVw, 1, nn_tlu_pobs)  ) )
      !
      print *, " " 
      print *, " " 
      do i = 1, nn_tlu_pobs
         !
         print *, (delta(i, j), j = 1, nn_tlu_pobs) 
         !
      end do

      stop
      !      
   END SUBROUTINE check_spmodes
   ! [check_spmodes]


END MODULE tlusvd


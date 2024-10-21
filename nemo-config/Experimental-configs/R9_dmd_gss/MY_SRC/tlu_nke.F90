!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlunke
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2023 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!
!!------------------------------------------------------------------------------
MODULE tlunke
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE tlu              ! Initialization of stochastic structures
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_nke_init
   PUBLIC tlu_nke
   PUBLIC tlu_amp
   ! [public_sub]
   !
   INCLUDE 'mpif.h'
   !
   INTEGER,  PARAMETER                          ::   nsamples = 10     !< @private number of samples in the ensemble
   ! [tlu_noise]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  unke   !< @public Modified advection: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  vnke   !< @public Modified advection: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  wnke   !< @public Modified advection: z component
   ! [tlu_noise]
   !
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: var_sav
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: uisd_s   !< @public Ito-Stokes drift: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: visd_s   !< @public Ito-Stokes drift: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wisd_s   

   !
   REAL(wp),            PARAMETER :: vel = 0.005_wp                       ! prescribed amplitude of noise [m/s]
   REAL(wp)                       :: vel_nke
   REAL(wp)                       :: bia_nke
   REAL(wp)                       :: amp   
   REAL(wp), SAVE                 :: dom_siz
   !
   !! * Substitutions [TODO] enable substitutions for NEMO
! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------

CONTAINS

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> P. Derian, P. Chandramouli, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief   Initialize module variables
   !!
   !!  @details
   !!  Reads namelist, compute other parameters and allocate arrays.
   !! 
   !!
   !! 
   !! 
   !> @param[in]     ln_tlu: logical switch for location uncertainty
   !! @param[in]     jpi: the first dimension of the spatial arrays
   !! @param[in]     jpj: the first dimension of the spatial arrays
   !! @param[in]     jpk: the first dimension of the spatial arrays
   !! @param[in]     numout: writing unit 
   !! @result        Arrays allocated
   !! @warning  
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_nke_init 
   ! [tlu_nke_init]
   SUBROUTINE tlu_nke_init
      INTEGER(i4) :: ierr   ! namelist output, allocation statistics
      !
      !!------------------------------------------------------------------------
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_nke_init : Transport under Location Uncertainty, Diagnostic routine '
         WRITE(numout,*) '~~~~~~~~~~~~'
      END IF
      !
      !
      ! Allocate
      !------------------------------
      ALLOCATE(    unke(jpi,jpj,jpk), &
      &            vnke(jpi,jpj,jpk), &
      &            wnke(jpi,jpj,jpk), &
      &          uisd_s(jpi,jpj,jpk), &
      &          visd_s(jpi,jpj,jpk), &
      &          wisd_s(jpi,jpj,jpk), &
      &         var_sav(jpi,jpj,jpk,6), stat = ierr )  
      !
      ! Allocate
      !------------------------------
      !
      ! Initialise
      !------------------------------
      !
      !
      IF (ierr/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_nke_init(): allocation failed = ', ierr   ! [TODO] should be ctmp1? instead of numout
         !CALL ctl_stop( ctmp1 )   !< @todo enable
         STOP
      END IF
      !
      !
   END SUBROUTINE tlu_nke_init
   ! [tlu_nke_init]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_fields ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Compute the fields necessary for the method
   !!
   !> @details       Redirects to the particular initialization procedures that 
   !!                are required to each kind of noise that is chosen in the 
   !!                namelist. Calls the specific initialization routine of the 
   !!                method.
   !!  
   !> @note          Uses the modules tlupod, tludmd, tlupso, tluwlt, tlugss
   !!
   !> @param[in]     ln_tlu_pod (as a global variable)
   !> @param[in]     ln_tlu_dmd (as a global variable)
   !> @param[in]     ln_tlu_pso (as a global variable)
   !> @param[in]     ln_tlu_wlt (as a global variable)
   !> @param[in]     ln_tlu_gss (as a global variable)
   !! 
   !! @result        All the structures are initialized
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_fields
   ! [tlu_fields]
   SUBROUTINE tlu_nke( kt )
      USE tlurnd
      USE tlupod
      USE tludmd
      USE tlupso
      USE tluwlt
      USE tlugss
      USE tluprj
      !
      INTEGER,  INTENT(in   ) :: kt         ! ocean time-step index
      INTEGER                 :: ierr, myID
      INTEGER                 :: j
      !
      IF( ln_timing ) CALL timing_start('tlu_nke') ![NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_nke : TLU fields construction'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '

         CALL tlu_nke_init
         !
         var_sav = var_ten
         !
         IF ( ln_pyp ) THEN
            uisd_s = uisd_0
            visd_s = visd_0
            wisd_s = wisd_0
         ELSE
            uisd_s = uisd_n
            visd_s = visd_n
            wisd_s = wisd_n
         END IF
         !
         unke = 1._wp
         vnke = 1._wp
         wnke = 1._wp
         !
         CALL int3d_ene3c( unke, vnke, wnke, 'dom', dom_siz, 1._wp )
         !
      ENDIF
      !
      unke = 0._wp
      vnke = 0._wp
      wnke = 0._wp
      !
      DO j = 1, nsamples
         !
         ! Initialize the noise
         !
         unoi = 0._wp
         vnoi = 0._wp
         wnoi = 0._wp      
         !
         ! POD-Based noise
         !
         IF ( ln_tlu_pod ) THEN
            ! 
            ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
            !
            CALL tlu_tcoef_pod( kt, brwn_rv )
            CALL tlu_noi( pod_xmd, pod_ymd, nn_tlu_nmod, brwn_rv, .FALSE. )
            !       
         END IF
         !
         ! DMD-Based noise
         !
         IF ( ln_tlu_dmd ) THEN
            ! 
            ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
            !
            CALL tlu_tcoef_dmd( kt,   rtime_coef,    itime_coef )
            CALL tlu_noi( dmd_rxmd_r, dmd_rymd_r, nn_tlu_nmod_r, rtime_coef,  .FALSE. )
            CALL tlu_noi( dmd_ixmd_r, dmd_iymd_r, nn_tlu_nmod_r, itime_coef,  .FALSE. )
            !   
         END IF
         !
         ! SVD-Based noise
         !
         IF ( ln_tlu_pso ) THEN
            !
            ! Generate the random modes
            !
            CALL gen_vel_mod( kt )
            !
            ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
            !
            CALL Standard_Gaussian( nn_tlu_pobs,  gaus_rv ) !, allones = .FALSE. )  
            CALL MPI_BCAST(  gaus_rv,   nn_tlu_pobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
            CALL tlu_noi( pso_xmd, pso_ymd, nn_tlu_pobs, gaus_rv, .FALSE. )
            !
         END IF      
         !
         ! Wavelet Based noise
         !
         IF ( ln_tlu_wlt ) THEN
            !
            ! Generate the random modes
            !
            CALL gen_wlt_mod( kt )
            !
            unoi = unoi + wlx_noi * umask
            vnoi = vnoi + wly_noi * vmask
            !
         END IF 
         !
         ! Gaussian white noise
         !
         IF ( ln_tlu_gss ) THEN
            !
            ! Generate the random modes
            !
            CALL gen_white_noise( kt )
            !
            unoi = unoi + u_wht * umask
            vnoi = vnoi + v_wht * vmask
            !          
         END IF   
         !
         ! Compute vertical components
         !
         CALL build_vmodes(wnoi, unoi, vnoi)
         !
         unke = unke + ( unoi**2 ) 
         vnke = vnke + ( vnoi**2 )
         wnke = wnke + ( wnoi**2 )
         !
      END DO
      !
      ! Compute expectation by averaging
      !
      unke = unke / nsamples
      vnke = vnke / nsamples
      wnke = wnke / nsamples
      !
      CALL int3d_ene3c( unke, vnke, wnke, 'Unke', vel_nke, 0.5_wp )
      !
      amp = vel * SQRT( dom_siz / vel_nke )

     ! call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
     ! IF ( myID == 0 ) THEN
     !    print '(A36, E16.7)', TRIM('       ')//'    final amplitude value:', amp
     ! END IF      
      !
      IF( ln_timing )  CALL timing_stop('tlu_nke')   ! [NEMO] check
      !
   END SUBROUTINE tlu_nke


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_noi ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Computes one realization of the noise
   !!
   !> @details       Computes one realization of the noise as the sum of a set of 
   !!                noise modes in each direction, x_mod, y_mod and z_mod, times
   !!                a given set of temporal coefficients time_coeff.
   !!                \f{align*}{
   !!                  \sigma\mathrm{d}\mathbf{B}_{t} = \sum_{k=1}^{N} \phi_{k}(x)a_{k}(t) 
   !!                \f}
   !!                where $\phi_{k}(x)$ is represented by x_mod, y_mod, z_mod,
   !!                $a_{k}(t)$  is represented by time_coeff, $N$ is nn_nmod.
   !!                The results is:
   !!                \f{align*}{
   !!                  \texttt{unoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{x} \\
   !!                  \texttt{vnoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{y} \\
   !!                  \texttt{wnoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{z} 
   !!               \f}
   !!  
   !> @note          It may or may not communicate the halo to other processors
   !!                depending on the value of lbc_comm
   !!
   !> @param[in]     x_mod 
   !> @param[in]     y_mod 
   !> @param[in]     z_mod 
   !> @param[in]     nn_nmod 
   !> @param[in]     time_coeff 
   !> @param[in]     lbc_comm 
   !! 
   !! @result        Adds to the noise field the state given by the input modes
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_noi
   ! [tlu_noi]
   SUBROUTINE tlu_amp( kt )
     !
     INTEGER,  INTENT(in   ) :: kt         ! ocean time-step index
      INTEGER                 :: ierr, myID
     !
     IF( kt == nit000 + dt_delay )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_nke : TLU fields construction'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      END IF
      !      
      !
      unoi = amp * unoi
      vnoi = amp * vnoi
      wnoi = amp * wnoi
      !
      var_ten = (amp**2) * var_sav 
      !
      IF ( ln_pyp ) THEN
         uisd_0 = (amp**2) * uisd_s
         visd_0 = (amp**2) * visd_s
         wisd_0 = (amp**2) * wisd_s
      ELSE
         uisd_n = (amp**2) * uisd_s
         visd_n = (amp**2) * visd_s
         wisd_n = (amp**2) * wisd_s
      END IF
      !
      IF (ln_tlu_bia) THEN
         !
         CALL int3d_ene3c( ubia_n**2, vbia_n**2, wbia_n**2, 'dom', bia_nke, 0.5_wp )
         !
         ubia_n = ubia_n * SQRT( dom_siz / bia_nke ) * 0.0005_wp
         vbia_n = vbia_n * SQRT( dom_siz / bia_nke ) * 0.0005_wp
         wbia_n = wbia_n * SQRT( dom_siz / bia_nke ) * 0.0005_wp
         !
         !call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
         !IF ( myID == 0 ) THEN
         !   print '(A36, E16.7)', TRIM('       ')//'    final amplitude value:', amp
         !END IF  
         !
      END IF
      !
      !
   END SUBROUTINE tlu_amp
   ! [tlu_noi]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_noi ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Computes one realization of the noise
   !!
   !> @details       Computes one realization of the noise as the sum of a set of 
   !!                noise modes in each direction, x_mod, y_mod and z_mod, times
   !!                a given set of temporal coefficients time_coeff.
   !!                \f{align*}{
   !!                  \sigma\mathrm{d}\mathbf{B}_{t} = \sum_{k=1}^{N} \phi_{k}(x)a_{k}(t) 
   !!                \f}
   !!                where $\phi_{k}(x)$ is represented by x_mod, y_mod, z_mod,
   !!                $a_{k}(t)$  is represented by time_coeff, $N$ is nn_nmod.
   !!                The results is:
   !!                \f{align*}{
   !!                  \texttt{unoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{x} \\
   !!                  \texttt{vnoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{y} \\
   !!                  \texttt{wnoi} & += \sigma\mathrm{d}\mathbf{B}_{t}^{z} 
   !!               \f}
   !!  
   !> @note          It may or may not communicate the halo to other processors
   !!                depending on the value of lbc_comm
   !!
   !> @param[in]     x_mod 
   !> @param[in]     y_mod 
   !> @param[in]     z_mod 
   !> @param[in]     nn_nmod 
   !> @param[in]     time_coeff 
   !> @param[in]     lbc_comm 
   !! 
   !! @result        Adds to the noise field the state given by the input modes
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_noi
   ! [tlu_noi]
   SUBROUTINE tlu_noi( x_mod, y_mod, nn_nmod, time_coeff, lbc_comm )
      !      
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: x_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: y_mod
      REAL(wp), DIMENSION(:),     INTENT(in   ) :: time_coeff
      INTEGER,                    INTENT(in   ) :: nn_nmod
      LOGICAL,                    INTENT(in   ) :: lbc_comm
      INTEGER                                   :: jm, m_idx      
      !
      IF( ln_timing ) CALL timing_start('tlu_noi') ![NEMO] check
      !      
      DO jm = 1, nn_nmod
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
         !
         unoi(:,:,:) = unoi(:,:,:) + x_mod(:,:,m_idx + 1 : m_idx + jpk ) * umask * time_coeff(jm)
         vnoi(:,:,:) = vnoi(:,:,:) + y_mod(:,:,m_idx + 1 : m_idx + jpk ) * vmask * time_coeff(jm)
         !
      ENDDO
      !      
      ! Compute vertical velocity
      !
      CALL build_vmodes( wnoi, unoi, vnoi)      
      !
      ! Lateral boundary condition transfer across nodes
      !
      IF ( lbc_comm ) THEN
         !
         CALL lbc_lnk_multi( 'tlu_noi', unoi , 'U', -1., vnoi , 'V', -1.,  wnoi , 'T', -1.  )
         !
      END IF
      !
   END SUBROUTINE tlu_noi
   ! [tlu_noi]


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
      wout = 0._wp
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         wout( 2:jpi-1, 2:jpj-1, jk) =     wout( 2:jpi-1, 2:jpj-1, jk + 1) -   &
         !
         ! ( dx( dy*dz*u_in )
         &                           (      e2u( 2:jpi-1, 2:jpj-1 )            &
         &                           *    e3u_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                           *      uin( 2:jpi-1, 2:jpj-1, jk )        &
         &                           -      e2u( 1:jpi-2, 2:jpj-1 )            &
         &                           *    e3u_n( 1:jpi-2, 2:jpj-1, jk )        &
         &                           *      uin( 1:jpi-2, 2:jpj-1, jk )        &
         !   dy ( dx*dz*v_isd )
         &                           +      e1v( 2:jpi-1, 2:jpj-1 )            &
         &                           *    e3v_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                           *      vin( 2:jpi-1, 2:jpj-1, jk )        &
         &                           -      e1v( 2:jpi-1, 1:jpj-2 )            &
         &                           *    e3v_n( 2:jpi-1, 1:jpj-2, jk )        &
         &                           *      vin( 2:jpi-1, 1:jpj-2, jk ) )      &
         !   .../dV ) * dz
         &                           *    tmask( 2:jpi-1, 2:jpj-1, jk)         & 
         &                           * r1_e1e2t( 2:jpi-1, 2:jpj-1 )
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
   END SUBROUTINE build_vmodes   
   ! [tlu_wzvcmp]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE int3d_ene3c ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes a 2-dimensional "energy-like" integral of a 
   !!                3-components field, i.e.
   !!
   !!                    \int\int u_{x}^{2} + u_{y}^{2} + u_{z}^{2} dxdy
   !!                
   !!                
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     cx: x-component 
   !> @param[in]     cy: y-component
   !> @param[in]     cz: z-component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this int3d_ene3c
   ! [int3d_ene3c]
   SUBROUTINE int3d_ene3c( cx, cy, cz, fieldname, val, norm_const )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::  cx, cy, cz    ! components 
      CHARACTER(len = * )         , INTENT(in   ) ::  fieldname     ! string to print
      REAL(wp)                    , INTENT(  out) ::  val           ! components 
      REAL(wp)                                    ::  norm_const    ! components 
      !
      INTEGER                                     ::  ji, jj, jk    ! dummy loop indices
      REAL(wp)                                    ::  sub_int       ! 2D workspace
      REAL(wp)                                    ::  dom_int       ! 2D workspace
      REAL(wp), DIMENSION(mppsize)                ::  zw1
      !
      INTEGER                                     ::  ierr, myID

      !!----------------------------------------------------------------------
      !
      sub_int = 0._wp
      !
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               sub_int = sub_int + ( e1u(ji,jj) * e2u(ji,jj) * e3u_n(ji,jj,jk) * cx(ji,jj) &
                               &   + e1v(ji,jj) * e2v(ji,jj) * e3v_n(ji,jj,jk) * cy(ji,jj) &
                               &   + e1t(ji,jj) * e2t(ji,jj) * e3w_n(ji,jj,jk) * cz(ji,jj) ) * norm_const

            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      

      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLGATHER( sub_int, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)
      val = SUM(zw1)

      !call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
      !IF ( myID == 0 ) THEN
      !   print *, TRIM(fieldname)//' 3D domain integral value:', val
      !END IF

   END SUBROUTINE int3d_ene3c
   ! [int3d_ene3c]


END MODULE tlunke





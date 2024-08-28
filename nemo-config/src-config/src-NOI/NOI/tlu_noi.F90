!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_noi
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                P. Derian, P. Chandramouli, F. Tucciarone
!!
!> @version
!!                2017 -  ( P. DERIAN       )  Original code <br>
!!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
!!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Contains the initialization of the noise and the procedure-independent
!!                routines for noise generation
!! 
!! 
!! @par           Procedure specifics      
!> @details       Each routine in this module is either independent of the noise 
!!                chosen (viz. the Ito-Stokes drift) or discriminates between the 
!!                noise models and refers to the relative subroutines.
!!
!> @snippet       this tlu_noise
!!                \f{align*}{
!!                  \texttt{unoi} &\leftarrow d\boldsymbol{\eta}_{t,x} \\
!!                  \texttt{vnoi} &\leftarrow d\boldsymbol{\eta}_{t,y} \\
!!                  \texttt{wnoi} &\leftarrow d\boldsymbol{\eta}_{t,z} 
!!               \f}
!!
!> @snippet      this tlu_bias
!!               \f{align*}{
!!                  \texttt{ubias} &\leftarrow d\boldsymbol{\eta}_{t,x}^{b} \\
!!                  \texttt{vbias} &\leftarrow d\boldsymbol{\eta}_{t,y}^{b} \\
!!                  \texttt{wbias} &\leftarrow d\boldsymbol{\eta}_{t,z}^{b} 
!!               \f}
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
!!
!> @par           Other modules dependencies
!> @snippet       this mod_dep
!!
!> @par           Public variables
!> @snippet       this public_vars
!!
!> @par           Public routines
!> @snippet       this public_sub
!!
!> @todo          read namelist from _cfg, not only _ref
!! @todo          write error messages to appropriate unit
!!
!!------------------------------------------------------------------------------
MODULE tlunoi
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
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC tlu_noi_init    ! called by tlu_init in tlu.f90
   PUBLIC tlu_fields      ! called by step.f90
   PUBLIC tlu_noi         ! called by tlu_fields and noise_energy in tlu_dgns.F90
   !
   INCLUDE 'mpif.h'
   !
   REAL(wp),    PUBLIC    :: cof_tke     
   REAL(wp)               :: KE_ratio_b, KE_ratio_n 
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_noi_init ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Initialize the noise model chosen
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
   !! @result        All the strictures are initialized
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_noi_init
   ! [tlu_noi_init]
   SUBROUTINE tlu_noi_init
      USE tlupod
      USE tludmd
      USE tlupso
      USE tluwlt
      USE tlugss
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_noi_init : TLU, Proper orthogonal decomposition  noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*)
      END IF
      !
      IF ( ln_tlu_pod ) CALL tlu_init_pod ! POD-Based noise
      IF ( ln_tlu_dmd ) CALL tlu_init_dmd ! DMD-Based noise
      IF ( ln_tlu_pso ) CALL tlu_init_pso ! Pseudo-Observations model
      IF ( ln_tlu_wlt ) CALL tlu_init_wlt ! Pseudo-Observations model
      IF ( ln_tlu_gss ) CALL tlu_init_gss ! White Gaussian noise
      !
      CALL tlu_arrays_init( nit000 )
      IF ( dt_delay .ne. 0 ) ld_tlu = .false.
      !
   END SUBROUTINE tlu_noi_init
   ! [tlu_noi_init]


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
   SUBROUTINE tlu_fields( kt )
      USE tlurnd
      USE tlupod
      USE tludmd
      USE tlupso
      USE tluwlt
      USE tlugss
      USE tluprj
      USE tlunke
      !
      INTEGER, INTENT(in   ) :: kt         ! ocean time-step index
      INTEGER                :: ierr
      !
      IF( ln_timing ) CALL timing_start('tlu_noi') ![NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_fields : TLU fields construction'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF
      !
      CALL tlu_arrays_init( kt ) 
      !
      ! POD-Based noise
      !
      IF ( ln_tlu_pod ) THEN
         ! 
         IF( (kt == nit000 + dt_delay ) ) THEN      ! Stationary assumption here
            !
            ! Compute the variance tensor a = sum_{k} \phi_{k}^{0}\phi_{k}^{0} * dt and Ito-Stokes drift
            !
            !    tlu_var(   x_mod,   y_mod,   z_mod,     nn_nmod,     dt, lbc_comm )
            CALL tlu_var( pod_xmd, pod_ymd, pod_zmd, nn_tlu_nmod, rn_rdt,   .TRUE. )
            CALL tlu_isd( uisd_n, visd_n, wisd_n )
            !         
            ! Transfer initial fields into now fields for later projection
            !
            IF ( ln_pyp ) CALL tlu_tnsf_isd  
            !
         END IF
         !
         ! Create one realisation of bias \mu_{t}
         !
         IF ( ln_tlu_bia ) CALL tlu_bia_pod ( kt )         
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         CALL tlu_tcoef_pod( kt, brwn_rv )
         CALL tlu_noi( pod_xmd, pod_ymd, nn_tlu_nmod, brwn_rv, .TRUE. )
         !       
      END IF
      !
      ! DMD-Based noise
      !
      IF ( ln_tlu_dmd ) THEN
         ! 
         IF( (kt == nit000 + dt_delay ) .and. ( .not. only_dmd_bias ) ) THEN      ! Stationary assumption here
            !
            ! Compute the variance tensor a = sum_{k} \phi_{k}^{0}\phi_{k}^{0} * dt and Ito-Stokes drift
            !
            !    tlu_var(      x_mod,      y_mod,      z_mod,       nn_nmod,     dt, lbc_comm )      
            CALL tlu_var( dmd_rxmd_r, dmd_rymd_r, dmd_rzmd_r, nn_tlu_nmod_r, rn_rdt,  .FALSE. )
            CALL tlu_var( dmd_ixmd_r, dmd_iymd_r, dmd_izmd_r, nn_tlu_nmod_r, rn_rdt,   .TRUE. )
            CALL tlu_isd( uisd_n, visd_n, wisd_n )
            !         
            ! Transfer initial fields into now fields for projection
            !
            IF ( ln_pyp ) CALL tlu_tnsf_isd  
            !
         END IF
         !
         ! Create one realisation of bias \mu_{t}
         !
         IF (ln_tlu_bia) CALL tlu_bia_dmd ( kt, omega_c )
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         IF ( .not. only_dmd_bias ) THEN
           !
           CALL tlu_tcoef_dmd( kt,   rtime_coef,    itime_coef )
           CALL tlu_noi( dmd_rxmd_r, dmd_rymd_r, nn_tlu_nmod_r, rtime_coef,  .FALSE. )
           CALL tlu_noi( dmd_ixmd_r, dmd_iymd_r, nn_tlu_nmod_r, itime_coef,   .TRUE. )
           !
         END IF
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
         ! Compute the variance tensor a = sum_{k} \phi_{k}^{0}\phi_{k}^{0} * dt and Ito-Stokes drift
         !
         CALL tlu_var( pso_xmd, pso_ymd, pso_zmd, nn_tlu_pobs, rn_rdt,   .TRUE. )  
         CALL tlu_isd( uisd_n, visd_n, wisd_n )
         !         
         ! Transfer initial fields into now fields for projection
         !
         IF ( ln_pyp ) CALL tlu_tnsf_isd  
         !
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         CALL Standard_Gaussian( nn_tlu_pobs,  gaus_rv ) !, allones = .FALSE. )  
         CALL MPI_BCAST(  gaus_rv,   nn_tlu_pobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
         CALL tlu_noi( pso_xmd, pso_ymd, nn_tlu_pobs,    gaus_rv,   .TRUE. )
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
         unoi = wlx_noi * umask
         vnoi = wly_noi * vmask
         !
      END IF
      !
      ! Compute vertical components
      !
      CALL build_vmodes(wnoi, unoi, vnoi)
      IF (ln_tlu_bia) CALL build_vmodes(wbia_n, ubia_n, vbia_n)
      !
      ! Projection on isopycnals
      !
      IF ( ln_pyp ) THEN
         !
         CALL  tlu_vel_proj( uisd_0, visd_0, wisd_0, uisd_n, visd_n, wisd_n )
         CALL  tlu_vel_proj(   unoi,   vnoi,   wnoi,   unoi,   vnoi,   wnoi )
         !    
         IF (ln_tlu_bia) CALL tlu_vel_proj( ubia_n, vbia_n, wbia_n, ubia_n, vbia_n, wbia_n )
         IF (ln_tlu_bia) wbia_n(:,:,1) = 0._wp
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
         vnoi = vnoi * v_wht * vmask
         !
         CALL build_vmodes(w_wht, u_wht, v_wht)
         wnoi = wnoi + w_wht * wmask
         
         !          
         IF ( ln_pyp ) CALL tlu_tnsf_isd  
         !
      END IF   
      ! 
      ! Rescale the noise by some manipulation
      ! dBt = C * (      sum_{k} \phi_{k} \xi_{k} )
      !   a = C * ( dt * sum_{k} \phi_{k}\phi_{k} )
      !
      IF (ln_tlu_nke) CALL tlu_nke( kt )
      IF (ln_tlu_nke) CALL tlu_amp( kt )
      !
      CALL tlu_norm_mod ( kt, .TRUE. )
      !
      IF(  (kt == nit000 + dt_delay) .AND. (.NOT. ln_pyp) ) THEN 
         !
         ! Print variance tensor (a) at time t0 i.e. stationary
         !
         CALL iom_put( 'var_uu_once', var_ten(:,:,:,1) )  
         CALL iom_put( 'var_vv_once', var_ten(:,:,:,2) )  
         CALL iom_put( 'var_ww_once', var_ten(:,:,:,3) )  
         CALL iom_put( 'var_uv_once', var_ten(:,:,:,4) )  
         CALL iom_put( 'var_uw_once', var_ten(:,:,:,5) )  
         CALL iom_put( 'var_vw_once', var_ten(:,:,:,6) ) 
         !
         ! Print Ito-Stokes Drift at time t0 i.e. stationary
         !
         CALL iom_put( 'U_isd_once', uisd_n )  
         CALL iom_put( 'V_isd_once', visd_n )  
         CALL iom_put( 'W_isd_once', wisd_n )   
         !
      ELSE IF ( ln_pyp ) THEN
         !
         ! Print variance tensor (a) at each time
         !
         CALL iom_put( 'var_uu', var_ten(:,:,:,1) )  
         CALL iom_put( 'var_vv', var_ten(:,:,:,2) )  
         CALL iom_put( 'var_ww', var_ten(:,:,:,3) )  
         CALL iom_put( 'var_uv', var_ten(:,:,:,4) )  
         CALL iom_put( 'var_uw', var_ten(:,:,:,5) )  
         CALL iom_put( 'var_vw', var_ten(:,:,:,6) ) 
         !
         ! Print Ito-Stokes Drift
         !
         CALL iom_put( 'U_isd', uisd_n )  
         CALL iom_put( 'V_isd', visd_n )  
         CALL iom_put( 'W_isd', wisd_n )   
         !
      END IF
      !
      ! Print noise fields
      !
      CALL iom_put( "unoi", unoi )
      CALL iom_put( "vnoi", vnoi )
      CALL iom_put( "wnoi", wnoi )
      !
      IF (ln_tlu_bia) THEN
         !
         ! Print bias fields
         !
         CALL iom_put( "ubia_n", ubia_n )
         CALL iom_put( "vbia_n", vbia_n )
         CALL iom_put( "wbia_n", wbia_n )
         !
         ubia_n = 0.75 * ubia_n
         vbia_n = 0.75 * vbia_n
         wbia_n = 0.75 * wbia_n
         !
      END IF
      !
      IF( ln_timing )  CALL timing_stop('tlu_noi')   ! [NEMO] check
      !
   END SUBROUTINE tlu_fields
   ! [tlu_fields]   

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_arrays_init ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Initialize to zero all "now" fields
   !!
   !> @param[in]     kt : ocean integer time step
   !! 
   !> @result        All fields that modify the dynamics are set to zero
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_arrays_init
   ! [tlu_arrays_init]
   SUBROUTINE tlu_arrays_init( kt )
      INTEGER, INTENT(in   ) :: kt         ! ocean time-step index
      !
      unoi = 0._wp
      vnoi = 0._wp
      wnoi = 0._wp
      !
      IF (kt == nit000 + dt_delay)  uisd_n = 0._wp
      IF (kt == nit000 + dt_delay)  visd_n = 0._wp
      IF (kt == nit000 + dt_delay)  wisd_n = 0._wp
      !
      IF (kt == nit000 + dt_delay) var_ten = 0._wp
      !
      IF (ln_tlu_bia)    ubia_n = 0._wp
      IF (ln_tlu_bia)    vbia_n = 0._wp
      IF (ln_tlu_bia)    wbia_n = 0._wp
      !
   END SUBROUTINE tlu_arrays_init
   ! [tlu_arrays_init]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tnsf_isd ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         transfers the now stokes drift to the initial field
   !!
   !> @result        Now stokes drift transferred to the initial field
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_tnsf_isd
   ! [tlu_tnsf_isd]
   SUBROUTINE tlu_tnsf_isd
      !
      uisd_0 = uisd_n
      visd_0 = visd_n
      wisd_0 = wisd_n
      !
   END SUBROUTINE tlu_tnsf_isd
   ! [tlu_tnsf_isd]


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
         CALL lbc_lnk_multi( 'tlu_noi', unoi , 'U', -1., vnoi , 'V', -1. )
         !
      END IF
      !
   END SUBROUTINE tlu_noi
   ! [tlu_noi]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_var ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2023 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes the variance
   !!
   !> @details       Computes the variance starting from the noise modes in each 
   !!                direction, x_mod, y_mod and z_mod, times the timestep dt
   !!                needed for units consistency
   !!                \f{align*}{
   !!                  a(x) = \sum_{k=1}^{N} \phi_{k}(x)\phi^{T}_{k} \Delta t
   !!                \f}
   !!                where $\phi_{k}(x)$ is represented by x_mod, y_mod, z_mod,
   !!                $\Delta t$  is represented by dt, $N$ is nn_nmod.
   !!  
   !> @note          It may or may not communicate the halo to other processors
   !!                depending on the value of lbc_comm
   !!
   !> @param[in]     x_mod 
   !> @param[in]     y_mod 
   !> @param[in]     z_mod 
   !> @param[in]     nn_nmod 
   !> @param[in]     dt
   !> @param[in]     lbc_comm 
   !! 
   !! @result        Adds to the variance field the state given by the input modes
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_var
   ! [tlu_var]
   SUBROUTINE tlu_var( x_mod, y_mod, z_mod, nn_nmod, dt, lbc_comm )
      REAL(wp),                   INTENT(in  ) :: dt        ! Time interval
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: x_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: y_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: z_mod

      INTEGER,                    INTENT(in  ) :: nn_nmod
      LOGICAL,                    INTENT(in  ) :: lbc_comm

      INTEGER                                  :: jm, m_idx

      INTEGER, PARAMETER                       :: ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      DO jm = 1, nn_nmod
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
         !
         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia11) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia11)                              &
         !
         &                                               + ( x_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2   & 
         &                                                 + x_mod( 1:jpi-2, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia22) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia22)                              &
         !
         &                                               + ( y_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2   & 
         &                                                 + y_mod( 2:jpi-1, 1:jpj-2, m_idx + 1 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia33) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia33)                              &
         !
         &                                               + ( z_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 1 )**2   & 
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia12) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia12)                             &
         !
         &                                               + ( x_mod( 2:jpi-1, 3:jpj-1, m_idx + 1 : m_idx + jpk     )      &
         &                                                 + x_mod( 2:jpi-1, 2:jpj-2, m_idx + 1 : m_idx + jpk     ) )    &
         &                                               * ( y_mod( 3:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )      &
         &                                                 + y_mod( 2:jpi-2, 2:jpj-1, m_idx + 1 : m_idx + jpk     ) )    &
         &                                               * 0.25_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia13) = var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia13)                             &
         !
         &                                               + ( x_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 2 )      &
         &                                                 + x_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * ( z_mod( 3:jpi  , 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 )      &
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * 0.25_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia23) = var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia23)                             &
         !
         &                                               + ( y_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 2 )      &
         &                                                 + y_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * ( z_mod( 2:jpi-1, 3:jpj  , m_idx + 2 : m_idx + jpk - 1 )      &
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * 0.25_wp * dt
         !
      END DO
      !
      ! Mask the components
      !
      var_ten(:,:,:,ia11) = var_ten(:,:,:,ia11) * tmask
      var_ten(:,:,:,ia22) = var_ten(:,:,:,ia22) * tmask
      var_ten(:,:,:,ia33) = var_ten(:,:,:,ia33) * tmask
      var_ten(:,:,:,ia12) = var_ten(:,:,:,ia12) * fmask
      var_ten(:,:,:,ia13) = var_ten(:,:,:,ia13) * wumask
      var_ten(:,:,:,ia23) = var_ten(:,:,:,ia23) * wvmask
      !
      ! Lateral boundary condition transfer across nodes
      !
      IF ( lbc_comm ) THEN
         !
         CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,ia11) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia22) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia33) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia12) , 'F', -1.,   &
         &                                  var_ten(:,:,:,ia13) , 'U', -1.,   &
         &                                  var_ten(:,:,:,ia23) , 'V', -1.    )
         !
      END IF
      !
   END SUBROUTINE tlu_var
   ! [tlu_var]


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
      !
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'build_vmodes', wout , 'T', -1.  )
      !
   END SUBROUTINE build_vmodes   
   ! [tlu_wzvcmp]
   

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
   !> @brief         Computes the Ito-Stokes drift
   !!
   !> @details       Computes the Ito-Stokes drift starting from the variance
   !!                direction, x_mod, y_mod and z_mod, times the timestep dt
   !!                needed for units consistency
   !!                \f{align*}{
   !!                  u_{s} = \frac{1}{2} \nabla \cdot a
   !!                \f}
   !!  
   !> @param[in]     var_ten (as a global variable) 
   !! 
   !! @result        Computes Ito-Stokes drift
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_isd
   ! [tlu_isd]
   SUBROUTINE tlu_isd( pcu, pcv, pcw )
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcu                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcv                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcw                   ! before field (so its euler)
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)         ::   flux_var 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3u   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3v   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3w        
      INTEGER                                           ::   ierr                
      !
      INTEGER                                           ::   jk                    ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia21 = ndiffidx(2,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia31 = ndiffidx(3,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia32 = ndiffidx(3,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      ALLOCATE( flux_var(jpi,jpj,jpk,9), stat=ierr )

      ALLOCATE( r1_e1e2e3u(jpi,jpj,jpk), &
      &         r1_e1e2e3v(jpi,jpj,jpk), &
      &         r1_e1e2e3w(jpi,jpj,jpk), stat=ierr )

      flux_var(:,:,:,:) = 0._wp
      flux_var(:,:,:,3) = var_ten(:,:,:,ia31) * 0.5_wp
      flux_var(:,:,:,6) = var_ten(:,:,:,ia32) * 0.5_wp
      flux_var(:,:,:,9) = var_ten(:,:,:,ia33) * 0.5_wp
      !                     
      ! Divergence is by column
      !
      DO jk = 1, jpkm1
         !
         r1_e1e2e3u(:,:,jk) = r1_e1e2u(:,:) * umask(:,:,jk) / e3u_n(:,:,jk) 
         r1_e1e2e3v(:,:,jk) = r1_e1e2v(:,:) * vmask(:,:,jk) / e3v_n(:,:,jk) 
         r1_e1e2e3w(:,:,jk) = r1_e1e2t(:,:) * wmask(:,:,jk) / e3w_n(:,:,jk)
         !
         flux_var(:,:,jk,1) = var_ten(:,:,jk,ia11) * e2t(:,:) * e3t_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,2) = var_ten(:,:,jk,ia21) * e1f(:,:) * e3f_n(:,:,jk) * 0.5_wp
         !
         flux_var(:,:,jk,4) = var_ten(:,:,jk,ia12) * e2f(:,:) * e3f_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,5) = var_ten(:,:,jk,ia22) * e1t(:,:) * e3t_n(:,:,jk) * 0.5_wp
         !
         flux_var(:,:,jk,7) = var_ten(:,:,jk,ia13) * e2u(:,:) * e3uw_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,8) = var_ten(:,:,jk,ia23) * e1v(:,:) * e3vw_n(:,:,jk) * 0.5_wp
         !
         !
         !
         pcu(2:jpi -1, 2:jpj -1, jk   ) =        pcu( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3u( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_11)
         &                              * ( flux_var( 3:jpi   , 2:jpj -1, jk   , 1)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk   , 1)    &
         !  d_y(a_21)
         &                                + flux_var( 2:jpi -1, 2:jpj -1, jk   , 2)    &
         &                                - flux_var( 2:jpi -1, 1:jpj -2, jk   , 2) )  &
         !  d_z(a_31)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 3)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 3) )  &
         &                               /      e3u_n( 2:jpi -1, 2:jpj -1, jk   )
         !
         !
         !
         pcv(2:jpi -1, 2:jpj -1, jk   ) =        pcv( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3v( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_12)
         &                              * ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 4)    &
         &                                - flux_var( 1:jpi -2, 2:jpj -1, jk   , 4)    &
         !  d_y(a_22)
         &                                + flux_var( 2:jpi -1, 3:jpj -1, jk   , 5)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -2, jk   , 5) )  &
         !  d_z(a_32)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 6)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 6) )  &
         &                              /      e3v_n( 2:jpi -1, 2:jpj -1, jk   )
         !
         !
         !
         pcw(2:jpi -1, 2:jpj -1, jk   ) =        pcw( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3w( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_13)
         &                              * ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 7)    &
         &                                - flux_var( 1:jpi -2, 2:jpj -1, jk   , 7)    &
         !  d_y(a_23)
         &                                + flux_var( 2:jpi -1, 2:jpj -1, jk   , 8)    &
         &                                - flux_var( 2:jpi -1, 1:jpj -2, jk   , 8) )

         pcw(2:jpi -1, 2:jpj -1, jk +1) =        pcw( 2:jpi -1, 2:jpj -1, jk +1)       &
         !  d_z(a_33)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 9)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 9) )  &
         &                              *      wmask( 2:jpi -1, 2:jpj -1, jk +1)       &
         &                              /      e3w_n( 2:jpi -1, 2:jpj -1, jk +1)
         !
      END DO
      !
      ! Top boundary condition
      !
      pcw(2:jpi -1, 2:jpj -1, 1) =      pcw( 2:jpi -1, 2:jpj -1, 1)                &
      !  d_z(a_33)|_1
      &                          - flux_var( 2:jpi -1, 2:jpj -1, 1, 9)             &
      &                          *    wmask( 2:jpi -1, 2:jpj -1, 1)                &
      &                          /    e3w_n( 2:jpi -1, 2:jpj -1, 1)
      !
      ! Bottom boundary condition
      !
      pcw(:,:,jpk) = 0._wp
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_isd', pcu , 'U', -1., pcv , 'V', -1., pcw , 'T', -1.  )
      !
   END SUBROUTINE tlu_isd
   ! [tlu_isd]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_norm_mod ***
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
   ! @snippet this tlu_norm_mod
   ! [tlu_norm_mod]
   SUBROUTINE tlu_norm_mod ( kt , prt_noi )
      USE daymod
      INTEGER                           , INTENT(in   ) ::   kt         ! ocean time-step index
      LOGICAL                           , INTENT(in   ) ::   prt_noi    ! ocean time-step index

      REAL(wp), DIMENSION(jpi,jpj,jpk)                  :: zwke_noi, zwke_vel
      REAL(wp), DIMENSION(mppsize)                      :: zw1, zw2
      REAL(wp), DIMENSION(12,mppsize)                   :: zw3
      REAL(wp)                                          :: tke_vel_sub
      REAL(wp)                                          :: tke_noi_sub
      REAL(wp)                                          :: tke_noi
      REAL(wp)                                          :: tke_vel
      REAL(wp), DIMENSION(12) ::    zmaxl,    zmaxg
      REAL(wp), DIMENSION(12) ::   zmeanl,   zmeang

      INTEGER                                           :: i, ierr
      !!------------------------------------------------------------------------
      !

      zmeanl = 0._wp
      zmaxl = 0._wp


!      zwke_noi = 0._wp
!      zwke_noi = ( unoi**2 ) * (spread(e1e2u,3,jpk)*e3u_n) + &
!               & ( vnoi**2 ) * (spread(e1e2v,3,jpk)*e3v_n) + &
!               & ( wnoi**2 ) * (spread(e1e2t,3,jpk)*e3w_n) 

!      tke_noi_sub = SUM( zwke_noi(2:jpim1, 2:jpjm1, 1:jpkm1) )

!      zwke_vel = 0._wp
!      zwke_vel = ( ubia_n**2 ) * (spread(e1e2u,3,jpk)*e3u_n) + &
!               & ( vbia_n**2 ) * (spread(e1e2v,3,jpk)*e3v_n) + &
!               & ( wbia_n**2 ) * (spread(e1e2t,3,jpk)*e3w_n)  

!      tke_vel_sub = SUM( zwke_vel(2:jpim1, 2:jpjm1, 1:jpkm1) )

      zmeanl(1) = SUM( unoi(2:jpim1, 2:jpjm1, 1:jpkm1) )
      zmeanl(2) = SUM( vnoi(2:jpim1, 2:jpjm1, 1:jpkm1) )
      zmeanl(3) = SUM( wnoi(2:jpim1, 2:jpjm1, 1:jpkm1) )

      zmeanl(4) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,1) )
      zmeanl(5) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,2) )
      zmeanl(6) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,3) )

      zmeanl(7) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,4) )
      zmeanl(8) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,5) )
      zmeanl(9) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,6) )

      IF (ln_tlu_bia) THEN 
         zmeanl(10) = SUM( ubia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
         zmeanl(11) = SUM( vbia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
         zmeanl(12) = SUM( wbia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
      ELSE
         zmeanl(10) = 0._wp 
         zmeanl(11) = 0._wp 
         zmeanl(12) = 0._wp
      END IF

      call MPI_BARRIER(mpi_comm_oce, ierr)
!      call MPI_ALLGATHER( tke_noi_sub, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
!      call MPI_ALLGATHER( tke_vel_sub, 1, mpi_double_precision, zw2, 1, mpi_double_precision, mpi_comm_oce, ierr)
      
      call MPI_ALLREDUCE(zmeanl, zmeang, 12, mpi_double_precision, MPI_SUM, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)

      zmeang = zmeang / ( mppsize * jpim1 * jpjm1 * jpkm1 )
!      tke_noi = SUM(zw1)
!      tke_vel = SUM(zw2)

      KE_ratio_n = 0.2 !* SQRT( tke_vel / tke_noi )

      unoi = unoi! * KE_ratio_n
      vnoi = vnoi! * KE_ratio_n
      wnoi = wnoi! * KE_ratio_n

      var_ten(:,:,:,:) = var_ten(:,:,:,:)! * (KE_ratio_n**2) / (KE_ratio_b**2)
!      KE_ratio_b = KE_ratio_n
      IF (prt_noi) THEN
        zmaxl(1) = maxval(unoi); zmaxl(2) = maxval(vnoi); zmaxl(3) = maxval(wnoi);
        zmaxl(4) = maxval(var_ten(:,:,:,1)); zmaxl(5) = maxval(var_ten(:,:,:,2)); zmaxl(6) = maxval(var_ten(:,:,:,3));
        zmaxl(7) = maxval(var_ten(:,:,:,4)); zmaxl(8) = maxval(var_ten(:,:,:,5)); zmaxl(9) = maxval(var_ten(:,:,:,6));
        IF (ln_tlu_bia) THEN 
           zmaxl(10) = maxval(ubia_n); zmaxl(11) = maxval(vbia_n); zmaxl(12) = maxval(wbia_n);
        END IF
    
        call MPI_BARRIER(mpi_comm_oce, ierr)
        call MPI_REDUCE(zmaxl, zmaxg, 12, mpi_double_precision, MPI_MAX, 0, mpi_comm_oce, ierr)

        IF (lwp .AND. mod(kt,500) .eq. 0) THEN
            print *, ""
            print '(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i3.3)', &
                  &    '======================>> time-step =', kt,'      DATE Y/M/D = ', nyear, '/', nmonth, '/', nday, '      nday_year = ', nday_year
            print *, '--------------------------------------------------------------------------------------------- Maximum values'
            print *, '     Max noise: ',  zmaxg( 1), zmaxg( 2), zmaxg( 3)
            print *, ' Max trace var: ',  zmaxg( 4), zmaxg( 5), zmaxg( 6)
            print *, ' Max cross var: ',  zmaxg( 7), zmaxg( 8), zmaxg( 9)
            print *, '      Max bias: ',  zmaxg(10), zmaxg(11), zmaxg(12)
            print *, '------------------------------------------------------------------------------------------------ Mean values'
            print *, '    Mean noise: ', zmeang( 1),zmeang( 2),zmeang( 3)
            print *, 'Mean trace var: ', zmeang( 4),zmeang( 5),zmeang( 6)
            print *, 'Mean cross var: ', zmeang( 7),zmeang( 8),zmeang( 9)
            print *, '     Mean bias: ', zmeang(10),zmeang(11),zmeang(12)
            print *, '------------------------------------------------------------------------------------------ Ratio of energies'
            print *, 'KEbias/KEnoise: '!, SQRT( tke_vel / tke_noi )
            print *, ""
        ENDIF
      ENDIF

   END SUBROUTINE tlu_norm_mod
   ! [tlu_norm_mod]

END MODULE tlunoi

!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_step
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
!> @brief         Transport under Location Uncertainty module.
!! 
!! 
!! @par           Procedure specifics      
!> @details       Defines variables and dependencies for the stochastic parametrization of 
!!                Navier-Stokes equations. 
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
MODULE tlustep
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE iom
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
!  USE mockup_nemo      ! [NEMO] to be replaced by relevant actual NEMO modules:  wrk_nemo
!  USE wrk_nemo         ! Memory allocation
   USE tlu
   USE tlusvd
   USE tluhdyn          ! Horizontal dynamics contribution
   USE tlubias
   USE tluzdyn          ! Vertical dynamics contribution
   USE tluhzdyn         ! Mixed Horizontal-Vertical dynamics contribution
   USE tlusbc           ! Surface boundary conditions
   USE tlustp           ! Stochastic pressure
   USE tlucmp           ! Stochastic compressibility
   USE tlutraadv        ! Tracer advection
   USE tlutradiff       ! Tracer diffusion 
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_bcdyn     ! Called by step.f90
   PUBLIC tlu_tradyn    ! Called by step.f90
   ! [public_sub]  
   !

CONTAINS


   !!--------------------------------------------------------------------------------------------------------------------------- 
   !!            ***  ROUTINE tlu_bcdyn ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Computation of location uncertainty terms in the BAROCLINIC time step 
   !!
   !> @par           Code specifics 
   !> @details       The code makes a sequence of calls to the single routines that implement the dynamic parts
   !! @snippet       this Advection  
   !!                @ref tlu_hadv_noi for the advection of horizontal components, 
   !!                @ref tlu_zadv_noi for the advection of vertical components,
   !! @snippet       this Stochastic drift
   !!                @ref tlu_hhdrift for the drift of horizontal components, 
   !!                @ref tlu_zzdrift for the drift of vertical components,
   !! @snippet       this Stochastic diffusion
   !!                @ref tlu_hhdiff for the diffusion of horizontal components, 
   !!                @ref tlu_zzdiff for the diffusion of vertical components,
   !! @snippet       this Stochastic surface boundary conditions
   !!                @ref tlu_zadv_sbc for the advection of horizontal components with wind influence,
   !!                @ref    tlu_zzsbc for the drift-diffusion of vertical components with wind influence, 
   !! @snippet       this Stochastic pressure model
   !!                @ref tlu_stpdyn for the stochastic pressure dynamics
   !! 
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !! 
   !! @result        Updated (ua,va) with all the Location Uncertainty terms
   !!  
   !!
   !!---------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_bcdyn 
   ! [tlu_bcdyn]
   SUBROUTINE tlu_bcdyn(kt)
      INTEGER,                          INTENT(in)     ::  kt      
      REAL(wp), DIMENSION(jpi,jpj,jpk)                 ::  unSC, ubSC   
      REAL(wp), DIMENSION(jpi,jpj,jpk)                 ::  vnSC, vbSC    
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_bcdyn')   ! [NEMO] check

      unSC = un * r1_SC
      vnSC = vn * r1_SC
      ubSC = ub * r1_SC
      vbSC = vb * r1_SC

      ! 
      ! [Advection]
      CALL tlu_hadv_noi( kt, 'U', 1 , unSC, ua )
      CALL tlu_zadv_noi( kt, 'U', 1 , unSC, ua )
      !
      CALL tlu_hadv_noi( kt, 'V', 2 , vnSC, va )
      CALL tlu_zadv_noi( kt, 'V', 2 , vnSC, va )
      ! [Advection]
 
      ! [Bias advection]
      IF (ln_tlu_bia) THEN
         CALL tlu_adv_bias( kt, 'U', 1 , un, ua, biaSIGN )
         !
         CALL tlu_adv_bias( kt, 'V', 2 , vn, va, biaSIGN )
      END IF
      ! [Bias advection]

      ! [Stochastic drift]
      CALL  tlu_hhdrift( kt, 'U', 1 , unSC, ua )
      CALL  tlu_hzdrift( kt, 'U', 1 , unSC, ua )
      CALL  tlu_zzdrift( kt, 'U', 1 , unSC, ua )
      !
      CALL  tlu_hhdrift( kt, 'V', 2 , vnSC, va )
      CALL  tlu_hzdrift( kt, 'V', 2 , vnSC, va )
      CALL  tlu_zzdrift( kt, 'V', 2 , vnSC, va )
      ! [Stochastic drift]

      ! [Stochastic diffusion]  ! b-fields because of euler scheme
      CALL   tlu_hhdiff( kt, 'U', 1 , ubSC, ua )
      CALL   tlu_hzdiff( kt, 'U', 1 , ubSC, ua )
      CALL   tlu_zzdiff( kt, 'U', 1 , ubSC, ua )
      !
      CALL   tlu_hhdiff( kt, 'V', 2 , vbSC, va )
      CALL   tlu_hzdiff( kt, 'V', 2 , vbSC, va )
      CALL   tlu_zzdiff( kt, 'V', 2 , vbSC, va )
      ! [Stochastic diffusion]

      ! [Stochastic surface boundary conditions]
      CALL tlu_zadv_sbc( kt, 'U', 1 , un, ua )
      CALL    tlu_zzsbc( kt, 'U', 1 , ub, un, ua )
      !
      CALL tlu_zadv_sbc( kt, 'V', 2 , vn, va )
      CALL    tlu_zzsbc( kt, 'V', 2 , vb, vn, va )
      ! [Stochastic surface boundary conditions]

      ! [Stochastic compressibility]
      CALL  tlu_hhcmp( kt, 'U', 1 , unSC, ua )
      CALL  tlu_hzcmp( kt, 'U', 1 , unSC, ua )
      CALL  tlu_zzcmp( kt, 'U', 1 , unSC, ua )
      CALL tlu_cmpsbc( kt, 'U', 1 , unSC, ua )
      !
      CALL  tlu_hhcmp( kt, 'V', 2 , vnSC, va )
      CALL  tlu_hzcmp( kt, 'V', 2 , vnSC, va )
      CALL  tlu_zzcmp( kt, 'V', 2 , vnSC, va )
      CALL tlu_cmpsbc( kt, 'V', 2 , vnSC, va )
      ! [Stochastic compressibility]

      ! [Stochastic pressure model]
      !    Ge - Geostrophic Balance noise-stp
      !    QG - Quasi-Geostrophic noise (not implemented)
      CALL   tlu_stpdyn( kt, 'Ge' )
      ! [Stochastic pressure model]


      !
      IF( ln_timing ) CALL timing_stop('tlu_dyn')   ! [NEMO] check
      !
   END SUBROUTINE tlu_bcdyn
   ! [tlu_bcdyn]

   !!--------------------------------------------------------------------------------------------------------------------------- 
   !!            ***  ROUTINE tlu_tradyn ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Computation of location uncertainty terms in the BAROCLINIC time step 
   !!
   !> @par           Code specifics 
   !> @details       The code makes a sequence of calls to the single routines that implement the dynamic parts
   !! @snippet       this Tracer Advection  
   !!                @ref tlu_tra_adv
   !! @snippet       this Bias Tracer advection
   !!                @ref tlu_tra_adv
   !! @snippet       this Tracer Stochastic diffusion
   !!                @ref tlu_trahhdiff for the diffusion of horizontal components,
   !!                @ref tlu_trahzdiff for the diffusion of horizontal components, 
   !!                @ref tlu_trazzdiff for the diffusion of horizontal components,  
   !! 
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !! 
   !! @result        Updated (ua,va) with all the Location Uncertainty terms
   !!  
   !!
   !!---------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_tradyn 
   ! [tlu_tradyn]
   SUBROUTINE tlu_tradyn(kt)
      INTEGER, INTENT(in)                     ::   kt       
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_tradyn')   ! [NEMO] check

      ! 
      ! [Tracer Advection]
      CALL  tlu_tra_adv( kt, unoi, vnoi, wnoi)
      CALL  tlu_tra_adv( kt, -1._wp * uisd_n, &
                       &     -1._wp * visd_n, &
                       &     -1._wp * wisd_n  )
      ! [Tracer Advection]
 
      ! [Bias Tracer advection]
      IF (ln_tlu_bia) THEN
         CALL  tlu_tra_adv( kt, -1._wp * biaSIGN * ubia, &
                          &     -1._wp * biaSIGN * vbia, &
                          &     -1._wp * biaSIGN * wbia  )
      END IF
      ! [Bias Tracer advection]

      ! [Tracer Stochastic diffusion] ! b-fields because of euler scheme
      CALL tlu_trahhdiff( kt, tsb, tsa )
!     CALL tlu_trahzdiff( kt, tsb, tsa )
      CALL tlu_trazzdiff( kt, tsb, tsa )
      ! [Tracer Stochastic diffusion]
      !
      IF( ln_timing ) CALL timing_stop('tlu_tradyn')   ! [NEMO] check
      !
   END SUBROUTINE tlu_tradyn
   ! [tlu_tradyn]

END MODULE tlustep




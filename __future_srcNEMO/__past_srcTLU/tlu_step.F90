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

   ! General definition 
   USE tlu

   ! Noise Definition 
   USE tlunoi

   ! Momentum dynamics
   USE tluhdifDYN
   USE tluzdifDYN

!   USE tlusbcDYN           ! Surface boundary conditions

   ! Stochastic compressibility model
!   USE tlucmp           ! Stochastic compressibility

   ! Stochastic Pressure model
   USE tlustp           ! Stochastic pressure

   ! Tracer dynamics
   USE tluhdifTRA       ! Tracer horizontal diffusion
   USE tluzdfTRA        ! Tracer vertical diffusion

!   USE tluadv_ubs
!   USE tluadv_cen2
   ! 
 
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
   !! @snippet       this Ito-Stokes Advection
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
   SUBROUTINE tlu_bcdyn(kt, opt)
      USE dynadv           ! For n_dynadv
      !
      INTEGER,          INTENT(in   )  ::  kt      
      CHARACTER(len=3), INTENT(in   )  ::  opt    
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_bcdyn')   ! [NEMO] check

      SELECT CASE( opt )     

         !> @par
         ! [Stochastic diffusion] ! b-fields because of euler scheme            
         CASE ( 'hdf' )
             CALL   tlu_hhdiff( kt, 'U', 1 , ub, ua )
             CALL   tlu_hhdiff( kt, 'V', 2 , vb, va )
                   
         CASE ( 'cdf' )    
!!             CALL   tlu_hzdiff( kt, 'U', 1 , ub, ua, 1._wp )
!!             CALL   tlu_hzdiff( kt, 'V', 2 , vb, va, 1._wp ) 

         CASE ( 'zdf' )
              CALL   tlu_dyn_zdf ( kt )
         ! [Stochastic diffusion]
         !



         ! [Stochastic surface boundary conditions]
         CASE ( 'sbc' )
         ! [Stochastic surface boundary conditions]

         ! [Stochastic compressibility]
         CASE ( 'cmp' )
         ! [Stochastic compressibility]

         ! [Stochastic pressure model]
         CASE ( 'stp' )
            !    Ge - Geostrophic Balance noise-stp
            !    QG - Quasi-Geostrophic noise (not implemented)
            CALL   tlu_stpdyn( kt, 'Ge' )

         ! [Stochastic pressure model]

         CASE DEFAULT 
            CALL ctl_stop('STOP','tlu_tradyn: wrong value for opt'  )
      END SELECT
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
   !!                @ref tlu_trahzdiff for the diffusion of mixed components, 
   !!                @ref tlu_trazzdiff for the diffusion of vertical components,  
   !! 
   !! 
   !> @param[in]      kt: ocean time-step integer index
   !> @param[in]     opt: logical switch to choose different components in different points of step.F90
   !! 
   !! @result        Updated (ua,va) with all the Location Uncertainty terms
   !!  
   !!
   !!---------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_tradyn 
   ! [tlu_tradyn]
   SUBROUTINE tlu_tradyn(kt, opt)
!      USE tludia, ONLY :  Rt_NoD_hh,  Rt_NoD_zz
      INTEGER,          INTENT(in   )  ::  kt      
      CHARACTER(len=3), INTENT(in   )  ::  opt 
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_tradyn')   ! [NEMO] check

      SELECT CASE( opt )     

         !> @par
         ! [Tracer Stochastic diffusion] ! b-fields because of euler scheme            
         CASE ( 'hdf' )
             CALL tlu_trahhdiff( kt, tsb, tsa, (/ Sc_T**2, Sc_S**2 /) ) ! (/ Rt_NoD_hh(1), Rt_NoD_hh(2) /) )
      
         CASE ( 'cdf' )
             CALL tlu_trahzdiff( kt, tsb, tsa, (/ Sc_T**2, Sc_S**2 /) ) !  (/ Rt_NoD_hh(1), Rt_NoD_hh(2) /) )
      
         CASE ( 'zdf' )
             CALL tlu_tra_zdf( kt, (/ Sc_T**2, Sc_S**2 /) ) !   (/ Rt_NoD_zz(1), Rt_NoD_zz(2) /) )
      
         ! [Tracer Stochastic diffusion]
         !
         CASE DEFAULT 
            CALL ctl_stop('STOP','tlu_tradyn: wrong value for opt'  )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_tradyn')   ! [NEMO] check
      !
   END SUBROUTINE tlu_tradyn
   ! [tlu_tradyn]

END MODULE tlustep




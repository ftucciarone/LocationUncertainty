!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlusbc
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2017 -  ( F. TUCCIARONE   )  Original code and documentation <br>
!! 
!> @brief         Transport under Location Uncertainty dynamics components.
!! 
!! @par           Procedure specifics      
!> @details       Implements the modification to the standard NEMO dynamics. The transport 
!!                operator implemented is in its incompressible form, that is  
!!                
!!          
!!
!> @par           Code specifics
!!                The public routine @ref tlu_dyn performs the following calls to compute all the terms:
!!                - The advection term @f$ \mathcal{A}_{\boldsymbol{u}} @f$ is computed in @ref tlu_adv_noi_cmp
!!                - The coriolis term @f$ \mathcal{C}_{\boldsymbol{u}} @f$ is computed in @ref tlu_cor_noi
!!                - The drift-diffusion term @f$ \mathcal{A}_{\boldsymbol{u}} @f$ is computed in several rouines to ease 
!!                the calculations. Splitting the tensor @f$\boldsymbol{a}@f$ into the two components
!!               
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
!!
!> @par           Other modules dependencies
!> @snippet       this mod_dep
!!
!> @par           Public routines
!> @snippet       this public_sub
!!
!> @todo          read namelist from _cfg, not only _ref
!! @todo          write error messages to appropriate unit
!!
!> @warning       Use of `dynvor` sketchy
!! @warning       Several errors in the positioning of the scale factors  
!! @warning       No implementation of the @f$ \boldsymbol{u}\nabla\cdot\nabla\cdot\boldsymbol{a} @f$ term
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tlusbcTRA
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE sbc_oce          ! surface boundary condition: ocean
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE tlu              ! Initialization of stochastic structures
   ! [mod_dep]
   !
   IMPLICIT NONE   
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_traadv_sbc  ! Called by tlu_tra_adv_noi in tlu_advTRA.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS


   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_traadv_sbc ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computation of the advection of velocity due to the noise
   !!
   !> @par           Code specifics 
   !> @details       The advection term due to the noise velocity is written explicitly as:
   !!
   !> @param[in]     kt: ocean time-step integer indexi
   !> @param[in]     kcmp: Define velocity type for calculation
   !> @param[in]     cdtype: =U, V or W (component indicator)
   !> @param[inout]  pcn, pca:  modified advection term for this component
   !! 
   !! @result        Update (ua,va) with the noise advection term trend
   !!               
   !!--------------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_traadv_sbc 
   ! [tlu_traadv_sbc]
   SUBROUTINE tlu_traadv_sbc( kt, idtrc, surface_layer)
      INTEGER                     , INTENT(in   ) ::   kt                     ! ocean time-step index
      INTEGER                     , INTENT(in   ) ::   idtrc                  ! Define velocity type for calculation
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   surface_layer          ! modified advection term for this component
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_traadv_sbc')   ! [NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_traadv_sbc : surface boundary condition tracer number', idtrc, '  (1 = temp, 2 = salt)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      SELECT CASE( idtrc )     !Select velocity type to interpolate sigma dbt in to corresponding mesh point for advection calculation
         !> @par               U Velocity
         CASE ( jp_tem )
            surface_layer = 0._wp
         !> @par               V Velocity
         CASE ( jp_sal )
            surface_layer = 0._wp
         CASE DEFAULT
            CALL ctl_stop('STOP','tlu_traadv_sbc: wrong value for idtrc' )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_traadv_sbc')   ! [NEMO] check
      !
   END SUBROUTINE tlu_traadv_sbc
   ! [tlu_traadv_sbc]
 

END MODULE tlusbcTRA


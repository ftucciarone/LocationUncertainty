!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlutrasbc
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
!> @brief         Transport under Location Uncertainty dynamics components.
!! 
!! @par           Procedure specifics      
!> @details       Implements the modification to the standard NEMO dynamics. The transport 
!!                operator implemented is in its incompressible form, that is  
!!          
!!
!> @par           Code specifics
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
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tlutrasbc
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE tlu              ! Initialization of stochastic structures
   ! [mod_dep]
   !
   IMPLICIT NONE   
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_tradflux      ! Called by 
   ! [public_sub]
   !

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !
CONTAINS



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tradflux ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes the drift-diffusion term @f$ \dfrac{1}{2} \nabla\cdot\nabla\cdot\left[( \boldsymbol{a} - 
   !!                \text{Dg}(\boldsymbol{a}) )\boldsymbol{u}\right]@f$
   !!
   !> @details
   !!                \f{align*}{
   !!                \mathcal{D}_{\boldsymbol{u}} & = - \dfrac{1}{2}\dfrac{1}{e_{1}e_{2}e_{3}}\dfrac{\partial}{\partial x}
   !!                        \left[\dfrac{1}{e_{1}}
   !!                        \dfrac{\partial \left(e_{1}e_{3}a_{21} u\right) }{\partial y}\right]
   !!                      - \dfrac{1}{2}\dfrac{1}{e_{1}e_{2}e_{3}}\dfrac{\partial}{\partial y}\left[\dfrac{1}{e_{2}}
   !!                        \dfrac{\partial \left(e_{2}e_{3}a_{12} u\right) }{\partial x}\right]
   !!                \f} 
   !!                
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     kcmp: Define velocity type for calculation
   !> @param[in]     cdtype: Define velocity type   !> @param[in]     pcmask: mask for vertical component
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
   ! @snippet this tlu_tradflux
   ! [tlu_tradflux]
   SUBROUTINE tlu_tradflux( kt, jn, ptb, dflux)
      INTEGER                     , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                     , INTENT(in   ) ::   jn                    ! tracer index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   dflux                 ! flux at the top layer
      !
      INTEGER                                     ::   ji, jj, jk            ! dummy loop indices
      !
      INTEGER, PARAMETER                          ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_tradflux')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_tradflux : '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      dflux = 0._wp
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
!            dflux(ji,jj) = 0.5_wp * var_ten(ji,jj,1 ,ia33) * ( 0._wp - ptb(ji,jj) ) / ( e3w_n(ji,jj,1) )

!            dflux(ji,jj) = 0.5_wp * 0.001_wp * ( 0._wp - ptb(ji,jj) ) / ( e3w_n(ji,jj,1) )
         END DO
      END DO
      !
   END SUBROUTINE tlu_tradflux
   ! [tlu_tradflux]

END MODULE tlutrasbc

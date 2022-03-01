!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tluhdyn
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
MODULE tlutradiff
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
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
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_trahhdiff      ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_trahzdiff      ! Called by 
   PUBLIC tlu_trazzdiff      ! Called by 
   ! [public_sub]
   !

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trahhdiff ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes the Horizontal components of the Stochastic Diffusion
   !!                 @f$ \dfrac{1}{2} \nabla\cdot\nabla\cdot\left[\boldsymbol{a}_{H} - 
   !!                 \boldsymbol{u}_{H}\right]@f$
   !!
   !> @details
   !!                \f{align*}{
   !!                \mathcal{D}_{\boldsymbol{u}}^{hh} & = - \dfrac{1}{2}\dfrac{1}{e_{1}e_{2}e_{3}}\dfrac{\partial}{\partial x}
   !!                        \left[\dfrac{1}{e_{1}}
   !!                        \dfrac{\partial \left(e_{1}e_{3}a_{21} u\right) }{\partial y}\right]
   !!                      - \dfrac{1}{2}\dfrac{1}{e_{1}e_{2}e_{3}}\dfrac{\partial}{\partial y}\left[\dfrac{1}{e_{2}}
   !!                        \dfrac{\partial \left(e_{2}e_{3}a_{12} u\right) }{\partial x}\right]
   !!                \f} 
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
   ! @snippet this tlu_trahhdiff
   ! [tlu_trahhdiff]
   SUBROUTINE tlu_trahhdiff( kt, ptb, pta)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw0                   ! accumulation array
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var11             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var12             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var21             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var22             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztuu                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztuv                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztvu                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztvv                  ! workspace array

      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_trahhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trahhdiff : Horizontal diffusion components on velocity type '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF


      ALLOCATE( int_var11(jpi,jpj,jpk), &
              & int_var12(jpi,jpj,jpk), &
              & int_var21(jpi,jpj,jpk), &
              & int_var22(jpi,jpj,jpk)  )

      ALLOCATE( ztuu(jpi,jpj,jpk), &
              & ztuv(jpi,jpj,jpk), &
              & ztvu(jpi,jpj,jpk), &
              & ztvv(jpi,jpj,jpk)  )

      ALLOCATE( zw0(jpi,jpj,jpk,2) )



      DO jk = 1, jpkm1              !==  Interpolation of variance tensor  ==!
         DO jj = 2, jpjm1
            DO ji =  2, jpim1
               int_var11(ji,jj,jk) = ( var_ten(ji+1,jj  ,jk,ia11) + var_ten(ji  ,jj  ,jk,ia11) ) * 0.5_wp
               int_var12(ji,jj,jk) = ( var_ten(ji  ,jj  ,jk,ia12) + var_ten(ji  ,jj-1,jk,ia12) ) * 0.5_wp
               int_var21(ji,jj,jk) = ( var_ten(ji  ,jj  ,jk,ia12) + var_ten(ji-1,jj  ,jk,ia12) ) * 0.5_wp
               int_var22(ji,jj,jk) = ( var_ten(ji  ,jj+1,jk,ia22) + var_ten(ji  ,jj  ,jk,ia22) ) * 0.5_wp
            END DO
         END DO
      END DO  
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_trahhdiff', int_var11 , 'U', 1.,   &
                        &                  int_var12 , 'V', 1.,   &
                        &                  int_var21 , 'U', 1.,   &
                        &                  int_var22 , 'V', 1.    )
      !
      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !    
         !                               
         DO jk = 1, jpkm1              !==  First derivative (gradient)  ==!
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ztuu(ji,jj,jk) = int_var11(ji,jj,jk) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztuv(ji,jj,jk) = int_var12(ji,jj,jk) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztvu(ji,jj,jk) = int_var21(ji,jj,jk) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztvv(ji,jj,jk) = int_var22(ji,jj,jk) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
         END DO  
         !
         zw0 = 0._wp
         !
         !
         DO jk = 1, jpkm1              !==  Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !
                  ! Diagonal terms of the diffusion
                  !
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  ztuu(ji,jj,jk) - ztuu(ji-1,jj,jk)     &
                     &                                   + ztvv(ji,jj,jk) - ztvv(ji,jj-1,jk) )   &
                     &                                / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
                  !
                  ! Extra-diagonal terms (need interpolation later)
                  !
                  zw0(ji,jj,jk,jn) = zw0(ji,jj,jk,jn) + (  ztuv(ji,jj,jk) - ztuv(ji-1,jj,jk)     &
                     &                                   + ztvu(ji,jj,jk) - ztvu(ji,jj-1,jk) )   &
                     &                                / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )

               END DO
            END DO
         END DO  
         !                          ! ==================
      END DO                        ! end of tracer loop
      !                             ! ==================

      !
      CALL lbc_lnk( 'tlu_trahhdiff', zw0, 'F', 1. )
      !

      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !    
         !
         DO jk = 1, jpkm1              !==  Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !
                  ! Extra-diagonal terms 
                  !
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  zw0(ji  ,jj  ,jk,jn) + zw0(ji-1,jj  ,jk,jn)     &
                     &                                   + zw0(ji  ,jj-1,jk,jn) + zw0(ji-1,jj-1,jk,jn) ) * 0.25_wp 
                  !
               END DO
            END DO
         END DO  
         !                          ! ==================
      END DO                        ! end of tracer loop
      !                             ! ==================

      DEALLOCATE( int_var11, int_var12, int_var21, int_var22 )

      DEALLOCATE( ztuu, ztuv, ztvu, ztvv )

      DEALLOCATE( zw0 )

   END SUBROUTINE tlu_trahhdiff
   ! [tlu_trahhdiff]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trahzdiff ***
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
   ! @snippet this tlu_trahzdiff
   ! [tlu_trahzdiff]
   SUBROUTINE tlu_trahzdiff( kt, ptb, pta)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   zw0                   ! accumulation array
      REAL(wp)                                             ::   zwb , zwt             ! local scalars (FD bottom top)
      REAL(wp)                                             ::   zwl , zwr             ! local scalars (FD left right)
      REAL(wp)                                             ::   zwbl, zwbr            ! local scalars (FD bottom-left bottom-right)
      REAL(wp)                                             ::   zwtl, zwtr            ! local scalars (FD top-left top-right)
      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   zacni                 ! Average in i of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   zacnj                 ! Average in j of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   zacnik                ! Average in ik of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   zacnjk                ! Average in jk of pcn

      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_trahzdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trahzdiff : Mixed Horizontal-Vertical diffusion components on velocity type '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF



      DO jn = 1, jpts            !==  loop over the tracers  ==!
         
         
         
         
         
          
         
         
         
         
         
         
         
         
                                                           
      END DO

   END SUBROUTINE tlu_trahzdiff
   ! [tlu_trahzdiff]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trazzdiff ***
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
   ! @snippet this tlu_trazzdiff
   ! [tlu_trazzdiff]
   SUBROUTINE tlu_trazzdiff( kt, ptb, pta)
      !
      USE tlutrasbc, ONLY : tlu_tradflux     
      !
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var33             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztww                  ! workspace array
      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_trazzdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trazzdiff : '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF

      ALLOCATE(ztww(jpi,jpj,jpk), int_var33(jpi,jpj,jpk) )

      DO jj = 1, jpj               !==  Interpolation of variance tensor  ==!
         DO ji = 1, jpi
            int_var33(ji,jj,1) = ( var_ten(ji,jj,1  ,ia33) ) * 0.5_wp
         END DO
      END DO

      DO jk = 2, jpkm1             !==  Interpolation of variance tensor  ==!  
         DO jj = 1, jpj
            DO ji = 1, jpi
               int_var33(ji,jj,jk) = ( var_ten(ji,jj,jk  ,ia33) + var_ten(ji,jj,jk-1,ia33) ) * 0.5_wp
            END DO
         END DO
      END DO  

      DO jn = 1, jpts              !==  loop over the tracers  ==!
         !
         ztww = 0._wp   
         CALL tlu_tradflux( kt, jn, ptb(:,:,1,jn), ztww(:,:,1)  )                  
         !
         DO jk = 2, jpk             !==  First derivative (gradient)  ==!
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztww(ji,jj,jk) = int_var33(ji,jj,jk) * ( ptb(ji,jj,jk-1,jn) - ptb(ji,jj,jk,jn) )   &
                      &                                / ( e3w_n(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
         ztww(:,:,jpk) = 0._wp      ! Avoid bottom boudary fluxes
         
         DO jk = 1, jpkm1           !==  Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  ! Vertical term of the diffusion
                  !
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  ztww(ji,jj,jk) - ztww(ji,jj,jk+1) ) * 0.5_wp  &
                     &                                / ( e3t_n(ji,jj,jk) )
                  !
               END DO
            END DO
         END DO

      END DO

      DEALLOCATE( ztww, int_var33 )

   END SUBROUTINE tlu_trazzdiff
   ! [tlu_trazzdiff]

END MODULE tlutradiff

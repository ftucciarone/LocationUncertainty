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
!> @warning       Sketchy use of `dynvor`: why is it needed?
!! @warning       BIGGEST WARN OF THEM ALL: THIS CODE IS NOT TESTED  
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tluhdifDYN
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
   PUBLIC tlu_hhdiff       ! Called by tlu_bcdyn in tlu.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hhdiff ***
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
   ! @snippet this tlu_hhdiff
   ! [tlu_hhdiff]
   SUBROUTINE tlu_hhdiff( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zcnt, zcnf            ! averaged component on t and f grid
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   dxuv, dyuv, dzuv      ! averaged component on t and f grid
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw1, zw2, zw3         ! workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwb, zwa              ! local scalars (FD before after)
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia21 = ndiffidx(2,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(1,3)  ! Mapping the indeces of a

      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_hhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hhdiff : Horizontal diffusion components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))
      ALLOCATE(zw1(jpi,jpj,jpk), zw2(jpi,jpj,jpk), zw3(jpi,jpj,jpk))


      ! Velocity averages
      zcnt = 0._wp     ! For diagonal components (at T-point)
      zcnf = 0._wp     ! For extra-diagonal components (at f-point)

      ! Support vectors
      zw0 = 0._wp
      zw1 = 0._wp
      zw2 = 0._wp
      zw3 = 0._wp

      dxuv = 0._wp
      dyuv = 0._wp
      dzuv = 0._wp

      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 
            !
            ! Computation of the Horizontal Gradient of U
            !
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 2, jpi  
                     dxuv(ji,jj,jk) = ( pcn(ji,jj  ,jk) - pcn(ji-1,jj,jk) ) * tmask(ji,jj,jk) * r1_e1t(ji,jj) 
                     dyuv(ji,jj,jk) = ( pcn(ji,jj+1,jk) - pcn(ji  ,jj,jk) ) * fmask(ji,jj,jk) * r1_e2f(ji,jj)  
                  END DO
               END DO
            END DO
            !
            ! Computation of the vertical derivative of U
            !
            DO jk = 2, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi   
                     dzuv(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) - pcn(ji,jj,jk) ) * wumask(ji,jj,jk) / e3uw_n(ji,jj,jk) 
                  END DO
               END DO
            END DO
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hhdiff', dxuv , 'T', 1.,     &
                              &               dyuv , 'F', 1.      )
            !
            ! Computation of the X-directed diffusive flux of U
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     !
                     ! a_{11}\partial_{x} u is at T-point
                     ! 
                     zw0(ji,jj,jk) = dxuv(ji,jj,jk) * var_ten(ji,jj,jk,ia11) * e2t(ji,jj) * e3t_n(ji,jj,jk)
                     !
                     ! a_{12}\partial_{y} u is at f-point, needs 4-point interpolation
                     !
                     zw1(ji,jj,jk) = ( dyuv(ji  ,jj  ,jk) * var_ten(ji  ,jj  ,jk,ia21) * e2f(ji  ,jj  ) * e3f_n(ji  ,jj  ,jk) +  &
                                   &   dyuv(ji-1,jj  ,jk) * var_ten(ji-1,jj  ,jk,ia21) * e2f(ji-1,jj  ) * e3f_n(ji-1,jj  ,jk) +  &
                                   &   dyuv(ji-1,jj-1,jk) * var_ten(ji-1,jj-1,jk,ia21) * e2f(ji-1,jj-1) * e3f_n(ji-1,jj-1,jk) +  &
                                   &   dyuv(ji  ,jj-1,jk) * var_ten(ji  ,jj-1,jk,ia21) * e2f(ji  ,jj-1) * e3f_n(ji  ,jj-1,jk) )  &
                                   & * 0.25_wp
                     !
                     ! Cumulate on zw0
                     !
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + zw1(ji,jj,jk)
                     !
                     ! a_{13}\partial_{z} u is at uw-point, needs 2-points interpolation in x and 2-points interpolation in z (later)
                     !
                     zw2(ji,jj,jk) = ( dzuv(ji  ,jj  ,jk) * var_ten(ji  ,jj  ,jk,ia21) * e2f(ji  ,jj) * e3f_n(ji  ,jj  ,jk) +  &
                                   &   dzuv(ji-1,jj  ,jk) * var_ten(ji-1,jj  ,jk,ia21) * e2f(ji-1,jj) * e3f_n(ji-1,jj  ,jk) )  &
                                   & * 0.5_wp 
                  END DO
               END DO
            END DO
            !
            ! Surface boundary condition on the shear of U
            !
            zw2(:,:,1) = 0._wp
            !
            ! Finalisation of the X-directed diffusive flux of U with vertical interpolation of zw2
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + ( zw2(ji,jj,jk) + zw2(ji,jj,jk+1) ) * 0.5_wp
                  END DO 
               END DO
            END DO
            !
            ! Reinitialise support vectors
            !  
            zw1 = 0._wp
            zw2 = 0._wp
            zw3 = 0._wp
            !
            ! Computation of the Y-directed diffusive flux of U
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     !
                     ! \partial_{x} u is at T-point, needs 4-points interpolation
                     !
                     zw1(ji,jj,jk) = (    dxuv(ji+1,jj+1,jk) * e1t(ji+1,jj+1) * e3t_n(ji+1,jj+1,jk) +  &
                                   &      dxuv(ji  ,jj+1,jk) * e1t(ji  ,jj+1) * e3t_n(ji  ,jj+1,jk) +  &
                                   &      dxuv(ji  ,jj  ,jk) * e1t(ji  ,jj  ) * e3t_n(ji  ,jj  ,jk) +  &
                                   &      dxuv(ji+1,jj  ,jk) * e1t(ji+1,jj  ) * e3t_n(ji+1,jj  ,jk) )  &
                                   & * var_ten(ji  ,jj  ,jk,ia21) * 0.25_wp
                     !
                     ! a_{22} is at T-point, needs 4-points interpolation
                     ! 
                     zw2(ji,jj,jk) = ( var_ten(ji+1,jj+1,jk,ia22) + var_ten(ji  ,jj+1,jk,ia22) + &
                                   &   var_ten(ji  ,jj  ,jk,ia22) + var_ten(ji+1,jj  ,jk,ia22) ) & 
                                   & *    dyuv(ji,jj,jk) * e1f(ji,jj) * e3f_n(ji,jj,jk) * 0.25_wp
                     !
                     ! Cumulate on zw1
                     !
                     zw1(ji,jj,jk) = zw1(ji,jj,jk) + zw2(ji,jj,jk)
                     !
                     ! a_{23} is at vw-point, \partial_{z} u is at uw-point. Double 2-points interpolation in x and y now
                     !
                     zw3(ji,jj,jk) = (      dzuv(ji  ,jj+1,jk) * e1u(ji,jj+1) * e3uw_n(ji,jj+1,jk) +  &
                                   &        dzuv(ji  ,jj  ,jk) * e1u(ji,jj  ) * e3uw_n(ji,jj  ,jk) )  &
                                   & * ( var_ten(ji+1,jj  ,jk,ia23) + var_ten(ji,jj,jk,ia23) ) * 0.25_wp 
                  END DO
               END DO
            END DO
            !
            ! Surface boundary condition on the shear of u
            !
            zw3(:,:,1) = 0._wp
            !
            ! Finalisation of the Y-directed diffusive flux of U with vertical interpolation of zw3
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     zw1(ji,jj,jk) = zw1(ji,jj,jk) + ( zw3(ji,jj,jk) + zw3(ji,jj,jk+1) ) * 0.5_wp
                  END DO 
               END DO
            END DO
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hhdiff', zw0 , 'T', 1., &
                              &               zw1 , 'F', 1.  )
            !
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     !
                     ! Computing the divergence of the fluxes, centred in U-point 
                     ! 
                     zw2(ji,jj,jk) = ( zw0(ji+1,jj,jk) - zw0(ji,jj  ,jk) ) / ( e1e2u(ji,jj) * e3u_n(ji,jj,jk) )
                     zw3(ji,jj,jk) = ( zw1(ji  ,jj,jk) - zw1(ji,jj-1,jk) ) / ( e1e2u(ji,jj) * e3u_n(ji,jj,jk) )
                     !
                     ! Updating the component
                     !
                     pca(ji,jj,jk) = pca(ji,jj,jk)  + 0.5_wp * umask(ji,jj,jk) * ( zw2(ji,jj,jk) + zw3(ji,jj,jk) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 
            !
            ! Computation of the Horizontal Gradient of U
            !
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 2, jpi  
                     dxuv(ji,jj,jk) = ( pcn(ji,jj  ,jk) - pcn(ji-1,jj,jk) ) * fmask(ji,jj,jk) * r1_e1t(ji,jj) 
                     dyuv(ji,jj,jk) = ( pcn(ji,jj+1,jk) - pcn(ji  ,jj,jk) ) * tmask(ji,jj,jk) * r1_e2t(ji,jj)  
                  END DO
               END DO
            END DO
            !
            ! Computation of the vertical derivative of U
            !
            DO jk = 2, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi   
                     dzuv(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) - pcn(ji,jj,jk) ) * wvmask(ji,jj,jk) / e3vw_n(ji,jj,jk) 
                  END DO
               END DO
            END DO
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hhdiff', dxuv , 'F', 1.,     &
                              &               dyuv , 'T', 1.      )
            !
            ! Computation of the X-directed diffusive flux of V
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     !
                     ! a_{11} is at T-point, needs 4-points interpolation
                     ! 
                     zw0(ji,jj,jk) = ( var_ten(ji+1,jj+1,jk,ia11) + var_ten(ji  ,jj+1,jk,ia11) + &
                                   &   var_ten(ji  ,jj  ,jk,ia11) + var_ten(ji+1,jj  ,jk,ia11) ) & 
                                   & *    dxuv(ji,jj,jk) * e2f(ji,jj) * e3f_n(ji,jj,jk) * 0.25_wp
                     !
                     ! \partial_{x} u is at T-point, needs 4-points interpolation
                     !
                     zw1(ji,jj,jk) = (    dyuv(ji+1,jj+1,jk) * e2t(ji+1,jj+1) * e3t_n(ji+1,jj+1,jk) +  &
                                   &      dyuv(ji  ,jj+1,jk) * e2t(ji  ,jj+1) * e3t_n(ji  ,jj+1,jk) +  &
                                   &      dyuv(ji  ,jj  ,jk) * e2t(ji  ,jj  ) * e3t_n(ji  ,jj  ,jk) +  &
                                   &      dyuv(ji+1,jj  ,jk) * e2t(ji+1,jj  ) * e3t_n(ji+1,jj  ,jk) )  &
                                   & * var_ten(ji  ,jj  ,jk,ia12) * 0.25_wp
                     !
                     ! Cumulate on zw1
                     !
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + zw1(ji,jj,jk)
                     !
                     ! a_{23} is at vw-point, \partial_{z} u is at uw-point. Double 2-points interpolation in x and y now
                     !
                     zw2(ji,jj,jk) = (      dzuv(ji+1,jj  ,jk) * e1u(ji+1,jj) * e3uw_n(ji+1,jj,jk) +  &
                                   &        dzuv(ji  ,jj  ,jk) * e1u(ji  ,jj) * e3uw_n(ji  ,jj,jk) )  &
                                   & * ( var_ten(ji  ,jj+1,jk,ia23) + var_ten(ji,jj,jk,ia23) ) * 0.25_wp 
                  END DO
               END DO
            END DO
            !
            ! Surface boundary condition on the shear of U
            !
            zw2(:,:,1) = 0._wp
            !
            ! Finalisation of the X-directed diffusive flux of U with vertical interpolation of zw2
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + ( zw2(ji,jj,jk) + zw2(ji,jj,jk+1) ) * 0.5_wp
                  END DO 
               END DO
            END DO
            !
            ! Reinitialise support vectors
            !  
            zw1 = 0._wp
            zw2 = 0._wp
            zw3 = 0._wp
            !
            ! Computation of the Y-directed diffusive flux of V
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     !
                     ! a_{21}\partial_{y} v is at f-point, needs 4-point interpolation
                     !
                     zw1(ji,jj,jk) = ( dxuv(ji  ,jj  ,jk) * var_ten(ji  ,jj  ,jk,ia21) * e1f(ji  ,jj  ) * e3f_n(ji  ,jj  ,jk) +  &
                                   &   dxuv(ji-1,jj  ,jk) * var_ten(ji-1,jj  ,jk,ia21) * e1f(ji-1,jj  ) * e3f_n(ji-1,jj  ,jk) +  &
                                   &   dxuv(ji-1,jj-1,jk) * var_ten(ji-1,jj-1,jk,ia21) * e1f(ji-1,jj-1) * e3f_n(ji-1,jj-1,jk) +  &
                                   &   dxuv(ji  ,jj-1,jk) * var_ten(ji  ,jj-1,jk,ia21) * e1f(ji  ,jj-1) * e3f_n(ji  ,jj-1,jk) )  &
                                   & * 0.25_wp
                     !
                     ! a_{22}\partial_{x} v is at T-point
                     ! 
                     zw1(ji,jj,jk) = zw1(ji,jj,jk) + dxuv(ji,jj,jk) * var_ten(ji,jj,jk,ia22) * e1t(ji,jj) * e3t_n(ji,jj,jk)

                     !
                     ! Cumulate on zw0
                     !
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + zw1(ji,jj,jk)
                     !
                     ! a_{13}\partial_{z} u is at uw-point, needs 2-points interpolation in x and 2-points interpolation in z (later)
                     !
                     zw2(ji,jj,jk) = ( dzuv(ji  ,jj  ,jk) * var_ten(ji  ,jj  ,jk,ia23) * e1v(ji  ,jj  ) * e3f_n(ji  ,jj  ,jk) +  &
                                   &   dzuv(ji  ,jj-1,jk) * var_ten(ji  ,jj-1,jk,ia23) * e1v(ji  ,jj-1) * e3f_n(ji  ,jj-1,jk) )  &
                                   & * 0.5_wp 
                  END DO
               END DO
            END DO
            !
            ! Surface boundary condition on the shear of U
            !
            zw2(:,:,1) = 0._wp
            !
            ! Finalisation of the X-directed diffusive flux of U with vertical interpolation of zw2
            !
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi  
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + ( zw2(ji,jj,jk) + zw2(ji,jj,jk+1) ) * 0.5_wp
                  END DO 
               END DO
            END DO
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hhdiff', zw0 , 'F', 1., &
                              &               zw1 , 'T', 1.  )

            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     !
                     ! Computing the divergence of the fluxes, centred in V-point 
                     ! 
                     zw2(ji,jj,jk) = ( zw0(ji,jj  ,jk) - zw0(ji-1,jj,jk) ) / ( e1e2v(ji,jj) * e3v_n(ji,jj,jk) )
                     zw3(ji,jj,jk) = ( zw1(ji,jj+1,jk) - zw1(ji  ,jj,jk) ) / ( e1e2v(ji,jj) * e3v_n(ji,jj,jk) )
                     !
                     ! Updating the component
                     !
                     pca(ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  )  + 0.5_wp * vmask(ji,jj,jk) * ( zw2(ji,jj,jk) + zw3(ji,jj,jk) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_hhdiff: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0, zw1, zw2)

   END SUBROUTINE tlu_hhdiff
   ! [tlu_hhdiff]


END MODULE tluhdifDYN

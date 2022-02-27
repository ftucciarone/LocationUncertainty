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
MODULE tlusbc
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
   PUBLIC tlu_zadv_sbc  ! Called by step.f90
   PUBLIC tlu_zzsbc     ! CALLED BY STEP.F90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS


   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zadv_sbc ***
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
   !> @snippet this tlu_zadv_sbc 
   ! [tlu_zadv_sbc]
   SUBROUTINE tlu_zadv_sbc( kt, cdtype, kcmp, pcn, pca)
      INTEGER                         , INTENT(in   ) ::   kt                     ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kcmp                   ! Define velocity type for calculation
      CHARACTER(len=1)                , INTENT(in   ) ::   cdtype                 ! =U, V or W (component indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pcn, pca               ! modified advection term for this component
      !
      INTEGER                                         ::   ji, jj                 ! dummy loop indices
      REAL(wp), POINTER, DIMENSION(:,:,:)             ::   zwa0                   ! workspace
      REAL(wp)                                        ::   zwt, zwb               ! top and bottom components for finite difference
      REAL(wp)                                        ::   ztau                   ! wind stress on u and v points

      REAL(wp), PARAMETER                             ::   taurnm = 1             ! renormalization constant for wind stress
      !
      !!----------------------------------------------------------------------
      !
!      IF( ln_timing ) CALL timing_start('tlu_zadv_sbc')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zadv_sbc : ', cdtype, ' component of surface layer noise advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( zwa0(jpi,jpj,2))
      !

      SELECT CASE( kcmp )     !Select velocity type to interpolate sigma dbt in to corresponding mesh point for advection calculation

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity        
            zwa0 = 0._wp
            DO jj = 1, jpj
               DO ji = 1, jpim1
                  !  Interpolation of the wind stress on the advection grid (originally on T-grid)
                  ztau = 0.5_wp * ( utau(ji+1,jj  ) + utau(ji,jj  ) ) * umask(ji,jj,1) * taurnm
                  ! Interpolation from w -> uw and adding scale factor for the first two levels)
                  zwa0(ji,jj,1) = 0.5_wp * ( wnoi(ji+1,jj,1) + wnoi(ji  ,jj,1) ) * wumask(ji,jj,1) / e3uw_n(ji,jj,1) 
                  zwa0(ji,jj,2) = 0.5_wp * ( wnoi(ji+1,jj,2) + wnoi(ji  ,jj,2) ) * wumask(ji,jj,2) / e3uw_n(ji,jj,2) 
                  !
                  ! Vertical derivative to the first layer (Vertical average included)
                  !
                  !                   Top        (   lvl0 -  lvl1   )
                  zwt =  (     ztau     - pcn(ji,jj,1) ) * wumask(ji,jj,1) * zwa0(ji,jj,1)
                  !                  Below       (   lvl1  - lvl2   )
                  zwb =  ( pcn(ji,jj,1) - pcn(ji,jj,2) ) * wumask(ji,jj,2) * zwa0(ji,jj,2)
                  ! Average in z;   0.5 * ( Top + Below )
                  pca(ji,jj,1) = pca(ji,jj,1) - 0.5_wp * ( zwt + zwb ) * umask(ji,jj,1)
               END DO
            END DO
            
         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity    
            zwa0 = 0._wp
            DO jj = 1, jpjm1
               DO ji = 1, jpi
                  !  Interpolation of the wind stress on the advection grid (originally on T-grid)
                  ztau = 0.5_wp * (vtau(ji  ,jj+1) + vtau(ji  ,jj  ) ) * vmask(ji,jj,1) * taurnm
                  ! Interpolation from w -> vw and adding scale factor for the first two levels)
                  zwa0(ji,jj,1) = 0.5_wp * ( wnoi(ji,jj+1,1) + wnoi(ji,jj,1) ) * wvmask(ji,jj,1) / e3vw_n(ji,jj,1) 
                  zwa0(ji,jj,2) = 0.5_wp * ( wnoi(ji,jj+1,2) + wnoi(ji,jj,2) ) * wvmask(ji,jj,2) / e3vw_n(ji,jj,2) 
                  !
                  ! Vertical derivative to the first layer (Vertical average included)
                  !
                  !                   Top        (   lvl0 -  lvl1   )
                  zwt =  (     ztau     - pcn(ji,jj,1) ) * wvmask(ji,jj,1) * zwa0(ji,jj,1)
                  !                  Below       (   lvl1  - lvl2   )
                  zwb =  ( pcn(ji,jj,1) - pcn(ji,jj,2) ) * wvmask(ji,jj,2) * zwa0(ji,jj,2)
                  ! Average in z;   0.5 * ( Top + Below )
                  pca(ji,jj,1) = pca(ji,jj,1) - 0.5_wp * ( zwt + zwb ) * vmask(ji,jj,1)
               END DO
            END DO

         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','tlu_zadv_sbc: wrong value for kcmp'  )
      END SELECT


      DEALLOCATE(zwa0)
      !
      IF( ln_timing ) CALL timing_stop('tlu_zadv_sbc')   ! [NEMO] check
      !
   END SUBROUTINE tlu_zadv_sbc
   ! [tlu_zadv_sbc]




   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zzsbc ***
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
   ! @snippet this tlu_zzsbc
   ! [tlu_zzsbc]
   SUBROUTINE tlu_zzsbc( kt, cdtype, kcmp , pcb, pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcb                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! now field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp), DIMENSION(jpi,jpj,2)                    ::   zcnt_b                ! averaged before-component on t
      REAL(wp), DIMENSION(jpi,jpj,2)                    ::   zcnt_n                ! averaged now-component on t

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zzsbc')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zzsbc : ', cdtype, ' surface layer component of vertical diffusion (wind)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity

            ! Interpolation of the velocity component
            DO jk = 1, 2
               DO jj = 1, jpj
                  DO ji = 2, jpi   
                     zcnt_b(ji,jj,jk) = 0.5_wp * ( pcb(ji,jj,jk) + pcb(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                     zcnt_n(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO

            DO jj = 1, jpj
               DO ji = 2, jpi  
                     !
                     ! Surface boundary contition, Diffusive condition (au)
                     !
                     ! Surface condition;             Top      (   (k+1) -   k      )
                     zwt =  ( var_tau(ji,jj) * utau_b(ji,jj) - var_ten(ji,jj,1,ia33) * zcnt_b(ji,jj,1) ) &
                         &  * wmask(ji,jj,1) / e3w_n(ji,jj,1)
                     ! First layer;                  Below      (     k   - (k-1)    )
                     zwb =  ( var_ten(ji,jj,1,ia33) * zcnt_b(ji,jj,1) - var_ten(ji,jj,2,ia33) * zcnt_b(ji,jj,2) ) &
                         &  * wmask(ji,jj,2) / e3w_n(ji,jj,2)
                     ! Double derivative in z;     Top - Below  ( (k+1) - 2k + (k-1) )
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     zw0(ji,jj,1) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,1) / e3t_n(ji,jj,1)
                     !
                     ! Surface boundary contition, Drift condition (a)
                     !
                     !            Top      (   (k+1) -   k      )
                     zwt =  ( var_tau(ji,jj) - var_ten(ji,jj,1,ia33) ) &
                         &  * wmask(ji,jj,1) / e3w_n(ji,jj,1)
                     ! First layer;                  Below      (     k   - (k-1)    )
                     zwb =  ( var_ten(ji,jj,1,ia33) - var_ten(ji,jj,2,ia33) ) &
                         &  * wmask(ji,jj,2) / e3w_n(ji,jj,2)
                     ! Double derivative in z;     Top - Below  ( (k+1) - 2k + (k-1) )
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     zw0(ji,jj,1) = zw0(ji,jj,1) - 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,1) / e3t_n(ji,jj,1)
               END DO
            END DO

            DO jj = 1, jpj
               DO ji = 1, jpim1
                  ! Averaging and update of the velocity component
                  pca (ji  ,jj  ,1  ) = pca(ji  ,jj  ,1  ) + 0.5_wp * umask(ji,jj,1)   &
                                      & * ( zw0 (ji+1,jj  ,1  ) + zw0 (ji  ,jj  ,1  ) )
               END DO  
            END DO


         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 

            ! Interpolation of the velocity component
            DO jk = 1, 2
               DO jj = 2, jpj
                  DO ji = 1, jpi   
                     zcnt_b(ji,jj,jk) = 0.5_wp * ( pcb(ji,jj,jk) + pcb(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                     zcnt_n(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO

            DO jj = 1, jpj
               DO ji = 2, jpi   
                     !
                     ! Surface boundary contition, Diffusive condition (au)
                     !
                     ! Surface condition;             Top      (   (k+1) -   k      )
                     zwt =  ( var_tau(ji,jj) * vtau_b(ji,jj) - var_ten(ji,jj,1,ia33) * zcnt_b(ji,jj,1) ) &
                         &  * wmask(ji,jj,1) / e3w_n(ji,jj,1)
                     ! First layer;                  Below      (    k   - (k-1)    )
                     zwb =  ( var_ten(ji,jj,1,ia33) * zcnt_b(ji,jj,1) - var_ten(ji,jj,2,ia33) * zcnt_b(ji,jj,2) ) &
                         &  * wmask(ji,jj,2) / e3w_n(ji,jj,2)
                     ! Double derivative in z;     Top - Below  ( (k+1) - 2k + (k-1) )
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     zw0(ji,jj,1) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,1) / e3t_n(ji,jj,1)
                     !
                     ! Surface boundary contition, Drift condition (a)
                     !
                     ! Surface condition;             Top      (   (k+1) -   k      )
                     zwt =  ( var_tau(ji,jj) - var_ten(ji,jj,1,ia33) )  &
                         &  * wmask(ji,jj,1) / e3w_n(ji,jj,1)
                     ! First layer;                  Below      (     k   - (k-1)    )
                     zwb =  ( var_ten(ji,jj,1,ia33) - var_ten(ji,jj,2,ia33) ) &
                         &  * wmask(ji,jj,2) / e3w_n(ji,jj,2)
                     ! Double derivative in z;     Top - Below  ( (k+1) - 2k + (k-1) )
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     zw0(ji,jj,1) = zw0(ji,jj,1) - 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,1) / e3t_n(ji,jj,1)
               END DO
            END DO
           
            DO jj = 1, jpjm1
               DO ji = 1, jpi
                  ! Averaging and update of the velocity component
                  pca (ji  ,jj  ,1  ) = pca(ji  ,jj  ,1  ) + 0.5_wp * vmask(ji,jj,1)   &
                                      & * ( zw0 (ji  ,jj+1,1  ) + zw0 (ji  ,jj  ,1  ) )
               END DO  
            END DO


         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_zzsbc: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0)
      !
      IF( ln_timing ) CALL timing_stop('tlu_zzsbc')   ! [NEMO] check
      !
   END SUBROUTINE tlu_zzsbc
   ! [tlu_zzsbc]
 





END MODULE tlusbc


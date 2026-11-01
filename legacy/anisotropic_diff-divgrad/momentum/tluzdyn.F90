!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tluzdyn
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
!> @details         
!!
!> @par           Code specifics
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
!!
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tluzdyn
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
   PUBLIC tlu_zadv_noi  ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_zzdrift   ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_zzdiff    ! Called by tlu_bcdyn in tlu.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS


   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zadv_noi ***
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
   !> @snippet this tlu_zadv_noi 
   ! [tlu_zadv_noi]
   SUBROUTINE tlu_zadv_noi( kt, cdtype, kcmp, pcn, pca)
      INTEGER                         , INTENT(in   ) ::   kt                     ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kcmp                   ! Define velocity type for calculation
      CHARACTER(len=1)                , INTENT(in   ) ::   cdtype                 ! =U, V or W (component indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pcn, pca               ! modified advection term for this component
      !
      INTEGER                                         ::   ji, jj, jk             ! dummy loop indices
      REAL(wp)                                        ::   zwb, zwt               ! FD scalars (bottom, top)        
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zadv_noi')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zadv_noi :', cdtype, ' component of vertical-noise advection (no surface layer)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !

      SELECT CASE( kcmp )     

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity
            !
            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji= 1, jpim1
                     !                   Top        (   (k+1) -   k      )
                     zwt =  0.5_wp * ( wnoi(ji+1,jj,jk-1) + wnoi(ji,jj,jk-1) ) * wumask(ji,jj,jk-1) &
                         &         * (  pcn(ji  ,jj,jk-1) -  pcn(ji,jj,jk  ) ) / e3uw_n(ji,jj,jk-1)
                     !                  Below       (     k   - (k-1)    )
                     zwb =  0.5_wp * ( wnoi(ji+1,jj,jk  ) + wnoi(ji,jj,jk  ) ) * wumask(ji,jj,jk  ) &
                         &         * (  pcn(ji  ,jj,jk  ) -  pcn(ji,jj,jk+1) ) / e3uw_n(ji,jj,jk  )
                     ! Average in z;   0.5 * ( Top + Below )
                     pca(ji,jj,jk) = pca(ji,jj,jk) - 0.5_wp * ( zwt + zwb ) * umask(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity
            !
            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji= 1, jpi
                     !                   Top        (   (k+1) -   k      )
                     zwt =  0.5_wp * ( wnoi(ji,jj+1,jk-1) + wnoi(ji,jj,jk-1) ) * wvmask(ji,jj,jk-1) &
                         &         * (  pcn(ji,jj  ,jk-1) -  pcn(ji,jj,jk  ) ) / e3vw_n(ji,jj,jk-1) 
                     !                  Below       (     k   - (k-1)    )
                     zwb =  0.5_wp * ( wnoi(ji,jj+1,jk  ) + wnoi(ji,jj,jk  ) ) * wvmask(ji,jj,jk) &
                         &         * (  pcn(ji,jj  ,jk  ) -  pcn(ji,jj,jk+1) ) / e3vw_n(ji,jj,jk  ) 
                     ! Average in z;   0.5 * ( Top + Below )
                     pca(ji,jj,jk) = pca(ji,jj,jk) - 0.5_wp * ( zwt + zwb ) * vmask(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','tlu_zadv_noi: wrong value for kcmp'  )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_zadv_noi')   ! [NEMO] check
      !
   END SUBROUTINE tlu_zadv_noi
   ! [tlu_zadv_noi]




   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zzdiff ***
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
   ! @snippet this tlu_zzdiff
   ! [tlu_zzdiff]
   SUBROUTINE tlu_zzdiff( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zcnt                  ! averaged component on t grid
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zzdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zzdiff : ', cdtype, ' component of vertical diffusion (no surface layer)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Velocity averages
      zcnt = 0._wp     ! For diagonal components (at T-point)

      ! Support vectors
      zw0 = 0._wp


      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity

            ! Interpolation of the velocity component
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 2, jpi   
                     zcnt(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO

            DO jj = 1,jpj
               DO ji = 2,jpi 
                  !                                ! ================
                  DO jk = 2, jpkm1                 ! Horizontal slab
                     !                             ! ================
                     !
                     !                               Top        (   (k+1) -   k      )
                     zwt =   ( var_ten(ji,jj,jk-1,ia33) * zcnt(ji,jj,jk-1) -   & 
                         &     var_ten(ji,jj,jk  ,ia33) * zcnt(ji,jj,jk  ) ) * &
                         &     wmask(ji,jj,jk) / e3w_n(ji,jj,jk  )
                     !                              Below       (     k   - (k-1)    )
                     zwb =   ( var_ten(ji,jj,jk  ,ia33) * zcnt(ji,jj,jk  ) -   & 
                         &     var_ten(ji,jj,jk+1,ia33) * zcnt(ji,jj,jk+1) ) * &
                         &     wmask(ji,jj,jk+1) / e3w_n(ji,jj,jk+1)
                     ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
                     zw0(ji,jj,jk) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                     !
                     !                             ! ================
                  END DO                           !   End of slab
                  !                                ! ================
               END DO
            END DO

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji = 1, jpim1
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + umask(ji,jj,jk) * 0.5_wp  &
                                        & * ( zw0 (ji+1,jj  ,jk  ) + zw0 (ji  ,jj  ,jk  ) ) 
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 

            ! Interpolation of the velocity component
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 1, jpi   
                     zcnt(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji,jj-1,jk) ) * tmask(ji,jj,jk) 
                  END DO
               END DO
            END DO

            DO jj = 1,jpj
               DO ji = 2,jpi 
                  !                                ! ================
                  DO jk = 2, jpkm1                 ! Horizontal slab
                     !                             ! ================
                     !
                     !                               Top        (   (k+1) -   k      )
                     zwt =   ( var_ten(ji,jj,jk-1,ia33) * zcnt(ji,jj,jk-1) -   & 
                         &     var_ten(ji,jj,jk  ,ia33) * zcnt(ji,jj,jk  ) ) * &
                         &     wmask(ji,jj,jk) / e3w_n(ji,jj,jk  )
                     !                              Below       (     k   - (k-1)    )
                     zwb =   ( var_ten(ji,jj,jk  ,ia33) * zcnt(ji,jj,jk  ) -   & 
                         &     var_ten(ji,jj,jk+1,ia33) * zcnt(ji,jj,jk+1) ) * &
                         &     wmask(ji,jj,jk+1)  / e3w_n(ji,jj,jk+1)
                     ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
                     zw0(ji,jj,jk) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                     !
                     !                             ! ================
                  END DO                           !   End of slab
                  !                                ! ================
               END DO
            END DO

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + vmask(ji,jj,jk) * 0.5_wp *  &
                                        & ( zw0 (ji  ,jj+1,jk  ) + zw0 (ji  ,jj  ,jk  ) ) 
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            ! 

         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_zzdiff: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0)


   END SUBROUTINE tlu_zzdiff
   ! [tlu_zzdiff]



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zzdrift ***
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
   ! @snippet this tlu_zzdrift
   ! [tlu_zzdrift]
   SUBROUTINE tlu_zzdrift( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zzdrift')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zzdrift : ', cdtype, ' component of vertical drift (no surface layer)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Support vectors
      zw0 = 0._wp

      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity

            DO jj = 1,jpj
               DO ji = 1,jpi 
                  !                                ! ================
                  DO jk = 2, jpkm1                 ! Horizontal slab
                     !                             ! ================
                     !
                     !                               Top        (   (k+1) -   k      )
                     zwt =   ( var_ten(ji,jj,jk-1,ia33) - var_ten(ji,jj,jk  ,ia33) ) * &
                         &     wmask(ji,jj,jk)  / e3w_n(ji,jj,jk  )
                     !                              Below       (     k   - (k-1)    )
                     zwb =   ( var_ten(ji,jj,jk  ,ia33) - var_ten(ji,jj,jk+1,ia33) ) * &
                         &     wmask(ji,jj,jk+1)  / e3w_n(ji,jj,jk+1)
                     ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
                     zw0(ji,jj,jk) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                     !
                     !                             ! ================
                  END DO                           !   End of slab
                  !                                ! ================
               END DO
            END DO

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji = 1, jpim1
                     ! Averaging and update of the velocity component
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) - 0.5_wp * umask(ji,jj,jk) *      &
                                        &   pcn(ji  ,jj  ,jk  ) * ( zw0 (ji+1,jj  ,jk  )          &
                                        &                         + zw0 (ji  ,jj  ,jk  ) ) 
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 

            DO jj = 1,jpj
               DO ji = 1,jpi 
                  !                                ! ================
                  DO jk = 2, jpkm1                 ! Horizontal slab
                     !                             ! ================
                     !
                     !                               Top        (   (k+1) -   k      )
                     zwt =   ( var_ten(ji,jj,jk-1,ia33) - var_ten(ji,jj,jk  ,ia33) ) * &
                         &     wmask(ji,jj,jk) / e3w_n(ji,jj,jk  )
                     !                              Below       (     k   - (k-1)    )
                     zwb =   ( var_ten(ji,jj,jk  ,ia33) - var_ten(ji,jj,jk+1,ia33) ) * &
                         &     wmask(ji,jj,jk+1) / e3w_n(ji,jj,jk+1)
                     ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
                     zw0(ji,jj,jk) = + 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                     !
                     !                             ! ================
                  END DO                           !   End of slab
                  !                                ! ================
               END DO
            END DO

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     ! Averaging and update of the velocity component
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) - 0.5_wp * vmask(ji,jj,jk) *      &
                                        &   pcn(ji  ,jj  ,jk  ) * ( zw0 (ji  ,jj+1,jk  )          &
                                        &                         + zw0 (ji  ,jj  ,jk  ) ) 
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            ! 


         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_zzdrift: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0)

   END SUBROUTINE tlu_zzdrift
   ! [tlu_zzdrift]


END MODULE tluzdyn

!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlucmp
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
!> @brief         
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!
!!------------------------------------------------------------------------------
MODULE tlucmp
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
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_hhcmp      ! Called by tlu_bcdyn in tlu.F90
   PUBLIC tlu_hzcmp      ! Called by tlu_bcdyn in tlu.F90
   PUBLIC tlu_zzcmp      ! Called by tlu_bcdyn in tlu.F90
   PUBLIC tlu_cmpsbc     ! Called by tlu_bcdyn in tlu.F90
   ! [public_sub]  
   !
   !! * Substitutions [TODO] enable substitutions for NEMO
! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------

CONTAINS

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hhcmp ***
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
   ! @snippet this tlu_hhcmp
   ! [tlu_hhcmp]
   SUBROUTINE tlu_hhcmp( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp)                                          ::   rmu_3                 ! Molecular viscosity divided by 3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw1, zw2              ! workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwb, zwa              ! local scalars (FD before after)
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_hhcmp')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hhcmp : Horizontal compressibility components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))
      ALLOCATE(zw1(jpi,jpj,jpk), zw2(jpi,jpj,jpk))

      ! Molecular viscosity divided by 3
      rmu_3 = 1.1e-6 / 3

      ! Support vectors
      zw0 = 0._wp
      zw1 = 0._wp
      zw2 = 0._wp

      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 
            !
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               !
               ! First component: Double x derivative using a centered scheme inside comp. domain 
               ! No halo is computed here!!!
               DO jj = 1,jpj
                  DO ji = 2,jpim1 
                     ! Before (     i   - (i-1)   )
                     zwb =  ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11)   &
                          & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia11) ) &
                          & * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk)
                     ! After  (   (i+1) -   i     )
                     zwa =  ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11)   &
                          & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) ) &
                          & * r1_e1u(ji  ,jj) * umask(ji  ,jj,jk)
                     ! Double derivative in x
                     zw0(ji,jj,jk) = zw0(ji,jj,jk)  + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zwa  - zwb ) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               ! scheme check'd with quad3D(zcnt) before jk loop. 23/06/2021
               !
               !
               ! Second component: Double y derivative
               DO jj = 2,jpjm1
                  DO ji = 1,jpi
                     ! Before (     j   - (j-1)   )
                     zwb = ( e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11)   &
                         & - e1t(ji,jj-1) * e3t_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia11) ) &
                         & * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk)
                     ! After (   (j+1) -   j     )
                     zwa = ( e1t(ji,jj+1) * e3t_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia11)   &
                         & - e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) ) &
                         & * r1_e2v(ji,jj  ) * vmask(ji,jj  ,jk)
                     ! Double derivative in x
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zwa  - zwb ) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               ! scheme check'd with quad3D(zcnt) before jk loop. 23/06/2021
               !
               !
               ! Third component: mixed derivatives 
               DO jj = 2, jpj
                  DO ji = 1, jpi
                     !  Derivation in j (inner divergence) leads to u-points
                     zw1(ji,jj,jk) = (  e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia12)   &
                                   & -  e1f(ji,jj-1) * e3f_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia12) ) &
                                   & * umask(ji,jj,jk) * r1_e1u(ji,jj)
                  END DO  
               END DO
               DO jj = 1, jpj
                  DO ji = 2, jpi
                     !  Derivation in i (inner divergence) leads to v-points
                     zw2(ji,jj,jk) = (  e2f(ji  ,jj) * e3f_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia12)   &
                                   & -  e2f(ji-1,jj) * e3f_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia12) ) &
                                   & * vmask(ji,jj,jk) * r1_e2v(ji,jj)
                  END DO  
               END DO
               DO jj = 2, jpj
                  DO ji = 2, jpi
                     ! Derivation in i (outer divergence) leads to T-points
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zw1(ji  ,jj  ,jk  ) - zw1(ji-1,jj  ,jk  ) ) / e3t_n(ji,jj,jk)
                  END DO  
               END DO
               DO jj = 2, jpj
                  DO ji = 2, jpi
                     ! Derivation in j (outer divergence) leads to T-points
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zw2(ji  ,jj  ,jk  ) - zw2(ji  ,jj-1,jk  ) ) / e3t_n(ji,jj,jk)
                  END DO  
               END DO
               ! scheme check'd with biLin3D(zcnf) before jk loop. 23/06/2021

               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
     
            CALL lbc_lnk_multi( 'tlu_hhdrift', zw0 , 'T', 1.)

            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji = 1, jpim1
                     ! Gradient (in x) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + umask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji+1,jj  ,jk  ) - zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            ! 


         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 
            !
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               !
               ! First component: Double x derivative using a centered scheme inside comp. domain 
               ! No halo is computed here!!!
               DO jj = 1,jpj
                  DO ji = 2,jpim1 
                     ! Before (     i   - (i-1)   )
                     zwb =  ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11)   &
                          & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia11) ) &
                          & * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk)
                     ! After  (   (i+1) -   i     )
                     zwa =  ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11)   &
                          & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) ) &
                          & * r1_e1u(ji  ,jj) * umask(ji  ,jj,jk)
                     ! Double derivative in x
                     zw0(ji,jj,jk) = zw0(ji,jj,jk)  + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zwa  - zwb ) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               ! scheme check'd with quad3D(zcnt) before jk loop. 23/06/2021
               !
               !
               ! Second component: Double y derivative
               DO jj = 2,jpjm1
                  DO ji = 1,jpi
                     ! Before (     j   - (j-1)   )
                     zwb = ( e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11)   &
                         & - e1t(ji,jj-1) * e3t_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia11) ) &
                         & * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk)
                     ! After (   (j+1) -   j     )
                     zwa = ( e1t(ji,jj+1) * e3t_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia11)   &
                         & - e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) ) &
                         & * r1_e2v(ji,jj  ) * vmask(ji,jj  ,jk)
                     ! Double derivative in x
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zwa  - zwb ) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               ! scheme check'd with quad3D(zcnt) before jk loop. 23/06/2021
               !
               !
               ! Third component: mixed derivatives 
               DO jj = 2, jpj
                  DO ji = 1, jpi
                     !  Derivation in j (inner divergence) leads to u-points
                     zw1(ji,jj,jk) = (  e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia12)   &
                                   & -  e1f(ji,jj-1) * e3f_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia12) ) &
                                   & * umask(ji,jj,jk) * r1_e1u(ji,jj)
                  END DO  
               END DO
               DO jj = 1, jpj
                  DO ji = 2, jpi
                     !  Derivation in i (inner divergence) leads to v-points
                     zw2(ji,jj,jk) = (  e2f(ji  ,jj) * e3f_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia12)   &
                                   & -  e2f(ji-1,jj) * e3f_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia12) ) &
                                   & * vmask(ji,jj,jk) * r1_e2v(ji,jj)
                  END DO  
               END DO
               DO jj = 2, jpj
                  DO ji = 2, jpi
                     ! Derivation in i (outer divergence) leads to T-points
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zw1(ji  ,jj  ,jk  ) - zw1(ji-1,jj  ,jk  ) ) / e3t_n(ji,jj,jk)
                  END DO  
               END DO
               DO jj = 2, jpj
                  DO ji = 2, jpi
                     ! Derivation in j (outer divergence) leads to T-points
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) &
                                   & * ( zw2(ji  ,jj  ,jk  ) - zw2(ji  ,jj-1,jk  ) ) / e3t_n(ji,jj,jk)
                  END DO  
               END DO
               ! scheme check'd with biLin3D(zcnf) before jk loop. 23/06/2021

               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
     
            CALL lbc_lnk_multi( 'tlu_hhdrift', zw0 , 'T', 1.)

            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     ! Gradient (in y) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + vmask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji  ,jj+1,jk  ) - zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_hhcmp: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0, zw1, zw2)

   END SUBROUTINE tlu_hhcmp
   ! [tlu_hhcmp]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hzcmp ***
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
   ! @snippet this tlu_hzcmp
   ! [tlu_hzcmp]
   SUBROUTINE tlu_hzcmp( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp)                                          ::   rmu_3                 ! Molecular viscosity divided by 3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwb , zwt             ! local scalars (FD bottom top)
      REAL(wp)                                          ::   zwl , zwr             ! local scalars (FD left right)
      REAL(wp)                                          ::   zwbl, zwbr            ! local scalars (FD bottom-left bottom-right)
      REAL(wp)                                          ::   zwtl, zwtr            ! local scalars (FD top-left top-right)
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zacnk                 ! Average in k of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zacnij                ! Average in i or j of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zacnijk               ! Average in ijk of pcn
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_hzcmp')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hzcmp : Mixed Horizontal-Vertical compressibility components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Molecular viscosity divided by 3
      rmu_3 = 1.1e-6 / 3

      ! Velocity averages
      zacnk = 0._wp      ! average in k
      zacnij = 0._wp     ! average in i and j
      zacnijk = 0._wp    ! average in ijk simultaneously

      ! Support vectors
      zw0 = 0._wp


      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 


            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (i  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia13)
                     !              Below-right       (   (i  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia13)
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji  ,jj,jk)
                     !               Top-left         (   (i-1,k-1)   )
                     zwtl = var_ten(ji-1,jj  ,jk-1,ia13)
                     !              Below-left        (   (i-1,k  )   )
                     zwbl = var_ten(ji-1,jj  ,jk  ,ia13)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwl = e2u(ji-1,jj  ) * (zwtl - zwbl) * umask(ji-1,jj,jk)
                     ! Right-Left derivative (in x)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia23) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia23) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e1v(ji  ,jj  ) * (zwtr - zwbr) * vmask(ji,jj  ,jk)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = var_ten(ji  ,jj-1,jk-1,ia23)
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = var_ten(ji  ,jj-1,jk  ,ia23)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwl = e1v(ji  ,jj-1) * (zwtl - zwbl) * vmask(ji,jj-1,jk)
                     ! Right-Left derivative (in y)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia13)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk-1) * var_ten(ji-1,jj  ,jk-1,ia13) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia13) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk  ) * var_ten(ji-1,jj  ,jk  ,ia13)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Right-Left derivative (in y)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 1, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia23)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk-1) * var_ten(ji  ,jj-1,jk-1,ia23)
                     ! Top derivative in y;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia23) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk  ) * var_ten(ji  ,jj-1,jk  ,ia23)
                     ! Bottom derivative in y;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Derivative in z
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !
               !
               DO jj = 1, jpj
                  DO ji = 1, jpim1
                     ! Gradient (in x) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + umask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji+1,jj  ,jk  ) - zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================



         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia23)
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia23) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e1v(ji  ,jj  ) * (zwtr - zwbr) * vmask(ji,jj  ,jk)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = var_ten(ji  ,jj-1,jk-1,ia23) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = var_ten(ji  ,jj-1,jk  ,ia23)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwl = e1v(ji  ,jj-1) * (zwtl - zwbl) * vmask(ji,jj-1,jk)
                     ! Right-Left derivative (in y)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (i  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia13) 
                     !              Below-right       (   (i  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia13)  
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji  ,jj,jk)
                     !               Top-left         (   (i-1,k-1)   )
                     zwtl = var_ten(ji-1,jj  ,jk-1,ia13) 
                     !              Below-left        (   (i-1,k  )   )
                     zwbl = var_ten(ji-1,jj  ,jk  ,ia13)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwl = e2u(ji-1,jj  ) * (zwtl - zwbl) * umask(ji-1,jj,jk)
                     ! Right-Left derivative (in x)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia23)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk-1) * var_ten(ji  ,jj-1,jk-1,ia23) 
                     ! Top derivative in y;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia23) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk  ) * var_ten(ji  ,jj-1,jk  ,ia23)
                     ! Bottom derivative in y;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Derivative in z
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia13)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk-1) * var_ten(ji-1,jj  ,jk-1,ia13) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia13) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk  ) * var_ten(ji-1,jj  ,jk  ,ia13)
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Right-Left derivative (in y)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !
               !
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     ! Gradient (in y) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + vmask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji  ,jj+1,jk  ) - zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         CASE DEFAULT          
            CALL ctl_stop('STOP','tlu_hzcmp: wrong value for kcmp'  )
      END SELECT
      !

   END SUBROUTINE tlu_hzcmp
   ! [tlu_hzcmp]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_zzcmp ***
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
   ! @snippet this tlu_zzcmp
   ! [tlu_zzcmp]
   SUBROUTINE tlu_zzcmp( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp)                                          ::   rmu_3                 ! Molecular viscosity divided by 3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zzcmp')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_zzcmp : ', cdtype, ' component of vertical compressibility (no surface layer)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Molecular viscosity divided by 3
      rmu_3 = 1.1e-6 / 3

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
                     ! Gradient (in x) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + umask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji+1,jj  ,jk  ) - zw0 (ji  ,jj  ,jk  ) ) 
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
                     ! Gradient (in y) and update of the velocity component
                     ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + vmask(ji,jj,jk) * rmu_3        &
                                        & * ( zw0 (ji  ,jj+1,jk  ) - zw0 (ji  ,jj  ,jk  ) ) 
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            ! 


         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_zzcmp: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0)

   END SUBROUTINE tlu_zzcmp
   ! [tlu_zzcmp]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_cmpsbc ***
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
   ! @snippet this tlu_cmpsbc
   ! [tlu_cmpsbc]
   SUBROUTINE tlu_cmpsbc( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
      REAL(wp)                                          ::   rmu_3                 ! Molecular viscosity divided by 3
      REAL(wp), DIMENSION(jpi,jpj,2)                    ::   zcnt                  ! averaged component on t
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   zw0                   ! accumulation array
      REAL(wp)                                          ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_cmpsbc')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_cmpsbc : ', cdtype, ' surface layer component of compressibility (wind)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Molecular viscosity divided by 3
      rmu_3 = 1.1e-6 / 3

      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity

            DO jj = 1, jpj
               DO ji = 2, jpi  
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
                     zw0(ji,jj,1) = zw0(ji,jj,1) - 0.5_wp * ( zwt  - zwb ) * tmask(ji,jj,jk) / e3t_n(ji,jj,1)
               END DO
            END DO

            DO jj = 1, jpj
               DO ji = 1, jpim1
                  ! Gradient (in x) and update of the velocity component
                  ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                  pca (ji  ,jj  ,1  ) = pca(ji  ,jj  ,1  ) + umask(ji,jj,1) * rmu_3        &
                                    & * ( zw0 (ji+1,jj  ,1  ) - zw0 (ji  ,jj  ,1  ) )
               END DO  
            END DO



         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 

            DO jj = 1, jpj
               DO ji = 2, jpi   
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
                  ! Gradient (in y) and update of the velocity component
                  ! This is a compressibility component added to the RIGHT-HAND SIDE, so positive
                  pca (ji  ,jj  ,1  ) = pca(ji  ,jj  ,1  ) + vmask(ji,jj,1) * rmu_3        &
                                    & * ( zw0 (ji  ,jj+1,1  ) - zw0 (ji  ,jj  ,1  ) )
               END DO  
            END DO


         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_cmpsbc: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0)
      !
      IF( ln_timing ) CALL timing_stop('tlu_cmpsbc')   ! [NEMO] check
      !
   END SUBROUTINE tlu_cmpsbc
   ! [tlu_cmpsbc]
 





END MODULE tlucmp




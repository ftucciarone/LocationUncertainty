!! !!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tluhzdyn
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
!!          
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
MODULE tluhzdyn
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
!   USE tlusvd           ! Initialization of the noise field
!   USE dynvor           ! Computation of vorticity (see warnings)
   ! [mod_dep]
   !
   IMPLICIT NONE   
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_hzdiff    ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_hzdrift   ! Called by tlu_bcdyn in tlu.F90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hzdiff ***
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
   ! @snippet this tlu_hzdiff
   ! [tlu_hzdiff]
   SUBROUTINE tlu_hzdiff( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
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
      IF( ln_timing ) CALL timing_start('tlu_hzdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hzdiff : Mixed Horizontal-Vertical diffusion components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      ! Velocity averages
      zacnk = 0._wp      ! average in k
      zacnij = 0._wp     ! average in i and j
      zacnijk = 0._wp    ! average in ijk simultaneously

      ! Support vectors
      zw0 = 0._wp


      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 


            ! First of all I need to take the velocity component to the other grid
            ! performing thus u -> v 
            !
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 2, jpi
                     zwl = ( pcn(ji-1,jj+1,jk  ) + pcn(ji-1,jj  ,jk  ) ) * fmask(ji-1,jj,jk)
                     zwr = ( pcn(ji  ,jj+1,jk  ) + pcn(ji  ,jj  ,jk  ) ) * fmask(ji  ,jj,jk)
                     zacnij(ji,jj,jk) = 0.25_wp * ( zwl + zwr ) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO

            DO jk = 2, jpk
               DO jj = 1, jpjm1
                  DO ji = 2, jpi
                     zacnk(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji,jj,jk-1) ) * wumask(ji,jj,jk)
                     zacnijk(ji,jj,jk) = 0.5_wp * ( zacnij(ji,jj,jk) + zacnij(ji,jj,jk-1) ) *wvmask(ji,jj,jk)
                  END DO
               END DO
            END DO

            ! Both averages checkd
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hzdiff', zacnk , 'U', -1., zacnijk , 'V', -1. )


            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (i  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia13) * zacnk(ji  ,jj  ,jk-1)
                     !              Below-right       (   (i  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia13) * zacnk(ji  ,jj  ,jk  ) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji  ,jj,jk)
                     !               Top-left         (   (i-1,k-1)   )
                     zwtl = var_ten(ji-1,jj  ,jk-1,ia13) * zacnk(ji-1,jj  ,jk-1) 
                     !              Below-left        (   (i-1,k  )   )
                     zwbl = var_ten(ji-1,jj  ,jk  ,ia13) * zacnk(ji-1,jj  ,jk  )
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
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia23) * zacnijk(ji  ,jj  ,jk-1)
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia23) * zacnijk(ji  ,jj  ,jk  ) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e1v(ji  ,jj  ) * (zwtr - zwbr) * vmask(ji,jj  ,jk)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = var_ten(ji  ,jj-1,jk-1,ia23) * zacnijk(ji  ,jj-1,jk-1) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = var_ten(ji  ,jj-1,jk  ,ia23) * zacnijk(ji  ,jj-1,jk  )
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
                     zwtr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia13) * zacnk(ji  ,jj  ,jk-1)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk-1) * var_ten(ji-1,jj  ,jk-1,ia13) * zacnk(ji-1,jj  ,jk-1) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia13) * zacnk(ji  ,jj  ,jk  ) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk  ) * var_ten(ji-1,jj  ,jk  ,ia13) * zacnk(ji-1,jj  ,jk  )
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
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia23)*zacnijk(ji  ,jj  ,jk-1)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk-1) * var_ten(ji  ,jj-1,jk-1,ia23)*zacnijk(ji  ,jj-1,jk-1) 
                     ! Top derivative in y;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia23)*zacnijk(ji  ,jj  ,jk  ) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk  ) * var_ten(ji  ,jj-1,jk  ,ia23)*zacnijk(ji  ,jj-1,jk  )
                     ! Bottom derivative in y;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Derivative in z
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !
               !
               DO jj = 2, jpj
                  DO ji = 2, jpim1
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + 0.5_wp * umask(ji,jj,jk)        &
                                        & * ( zw0 (ji+1,jj  ,jk  ) + zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================



         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 


            ! First of all I need to take the velocity component to the other grid
            ! performing thus v -> u 
            !
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 2, jpi
                     zwl = pcn(ji  ,jj  ,jk  ) + pcn(ji-1,jj  ,jk  ) * fmask(ji,jj  ,jk)
                     zwr = pcn(ji  ,jj+1,jk  ) + pcn(ji-1,jj+1,jk  ) * fmask(ji,jj+1,jk)
                     zacnij(ji,jj,jk) = 0.25_wp * ( zwl + zwr ) * umask(ji,jj,jk)
                  END DO
               END DO
            END DO

            DO jk = 2, jpk
               DO jj = 1, jpjm1
                  DO ji = 2, jpi
                     zacnk(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji,jj,jk-1) ) * wvmask(ji,jj,jk)
                     zacnijk(ji,jj,jk) = 0.5_wp * ( zacnij(ji,jj,jk) + zacnij(ji,jj,jk-1) ) * wumask(ji,jj,jk)
                  END DO
               END DO
            END DO

            ! Both averages checkd
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hzdiff', zacnk , 'V', -1., zacnijk , 'U', -1. )


            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 2, jpj
                  DO ji= 2, jpi
                     !               Top-right        (   (j  ,k-1)   )
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia23) * zacnk(ji  ,jj  ,jk-1)
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia23) * zacnk(ji  ,jj  ,jk  ) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e1v(ji  ,jj  ) * (zwtr - zwbr) * vmask(ji,jj  ,jk)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = var_ten(ji  ,jj-1,jk-1,ia23) * zacnk(ji  ,jj-1,jk-1) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = var_ten(ji  ,jj-1,jk  ,ia23) * zacnk(ji  ,jj-1,jk  )
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
                     zwtr = var_ten(ji  ,jj  ,jk-1,ia13) * zacnijk(ji  ,jj  ,jk-1)
                     !              Below-right       (   (i  ,k  )   )
                     zwbr = var_ten(ji  ,jj  ,jk  ,ia13) * zacnijk(ji  ,jj  ,jk  ) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji,jj,jk)
                     !               Top-left         (   (i-1,k-1)   )
                     zwtl = var_ten(ji-1,jj  ,jk-1,ia13) * zacnijk(ji-1,jj  ,jk-1) 
                     !              Below-left        (   (i-1,k  )   )
                     zwbl = var_ten(ji-1,jj  ,jk  ,ia13) * zacnijk(ji-1,jj  ,jk  )
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
                     zwtr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia23)*zacnk(ji  ,jj  ,jk-1)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk-1) * var_ten(ji  ,jj-1,jk-1,ia23)*zacnk(ji  ,jj-1,jk-1) 
                     ! Top derivative in y;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia23)*zacnk(ji  ,jj  ,jk  ) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk  ) * var_ten(ji  ,jj-1,jk  ,ia23)*zacnk(ji  ,jj-1,jk  )
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
                     zwtr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia13)*zacnijk(ji  ,jj  ,jk-1)
                     !               Top-left         (   (j-1,k-1)   )
                     zwtl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk-1) * var_ten(ji-1,jj  ,jk-1,ia13)*zacnijk(ji-1,jj  ,jk-1) 
                     ! Right derivative in z;      ( Top-Right - Below-Right )
                     zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                     !              Below-right       (   (j  ,k  )   )
                     zwbr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia13)*zacnijk(ji  ,jj  ,jk  ) 
                     !              Below-left        (   (j-1,k  )   )
                     zwbl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk  ) * var_ten(ji-1,jj  ,jk  ,ia13)*zacnijk(ji-1,jj  ,jk  )
                     ! Left derivative in z;       ( Top-Left - Below-Left )
                     zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                     ! Right-Left derivative (in y)
                     zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
                  END DO
               END DO
               !
               !
               DO jj = 2, jpjm1
                  DO ji = 2, jpi
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  )  + 0.5_wp * vmask(ji,jj,jk)   &
                                        & * ( zw0 (ji  ,jj+1,jk  ) + zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         CASE DEFAULT          
            CALL ctl_stop('STOP','tlu_hzdiff: wrong value for kcmp'  )
      END SELECT
      !

   END SUBROUTINE tlu_hzdiff
   ! [tlu_hzdiff]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hzdrift ***
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
   ! @snippet this tlu_hzdrift
   ! [tlu_hzdrift]
   SUBROUTINE tlu_hzdrift( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
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
      IF( ln_timing ) CALL timing_start('tlu_hzdrift')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hzdrift : Mixed Horizontal-Vertical drift components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

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
               !
               !
               DO jj = 2, jpj
                  DO ji = 2, jpim1
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
                     zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji,jj,jk)
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
               DO jj = 2, jpjm1
                  DO ji = 2, jpi
                     ! Averaging and update of the velocity component
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     pca (ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  )  - 0.5_wp * vmask(ji,jj,jk) *     &
                                        &   pcn(ji  ,jj  ,jk  ) * ( zw0 (ji  ,jj+1,jk  )          &
                                        &                         + zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         CASE DEFAULT          
            CALL ctl_stop('STOP','tlu_hzdrift: wrong value for kcmp'  )
      END SELECT
      !

   END SUBROUTINE tlu_hzdrift
   ! [tlu_hzdrift]







END MODULE tluhzdyn


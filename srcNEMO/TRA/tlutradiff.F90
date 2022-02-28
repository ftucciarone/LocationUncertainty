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
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
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
   SUBROUTINE tlu_trahhdiff( kt, ptn, pta)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   zcnf                  ! averaged component on f grid
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   zw1, zw2              ! workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   zw0                   ! accumulation array
      REAL(wp)                                             ::   zwb, zwa              ! local scalars (FD before after)
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
      IF( ln_timing ) CALL timing_start('tlu_hhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trahhdiff : Horizontal diffusion components on velocity type '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))
      ALLOCATE(zw1(jpi,jpj,jpk), zw2(jpi,jpj,jpk))

      !

      DO jn = 1, jpts            !==  loop over the tracers  ==!

         ! Tracer averages
         zcnf = 0._wp     ! For extra-diagonal components (at f-point)

         ! Support vectors
         zw0 = 0._wp
         zw1 = 0._wp
         zw2 = 0._wp
         !
         ! Interpolation of the Tracer component
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1 
                  zcnf(ji,jj,jk) = 0.5_wp * ( ptn(ji  ,jj  ,jk  ,jn) + ptn(ji+1,jj  ,jk  ,jn) + &
                                   &          ptn(ji  ,jj+1,jk  ,jn) + ptn(ji+1,jj+1,jk  ,jn) ) * fmask(ji,jj,jk)  
               END DO
            END DO
         END DO
         ! Both averages checkd
         !
         ! Lateral boundary condition transfer across nodes
         CALL lbc_lnk_multi( 'tlu_trahhdiff', zcnf, 'F', 1. )
         !                                ! ================
         DO jk = 1, jpkm1                 ! Horizontal slab
            !                             ! ================
            !
            ! First component: Double x derivative using a centered scheme inside comp. domain 
            ! No halo is computed here!!!
            DO jj = 1,jpj
               DO ji = 2,jpim1 
                  ! Before (     i   - (i-1)   )
                  zwb =  ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * ptn(ji  ,jj,jk,jn)   &
                       & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia11) * ptn(ji-1,jj,jk,jn) ) &
                       & * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk)
                  ! After  (   (i+1) -   i     )
                  zwa =  ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11) * ptn(ji+1,jj,jk,jn)   &
                       & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * ptn(ji  ,jj,jk,jn) ) &
                       & * r1_e1u(ji  ,jj) * umask(ji  ,jj,jk)
                  ! Double derivative in x
                  zw0(ji,jj,jk) = zw0(ji,jj,jk)  + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) * (SC**2) &
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
                  zwb = ( e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * ptn(ji,jj  ,jk,jn)   &
                      & - e1t(ji,jj-1) * e3t_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia11) * ptn(ji,jj-1,jk,jn) ) &
                      & * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk)
                  ! After (   (j+1) -   j     )
                  zwa = ( e1t(ji,jj+1) * e3t_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia11) * ptn(ji,jj+1,jk,jn)   &
                      & - e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * ptn(ji,jj  ,jk,jn) ) &
                      & * r1_e2v(ji,jj  ) * vmask(ji,jj  ,jk)
                  ! Double derivative in x
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) * (SC**2)  &
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
                  zw1(ji,jj,jk) = (  e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia12) * zcnf(ji,jj  ,jk)   &
                                & -  e1f(ji,jj-1) * e3f_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia12) * zcnf(ji,jj-1,jk) ) &
                                & * umask(ji,jj,jk) * r1_e1u(ji,jj)
               END DO  
            END DO
            DO jj = 1, jpj
               DO ji = 2, jpi
                  !  Derivation in i (inner divergence) leads to v-points
                  zw2(ji,jj,jk) = (  e2f(ji  ,jj) * e3f_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia12) * zcnf(ji  ,jj,jk)   &
                                & -  e2f(ji-1,jj) * e3f_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia12) * zcnf(ji-1,jj,jk) ) &
                                & * vmask(ji,jj,jk) * r1_e2v(ji,jj)
               END DO  
            END DO
            DO jj = 2, jpj
               DO ji = 2, jpi
                  !  Derivation in i (outer divergence) leads to T-points
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) * (SC**2)  &
                                & * ( zw1(ji  ,jj  ,jk  ) - zw1(ji-1,jj  ,jk  ) ) / e3t_n(ji,jj,jk)
               END DO  
            END DO
            DO jj = 2, jpj
               DO ji = 2, jpi
                  !  Derivation in j (outer divergence) leads to T-points
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * tmask(ji,jj,jk) * (SC**2)  &
                                & * ( zw2(ji  ,jj  ,jk  ) - zw2(ji  ,jj-1,jk  ) ) / e3t_n(ji,jj,jk)
               END DO  
            END DO
            ! scheme check'd with biLin3D(zcnf) before jk loop. 23/06/2021

            !                             ! ================
         END DO                           !   End of slab
         !                                ! ================

         CALL lbc_lnk_multi( 'tlu_trahhdiff', zw0 , 'T', 1.)

         !                                ! ================
         DO jk = 1, jpkm1                 ! Horizontal slab
            !                             ! ================
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Averaging and update of the velocity component
                  ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zw0 (ji,jj,jk)
               END DO
            END DO
            !                             ! ================
         END DO                           !   End of slab
         !                                ! ================
      END DO
      !
      DEALLOCATE(zw0, zw1, zw2)

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
   SUBROUTINE tlu_trahzdiff( kt, ptn, pta)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptn                   ! before field
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


      ALLOCATE(zw0(jpi,jpj,jpk))


      DO jn = 1, jpts            !==  loop over the tracers  ==!

         ! Velocity averages
         zacni = 0._wp     ! average in i
         zacnj = 0._wp     ! average in j
         zacnik = 0._wp    ! average in ik simultaneously
         zacnjk = 0._wp    ! average in jk simultaneously

         ! Support vectors
         zw0 = 0._wp

         ! First of all I need to take the tracer to uw and vw grids:
         ! performing T -> v 
         ! performing T -> u
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zacni(ji,jj,jk) = 0.5_wp * ( ptn(ji  ,jj  ,jk  ,jn) + ptn(ji+1,jj  ,jk  ,jn) ) * umask(ji,jj,jk)
                  zacnj(ji,jj,jk) = 0.5_wp * ( ptn(ji  ,jj  ,jk  ,jn) + ptn(ji  ,jj+1,jk  ,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO

         DO jk = 2, jpk
            DO jj = 1, jpjm1
               DO ji = 2, jpim1
                  zacnik(ji,jj,jk) = 0.5_wp * ( zacni(ji,jj,jk) + zacni(ji,jj,jk-1) ) * wumask(ji,jj,jk)
                  zacnjk(ji,jj,jk) = 0.5_wp * ( zacnj(ji,jj,jk) + zacnj(ji,jj,jk-1) ) * wvmask(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Both averages checkd
         !
         ! Lateral boundary condition transfer across nodes
         !
         CALL lbc_lnk_multi( 'tlu_hzdiff', zacnik, 'U', -1., zacnjk, 'V', -1. )


         !                                ! ================
         DO jk = 2, jpkm1                 ! Horizontal slab
            !                             ! ================
            DO jj = 2, jpj
               DO ji= 2, jpi
                  !               Top-right        (   (i  ,k-1)   )
                  zwtr = var_ten(ji  ,jj  ,jk-1,ia13) * zacnik(ji  ,jj  ,jk-1)
                  !              Below-right       (   (i  ,k  )   )
                  zwbr = var_ten(ji  ,jj  ,jk  ,ia13) * zacnik(ji  ,jj  ,jk  ) 
                  ! Right derivative in z;      ( Top-Right - Below-Right )
                  zwr = e2u(ji  ,jj  ) * (zwtr - zwbr) * umask(ji  ,jj,jk)
                  !               Top-left         (   (i-1,k-1)   )
                  zwtl = var_ten(ji-1,jj  ,jk-1,ia13) * zacnik(ji-1,jj  ,jk-1) 
                  !              Below-left        (   (i-1,k  )   )
                  zwbl = var_ten(ji-1,jj  ,jk  ,ia13) * zacnik(ji-1,jj  ,jk  )
                  ! Left derivative in z;       ( Top-Left - Below-Left )
                  zwl = e2u(ji-1,jj  ) * (zwtl - zwbl) * umask(ji-1,jj,jk)
                  ! Right-Left derivative (in x)
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) * (SC**2)  / e3t_n(ji,jj,jk)
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
                  zwtr = var_ten(ji  ,jj  ,jk-1,ia23) * zacnjk(ji  ,jj  ,jk-1)
                  !              Below-right       (   (j  ,k  )   )
                  zwbr = var_ten(ji  ,jj  ,jk  ,ia23) * zacnjk(ji  ,jj  ,jk  ) 
                  ! Right derivative in z;      ( Top-Right - Below-Right )
                  zwr = e1v(ji  ,jj  ) * (zwtr - zwbr) * vmask(ji,jj  ,jk)
                  !               Top-left         (   (j-1,k-1)   )
                  zwtl = var_ten(ji  ,jj-1,jk-1,ia23) * zacnjk(ji  ,jj-1,jk-1) 
                  !              Below-left        (   (j-1,k  )   )
                  zwbl = var_ten(ji  ,jj-1,jk  ,ia23) * zacnjk(ji  ,jj-1,jk  )
                  ! Left derivative in z;       ( Top-Left - Below-Left )
                  zwl = e1v(ji  ,jj-1) * (zwtl - zwbl) * vmask(ji,jj-1,jk)
                  ! Right-Left derivative (in y)
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * r1_e1e2t(ji,jj) * &
                                   & ( zwr - zwl ) * tmask(ji,jj,jk) * (SC**2)  / e3t_n(ji,jj,jk)
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
                  zwtr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia13) * zacnik(ji  ,jj  ,jk-1)
                  !               Top-left         (   (j-1,k-1)   )
                  zwtl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk-1) * var_ten(ji-1,jj  ,jk-1,ia13) * zacnik(ji-1,jj  ,jk-1) 
                  ! Right derivative in z;      ( Top-Right - Below-Right )
                  zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                  !              Below-right       (   (j  ,k  )   )
                  zwbr = e2u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia13) * zacnik(ji  ,jj  ,jk  ) 
                  !              Below-left        (   (j-1,k  )   )
                  zwbl = e2u(ji-1,jj  ) * e3uw_n(ji-1,jj  ,jk  ) * var_ten(ji-1,jj  ,jk  ,ia13) * zacnik(ji-1,jj  ,jk  )
                  ! Left derivative in z;       ( Top-Left - Below-Left )
                  zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                  ! Right-Left derivative (in y)
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * (SC**2) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
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
                  zwtr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk-1) * var_ten(ji  ,jj  ,jk-1,ia23)*zacnjk(ji  ,jj  ,jk-1)
                  !               Top-left         (   (j-1,k-1)   )
                  zwtl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk-1) * var_ten(ji  ,jj-1,jk-1,ia23)*zacnjk(ji  ,jj-1,jk-1) 
                  ! Top derivative in y;      ( Top-Right - Below-Right )
                  zwt = r1_e1e2t(ji  ,jj  ) * ( zwtr - zwtl ) * wmask(ji,jj,jk-1) / e3w_n(ji  ,jj  ,jk-1) 
                  !              Below-right       (   (j  ,k  )   )
                  zwbr = e1u(ji  ,jj  ) * e3uw_n(ji  ,jj  ,jk  ) * var_ten(ji  ,jj  ,jk  ,ia23)*zacnjk(ji  ,jj  ,jk  ) 
                  !              Below-left        (   (j-1,k  )   )
                  zwbl = e1u(ji  ,jj-1) * e3uw_n(ji  ,jj-1,jk  ) * var_ten(ji  ,jj-1,jk  ,ia23)*zacnjk(ji  ,jj-1,jk  )
                  ! Bottom derivative in y;       ( Top-Left - Below-Left )
                  zwb = r1_e1e2t(ji  ,jj  ) * ( zwbr - zwbl ) * wmask(ji,jj,jk  ) / e3w_n(ji  ,jj  ,jk  )
                  ! Derivative in z
                  zw0(ji,jj,jk) = zw0(ji,jj,jk) + 0.5_wp * (zwt - zwb) * (SC**2) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk) 
               END DO
            END DO
            !
            !
            DO jj = 2, jpj
               DO ji = 2, jpim1
                  ! Averaging and update of the velocity component
                  ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                  pta (ji  ,jj  ,jk  ,jn) = pta(ji  ,jj  ,jk  ,jn) + zw0 (ji+1,jj  ,jk  )
               END DO  
            END DO
            !                             ! ================
         END DO                           !   End of slab
         !                                ! ================
      END DO
      !
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
   SUBROUTINE tlu_trazzdiff( kt, ptn, pta)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   zw0                   ! accumulation array
      REAL(wp)                                             ::   zwt, zwb              ! local scalars (FD top bottom)
      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zzdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trazzdiff : '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))

      DO jn = 1, jpts            !==  loop over the tracers  ==!
         ! Support vectors
         zw0 = 0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               !
               ! ==========================
               ! Surface boundary condition 
               ! ==========================
               !
               !                               Top        (   (k+1) -   k      )
               zwt =     var_ten(ji,jj,1  ,ia33) * ptn(ji,jj,1  ,jn) * (-1._wp) * &
                   &     wmask(ji,jj,1) / e3w_n(ji,jj,1  )
               !                              Below       (     k   - (k-1)    )
               zwb =   ( var_ten(ji,jj,1,ia33) * ptn(ji,jj,1,jn) -   & 
                   &     var_ten(ji,jj,2,ia33) * ptn(ji,jj,2,jn) ) * &
                   &     wmask(ji,jj,2) / e3w_n(ji,jj,2)
               ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
               zw0(ji,jj,1) = + 0.5_wp * ( zwt  - zwb ) * (SC**2) * tmask(ji,jj,1) / e3t_n(ji,jj,1)
               !
               !                          ! ================
               DO jk = 2, jpkm1           ! Horizontal slab
                  !                       ! ================
                  !
                  !                               Top        (   (k+1) -   k      )
                  zwt =   ( var_ten(ji,jj,jk-1,ia33) * ptn(ji,jj,jk-1,jn) -   & 
                      &     var_ten(ji,jj,jk  ,ia33) * ptn(ji,jj,jk  ,jn) ) * &
                      &     wmask(ji,jj,jk) / e3w_n(ji,jj,jk  )
                  !                              Below       (     k   - (k-1)    )
                  zwb =   ( var_ten(ji,jj,jk  ,ia33) * ptn(ji,jj,jk  ,jn) -   & 
                      &     var_ten(ji,jj,jk+1,ia33) * ptn(ji,jj,jk+1,jn) ) * &
                      &     wmask(ji,jj,jk+1) / e3w_n(ji,jj,jk+1)
                  ! Double derivative in z;    Top - Below   ( (k+1) - 2k + (k-1) )
                  zw0(ji,jj,jk) = + 0.5_wp * ( zwt  - zwb ) * (SC**2) * tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
                  !
                  !                       ! ================
               END DO                     !   End of slab
               !                          ! ================
            END DO
         END DO
         !                                ! ================
         DO jk = 1, jpkm1                 ! Horizontal slab   replace 1 with 2 in case of no sbc
            !                             ! ================
            DO jj = 1, jpj
               DO ji = 1, jpim1
                  ! Averaging and update of the velocity component
                  ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                  pta (ji  ,jj  ,jk  ,jn) = pta(ji  ,jj  ,jk  ,jn) + zw0 (ji  ,jj  ,jk)
               END DO  
            END DO
            !                             ! ================
         END DO                           !   End of slab
         !                                ! ================
      END DO

      DEALLOCATE(zw0)


   END SUBROUTINE tlu_trazzdiff
   ! [tlu_trazzdiff]

END MODULE tlutradiff

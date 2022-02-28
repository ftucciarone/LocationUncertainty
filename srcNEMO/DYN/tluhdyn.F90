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
MODULE tluhdyn
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
   PUBLIC tlu_hadv_noi     ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_hhdrift      ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_hhdiff       ! Called by tlu_bcdyn in tlu.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS

   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hadv_noi ***
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
   !> @details       
   !!
   !> @param[in]     kt: ocean time-step integer indexi
   !> @param[in]     kcmp: Define velocity type for calculation
   !> @param[in]     cdtype: =U, V or W (component indicator)
   !> @param[inout]  pcn, pca:  modified advection term for this component
   !! 
   !! @result        Update (ua,va) with the noise advection term trend
   !!               
   !!--------------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_hadv_noi 
   ! [tlu_hadv_noi]
   SUBROUTINE tlu_hadv_noi( kt, cdtype, kcmp, pcn, pca)
      INTEGER                         , INTENT(in   ) ::   kt                     ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kcmp                   ! Define velocity type for calculation
      CHARACTER(len=1)                , INTENT(in   ) ::   cdtype                 ! =U, V or W (component indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pcn, pca               ! modified advection term for this component
      !
      INTEGER                                         ::   ji, jj, jk             ! dummy loop indices
      REAL(wp)                                        ::   zwb, zwt, zwl, zwr     ! FD scalars (bottom, top, left, right)
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_hadv_noi')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hadv_noi : ', cdtype, ' component of horizontal-noise advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !

      SELECT CASE( kcmp )     !Select velocity type to interpolate sigma dbt in to corresponding mesh point for advection calculation

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 
            !
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji= 2, jpim1
                     ! Interpolation from U -> T and adding scale factor
                     zwr =  0.5_wp * ( unoi(ji+1,jj,jk) + unoi(ji  ,jj,jk) ) * r1_e1t(ji+1,jj)        &
                         &         * (  pcn(ji+1,jj,jk) -  pcn(ji  ,jj,jk) ) *  tmask(ji+1,jj,jk)
                     !
                     zwl =  0.5_wp * ( unoi(ji  ,jj,jk) + unoi(ji-1,jj,jk) ) * r1_e1t(ji  ,jj)        &
                         &         * (  pcn(ji  ,jj,jk) -  pcn(ji-1,jj,jk) ) *  tmask(ji  ,jj,jk)
                     !          
                     pca(ji,jj,jk) = pca(ji,jj,jk) - 0.5_wp * ( zwl + zwr  ) * umask(ji,jj,jk)
                  ENDDO
               ENDDO
               DO jj = 2, jpjm1
                  DO ji= 1, jpim1
                     ! Interpolation from V -> F and adding scale factor
                     zwt =  0.5_wp * ( vnoi(ji+1,jj  ,jk) + vnoi(ji  ,jj  ,jk) ) * r1_e2f(ji,jj  )    &
                         &         * (  pcn(ji  ,jj+1,jk) -  pcn(ji  ,jj  ,jk) ) *  fmask(ji,jj  ,jk)
                     !
                     zwb =  0.5_wp * ( vnoi(ji+1,jj-1,jk) + vnoi(ji  ,jj-1,jk) ) * r1_e2f(ji,jj-1)    &
                         &         * (  pcn(ji  ,jj  ,jk) -  pcn(ji  ,jj-1,jk) ) *  fmask(ji,jj-1,jk)
                     ! 
                     pca(ji,jj,jk) = pca(ji,jj,jk) - 0.5_wp * ( zwt + zwb ) * umask(ji,jj,jk)
                  ENDDO
               ENDDO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 
            !
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji= 2, jpim1   !
                     ! Interpolation from U -> F and adding scale factor
                     zwr = 0.5_wp * ( unoi(ji  ,jj+1,jk) + unoi(ji  ,jj,jk) ) * r1_e1f(ji  ,jj)       &
                         &        * (  pcn(ji+1,jj  ,jk) -  pcn(ji  ,jj,jk) ) *  fmask(ji  ,jj,jk)
                     !
                     zwl = 0.5_wp * ( unoi(ji-1,jj+1,jk) + unoi(ji-1,jj,jk) ) * r1_e1f(ji-1,jj)       & 
                         &        * (  pcn(ji  ,jj  ,jk) -  pcn(ji-1,jj,jk) ) *  fmask(ji-1,jj,jk)
                     ! 
                     pca(ji,jj,jk) = pca(ji,jj,jk) - 0.5_wp * ( zwr + zwl ) * vmask(ji,jj,jk)
                  ENDDO
               ENDDO
               DO jj = 2, jpjm1
                  DO ji= 1, jpi   !
                     ! Interpolation from V -> T and adding scale factor
                     zwt = 0.5_wp * ( vnoi(ji,jj+1,jk) + vnoi(ji,jj  ,jk) ) * r1_e2t(ji,jj+1)         &
                         &        * (  pcn(ji,jj+1,jk) -  pcn(ji,jj  ,jk) ) *  tmask(ji,jj+1,jk)
                     !
                     zwb = 0.5_wp * ( vnoi(ji,jj  ,jk) + vnoi(ji,jj-1,jk) ) * r1_e2t(ji,jj  )         &
                         &        * (  pcn(ji,jj  ,jk) -  pcn(ji,jj-1,jk) ) *  tmask(ji,jj  ,jk) 
                     ! 
                     pca(ji,jj,jk) =  pca(ji,jj,jk) - 0.5_wp * ( zwt + zwb ) * vmask(ji,jj,jk)
                  ENDDO
               ENDDO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','tlu_hadv_noi: wrong value for kcmp'  )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_hadv_noi')   ! [NEMO] check
      !
   END SUBROUTINE tlu_hadv_noi
   ! [tlu_hadv_noi]




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
      IF( ln_timing ) CALL timing_start('tlu_hhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hhdiff : Horizontal diffusion components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))
      ALLOCATE(zw1(jpi,jpj,jpk), zw2(jpi,jpj,jpk))


      ! Velocity averages
      zcnt = 0._wp     ! For diagonal components (at T-point)
      zcnf = 0._wp     ! For extra-diagonal components (at f-point)

      ! Support vectors
      zw0 = 0._wp
      zw1 = 0._wp
      zw2 = 0._wp


      SELECT CASE( kcmp )

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity 
            !
            ! Interpolation of the velocity component
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 2, jpi   
                     zcnt(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji-1,jj,jk) ) * tmask(ji,jj,jk)
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, jpi   
                     zcnf(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj+1,jk) + pcn(ji,jj,jk) ) * fmask(ji,jj,jk)  
                  END DO
               END DO
            END DO
            ! Both averages checkd
            !
            ! Lateral boundary condition transfer across nodes
            CALL lbc_lnk_multi( 'tlu_hhdiff', zcnt , 'T', 1., zcnf , 'F', 1. )
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
                     zwb =  ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * zcnt(ji  ,jj,jk)   &
                          & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia11) * zcnt(ji-1,jj,jk) ) &
                          & * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk)
                     ! After  (   (i+1) -   i     )
                     zwa =  ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11) * zcnt(ji+1,jj,jk)   &
                          & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * zcnt(ji  ,jj,jk) ) &
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
                     zwb = ( e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * zcnt(ji,jj  ,jk)   &
                         & - e1t(ji,jj-1) * e3t_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia11) * zcnt(ji,jj-1,jk) ) &
                         & * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk)
                     ! After (   (j+1) -   j     )
                     zwa = ( e1t(ji,jj+1) * e3t_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia11) * zcnt(ji,jj+1,jk)   &
                         & - e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * zcnt(ji,jj  ,jk) ) &
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
     
            CALL lbc_lnk_multi( 'tlu_hhdiff', zw0 , 'T', 1.)

            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji = 1, jpim1
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca(ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) + 0.5_wp * umask(ji,jj,jk)        &
                                       & * ( zw0 (ji+1,jj  ,jk  ) + zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity 
            !
            ! Interpolation of the velocity component
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpim1  
                     zcnf(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji+1,jj,jk) ) * fmask(ji,jj,jk)
                  END DO
               END DO
               DO jj = 2, jpj
                  DO ji = 1, jpi   
                     zcnt(ji,jj,jk) = 0.5_wp * ( pcn(ji,jj,jk) + pcn(ji,jj-1,jk) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            ! Both averages checkd
            !
            ! Lateral boundary condition transfer across nodes
            !
            CALL lbc_lnk_multi( 'tlu_hhdiff', zcnt , 'T', 1., zcnf , 'F', 1. )
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
                     zwb =  ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * zcnt(ji  ,jj,jk)   &
                          & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia11) * zcnt(ji-1,jj,jk) ) &
                          & * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk)
                     ! After  (   (i+1) -   i     )
                     zwa =  ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11) * zcnt(ji+1,jj,jk)   &
                          & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) * zcnt(ji  ,jj,jk) ) &
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
                     zwb = ( e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * zcnt(ji,jj  ,jk)   &
                         & - e1t(ji,jj-1) * e3t_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia11) * zcnt(ji,jj-1,jk) ) &
                         & * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk)
                     ! After (   (j+1) -   j     )
                     zwa = ( e1t(ji,jj+1) * e3t_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia11) * zcnt(ji,jj+1,jk)   &
                         & - e1t(ji,jj  ) * e3t_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia11) * zcnt(ji,jj  ,jk) ) &
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
     
            CALL lbc_lnk_multi( 'tlu_hhdiff', zw0 , 'T', 1.)

            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     ! Averaging and update of the velocity component
                     ! This is a diffusion component added to the RIGHT-HAND SIDE, so positive
                     pca(ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  )  + 0.5_wp * vmask(ji,jj,jk)   &
                                       & * ( zw0 (ji  ,jj+1,jk  ) + zw0 (ji  ,jj  ,jk  ) )
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




   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_hhdrift ***
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
   ! @snippet this tlu_hhdrift
   ! [tlu_hhdrift]
   SUBROUTINE tlu_hhdrift( kt, cdtype, kcmp , pcn, pca)
      INTEGER                           , INTENT(in   ) ::   kt                    ! ocean time-step index
      INTEGER                           , INTENT(in   ) ::   kcmp                  ! Define velocity type for calculation
      CHARACTER(len=1)                  , INTENT(in   ) ::   cdtype                ! =U or V (field indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(in   ) ::   pcn                   ! before field
      REAL(wp), DIMENSION(jpi,jpj,jpk)  , INTENT(inout) ::   pca                   ! field trend
      !
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
      IF( ln_timing ) CALL timing_start('tlu_hhdrift')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_hhdrift : Horizontal drift components on velocity type ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF


      ALLOCATE(zw0(jpi,jpj,jpk))
      ALLOCATE(zw1(jpi,jpj,jpk), zw2(jpi,jpj,jpk))


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
                     ! Averaging and update of the velocity component
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     pca(ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  ) - 0.5_wp * umask(ji,jj,jk)        &
                                       & * pcn(ji  ,jj  ,jk  ) * ( zw0 (ji+1,jj  ,jk  )          &
                                       &                         + zw0 (ji  ,jj  ,jk  ) )
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
                     ! Averaging and update of the velocity component
                     ! This is a drift component added to the RIGHT-HAND SIDE, so negative
                     pca(ji  ,jj  ,jk  ) = pca(ji  ,jj  ,jk  )  - 0.5_wp * vmask(ji,jj,jk)       &
                                       & * pcn(ji  ,jj  ,jk  ) * ( zw0 (ji  ,jj+1,jk  )          &
                                       &                         + zw0 (ji  ,jj  ,jk  ) )
                  END DO  
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         CASE DEFAULT          
           CALL ctl_stop('STOP','tlu_hhdrift: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zw0, zw1, zw2)

   END SUBROUTINE tlu_hhdrift
   ! [tlu_hhdrift]


END MODULE tluhdyn

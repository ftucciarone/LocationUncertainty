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
MODULE tlubias
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
   PUBLIC tlu_adv_bias  ! Called by tlu_bcdyn in tlu_step.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS


   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_adv_bias ***
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
   !> @snippet this tlu_adv_bias 
   ! [tlu_adv_bias]
   SUBROUTINE tlu_adv_bias( kt, cdtype, kcmp, pcn, pca, numsign)
      INTEGER                         , INTENT( in  ) ::   kt                     ! ocean time-step index
      INTEGER                         , INTENT( in  ) ::   kcmp                   ! Define velocity type for calculation
      CHARACTER(len=1)                , INTENT( in  ) ::   cdtype                 ! =U, V or W (component indicator)
      REAL(wp)                        , INTENT( in  ) ::   numsign                ! Bias sign
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pcn, pca               ! modified advection term for this component
      !
      INTEGER                                         ::   ji, jj, jk             ! dummy loop indices
      REAL(wp)                                        ::   zwb, zwt               ! FD scalars (bottom, top)     
      REAL(wp)                                        ::   zwl, zwr               ! FD scalars (left, right)           
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_zadv_noi')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_adv_bias :', cdtype, ' component of vertical-noise advection (no surface layer)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !

      SELECT CASE( kcmp )     

         !> @par               U Velocity
         CASE ( np_ucmp ) !For U velocity
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji= 2, jpim1
                     ! Interpolation from U -> T and adding scale factor
                     zwr =  0.5_wp * ( ubia(ji+1,jj,jk) + ubia(ji  ,jj,jk) ) * r1_e1t(ji+1,jj)        &
                         &         * (  pcn(ji+1,jj,jk) -  pcn(ji  ,jj,jk) ) *  tmask(ji+1,jj,jk)
                     !
                     zwl =  0.5_wp * ( ubia(ji  ,jj,jk) + ubia(ji-1,jj,jk) ) * r1_e1t(ji  ,jj)        &
                         &         * (  pcn(ji  ,jj,jk) -  pcn(ji-1,jj,jk) ) *  tmask(ji  ,jj,jk)
                     !          
                     pca(ji,jj,jk) = pca(ji,jj,jk) + 0.5_wp * numsign * ( zwl + zwr  ) * umask(ji,jj,jk)
                  ENDDO
               ENDDO
               DO jj = 2, jpjm1
                  DO ji= 1, jpim1
                     ! Interpolation from V -> F and adding scale factor
                     zwt =  0.5_wp * ( vbia(ji+1,jj  ,jk) + vbia(ji  ,jj  ,jk) ) * r1_e2f(ji,jj  )    &
                         &         * (  pcn(ji  ,jj+1,jk) -  pcn(ji  ,jj  ,jk) ) *  fmask(ji,jj  ,jk)
                     !
                     zwb =  0.5_wp * ( vbia(ji+1,jj-1,jk) + vbia(ji  ,jj-1,jk) ) * r1_e2f(ji,jj-1)    &
                         &         * (  pcn(ji  ,jj  ,jk) -  pcn(ji  ,jj-1,jk) ) *  fmask(ji,jj-1,jk)
                     ! 
                     pca(ji,jj,jk) = pca(ji,jj,jk) + 0.5_wp * numsign * ( zwt + zwb ) * umask(ji,jj,jk)
                  ENDDO
               ENDDO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            !
            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpj
                  DO ji= 1, jpim1
                     !                   Top        (   (k+1) -   k      )
                     zwt =  0.5_wp * ( wbia(ji+1,jj,jk-1) + wbia(ji,jj,jk-1) ) * wumask(ji,jj,jk-1) &
                         &         * (  pcn(ji  ,jj,jk-1) -  pcn(ji,jj,jk  ) ) / e3uw_n(ji,jj,jk-1)
                     !                  Below       (     k   - (k-1)    )
                     zwb =  0.5_wp * ( wbia(ji+1,jj,jk  ) + wbia(ji,jj,jk  ) ) * wumask(ji,jj,jk  ) &
                         &         * (  pcn(ji  ,jj,jk  ) -  pcn(ji,jj,jk+1) ) / e3uw_n(ji,jj,jk  )
                     ! Average in z;   0.5 * ( Top + Below )
                     pca(ji,jj,jk) = pca(ji,jj,jk) + 0.5_wp * numsign * ( zwt + zwb ) * umask(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================

         !> @par               V Velocity
         CASE ( np_vcmp ) !For V velocity
            !                                ! ================
            DO jk = 1, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji= 2, jpim1   !
                     ! Interpolation from U -> F and adding scale factor
                     zwr = 0.5_wp * ( ubia(ji  ,jj+1,jk) + ubia(ji  ,jj,jk) ) * r1_e1f(ji  ,jj)       &
                         &        * (  pcn(ji+1,jj  ,jk) -  pcn(ji  ,jj,jk) ) *  fmask(ji  ,jj,jk)
                     !
                     zwl = 0.5_wp * ( ubia(ji-1,jj+1,jk) + ubia(ji-1,jj,jk) ) * r1_e1f(ji-1,jj)       & 
                         &        * (  pcn(ji  ,jj  ,jk) -  pcn(ji-1,jj,jk) ) *  fmask(ji-1,jj,jk)
                     ! 
                     pca(ji,jj,jk) = pca(ji,jj,jk) + 0.5_wp * numsign * ( zwr + zwl ) * vmask(ji,jj,jk)
                  ENDDO
               ENDDO
               DO jj = 2, jpjm1
                  DO ji= 1, jpi   !
                     ! Interpolation from V -> T and adding scale factor
                     zwt = 0.5_wp * ( vbia(ji,jj+1,jk) + vbia(ji,jj  ,jk) ) * r1_e2t(ji,jj+1)         &
                         &        * (  pcn(ji,jj+1,jk) -  pcn(ji,jj  ,jk) ) *  tmask(ji,jj+1,jk)
                     !
                     zwb = 0.5_wp * ( vbia(ji,jj  ,jk) + vbia(ji,jj-1,jk) ) * r1_e2t(ji,jj  )         &
                         &        * (  pcn(ji,jj  ,jk) -  pcn(ji,jj-1,jk) ) *  tmask(ji,jj  ,jk) 
                     ! 
                     pca(ji,jj,jk) =  pca(ji,jj,jk) + 0.5_wp * numsign * ( zwt + zwb ) * vmask(ji,jj,jk)
                  ENDDO
               ENDDO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================
            !
            !                                ! ================
            DO jk = 2, jpkm1                 ! Horizontal slab
               !                             ! ================
               DO jj = 1, jpjm1
                  DO ji= 1, jpi
                     !                   Top        (   (k+1) -   k      )
                     zwt =  0.5_wp * ( wbia(ji,jj+1,jk-1) + wbia(ji,jj,jk-1) ) * wvmask(ji,jj,jk-1) &
                         &         * (  pcn(ji,jj  ,jk-1) -  pcn(ji,jj,jk  ) ) / e3vw_n(ji,jj,jk-1) 
                     !                  Below       (     k   - (k-1)    )
                     zwb =  0.5_wp * ( wbia(ji,jj+1,jk  ) + wbia(ji,jj,jk  ) ) * wvmask(ji,jj,jk) &
                         &         * (  pcn(ji,jj  ,jk  ) -  pcn(ji,jj,jk+1) ) / e3vw_n(ji,jj,jk  ) 
                     ! Average in z;   0.5 * ( Top + Below )
                     pca(ji,jj,jk) = pca(ji,jj,jk) + 0.5_wp * numsign * ( zwt + zwb ) * vmask(ji,jj,jk)
                  END DO
               END DO
               !                             ! ================
            END DO                           !   End of slab
            !                                ! ================


         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','tlu_adv_bias: wrong value for kcmp'  )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_adv_bias')   ! [NEMO] check
      !
   END SUBROUTINE tlu_adv_bias
   ! [tlu_adv_bias]


END MODULE tlubias

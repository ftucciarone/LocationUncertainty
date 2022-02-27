!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tluwzv
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
MODULE tluwzv   
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE tlu              ! Initialization of stochastic structures
   ! [mod_dep]

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tlu_isd     ! called by step.F90
   PUBLIC   tlu_wzvcmp  ! called by step.F90

CONTAINS

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_isd ***
   !!
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
   !> @details
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
   ! @snippet this tlu_isd
   ! [tlu_isd]
   SUBROUTINE tlu_isd
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      uisd_n = 0._wp
      visd_n = 0._wp
      wisd_n = 0._wp
      !                                ! ================
      DO jk = 1, jpkm1                 ! Horizontal slab
         !                             ! ================
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               ! First component
               !
               ! dx(a11)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) + r1_e1e2u(ji,jj) *                           & 
                              & ( e2t(ji+1,jj) * e3t_n(ji+1,jj,jk) * var_ten(ji+1,jj,jk,ia11)    &
                              & - e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia11) )  &
                              & * umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               ! dy(a21)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) + r1_e1e2u(ji,jj) *                           & 
                              & ( e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia12)    &
                              & - e1f(ji,jj-1) * e3f_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia12) )  &
                              & * umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               ! dz(a31)
               uisd_n(ji,jj,jk) = uisd_n(ji,jj,jk) +                                             & 
                              & ( var_ten(ji,jj,jk  ,ia13) - var_ten(ji,jj,jk+1,ia13) ) *        &
                              &   umask(ji,jj,jk) / e3u_n(ji,jj,jk)
               !
               ! Second component
               !
               ! dx(a12)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) + r1_e1e2v(ji,jj) *                           & 
                              & ( e2t(ji  ,jj) * e3t_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia12)    &
                              & - e2t(ji-1,jj) * e3t_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia12) )  &
                              & * vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               ! dy(a22)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) + r1_e1e2v(ji,jj) *                           & 
                              & ( e1f(ji,jj+1) * e3f_n(ji,jj+1,jk) * var_ten(ji,jj+1,jk,ia22)    &
                              & - e1f(ji,jj  ) * e3f_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia22) )  &
                              & * vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               ! dz(a32)
               visd_n(ji,jj,jk) = visd_n(ji,jj,jk) +                                             & 
                              & ( var_ten(ji,jj,jk  ,ia23) - var_ten(ji,jj,jk+1,ia23) ) *        &
                              &   vmask(ji,jj,jk) / e3v_n(ji,jj,jk)
               !
               ! Third component
               !
               ! dx(a13)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) + r1_e1e2t(ji,jj) *                           & 
                              & ( e2t(ji  ,jj) * e3uw_n(ji  ,jj,jk) * var_ten(ji  ,jj,jk,ia13)   &
                              & - e2t(ji-1,jj) * e3uw_n(ji-1,jj,jk) * var_ten(ji-1,jj,jk,ia13) ) &
                              & * wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
               ! dy(a23)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) + r1_e1e2t(ji,jj) *                           & 
                              & ( e1t(ji,jj  ) * e3vw_n(ji,jj  ,jk) * var_ten(ji,jj  ,jk,ia23)   &
                              & - e1t(ji,jj-1) * e3vw_n(ji,jj-1,jk) * var_ten(ji,jj-1,jk,ia23) ) &
                              & * wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            !
            ! Third component Surface boundary condition
            !
            ! dz(a33)
            wisd_n(ji,jj, 1) =  wisd_n(ji,jj,1)                                                  & 
                           & - var_ten(ji,jj,1,ia33) * wmask(ji,jj,1) / e3w_n(ji,jj,1)
         END DO
      END DO

      !                                ! ================
      DO jk = 2, jpkm1                 ! Horizontal slab
         !                             ! ================
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               ! Third component
               !
               ! dz(a33)
               wisd_n(ji,jj,jk) = wisd_n(ji,jj,jk) +                                            & 
                              & ( var_ten(ji,jj,jk  ,ia33) - var_ten(ji,jj,jk+1,ia33) ) *       &
                              &   wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_hhdiff', uisd_n , 'U', 1., visd_n , 'V', 1., wisd_n , 'T', 1.  )
      !
      ! Print Ito-Stokes Drift
      !
      IF( kt == nit000 )  THEN
         CALL iom_put( 'U_isd', uisd_n )  
         CALL iom_put( 'V_isd', visd_n )  
         CALL iom_put( 'W_isd', wisd_n )  
      END IF
      !
   END SUBROUTINE tlu_isd
   ! [tlu_isd]

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_wzvcmp ***
   !!
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
   !> @details
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
   ! @snippet this tlu_wzvcmp
   ! [tlu_wzvcmp]
   SUBROUTINE tlu_wzvcmp
      REAL(wp), DIMENSION(jpi,jpj,jpk)                  ::   zwz         ! 2D workspace
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !
      !!----------------------------------------------------------------------
      !
      zwz = 0._wp
      tlu_wcorr = 0._wp
      !                                ! ================
      DO jk = 1, jpkm1                 ! Horizontal slab
         !                             ! ================
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               !
               ! First component
               !
               ! dx( udiva )
               zwz(ji,jj,jk) =   zwz(ji,jj,jk) + r1_e1e2t(ji,jj)  *                           & 
                             & ( e2u(ji+1,jj) * e3u_n(ji+1,jj,jk) * uisd_n(ji+1,jj,jk)        &
                             & - e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * uisd_n(ji  ,jj,jk) ) *    &
                             &   tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
               !
               ! Second component
               !
               ! dy( vdiva )
               zwz(ji,jj,jk) =   zwz(ji,jj,jk) + r1_e1e2t(ji,jj)  *                           & 
                             & ( e1v(ji,jj+1) * e3v_n(ji,jj+1,jk) * visd_n(ji,jj+1,jk)        &
                             & - e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * visd_n(ji,jj  ,jk) ) *    &
                             &   tmask(ji,jj,jk) / e3t_n(ji,jj,jk)
               !
               ! Third component
               !
               ! dz( wdiva )
               zwz(ji,jj,jk) =   zwz(ji,jj,jk) +                                  & 
                             & ( wisd_n(ji,jj,jk  ) - wisd_n(ji,jj,jk+1) )  *     &
                             &   wmask(ji,jj,jk) / e3w_n(ji,jj,jk)
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      DO jk = jpkm1, 1, -1
         tlu_wcorr(:,:,jk) = tlu_wcorr(:,:,jk+1) + 0.5_wp * e3t_n(:,:,jk) * zwz(:,:,jk) 
         ! computation of w
         wn(:,:,jk) = wn(:,:,jk) - tlu_wcorr(:,:,jk)
      END DO
      !
   END SUBROUTINE tlu_wzvcmp
   ! [tlu_wzvcmp]


END MODULE tluwzv


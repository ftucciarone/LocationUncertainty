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

   PUBLIC   tlu_wzvcmp  ! called by step.F90

CONTAINS


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
            END DO
         END DO
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      DO jk = jpkm1, 1, -1
         tlu_wcorr(:,:,jk) = tlu_wcorr(:,:,jk+1) + e3t_n(:,:,jk) * zwz(:,:,jk) 
         ! computation of w
         wn(:,:,jk) = wn(:,:,jk) + tlu_wcorr(:,:,jk) + wisd_n(:,:,jk)
      END DO
      !
   END SUBROUTINE tlu_wzvcmp
   ! [tlu_wzvcmp]


END MODULE tluwzv


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
MODULE tluwzvDYN   
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
      !
      INTEGER                                           ::   ji, jj, jk            ! dummy loop indices
      !
      !!----------------------------------------------------------------------
      !
      tlu_wcorr = 0._wp
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         tlu_wcorr( 2:jpi-1, 2:jpj-1, jk) = tlu_wcorr( 2:jpi-1, 2:jpj-1, jk + 1) +    &

         &                                  (                                         &
         !  dx( dy*dz*u_isd )
         &                                    (    e2u( 2:jpi-1, 2:jpj-1 )            &
         &                                  *    e3u_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                                  *   uisd_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                                  -      e2u( 1:jpi-2, 2:jpj-1 )            &
         &                                  *    e3u_n( 1:jpi-2, 2:jpj-1, jk )        &
         &                                  *   uisd_n( 1:jpi-2, 2:jpj-1, jk )        &
         !  dy ( dx*dz*v_isd )
         &                                  +      e1v( 2:jpi-1, 2:jpj-1 )            &
         &                                  *    e3v_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                                  *   visd_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                                  -      e1v( 2:jpi-1, 1:jpj-2 )            &
         &                                  *    e3v_n( 2:jpi-1, 1:jpj-2, jk )        &
         &                                  *   visd_n( 2:jpi-1, 1:jpj-2, jk ) )      &
         !  .../dV
         &                                  *    tmask( 2:jpi-1, 2:jpj-1, jk)         & 
         &                                  * r1_e1e2t( 2:jpi-1, 2:jpj-1 )            & 
         &                                  /    e3t_n( 2:jpi-1, 2:jpj-1, jk)         &
         !  + dz( w_isd )
         &                                  + ( wisd_n( 2:jpi-1, 2:jpj-1, jk )        &
         &                                  -   wisd_n( 2:jpi-1, 2:jpj-1, jk + 1 ) )  &
         &                                  *    tmask( 2:jpi-1, 2:jpj-1, jk)         & 
         &                                  /    e3t_n( 2:jpi-1, 2:jpj-1, jk) )       &
         !
         &                                     * e3t_n( 2:jpi-1, 2:jpj-1, jk )
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_wzvcmp', tlu_wcorr , 'T', -1.  )
      wn = wn - tlu_wcorr
      !
   END SUBROUTINE tlu_wzvcmp   
   ! [tlu_wzvcmp]


END MODULE tluwzvDYN


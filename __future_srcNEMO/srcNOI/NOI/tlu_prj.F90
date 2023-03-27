!!===========================================================================
!!                       ***  MODULE  tlusvd  ***
!! Transport under Location Uncertainty: stochastic noise SVD
!!                       noise generation based on pseudo-observations and SVD.
!!===========================================================================
!! History : 0.0  !  2017-..  (P. DERIAN)  Original code
!!
!! [TODO]    - read namelist from _cfg, not only _ref
!!           - write error messages to appropriate unit
!!           - use NEMO's work arrays?
!!---------------------------------------------------------------------------
MODULE tluprj
   ! [mod_dep]
   USE par_kind           ! data types defined in par_kind module
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce            ! ocean dynamics parameters
   USE oce                ! ocean dynamics and active tracers
   USE dom_oce            ! ocean space and time domain
   USE lib_mpp            ! MPP library
   USE timing             ! Timing
   USE lbclnk             ! ocean lateral boundary conditions (or mpp link)
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_vel_proj    !
   PUBLIC tlu_var_proj
   !
   INCLUDE 'mpif.h'
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_vel_proj ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_vel_proj
   ! [tlu_vel_proj]
   SUBROUTINE tlu_vel_proj( uin, vin, win, pcu, pcv, pcw )
      USE ldfslp, ONLY : uslp, vslp, wslpi, wslpj
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE Create one realisation of noise  ***
      !!
      !! ** Purpose :   Formulate noise from spatial modes and time coeff/random numbers
      !!
      !! ** Method  :   Project spatial modes on corresponding temporal coefficient or
      !!                on gaussian random number as provided by input.
      !!                The noise is divided by the time step because later on in the 
      !!                NEMO deterministic kernel it will multiplied by rn_dt
      !!------------------------------------------------------------------------
      INTEGER                                        :: ji, jj, jk, ierr 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   uin                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   vin                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   win                   ! before field (so its euler)

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwa1,   zwa2,   zwa3

      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pcu                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pcv                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pcw                   ! before field (so its euler)
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)         ::   zwr_1, zwr_2, zur_1, zvr_2
      !
      ! Allocate dummy matrices
      !
      ALLOCATE( zwa1(jpi,jpj,jpk), &
      &         zwa2(jpi,jpj,jpk), &
      &         zwa3(jpi,jpj,jpk),  stat=ierr ) 
      !
      ! Safety initialization of output
      !
      zwa1 = 0._wp
      zwa2 = 0._wp
      zwa3 = 0._wp
      !
      ! Allocate dummy matrices
      !
      ALLOCATE( zwr_1(jpi,jpj,jpk), &
      &         zwr_2(jpi,jpj,jpk), &
      &         zur_1(jpi,jpj,jpk), &
      &         zvr_2(jpi,jpj,jpk),  stat=ierr ) 
      !
      ! Initialize Support Arrays
      !
      zwr_1 = wslpi * win
      zwr_2 = wslpj * win
      zur_1 =  uslp * uin
      zvr_2 =  vslp * vin
      !
      zwa3 = zwr_1 + zwr_2
      !
      zwa1(1:jpim1,2:jpj-1,1:jpk-1) =       uin(1:jpim1,1:jpjm1,1:jpkm1) -                                &
      &                                 ( zwr_1(2:jpi  ,1:jpj  ,1:jpkm1) + zwr_1(1:jpim1,1:jpj  ,1:jpkm1) &
      &                                 + zwr_1(2:jpi  ,1:jpj  ,2:jpk  ) + zwr_1(1:jpim1,1:jpj  ,2:jpk  ) ) * 0.25_wp
      !
      zwa2(2:jpi-1,2:jpj-1,1:jpk-1) =       vin(1:jpi  ,1:jpjm1,1:jpkm1) -                                &
      &                                 ( zwr_2(1:jpi  ,2:jpj  ,1:jpkm1) + zwr_2(1:jpi  ,1:jpjm1,1:jpkm1) &
      &                                 + zwr_2(1:jpi  ,2:jpj  ,2:jpk  ) + zwr_2(1:jpi  ,1:jpjm1,2:jpk  ) ) * 0.25_wp
      !
      zwa3(2:jpi-1,2:jpj-1,2:jpk-1) =     zwr_1(2:jpi-1,2:jpj-1,2:jpk-1) * wslpi(2:jpi-1,2:jpj-1,2:jpk-1) &
      &                                 + zwr_2(2:jpi-1,2:jpj-1,2:jpk-1) * wslpj(2:jpi-1,2:jpj-1,2:jpk-1) &
      !                              
      !                               Interpolation of (r_1 * u) in  and z
      !                              
      &                               - ( zur_1(2:jpi-1,2:jpj-1,1:jpk-2) + zur_1(1:jpi-2,2:jpj-1,1:jpk-2) &
      &                                 + zur_1(2:jpi-1,2:jpj-1,2:jpk-1) + zur_1(1:jpi-2,2:jpj-1,2:jpk-1) &
      !                              
      !                               Interpolation of (r_2 * v) in  and z
      !                              
      &                                 + zvr_2(2:jpi-1,2:jpj-1,1:jpk-2) + zvr_2(2:jpi-1,1:jpj-2,1:jpk-1) &
      &                                 + zvr_2(2:jpi-1,2:jpj-1,2:jpk-1) + zvr_2(2:jpi-1,1:jpj-2,2:jpk-1) ) * 0.25_wp 
      !
      ! Re-imposition of bottom boundary no slip condition
      !
      zwa3(:,:, 1 ) = 0._wp
      zwa3(:,:,jpk) = 0._wp
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_noi_proj', zwa1 , 'U', -1., zwa2 , 'V', -1., zwa3 , 'T', -1.  )
      !
      ! Mask the noise to enforce No-Flux boundary conditions 
      !
      pcu = zwa1 * umask
      pcv = zwa2 * vmask
      pcw = zwa3 * wmask
      !
      DEALLOCATE( zwa1, zwa2, zwa3, zwr_1, zwr_2, zur_1, zvr_2 ) 
      !
   END SUBROUTINE tlu_vel_proj
   ! [tlu_vel_proj]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_var_proj ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this tlu_var_proj
   ! [tlu_var_proj]
   SUBROUTINE tlu_var_proj(  )
      USE ldfslp, ONLY : uslp, vslp, wslpi, wslpj
      !
   END SUBROUTINE tlu_var_proj 
   ! [tlu_var_proj]

END MODULE tluprj


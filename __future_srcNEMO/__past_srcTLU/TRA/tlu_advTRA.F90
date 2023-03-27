!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tluhdifTRA
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!!
!> @authors
!>                P. Derian, P. Chandramouli, F. Tucciarone
!!
!> @version
!!                2017 -  ( P. DERIAN       )  Original code <br>
!!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
!!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Transport under Location Uncertainty: TRACER diffusion.
!! 
!! @par           Procedure specifics      
!>
!!                \f{align*}{
!!                {\color{gray}  \mathrm{d}_{t} \phi } + &
!!                {\color{gray}  \nabla\cdot\left\lbrace \left[\boldsymbol{u}-\boldsymbol{u}_{s}+\boldsymbol{\mu}_{t}\right)\phi\mathrm{d}t
!!                                 +  \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}\phi \right\rbrace } - 
!!                {\color{gray}  \dfrac{1}{2}\nabla\cdot\left( \boldsymbol{a}\nabla\phi \right)\mathrm{d}t } \\ &+
!!                {\color{black} \boldsymbol{f}\times } \left( 
!!                {\color{gray}  \boldsymbol{u}\mathrm{d}t + }
!!                {\color{black} \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }
!!                              \right) = 
!!                {\color{gray}  - \nabla p\mathrm{d}t }
!!                {\color{black} - \nabla\mathrm{d}p_{t}^{\sigma} } 
!!                {\color{gray}  + \mathrm{D}^{\phi} }
!!                {\color{gray}  + \mathrm{F}^{\phi} } 
!!                       
!!                \f}

!> @par           Code specifics
!!                This module is called by tlu_bcdyn in tlu.f90. The terms @f$ F_{x} @f$ and @f$ F_{y} @f$ are computed in \ref tlu_trahhdiff,
!!                while the term @f$ F_{z} @f$ (which does not contain the purely vertical diffusion term) is computed in \ref tlu_trahzdiff,
!!                so that in terms of NEMO operators it becomes
!!                \f{align*}{
!!                \nabla\cdot\boldsymbol{F}_{T} = \underbrace{\dfrac{1}{e_{1t}e_{2t}e_{3t}}\left\lbrace 
!!                       \delta_{i}\left[ F_{T,x} \right]
!!                     + \delta_{j}\left[ F_{T,y} \right] 
!!                       \right\rbrace}_{\texttt{tlu\_trahhdiff}}
!!                     + \underbrace{\dfrac{1}{e_{3t}}\delta_{k}\left[ F_{T,z} \right]}_{\texttt{tlu\_trahzdiff}}
!!                \f}
 
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
!> @todo          Write error messages to appropriate unit
!!
!!-----------------------------------------------------------------------------------------------------------------------------------

MODULE tluadvTRA
   !!===========================================================================
   !!                       ***  MODULE  tlutra  ***
   !! Transport under Location Uncertainty: Tracer contributions of the stochastic
   !! parametrization of N-S equations.
   !!
   !!===========================================================================
   !! History :         (P. CHANDRAMOULI)  Original code
   !!
   !! [TODO]    - write error messages to appropriate unit
   !!
   !!---------------------------------------------------------------------------
    USE par_kind         ! data types defined in par_kind module
    USE in_out_manager   ! I/O manager
    USE iom              ! I/O manager library
    USE par_oce
    USE oce              ! ocean dynamics and active tracers
    USE dom_oce          ! ocean space and time domain
    USE lib_mpp         ! MPP library
    USE timing           ! Timing
    USE tlu
   !
   IMPLICIT NONE   ! turn off implicit variable declaration
   PRIVATE         ! and make stuff private by default
   !
   PUBLIC tlu_tra_adv_noi               ! Called by tlu_tra_adv in tlu_fluxTRA.f90
   !
   !! * Substitutions [TODO] enable substitutions for NEMO
! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------
CONTAINS

   SUBROUTINE tlu_tra_adv_noi( kt, cdtype, ptn, pta, uadv, vadv, wadv, r_Scale)
      USE tlusbcTRA
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tlu_adv_sto  ***
      !!
      !! ** Purpose :   Compute the stochastic advection velocity -0.5*div(a)'
      !!                for the given component.
      !!
      !! ** Method  :   Set pc to -0.5*divc
      !!                where divc = 1/(e1c e2c)*(di[e2c a1] + dj[e1c a2])
      !!                           + 1/(e3c)*dk[a3]
      !!                and a1, a2, a3 are the appropriate components of the
      !!                stochastic diffusion tensor var_ten.
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER,                               INTENT(in   ) ::   kt             ! ocean time-step index
      CHARACTER(len=1),                      INTENT(in   ) ::   cdtype         ! =U, V or W (component indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk),      INTENT(in   ) ::   uadv, vadv, wadv !The components advecting the velocity field
      REAL(wp), DIMENSION(jpts)            , INTENT(in   ) ::   r_Scale               ! Schmidt number      
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   ptn, pta             ! modified advection term for this component
      !
      INTEGER  ::   ji, jj, jk, jn                               ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk)                     ::   div           ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwa1, zwa2, zwa3, zwpw   ! workspace
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_tra_adv_noi')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'tlu_tra_adv_noi : Noise adv. for tracer'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF
      !

      ALLOCATE( zwpw(jpi,jpj,jpk),   &
              & zwa1(jpi,jpj,jpk),   &
              & zwa2(jpi,jpj,jpk),   &
              & zwa3(jpi,jpj,jpk)    )
      !
      zwa1 = 0._wp
      zwa2 = 0._wp
      zwa3 = 0._wp
      !
      zwa1(:,:,1:jpkm1) = uadv(:,:,1:jpkm1) / spread(e1u,3,jpkm1)        !adding scale factor
      zwa2(:,:,1:jpkm1) = vadv(:,:,1:jpkm1) / spread(e2v,3,jpkm1)        !adding scale factor
      zwa3(:,:,2:jpkm1) = wadv(:,:,2:jpkm1) / e3w_n(:,:,2:jpkm1)           !adding scale factor
      !
      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !  
         !
         ! Calculate vertical advection separately due to surface issues in direct calculation
         !
         zwpw = 0._wp
         zwpw(:,:,2:jpkm1) = zwa3(:,:,2:jpkm1) * ( ptn(:,:,1:jpkm1-1,jn) - ptn(:,:,2:jpkm1,jn) ) &
         &                 + zwa3(:,:,3:jpk  ) * ( ptn(:,:,2:jpkm1  ,jn) - ptn(:,:,3:jpk  ,jn) )
	 !
	 ! Surface boundary conditions
	 !
	 CALL tlu_traadv_sbc( kt, jn, zwpw(:,:,1) )
         !
         !
         ! Add advection by sigma dbt to the after trends of momentum through direction calculation in horizontal and pre-calculated in vertical.
         !
         pta(2:jpim1,2:jpjm1,1:jpkm1,jn) = pta(2:jpim1,2:jpjm1,1:jpkm1,jn) - 0.5_wp * ( &
         !
         ! x-Directed gradient with interpolation in x
         !
         &                                 zwa1(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( ptn(3:jpi  ,2:jpjm1,1:jpkm1,jn) - ptn(2:jpim1  ,2:jpjm1,1:jpkm1,jn) ) &
         &                               + zwa1(1:jpim1-1,2:jpjm1  ,1:jpkm1) * ( ptn(2:jpim1,2:jpjm1,1:jpkm1,jn) - ptn(1:jpim1-1,2:jpjm1,1:jpkm1,jn) ) &
         !
         ! y-Directed gradient with interpolation in y
         !
         &                               + zwa2(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( ptn(2:jpim1,3:jpj  ,1:jpkm1,jn) - ptn(2:jpim1,2:jpjm1  ,1:jpkm1,jn) ) &
         &                               + zwa2(2:jpim1  ,1:jpjm1-1,1:jpkm1) * ( ptn(2:jpim1,2:jpjm1,1:jpkm1,jn) - ptn(2:jpim1,1:jpjm1-1,1:jpkm1,jn) ) &
         !
         ! z-Directed gradient
         !
         &                               + zwpw(2:jpim1,2:jpjm1,1:jpkm1) ) * tmask(2:jpim1,2:jpjm1,1:jpkm1) / r_Scale(jn)
         !
      END DO
      !
      !
      zwa1 = 0._wp
      zwa2 = 0._wp
      zwa3 = 0._wp
      !
      zwa1(:,:,1:jpkm1) = uadv(:,:,1:jpkm1) * spread(  e2u,3,jpkm1) * e3u_n(:,:,1:jpkm1)         !adding scale factor
      zwa2(:,:,1:jpkm1) = vadv(:,:,1:jpkm1) * spread(  e1v,3,jpkm1) * e3v_n(:,:,1:jpkm1)         !adding scale factor
      zwa3(:,:,2:jpkm1) = wadv(:,:,2:jpkm1) * spread(e1e2t,3,jpkm1)                          !adding scale factor
      !
      div = 0._wp
      !
      div(2:jpim1,2:jpjm1,1:jpkm1) = ( zwa1(2:jpim1,2:jpjm1,1:jpkm1  ) - zwa1(1:jpim1-1,2:jpjm1  ,1:jpkm1) &
      &                              + zwa2(2:jpim1,2:jpjm1,1:jpkm1  ) - zwa2(2:jpim1  ,1:jpjm1-1,1:jpkm1) &
      &                              + zwa3(2:jpim1,2:jpjm1,1:jpkm1-1) - zwa3(2:jpim1  ,2:jpjm1  ,2:jpkm1) )
      !
      div(2:jpim1,2:jpjm1,1:jpkm1) = div(2:jpim1,2:jpjm1,1:jpkm1)  * spread(r1_e1e2t(2:jpim1,2:jpjm1),3,jpk-1) / e3t_n(2:jpim1,2:jpjm1,1:jpkm1) 

      CALL iom_put( "tlu_noidiv", div )

      DEALLOCATE(zwa1, zwa2, zwa3, zwpw)
      !
      IF( ln_timing ) CALL timing_stop('tlu_tra_adv_noi')   ! [NEMO] check
      !
   END SUBROUTINE tlu_tra_adv_noi

END MODULE tluadvTRA


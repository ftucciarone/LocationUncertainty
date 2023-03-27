!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlufluxTRA
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Transport under Location Uncertainty: TRACER advection.
!! 
!! @par           Procedure specifics      
!>                Modifies the standard NEMO module traadv.F90 to include the stochastic drifting terms:
!!                \f{align*}{
!!                {\color{gray} \mathrm{d}_{t} \theta } + &
!!                {\color{black} \nabla\cdot\left\lbrace \left[\boldsymbol{u}-\boldsymbol{u}_{s}+\boldsymbol{\mu}_{t}\right)\theta\mathrm{d}t
!!                                 +  \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}\theta \right\rbrace } - 
!!                {\color{gray} \dfrac{1}{2}\nabla\cdot\left( \boldsymbol{a}\nabla\theta \right)\mathrm{d}t } =
!!                {\color{gray} \mathrm{D}^{\theta} } +
!!                {\color{gray} \mathrm{F}^{\theta} } 
!!                       
!!                \f}

!> @par           Code specifics
!!                
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
!!
!> @par           Original History of the module
!> @snippet       this Original_History
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
   !
   ! [Original_History]
   !!==============================================================================
   !!                       ***  MODULE  tlufluxTRA  ***
   !! Ocean active tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-09  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            3.6  !  2011-06  (G. Madec)  Addition of Mixed Layer Eddy parameterisation
   !!            3.7  !  2014-05  (G. Madec)  Add 2nd/4th order cases for CEN and FCT schemes 
   !!             -   !  2014-12  (G. Madec) suppression of cross land advection option
   !!            3.6  !  2015-06  (E. Clementi) Addition of Stokes drift in case of wave coupling
   !!----------------------------------------------------------------------
   ! [Original_History]
   !
   ! [mod_dep]
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE domvvl         ! variable vertical scale factors
   USE sbcwave        ! wave module
   USE sbc_oce        ! surface boundary condition: ocean
   USE traadv_cen     ! centered scheme            (tra_adv_cen  routine)
   USE traadv_fct     ! FCT      scheme            (tra_adv_fct  routine)
   USE traadv_mus     ! MUSCL    scheme            (tra_adv_mus  routine)
   USE traadv_ubs     ! UBS      scheme            (tra_adv_ubs  routine)
   USE traadv_qck     ! QUICKEST scheme            (tra_adv_qck  routine)
   USE tramle         ! Mixed Layer Eddy transport (tra_mle_trp  routine)
   USE ldftra         ! Eddy Induced transport     (ldf_eiv_trp  routine)
   USE ldfslp         ! Lateral diffusion: slopes of neutral surfaces
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE diaptr         ! Poleward heat transport 
   !
   USE tlu
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O module
   USE prtctl         ! Print control
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE
   !
   ! [public_sub]
   !!-------------------------------------------------------------------------
   !!   tlu_tra_adv       : compute ocean tracer advection trend
   !!   tlu_tra_adv_init  : control the different options of advection scheme
   !!-------------------------------------------------------------------------
   !
   PUBLIC   tlu_tra_adv        ! called by step.F90
   PUBLIC   tlu_tra_adv_init   ! called by nemogcm.F90
   ! [public_sub]
   !                            !!* Namelist namtra_adv *
   LOGICAL ::   ln_traadv_OFF    ! no advection on T and S
   LOGICAL ::   ln_traadv_cen    ! centered scheme flag
   INTEGER ::      nn_cen_h, nn_cen_v   ! =2/4 : horizontal and vertical choices of the order of CEN scheme
   LOGICAL ::   ln_traadv_fct    ! FCT scheme flag
   INTEGER ::      nn_fct_h, nn_fct_v   ! =2/4 : horizontal and vertical choices of the order of FCT scheme
   LOGICAL ::   ln_traadv_mus    ! MUSCL scheme flag
   LOGICAL ::      ln_mus_ups           ! use upstream scheme in vicinity of river mouths
   LOGICAL ::   ln_traadv_ubs    ! UBS scheme flag
   INTEGER ::      nn_ubs_v             ! =2/4 : vertical choice of the order of UBS scheme
   LOGICAL ::   ln_traadv_qck    ! QUICKEST scheme flag

   INTEGER ::   nadv             ! choice of the type of advection scheme
   !                             ! associated indices:
   INTEGER, PARAMETER ::   np_NO_adv  = 0   ! no T-S advection
   INTEGER, PARAMETER ::   np_CEN     = 1   ! 2nd/4th order centered scheme
   INTEGER, PARAMETER ::   np_FCT     = 2   ! 2nd/4th order Flux Corrected Transport scheme
   INTEGER, PARAMETER ::   np_MUS     = 3   ! MUSCL scheme
   INTEGER, PARAMETER ::   np_UBS     = 4   ! 3rd order Upstream Biased Scheme
   INTEGER, PARAMETER ::   np_QCK     = 5   ! QUICK scheme
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: traadv.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tlu_tra_nonCons( kt, ptn, pta, uadv, vadv, wadv, r_Scale)
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
      REAL(wp), DIMENSION(jpi,jpj,jpk),      INTENT(in   ) ::   uadv, vadv, wadv !The components advecting the velocity field
      REAL(wp), DIMENSION(jpts)            , INTENT(in   ) ::   r_Scale               ! Schmidt number      
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptn
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pta             ! modified advection term for this component
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
      ALLOCATE( zwpw(jpi,jpj,jpk), &
      &         zwa1(jpi,jpj,jpk), &
      &         zwa2(jpi,jpj,jpk), &
      &         zwa3(jpi,jpj,jpk)  )
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
         &                                 zwa1(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( ptn(3:jpi  ,2:jpjm1,1:jpkm1,jn) - ptn(2:jpi-1,2:jpjm1,1:jpkm1,jn) ) &
         &                               + zwa1(1:jpim1-1,2:jpjm1  ,1:jpkm1) * ( ptn(2:jpi-1,2:jpjm1,1:jpkm1,jn) - ptn(1:jpi-2,2:jpjm1,1:jpkm1,jn) ) &
         !
         ! y-Directed gradient with interpolation in y
         !
         &                               + zwa2(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( ptn(2:jpim1,3:jpj  ,1:jpkm1,jn) - ptn(2:jpim1,2:jpj-1,1:jpkm1,jn) ) &
         &                               + zwa2(2:jpim1  ,1:jpjm1-1,1:jpkm1) * ( ptn(2:jpim1,2:jpj-1,1:jpkm1,jn) - ptn(2:jpim1,1:jpj-2,1:jpkm1,jn) ) &
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
   END SUBROUTINE tlu_tra_nonCons

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trahhdiff ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Compute the ocean tracer advection trend including the stochastic
   !!                contributions
   !!
   !> @par           Procedure specifics
   !>                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !! 
   !> @par           Code specifics               
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!
   !!
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     ptb: Tracer before field (Euler scheme)
   !> @param[inout]  pta: Tracer trend to be updated
   !> @param[in]     r_Scale: Scaling factor for diffusion
   !! 
   !! @result        Update (pta) with the advection term following nadv
   !!
   !> @par           Diagnostic
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !  @note          
   !  @todo              
   !!
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_tra_adv
   ! [tlu_tra_adv]
   SUBROUTINE tlu_tra_adv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   jk   ! dummy loop index
      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zun, zvn, zwn   ! 3D workspace


      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zuADV, zvADV, zwADV   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zunoi, zvnoi, zwnoi   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zuisd, zvisd, zwisd   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zubia, zvbia, zwbia   ! 3D workspace

      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdt, ztrds

      
      REAL(wp), DIMENSION(jpi,jpj,jpk,2)      :: noitrnd, isdtrnd, biatrnd, utrnd   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk,2)      :: tsn_T, tsn_S
      REAL(wp), DIMENSION(jpi,jpj,jpk,2)      :: noi_gradtsn
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tlu_tra_adv')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dt =         rdt   ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dt = 2._wp * rdt   ! at nit000 or nit000+1 (Leapfrog)
      ENDIF
      !
      !                                         !==  effective transport  ==!
      zun(:,:,jpk) = 0._wp
      zvn(:,:,jpk) = 0._wp
      zwn(:,:,jpk) = 0._wp
      ! 
      IF( ln_wave .AND. ln_sdw )  THEN
         DO jk = 1, jpkm1                                                       ! eulerian transport + Stokes Drift
            zun(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( un(:,:,jk) + usd(:,:,jk) )
            zvn(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( vn(:,:,jk) + vsd(:,:,jk) )
            zwn(:,:,jk) = e1e2t(:,:)                 * ( wn(:,:,jk) + wsd(:,:,jk) )
         END DO
      ELSE
         DO jk = 1, jpkm1
            zun(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( un(:,:,jk) )               ! eulerian transport only
            zvn(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( vn(:,:,jk) )
            zwn(:,:,jk) = e1e2t(:,:)                 * ( wn(:,:,jk) )
         END DO
      ENDIF
      ! 
      DO jk = 1, jpkm1
         !
         zuisd(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( uisd_n(:,:,jk) )
         zvisd(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( visd_n(:,:,jk) )
         zwisd(:,:,jk) = e1e2t(:,:)                 * ( wisd_n(:,:,jk) )
         !
         zunoi(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( unoi(:,:,jk) )
         zvnoi(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( vnoi(:,:,jk) )
         zwnoi(:,:,jk) = e1e2t(:,:)                 * ( wnoi(:,:,jk) )
         !
         IF ( ln_tlu_bia ) THEN
            ! 
            zubia(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( ubia_n(:,:,jk) )
            zvbia(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( vbia_n(:,:,jk) )
            zwbia(:,:,jk) = e1e2t(:,:)                 * ( wbia_n(:,:,jk) )
            !
         END IF 
         !
      END DO 
      !
      utrnd = 0._wp  
      noitrnd = 0._wp
      isdtrnd = 0._wp
      !
      CALL tra_adv_cen    ( kt, nit000, 'TRA',   zun,   zvn,   zwn, tsn,   utrnd, jpts, 2, 2 )
      CALL tra_adv_cen    ( kt, nit000, 'TRA', zunoi, zvnoi, zwnoi, tsn, noitrnd, jpts, 2, 2 )
      CALL tra_adv_cen    ( kt, nit000, 'TRA', zuisd, zvisd, zwisd, tsn, isdtrnd, jpts, 2, 2 )
      !
      IF ( ln_tlu_bia ) THEN
         !
         biatrnd = 0._wp
         CALL tra_adv_cen    ( kt, nit000, 'TRA', zubia, zvbia, zwbia, tsn, biatrnd, jpts, 2, 2 )
         !
      END IF
      !
      !
      !
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN                                ! add z-tilde and/or vvl corrections
         zun(:,:,:) = zun(:,:,:) + un_td(:,:,:)
         zvn(:,:,:) = zvn(:,:,:) + vn_td(:,:,:)
      ENDIF
      !
      zun(:,:,jpk) = 0._wp                                                      ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp
      zwn(:,:,jpk) = 0._wp
      !
      IF( ln_ldfeiv .AND. .NOT. ln_traldf_triad )   &
         &              CALL ldf_eiv_trp( kt, nit000, zun, zvn, zwn, 'TRA' )   ! add the eiv transport (if necessary)
      !
      IF( ln_mle    )   CALL tra_mle_trp( kt, nit000, zun, zvn, zwn, 'TRA' )   ! add the mle transport (if necessary)
      !
      CALL iom_put( "uocetr_eff", zun )                                        ! output effective transport      
      CALL iom_put( "vocetr_eff", zvn )
      CALL iom_put( "wocetr_eff", zwn )
      !
!!gm ???
      IF( ln_diaptr )   CALL dia_ptr( zvn )                                    ! diagnose the effective MSF 
!!gm ???
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk), ztrds(jpi,jpj,jpk) )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
      CALL tlu_tra_nonCons( kt, tsn, tsa, unoi, vnoi, wnoi, (/ Sc_T, Sc_S /))
      !
      ! Diffusion as double advection
      !
      noi_gradtsn = 0._wp
      CALL tlu_tra_nonCons( kt,                  tsn, noi_gradtsn, unoi, vnoi, wnoi, (/ Sc_T, Sc_S /))
      CALL tlu_tra_nonCons( kt, 0.5_wp * noi_gradtsn,         tsa, unoi, vnoi, wnoi, (/ Sc_T, Sc_S /))
      ! 
      ! Duplicate tsn and tsb to scale them down for advection purposes
      tsn_T = tsn
      tsn_S = tsn
      tsn_T(:,:,:,jp_sal) = 0._wp
      tsn_S(:,:,:,jp_tem) = 0._wp
      !
      !
      SELECT CASE ( nadv )                      !==  compute advection trend and add it to general trend  ==!
      !
      CASE ( np_CEN )                                 ! Centered scheme : 2nd / 4th order
         !
         ! Scale Temperature
         !
         zuADV = zun - zuisd / Sc_T**2
         zvADV = zvn - zvisd / Sc_T**2
         zwADV = zwn - zwisd / Sc_T**2
         !
         IF ( ln_tlu_bia ) THEN
            !
            zuADV = zuADV + biaSIGN * zubia / Sc_T
            zvADV = zvADV + biaSIGN * zvbia / Sc_T
            zwADV = zwADV + biaSIGN * zwbia / Sc_T
            !
         END IF
	 !
         CALL tra_adv_cen    ( kt, nit000, 'TRA',       zuADV, zvADV, zwADV     , tsn_T, tsa, jpts, nn_cen_h, nn_cen_v )
         !
         !  Scale Salinity
         !
         zuADV = zun - zuisd / Sc_S**2
         zvADV = zvn - zvisd / Sc_S**2
         zwADV = zwn - zwisd / Sc_S**2
         !
         IF ( ln_tlu_bia ) THEN
            !
            zuADV = zuADV + biaSIGN * zubia / Sc_S
            zvADV = zvADV + biaSIGN * zvbia / Sc_S
            zwADV = zwADV + biaSIGN * zwbia / Sc_S
            !
         END IF
         !
         CALL tra_adv_cen    ( kt, nit000, 'TRA',       zuADV, zvADV, zwADV     , tsn_S, tsa, jpts, nn_cen_h, nn_cen_v )


      CASE ( np_FCT )                                 ! FCT scheme      : 2nd / 4th order
         CALL tra_adv_fct    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts, nn_fct_h, nn_fct_v )
      CASE ( np_MUS )                                 ! MUSCL
         CALL tra_adv_mus    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb,      tsa, jpts        , ln_mus_ups ) 
      CASE ( np_UBS )                                 ! UBS
         CALL tra_adv_ubs    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts        , nn_ubs_v   )
      CASE ( np_QCK )                                 ! QUICKEST
         CALL tra_adv_qck    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts                     )
      !
      END SELECT
      !
      !      
      CALL iom_put(   "tlu_uadv_T", utrnd(:,:,:,1) )
      CALL iom_put(   "tlu_uadv_S", utrnd(:,:,:,2) )

      CALL iom_put( "tlu_noiadv_T", noitrnd(:,:,:,1) )
      CALL iom_put( "tlu_noiadv_S", noitrnd(:,:,:,2) )

      CALL iom_put( "tlu_isdadv_T", isdtrnd(:,:,:,1) )
      CALL iom_put( "tlu_isdadv_S", isdtrnd(:,:,:,2) )
      !
      IF ( ln_tlu_bia ) THEN
         !
         CALL iom_put( "tlu_biaadv_T", biatrnd(:,:,:,1) )
         CALL iom_put( "tlu_biaadv_S", biatrnd(:,:,:,2) )
         ! 
      END IF
      !
      IF( l_trdtra )   THEN                      ! save the advective trends for further diagnostics
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = tsa(:,:,jk,jp_tem) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = tsa(:,:,jk,jp_sal) - ztrds(:,:,jk)
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_totad, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_totad, ztrds )
         DEALLOCATE( ztrdt, ztrds )
      ENDIF
      !                                              ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop( 'tlu_tra_adv' )
      !
   END SUBROUTINE tlu_tra_adv
   ! [tlu_tra_adv]

   SUBROUTINE tlu_tra_adv_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_init  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options for 
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! Local integers
      !
      NAMELIST/namtra_adv/ ln_traadv_OFF,                        &   ! No advection
         &                 ln_traadv_cen , nn_cen_h, nn_cen_v,   &   ! CEN
         &                 ln_traadv_fct , nn_fct_h, nn_fct_v,   &   ! FCT
         &                 ln_traadv_mus , ln_mus_ups,           &   ! MUSCL
         &                 ln_traadv_ubs ,           nn_ubs_v,   &   ! UBS
         &                 ln_traadv_qck                             ! QCK
      !!----------------------------------------------------------------------
      !
      !                                !==  Namelist  ==!
      REWIND( numnam_ref )                   ! Namelist namtra_adv in reference namelist : Tracer advection scheme
      READ  ( numnam_ref, namtra_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_adv in reference namelist', lwp )
      !
      REWIND( numnam_cfg )                   ! Namelist namtra_adv in configuration namelist : Tracer advection scheme
      READ  ( numnam_cfg, namtra_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_adv in configuration namelist', lwp )
      IF(lwm) WRITE( numond, namtra_adv )
      !
      IF(lwp) THEN                           ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_tra_adv_init : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      No advection on T & S                     ln_traadv_OFF = ', ln_traadv_OFF
         WRITE(numout,*) '      centered scheme                           ln_traadv_cen = ', ln_traadv_cen
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_cen_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_cen_v   = ', nn_fct_v
         WRITE(numout,*) '      Flux Corrected Transport scheme           ln_traadv_fct = ', ln_traadv_fct
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_fct_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_fct_v   = ', nn_fct_v
         WRITE(numout,*) '      MUSCL scheme                              ln_traadv_mus = ', ln_traadv_mus
         WRITE(numout,*) '            + upstream scheme near river mouths    ln_mus_ups = ', ln_mus_ups
         WRITE(numout,*) '      UBS scheme                                ln_traadv_ubs = ', ln_traadv_ubs
         WRITE(numout,*) '            vertical   2nd/4th order               nn_ubs_v   = ', nn_ubs_v
         WRITE(numout,*) '      QUICKEST scheme                           ln_traadv_qck = ', ln_traadv_qck
      ENDIF
      !
      !                                !==  Parameter control & set nadv ==!
      ioptio = 0                       
      IF( ln_traadv_OFF ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_NO_adv   ;   ENDIF
      IF( ln_traadv_cen ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_CEN      ;   ENDIF
      IF( ln_traadv_fct ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_FCT      ;   ENDIF
      IF( ln_traadv_mus ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_MUS      ;   ENDIF
      IF( ln_traadv_ubs ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_UBS      ;   ENDIF
      IF( ln_traadv_qck ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_QCK      ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'tra_adv_init: Choose ONE advection option in namelist namtra_adv' )
      !
      IF( ln_traadv_cen .AND. ( nn_cen_h /= 2 .AND. nn_cen_h /= 4 )   &          ! Centered
                        .AND. ( nn_cen_v /= 2 .AND. nn_cen_v /= 4 )   ) THEN
        CALL ctl_stop( 'tra_adv_init: CEN scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_fct .AND. ( nn_fct_h /= 2 .AND. nn_fct_h /= 4 )   &          ! FCT
                        .AND. ( nn_fct_v /= 2 .AND. nn_fct_v /= 4 )   ) THEN
        CALL ctl_stop( 'tra_adv_init: FCT scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_ubs .AND. ( nn_ubs_v /= 2 .AND. nn_ubs_v /= 4 )   ) THEN     ! UBS
        CALL ctl_stop( 'tra_adv_init: UBS scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_ubs .AND. nn_ubs_v == 4 ) THEN
         CALL ctl_warn( 'tra_adv_init: UBS scheme, only 2nd FCT scheme available on the vertical. It will be used' )
      ENDIF
      IF( ln_isfcav ) THEN                                                       ! ice-shelf cavities
         IF(  ln_traadv_cen .AND. nn_cen_v == 4    .OR.   &                            ! NO 4th order with ISF
            & ln_traadv_fct .AND. nn_fct_v == 4   )   CALL ctl_stop( 'tra_adv_init: 4th order COMPACT scheme not allowed with ISF' )
      ENDIF
      !
      !                                !==  Print the choice  ==!  
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE ( nadv )
         CASE( np_NO_adv  )   ;   WRITE(numout,*) '   ==>>>   NO T-S advection'
         CASE( np_CEN     )   ;   WRITE(numout,*) '   ==>>>   CEN      scheme is used. Horizontal order: ', nn_cen_h,   &
            &                                                                        ' Vertical   order: ', nn_cen_v
         CASE( np_FCT     )   ;   WRITE(numout,*) '   ==>>>   FCT      scheme is used. Horizontal order: ', nn_fct_h,   &
            &                                                                        ' Vertical   order: ', nn_fct_v
         CASE( np_MUS     )   ;   WRITE(numout,*) '   ==>>>   MUSCL    scheme is used'
         CASE( np_UBS     )   ;   WRITE(numout,*) '   ==>>>   UBS      scheme is used'
         CASE( np_QCK     )   ;   WRITE(numout,*) '   ==>>>   QUICKEST scheme is used'
         END SELECT
      ENDIF
      !
      CALL tra_mle_init            !== initialisation of the Mixed Layer Eddy parametrisation (MLE)  ==!
      !
   END SUBROUTINE tlu_tra_adv_init

  !!======================================================================
END MODULE tluadvTRA



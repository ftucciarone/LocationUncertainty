!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlustp
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
!> @brief         Transport under Location Uncertainty: CORIOLIS effect and STOCHASTIC PRESSURE correction
!! 
!! @par           Procedure specifics      
!>                Compute all the terms of the Coriolis effect and the Stochastic Pressure corrections. Designed to handle all the possible
!!                models for the stochastic pressure through a SELECT CASE statement.
!!                \f{align*}{
!!                {\color{gray}  \mathrm{d}_{t} \boldsymbol{u} } + &
!!                {\color{gray}  \nabla\cdot\left\lbrace \left[\boldsymbol{u}-\boldsymbol{u}_{s}+\boldsymbol{\mu}_{t}\right)\boldsymbol{u}\mathrm{d}t
!!                                 +  \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}\boldsymbol{u} \right\rbrace } - 
!!                {\color{gray}  \dfrac{1}{2}\nabla\cdot\left( \boldsymbol{a}\nabla\boldsymbol{u} \right)\mathrm{d}t } \\ &+
!!                {\color{black} \boldsymbol{f}\times } \left( 
!!                {\color{gray}  \boldsymbol{u}\mathrm{d}t + }
!!                {\color{black} \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }
!!                              \right) = 
!!                {\color{gray}  - \nabla p\mathrm{d}t }
!!                {\color{black} - \nabla\mathrm{d}p_{t}^{\sigma} } 
!!                {\color{gray}  + \mathrm{D}^{\boldsymbol{u}} }
!!                {\color{gray}  + \mathrm{F}^{\boldsymbol{u}} } 
!!                       
!!                \f}
!> @par           Code specifics
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
!> @todo          Write error messages to appropriate unit
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tluSTP
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE tlu              ! Initialization of stochastic structures
   USE dynvor           ! Computation of vorticity (see warnings)
   ! [mod_dep]
   !
   IMPLICIT NONE   
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_stpdyn    ! Called by step.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_stpdyn ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Handles the different choices for the Coriolis-Pressure model.
   !!
   !> @par           Procedure specifics
   !>                Handles the different choices for the Coriolis-Pressure model.
   !!
   !!
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     stpmod: 2-characters wide argument for model selection 
   !!
   !! 
   !> @par           Code specifics               
   !!                The input parameter \p stpmod governs the SELECT CASE statement. The possible choices are: <br>
   !!                <br>
   !!                <b> Geostrophic balance </b> <br>
   !!                The horizontal gradient of pressure and the vertically averaged coriols term are assumed to be 
   !!                in approximate balance, i.e
   !!                \f{align*}{
   !!                \boldsymbol{f}\times 
   !!                      \overline{ \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }^{\, z} \approx 
   !!                - \nabla\mathrm{d}p_{t}^{\sigma}, 
   !!                \f}
   !!                thus the routine \ref tlu_cor_noi implements the coriolis correction for a purely baroclinic noise, that is:
   !!                \f{align*}{
   !!                \mathrm{RHS}^{\boldsymbol{u}} =  \mathrm{RHS}^{\boldsymbol{u}}  +
   !!                \boldsymbol{f}\times  \left(
   !!                       \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} - 
   !!                       \overline{ \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }^{\, z} 
   !!                \right)  
   !!                \f}
   !>                @snippet this Coriolis_Geostrophic
   !!                <b> Quasi-Geostrophic Regime </b>  
   !!                There is no balance between the terms, thus the coriolis term are implemented completely as 
   !!                \f{align*}{
   !!                & \mathrm{RHS}^{\boldsymbol{u}} =  \mathrm{RHS}^{\boldsymbol{u}}  +
   !!                \boldsymbol{f}\times  \left(
   !!                       \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} - 
   !!                       \overline{ \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }^{\, z} 
   !!                \right),  \\
   !!                & \mathrm{RHS}^{\boldsymbol{u}} =  \mathrm{RHS}^{\boldsymbol{u}}  +
   !!                \boldsymbol{f}\times  \overline{ \dfrac{1}{2}\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} }^{\, z},
   !!                \f}
   !!                the Stochastic pressure is defined as
   !!                \f{align*}{
   !!                \mathrm{d}p_{t}^{\sigma}\left(x, y, z\right) = \int_{\eta_{b}}^{z} \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} 
   !!                                           \cdot \nabla w \, \mathrm{d}\zeta
   !!                \f}
   !!                and implemented in \ref stoch_pressure_correction as
   !!                \f{align*}{
   !!                \mathrm{RHS}^{\boldsymbol{u}} =  \mathrm{RHS}^{\boldsymbol{u}}  +
   !!                \nabla \left( 
   !!                \int_{\eta_{b}}^{z} \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t} 
   !!                                           \cdot \nabla w \, \mathrm{d}\zeta
   !!                \right)
   !!                \f}
   !!                The usual pressure component is modified as:
   !!                \f{align*}{
   !!                p^{\prime}\left(x, y, z\right) =\left.\displaystyle\frac{\nu}{3}\nabla\boldsymbol{\cdot} \mathbf{u}_{s}
   !!                                                \right\vert_{\eta_{b}}^{z} 
   !!                + \int_{\eta_{b}}^{z} b + \mathbf{u}_{s}\boldsymbol{\cdot}\nabla w
   !!                + \frac{1}{2}\nabla\boldsymbol{\cdot}\left( \boldsymbol{a}\nabla w \right)\,\mathrm{d}\zeta,
   !!                \f}
   !!                and implemented in \ref pressure_correction as
   !!                \f{align*}{
   !!                 \mathrm{RHS}^{\boldsymbol{u}} =  \mathrm{RHS}^{\boldsymbol{u}}  +
   !!                \nabla \left( 
   !!                       \left.\displaystyle\frac{\nu}{3}\nabla\boldsymbol{\cdot} \mathbf{u}_{s}
   !!                                                \right\vert_{\eta_{b}}^{z} 
   !!                + \int_{\eta_{b}}^{z} \mathbf{u}_{s}\boldsymbol{\cdot}\nabla w
   !!                + \frac{1}{2}\nabla\boldsymbol{\cdot}\left( \boldsymbol{a}\nabla w \right)\,\mathrm{d}\zeta,
   !!                \right)
   !!                \f}
   !>                @snippet this Coriolis_QuasiGeostrophic
   !!
   !! 
   !! @result        Calls different routines depending on the model choice, these routines update ua, va
   !!
   !  @par           Diagnostic
   !  @note          
   !> @todo          Implement The Full Stochastic Model    
   !!
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_stpdyn
   ! [tlu_stpdyn]
   SUBROUTINE tlu_stpdyn(kt, stpmod)
      INTEGER,                                  INTENT(in   ) :: kt                     ! ocean time-step index
      CHARACTER(len=2),                         INTENT(in   ) :: stpmod                 ! = G (Model indicator)   
      REAL(wp),         DIMENSION(jpi,jpj,jpk)                :: pcorr, pvn    ! now velocities
      REAL(wp),         DIMENSION(jpi,jpj,jpk)                :: stpgu, stpgv  ! now velocities
      REAL(wp),         DIMENSION(jpi,jpj,jpk)                :: utrnd, vtrnd  ! now velocities  
      ! 
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_stpdyn')   ! [NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_stpdyn : Stochastic Pressure model '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      SELECT CASE( stpmod )     

         ! [Coriolis_Geostrophic]
         CASE ( 'Ge' ) !For Geostrophic balance between noise and stp
            !
            ! Geostrophic balance can only be between two-dimensional quantities,
            ! hence we use the barotropic part of the noise
            CALL tlu_cor_noi( kt,             unoi,             vnoi, ua, va )
            CALL tlu_cor_noi( kt, biaSign * ubia_n, biaSign * vbia_n, ua, va )
            !
         ! [Coriolis_Geostrophic]
         !
         ! [Coriolis_QuasiGeostrophic]
         CASE ( 'QG' ) !For Quasi-Geostrophic noise 
            !
            utrnd = ua
            vtrnd = va
            !
            !
            ! Coriolis contributions
            CALL tlu_cor_noi( kt, unoi, vnoi, ua, va )
            !
            ! Print noise fields
            !
            CALL iom_put( "f_x_sigmav", ua-utrnd )
            CALL iom_put( "f_x_sigmau", va-vtrnd )
            !
            IF (ln_tlu_bia) CALL tlu_cor_noi( kt, biaSign * ubia_n, biaSign * vbia_n, ua, va )
            !
            ! Stochastic pressure model:
            !
            pvn = 0._wp
            pcorr = 0._wp
            !            Diffusion
            CALL tlu_stp_adv( kt,    wn, pcorr, unoi, vnoi, wnoi)
            CALL tlu_stp_adv( kt, pcorr, pcorr, unoi, vnoi, wnoi)
            pcorr = pcorr * 0.5_wp
            !
            !     Ito-Stokes drift
            CALL tlu_stp_adv( kt,    wn, pcorr, uisd_n, visd_n, wisd_n)
            !
            !                Noise 
            CALL tlu_stp_adv( kt,    wn, pcorr, unoi, vnoi, wnoi)
            !
            !  Vertical cumulation
            CALL tlu_vert_int( kt, pcorr, pvn)
            !
            !      Compressibility is missing: \nabla\cdot\boldsymbol{u}_s
            ! CALL tlu_compr( kt, press_corr)
            !
            !   Gradient of all the terms
            !
            CALL tlu_stp_grad( kt, pvn, stpgu, stpgv)
            !
            ! Print noise fields
            !
            CALL iom_put( "p_correct", pvn )
            CALL iom_put( "dpt_dx", stpgu )
            CALL iom_put( "dpt_dy", stpgv )
            !
            ua = ua - stpgu
            va = va - stpgv
            !
         ! [Coriolis_QuasiGeostrophic]
         !
         CASE DEFAULT 
            CALL ctl_stop('STOP','tlu_stpdyn: wrong value for stpmod'  )
      END SELECT
      !
      IF( ln_timing ) CALL timing_stop('tlu_stpdyn')   ! [NEMO] check
      !
   END SUBROUTINE tlu_stpdyn
   ! [tlu_stpdyn]


   !!--------------------------------------------------------------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_cor_noi ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Compute the Coriolis term due to a given noise 
   !!                @f$ \boldsymbol{f}\times\boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}@f$
   !!                and add it to the general trend of the momentum equation.
   !!
   !> @details       The coriolis term under Location uncertainty is defined as:
   !!                Computes the horizontal terms of the stochastic diffusion, using a finite volumes approach
   !!                for the diffusive fluxes @f$ F_{T,x} @f$ and @f$ F_{T,4} @f$ centred on a @f$ T- @f$ point 
   !!                \f{align*}{
   !!                   \mathrm{RHS}^{u} =  \mathrm{RHS}^{u}
   !!                      + f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y}
   !!                \f}
   !!                \f{align*}{
   !!                   \mathrm{RHS}^{v} =  \mathrm{RHS}^{v} 
   !!                      -f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}
   !!                \f}
   !> @par           Code specifics
   !!                The coriolis terms are computed following the standard Energy conserving scheme in NEMO 
   !!                \f{align*}{
   !!                   f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x} = 
   !!                    \frac{1}{2}\dfrac{1}{e_{1u}}\overline{\dfrac{f + \mathrm{metric}}{e_{3f}} \,
   !!                    \overline{e_{1v}e_{3v}v}^{\, i+1/2}
   !!                    }^{\, j},
   !!                   \qquad
   !!                   f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y} = 
   !!                    \frac{1}{2}\dfrac{1}{e_{2v}}\overline{\dfrac{f + \mathrm{metric}}{e_{3f}} \,
   !!                    \overline{e_{2u}e_{3u}u}^{\, j+1/2}
   !!                    }^{\, i},
   !!                \f}
   !!                where the metric terms represent corrections due to the spherical polar coordinates in 
   !!                the flux form.
   !! 
   !! 
   !> @param[in]     kt: ocean time-step integer indexi
   !> @param[in]     punoi, pvnoi: @f$x-@f$wise and @f$y-@f$wise noise
   !> @param[inout]  pua, pva: velocity trends updated
   !! 
   !! @result        Update (ua,va) with the now noise based coriolisa term trend
   !!
   !  @note   
   !  @todo              
   !!
   !!--------------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_cor_noi 
   ! [tlu_cor_noi]
   SUBROUTINE tlu_cor_noi( kt, uadv, vadv, pua, pva )
      INTEGER                         , INTENT(in   ) ::   kt                  ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk)                ::   uadv, vadv        ! Noise from now realisation
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva           ! total v-trend
      !
      INTEGER                                         ::   ji, jj, jk          ! dummy loop indices
      REAL(wp)                                        ::   zx1, zy1, zx2, zy2  ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)                    ::   zwx, zwy, zwz       ! 2D workspace
      REAL(wp)                                        ::   r1_2  = 0.50_wp    ! =1/2
      REAL(wp)                                        ::   r1_4  = 0.25_wp    ! =1/4
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_cor_noi : Noise based Coriolis term : energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !   
      zwz(:,:) = ff_f(:,:) 
      !      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         IF( ln_dynvor_msk ) THEN          !==  mask/unmask vorticity ==!
            DO jj = 1, jpj
               DO ji = 1, jpi   ! vector opt.
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END DO
            END DO
         ENDIF

         IF( ln_sco ) THEN      ! Generalized s-coordinate
            zwz(:,:) = zwz(:,:) / e3f_n(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u_n(:,:,jk) * uadv(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v_n(:,:,jk) * vadv(:,:,jk)
         ELSE
            zwx(:,:) = e2u(:,:) * uadv(:,:,jk) * umask(:,:,jk)
            zwy(:,:) = e1v(:,:) * vadv(:,:,jk) * vmask(:,:,jk)
         ENDIF
         !                                   !==  compute and add the Noise based Coriolis term trend  =!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zy1 = zwy(ji  ,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji  ,jj  ) + zwy(ji+1,jj  )
               zx1 = zwx(ji-1,jj  ) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj  ) + zwx(ji  ,jj+1)
               pua(ji,jj,jk) = pua(ji,jj,jk) + r1_4 * r1_e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               pva(ji,jj,jk) = pva(ji,jj,jk) - r1_4 * r1_e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE tlu_cor_noi
   ! [tlu_cor_noi]

   SUBROUTINE tlu_stp_adv( kt, pcn, pca, uadv, vadv, wadv)
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
      INTEGER                         , INTENT(in   ) ::   kt               ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   uadv, vadv, wadv ! The components advecting the velocity field
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pca         ! modified advection term for this component
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwe1, zwe2         ! workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwa1, zwa2, zwa3   ! workspace
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_stp_adv')   ! [NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'tlu_stp_adv : advection by the noise of vertical velocity'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( zwa1(jpi,jpj,jpk) , zwa2(jpi,jpj,jpk) , zwa3(jpi,jpj,jpk))
      ALLOCATE( zwe1(jpi,jpj,jpk) , zwe2(jpi,jpj,jpk) )
      !
      zwa1 = 0._wp
      zwa2 = 0._wp
      zwa3 = 0._wp
      !
      zwe1 = spread(r1_e1t,3,jpkm1)
      zwe2 = spread(r1_e2f,3,jpkm1) 
      !
      zwa1( 1:jpim1, 2:jpjm1, 2:jpkm1) = ( uadv(1:jpim1, 2:jpjm1, 1:jpk-1 )   &
      &                                  + uadv(1:jpim1, 2:jpjm1, 2:jpk   ) ) &
      &                                * (  pcn(2:jpi  , 2:jpjm1, 1:jpkm1 )   &
      &                                  -  pcn(1:jpi-1, 2:jpjm1, 1:jpkm1 ) ) &
      &                                *   zwe1(1:jpim1, 2:jpjm1, 1:jpkm1 ) * 0.5_wp
      !
      zwa2( 2:jpim1, 1:jpjm1, 2:jpkm1) = ( vadv(2:jpim1, 1:jpjm1, 1:jpk-1 )   &
      &                                  + vadv(2:jpim1, 1:jpjm1, 2:jpk   ) ) &
      &                                * (  pcn(2:jpim1, 2:jpj  , 1:jpkm1 )   &
      &                                  -  pcn(2:jpim1, 1:jpj-1, 1:jpkm1 ) ) &
      &                                *   zwe2(2:jpim1, 1:jpjm1, 1:jpkm1 ) * 0.5_wp
      !
      zwa3( 2:jpim1, 2:jpjm1, 2:jpkm1) = ( wadv(2:jpi-1, 1:jpjm1, 2:jpkm1 )   &
      &                                  + wadv(3:jpi  , 1:jpjm1, 2:jpkm1 ) ) &
      &                                * (  pcn(2:jpim1, 2:jpjm1, 1:jpkm1 )   &
      &                                  -  pcn(2:jpim1, 2:jpjm1, 2:jpk   ) ) &
      &                                / e3uw_n(2:jpim1, 1:jpjm1, 2:jpkm1 )
      !
      zwa3( :, :,   1) = 0._wp
      zwa3( :, :, jpk) = 0._wp
      !
      zwa3( :, :, : ) = zwa3( :, :, :) * 0.5_wp
      !
      pca( 2:jpim1, 2:jpjm1, 1:jpkm1 ) =    pca(2:jpim1, 2:jpjm1, 1:jpkm1 )   &
      !
      &                                - ( zwa1(2:jpi-1, 2:jpjm1, 1:jpkm1 )   &
      &                                  + zwa1(1:jpi-2, 2:jpjm1, 1:jpkm1 )   &
      !
      &                                  + zwa2(2:jpim1, 2:jpj-1, 1:jpkm1 )   &
      &                                  + zwa2(2:jpim1, 1:jpj-2, 1:jpkm1 )   &
      !
      &                                  + zwa3(2:jpim1, 2:jpjm1, 1:jpkm1 )   &
      &                                  + zwa3(2:jpim1, 2:jpjm1, 2:jpk   ) ) * 0.5_wp
      !

      !
      DEALLOCATE(zwa1, zwa2, zwa3)
      DEALLOCATE(zwe1, zwe2)
      !
      IF( ln_timing ) CALL timing_stop('tlu_stp_adv')   ! [NEMO] check
      !
   END SUBROUTINE tlu_stp_adv


   SUBROUTINE tlu_stp_grad( kt, pcn, pcu, pcv)
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
      INTEGER                         , INTENT(in   ) ::   kt               ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pcu, pcv         ! modified advection term for this component
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwe1, zwe2         ! workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwa1, zwa2, zwa3   ! workspace
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_stp_grad')   ! [NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'tlu_stp_grad : gradient in x and y of (global) stochastic pressure'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( zwe1(jpi,jpj,jpk) , zwe2(jpi,jpj,jpk) )
      !
      pcu = 0._wp
      pcv = 0._wp
      !
      zwe1 = spread(r1_e1u,3,jpkm1)
      zwe2 = spread(r1_e2v,3,jpkm1) 
      !
      pcu( 2:jpim1, 2:jpjm1, 1:jpkm1) = (  pcn(3:jpi  , 2:jpjm1, 1:jpkm1 )   &
      &                                 -  pcn(2:jpi-1, 2:jpjm1, 1:jpkm1 ) ) &
      &                               *   zwe1(1:jpim1, 2:jpjm1, 1:jpkm1 ) 
      !
      pcv( 2:jpim1, 2:jpjm1, 1:jpkm1) = (  pcn(2:jpim1, 3:jpj  , 1:jpkm1 )   &
      &                                 -  pcn(2:jpim1, 2:jpj-1, 1:jpkm1 ) ) &
      &                               *   zwe2(2:jpim1, 1:jpjm1, 1:jpkm1 )
      !
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_stp_grad', pcu , 'U', -1., pcv , 'V', -1.  )
      !
      DEALLOCATE(zwe1, zwe2)
      !
      IF( ln_timing ) CALL timing_stop('tlu_stp_grad')   ! [NEMO] check
      !
   END SUBROUTINE tlu_stp_grad


   SUBROUTINE tlu_vert_int( kt, field_in, field_out)
      !
      INTEGER                         , INTENT(in   ) :: kt           ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) :: field_in     ! modified advection term for this component
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) :: field_out    ! modified advection term for this component
      !
      INTEGER                                         :: jk           ! dummy loop indices
      !
      !!----------------------------------------------------------------------
      !
      field_out = 0._wp
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         field_out( 2:jpi-1, 2:jpj-1, jk) = field_out( 2:jpi-1, 2:jpj-1, jk + 1) + &
         &                                   field_in( 2:jpi-1, 2:jpj-1, jk    ) * e3t_n(2:jpi-1, 2:jpj-1,jk) 
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_vert_int', field_out , 'T', -1.  )
      !
   END SUBROUTINE tlu_vert_int 
   ! [tlu_wzvcmp]

END MODULE tlustp


!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlustp
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
!> @todo          read namelist from _cfg, not only _ref
!! @todo          write error messages to appropriate unit
!!
!> @warning       Sketchy use of `dynvor`: why is it needed?
!! @warning       BIGGEST WARN OF THEM ALL: THIS CODE IS NOT TESTED  
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tlustp
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
!   USE tlusvd           ! Initialization of the noise field
   USE dynvor           ! Computation of vorticity (see warnings)
   ! [mod_dep]
   !
   IMPLICIT NONE   
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_stpdyn    ! Called by tlu_bcdyn.f90
   ! [public_sub]
   !
   ! #  include "domzgr_substitute.h90"
   !!--------------------------------------------------------------------------------------------------------------------------------
CONTAINS

   !!--------------------------------------------------------------------------------------------------------------------------- 
   !!            ***  ROUTINE tlu_stpdyn ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         Computation of location uncertainty terms in the BAROCLINIC time step 
   !!
   !> @par           Code specifics 
   !> @details       The code makes a sequence of calls to the single routines that implement the dynamic parts
   !! @snippet       this Advection  
   !!                @ref tlu_hadv_noi for the advection of horizontal components, 
   !!                @ref tlu_zadv_noi for the advection of vertical components,
   !! @snippet       this Stochastic drift
   !!                @ref tlu_hhdrift for the drift of horizontal components, 
   !!                @ref tlu_zzdrift for the drift of vertical components,
   !! @snippet       this Stochastic diffusion
   !!                @ref tlu_hhdiff for the diffusion of horizontal components, 
   !!                @ref tlu_zzdiff for the diffusion of vertical components,
   !! @snippet       this Stochastic surface boundary conditions
   !!                @ref tlu_zadv_sbc for the advection of horizontal components with wind influence,
   !!                @ref    tlu_zzsbc for the drift-diffusion of vertical components with wind influence, 
   !! 
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !! 
   !! @result        Updated (ua,va) with all the Location Uncertainty terms
   !!  
   !!
   !!---------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_bcdyn 
   ! [tlu_bcdyn]
   SUBROUTINE tlu_stpdyn(kt, stpmod)
      INTEGER,          INTENT(in   ) ::   kt                     ! ocean time-step index
      CHARACTER(len=2), INTENT(in   ) ::   stpmod                 ! = G (Model indicator)    
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_stpdyn')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_stpdyn : Stochastic Pressure model '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      SELECT CASE( stpmod )     

         !> @par            Geostrophic balance between noise and stp
         CASE ( 'Ge' ) !For Geostrophic balance between noise and stp
            !
            ! Do nothing
            !
         !> @par            Quasi-Geostrophic noise
         CASE ( 'QG' ) !For Quasi-Geostrophic noise 
            !
            IF( kt == nit000 )  THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'tlu_stpdyn : QG-noise model not implemented yet, proceeding with Geostrophic'
               IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
            ENDIF
            !
            !
            ! Coriolis contributions
            ! CALL tlu_cor_noi( kt, unoi, vnoi, ua, va )
            !
            ! Stochastic pressure model
            ! CALL tba
            !
         CASE DEFAULT 
            CALL ctl_stop('STOP','tlu_stpdyn: wrong value for stpmod'  )
      END SELECT
      !
!      IF( ln_timing ) CALL timing_stop('tlu_stpdyn')   ! [NEMO] check
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
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Compute the Coriolis term due to @f$\boldsymbol{\sigma}d\boldsymbol{B}_{t}@f$
   !!                and add it to the general trend of the momentum equation.
   !!
   !> @details       The coriolis term under Location uncertainty is defined as: 
   !!                \f{align*}{
   !!                   \mathcal{C}_{\boldsymbol{u}} = 
   !!                   -\boldsymbol{f}\times\dfrac{1}{2}\boldsymbol{\sigma}dB_{t} = \left(
   !!                   \begin{array}{c}
   !!                      +f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y} \\
   !!                      -f(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}
   !!        	        \end{array}
   !!                   \right) = \left(
   !!                   \begin{array}{c}
   !!                      +\frac{1}{e_{1u}}[\frac{f}{e_{3f}}e_{1v}e_{3v}(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y}] \\
   !!                      -\frac{1}{e_{2v}}[\frac{f}{e_{3f}}e_{1u}e_{3u}(\frac{1}{2}\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}]
   !!                   \end{array}
   !!                   \right) 
   !!                \f}
   !!                The computation of this term follows the following procedure.    
   !!                \f{algorithm}{
   !!                   \caption{Computation of Coriolis term}
   !!                   \begin{algorithmic}
   !!                   \If{$s-$coordinate is employed (\texttt{ln\_sco = .TRUE.)}}
   !!                      \State Weight each layer by the corresponding scale factor $e_{3}$
   !!                         \begin{tabular}{l}
   !!                            \texttt{zwz} = $e_{3f}f$\\
   !!                            \texttt{zwx} = $e_{2u}e^{n}_{3u}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}$\\
   !!                            \texttt{zwy} = $e_{1v}e^{n}_{3v}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y}$
   !!                         \end{tabular}
   !!                   \Else
   !!                      \State  
   !!                      \begin{tabular}{l}
   !!                         $\texttt{zwx} = e_{2u}(\boldsymbol{\sigma}dB_{t})_{x}$\\
   !!                         $\texttt{zwy} = e_{1v}(\boldsymbol{\sigma}dB_{t})_{y}$
   !!                      \end{tabular}
   !!                   \EndIf
   !!                   \State \texttt{zx1} = $ 2\overline{e_{1v}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y}}^{i+1/2}$ at $j-1$
   !!                   \State \texttt{zx2} = $ 2\overline{e_{1v}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{y}}^{i+1/2}$ at $j$
   !!                   \State \texttt{zy1} = $ 2\overline{e_{2u}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}}^{i+1/2}$ at $j-1$
   !!                   \State \texttt{zy2} = $ 2\overline{e_{2u}(\boldsymbol{\sigma}d\boldsymbol{B}_{t})_{x}}^{i+1/2}$ at $j$
   !!                   \State \texttt{pua} = \texttt{pua} + $ \frac{1}{4}\frac{1}{e_{1u}}2\overline{\texttt{zwz*zy1}}^{j} $
   !!                   \State \texttt{pva} = \texttt{pva} - $\frac{1}{4}\frac{1}{e_{2v}}2\overline{\texttt{zwz*zy2}}^{i}$
   !!                   \end{algorithmic}
   !!                \f}
   !!                Notice that in case of a generalized @f$s-@f$coordinate the addition of the vertical scale factor @f$e_{3f}^{n}@f$
   !!                is required.
   !!
   !! 
   !! 
   !> @param[in]     kt: ocean time-step integer indexi
   !> @param[in]     punoi, pvnoi: @f$x-@f$wise and @f$y-@f$wise noise
   !> @param[inout]  pua, pva: velocity trends updated
   !! 
   !! @result        Update (ua,va) with the now noise based coriolisa term trend
   !!
   !! @note   
   !! @todo              
   !!
   !!--------------------------------------------------------------------------------------------------------------------------------
   !> @snippet this tlu_cor_noi 
   ! [tlu_cor_noi]
   SUBROUTINE tlu_cor_noi( kt, punoi, pvnoi, pua, pva )
      INTEGER                         , INTENT(in   ) ::   kt                  ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   punoi, pvnoi        ! Noise from now realisation
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva            ! total v-trend
      !
      INTEGER                                         ::   ji, jj, jk          ! dummy loop indices
      REAL(wp)                                        ::   zx1, zy1, zx2, zy2  ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)                    ::   zwx, zwy, zwz       ! 2D workspace
      REAL(wp)                                        ::   r1_2  = 0.50_wp    ! =1/2
      REAL(wp)                                        ::   r1_4  = 0.25_wp    ! =1/4
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_cor_noi : Noise based Coriolis term : energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         zwz(:,:) = ff_f(:,:) * r1_2         !   r1_2 comes from LU in rotating frame
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
            zwx(:,:) = e2u(:,:) * e3u_n(:,:,jk) * punoi(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v_n(:,:,jk) * pvnoi(:,:,jk)
         ELSE
            zwx(:,:) = e2u(:,:) * punoi(:,:,jk)
            zwy(:,:) = e1v(:,:) * pvnoi(:,:,jk)
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



END MODULE tlustp


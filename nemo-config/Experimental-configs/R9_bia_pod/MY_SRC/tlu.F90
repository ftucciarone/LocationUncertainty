!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu
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
!> @brief         Transport under Location Uncertainty module.
!! 
!! 
!! @par           Procedure specifics      
!> @details       Defines variables and dependencies for the stochastic parametrization of 
!!                Navier-Stokes equations. 
!!
!> @par           Code specifics
!!                This module is called by nemogcm.F90 in tlu.f90. It contains the definition of all the important variables for 
!!                Location Uncertanity implementation. These includes:
!!
!!                The baroclinic noise; 
!> @snippet       this tlu_noise
!!                \f{align*}{
!!                  \texttt{unoi} &\leftarrow \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}^{x} \\
!!                  \texttt{vnoi} &\leftarrow \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}^{y} \\
!!                  \texttt{wnoi} &\leftarrow \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}^{z} \\
!!                \f}
!!                The variance tensors, for the baroclinic and for the wind modes; 
!> @snippet       this tlu_variance
!!                \f{align*}{
!!                  \texttt{var\_ten}    &\leftarrow \boldsymbol{a}\\
!!                  \texttt{var\_tau}    &\leftarrow \boldsymbol{a}_{\tau} \\
!!                \f}
!!                The Bias, that is the Girsenov drift of the Noise;
!> @snippet       this tlu_bias
!!                \f{align*}{
!!                  \texttt{ubia\_0} &\leftarrow \boldsymbol{\mu}_{t=0}^{x} \\
!!                  \texttt{vbia\_0} &\leftarrow \boldsymbol{\mu}_{t=0}^{y} \\
!!                  \texttt{wbia\_0} &\leftarrow \boldsymbol{\mu}_{t=0}^{z} \\
!!                  \texttt{b\_ubia} &\leftarrow \overline{\boldsymbol{\mu}_{t=0}^{x}}^{z} \\
!!                  \texttt{b\_vbia} &\leftarrow \overline{\boldsymbol{\mu}_{t=0}^{y}}^{z}
!!                \f}
!!                The "now" components of the bias are useful in case of projection on isopycnals, that is
!!                \f{align*}{
!!                  \texttt{ubia\_n} &\leftarrow \boldsymbol{\mu}_{t}^{x} = \mathrm{P}(t)\boldsymbol{\mu}^{x}_{t=0} \\
!!                  \texttt{vbia\_n} &\leftarrow \boldsymbol{\mu}_{t}^{y} = \mathrm{P}(t)\boldsymbol{\mu}^{y}_{t=0}  \\
!!                  \texttt{wbia\_n} &\leftarrow \boldsymbol{\mu}_{t}^{z} = \mathrm{P}(t)\boldsymbol{\mu}^{z}_{t=0} 
!!                \f}
!!                The Ito-Stokes drift
!> @snippet       this tlu_itoStokes
!!                \f{align*}{
!!                  \texttt{uisd\_0} &\leftarrow \dfrac{1}{2}(\nabla\cdot\boldsymbol{a})^{x} \\
!!                  \texttt{visd\_0} &\leftarrow \dfrac{1}{2}(\nabla\cdot\boldsymbol{a})^{y} \\
!!                  \texttt{wisd\_0} &\leftarrow \dfrac{1}{2}(\nabla\cdot\boldsymbol{a})^{z} \\
!!                \f}
!!                The correction of the \f@ w \f@ velocity, stemming from the integration of the modified incompressibility
!!                condition
!!                \f{equation*}{
!!                  w(z) = w(0) - \int_{\eta_{b}}^{z} \nabla_{H}\cdot \boldsymbol{u} + \nabla\cdot\boldsymbol{u}_{s} \,\mathrm{d}\zeta
!!                \f}
!!                so that
!!                \f{equation*}{
!!                   \texttt{tlu\_wcorr} \leftarrow \int_{\eta_{b}}^{z} \nabla\cdot\boldsymbol{u}_{s} \,\mathrm{d}\zeta
!!                \f}
!> @snippet       this tlu_w_correction
!!
!!
!> @par           Other modules dependencies
!> @snippet       this mod_dep
!!
!> @par           Public variables
!> @snippet       this public_vars
!!
!> @par           Public routines
!> @snippet       this public_sub
!!
!! @note          Dynamics of ISD not implemented         
!!
!!------------------------------------------------------------------------------
MODULE tlu
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE iom
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_init      ! Called by nemogcm.f90
   ! [public_sub]  
   !
   ! [public_vars]
   LOGICAL,     PUBLIC                                        :: ln_tlu     = .FALSE.    !< @public Switch for Location Uncertainty (LU)
   LOGICAL,     PUBLIC                                        :: ln_pyp     = .TRUE.     !< @public Switch for projection on isopycnals
   LOGICAL,     PUBLIC                                        :: ln_tlu_nke = .FALSE.    !< @public Switch for ke rescaling 
   LOGICAL,     PUBLIC                                        :: ln_tlu_bia = .FALSE.    !< @public Switch for bias 
   LOGICAL,     PUBLIC                                        :: ln_tlu_pod = .FALSE.    !< @public Switch for POD based noise
   LOGICAL,     PUBLIC                                        :: ln_tlu_dmd = .FALSE.    !< @public Switch for DMD based noise
   LOGICAL,     PUBLIC                                        :: ln_tlu_pso = .FALSE.    !< @public Switch for Pseudo-Observation model
   LOGICAL,     PUBLIC                                        :: ln_tlu_wlt = .FALSE.    !< @public Switch for Daubechies Wavelet model
   LOGICAL,     PUBLIC                                        :: ln_tlu_gss = .FALSE.    !< @public Switch for Gaussian White noise
   !
   ! [tlu_noise]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::     unoi                !< @public Modified advection: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::     vnoi                !< @public Modified advection: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::     wnoi                !< @public Modified advection: z component
   ! [tlu_noise]
   !
   ! [tlu_variance]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  var_ten , var_ten_n    !< @public Noise diffusion tensor components
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  var_tau                !< @public wind Noise diffusion tensor components
   ! [tlu_variance]
   !
   ! [tlu_itoStokes]                                          !!  initial !   now    
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   uisd_0 ,   uisd_n     !< @public Ito-Stokes drift: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   visd_0 ,   visd_n     !< @public Ito-Stokes drift: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   wisd_0 ,   wisd_n     !< @public Ito-Stokes drift: z component
   ! [tlu_itoStokes]
   !
   ! [tlu_bias]                                               !!  initial !  after  ! 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ubia_0 ,   ubia_n     !< @public Modified advection: x component BIAS
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vbia_0 ,   vbia_n     !< @public Modified advection: y component BIAS 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   wbia_0 ,   wbia_n     !< @public Modified advection: z component BIAS
   REAL(wp),    PUBLIC,              SAVE                     ::  biaSIGN
   ! [tlu_bias]
   !
   ! [tlu_w_correction]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tlu_wcorr               !< @public Incompressibility correction
   ! [tlu_w_correction]
   !
   ! [tlu_trend]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   utnd                  !< @public LU trend
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vtnd                  !< @public LU trend (for removal in other parts)
   ! [tlu_trend]
   !
   !
   INTEGER(i4), PUBLIC, PARAMETER,         DIMENSION(3,3)     :: ndiffidx  = &             !: the indices of tensor var_ten 4th dim
                                                              &  RESHAPE((/ 1,4,5, 4,2,6, 5,6,3 /), (/ 3, 3 /))
   !
   !
   INTEGER,     PUBLIC            ::   dt_delay = 0 * 360 * 86400                       ! 
   !
   INTEGER,     PUBLIC, PARAMETER ::   np_ucmp = 1                                   ! to calculate U contributions
   INTEGER,     PUBLIC, PARAMETER ::   np_vcmp = 2                                   ! to calculate V contributions
   ! [public_vars]
   !
   ! [tlu_Schmidt]
   REAL(wp),    PUBLIC, PARAMETER ::   Sc_T = 1._wp   ! Schmidt number for Temperature
   REAL(wp),    PUBLIC, PARAMETER ::   Sc_S = 1._wp   ! Schmidt number for Salinity
   ! [tlu_Schmidt]
   !
   ! [tlu_MKE_prcntg]
   REAL(wp),    PUBLIC            ::   MKE_prcntg
   ! [tlu_MKE_prcntg]

   !! * Substitutions [TODO] enable substitutions for NEMO
   ! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------

CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_init ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Initialize module baroclinic variables
   !!
   !!
   !! 
   !! @result        Arrays allocated
   !!               
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_init 
   ! [tlu_init]
   SUBROUTINE tlu_init
      USE tluprj, ONLY : tlu_proj_init
      INTEGER(i4) :: ios, chk, ierr(5)   ! namelist output, allocation statistics
      !
      NAMELIST/namtlu/ ln_tlu,      &
      &                ln_pyp,      &
      &                ln_tlu_nke,  &
      &                ln_tlu_bia,  &
      &                ln_tlu_pod,  &
      &                ln_tlu_dmd,  &
      &                ln_tlu_pso,  &
      &                ln_tlu_wlt,  &
      &                ln_tlu_gss,  &
      &                MKE_prcntg
      !!------------------------------------------------------------------------
      !
      ! Read namelist: transport under location uncertainty parametrization
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !                                          ! [TODO] also read from numnam_cfg
      !
      ! Control print
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_init : Transport under Location Uncertainty '
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu'
         WRITE(numout,*) '      transport under location uncertainty       ln_tlu = ', ln_tlu
         WRITE(numout,*) '             isopycnal projection of noise       ln_pyp = ', ln_pyp
         WRITE(numout,*) '       noise rescaling based on kin. enrgy   ln_tlu_nke = ', ln_tlu_nke
         WRITE(numout,*) '                       Girsanov correction   ln_tlu_bia = ', ln_tlu_bia
         WRITE(numout,*)
         WRITE(numout,*) '   Reduced Order Model '
         WRITE(numout,*) '           Proper Orthogonal Decomposition   ln_tlu_pod = ', ln_tlu_pod
         WRITE(numout,*) '              Dynamical Mode Decomposition   ln_tlu_dmd = ', ln_tlu_dmd
         WRITE(numout,*) '        Data-Free Pseudo-Observation Model   ln_tlu_pso = ', ln_tlu_pso
         WRITE(numout,*) '        Data-Free Daubechies Wavelet Model   ln_tlu_wlt = ', ln_tlu_wlt
         WRITE(numout,*) '            Data-Free Gaussian White noise   ln_tlu_gss = ', ln_tlu_gss
         WRITE(numout,*)
         WRITE(numout,*) '   Schmidt numbers employed '
         WRITE(numout,*) '                              Temperature          Sc_T = ', Sc_T
         WRITE(numout,*) '                                 Salinity          Sc_S = ', Sc_S

      END IF
      !
      ierr = 0
      IF ( ln_tlu ) THEN      
         !
         ! Check multiple method 
         !
         chk = 0
         IF ( ln_tlu_pod ) chk = chk+1
         IF ( ln_tlu_dmd ) chk = chk+1
         IF ( ln_tlu_pso ) chk = chk+1
         IF ( ln_tlu_wlt ) chk = chk+1
         IF ( ln_tlu_gss ) chk = chk+1
         !
         IF ( chk .lt. 1) THEN
            WRITE(numout,*) ' E R R O R! No model chosen for noise formulation! ', ierr
            STOP
         ELSEIF ( chk .gt. 1) THEN
            WRITE(numout,*)      
            WRITE(numout,*) ' W A R N I N G! Multiple models chosen for noise formulation! ', ierr
            WRITE(numout,*)     
         END IF
         !
         ! Allocate current time variables (all dynamics refers to these variables)
         !
         ALLOCATE(      unoi(jpi,jpj,jpk),   &
         &              vnoi(jpi,jpj,jpk),   &
         &              wnoi(jpi,jpj,jpk),   &
         &            uisd_n(jpi,jpj,jpk),   &
         &            visd_n(jpi,jpj,jpk),   &
         &            wisd_n(jpi,jpj,jpk),   &
         &           var_ten(jpi,jpj,jpk,6), &
         &         tlu_wcorr(jpi,jpj,jpk),   &
         &              utnd(jpi,jpj,jpk),   &
         &              vtnd(jpi,jpj,jpk),   stat = ierr(1) )   
         !
         ! Initialization
         !
         unoi = 0._wp
         vnoi = 0._wp
         wnoi = 0._wp
         uisd_n = 0._wp
         visd_n = 0._wp
         wisd_n = 0._wp 
         var_ten = 0._wp
         utnd = 0._wp
         vtnd = 0._wp
         !
         IF ( ln_tlu_bia ) THEN
            !
            ! Allocate current time bias
            !
            ALLOCATE( ubia_n(jpi,jpj,jpk), &
            &         vbia_n(jpi,jpj,jpk), &
            &         wbia_n(jpi,jpj,jpk), stat = ierr(2) )  
            !
            ubia_n = 0._wp
            vbia_n = 0._wp
            wbia_n = 0._wp
            !
         END IF
         !
         IF ( ln_pyp ) THEN
            !
            CALL tlu_proj_init         
            !
            ! Allocate initial fields (starting point for projection)
            !
            ALLOCATE( uisd_0(jpi,jpj,jpk), &
            &         visd_0(jpi,jpj,jpk), &
            &         wisd_0(jpi,jpj,jpk), stat = ierr(3) )   
            !
            uisd_0 = 0._wp
            visd_0 = 0._wp
            wisd_0 = 0._wp 
            !
            IF ( ln_tlu_bia ) THEN
               !
               ! Allocate current time bias
               !
               ALLOCATE( ubia_0(jpi,jpj,jpk), &
               &         vbia_0(jpi,jpj,jpk), &
               &         wbia_0(jpi,jpj,jpk), stat = ierr(4) )  
               !
               ubia_0 = 0._wp
               vbia_0 = 0._wp
               wbia_0 = 0._wp
               !
            END IF
            !    
         END IF    
         !
         ! Set up the delay
         !
         dt_delay = dt_delay / rn_rdt
         !
         ! var_tau works as a top boundary conditon for the variance in the dynamics. 
         !
         ALLOCATE( var_tau(jpi,jpj), stat = ierr(5) )
         var_tau = 0._wp 
         !
         ! Allocation check
         !
         IF (SUM(ierr)/=0) THEN   
            WRITE(numout,*) ' tlu_init(): allocation failed = ', ierr
            STOP
         END IF
         !
      END IF
      !
   END SUBROUTINE tlu_init
   ! [tlu_init]

END MODULE tlu


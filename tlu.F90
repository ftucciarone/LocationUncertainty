!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu
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
!> @brief         Transport under Location Uncertainty module.
!! 
!! 
!! @par           Procedure specifics      
!> @details       Defines variables and dependencies for the stochastic parametrization of 
!!                Navier-Stokes equations. 
!!
!> @snippet       this tlu_noise
!!                \f{align*}{
!!                  \texttt{unoi} &\leftarrow d\boldsymbol{\eta}_{t,x} \\
!!                  \texttt{vnoi} &\leftarrow d\boldsymbol{\eta}_{t,y} \\
!!                  \texttt{wnoi} &\leftarrow d\boldsymbol{\eta}_{t,z} 
!!               \f}
!!
!> @snippet      this tlu_bias
!!               \f{align*}{
!!                  \texttt{ubias} &\leftarrow d\boldsymbol{\eta}_{t,x}^{b} \\
!!                  \texttt{vbias} &\leftarrow d\boldsymbol{\eta}_{t,y}^{b} \\
!!                  \texttt{wbias} &\leftarrow d\boldsymbol{\eta}_{t,z}^{b} 
!!               \f}
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
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
!> @todo          read namelist from _cfg, not only _ref
!! @todo          write error messages to appropriate unit
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
!  USE mockup_nemo      ! [NEMO] to be replaced by relevant actual NEMO modules:  wrk_nemo
!  USE wrk_nemo         ! Memory allocation
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
   !
   LOGICAL,     PUBLIC                                        :: ln_tlu = .FALSE.          !< @public Switch for Location Uncertainty (LU)
   LOGICAL,     PUBLIC                                        :: ln_isd = .FALSE.          !< @public Switch for Ito-Stokes drift dynamics
   LOGICAL,     PUBLIC                                        :: ln_wnd = .TRUE.           !< @public Switch for Wind information
   !
   ! [tlu_noise]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: unoi                      !< @public Modified advection: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: vnoi                      !< @public Modified advection: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wnoi                      !< @public Modified advection: z component
   ! [tlu_noise]
   !
   ! [tlu_itoStokes]                                          !! before ! now    ! after  
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: uisd_b , uisd_n , uisd_a  !< @public Ito-Stokes drift: x component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: visd_b , visd_n , visd_a  !< @public Ito-Stokes drift: y component
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wisd_b , wisd_n , wisd_a  !< @public Ito-Stokes drift: z component
   ! [tlu_itoStokes]
   !
   ! [tlu_bias]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: ubia                      !< @public Modified advection: x component BIAS
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: vbia                      !< @public Modified advection: y component BIAS 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wbia                      !< @public Modified advection: z component BIAS
   ! [tlu_bias]
   !
   ! [tlu_w_correction]
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tlu_wcorr                 !< @public Incompressibility correction
   ! [tlu_w_correction]
   !
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: var_tau                   !< @public Noise diffusion tensor components
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: var_ten                   !< @public Noise diffusion tensor components
   INTEGER(i4), PUBLIC, PARAMETER,         DIMENSION(3,3)     :: ndiffidx  = &             !: the indices of tensor var_ten 4th dim
                                                              &  RESHAPE((/ 1,4,5, 4,2,6, 5,6,3 /), (/ 3, 3 /))
   INTEGER,     PUBLIC, PARAMETER ::   np_ucmp = 1                                         ! to calculate U contributions
   INTEGER,     PUBLIC, PARAMETER ::   np_vcmp = 2                                         ! to calculate V contributions
   ! [public_vars]
   !
   ! [tlu_Schmidt]
   REAL(wp),    PUBLIC, PARAMETER ::   SC = 1
   REAL(wp),    PUBLIC, PARAMETER ::   r1_SC = 1._wp / sqrt(50._wp)
   ! [tlu_Schmidt]
   !
   ! [tlu_Schmidt]
   REAL(wp),    PUBLIC            ::   biaSIGN
   ! [tlu_Schmidt]


   !! * Substitutions [TODO] enable substitutions for NEMO
   ! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------

CONTAINS

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> P. Derian, P. Chandramouli, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief   Initialize module variables
   !!
   !!  @details
   !!  Reads namelist, compute other parameters and allocate arrays.
   !! 
   !!
   !! 
   !! 
   !> @param[in]     ln_tlu: logical switch for location uncertainty
   !! @param[in]     jpi: the first dimension of the spatial arrays
   !! @param[in]     jpj: the first dimension of the spatial arrays
   !! @param[in]     jpk: the first dimension of the spatial arrays
   !! @param[in]     numout: writing unit 
   !! @result        Arrays allocated
   !! @warning  
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_init 
   ! [tlu_init]
   SUBROUTINE tlu_init
      INTEGER(i4) :: ios, ierr(2)   ! namelist output, allocation statistics
      !
      NAMELIST/namtlu/ ln_tlu, ln_isd, ln_wnd
      !!------------------------------------------------------------------------
      !
      ! Read namelist: transport under location uncertainty parametrization
      !------------------------------
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !                                          ! [TODO] also read from numnam_cfg
      !
      ! Control print
      !------------------------------
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_init : Transport under Location Uncertainty '
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu'
         WRITE(numout,*) '      transport under location uncertainty   ln_tlu = ', ln_tlu
      END IF
      !
      IF ( ln_tlu ) THEN
         !
         ! Allocate
         !------------------------------
         ALLOCATE( unoi(jpi,jpj,jpk),      &
             &     vnoi(jpi,jpj,jpk),      &
             &     wnoi(jpi,jpj,jpk),      &
             &     var_tau(jpi,jpj),       &
             &     var_ten(jpi,jpj,jpk,6), &
             &     ubia(jpi,jpj,jpk),      &
             &     vbia(jpi,jpj,jpk),      &
             &     wbia(jpi,jpj,jpk),      &
             &     tlu_wcorr(jpi,jpj,jpk), stat=ierr(1) )    
         IF ( ln_isd ) THEN
            !
            ! Allocate BEFORE, NOW and AFTER fields
            !------------------------------
            ALLOCATE( uisd_b(jpi,jpj,jpk),      &
                &     visd_b(jpi,jpj,jpk),      &
                &     wisd_b(jpi,jpj,jpk),      &
                &     uisd_n(jpi,jpj,jpk),      &
                &     visd_n(jpi,jpj,jpk),      &
                &     wisd_n(jpi,jpj,jpk),      &
                &     uisd_a(jpi,jpj,jpk),      &
                &     visd_a(jpi,jpj,jpk),      &
                &     wisd_a(jpi,jpj,jpk), stat=ierr(2) )   
         ELSE
            !
            ! Allocate NOW fields
            !------------------------------
            ALLOCATE( uisd_n(jpi,jpj,jpk),      &
                &     visd_n(jpi,jpj,jpk),      &
                &     wisd_n(jpi,jpj,jpk), stat=ierr(2) )           
         END IF        
         !
         ! Allocation check
         !------------------------------   
         IF (SUM(ierr)/=0) THEN   
            WRITE(numout,*) ' tlu_init(): allocation failed = ', ierr(1), ierr(2)  
            STOP
         END IF
         !
         ! For now, var_tau better be zero
         !------------------------------   
         var_tau = 0._wp

      END IF
   END SUBROUTINE tlu_init
   ! [tlu_init]


END MODULE tlu




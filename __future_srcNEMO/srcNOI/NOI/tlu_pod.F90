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
MODULE tlupod
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
   USE tlu
   !   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_init_pod    ! called by tlu_init in tlu.f90
   !
   INCLUDE 'mpif.h'
   !
   !
   ! [tlu_spmodes] 
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spx_mod                 !> @public   U-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spy_mod                 !> @public   V-velocity modes
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: spz_mod                 !> @public   W-velocity modes
   ! [tlu_spmodes]
   !
   ! [tlu_randgauss] 
   REAL(wp),    PUBLIC, ALLOCATABLE,       DIMENSION(:)     :: brwn_rv                 !> @public   Gaussian Random Variables
   ! [tlu_randgauss]
   !
   ! [tlu_ken_modes]
   INTEGER(i4), PUBLIC                                      :: nn_tlu_nmod             !> @public   Number of ocean POD modes (static) 
   ! [tlu_ken_modes]
   !
   ! [tlu_kenspect]
   CHARACTER(lc), PUBLIC                                    :: nm_mod_nc               !> @public   Name of the NetCDF file containing modes
   ! [tlu_kenspect]

CONTAINS

   SUBROUTINE tlu_init_pod
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_svd_init  ***
      !!
      !! ** Purpose :   Initialize module variables.
      !!
      !! ** Method  :   Read namelist, compute other parameters and allocate arrays.
      !!
      !! ** Action  :   Arrays allocated.
      !!
      !!------------------------------------------------------------------------
      INTEGER(i4) :: ios, ierr(4), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Read namelist: transport under location uncertainty parametrization
      !
      NAMELIST/namtlu_pod/   nm_mod_nc,    & 
                         &   nn_tlu_nmod,  &
                         &   biaSIGN
      !
      ! Read namelist
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namtlu_pod, IOSTAT = ios ) ! [TODO] error management - cf module_example
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namtlu_pod, IOSTAT = ios ) ! [TODO] error management - cf module_example
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_pod_init : TLU, Proper orthogonal decomposition  noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '      data-based noise model            ln_tlu_pod = ', ln_tlu_pod
         WRITE(numout,*) '      number of modes employed         nn_tlu_nmod = ', nn_tlu_nmod
      END IF
      !
      ! Allocate Spatial modes
      !
      ierr = 0
      !
      ALLOCATE( spx_mod(jpi,jpj,nn_tlu_nmod * jpk), &
      &         spy_mod(jpi,jpj,nn_tlu_nmod * jpk), &
      &         spz_mod(jpi,jpj,nn_tlu_nmod * jpk),  stat=ierr(1) )   ! the singular values
      !
      !
      ! Allocate Brownian motion
      !
      ALLOCATE( brwn_rv(nn_tlu_nmod), stat=ierr(2) ) 
      !
      ! Read Baroclinic modes
      !
      CALL read_vel_mod( )


      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_pod_init(): allocation failed = ', ierr  
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !
      !
   END SUBROUTINE tlu_init_pod


   SUBROUTINE read_vel_mod( )
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE read_spa_mod  ***
      !!
      !! ** Purpose :   Read from file the spatial bases for forming noise
      !!
      !! ** Method :    Read from pre-computed file the bases.
      !!
      !! ** Action :    spx/y/z_mod contains the three component spatially organised
      !!                modes with the total number specified from the namelist
      !!
      !! ** Note:       Reading the variable using fortran read may not work!
      !!                For mode matching, make sure to normalize spatial basis
      !!                with singular value
      !!------------------------------------------------------------------------
      !
      INTEGER ::   jk, jkk, jm, m_idx, ierr ! dummy loop arguments
      INTEGER ::   endloop
      INTEGER ::   str_at

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwa1,   zwa2,   zwa3

      INTEGER                                 :: ucid,   vcid,   wcid
      CHARACTER (len=90)                      :: ubas,   vbas,   wbas

      REAL(wp)  :: val

      !
      ! Allocate dummy matrices
      !
      ALLOCATE( zwa1(jpi,jpj,jpk), &
      &         zwa2(jpi,jpj,jpk), &
      &         zwa3(jpi,jpj,jpk),  stat=ierr ) 

      IF ( ierr /= 0 ) THEN  
         WRITE(numout,*) ' read_spa_mod2(): allocation failed = ', ierr   
         STOP
      END IF  
      !
      ! Get netCDF file ID
      !
      str_at = INDEX(nm_mod_nc, '@')
      CALL iom_open(nm_mod_nc(1:str_at-1)//'u'//nm_mod_nc(str_at+1:), ucid )
      CALL iom_open(nm_mod_nc(1:str_at-1)//'v'//nm_mod_nc(str_at+1:), vcid )
      CALL iom_open(nm_mod_nc(1:str_at-1)//'w'//nm_mod_nc(str_at+1:), wcid )

      val = 1.0_wp
      endloop = 0 

      spx_mod = 0._wp
      spy_mod = 0._wp
      spz_mod = 0._wp
         
      DO jm = 1, nn_tlu_nmod
         !
         ! Initialize dummy matrices
         !
         zwa1 = 0._wp
         zwa2 = 0._wp
         zwa3 = 0._wp
         !
         ! Write variable name 
         !
         WRITE (ubas, '(a, i3.3)' )  "spat_basis_u_", jm
         WRITE (vbas, '(a, i3.3)' )  "spat_basis_v_", jm
         WRITE (wbas, '(a, i3.3)' )  "spat_basis_w_", jm
         !
         ! Read modes (in X, Y, Z)
         !
         CALL iom_get(ucid, jpdom_autoglo, ubas, zwa1, ldxios = lrxios)
         CALL iom_get(vcid, jpdom_autoglo, vbas, zwa2, ldxios = lrxios)
         CALL iom_get(wcid, jpdom_autoglo, wbas, zwa3, ldxios = lrxios)
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         !
         spx_mod(:,:,m_idx + 1 : m_idx + jpk ) = zwa1 * umask * val
         spy_mod(:,:,m_idx + 1 : m_idx + jpk ) = zwa2 * vmask * val
         spz_mod(:,:,m_idx + 1 : m_idx + jpk ) = zwa3 * wmask * val
         !
      ENDDO
      !
      ! Read bias (if flagged)
      !
      IF (ln_tlu_bia) THEN
         CALL iom_get(ucid,jpdom_autoglo,'u_mean',ubia_0, ldxios = lrxios)
         CALL iom_get(vcid,jpdom_autoglo,'v_mean',vbia_0, ldxios = lrxios)
         CALL iom_get(wcid,jpdom_autoglo,'w_mean',wbia_0, ldxios = lrxios)
         ! 
         ! If one rescales the noise, must rescale the bias if projected
         !
         ubia_0 = ubia_0 * umask
         vbia_0 = vbia_0 * vmask
         wbia_0 = wbia_0 * wmask
      END IF
      !
      ! Close netCDF files
      !
      CALL iom_close(ucid)
      CALL iom_close(vcid)
      CALL iom_close(wcid)
      !
      ! Deallocate dummy matrices
      !
      DEALLOCATE(zwa1, zwa2, zwa3)
      !
   END SUBROUTINE read_vel_mod

END MODULE tlupod

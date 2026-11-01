!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: task_data
!!   Coarse Graining procedure to create successful noise structures in the LU 
!!   module
!!   
!
!> @authors
!>                L. Li
!!
!> @version
!!                2020 -  ( L. LI           )  .f77 original code <br>
!!                2021 -  ( F. TUCCIARONE   )  .F90 adaptation
!! 
!> @brief         Allocate state variables        
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!
!!------------------------------------------------------------------------------
MODULE HR_fields
   !
   USE class_precision
   !  
   IMPLICIT NONE
   PUBLIC
   !
   !  Array dimensions    
   !  ----------------
   INTEGER, SAVE :: nxHR  ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER, SAVE :: nyHR  ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER, SAVE :: nloHR   ! number of oceanic layers
   INTEGER, SAVE :: ntoHR   ! number of snapshots (time counter)
   !
   REAL(KIND = wp) :: HRgridsize
   !
   !  Input fields    
   !  ------------
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   latHR  ! High resolution latitude
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   lonHR  ! High resolution longitude
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: b_levHR  ! High resolution bottom level height
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: t_levHR  ! High resolution top level height
   ! 
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   uHR_b  ! High resolution zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   vHR_b  ! High resolution meridional velocity  !
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     uHR  ! High resolution zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     vHR  ! High resolution meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     wHR  ! High resolution vertical velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   sshHR  ! High resolution vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanuHR  ! High resolution mean zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanvHR  ! High resolution mean meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanwHR  ! High resolution mean vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctuHR  ! High resolution fluctuation of zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctvHR  ! High resolution fluctuation of meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctwHR  ! High resolution fluctuation of vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     tHR  ! High resolution temperature
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     sHR  ! High resolution salinity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1uHR  ! High resolution dx scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1vHR  ! High resolution dx scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1tHR  ! High resolution dx scale factor at t point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2uHR  ! High resolution dy scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2vHR  ! High resolution dy scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2tHR  ! High resolution dy scale factor at t point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3uHR  ! High resolution dz scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3vHR  ! High resolution dz scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3wHR  ! High resolution dz scale factor at w point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3tHR  ! High resolution dz scale factor at t point
   !
   !  Filtered fields   
   !  ---------------
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: uf1_HR, uf2_HR, uband_HR  ! filtered zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: vf1_HR, vf2_HR, vband_HR  ! filtered meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: wf1_HR, wf2_HR, wband_HR  ! filtered vertical velocity
   !
   !  NetCDF identifiers    
   !  ------------------
   INTEGER, PUBLIC :: HRoutID

   INTEGER, PUBLIC ::  HRlatID
   INTEGER, PUBLIC ::  HRlonID

   INTEGER, PUBLIC ::   HRtvID
   INTEGER, PUBLIC ::   HRzvID
   INTEGER, PUBLIC ::   HRyvID
   INTEGER, PUBLIC ::   HRxvID 

   INTEGER, PUBLIC :: uf1ID_HR, uf2ID_HR, ubID_HR
   INTEGER, PUBLIC :: vf1ID_HR, vf2ID_HR, vbID_HR
   INTEGER, PUBLIC :: wf1ID_HR, wf2ID_HR, wbID_HR
   !
   !  Public routines
   !  ---------------
   PUBLIC ::  HR_task_init
   PUBLIC ::  HR_state_init
   !
CONTAINS
 
   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Initialize identifiers of input variables and read stationary
   !!            states 
   !!
   !!  @details
   !!  Reads namelist, compute other parameters and allocate arrays.
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE HR_task_init (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt   ! Parent subroutine name (for EH)
      INTEGER                            :: ierr(3)

      IF ( ( TRIM(opt) /= 'bar' ) .AND. &
         & ( TRIM(opt) /= 'vel' ) .AND. &
         & ( TRIM(opt) /= 'avg' ) .AND. &
         & ( TRIM(opt) /= 'flc' ) .AND. &
         & ( TRIM(opt) /= 'dom' ) .AND. &
         & ( TRIM(opt) /= 'all' )       ) THEN
         print *, 'Input string must be either bar, vel, avg, flc, dom or all'
         STOP
      ENDIF
      !

      IF (( TRIM(opt) == 'bar' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                          &
            &      uHR_b(nxHR,nyHR),        &
            &      vHR_b(nxHR,nyHR),        &
            &                                STAT=ierr(1) )
      END IF

      IF (( TRIM(opt) == 'vel' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                          &
            &        uHR(nxHR,nyHR,nloHR),  &
            &        vHR(nxHR,nyHR,nloHR),  &          
            &        wHR(nxHR,nyHR,nloHR),  &        
            &      sshHR(nxHR,nyHR),        &
            &                                STAT=ierr(1) )
      END IF

      IF (( TRIM(opt) == 'avg' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                          &
            &  meanuHR(nxHR,nyHR,nloHR),    &
            &  meanvHR(nxHR,nyHR,nloHR),    &          
            &  meanwHR(nxHR,nyHR,nloHR),    &
            &                                STAT=ierr(1) )
      END IF

      IF (( TRIM(opt) == 'flc' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                          &
            &  flctuHR(nxHR,nyHR,nloHR),    &
            &  flctvHR(nxHR,nyHR,nloHR),    &          
            &  flctwHR(nxHR,nyHR,nloHR),    &
            &                                STAT=ierr(1) )
      END IF
      !
      IF (( TRIM(opt) == 'dom' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                          &
            &      e1uHR(nxHR,nyHR),        &
            &      e1vHR(nxHR,nyHR),        &
            &      e1tHR(nxHR,nyHR),        &
            &      e2uHR(nxHR,nyHR),        & 
            &      e2vHR(nxHR,nyHR),        &     
            &      e2tHR(nxHR,nyHR),        &    
            &      e3uHR(nxHR,nyHR,nloHR),  &
            &      e3vHR(nxHR,nyHR,nloHR),  &
            &      e3wHR(nxHR,nyHR,nloHR),  &
            &      e3tHR(nxHR,nyHR,nloHR),  &
            &                                STAT=ierr(2) )
            !     
         ALLOCATE(                        &
            &        latHR(nxHR,nyHR),    &
            &        lonHR(nxHR,nyHR),    &
            &      b_levHR(nxHR,nyHR),    &
            &      t_levHR(nxHR,nyHR),    &
            &                                STAT=ierr(3) )
      END IF
      !
      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in task_initHR()'
      END IF
   END SUBROUTINE HR_task_init
   

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Initialize identifiers of input variables and read stationary
   !!            states 
   !!
   !!  @details
   !!  Reads namelist, compute other parameters and allocate arrays.
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE HR_state_init (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt
      INTEGER :: ierr(1)

      IF (( TRIM(opt) /= 'band' ) .AND. ( TRIM(opt) /= '' ) ) THEN
         print *, 'Input string must be either band or empty'
         STOP
      ENDIF
      !
      ALLOCATE( uf1_HR(nxHR,nyHR,nloHR),    &
         &      vf1_HR(nxHR,nyHR,nloHR),    &          
         &      wf1_HR(nxHR,nyHR,nloHR),    &
         &      uf2_HR(nxHR,nyHR,nloHR),    &
         &      vf2_HR(nxHR,nyHR,nloHR),    &          
         &      wf2_HR(nxHR,nyHR,nloHR),    &
         &    uband_HR(nxHR,nyHR,nloHR),    &
         &    vband_HR(nxHR,nyHR,nloHR),    &
         &    wband_HR(nxHR,nyHR,nloHR),    &
         &                             STAT=ierr(1) )

      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in state_init()'
      END IF
   END SUBROUTINE HR_state_init

END MODULE HR_fields   


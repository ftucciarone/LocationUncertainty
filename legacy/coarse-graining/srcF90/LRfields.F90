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
MODULE LR_fields
   !
   USE class_precision
   !  
   IMPLICIT NONE
   PUBLIC
   !
   !  Array dimensions    
   !  ----------------
   INTEGER, SAVE :: nxLR   ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER, SAVE :: nyLR   ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER, SAVE :: nloLR  ! number of oceanic layers
   INTEGER, SAVE :: ntoLR  ! number of snapshots (time counter)
   !
   !  Input fields    
   !  ----------
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   latLR  ! Low resolution latitude
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   lonLR  ! Low resolution longitude
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: b_levLR  ! Low resolution bottom level height
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: t_levLR  ! Low resolution top level height
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   uLR_b  ! Low resolution zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   vLR_b  ! Low resolution meridional velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     uLR  ! Low resolution zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     vLR  ! Low resolution meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     wLR  ! Low resolution vertical velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   sshLR  ! Low resolution sea surface height
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanuLR  ! Low resolution mean zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanvLR  ! Low resolution mean meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanwLR  ! Low resolution mean vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctuLR  ! Low resolution fluctuation of zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctvLR  ! Low resolution fluctuation of meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: flctwLR  ! Low resolution fluctuation of vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     tLR  ! Low resolution temperature
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     sLR  ! Low resolution salinity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1uLR  ! Low resolution dx scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1vLR  ! Low resolution dx scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1tLR  ! Low resolution dx scale factor at t point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2uLR  ! Low resolution dy scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2vLR  ! Low resolution dy scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2tLR  ! Low resolution dy scale factor at t point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3uLR  ! Low resolution dz scale factor at u point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3vLR  ! Low resolution dz scale factor at v point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3wLR  ! Low resolution dz scale factor at w point
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3tLR  ! Low resolution dz scale factor at t point
   !
   !  Filtered fields   
   !  ---------------
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   uf1_LR,   uf2_LR,   uband_LR  ! filtered zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   vf1_LR,   vf2_LR,   vband_LR  ! filtered meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   wf1_LR,   wf2_LR,   wband_LR  ! filtered vertical velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: sshf1_LR, sshf2_LR, sshband_LR  ! filtered vertical velocity
   !
   !  NetCDF identifiers    
   !  ------------------
   INTEGER, PUBLIC ::  LRoutID

   INTEGER, PUBLIC ::  LRlatID
   INTEGER, PUBLIC ::  LRlonID

   INTEGER, PUBLIC ::   LRtvID
   INTEGER, PUBLIC ::   LRzvID
   INTEGER, PUBLIC ::   LRyvID
   INTEGER, PUBLIC ::   LRxvID
 
   INTEGER, PUBLIC ::   uID_LR
   INTEGER, PUBLIC ::   vID_LR
   INTEGER, PUBLIC ::   wID_LR
   INTEGER, PUBLIC :: sshID_LR

   INTEGER, PUBLIC :: u_btID_LR
   INTEGER, PUBLIC :: v_btID_LR

   INTEGER, PUBLIC :: meanuID_LR
   INTEGER, PUBLIC :: meanvID_LR
   INTEGER, PUBLIC :: meanwID_LR

   INTEGER, PUBLIC ::   uf1ID_LR,   uf2ID_LR,   ubID_LR
   INTEGER, PUBLIC ::   vf1ID_LR,   vf2ID_LR,   vbID_LR
   INTEGER, PUBLIC ::   wf1ID_LR,   wf2ID_LR,   wbID_LR
   INTEGER, PUBLIC :: sshf1ID_LR, sshf2ID_LR, sshbID_LR
   !
   !  Public routines
   !  ---------------
   PUBLIC ::  LR_task_init
   PUBLIC ::  LR_state_init
   !
CONTAINS
 

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE LR_task_init (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt   ! Parent subroutine name (for EH)
      INTEGER                            :: ierr(3)

      IF (( TRIM(opt) /= 'bar' ) .AND. ( TRIM(opt) /= 'vel' ) .AND. ( TRIM(opt) /= 'dom' ) .AND. ( TRIM(opt) /= 'all' ) ) THEN
         print *, 'Input string must be either bar, vel, dom or all'
         STOP
      ENDIF
      !
      IF (( TRIM(opt) == 'bar' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                             &
            &      uLR_b(nxLR,nyLR),           &
            &      vLR_b(nxLR,nyLR),           &
            &                                STAT=ierr(1) )
      END IF
      !
      IF (( TRIM(opt) == 'vel' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                             &
            &      uLR(nxLR,nyLR,nloLR),       &
            &      vLR(nxLR,nyLR,nloLR),       &
            &      wLR(nxLR,nyLR,nloLR),       &
            &    sshLR(nxLR,nyLR),             &
            &  meanuLR(nxLR,nyLR,nloLR),       &
            &  meanvLR(nxLR,nyLR,nloLR),       &
            &  meanwLR(nxLR,nyLR,nloLR),       &
            &                                STAT=ierr(1) )
      END IF
      !
      IF (( TRIM(opt) == 'dom' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                             &
            &      e1uLR(nxLR,nyLR),           &
            &      e1vLR(nxLR,nyLR),           &
            &      e1tLR(nxLR,nyLR),           &
            &      e2uLR(nxLR,nyLR),           &
            &      e2vLR(nxLR,nyLR),           &
            &      e2tLR(nxLR,nyLR),           &    
            &      e3uLR(nxLR,nyLR,nloLR),     &
            &      e3vLR(nxLR,nyLR,nloLR),     &
            &      e3wLR(nxLR,nyLR,nloLR),     &
            &      e3tLR(nxLR,nyLR,nloLR),     &
            &                                STAT=ierr(2) )
            !
         ALLOCATE(                           &
            &        latLR(nxLR,nyLR),       &
            &        lonLR(nxLR,nyLR),       &
            &      b_levLR(nxLR,nyLR),       &
            &      t_levLR(nxLR,nyLR),       &
            &                                STAT=ierr(3) )
      END IF
      !
      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in task_initLR()'
      END IF
   END SUBROUTINE LR_task_init

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE LR_state_init (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt
      INTEGER :: ierr(1)

      IF (( TRIM(opt) /= 'band' ) .AND. ( TRIM(opt) /= '' ) ) THEN
         print *, 'Input string must be either band or empty'
         STOP
      ENDIF
      !
      ALLOCATE(     uf1_LR(nxLR,nyLR,nloLR), &
         &          vf1_LR(nxLR,nyLR,nloLR), &          
         &          wf1_LR(nxLR,nyLR,nloLR), &          
         &        sshf1_LR(nxLR,nyLR),       &
         &          uf2_LR(nxLR,nyLR,nloLR), &
         &          vf2_LR(nxLR,nyLR,nloLR), &          
         &          wf2_LR(nxLR,nyLR,nloLR), &          
         &        sshf2_LR(nxLR,nyLR),       &
         &        uband_LR(nxLR,nyLR,nloLR), &
         &        vband_LR(nxLR,nyLR,nloLR), &
         &        wband_LR(nxLR,nyLR,nloLR), &
         &      sshband_LR(nxLR,nyLR),       &
         &                             STAT=ierr(1) )

      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in state_init()'
      END IF
   END SUBROUTINE LR_state_init

END MODULE LR_fields   


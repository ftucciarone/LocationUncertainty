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
MODULE task_data
   !
   USE class_precision
   !  
   IMPLICIT NONE
   PUBLIC
   !
   !  Array dimensions    
   !  ----------------
   INTEGER, SAVE :: nxR27  ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER, SAVE :: nyR27  ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER, SAVE :: nxR3   ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER, SAVE :: nyR3   ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER, SAVE :: nlo    ! number of oceanic layers
   INTEGER, SAVE :: nto    ! number of snapshots (time counter)
   !
   PUBLIC ::  task_initR3, task_initR27
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: uf, ur  ! coarse-grained and residual zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: vf, vr  ! coarse-grained and residual meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: wf, wr  ! coarse-grained and residual vertical velocity
   !  Additional netCDF identifiers
   INTEGER, PUBLIC :: ufID, urID
   INTEGER, PUBLIC :: vfID, vrID
   INTEGER, PUBLIC :: wfID, wrID


   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   latR27,   lonR27
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::    latR3,    lonR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: b_levR27, t_levR27
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::  b_levR3,  t_levR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     uR27,      uR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     vR27,      vR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     wR27,      wR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     tR27,      tR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     sR27,      sR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1uR27,    e1uR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1vR27,    e1vR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e1tR27,    e1tR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2uR27,    e2uR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2vR27,    e2vR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   e2tR27,    e2tR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3uR27,    e3uR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3vR27,    e3vR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3wR27,    e3wR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::   e3tR27,    e3tR3


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
   SUBROUTINE task_initR27 (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt   ! Parent subroutine name (for EH)
      INTEGER                            :: ierr(3)

      IF (( TRIM(opt) /= 'vel' ) .AND. ( TRIM(opt) /= 'dom' ) .AND. ( TRIM(opt) /= 'all' ) ) THEN
         print *, 'Input string must be either vel, dom or all'
         STOP
      ENDIF
      !
      IF (( TRIM(opt) == 'vel' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                           &
            &      uR27(nxR27,nyR27,nlo),    &
            &      vR27(nxR27,nyR27,nlo),    &          
            &      wR27(nxR27,nyR27,nlo),    &
            &                                STAT=ierr(1) )
      END IF
      !
      IF (( TRIM(opt) == 'dom' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                           &
            &      e1uR27(nxR27,nyR27),      &
            &      e1vR27(nxR27,nyR27),      &
            &      e1tR27(nxR27,nyR27),      &
            &      e2uR27(nxR27,nyR27),      & 
            &      e2vR27(nxR27,nyR27),      &     
            &      e2tR27(nxR27,nyR27),      &    
            &      e3uR27(nxR27,nyR27,nlo),  &
            &      e3vR27(nxR27,nyR27,nlo),  &
            &      e3wR27(nxR27,nyR27,nlo),  &
            &      e3tR27(nxR27,nyR27,nlo),  &
            &                                STAT=ierr(2) )
            !     
         ALLOCATE(                           &
            &      latR27(nxR27,nyR27),      &
            &      lonR27(nxR27,nyR27),      &
            &      b_levR27(nxR27,nyR27),    &
            &      t_levR27(nxR27,nyR27),    &
            &                                STAT=ierr(3) )
      END IF
      !
      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in task_initR27()'
      END IF
   END SUBROUTINE task_initR27
   
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
   SUBROUTINE task_initR3 (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt   ! Parent subroutine name (for EH)
      INTEGER                            :: ierr(3)

      IF (( TRIM(opt) /= 'vel' ) .AND. ( TRIM(opt) /= 'dom' ) .AND. ( TRIM(opt) /= 'all' ) ) THEN
         print *, 'Input string must be either vel, dom or all'
         STOP
      ENDIF
      !
      IF (( TRIM(opt) == 'vel' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                           &
            &      uR3(nxR3,nyR3,nlo),       &
            &      vR3(nxR3,nyR3,nlo),       &
            &      wR3(nxR3,nyR3,nlo),       &
            &                                STAT=ierr(1) )
      END IF
      !
      IF (( TRIM(opt) == 'dom' ) .OR. ( TRIM(opt) == 'all' )) THEN
         ALLOCATE(                           &
            &      e1uR3(nxR3,nyR3),         &
            &      e1vR3(nxR3,nyR3),         &
            &      e1tR3(nxR3,nyR3),         &
            &      e2uR3(nxR3,nyR3),         &
            &      e2vR3(nxR3,nyR3),         &
            &      e2tR3(nxR3,nyR3),         &    
            &      e3uR3(nxR3,nyR3,nlo),     &
            &      e3vR3(nxR3,nyR3,nlo),     &
            &      e3wR3(nxR3,nyR3,nlo),     &
            &      e3tR3(nxR3,nyR3,nlo),     &
            &                                STAT=ierr(2) )
            !
         ALLOCATE(                           &
            &      latR3(nxR3,nyR3),         &
            &      lonR3(nxR3,nyR3),         &
            &      b_levR3(nxR3,nyR3),       &
            &      t_levR3(nxR3,nyR3),       &
            &                                STAT=ierr(3) )
      END IF
      !
      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in task_initR3()'
      END IF
   END SUBROUTINE task_initR3

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
   SUBROUTINE state_init ()
      INTEGER :: ierr(1)

      ALLOCATE( uf(nxR3,nyR3,nlo),    &
         &      vf(nxR3,nyR3,nlo),    &          
         &      wf(nxR3,nyR3,nlo),    &
         &      ur(nxR3,nyR3,nlo),    &
         &      vr(nxR3,nyR3,nlo),    &
         &      wr(nxR3,nyR3,nlo),  STAT=ierr(1) )

      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in state_init()'
      END IF
   END SUBROUTINE state_init

END MODULE task_data   
!      
!**************************************************************************

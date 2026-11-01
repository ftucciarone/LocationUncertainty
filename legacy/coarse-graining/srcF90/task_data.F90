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
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: uf_b, ur_b  ! coarse-grained and residual barotropic zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: vf_b, vr_b  ! coarse-grained and residual barotropic meridional velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: uf1, uf2, uband  ! coarse-grained and residual zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: vf1, vf2, vband  ! coarse-grained and residual meridional velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: wf1, wf2, wband  ! coarse-grained and residual vertical velocity
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: uf1_b, uf2_b, uband_b  ! coarse-grained and residual barotropic zonal velocity
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: vf1_b, vf2_b, vband_b  ! coarse-grained and residual barotropic meridional velocity


   !  Additional netCDF identifiers
   INTEGER, PUBLIC :: ufID, urID
   INTEGER, PUBLIC :: vfID, vrID
   INTEGER, PUBLIC :: wfID, wrID
   !  Additional netCDF identifiers
   INTEGER, PUBLIC :: uf_bID, ur_bID
   INTEGER, PUBLIC :: vf_bID, vr_bID

   !  Additional netCDF identifiers
   INTEGER, PUBLIC :: uf1ID, uf2ID, ubID
   INTEGER, PUBLIC :: vf1ID, vf2ID, vbID
   INTEGER, PUBLIC :: wf1ID, wf2ID, wbID
   !  Additional netCDF identifiers
   INTEGER, PUBLIC :: uf1_bID, uf2_bID, ub_bID
   INTEGER, PUBLIC :: vf1_bID, vf2_bID, vb_bID


   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   latR27,   lonR27
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::    latR3,    lonR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: b_levR27, t_levR27
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::  b_levR3,  t_levR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     uR27,      uR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     vR27,      vR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     wR27,      wR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanuR27,  meanuR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanvR27,  meanvR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: meanwR27,  meanwR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     tR27,      tR3
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) ::     sR27,      sR3

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   uR27_b,    uR3_b
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   ::   vR27_b,    vR3_b

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
            &  meanuR27(nxR27,nyR27,nlo),    &
            &  meanvR27(nxR27,nyR27,nlo),    &          
            &  meanwR27(nxR27,nyR27,nlo),    &
            &    uR27_b(nxR27,nyR27),        &          
            &    vR27_b(nxR27,nyR27),        &
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
            &  meanuR3(nxR3,nyR3,nlo),       &
            &  meanvR3(nxR3,nyR3,nlo),       &
            &  meanwR3(nxR3,nyR3,nlo),       &
            &    uR3_b(nxR3,nyR3),           &
            &    vR3_b(nxR3,nyR3),           &
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
   SUBROUTINE state_init (opt)

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: opt
      INTEGER :: ierr(1)

      IF (( TRIM(opt) /= 'band' ) .AND. ( TRIM(opt) /= '' ) ) THEN
         print *, 'Input string must be either band or empty'
         STOP
      ENDIF
      !
      IF ( TRIM(opt) == 'band' ) THEN
         ALLOCATE( uf1(nxR3,nyR3,nlo),    &
            &      vf1(nxR3,nyR3,nlo),    &          
            &      wf1(nxR3,nyR3,nlo),    &
            &    uf1_b(nxR3,nyR3),        &          
            &    vf1_b(nxR3,nyR3),        &
            &      uf2(nxR3,nyR3,nlo),    &
            &      vf2(nxR3,nyR3,nlo),    &          
            &      wf2(nxR3,nyR3,nlo),    &
            &    uf2_b(nxR3,nyR3),        &          
            &    vf2_b(nxR3,nyR3),        &
            &    uband(nxR3,nyR3,nlo),    &
            &    vband(nxR3,nyR3,nlo),    &
            &    wband(nxR3,nyR3,nlo),    &
            &  uband_b(nxR3,nyR3),        &          
            &  vband_b(nxR3,nyR3),        &
            &                             STAT=ierr(1) )
      ELSE
         ALLOCATE( uf(nxR3,nyR3,nlo),    &
            &      vf(nxR3,nyR3,nlo),    &          
            &      wf(nxR3,nyR3,nlo),    &
            &    uf_b(nxR3,nyR3),        &          
            &    vf_b(nxR3,nyR3),        &
            &      ur(nxR3,nyR3,nlo),    &
            &      vr(nxR3,nyR3,nlo),    &
            &      wr(nxR3,nyR3,nlo),    &
            &    ur_b(nxR3,nyR3),        &          
            &    vr_b(nxR3,nyR3),        &
            &                             STAT=ierr(1) )
      END IF

      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in state_init()'
      END IF
   END SUBROUTINE state_init

END MODULE task_data   
!      
!**************************************************************************


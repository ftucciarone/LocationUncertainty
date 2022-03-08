!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: state
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
MODULE state
   !
   USE class_precision
   USE param
   !  
   IMPLICIT NONE
   !
   PUBLIC ::  state_init
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: uf, ur          ! coarse-grained and residual zonal velocity (m/s)
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: vf, vr          ! coarse-grained and residual meridional velocity (m/s)
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: wf, wr          ! coarse-grained and residual vertical velocity (m/s)
   !
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: uwindf, uwindr  ! coarse-grained and residual zonal wind stress (N/m^2)
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: vwindf, vwindr  ! coarse-grained and residual meridional velocity (N/m^2)

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: lat,lon
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: b_lev, t_lev

   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: e1v,e2u
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:)   :: e1t,e2t
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: e3u,e3v
   REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:) :: e3t 

 
!  TODO: [ Optional ]
!  Maybe you could complet the output data by defining properly
!  the NEMO grid: read the fine one from input data and subsample
!  them on coarse-grid, but these are only optional! 
!  We can do this latter :)
!  Here are my old version on a plane (much easier):
!  real xpo(nxpo),ypo(nypo),zo(nlo),tyrs(nto),xpoc(nxpoc),ypoc(nypoc)
!
!  xpo and ypo are fine-grid axis (km)
!  zo is vertical axis (km)
!  tyrs is time axis (years) 
!  xpoc and ypoc are coarse-grid axis (km)
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
   SUBROUTINE state_init ()
      INTEGER :: ierr(3)
      ALLOCATE( uf(nxpoc,nypoc,nlo),    &
         &      vf(nxpoc,nypoc,nlo),    &          
         &      wf(nxpoc,nypoc,nlo),    &
         &      ur(nxpoc,nypoc,nlo),    &
         &      vr(nxpoc,nypoc,nlo),    &
         &      wr(nxpoc,nypoc,nlo),    &    
         &      uwindf(nxpoc,nypoc),    &
         &      uwindr(nxpoc,nypoc),    &
         &      vwindf(nxpoc,nypoc),    &
         &      vwindr(nxpoc,nypoc),  STAT=ierr(1) )
         !     
      ALLOCATE( e2u(nxpoc,nypoc)    ,   &        
         &      e1v(nxpoc,nypoc)    ,   &
         &      e1t(nxpoc,nypoc)    ,   &
         &      e2t(nxpoc,nypoc)    ,   &
         &      e3u(nxpoc,nypoc,nlo),   &
         &      e3v(nxpoc,nypoc,nlo),   &  
         &      e3t(nxpoc,nypoc,nlo), STAT=ierr(2) )

      ALLOCATE( lat (nxpoc,nypoc) ,     &
         &      lon (nxpoc,nypoc) ,     &
         &      b_lev(nxpoc,nypoc),     &
         &      t_lev(nxpoc,nypoc),   STAT=ierr(3) )

      IF( sum(ierr) .ne. 0 ) THEN
         print *, 'Allocation error in state_init()'
      END IF
   END SUBROUTINE state_init
   

END MODULE state   
!      
!**************************************************************************

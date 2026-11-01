!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: param
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
!> @brief         Set the parameters of the coarse graining procedure       
!! 
!! 
!! @par           Procedure specifics      
!> @details       Defines the name of the input files, the names of the output,
!!                the name of the variables and the array dimensions
!!
!!------------------------------------------------------------------------------
MODULE param

   USE class_precision

   IMPLICIT NONE

   PUBLIC
   SAVE

!  I/O directory
!  -------------
   CHARACTER (len = lc)   :: iodir = ''  ! path of I/O data
   CHARACTER (len = lc)   :: trcnam  ! input NetCDF file of tracers
   CHARACTER (len = lc)   :: uncnam  ! input NetCDF file of zonal velocity
   CHARACTER (len = lc)   :: vncnam  ! input NetCDF file of meridional velocity
   CHARACTER (len = lc)   :: wncnam  ! input NetCDF file of vertical velocity
   CHARACTER (len = lc)   :: tggrid  ! input NetCDF file of target grid
   CHARACTER (len = lc)   :: outnam  ! output NetCDF file

!  Array dimensions    
!  ----------------
   INTEGER :: nxto  ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER :: nyto  ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER :: nlo   ! umber of oceanic layers
   INTEGER :: nto   ! number of snapshots (time counter)
   INTEGER :: nso   ! factor of both nxto and nyto (Ex. nso = 9 for data from R27 to R3) 
   INTEGER :: co    ! 

!  Derived grid parameters 
!  -------------------------------------- 
   INTEGER :: nxpo  ! = nxto + 2    oceanic complete grid with boundaries in x 
   INTEGER :: nypo  ! = nyto + 2    oceanic complete grid with boundaries in y
   INTEGER :: nxtoc ! = nxto/nso    coarse inner gridcells in x
   INTEGER :: nytoc ! = nyto/nso    coarse inner gridcells in y
   INTEGER :: nxpoc ! = nxtoc + 2   complete coarse grid points in x
   INTEGER :: nypoc ! = nytoc + 2   complete coarse grid points in y

   LOGICAL, PARAMETER :: rescale = .FALSE.


END MODULE param    

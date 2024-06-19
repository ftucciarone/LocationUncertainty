!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_newEnv
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Brief resume
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!                
!!                
!> @snippet       this tlu_noise
!!
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the second dimension of the spatial arrays
!! @param[in]     jpk: the third dimension of the spatial arrays
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
!> @todo          
!! @todo          
!!
!!------------------------------------------------------------------------------
MODULE tlunewEnv
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
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC tlu_newsub      ! Indicate which subroutine calls it
   !
   INCLUDE 'mpif.h'
   !
   REAL(wp),    PUBLIC    :: mod_var
   !
CONTAINS

   SUBROUTINE newsub
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_newsub  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!------------------------------------------------------------------------
      INTEGER(i4) :: ios, ierr(4), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_newsub : TLU, brief description to be printed in ocean.output' 
         WRITE(numout,*) '~~~~~~~~~~~'
      END IF
      ! 
      
      
      
      
      !
   END SUBROUTINE tlu_newsub


END MODULE tlunewEnv








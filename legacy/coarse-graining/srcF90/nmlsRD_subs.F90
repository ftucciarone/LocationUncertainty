!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: nmlsRD_subs
!!    
!!   
!!   
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2021 -  ( F. TUCCIARONE   )  .F90 original code <br>
!! 
!> @brief                 
!! 
!! 
!! @par                 
!> @details       
!!
!!------------------------------------------------------------------------------
MODULE nmlsRD_subs
  
   USE class_precision
   !
   IMPLICIT NONE

   INTEGER, PARAMETER :: stderr = 6 
   !
   PRIVATE

!  Subroutines      
   PUBLIC :: open_inputfile
   PUBLIC :: close_inputfile

CONTAINS

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    
   !!             
   !!
   !!  @details
   !!  
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE open_inputfile(file_path, file_unit, iostat)
      !! Check whether file exists, with consitent error message
      !! return the file unit
      CHARACTER(len = *),  INTENT( in  )  :: file_path
      INTEGER,              INTENT( out )  :: file_unit, iostat

      INQUIRE (FILe=file_path, IOSTAT=iostat)
      IF (iostat /= 0) THEN
         WRITE (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
      END IF
      OPEN (ACTION='read', FILE=file_path, IOSTAT=iostat, NEWUNIT=file_unit)
   END SUBROUTINE open_inputfile


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    
   !!             
   !!
   !!  @details
   !!  
   !! 
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   SUBROUTINE close_inputfile(file_path, file_unit, iostat)
      !! Check the reading was OK
      !! return error line IF not
      !! close the unit
      CHARACTER(len = lc),  INTENT( in  )  :: file_path
      CHARACTER(len=1000)                  :: line
      INTEGER,              INTENT( in  )  :: file_unit, iostat

      IF (iostat /= 0) THEN
         WRITE (stderr, '(2a)'      ) 'Error reading file :"', TRIM(file_path)
         WRITE (stderr, '(a, i0)'   ) 'iostat was:"', iostat
         BACKSPACE(file_unit)
         READ  (file_unit, FMT='(A)')  line
         WRITE (stderr, '(A)'       ) 'Invalid line : '//TRIM(line)
      END IF
      CLOSE (file_unit)   
   END SUBROUTINE close_inputfile

END MODULE nmlsRD_subs

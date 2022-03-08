!! ------------------------------------------------------------------------------
!! 
!!                               Module: class_precision
!!    Defines the precision of the code
!!    
!
!> @authors
!>                Carlo Janna, Massimiliano Ferronato and Nicola Castelletto (uniPD)
!!
!> @version
!!                2010 -  ( C.Janna, M.Ferronato, N.Castelletto ) original code <br>
!!                2020 -  ( F.Tucciarone ) adaptation to NEMO standards
!!
!> @brief         Set data precision codes        
!! 
!!
!!------------------------------------------------------------------------------

MODULE class_precision

   IMPLICIT NONE 

   !                                                                     !!** Floating point **
   INTEGER, PUBLIC, PARAMETER ::   single = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   double = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = single                           !: working precision

   !
   INTEGER, PUBLIC, PARAMETER ::   io_rp = single                        !: I/O reading precision
   INTEGER, PUBLIC, PARAMETER ::   io_wp = single                        !: I/O working precision

   !                                                                     !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   i4 = SELECTED_INT_KIND( 9)            !: single precision (integer 4)
   INTEGER, PUBLIC, PARAMETER ::   i8 = SELECTED_INT_KIND(14)            !: double precision (integer 8)
   
   !                                                                     !!** Integer **
   INTEGER, PUBLIC, PARAMETER ::   lc = 256                              !: Lenght of Character strings

END MODULE class_precision

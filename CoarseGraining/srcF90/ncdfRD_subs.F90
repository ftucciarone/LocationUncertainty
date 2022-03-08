!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: ncdfRD_subs
!!   This module is designed to contain all Ready-Made routines to read several
!!   data tyoes from a given netCDF tile.
!!   
!
!> @authors
!>                F. TUCCIARONE 
!!
!> @version
!!                2021 -  ( F. TUCCIARONE   )  Original code
!! 
!> @brief         Routines to read netCDF files        
!! 
!! 
!! @par           Procedure specifics      
!> @details       This module is used to read netCDF files in a standardised way.
!!                Routines names usually refer to the operation (read), the dimension
!!                of the required field (1-2-3D), the dimension of the reading target
!!                (f3D, f4D, from 3 or 4 dimensional fields) anbd the kind of variable.
!!                Standard arguments are: <br>
!!
!> @param[in]     subnam: program that executes the call (for error handling, EH) 
!! @param[in]     varnam: variable to inquire
!!                dimnam: dimension to inquire
!! @param[in]     ncID:   ID of the NetCDF file
!! @param[out]    varOUT: output variable
!!                dimOUT: output dimension
!!
!!                See details and examples of netCDF at <br>
!!     https://www.unidata.ucar.edu/software/netcdf/docs-fortran/nc_f77_interface_guide.html
!!
!!------------------------------------------------------------------------------
MODULE ncdfRD_subs
  
   USE class_precision
   !
   IMPLICIT NONE
   !
   !  I/O files
   !  -------------
   CHARACTER (LEN = lc), PUBLIC, SAVE ::  indir      !< @public path of input data
   CHARACTER (LEN = lc), PUBLIC, SAVE :: uncnam      !< @public input NetCDF file of zonal velocity ( U-grid)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: vncnam      !< @public input NetCDF file of meridional velocity ( V-grid)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: wncnam      !< @public input NetCDF file of vertical velocity ( W-grid)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: trcnam      !< @public input NetCDF file of tracer ( T-grid)
   CHARACTER (LEN = lc), PUBLIC, SAVE ::  dmdir      !< @public path of input data
   CHARACTER (LEN = lc), PUBLIC, SAVE :: dmgrid      !< @public input NetCDF file of domain grid
   !  Additional netCDF files
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile1     !< @public additional input NetCDF file (from 1 to 4)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile2
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile3
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile4
   !
   !  Dimension names 
   !  ---------------------
   CHARACTER (LEN = lc), PUBLIC, SAVE :: nm_x        !< @public longitude, I dimension
   CHARACTER (LEN = lc), PUBLIC, SAVE :: nm_y        !< @public latitude,  J dimension
   CHARACTER (LEN = lc), PUBLIC, SAVE :: nm_z1       !< @public depth, z dimension
   CHARACTER (LEN = lc), PUBLIC, SAVE :: nm_z2       !< @public levels, z dimension
   CHARACTER (LEN = lc), PUBLIC, SAVE :: nm_t        !< @public time dimension
   !
   !  Storage for identifiers
   !  -------------
   INTEGER, PUBLIC, SAVE :: uncID                    !< @public ID for U-grid file
   INTEGER, PUBLIC, SAVE :: vncID                    !< @public ID for V-grid file
   INTEGER, PUBLIC, SAVE :: wncID                    !< @public ID for W-grid file
   INTEGER, PUBLIC, SAVE :: trcID                    !< @public ID for T-grid file
   INTEGER, PUBLIC, SAVE :: domID                    !< @public ID for domain file
   !  Additional netCDF identifiers
   INTEGER, PUBLIC, SAVE :: ncID1                    !< @public ID for additional netCDF files (from 1 to 4)
   INTEGER, PUBLIC, SAVE :: ncID2
   INTEGER, PUBLIC, SAVE :: ncID3
   INTEGER, PUBLIC, SAVE :: ncID4
   !
   !
   INTEGER, PARAMETER    :: stderr = 6 
   PRIVATE
   !
   !  ID inquiring
   !  ---------------------
   PUBLIC  :: readFileID
   !
   !  Dimension reading
   !
   PUBLIC  :: readDim
   !
   !  1D subroutines
   !  ---------------------
   PUBLIC  :: read1D_REALvar
   !
   !  2D subroutines
   !  ---------------------
   PUBLIC  :: read2D_REALvar
   PUBLIC  :: read2Df3D_REALvar
   PUBLIC  :: read2Df4D_REALvar
   !
   !  3D subroutines
   !  ---------------------
   PUBLIC  :: read3D_REALvar
   PUBLIC  :: read3Df4D_REALvar
   !
   !  4D subroutines
   !  ---------------------
   PUBLIC  :: read4D_REALvar
   !
   !  Miscellaneous
   !  ---------------------
   PRIVATE :: handle_err
   PUBLIC  :: rdnc_namelist

CONTAINS


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads the file ID
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
   SUBROUTINE readFileID (subnam, iodir, ncnam, ncID)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: iodir    ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: ncnam    ! Parent subroutine name (for EH)
      INTEGER,             INTENT(  OUT) :: ncID     ! NetCDF file ID 

      ! Internal  
      INTEGER :: ncstat

      ! Open netCDF file 
      ncstat = nf_open (TRIM(iodir)//TRIM(ncnam), nf_nowrite, ncID) 
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//subnam//'] '//ncnam)

   END SUBROUTINE readFileID


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads the file ID
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
   SUBROUTINE readDim (subnam, dimnam, ncID, dimOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: dimnam   ! Dimension to be read name
      INTEGER,             INTENT(IN   ) :: ncID     ! NetCDF file ID 
      INTEGER,             INTENT(  OUT) :: dimOUT   ! NetCDF file ID
      INTEGER                            :: dimID    ! Dimension ID   

      ! Internal  
      INTEGER :: ncstat  

      ! Check existence of dimension
      ncstat = nf_inq_dimid (ncID, TRIM(dimnam), dimID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(dimnam))
      ! Read dimension 
      ncstat = nf_inq_dimlen (ncID, dimID , dimOUT )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(dimnam))

   END SUBROUTINE readDim

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 2D variable 
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
   SUBROUTINE read1D_REALvar (subnam, varnam, ncID, varOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )               :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )               :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )               :: ncID     ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:) :: varOUT   ! NetCDF file ID 
      INTEGER                                          :: varID    ! Variable ID

      ! Internal  
      INTEGER :: ncstat

      ! Check existence
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read variable
      ncstat = nf_get_var_real (ncID, varID, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read1D_REALvar


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 2D variable 
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
   SUBROUTINE read2D_REALvar (subnam, varnam, ncID, varOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                 :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                 :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                 :: ncID     ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                            :: varID    ! Variable ID 

      ! Internal  
      INTEGER :: ncstat

      ! Check existence
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read variable
      ncstat = nf_get_var_real (ncID, varID, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read2D_REALvar


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable from a 4D array
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
   SUBROUTINE read2Df3D_REALvar (subnam, varnam, ncID, varOUT, idz)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                 :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                 :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                 :: ncID     ! NetCDF file ID 
      INTEGER,             INTENT(IN   )                 :: idz      ! Index of sub-slicing
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                            :: varID    ! Variable ID   

      ! Internal  
      INTEGER :: ncstat  
      INTEGER :: starts(4) = 1, counts(4) = 1
      INTEGER :: xID      ! x-dimension ID
      INTEGER :: yID      ! y-dimension ID

!     Read eddy-resolving snapshot
      starts(3) = idz

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, TRIM(nm_x), xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, TRIM(nm_y), yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read dimension x
      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read dimension y
      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
     
      ! Read variable
      ncstat = nf_get_vara_real (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read2Df3D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable from a 4D array
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
   SUBROUTINE read2Df4D_REALvar (subnam, varnam, ncID, varOUT, idz, idt)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                 :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                 :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                 :: ncID     ! NetCDF file ID 
      INTEGER,             INTENT(IN   )                 :: idz, idt ! Index of sub-slicing
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                            :: varID    ! Variable ID   

      ! Internal  
      INTEGER :: ncstat  
      INTEGER :: starts(4) = 1, counts(4) = 1
      INTEGER :: xID      ! x-dimension ID
      INTEGER :: yID      ! y-dimension ID

!     Read eddy-resolving snapshot
      starts(3) = idz
      starts(4) = idt

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, TRIM(nm_x), xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, TRIM(nm_y), yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read dimension x
      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Read dimension y
      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
     
      ! Read variable
      ncstat = nf_get_vara_real (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read2Df4D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable 
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
   SUBROUTINE read3D_REALvar (subnam, varnam, ncID, varOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                   :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                   :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                   :: ncID     ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                              :: varID    ! Variable ID    

      ! Internal  
      INTEGER                            :: ncstat  
      INTEGER                            :: starts(3) = 1, counts(3) = 1
      CHARACTER (LEN = lc)               :: xID      ! x-dimension ID
      CHARACTER (LEN = lc)               :: yID      ! y-dimension ID
      CHARACTER (LEN = lc)               :: zID      ! z-dimension ID

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, TRIM(nm_x), xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, TRIM(nm_y), yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Check existence of dimension z
      ncstat = nf_inq_dimid (ncID, TRIM(nm_z1), zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'t', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'u', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'v', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'w', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z2),      zID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))

      ! Read dimension x
      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Read dimension y
      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Read dimension z
      ncstat = nf_inq_dimlen (ncID, zID , counts(3) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))

      ! Read variable
      ncstat = nf_get_vara_real (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read3D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable from a 4D array
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
   SUBROUTINE read3Df4D_REALvar (subnam, varnam, ncID, varOUT, idt)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                   :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                   :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                   :: ncID     ! NetCDF file ID 
      INTEGER,             INTENT(IN   )                   :: idt      ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                              :: varID    ! Variable ID    

      ! Internal  
      INTEGER                            :: ncstat  
      INTEGER                            :: starts(4) = 1, counts(4) = 1
      CHARACTER (LEN = lc)               :: xID      ! x-dimension ID
      CHARACTER (LEN = lc)               :: yID      ! y-dimension ID
      CHARACTER (LEN = lc)               :: zID      ! z-dimension ID

      ! Set time counter
      starts(4) = idt

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, TRIM(nm_x), xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, TRIM(nm_y), yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Check existence of dimension z
      ncstat = nf_inq_dimid (ncID, TRIM(nm_z1), zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'t', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'u', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'v', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'w', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z2),      zID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))

      ! Read dimension x
      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Read dimension y
      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Read dimension z
      ncstat = nf_inq_dimlen (ncID, zID , counts(3) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))

      ! Read variable
      ncstat = nf_get_vara_real (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read3Df4D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Reads a 3D variable from a 4D array
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
   SUBROUTINE read4D_REALvar (subnam, varnam, ncID, varOUT)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   )                     :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   )                     :: varnam   ! Variable to be read name
      INTEGER,             INTENT(IN   )                     :: ncID     ! NetCDF file ID 
      REAL (KIND = wp),    INTENT(  OUT), DIMENSION(:,:,:,:) :: varOUT   ! NetCDF file ID
      INTEGER                                                :: varID    ! Variable ID    

      ! Internal  
      INTEGER                            :: ncstat  
      INTEGER                            :: starts(4) = 1, counts(4) = 1
      CHARACTER (LEN = lc)               :: xID      ! x-dimension ID
      CHARACTER (LEN = lc)               :: yID      ! y-dimension ID
      CHARACTER (LEN = lc)               :: zID      ! z-dimension ID
      CHARACTER (LEN = lc)               :: tID      ! t-dimension ID

      ! Check existence of variable
      ncstat = nf_inq_varid (ncID, TRIM(varnam), varID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))
      ! Check existence of dimension x
      ncstat = nf_inq_dimid (ncID, TRIM(nm_x), xID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Check existence of dimension y
      ncstat = nf_inq_dimid (ncID, TRIM(nm_y), yID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Check existence of dimension z
      ncstat = nf_inq_dimid (ncID, TRIM(nm_z1), zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'t', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'u', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'v', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z1)//'w', zID)
      IF ( ncstat .ne. NF_NOERR ) ncstat = nf_inq_dimid (ncID, TRIM(nm_z2),      zID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))

      ! Check existence of dimension t
      ncstat = nf_inq_dimid (ncID, TRIM(nm_t), tID)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_t))
      ! Read dimension x
      ncstat = nf_inq_dimlen (ncID, xID , counts(1) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_x))
      ! Read dimension y
      ncstat = nf_inq_dimlen (ncID, yID , counts(2) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_y))
      ! Read dimension z
      ncstat = nf_inq_dimlen (ncID, zID , counts(3) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_z1))
      ! Read dimension t
      ncstat = nf_inq_dimlen (ncID, tID , counts(4) )
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(nm_t))

      ! Read variable
      ncstat = nf_get_vara_real (ncID, varID, starts, counts, varOUT)
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//TRIM(subnam)//'] '// TRIM(varnam))

   END SUBROUTINE read4D_REALvar

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
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
   SUBROUTINE rdnc_namelist (nmlIN)

!     Subroutine arguments
      CHARACTER (LEN = *   ), INTENT(IN   ), OPTIONAL :: nmlIN
      CHARACTER (LEN = 1000)                          :: line
      INTEGER                                         :: file_unit
      INTEGER                                         :: iostat

      ! ------------------------------------------------------------------------
      !   Reading namelist
      ! ------------------------------------------------------------------------
      NAMELIST /rdnc_param/ indir, trcnam, uncnam, vncnam, wncnam, dmdir, dmgrid
      NAMELIST  /dim_param/ nm_x, nm_y, nm_z1, nm_z2, nm_t


      INQUIRE (FILE = nmlIN, IOSTAT = iostat)
      IF (iostat /= 0) THEN
         WRITE (stderr, '(3a)') 'Error: file "', TRIM(nmlIN), '" not found!'
      END IF
      OPEN (ACTION='read', FILE = nmlIN, IOSTAT = iostat, NEWUNIT = file_unit)

      READ (NML = rdnc_param, IOSTAT = iostat, UNIT = file_unit)
      READ (NML =  dim_param, IOSTAT = iostat, UNIT = file_unit)

      print *
      print *, '  Input directory:'
      print *, TRIM(indir)
      print * 
      print *, '      Input files:'
      print *, 'T-grid:  ', TRIM(trcnam)
      print *, 'U-grid:  ', TRIM(uncnam)
      print *, 'V-grid:  ', TRIM(vncnam)
      print *, 'W-grid:  ', TRIM(wncnam)
      print * 
      print *, ' Domain directory:'
      print *, TRIM(dmdir)
      print * 
      print *, '  Domain filename:'
      print *, TRIM(dmgrid)

      IF (iostat /= 0) THEN
         WRITE (stderr, '(2a)'      ) 'Error reading file :"', TRIM(nmlIN)
         WRITE (stderr, '(a, i0)'   ) 'iostat was:"', iostat
         BACKSPACE(file_unit)
         READ  (file_unit, FMT='(A)')  line
         WRITE (stderr, '(A)'       ) 'Invalid line : '//TRIM(line)
      END IF
      CLOSE (file_unit)  

   END SUBROUTINE rdnc_namelist

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
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
   SUBROUTINE handle_err (ncstat, fromst)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Subroutine arguments
      INTEGER,             INTENT(IN   )           :: ncstat
      CHARACTER (LEN = *), INTENT(IN   ), OPTIONAL :: fromst

!     fromst is an optional string indicating where the CALL came
!     from that caused the netCDF problem (e.g. subroutine name).

!     Routine which interprets errors from netCDF output functions,
!     prints them to standard output and then kills the whole run.
      IF ( ncstat  .ne.  NF_NOERR ) THEN
         IF ( PRESENT(fromst) ) THEN
            PRINT *, TRIM(fromst)//'  '//TRIM( nf_strerror(ncstat) )
         ELSE
            PRINT *, TRIM( nf_strerror(ncstat) )
         ENDIF
         STOP 'netCDF:: STOPPED'
      ENDIF

   END SUBROUTINE handle_err

END MODULE ncdfRD_subs


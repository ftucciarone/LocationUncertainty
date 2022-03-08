!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: ncdfWR_subs
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
!> @brief         Input-output routines using netCDF files        
!! 
!! 
!! @par           Procedure specifics      
!> @details       This module is used to read and save netCDF files for performing
!!                the program main.F. See details and examples of netCDF at <br>
!!     https://www.unidata.ucar.edu/software/netcdf/docs-fortran/nc_f77_interface_guide.html
!!
!!------------------------------------------------------------------------------
MODULE ncdfWR_subs
  
   USE class_precision
   !
   IMPLICIT NONE
   !
   !  I/O files
   !  -------------
   CHARACTER (LEN = lc), PUBLIC, SAVE :: outdir      !< @public path of output data
   CHARACTER (LEN = lc), PUBLIC, SAVE :: outnam      !< @public output NetCDF file
   !  Additional netCDF files
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile1     !< @public additional input NetCDF file (from 1 to 4)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile2
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile3
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile4
   !
   !  Storage for identifiers
   !  -------------
   INTEGER, PUBLIC, SAVE :: outID                    !< @public ID for W-grid file
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
   !  Miscellaneous
   !  ---------------------
   PRIVATE :: handle_err, save_init
   PUBLIC  :: wrnc_namelist, wrnc_init, save_out

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
   SUBROUTINE save_init

!     Modules   
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE task_data
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'save_init'
      integer xpdim,ypdim,tid,xpid,ypid,zid,dims(4),starts(2),counts(2)

      integer ncstat

      ! Variables IDs
      INTEGER                            ::  tvID,  zvID, yvID, xvID   
      INTEGER                            :: latID, lonID

      dims = 0

      ! Define dimensions
      ncstat = nf_def_dim (outID,  nm_t,  nto, dims(4))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (outID, nm_z2,  nlo, dims(3))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (outID,  nm_y, nyR3, dims(2))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (outID,  nm_x, nxR3, dims(1))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define Time with NEMO convention
      ncstat = nf_def_var (outID, cn_vtimec, NF_DOUBLE, 1, dims(4), tvID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, tvID, 'long_name', 9, 'Time axis')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, tvID, 'units', 5, 'years')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define geographical variables 
      ncstat = nf_def_var (outID, cn_vlat2d, NF_FLOAT, 2, dims(1:2), latID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, latID, 'long_name', 8, 'Latitude')
      ncstat = nf_put_att_text (outID, latID, 'units', 13, 'degrees_north')
 
      ncstat = nf_def_var (outID, cn_vlon2d, NF_FLOAT, 2, dims(1:2), lonID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, lonID, 'long_name', 9, 'Longitude')
      ncstat = nf_put_att_text (outID, lonID, 'units', 12, 'degrees_east')

      ! Define variable of interest
      ncstat = nf_def_var (outID, 'uf', NF_FLOAT, 4, dims, ufID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, ufID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, ufID, 'long_name', 29,         &
     &                          'Coarse-grained zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (outID, 'ur', NF_FLOAT, 4, dims, urID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, urID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, urID, 'long_name', 23,         &
     &                          'Residual zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      
      ncstat = nf_def_var (outID, 'vf', NF_FLOAT, 4, dims, vfID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, vfID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, vfID, 'long_name', 34,         &
     &                          'Coarse-grained meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (outID, 'vr', NF_FLOAT, 4, dims, vrID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, vrID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, vrID, 'long_name', 28,         &
     &                          'Residual meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (outID, 'wf', NF_FLOAT, 4, dims, wfID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, wfID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, wfID, 'long_name', 32,         &
     &                          'Coarse-grained vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (outID, 'wr', NF_FLOAT, 4, dims, wrID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, wrID, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (outID, wrID, 'long_name', 26,         &
     &                          'Residual vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

!     Leave definition mode and entering data mode
      ncstat = nf_enddef (outID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      counts(1) = nxR3
      counts(2) = nyR3

!      ncstat = nf_put_vara_real (outID, latID, starts, counts, lat)
!      ncstat = nf_put_vara_real (outID, lonID, starts, counts, lon)

      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)


   END SUBROUTINE save_init

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
   SUBROUTINE save_out (idt)

!     Write out one data at time of index 'idt'

!     Modules      
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE task_data
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      INTEGER, INTENT(IN) :: idt

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'save_out'
      ! Internal  
      INTEGER                        :: ncstat
      INTEGER                        :: starts(4), counts(4)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      starts(4) = idt
      counts(1) = nxR3
      counts(2) = nyR3
      counts(3) = nlo
      counts(4) = 1

      ncstat = nf_put_vara_real (outID, ufID, starts, counts, uf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (outID, urID, starts, counts, ur)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (outID, vfID, starts, counts, vf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (outID, vrID, starts, counts, vr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (outID, wfID, starts, counts, wf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (outID, wrID, starts, counts, wr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE save_out


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
   SUBROUTINE wrnc_init (subnam, ncID)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: subnam   ! Parent subroutine name (for EH)
      INTEGER,             INTENT(  OUT) :: ncID     ! NetCDF file ID 

      ! Internal  
      INTEGER :: ncstat

      ncstat = nf_create (TRIM(outdir)//TRIM(outnam), &
             &            ior(nf_clobber,nf_64bit_offset), ncID)
      ! This is used to create 'big-data'
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//subnam//'] '//outnam)

      CALL save_init

   END SUBROUTINE wrnc_init


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
   SUBROUTINE wrnc_namelist (nmlIN)

!     Subroutine arguments
      CHARACTER (LEN = *   ), INTENT(IN   ), OPTIONAL :: nmlIN
      CHARACTER (LEN = 1000)                          :: line
      INTEGER                                         :: file_unit
      INTEGER                                         :: iostat

      ! ------------------------------------------------------------------------
      !   Reading namelist
      ! ------------------------------------------------------------------------
      NAMELIST /wrnc_param/ outdir, outnam


      INQUIRE (FILE = nmlIN, IOSTAT = iostat)
      IF (iostat /= 0) THEN
         WRITE (stderr, '(3a)') 'Error: file "', TRIM(nmlIN), '" not found!'
      END IF
      OPEN (ACTION='read', FILE = nmlIN, IOSTAT = iostat, NEWUNIT = file_unit)

      READ (NML = wrnc_param, IOSTAT = iostat, UNIT = file_unit)

      print *
      print *, ' Output directory:'
      print *, TRIM(outdir)
      print * 
      print *, '     Output files:'
      print *, TRIM(outnam)


      IF (iostat /= 0) THEN
         WRITE (stderr, '(2a)'      ) 'Error reading file :"', TRIM(nmlIN)
         WRITE (stderr, '(a, i0)'   ) 'iostat was:"', iostat
         BACKSPACE(file_unit)
         READ  (file_unit, FMT='(A)')  line
         WRITE (stderr, '(A)'       ) 'Invalid line : '//TRIM(line)
      END IF
      CLOSE (file_unit)  

   END SUBROUTINE wrnc_namelist


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
   SUBROUTINE handle_err (ncstat, fromst)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Subroutine arguments
      integer, INTENT(IN) :: ncstat
      character (len=*), INTENT(IN), OPTIONAL :: fromst

!     fromst is an optional string indicating where the call came
!     from that caused the netCDF problem (e.g. subroutine name).

!     Routine which interprets errors from netCDF output functions,
!     prints them to standard output and then kills the whole run.
      if ( ncstat .ne. NF_NOERR ) then
        if ( present(fromst) ) then
          print *, trim(fromst)//'  '//trim( nf_strerror(ncstat) )
         else
          print *, trim( nf_strerror(ncstat) )
        endif
        stop 'netCDF:: STOPPED'
      endif

   END SUBROUTINE handle_err

END MODULE ncdfWR_subs


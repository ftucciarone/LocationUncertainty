!,! ------------------------------------------------------------------------------
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
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile1wr     !< @public additional input NetCDF file (from 1 to 4)
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile2wr
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile3wr
   CHARACTER (LEN = lc), PUBLIC, SAVE :: ncfile4wr
   !
   !  Storage for identifiers
   !  -------------
   INTEGER, PUBLIC, SAVE :: outID                    !< @public ID for W-grid file
   !  Additional netCDF identifiers
   INTEGER, PUBLIC, SAVE :: ncID1wr                   !< @public ID for additional netCDF files (from 1 to 4)
   INTEGER, PUBLIC, SAVE :: ncID2wr
   INTEGER, PUBLIC, SAVE :: ncID3wr
   INTEGER, PUBLIC, SAVE :: ncID4wr
   !
   !
   INTEGER, PARAMETER    :: stderr = 6 
   PRIVATE
   !
   !  Miscellaneous
   !  ---------------------
   PRIVATE :: handle_err
   PUBLIC  :: HR_save_init, LR_save_init 
   PUBLIC  :: HR_save_out, LR_save_out
   PUBLIC  :: HR_save_onetime_out, LR_save_onetime_out
   PUBLIC  :: wrnc_namelist, wrnc_init 

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
   SUBROUTINE LR_save_init (closefile, time)

!     Modules   
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE LR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'
      ! I/O arguments
!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'LR_save_init'
      integer dims(4),starts(2),counts(2), tria(4)
      integer time
      integer ncstat
      logical closefile

      dims = 0
    
      ! Define dimensions
      ncstat = nf_def_dim (LRoutID,  nm_x,         nxLR, dims(1))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (LRoutID,  nm_y,         nyLR, dims(2))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (LRoutID, nm_z2,        nloLR, dims(3))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (LRoutID,  nm_t, time, dims(4))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define Time with NEMO convention
      ncstat = nf_def_var (LRoutID, cn_vtimec, NF_DOUBLE, 1, dims(4), LRtvID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, LRtvID, 'long_name',  9, 'Time axis')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, LRtvID,     'units',  5, 'years')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define geographical variables 
      ncstat = nf_def_var (LRoutID, cn_vlat2d, NF_FLOAT, 2, dims(1:2), LRlatID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, LRlatID, 'long_name',  8, 'Latitude')
      ncstat = nf_put_att_text (LRoutID, LRlatID,     'units', 13, 'degrees_north')
 
      ncstat = nf_def_var (LRoutID, cn_vlon2d, NF_FLOAT, 2, dims(1:2), LRlonID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, LRlonID, 'long_name',  9, 'Longitude')
      ncstat = nf_put_att_text (LRoutID, LRlonID,     'units', 12, 'degrees_east')

      ! Define variable of interest: time average
      ncstat = nf_def_var (LRoutID, 'umean', NF_FLOAT, 3, dims(1:3), meanuID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )
      ncstat = nf_put_att_text (LRoutID, meanuID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )
      ncstat = nf_put_att_text (LRoutID, meanuID_LR, 'long_name', 30,         &
      &                          'Time average of zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )

      ncstat = nf_def_var (LRoutID, 'vmean', NF_FLOAT, 3, dims(1:3), meanvID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )
      ncstat = nf_put_att_text (LRoutID, meanvID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )
      ncstat = nf_put_att_text (LRoutID, meanvID_LR, 'long_name', 35,         &
      &                          'Time average of meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )

      ncstat = nf_def_var (LRoutID, 'wmean', NF_FLOAT, 3, dims(1:3), meanwID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )
      ncstat = nf_put_att_text (LRoutID, meanwID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )
      ncstat = nf_put_att_text (LRoutID, meanwID_LR, 'long_name', 33,         &
      &                          'Time average of vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )

      ! Define variable of interest
      ncstat = nf_def_var (LRoutID, 'uf1', NF_FLOAT, 4, dims, uf1ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uf1ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uf1ID_LR, 'long_name', 29,         &
      &                          'First filtered zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'uf2', NF_FLOAT, 4, dims, uf2ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uf2ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uf2ID_LR, 'long_name', 30,         &
      &                          'Second filtered zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'uband', NF_FLOAT, 4, dims, ubID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, ubID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, ubID_LR, 'long_name', 35,         &
      &                          'Bandpassed (f1 - f2) zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'vf1', NF_FLOAT, 4, dims, vf1ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vf1ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vf1ID_LR, 'long_name', 34,         &
      &                          'First filtered meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'vf2', NF_FLOAT, 4, dims, vf2ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vf2ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vf2ID_LR, 'long_name', 35,         &
      &                          'Second filtered meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'vband', NF_FLOAT, 4, dims, vbID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vbID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vbID_LR, 'long_name', 40,         &
      &                          'Bandpassed (f1 - f2) meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'wf1', NF_FLOAT, 4, dims, wf1ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wf1ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wf1ID_LR, 'long_name', 32,         &
      &                          'First filtered vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_def_var (LRoutID, 'wf2', NF_FLOAT, 4, dims, wf2ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wf2ID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wf2ID_LR, 'long_name', 33,         &
      &                          'Second filtered vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'wband', NF_FLOAT, 4, dims, wbID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wbID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wbID_LR, 'long_name', 38,         &
      &                          'Bandpassed (f1 - f2) vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
         
      ! Define variable of interest
      ncstat = nf_def_var (LRoutID, 'uds', NF_FLOAT, 4, dims, uID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uID_LR, 'long_name', 26,         &
      &                          'Downsampled zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'vds', NF_FLOAT, 4, dims, vID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vID_LR, 'long_name', 31,         &
      &                          'Downsampled meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'wds', NF_FLOAT, 4, dims, wID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, wID_LR, 'long_name', 29,         &
      &                          'Downsampled vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define variable of interest
      ncstat = nf_def_var (LRoutID, 'sshf1', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), sshf1ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshf1ID_LR, 'units', 1, 'm')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshf1ID_LR, 'long_name', 33,         &
      &                          'First filtered sea surface height')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'sshf2', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), sshf2ID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshf2ID_LR, 'units', 1, 'm')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshf2ID_LR, 'long_name', 34,         &
      &                          'Second filtered sea surface height')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'sshband', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), sshbID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshbID_LR, 'units', 1, 'm')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshbID_LR, 'long_name', 39,         &
      &                          'Bandpassed (f1 - f2) sea surface height')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'sshds', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), sshID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshID_LR, 'units', 1, 'm')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, sshID_LR, 'long_name', 30,         &
      &                          'Downsampled sea surface height')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ! Define variable of interest
      ncstat = nf_def_var (LRoutID, 'u_b', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), u_btID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, uID_LR, 'long_name', 37,         &
      &                          'Downsampled barotropic zonal velocity')

      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (LRoutID, 'v_b', NF_FLOAT, 3, (/ dims(1:2), dims(4) /), v_btID_LR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vID_LR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (LRoutID, vID_LR, 'long_name', 42,         &
      &                          'Downsampled barotropic meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
!     Leave definition mode and entering data mode
      if ( closefile ) then
         ncstat = nf_enddef (LRoutID)
         if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      end if

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      counts(1) = nxLR
      counts(2) = nyLR

!      ncstat = nf_put_vara_real (outID, lonID, starts, counts, lon)

      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE LR_save_init




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
   SUBROUTINE HR_save_init (closefile)

!     Modules   
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE HR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'
      ! I/O arguments
      logical closefile      
!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'HR_save_init'
      integer dims(4),starts(2),counts(2)

      integer ncstat


      dims = 0


      ! Define dimensions
      ncstat = nf_def_dim (HRoutID,  nm_t, ntoHR, dims(4))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(nm_t) )
      ncstat = nf_def_dim (HRoutID, nm_z2, nloHR, dims(3))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(nm_z2) )
      ncstat = nf_def_dim (HRoutID,  nm_y,  nyHR, dims(2))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(nm_y) )
      ncstat = nf_def_dim (HRoutID,  nm_x,  nxHR, dims(1))
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(nm_x) )

      ! Define Time with NEMO convention
      ncstat = nf_def_var (HRoutID, cn_vtimec, NF_DOUBLE, 1, dims(4), HRtvID)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(cn_vtimec) )
      ncstat = nf_put_att_text (HRoutID, HRtvID, 'long_name', 9, 'Time axis')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(cn_vtimec) )
      ncstat = nf_put_att_text (HRoutID, HRtvID, 'units', 5, 'years')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM(cn_vtimec) )

      ! Define variable of interest: time average
      ncstat = nf_def_var (HRoutID, 'umean', NF_FLOAT, 3, dims(1:3), ubID_HR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )
      ncstat = nf_put_att_text (HRoutID, ubID_HR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )
      ncstat = nf_put_att_text (HRoutID, ubID_HR, 'long_name', 30,         &
      &                          'Time average of zonal velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('umean') )

      ncstat = nf_def_var (HRoutID, 'vmean', NF_FLOAT, 3, dims(1:3), vbID_HR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )
      ncstat = nf_put_att_text (HRoutID, vbID_HR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )
      ncstat = nf_put_att_text (HRoutID, vbID_HR, 'long_name', 35,         &
      &                          'Time average of meridional velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vmean') )

      ncstat = nf_def_var (HRoutID, 'wmean', NF_FLOAT, 3, dims(1:3), wbID_HR)
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )
      ncstat = nf_put_att_text (HRoutID, wbID_HR, 'units', 3, 'm/s')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )
      ncstat = nf_put_att_text (HRoutID, wbID_HR, 'long_name', 33,         &
      &                          'Time average of vertical velocity')
      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wmean') )

!      ! Define variable of interest: fluctuation
!      ncstat = nf_def_var (HRoutID, 'uflct', NF_FLOAT, 4, dims, uf1ID_HR)
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('uflct') )
!      ncstat = nf_put_att_text (HRoutID, uf1ID_HR, 'units', 3, 'm/s')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('uflct') )
!      ncstat = nf_put_att_text (HRoutID, uf1ID_HR, 'long_name', 32,        &
!      &                          'Time fluctuations zonal velocity')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('uflct') )
!
!      ncstat = nf_def_var (HRoutID, 'vflct', NF_FLOAT, 4, dims, vf1ID_HR)
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vflct') )
!      ncstat = nf_put_att_text (HRoutID, vf1ID_HR, 'units', 3, 'm/s')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vflct') )
!      ncstat = nf_put_att_text (HRoutID, vf1ID_HR, 'long_name', 37,         &
!      &                          'Time fluctuations meridional velocity')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('vflct') )
!
!      ncstat = nf_def_var (HRoutID, 'wflct', NF_FLOAT, 4, dims, wf1ID_HR)
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wflct') )
!      ncstat = nf_put_att_text (HRoutID, wf1ID_HR, 'units', 3, 'm/s')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wflct') )
!      ncstat = nf_put_att_text (HRoutID, wf1ID_HR, 'long_name', 35,         &
!      &                          'Time fluctuations vertical velocity')
!      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam//'::'//TRIM('wflct') )
!


!     Leave definition mode and entering data mode
      if ( closefile ) then
         ncstat = nf_enddef (HRoutID)
         if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)
      end if

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      counts(1) = nxHR
      counts(2) = nyHR

!      ncstat = nf_put_vara_real (outID, lonID, starts, counts, lon)

      if ( ncstat .ne. NF_NOERR ) call handle_err (ncstat, subnam)


   END SUBROUTINE HR_save_init







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
   SUBROUTINE LR_save_out (idt)

!     Write out one data at time of index 'idt'

!     Modules      
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE LR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      INTEGER, INTENT(IN) :: idt

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'LR_save_out'
      ! Internal  
      INTEGER                        :: ncstat
      INTEGER                        :: starts(4), counts(4)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      starts(4) = idt
      counts(1) = nxLR
      counts(2) = nyLR
      counts(3) = nloLR
      counts(4) = 1

      ncstat = nf_put_vara_real (LRoutID, uf1ID_LR, starts, counts, uf1_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, uf2ID_LR, starts, counts, uf2_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, ubID_LR, starts, counts, uband_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_put_vara_real (LRoutID, vf1ID_LR, starts, counts, vf1_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, vf2ID_LR, starts, counts, vf2_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, vbID_LR, starts, counts, vband_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_put_vara_real (LRoutID, wf1ID_LR, starts, counts, wf1_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, wf2ID_LR, starts, counts, wf2_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, wbID_LR, starts, counts, wband_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_put_vara_real (LRoutID, v_btID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/), vLR_b)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, u_btID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/), uLR_b)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_put_vara_real (LRoutID, sshf1ID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/),  sshf1_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, sshf2ID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/),  sshf2_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, sshbID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/), sshband_LR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      ncstat = nf_put_vara_real (LRoutID, uID_LR, starts, counts, uLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, vID_LR, starts, counts, vLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, wID_LR, starts, counts, wLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (LRoutID, sshID_LR, (/starts(1:2), starts(4)/), (/counts(1:2), counts(4)/),  sshLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE LR_save_out


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
   SUBROUTINE HR_save_out (idt)

!     Write out one data at time of index 'idt'

!     Modules      
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE HR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      INTEGER, INTENT(IN) :: idt

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'HR_save_out'
      ! Internal  
      INTEGER                        :: ncstat
      INTEGER                        :: starts(4), counts(4)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      starts(4) = idt
      counts(1) = nxHR
      counts(2) = nyHR
      counts(3) = nloHR
      counts(4) = 1

      ncstat = nf_put_vara_real (HRoutID, uf1ID_HR, starts, counts, uf1_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, uf2ID_HR, starts, counts, uf2_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, ubID_HR, starts, counts, uband_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (HRoutID, vf1ID_HR, starts, counts, vf1_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, vf2ID_HR, starts, counts, vf2_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, vbID_HR, starts, counts, vband_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (HRoutID, wf1ID_HR, starts, counts, wf1_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, wf2ID_HR, starts, counts, wf2_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (HRoutID, wbID_HR, starts, counts, wband_HR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE HR_save_out

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
   SUBROUTINE HR_save_onetime_out (   )

!     Write out one data at time of index 'idt'

!     Modules      
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE HR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'HR_save_onetime_out'
      ! Internal  
      INTEGER                        :: ncstat
      INTEGER                        :: starts(3), counts(3)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      counts(1) = nxHR
      counts(2) = nyHR
      counts(3) = nloHR

      ncstat = nf_put_vara_real (HRoutID, ubID_HR, starts, counts, meanuHR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (HRoutID, vbID_HR, starts, counts, meanvHR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (HRoutID, wbID_HR, starts, counts, meanwHR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE HR_save_onetime_out

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
   SUBROUTINE LR_save_onetime_out (   )

!     Write out one data at time of index 'idt'

!     Modules      
      USE ncdfRD_subs, ONLY: nm_t, nm_z2, nm_y, nm_x
      USE LR_fields
      USE NEMO_nmconv

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments

!     Local parameters
      CHARACTER (LEN = *), PARAMETER :: subnam = 'LR_save_onetime_out'
      ! Internal  
      INTEGER                        :: ncstat
      INTEGER                        :: starts(3), counts(3)

!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      counts(1) = nxLR
      counts(2) = nyLR
      counts(3) = nloLR

      ncstat = nf_put_vara_real (LRoutID, meanuID_LR, starts, counts, meanuLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (LRoutID, meanvID_LR, starts, counts, meanvLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (LRoutID, meanwID_LR, starts, counts, meanwLR)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE LR_save_onetime_out

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
   SUBROUTINE wrnc_init (subnam, initag, ncID)

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      ! I/O arguments
      CHARACTER (LEN = *), INTENT(IN   ) :: subnam   ! Parent subroutine name (for EH)
      CHARACTER (LEN = *), INTENT(IN   ) :: initag   ! Parent subroutine name (for EH)      
      INTEGER,             INTENT(  OUT) :: ncID     ! NetCDF file ID 

      ! Internal  
      INTEGER :: ncstat, omode

      ncstat = nf_create (TRIM(outdir)//TRIM(initag)//TRIM(outnam), &
             &            nf_64bit_offset, ncID)

     ! omode = nf_set_fill(ncID, nf_fill)
      ! This is used to create 'big-data'
      IF ( ncstat .ne. NF_NOERR ) CALL handle_err (ncstat, '['//subnam//'] '//outnam)

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
   SUBROUTINE wrnc_namelist (nmlIN, suffix)

!     Subroutine arguments
      CHARACTER (LEN = *   ), INTENT(IN   )           :: nmlIN
      CHARACTER (LEN = *   ), INTENT(IN   )           :: suffix
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
      outnam = TRIM(outnam)//TRIM(suffix)//".nc"

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





!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: ionc_subs
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
MODULE ionc_subs
  
   USE class_precision
   USE state 
  !
   IMPLICIT NONE
   !
   PRIVATE

!  Subroutines      
   PUBLIC :: read_ini, save_ini, read_out, save_out, handle_err, read_out_sl, read_out_w

!  Storage for identifiers
   INTEGER, PUBLIC,  SAVE :: trcid
   INTEGER, PUBLIC,  SAVE :: uncid
   INTEGER, PUBLIC,  SAVE :: vncid
   INTEGER, PUBLIC,  SAVE :: wncid
   INTEGER, PUBLIC,  SAVE :: oncid   ! ID for files

   INTEGER, PUBLIC,  SAVE :: domID   ! nc ID for domain configuration file


   INTEGER, PRIVATE, SAVE ::    uid,    vid,  wid
   INTEGER, PRIVATE, SAVE :: utauid, vtauid

   INTEGER, PRIVATE, SAVE :: uwndfid, vwndfid
   INTEGER, PRIVATE, SAVE :: uwndrid, vwndrid
   INTEGER, PRIVATE, SAVE ::    ufid,    vfid, wfid
   INTEGER, PRIVATE, SAVE ::    urid,    vrid, wrid  ! ID for variables
   
   INTEGER, PRIVATE, SAVE ::  e1v_ID,  e3v_ID        ! nc IDs for domain variables
   INTEGER, PRIVATE, SAVE ::  e2u_ID,  e3u_ID        ! nc IDs for domain variables
   INTEGER, PRIVATE, SAVE ::  e1t_ID,  e2t_ID        ! nc IDs for domain variables
   INTEGER, PRIVATE, SAVE ::           e3t_ID        ! nc IDs for domain variables
   INTEGER, PRIVATE, SAVE ::  lat_ID,  lon_ID        ! nc IDs for domain variables






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
   SUBROUTINE read_ini
  
      USE param, ONLY : nxpo ,nypo ,nlo, nto    
!!    USE state, ONLY : tyrs,xpo,ypo,zo


      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Local parameters
      CHARACTER (len=*), PARAMETER :: subnam = 'read_ini'
      CHARACTER( len = lc )        :: varnam
      INTEGER ncstat,varid


      !  Velocities
      !  -------------------------------------------------------------
      varnam = 'vozocrtx'
      ncstat = nf_inq_varid (uncid, TRIM(varnam), uid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'vomecrty'
      ncstat = nf_inq_varid (vncid, TRIM(varnam), vid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'vovecrtz'
      ncstat = nf_inq_varid (wncid, TRIM(varnam), wid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)


      !  Wind stress
      !  -------------------------------------------------------------
      varnam = 'sozotaux'
      ncstat = nf_inq_varid (uncid, TRIM(varnam), utauid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'sometauy'
      ncstat = nf_inq_varid (vncid, TRIM(varnam), vtauid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)


      !  Scale factors
      !  -------------------------------------------------------------
      varnam = 'e1v'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e1v_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e1v_ID, e1v)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e3v_0'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e3v_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e3v_ID, e3v)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e2u'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e2u_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e2u_ID, e2u)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e3u_0'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e3u_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e3u_ID, e3u)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e1t'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e1t_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e1t_ID, e1t)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e2t'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e2t_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e2t_ID, e2t)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'e3t_0'
      ncstat = nf_inq_varid (domID, TRIM(varnam), e3t_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, e3t_ID, e3t)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)


      !  Geographical coordinats
      !  -------------------------------------------------------------
      varnam = 'nav_lat'
      ncstat = nf_inq_varid (domID, TRIM(varnam), lat_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, lat_ID, lat)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)

      varnam = 'nav_lon'
      ncstat = nf_inq_varid (domID, TRIM(varnam), lon_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)
      ncstat = nf_get_var_real (domID, lon_ID, lon)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, '['//subnam//'] '//varnam)



!!TODO: [ Optional ] 
!!      Read time axis and domains of NEMO for latter
!!      Here is my version:

!     Read time axis (yrs)
!!      ncstat = nf_inq_varid (incid, 'time', varid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_get_var (incid, varid, tyrs)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!     Read p grid x-axis (km) 
!!      ncstat = nf_inq_varid (incid, 'xp', varid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_get_var (incid, varid, xpo)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!     Read p-grid y-axis (km)     
!!      ncstat = nf_inq_varid (incid, 'yp', varid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_get_var (incid, varid, ypo)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
     
!     Read z-axis (km)     
!!      ncstat = nf_inq_varid (incid, 'z', varid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_get_var (incid, varid, zo)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

   END SUBROUTINE read_ini

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
   SUBROUTINE read_out (u, v, idt, idz)

!     Read input data at one time of index 'idt' 

!     Modules
      USE param, ONLY : nxpo,nypo,nlo    

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpo,nypo) :: u,v
      INTEGER, INTENT(IN) :: idt,idz

!     Local parameters
      CHARACTER (len=*), PARAMETER :: subnam = 'read_out'
      INTEGER ncstat,starts(4),counts(4)

!     Read eddy-resolving snapshot
      starts(1) = 1
      starts(2) = 1
      starts(3) = idz
      starts(4) = idt
      counts(1) = nxpo
      counts(2) = nypo
      counts(3) = 1
      counts(4) = 1

      !  Velocities
      !  -------------------------------------------------------------
      ncstat = nf_get_vara_real (uncid, uid, starts, counts, u)
      IF ( ncstat.ne.NF_NOERR ) CALL handle_err (ncstat, subnam)
      ncstat = nf_get_vara_real (vncid, vid, starts, counts, v)
      IF ( ncstat.ne.NF_NOERR ) CALL handle_err (ncstat, subnam)

   END SUBROUTINE read_out

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
   SUBROUTINE read_out_w (w, idt, idz)

!     Read input data at one time of index 'idt' 

!     Modules
      USE param, ONLY : nxpo,nypo,nlo    

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpo,nypo) :: w
      INTEGER, INTENT(IN) :: idt,idz

!     Local parameters
      CHARACTER (len=*), PARAMETER :: subnam = 'read_out_w'
      INTEGER ncstat,starts(4),counts(4)

!     Read eddy-resolving snapshot
      starts(1) = 1
      starts(2) = 1
      starts(3) = idz
      starts(4) = idt
      counts(1) = nxpo
      counts(2) = nypo
      counts(3) = 1
      counts(4) = 1

      !  Velocities
      !  -------------------------------------------------------------
      ncstat = nf_get_vara_real (wncid, wid, starts, counts, w)
      IF ( ncstat.ne.NF_NOERR ) CALL handle_err (ncstat, subnam)
   END SUBROUTINE read_out_w

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
   SUBROUTINE read_out_sl (utau, vtau, idt)

!     Read input data at one time of index 'idt' 

!     Modules
      USE param, ONLY : nxpo,nypo,nlo    

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpo,nypo) :: utau,vtau
      INTEGER, INTENT(IN) :: idt

!     Local parameters
      CHARACTER (len=*), PARAMETER :: subnam = 'read_out'
      INTEGER ncstat,starts(3),counts(3)

!     Read eddy-resolving snapshot
      starts(1) = 1
      starts(2) = 1
      starts(3) = idt
      counts(1) = nxpo
      counts(2) = nypo
      counts(3) = 1

      !  Wind stress
      !  -------------------------------------------------------------
      ncstat = nf_get_vara_real (uncid, utauid, starts, counts, utau)
      IF ( ncstat.ne.NF_NOERR ) CALL handle_err (ncstat, subnam)
      ncstat = nf_get_vara_real (vncid, vtauid, starts, counts, vtau)
      IF ( ncstat.ne.NF_NOERR ) CALL handle_err (ncstat, subnam)

   END SUBROUTINE read_out_sl



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
   SUBROUTINE save_ini

!     Modules      
      USE param, ONLY : nto,nxpoc,nypoc,nxtoc,nytoc,nlo    
!!      USE state, ONLY : xpoc,ypoc,xtoc,ytoc,zo,tyrs 

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     Local parameters
      character (len=*), parameter :: subnam = 'save_ini'
      integer ncstat,tdim,xpdim,ypdim,zdim,tid,xpid,ypid,zid,dims(4),starts(2),counts(2)
    
      INTEGER latid, lonid

!     Define dimensions
      ncstat = nf_def_dim (oncid, 'time', nto, tdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (oncid, 'xp', nxpoc, xpdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (oncid, 'yp', nypoc, ypdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_def_dim (oncid, 'z', nlo, zdim)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

!!TODO: [ Optional ] 
!!      Create output time axis and domains of NEMO for latter
!!      Here is my version:

!     Define 1D variables with attributes
!!      ncstat = nf_def_var (oncid, 'time', NF_DOUBLE, 1, tdim, tid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, tid, 'long_name', 9, 'Time axis')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, tid, 'units', 5, 'years')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      
!!      ncstat = nf_def_var (oncid, 'xp', NF_DOUBLE, 1, xpdim, xpid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, xpid, 'units', 2, 'km')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, xpid, 'long_name', 15, 
!!     &                          'X axis (p-grid)')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      
!!      ncstat = nf_def_var (oncid, 'yp', NF_DOUBLE, 1, ypdim, ypid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, ypid, 'units', 2, 'km')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, ypid, 'long_name', 15, 
!!     &                          'Y axis (p-grid)')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      
!!     ncstat = nf_def_var (oncid, 'z', NF_DOUBLE, 1, zdim, zid)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, zid, 'units', 2, 'km')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_att_text (oncid, zid, 'long_name', 20,
!!     &                          'Mid-layer depth axis')
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
 
!     Define 4D variables with attributes
      dims(1) = xpdim 
      dims(2) = ypdim
      dims(3) = zdim
      dims(4) = tdim
   

!!---------- Define variables 
 
      ncstat = nf_def_var (oncid, 'nav_lat', NF_FLOAT, 2, dims(1:2), lat_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, latid, 'long_name', 8, 'Latitude')
      ncstat = nf_put_att_text (oncid, latid, 'units', 13, 'degrees_north')



 
      ncstat = nf_def_var (oncid, 'nav_lon', NF_FLOAT, 2, dims(1:2), lon_ID)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, lonid, 'long_name', 9, 'Longitude')
      ncstat = nf_put_att_text (oncid, lonid, 'units', 12, 'degrees_east')



!!----------- 
 
      ncstat = nf_def_var (oncid, 'uf', NF_FLOAT, 4, dims, ufid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, ufid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, ufid, 'long_name', 29,         &
     &                          'Coarse-grained zonal velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'ur', NF_FLOAT, 4, dims, urid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, urid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, urid, 'long_name', 23,         &
     &                          'Residual zonal velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      
      ncstat = nf_def_var (oncid, 'vf', NF_FLOAT, 4, dims, vfid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vfid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vfid, 'long_name', 34,         &
     &                          'Coarse-grained meridional velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'vr', NF_FLOAT, 4, dims, vrid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vrid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vrid, 'long_name', 28,         &
     &                          'Residual meridional velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'wf', NF_FLOAT, 4, dims, wfid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, wfid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, wfid, 'long_name', 32,         &
     &                          'Coarse-grained vertical velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'wr', NF_FLOAT, 4, dims, wrid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, wrid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, wrid, 'long_name', 26,         &
     &                          'Residual vertical velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

 
      ncstat = nf_def_var (oncid, 'uwndf', NF_FLOAT, 3, [dims(1:2),dims(4)], uwndfid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, uwndfid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, uwndfid, 'long_name', 34,         &
     &                          'Coarse-grained zonal wind velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'uwndr', NF_FLOAT, 3, [dims(1:2),dims(4)], uwndrid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, uwndrid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, uwndrid, 'long_name', 28,         &
     &                          'Residual zonal wind velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      
      ncstat = nf_def_var (oncid, 'vwndf', NF_FLOAT, 3, [dims(1:2),dims(4)], vwndfid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vwndfid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vwndfid, 'long_name', 39,         &
     &                          'Coarse-grained meridional wind velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_def_var (oncid, 'vwndr', NF_FLOAT, 3, [dims(1:2),dims(4)], vwndrid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vwndrid, 'units', 3, 'm/s')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_att_text (oncid, vwndrid, 'long_name', 33,         &
     &                          'Residual meridional wind velocity')
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)



!     Leave definition mode and entering data mode
      ncstat = nf_enddef (oncid)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      counts(1) = nxpoc
      counts(2) = nypoc

      ncstat = nf_put_vara_real (oncid, lat_ID, starts, counts, lat)
      ncstat = nf_put_vara_real (oncid, lon_ID, starts, counts, lon)


      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)








!!TODO: [ Optional ] 
!!      Write out time axis and domains of NEMO for latter
!!      Here is my version:

!     Write 1D data to variables
!!      ncstat = nf_put_vara_double (oncid, tid, 1, nto, tyrs)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_vara_double (oncid, xpid, 1, nxpoc, xpoc)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!      ncstat = nf_put_vara_double (oncid, ypid, 1, nypoc, ypoc)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
!!     ncstat = nf_put_vara_double (oncid, zid, 1, nlo, zo)
!!      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
   END SUBROUTINE save_ini



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
      USE param, ONLY : nxpoc,nypoc,nlo    
      USE state 

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

!     I/O arguments
      INTEGER, INTENT(IN) :: idt

!     Local parameters
      CHARACTER (len=*), PARAMETER :: subnam = 'save_out'
      INTEGER ncstat,starts(4),counts(4)
!     Store 4D variables
      starts(1) = 1
      starts(2) = 1
      starts(3) = 1
      starts(4) = idt
      counts(1) = nxpoc
      counts(2) = nypoc
      counts(3) = nlo
      counts(4) = 1
      ncstat = nf_put_vara_real (oncid, ufid, starts, counts, uf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (oncid, urid, starts, counts, ur)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (oncid, vfid, starts, counts, vf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (oncid, vrid, starts, counts, vr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (oncid, wfid, starts, counts, wf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (oncid, wrid, starts, counts, wr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)


      !  Wind velocity
      !  -------------------------------------------------------------
      ncstat = nf_put_vara_real (oncid, uwndfid, [1,1,idt], [nxpoc,nypoc,1], uwindf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (oncid, uwndrid, [1,1,idt], [nxpoc,nypoc,1], uwindr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      ncstat = nf_put_vara_real (oncid, vwndfid, [1,1,idt], [nxpoc,nypoc,1], vwindf)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)
      ncstat = nf_put_vara_real (oncid, vwndrid, [1,1,idt], [nxpoc,nypoc,1], vwindr)
      if ( ncstat.ne.NF_NOERR ) call handle_err (ncstat, subnam)

      END SUBROUTINE save_out

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
      if ( ncstat.ne.NF_NOERR ) then
        if ( present(fromst) ) then
          print *, trim(fromst)//'  '//trim( nf_strerror(ncstat) )
         else
          print *, trim( nf_strerror(ncstat) )
        endif
        stop 'netCDF:: STOPPED'
      endif

   END SUBROUTINE handle_err

END MODULE ionc_subs


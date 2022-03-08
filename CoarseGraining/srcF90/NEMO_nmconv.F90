!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: NEMO_nmconv
!!   Sets global variables necessary to work with the NetCDF files
!!   
!!   
!
!> @authors
!>                F. TUCCIARONE 
!!
!> @version
!!                2021 -  ( F. TUCCIARONE   )  .F90 adaptation
!! 
!> @brief         
!! 
!! 
!! @par           Procedure specifics      
!> @details       Defines the name of the input files, the names of the output,
!!                the name of the variables and the array dimensions
!!
!!------------------------------------------------------------------------------
MODULE NEMO_nmconv

   USE class_precision

   IMPLICIT NONE

   PUBLIC
   SAVE
   !
   !  Dimension name : cn_. [ 1 letter only ]
   !  -------------
   CHARACTER (LEN = lc) :: cn_x = 'x'                               ! longitude, I dimension
   CHARACTER (LEN = lc) :: cn_y = 'y'                               ! latitude,  J dimension
   CHARACTER (LEN = lc) :: cn_z = 'depth'                           ! depth, z dimension
   CHARACTER (LEN = lc) :: cn_t = 'time_counter'                    ! time dimension
   !
   !  Dimension variable
   !  -------------
   CHARACTER (LEN = lc) :: cn_vlon2d  = 'nav_lon'                   ! longitude
   CHARACTER (LEN = lc) :: cn_vlat2d  = 'nav_lat'                   ! latitude
   CHARACTER (LEN = lc) :: cn_vdeptht = 'deptht'                    ! depth
   CHARACTER (LEN = lc) :: cn_vdepthu = 'depthu'                    ! depth
   CHARACTER (LEN = lc) :: cn_vdepthv = 'depthv'                    ! depth
   CHARACTER (LEN = lc) :: cn_vdepthw = 'depthw'                    ! depth
   CHARACTER (LEN = lc) :: cn_vtimec  = 'time_counter'              ! time 
   CHARACTER (LEN = lc) :: cn_vlon1d  = 'lon'                       ! longitude 1d
   CHARACTER (LEN = lc) :: cn_vlat1d  = 'lat'                       ! latitude  1d
   !
   !  Metrics
   !  -------------
   CHARACTER (LEN = lc) :: cn_ve1t = 'e1t', cn_ve2t = 'e2t'         ! e.t
   CHARACTER (LEN = lc) :: cn_ve1u = 'e1u', cn_ve2u = 'e2u'         ! e.u
   CHARACTER (LEN = lc) :: cn_ve1v = 'e1v', cn_ve2v = 'e2v'         ! e.v
   CHARACTER (LEN = lc) :: cn_ve1f = 'e1f', cn_ve2f = 'e2f'         ! e.v
   CHARACTER (LEN = lc) :: cn_ve3t1d = 'e3t'                        ! e3  (1D). 
   CHARACTER (LEN = lc) :: cn_ve3w1d = 'e3w'                        ! e3. (1D). 
   CHARACTER (LEN = lc) :: cn_ve3t = 'e3t', cn_ve3w = 'e3w'         ! e3. (3D). 
   CHARACTER (LEN = lc) :: cn_ve3u = 'e3u', cn_ve3v = 'e3v'         ! e3.
   !
   !  VVL case
   !  -------------
   CHARACTER (LEN = lc) :: cn_ve3tvvl = 'e3t', cn_ve3wvvl = 'e3w'   ! e3. (3D). 
   CHARACTER (LEN = lc) :: cn_ve3uvvl = 'e3u', cn_ve3vvvl = 'e3v'   ! e3.
   CHARACTER (LEN = lc) :: cn_ve3t0 = 'e3t_0', cn_ve3w0 = 'e3w_0'   ! e3. (3D). (at rest)
   CHARACTER (LEN = lc) :: cn_ve3u0 = 'e3u_0', cn_ve3v0 = 'e3v_0'   ! e3.

   CHARACTER (LEN = lc) :: cn_vff = 'ff'

   CHARACTER (LEN = lc) :: cn_gdept = 'gdept', cn_gdepw = 'gdepw'   ! 1d dep variable
   CHARACTER (LEN = lc) :: cn_hdept = 'hdept', cn_hdepw = 'hdepw'   ! 2d dep variable

   CHARACTER (LEN = lc) :: cn_dept3d = 'gdept_0'                    ! initial dept 3D
   CHARACTER (LEN = lc) :: cn_depu3d = 'depu3d'                     ! Local depth U in broken line extraction
   CHARACTER (LEN = lc) :: cn_depw3d = 'depw3d'                     ! Local depth W in broken line extraction

   CHARACTER (LEN = lc) :: cn_glamt = 'glamt', cn_gphit = 'gphit'   ! glam gphi
   CHARACTER (LEN = lc) :: cn_glamu = 'glamu', cn_gphiu = 'gphiu'   ! glam gphi
   CHARACTER (LEN = lc) :: cn_glamv = 'glamv', cn_gphiv = 'gphiv'   ! glam gphi
   CHARACTER (LEN = lc) :: cn_glamf = 'glamf', cn_gphif = 'gphif'   ! glam gphi
   !
   !  Mask variables
   !  -------------
   CHARACTER (LEN = lc) :: cn_tmask = 'tmask', cn_umask = 'umask'   ! tmask, umask
   CHARACTER (LEN = lc) :: cn_vmask = 'vmask', cn_fmask = 'fmask'   ! vmask, fmask
   CHARACTER (LEN = lc) :: cn_tmaskutil = 'tmaskutil'               ! tmaskutil
   CHARACTER (LEN = lc) :: cn_polymask = 'polymask'                 ! polymask
   CHARACTER (LEN = lc) :: cn_tmaskatl = 'tmaskatl'                 ! atlantic mask in cn_fbasins
   CHARACTER (LEN = lc) :: cn_tmaskpac = 'tmaskpac'                 ! pacific mask in cn_fbasins
   CHARACTER (LEN = lc) :: cn_tmaskind = 'tmaskind'                 ! indian mask in cn_fbasins
   CHARACTER (LEN = lc) :: cn_tmaskant = 'tmaskant'                 ! austral mask in cn_fbasins
   CHARACTER (LEN = lc) :: cn_tmaskmed = 'tmaskmed'                 ! mediterranean mask in cn_fbasins
   !
   !  Generic mesh-mask file names  cn_f...
   !  -------------
   CHARACTER (LEN = lc) :: cn_fzgr = 'mesh_zgr.nc'
   CHARACTER (LEN = lc) :: cn_fe3t = 'mesh_zgr.nc'
   CHARACTER (LEN = lc) :: cn_fe3u = 'mesh_zgr.nc'
   CHARACTER (LEN = lc) :: cn_fe3v = 'mesh_zgr.nc'
   CHARACTER (LEN = lc) :: cn_fe3w = 'mesh_zgr.nc'
   CHARACTER (LEN = lc) :: cn_fhgr = 'mesh_hgr.nc'
   CHARACTER (LEN = lc) :: cn_fmsk = 'mask.nc'
   CHARACTER (LEN = lc) :: cn_fcoo = 'coordinates.nc'
   CHARACTER (LEN = lc) :: cn_fbasins = 'new_maskglo.nc'
   !
   !  Variable name  : cn_v... [ starts with cn_v ]
   !  -------------
   CHARACTER (LEN = lc) :: cn_votemper = 'votemper'                 ! temperature
   CHARACTER (LEN = lc) :: cn_vosaline = 'vosaline'                 ! salinity
   CHARACTER (LEN = lc) :: cn_vozocrtx = 'vozocrtx'                 ! zonal velocity
   CHARACTER (LEN = lc) :: cn_vomecrty = 'vomecrty'                 ! meridional velocity
   CHARACTER (LEN = lc) :: cn_vomeeivv = 'vomeeivv'                 ! meridional Eddy Induced Velocity
   CHARACTER (LEN = lc) :: cn_vovecrtz = 'vovecrtz'                 ! vertical velocity
   CHARACTER (LEN = lc) :: cn_sossheig = 'sossheig'                 ! Sea Surface Height
   CHARACTER (LEN = lc) :: cn_somxl010 = 'somxl010'                 ! Mixed layer depth (density criterion)
   CHARACTER (LEN = lc) :: cn_somxlt02 = 'somxlt02'                 ! Mixed layer depth (temperature criterion)
   CHARACTER (LEN = lc) :: cn_sozotaux = 'sozotaux'                 ! Zonal wind stress

END MODULE NEMO_nmconv    

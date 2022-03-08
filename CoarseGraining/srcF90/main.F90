!! ------------------------------------------------------------------------------
!! 
!!                               Program: main
!!   Template to open NetCDF files in Fortran90 with description of the operations
!!   
!!   
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2021 -  ( F. TUCCIARONE   )  .F90 adaptation
!! 
!> @brief         Filtering routines        
!! 
!! 
!! @par           Procedure specifics      
!> @details       This program is used to perform the coarse-graining procedure of
!!                eddy-resolving snapshots in order to construct the POD hereafter.
!!
!!------------------------------------------------------------------------------
PROGRAM main

!  Modules
   USE class_precision        ! Defines single and double precision
   USE nmlsRD_subs            ! Manages the namelist
   USE ncdfRD_subs            ! Reading netCDF routines
   USE ncdfWR_subs            ! Writing netCDF routines
   USE NEMO_nmconv            ! Defines NEMO conventions for naming

   USE task_data
   USE task_operations

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   !-------------------------------------------------------------------------------------
   !   Local variables
   !-------------------------------------------------------------------------------------
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:) :: gfk
   REAL(kind = wp )                              :: wtot 
   INTEGER                                       :: idt, idz
   INTEGER                                       :: ierr(3)

   !-------------------------------------------------------------------------------------
   !   netCDF variables 
   !-------------------------------------------------------------------------------------
   CHARACTER (len=*), PARAMETER                     :: subnam = 'main'
   INTEGER                                          :: status, iostat, file_unit

   !-------------------------------------------------------------------------------------
   !   openMP Variables
   !-------------------------------------------------------------------------------------
   !$ INTEGER :: nprocs, OMP_GET_NUM_PROCS, nthmax, OMP_GET_MAX_THREADS

   !-------------------------------------------------------------------------------------
   !   Reading domain configuration and dimensions
   !-------------------------------------------------------------------------------------
   CALL rdnc_namelist ('namelist_ioNetCDF.nml')     ! Read address files

   CALL readFileID (subnam, indir, uncnam, uncID)      ! R27 domain file
   CALL readDim (subnam,             nm_x, uncID, nxR27)
   CALL readDim (subnam,             nm_y, uncID, nyR27)
   CALL readDim (subnam, TRIM(nm_z1)//'u', uncID,   nlo)
   CALL readDim (subnam,             nm_t, uncID,   nto)

   CALL readFileID (subnam, dmdir, dmgrid, domID)      !  R3 domain file
   CALL readDim (subnam,             nm_x, domID,  nxR3)
   CALL readDim (subnam,             nm_y, domID,  nyR3)
   CALL readDim (subnam,      TRIM(nm_z2), domID,   nlo)

   print *
   print *, 'R27 Dimensions (nx,ny,nz): ', nxR27, nyR27, nlo
   print *, ' R3 Dimensions (nx,ny,nz): ',  nxR3,  nyR3, nlo
   print *, 'R27        Time axis (nt): ',   nto
 
   CALL task_initR27 ('vel')   
   CALL  task_initR3 ('all')
   CALL   state_init ()
 
   !-------------------------------------------------------------------------------------
   !   Initialize outputs
   !-------------------------------------------------------------------------------------
   ! Create netCDF file
   CALL wrnc_namelist ('namelist_ioNetCDF.nml')
   CALL wrnc_init (subnam, outID)

   !-------------------------------------------------------------------------------------
   !   Derived grid parameters (do not alter)
   !-------------------------------------------------------------------------------------
   nso   = 9          ! Resolutions Ratio 
   nxpo  = nxR27      ! oceanic complete grid with boundaries in x 
   nypo  = nyR27      ! oceanic complete grid with boundaries in y
   nxto  = nxpo - 2   ! oceanic grid without boundaries in x 
   nyto  = nypo - 2   ! oceanic grid without boundaries in x 
   nxtoc = nxto / nso ! coarse inner gridcells in x
   nytoc = nyto / nso ! coarse inner gridcells in y
   nxpoc = nxtoc + 2  ! complete coarse grid points with boundaries in x
   nypoc = nytoc + 2  ! complete coarse grid points with boundaries in y
   co    = 2

   ALLOCATE( gfk(2*co*nso+1,2*co*nso+1) )
  
   print *
   print *,            '     Generate Gaussian filter kernel'

   CALL gaussf (gfk, wtot, 2*co*nso)

   write(*,'(a,i2)')   '         Width of Gaussian filter = ', 2*co*nso
   write(*,'(a,f6.3)') '  Total weight of Gaussian filter = ', wtot


   !-------------------------------------------------------------------------------------
   !   Reading Velocity fields
   !-------------------------------------------------------------------------------------
   ! Open netCDF file 
   CALL readFileID (subnam, indir, trcnam, trcID)  ! T-grid
   CALL readFileID (subnam, indir, uncnam, uncID)  ! U-grid
   CALL readFileID (subnam, indir, vncnam, vncID)  ! V-grid
   CALL readFileID (subnam, indir, wncnam, wncID)  ! W-grid

   CALL    read2D_REALvar (subnam,  cn_ve2u, domID, e2uR3)
   CALL    read2D_REALvar (subnam,  cn_ve1v, domID, e1vR3)
   CALL    read2D_REALvar (subnam,  cn_ve1t, domID, e1tR3)
   CALL    read2D_REALvar (subnam,  cn_ve2t, domID, e2tR3)
   CALL read3Df4D_REALvar (subnam, cn_ve3u0, domID, e3uR3, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3v0, domID, e3vR3, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3t0, domID, e3tR3, 1)

   !-------------------------------------------------------------------------------------
   !   Examine openMP environment
   !-------------------------------------------------------------------------------------
   !$    nprocs = OMP_GET_NUM_PROCS()
   !$    nthmax = OMP_GET_MAX_THREADS()
   !$    write(*,*) ' ' 
   !$    write(*,*) ' OpenMP parallelism is activated'
   !$    write(*,'(a,i5)') '   No. of processors available = ', nprocs
   !$    write(*,'(a,i3)') ' Max. no. of threads available = ', nthmax

   !-------------------------------------------------------------------------------------
   !   Initialize inputs
   !-------------------------------------------------------------------------------------
   !
   !                                ! ================
   DO idt = 1, nto                  !  Time instant
      !                             ! ================
      !
      uf = 0._wp
      vf = 0._wp
      wf = 0._wp
      !
      print *, 'Snapshot No. = ', idt
      !
      !                             ! ================
      DO idz = 1, nlo               ! Horizontal slab
         !                          ! ================
         !
         ! Read one eddy-resolving snapshot at one layer
         CALL read2Df4D_REALvar (subnam, cn_vozocrtx, uncID, uR27(:,:,idz), idz, idt)
         CALL read2Df4D_REALvar (subnam, cn_vomecrty, vncID, vR27(:,:,idz), idz, idt)
         !
         call filt (uf(:,:,idz), ur(:,:,idz), uR27(:,:,idz),                   & 
              &     gfk, wtot, 2*co*nso, nxtoc-1, nytoc) 
         call filt (vf(:,:,idz), vr(:,:,idz), vR27(:,:,idz),                   &
              &     gfk, wtot, 2*co*nso, nxtoc, nytoc-1)
         !
         !                          ! ================
      END DO                        ! End of slab
      !                             ! ================
      !
      ! Derive vertical velocities
      CALL uv2w (wf, uf, vf, e1tR3, e2tR3, e1vR3, e2uR3, e3uR3, e3vR3, e3tR3, nxR3, nyR3)
      CALL uv2w (wr, ur, vr, e1tR3, e2tR3, e1vR3, e2uR3, e3uR3, e3vR3, e3tR3, nxR3, nyR3)
      !
      !
      !
      ! Save output file
      CALL save_out (idt) !! defined in 'ionc_subs.F'
      !
      !                             ! ================
   END DO                           !  End of time
   !                                ! ================
   print *
   print *, '*******************************************************'
   print *, '******************  END of program   ******************'
   print *, '*******************************************************'

END PROGRAM main    



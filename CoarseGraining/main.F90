!! ------------------------------------------------------------------------------
!! 
!!                               Program: main
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
!> @brief         Filtering routines        
!! 
!! 
!! @par           Procedure specifics      
!> @details       This program is used to perform the coarse-graining procedure of
!!                eddy-resolving snapshots in order to construct the POD hereafter.
!!
!!------------------------------------------------------------------------------
PROGRAM main

!    Modules
   USE class_precision
   USE ionml_subs
   USE param    
   USE state
   USE ionc_subs
   USE coarse_subs

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   !-------------------------------------------------------------------------------------
   !   Local variables
   !-------------------------------------------------------------------------------------
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:,:)  ::   uo,vo,wo  ! (m/s)
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:,:)  ::   wtilde    ! (m/s)

   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:)    ::   utau,vtau ! (N/m^2)  
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:)    ::   uwnd,vwnd ! (N/m^2)  

   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:)    ::   gfk
   REAL(kind = wp )                                 ::   wtot 
   INTEGER                                          ::   i,j

   !-------------------------------------------------------------------------------------
   !   netCDF variables 
   !-------------------------------------------------------------------------------------
   CHARACTER (len=*), PARAMETER                    ::   subnam = 'main'
   INTEGER                                         ::  status, iostat, file_unit

   !-------------------------------------------------------------------------------------
   !   openMP Variables
   !-------------------------------------------------------------------------------------
!$ INTEGER :: nprocs, OMP_GET_NUM_PROCS, nthmax, OMP_GET_MAX_THREADS

   !-------------------------------------------------------------------------------------
   !   Reading namelist
   !-------------------------------------------------------------------------------------
   NAMELIST/inout_param/ iodir, trcnam, uncnam, vncnam, wncnam, tggrid, outnam
   NAMELIST/array_param/ nxto, nyto, nlo, nto, nso, co 

   CALL open_inputfile('namelist_coarsegrain.nml', file_unit, iostat)
   READ (NML=inout_param, IOSTAT=iostat, UNIT=file_unit)
   READ (NML=array_param, IOSTAT=iostat, UNIT=file_unit)
   CALL close_inputfile('namelist_coarsegrain.nml', file_unit, iostat)
 
   !  Derived grid parameters (do not alter)
   !  -------------------------------------- 
   nxpo  = nxto + 2   ! oceanic complete grid with boundaries in x 
   nypo  = nyto + 2   ! oceanic complete grid with boundaries in y
   nxtoc = nxto/nso   ! coarse inner gridcells in x
   nytoc = nyto/nso   ! coarse inner gridcells in y
   nxpoc = nxtoc + 2  ! complete coarse grid points in x
   nypoc = nytoc + 2  ! complete coarse grid points in y


   CALL state_init()
   ALLOCATE( uo(nxpo,nypo,nlo),   &
        &    vo(nxpo,nypo,nlo),   &
        &    wo(nxpo,nypo,nlo),   &
        &    wtilde(nxpo,nypo,nlo),   &
        &    utau(nxpo,nypo),     &
        &    vtau(nxpo,nypo),     &
        &    uwnd(nxpo,nypo),     &
        &    vwnd(nxpo,nypo),     &
        &    gfk(2*co*nso+1,2*co*nso+1))

 
   print *, '*******************************************************'
   print *, '****  Coarse-graining of eddy-resolving snapshots   ***'
   print *, '*******************************************************'

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
   ! Open netCDF file for T-grid 
   status = nf_open (TRIM(iodir)//TRIM(trcnam), nf_nowrite, trcid) 
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//trcnam)
   ! Open netCDF file for U-grid 
   status = nf_open (TRIM(iodir)//TRIM(uncnam), nf_nowrite, uncid) 
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//uncnam)
   ! Open netCDF file for V-grid 
   status = nf_open (TRIM(iodir)//TRIM(vncnam), nf_nowrite, vncid)
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//vncnam)
   ! Open netCDF file for W-grid 
   status = nf_open (TRIM(iodir)//TRIM(wncnam), nf_nowrite, wncid) 
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//wncnam)

 
   ! This is to read the domain file
   status = nf_open (TRIM(iodir)//TRIM(tggrid), nf_nowrite, domID) 
   ! This will generate an identifier 'uncid' for input *.nc file
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//tggrid)
   CALL read_ini !! defined in 'ionc_subs.F'
   print *
   print *, 'Input file initialized'

   !-------------------------------------------------------------------------------------
   !   Construct coarse grid
   !-------------------------------------------------------------------------------------
   !   x axis of p grid      
!!      do i=1,nxpoc 
!!        xpoc(i) = xpo(1+(i-1)*nso)
!!      enddo
  !   y axis of p grid        
!!      do j=1,nypoc
!!        ypoc(j) = ypo(1+(j-1)*nso)
!!      enddo
   !-------------------------------------------------------------------------------------
   !   Initialize outputs
   !-------------------------------------------------------------------------------------
   ! Create netCDF file
   status = nf_create (TRIM(iodir)//TRIM(outnam), &
          &            ior(nf_clobber,nf_64bit_offset), oncid)
   ! This is used to create 'big-data'
   IF ( status.ne.NF_NOERR ) CALL handle_err (status, '['//subnam//'] '//outnam)
   CALL save_ini !! defined in 'ionc_subs.F'
   print *, 'Output file created'
   !-------------------------------------------------------------------------------------
   !   Coarse-graining of input snapshots 
   !-------------------------------------------------------------------------------------
  
   print *
   print *, 'Generate Gaussian filter kernel'
   call gaussf (gfk, wtot, 2*co*nso) !! defined in 'coarse_subs.F'
   write(*,'(a,i2)') '  Width of Gaussian filter = ', 2*co*nso
   write(*,'(a,f6.3)') '  Total weight of Gaussian filter = ', wtot

   print *
   print *, 'Coarse-graining procedure'
   print *, '========================='

   !
   !                                ! ================
   DO i = 1, nto                    !  Time instant
      !                             ! ================
      !
      uo = 0._wp
      vo = 0._wp
      !
      print *, 'Snapshot No. = ', i
      !
      !                             ! ================
      DO j = 1, nlo                 ! Horizontal slab
         !                          ! ================
         !
         ! Read one eddy-resolving snapshot at one layer
         call read_out (uo(:,:,j), vo(:,:,j), i, j) !! defined in 'ionc_subs.F'         !
         ! Coarse-graining of horizontal velocities (defined in 'coarse_subs.F')
         call filt (uf(:,:,j), ur(:,:,j), uo(:,:,j),                   & 
              &     gfk, wtot, 2*co*nso, nxtoc-1, nytoc) 
         call filt (vf(:,:,j), vr(:,:,j), vo(:,:,j),                   &
              &     gfk, wtot, 2*co*nso, nxtoc, nytoc-1)
         !
         !                          ! ================
      END DO                        ! End of slab
      !                             ! ================
      !
      ! Derive vertical velocities
      CALL uv2w (wf, uf, vf)
      CALL uv2w (wr, ur, vr)
      !
      ! Read wind stress 
      CALL read_out_sl (utau, vtau, i)
      !
      ! Compute wind velocity
      CALL tau2v (utau, uwnd)
      CALL tau2v (vtau, vwnd)
      !
      ! Coarse-graining of horizontal velocities (defined in 'coarse_subs.F')
      call filt (uwindf(:,:), uwindr(:,:), uwnd,                   & 
           &     gfk, wtot, 2*co*nso, nxtoc-1, nytoc) 
      call filt (vwindf(:,:), vwindr(:,:), vwnd,                   &
           &     gfk, wtot, 2*co*nso, nxtoc, nytoc-1)
      !
      ! Save output file
      CALL save_out (i) !! defined in 'ionc_subs.F'
      !
      !                             ! ================
   END DO                           !  End of time
   !                                ! ================

   !-------------------------------------------------------------------------------------
   !   Close file and free all the resources
   !-------------------------------------------------------------------------------------
   status = nf_close (uncid)
   if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
   status = nf_close (vncid)
   if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
   status = nf_close (oncid)
   if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

   print *
   print *, '*******************************************************'
   print *, '******************  END of program   ******************'
   print *, '*******************************************************'

END PROGRAM main    

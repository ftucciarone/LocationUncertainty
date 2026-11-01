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

   USE HR_fields
   USE LR_fields
   USE task_operations

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   !-------------------------------------------------------------------------------------
   !   Input/Output and Internal variables
   !-------------------------------------------------------------------------------------
   CHARACTER (len=*), PARAMETER                      :: subnam = 'main'
   CHARACTER (len=255)                               :: w1char, w2char, w3char, w4char, w5char
   CHARACTER (len=255)                               :: charName

   !-------------------------------------------------------------------------------------
   !   Local variables
   !-------------------------------------------------------------------------------------
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:)     :: gfk1, gfk2
   REAL(kind = wp )                                  :: wtot1, wtot2
   INTEGER                                           :: gwd1, gwd2 
   INTEGER                                           :: idt, idz, i, nt, time
   INTEGER                                           :: ierr(3)

   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:)     :: sume3u, sume3v
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:,:)   :: u_baro, v_baro, sshstore
   REAL(kind = wp ), ALLOCATABLE, DIMENSION(:,:,:,:) :: ustore, vstore, wstore

   LOGICAL                                           :: keRescale = .true.
   LOGICAL                                           :: taver, zaver
   REAL(kind = wp )                                  :: keHR, keLR

   !-------------------------------------------------------------------------------------
   !   netCDF variables 
   !-------------------------------------------------------------------------------------
   INTEGER                                           :: status, iostat, file_unit
   INTEGER, PARAMETER                                :: maxtime  = 690
   CHARACTER(len=255)                                :: cwd, fwd1, fwd2 

   CHARACTER (len=8)                                 :: subdir

   !-------------------------------------------------------------------------------------
   !   openMP Variables
   !-------------------------------------------------------------------------------------
   !$ INTEGER :: nprocs, OMP_GET_NUM_PROCS, nthmax, OMP_GET_MAX_THREADS


   !-------------------------------------------------------------------------------------
   !   Handles the command line parameters -w1, -w2, -units, -baroclinic and -timeavgrem
   !-------------------------------------------------------------------------------------

   IF(COMMAND_ARGUMENT_COUNT().NE.10)THEN
      print *, 'ERROR, two command-line arguments required'
      print *, 'SYNTAX:    ./CsGrain_BtBcSpt_ZeroMean.out -w1 int8 -w2 int8 -units char -baroclinic char -timeavgrem char '
      print *, 'SYNTAX:         with:'
      print *, 'SYNTAX:               -w1 int8      width of first filter'
      print *, 'SYNTAX:               -w2 int8      width of second filter'
      print *, 'SYNTAX:               -units        km/gridpoints   units of the width'
      print *, 'SYNTAX:               -baroclinic   yes/no, if yes removes the vertical average'
      print *, 'SYNTAX:               -timeavgrem   yes/no, if yes removes the HR time average '
      print *, 'STOPPING'
      STOP
   ENDIF
   CALL GET_COMMAND_ARGUMENT( 2,w1char)   !first, read the command-line values
   CALL GET_COMMAND_ARGUMENT( 4,w2char)
   CALL GET_COMMAND_ARGUMENT( 6,w3char)
   CALL GET_COMMAND_ARGUMENT( 8,w4char) 
   CALL GET_COMMAND_ARGUMENT(10,w5char)

   SELECT CASE (w4char)
   CASE('yes')
      charName = '_BtBcSpt'
      zaver = .true.
   CASE('no')
      charName = ''
      zaver = .false.
      CASE DEFAULT
      print *, 'SYNTAX ERROR in       -baroclinic   yes/no, if yes removes the vertical average'
      print *, 'STOPPING'
      STOP
   END SELECT

   SELECT CASE (w5char)
   CASE('yes')
      charName = TRIM(charName)//'_'//'ZeroMean'
      taver = .true.
   CASE('no')
      charName = TRIM(charName)//'_'//'WithMean'
      taver = .false.
   CASE DEFAULT
      print *, 'SYNTAX ERROR in       -timeavgrem   yes/no, if yes removes the HR time average'
      print *, 'STOPPING'
      STOP
   END SELECT


   !-------------------------------------------------------------------------------------
   !   Reading domain configuration and dimensions
   !-------------------------------------------------------------------------------------

   CALL rdnc_namelist ('namelist_ioNetCDF.nml')     ! Read address files

   CALL readFileID (subnam, dmdir, TRIM(dmgrid)//"_R27.nc", ncID1rd       ) ! HR domain file
   CALL readDim    (subnam,                           nm_x, ncID1rd,  nxHR)
   CALL readDim    (subnam,                           nm_y, ncID1rd,  nyHR)
   CALL readDim    (subnam,                      'nav_lev', ncID1rd, nloHR)

   CALL readFileID (subnam, indir,                 uncnam, uncID         ) ! HR domain file
   CALL readDim    (subnam,                          nm_t, uncID,   ntoHR)

   CALL readFileID (subnam, dmdir, TRIM(dmgrid)//"_R3.nc", domID       )    !  LR domain file
   CALL readDim    (subnam,                          nm_x, domID,  nxLR)
   CALL readDim    (subnam,                          nm_y, domID,  nyLR)
   CALL readDim    (subnam,                   TRIM(nm_z2), domID, nloLR)

   print *
   print *, 'HR Dimensions (nx,ny,nz): ', nxHR, nyHR, nloHR
   print *, 'LR Dimensions (nx,ny,nz): ', nxLR, nyLR, nloLR
   print *, 'HR        Time axis (nt): ', ntoHR
 
   CALL HR_task_init ('all')   
   CALL LR_task_init ('all')

   CALL HR_state_init ('band')
   CALL LR_state_init ('band')

   CALL    read2D_REALvar (subnam,  cn_ve1u, domID, e1uLR)
   CALL    read2D_REALvar (subnam,  cn_ve1v, domID, e1vLR)
   CALL    read2D_REALvar (subnam,  cn_ve1t, domID, e1tLR)
   CALL    read2D_REALvar (subnam,  cn_ve2u, domID, e2uLR)
   CALL    read2D_REALvar (subnam,  cn_ve2v, domID, e2vLR)
   CALL    read2D_REALvar (subnam,  cn_ve2t, domID, e2tLR)
   CALL read3Df4D_REALvar (subnam, cn_ve3u0, domID, e3uLR, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3v0, domID, e3vLR, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3t0, domID, e3tLR, 1)

   CALL    read2D_REALvar (subnam,  cn_ve1u, ncID1rd, e1uHR)
   CALL    read2D_REALvar (subnam,  cn_ve1v, ncID1rd, e1vHR)
   CALL    read2D_REALvar (subnam,  cn_ve1t, ncID1rd, e1tHR)
   CALL    read2D_REALvar (subnam,  cn_ve2u, ncID1rd, e2uHR)
   CALL    read2D_REALvar (subnam,  cn_ve2v, ncID1rd, e2vHR)
   CALL    read2D_REALvar (subnam,  cn_ve2t, ncID1rd, e2tHR)
   CALL read3Df4D_REALvar (subnam, cn_ve3u0, ncID1rd, e3uHR, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3v0, ncID1rd, e3vHR, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3t0, ncID1rd, e3tHR, 1)


   !-------------------------------------------------------------------------------------
   !   Derived grid parameters (do not alter)
   !-------------------------------------------------------------------------------------
   nso   = (nxHR - 2) / (nxLR - 2)   ! Resolutions Ratio (Command this through domain_cfg.nc)
   nxpo  = nxHR                      ! oceanic complete grid with boundaries in x 
   nypo  = nyHR                      ! oceanic complete grid with boundaries in y
   nxto  = nxpo - 2                  ! oceanic grid without boundaries in x 
   nyto  = nypo - 2                  ! oceanic grid without boundaries in x 
   nxtoc = nxto / nso                ! coarse inner gridcells in x
   nytoc = nyto / nso                ! coarse inner gridcells in y
   nxpoc = nxtoc + 2                 ! complete coarse grid points with boundaries in x
   nypoc = nytoc + 2                 ! complete coarse grid points with boundaries in y
!   co1   = 0
!   co2   = 2


   SELECT CASE(TRIM(w3char))

   CASE('gridpoints')                       ! if they are in gidpoints
      READ(w1char,*) co1                    ! convert them to INTEGER
      READ(w2char,*) co2
   CASE('km')                               ! Else, convert them in gridpoints
      READ(w1char,*) km1                    ! and then convert them to INTEGER
      READ(w2char,*) km2
      print*, "qui"
      HRgridsize = 0.5_wp * ( ( sum(e1tHR)/max(1,size(e1tHR)) + sum(e2tHR)/max(1,size(e2tHR)) ) ) / 1000
      co1 = int( km1 / (2 * HRgridsize) )
      co2 = int( km2 / (2 * HRgridsize) )
   CASE DEFAULT
      print *, 'ERROR, wrong command line argument given to units'
      print *, '         acceptable options are  km or gridpoints'
      print *, 'STOPPING'
      STOP
   END SELECT

   ! Then generate the filters
   gwd1 = 2 * nso * co1
   gwd2 = 2 * nso * co2



   ALLOCATE( gfk1(gwd1+1,gwd1+1) )
   ALLOCATE( gfk2(gwd2+1,gwd2+1) )
  
   print *
   print *,            '     Generate Gaussian filter kernel'

   CALL gaussf (gfk1, wtot1, gwd1)
   CALL gaussf (gfk2, wtot2, gwd2)

   write(*,'(a,i0)')   '         Width of Gaussian filter n.1= ', gwd1
   write(*,'(a,f6.3)') '  Total weight of Gaussian filter n.1= ', wtot1
   write(*,'(a,i0)')   '         Width of Gaussian filter n.2= ', gwd2
   write(*,'(a,f6.3)') '  Total weight of Gaussian filter n.2= ', wtot2


   !-------------------------------------------------------------------------------------
   !   Initialize outputs
   !-------------------------------------------------------------------------------------

   SELECT CASE(TRIM(w3char))

   CASE('gridpoints')                       ! if they are in gidpoints
      write(fwd1,'(i0)')  gwd1
      write(fwd2,'(i0)')  gwd2
      charName = TRIM(charName)//'_'//TRIM(fwd1)//'d1_'//TRIM(fwd2)//'d2'
      CALL wrnc_namelist ('namelist_ioNetCDF.nml', charName)
   CASE('km')                               ! Else, convert them in gridpoints
      charName = TRIM(charName)//'_'//TRIM(w1char)//'km_'//TRIM(w2char)//'km'
      CALL wrnc_namelist ('namelist_ioNetCDF.nml', charName)
   END SELECT

   IF (keRescale) THEN
      CALL wrnc_init (subnam, 'keHR_', HRoutID)
      CALL wrnc_init (subnam, 'keLR_', LRoutID)
   ELSE
      CALL wrnc_init (subnam, 'HR_', HRoutID)
      CALL wrnc_init (subnam, 'LR_', LRoutID)
   END IF


   CALL HR_save_init( .TRUE. )
   CALL LR_save_init( .TRUE., 690 )


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
   !   First loop on the time instants to evaluate the average
   !-------------------------------------------------------------------------------------
   !
   meanuHR = 0._wp
   meanvHR = 0._wp
   meanwHR = 0._wp  

   nt = 0
   time = 0 



   ALLOCATE(  sume3u(nxHR, nyHR), &
           &  sume3v(nxHR, nyHR)  )
   ALLOCATE(    ustore(nxHR, nyHR, nloHR, 720), &
           &    vstore(nxHR, nyHR, nloHR, 720), &
           &    wstore(nxHR, nyHR, nloHR, 720), &
           &  sshstore(nxHR, nyHR,        720)   )


   sume3u = SUM(e3uHR, DIM=3)
   sume3v = SUM(e3vHR, DIM=3)

   DO i = 0, 8, 2 
      !
      WRITE (subdir, "(I3,A1,I3,A1)") 100+i, "-", 102+i, "y"
      print *, "Opening from directory"
      print *, indir(1:25)//subdir//"/"
      !-------------------------------------------------------------------------------------
      !   Reading Velocity fields
      !-------------------------------------------------------------------------------------
      ! Open netCDF file 
      CALL readFileID (subnam, indir(1:25)//subdir//"/", uncnam, uncID)  ! U-grid
      CALL readFileID (subnam, indir(1:25)//subdir//"/", vncnam, vncID)  ! V-grid
      CALL readFileID (subnam, indir(1:25)//subdir//"/", trcnam, trcID)  ! T-grid
      !
      nt = nt + ntoHR
      !
      !
      !                                ! ================
      DO idt = 1, ntoHR                !  Time instant
         !                             ! ================
         !
         uHR_b = 0._wp
         vHR_b = 0._wp
         !
         time = time + 1
         !
         print *, 'Snapshot No. = ', idt, 'stored in: ', time
         !
         CALL read3Df4D_REALvar (subnam, cn_vozocrtx, uncID, ustore(:,:,:,time), idt)
         CALL read3Df4D_REALvar (subnam, cn_vomecrty, vncID, vstore(:,:,:,time), idt)
         !
         CALL read2Df3Dt_REALvar(subnam,cn_sossheig, trcID, sshstore(:,:,time), idt)
         !
         meanuHR = meanuHR + ustore(:,:,:,idt)
         meanvHR = meanvHR + vstore(:,:,:,idt)
         !
         !                             ! ================
      END DO                           !  End of time
      !                                ! ================
      !
      !                             ! ================
   END DO                           !  End of file
   !                                ! ================
   meanuHR = meanuHR / nt
   meanvHR = meanvHR / nt
   !
   CALL uv2w (meanwHR, meanuHR, meanvHR, e1tHR, e2tHR, e1vHR, e2uHR, e3uHR, e3vHR, e3tHR, nxHR, nyHR, nloHR)
   !
   !
   !
   CALL HR_save_onetime_out (   )


   !                             ! ================
   DO idz = 1, nloHR             !  Horizontal slab
      !                          ! ================
      CALL subsample (meanuLR(:,:,idz), meanuHR(:,:,idz), nxtoc, nytoc)
      CALL subsample (meanvLR(:,:,idz), meanvHR(:,:,idz), nxtoc, nytoc)  
      !                          ! ================
   END DO                        !  End of slab
   !                             ! ================
   !
   CALL uv2w (meanwLR, meanuLR, meanvLR, e1tLR, e2tLR, e1vLR, e2uLR, e3uLR, e3vLR, e3tLR, nxLR, nyLR, nloHR)
   !
   IF (keRescale) THEN
      !
      CALL  eke (keHR, meanuHR, meanvHR, meanwHR, e1uHR, e1vHR, e1tHR, e2uHR, e2vHR, e2tHR, e3uHR, e3vHR, e3wHR)
      CALL  eke (keLR, meanuLR, meanvLR, meanwLR, e1uLR, e1vLR, e1tLR, e2uLR, e2vLR, e2tLR, e3uLR, e3vLR, e3wLR)
      !
      meanuHR = meanuHR * ( keHR / keLR  )
      meanvHR = meanvHR * ( keHR / keLR  )
      meanwHR = meanwHR * ( keHR / keLR  )
      !
   END IF
   !
   CALL LR_save_onetime_out (   )
   !
   time = 0
   !
   !                                ! ================
   DO idt = 1, maxtime              !  Time instant
      !                             ! ================
      time = time + 1
      !
      flctuHR = 0._wp
      flctuHR = 0._wp
      flctuHR = 0._wp
      ! 
      uf1_LR = 0._wp
      vf1_LR = 0._wp
      wf1_LR = 0._wp
      sshf1_LR = 0._wp
      !
      uf2_LR = 0._wp
      vf2_LR = 0._wp
      wf2_LR = 0._wp
      sshf2_LR = 0._wp
      !
      uband_LR = 0._wp
      vband_LR = 0._wp
      wband_LR = 0._wp
      sshband_LR = 0._wp
      ! 
      print *, 'Snapshot No. = ', idt
      !
      IF ( taver ) THEN
         !
         ustore(:,:,:,idt) = ustore(:,:,:,idt) - meanuHR
         vstore(:,:,:,idt) = vstore(:,:,:,idt) - meanvHR
         !
      END IF
      !
      !
      uHR_b = SUM(ustore(:,:,:,time) * e3uHR, DIM=3) / sume3u
      vHR_b = SUM(vstore(:,:,:,time) * e3vHR, DIM=3) / sume3v
      !
      CALL subsample (uLR_b(:,:), uHR_b(:,:), nxtoc-1, nytoc)
      CALL subsample (vLR_b(:,:), vHR_b(:,:), nxtoc-1, nytoc)  
      ! 
      !
      IF ( zaver ) THEN
         !
         ustore(:,:,:,idt) = ustore(:,:,:,idt) - SPREAD( uHR_b(:,:), 3, nloHR)
         vstore(:,:,:,idt) = vstore(:,:,:,idt) - SPREAD( vHR_b(:,:), 3, nloHR)
         !
      END IF
      !
      ustore(:,:,nloHR,idt) = 0._wp
      vstore(:,:,nloHR,idt) = 0._wp
      !
 
      CALL subsample (sshLR(:,:), sshstore(:,:,time), nxtoc, nytoc)
      CALL DifferenceOfGaussians(sshband_LR(:,:), sshf1_LR(:,:), sshf2_LR(:,:), sshstore(:,:,idt), &
           &                  gfk1, wtot1, gwd1, gfk2, wtot2, gwd2, nxtoc, nytoc  )

      !                             ! ================
      DO idz = 1, nloHR             !  Horizontal slab
         !                          ! ================
         !
         CALL subsample (uLR(:,:,idz), ustore(:,:,idz,time), nxtoc-1, nytoc)
         CALL subsample (vLR(:,:,idz), vstore(:,:,idz,time), nxtoc-1, nytoc)  
         !
         CALL DifferenceOfGaussians(uband_LR(:,:,idz), uf1_LR(:,:,idz), uf2_LR(:,:,idz), ustore(:,:,idz,idt), &
              &                  gfk1, wtot1, gwd1, gfk2, wtot2, gwd2, nxtoc, nytoc )
         CALL DifferenceOfGaussians(vband_LR(:,:,idz), vf1_LR(:,:,idz), vf2_LR(:,:,idz), vstore(:,:,idz,idt), &
              &                  gfk1, wtot1, gwd1, gfk2, wtot2, gwd2, nxtoc, nytoc )
         !
         !
         !                          ! ================
      END DO                        !  End of slab
      !                             ! ================
      !
      CALL uv2w (     wLR,      uLR,      vLR, e1tLR, e2tLR, e1vLR, e2uLR, e3uLR, e3vLR, e3tLR, nxLR, nyLR, nloHR)
      !
      ! Derive vertical velocities
      CALL uv2w (  wf1_LR,   uf1_LR,   vf1_LR, e1tLR, e2tLR, e1vLR, e2uLR, e3uLR, e3vLR, e3tLR, nxLR, nyLR, nloHR)
      CALL uv2w (  wf2_LR,   uf2_LR,   vf2_LR, e1tLR, e2tLR, e1vLR, e2uLR, e3uLR, e3vLR, e3tLR, nxLR, nyLR, nloHR)
      CALL uv2w (wband_LR, uband_LR, vband_LR, e1tLR, e2tLR, e1vLR, e2uLR, e3uLR, e3vLR, e3tLR, nxLR, nyLR, nloHR)
      !
      CALL uv2w (wstore(:,:,:,idt), &
         &       ustore(:,:,:,idt), &
         &       vstore(:,:,:,idt),            e1tHR, e2tHR, e1vHR, e2uHR, e3uHR, e3vHR, e3tHR, nxHR, nyHR, nloHR)
      !
      !
      IF (keRescale) THEN
         !      
         CALL  eke (keHR, ustore(:,:,:,idt), &
                  &       vstore(:,:,:,idt), & 
                  &       wstore(:,:,:,idt), e1uHR, e1vHR, e1tHR, e2uHR, e2vHR, e2tHR, e3uHR, e3vHR, e3wHR)
         CALL  eke (keLR, uband_LR, vband_LR, wband_LR, e1uLR, e1vLR, e1tLR, e2uLR, e2vLR, e2tLR, e3uLR, e3vLR, e3wLR)
         !
         uband_LR = uband_LR * ( keHR / keLR  )
         vband_LR = vband_LR * ( keHR / keLR  )
         wband_LR = wband_LR * ( keHR / keLR  )
         !
      END IF
      !
      !
      ! Save output file
      CALL LR_save_out (time) !! defined in 'ionc_subs.F'
      !
      !                             ! ================
   END DO                           !  End of time
   !                                ! ================
   !
   !
   !
   !
   print *
   print *, '*******************************************************'
   print *, '******************  END of program   ******************'
   print *, '*******************************************************'

END PROGRAM main    







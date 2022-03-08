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

!    Modules
   USE class_precision        ! Defines single and double precision
   USE ionml_subs             ! Manages the namelist
   USE NEMO_nmconv            ! Defines NEMO conventions for naming
   USE state
   USE rdnc_subs

   IMPLICIT NONE
   INCLUDE 'netcdf.inc'

   !-------------------------------------------------------------------------------------
   !   Local variables
   !-------------------------------------------------------------------------------------
   REAL( kind = wp ), ALLOCATABLE, DIMENSION(:,:,:) :: uin,vin,win,wtilde     ! (m/s)
   INTEGER                                          :: i,j
      INTEGER :: ierr(3)
   !-------------------------------------------------------------------------------------
   !   netCDF variables 
   !-------------------------------------------------------------------------------------
   CHARACTER (len=*), PARAMETER                     :: subnam = 'main'
   INTEGER                                          :: status, iostat, file_unit

   !-------------------------------------------------------------------------------------
   !   Reading namelist
   !-------------------------------------------------------------------------------------
   NAMELIST/array_param/ nxto, nyto, nlo, nto, nso, co 



   CALL  rdnc_namelist ('namelist_inputNetCDF.nml')
   CALL open_inputfile('namelist_coarsegrain.nml', file_unit, iostat)
   READ (NML=array_param, IOSTAT=iostat, UNIT=file_unit)
   CALL close_inputfile('namelist_coarsegrain.nml', file_unit, iostat)
   !-------------------------------------------------------------------------------------
   !   Reading domain configuration and dimensions
   !-------------------------------------------------------------------------------------
   CALL readFileID (subnam, iodir, dmgrid, domID)  ! domain file
   CALL readDim (subnam,      cn_x, domID, nxpo)
   CALL readDim (subnam,      cn_y, domID, nypo)
   CALL readDim (subnam, 'nav_lev', domID,  nlo)

   print *, 'Dimensions (nx,ny,nz): ', nxpo, nypo, nlo
 
   ALLOCATE(    uin(nxpo,nypo,nlo),   &
        &       vin(nxpo,nypo,nlo),   &
        &       win(nxpo,nypo,nlo),   &
        &    wtilde(nxpo,nypo,nlo)    )

   !-------------------------------------------------------------------------------------
   !   Reading Velocity fields
   !-------------------------------------------------------------------------------------
   ! Open netCDF file 
   CALL readFileID (subnam, iodir, trcnam, trcID)  ! T-grid
   CALL readFileID (subnam, iodir, uncnam, uncID)  ! U-grid
   CALL readFileID (subnam, iodir, vncnam, vncID)  ! V-grid
   CALL readFileID (subnam, iodir, wncnam, wncID)  ! W-grid


   CALL readDim (subnam,      cn_t, wncID,  nto)
   print *, '       Time axis (nt): ', nto


   ALLOCATE( e2u(nxpo,nypo)    ,   &        
      &      e1v(nxpo,nypo)    ,   &
      &      e1t(nxpo,nypo)    ,   &
      &      e2t(nxpo,nypo)    ,   &
      &      e3u(nxpo,nypo,nlo),   &
      &      e3v(nxpo,nypo,nlo),   &
      &      e3t(nxpo,nypo,nlo), STAT=ierr(2) )


   CALL read2D_REALvar (subnam, cn_ve2u, domID, e2u)
   CALL read2D_REALvar (subnam, cn_ve1v, domID, e1v)
   CALL read2D_REALvar (subnam, cn_ve1t, domID, e1t)
   CALL read2D_REALvar (subnam, cn_ve2t, domID, e2t)
   CALL read3Df4D_REALvar (subnam, cn_ve3u0, domID, e3u, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3v0, domID, e3v, 1)
   CALL read3Df4D_REALvar (subnam, cn_ve3t0, domID, e3t, 1)




   print *, MAXVAL(ABS( e2u)), MAXVAL(ABS( e1v)), MAXVAL(ABS( e1t)), MAXVAL(ABS( e2t))
   print *, MAXVAL(ABS( e3t))



   !
   !                                ! ================
   DO i = 1, nto                    !  Time instant
      !                             ! ================
      !
      uin = 0._wp
      vin = 0._wp
      win = 0._wp
      wtilde = 0._wp
      !
      print *, 'Snapshot No. = ', i
      !
      CALL read3Df4D_REALvar (subnam, cn_vozocrtx, uncID, uin, i)
      CALL read3Df4D_REALvar (subnam, cn_vomecrty, vncID, vin, i)
      CALL read3Df4D_REALvar (subnam, cn_vovecrtz, wncID, win, i)

      !
      ! Derive vertical velocities
      CALL uv2w (wtilde, uin, vin)
      print *, MAXVAL(ABS( win )), MAXVAL(ABS( wtilde )), MAXVAL(ABS( win - wtilde))

      !
      !
      !                             ! ================
   END DO                           !  End of time
   !                                ! ================


   print *
   print *, '*******************************************************'
   print *, '******************  END of program   ******************'
   print *, '*******************************************************'

END PROGRAM main    





   SUBROUTINE uv2w (w, u, v)

!     Given horizontal velocities u and v, compute the vertical
!     component w using divergence-free condition.

!     Modules
      USE param
      USE state

      IMPLICIT NONE    
      
!     I/O arguments
      REAL(KIND = wp), INTENT(  in ), DIMENSION(nxpo,nypo,nlo) :: u,v
      REAL(KIND = wp), INTENT( out ), DIMENSION(nxpo,nypo,nlo) :: w
!     Local variables      
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)        :: hdivn
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)        :: tmask
      INTEGER                                               :: ierr
      INTEGER                                               :: ji,jj,jk


      ALLOCATE( hdivn(nxpo,nypo,nlo) ,  &
              & tmask(nxpo,nypo,nlo) ,  STAT=ierr )
      IF( ierr .ne. 0 )  print *, 'Allocation error in uv2w'

      w = 0._wp !! null boundaries
      hdivn = 0._wp
      tmask = 0._wp


      DO jk = 1, nlo-1                                      !==  Horizontal divergence  ==!
         DO jj = 2, nypo-1
            DO ji = 2, nxpo-1   ! vector opt.
               hdivn(ji,jj,jk) = (  e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk) * u(ji  ,jj  ,jk)       &
                  &               - e2u(ji-1,jj  ) * e3u(ji-1,jj  ,jk) * u(ji-1,jj  ,jk)       &
                  &               + e1v(ji  ,jj  ) * e3v(ji  ,jj  ,jk) * v(ji  ,jj  ,jk)       &
                  &               - e1v(ji  ,jj-1) * e3v(ji  ,jj-1,jk) * v(ji  ,jj-1,jk)   )   & 
                  &             / ( e3t(ji  ,jj  ,jk) * e1t(ji,jj) * e2t(ji,jj) )
         !      IF ( ( b_lev(ji,jj) <= jk ) .or. ( jk <= t_lev(ji,jj) ) ) THEN
         !         tmask(ji,jj,jk) = 1._wp
         !      END IF 

            END DO  
         END DO  
      END DO

      DO jk = nlo-1, 1, -1                       ! integrate from the bottom the hor. divergence
         ! computation of w
         w(:,:,jk) = w(:,:,jk+1) - ( e3t(:,:,jk) * hdivn(:,:,jk) )! * tmask(:,:,jk)
      END DO

   END SUBROUTINE uv2w

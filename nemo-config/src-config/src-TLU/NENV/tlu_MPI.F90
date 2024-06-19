!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_newEnv
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Brief resume
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!                
!!                
!> @snippet       this tlu_noise
!!
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the second dimension of the spatial arrays
!! @param[in]     jpk: the third dimension of the spatial arrays
!!
!> @par           Other modules dependencies
!> @snippet       this mod_dep
!!
!> @par           Public variables
!> @snippet       this public_vars
!!
!> @par           Public routines
!> @snippet       this public_sub
!!
!> @todo          
!! @todo          
!!
!!------------------------------------------------------------------------------
MODULE tluMPI
   ! [mod_dep]
   USE par_kind           ! data types defined in par_kind module
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce            ! ocean dynamics parameters
   USE oce                ! ocean dynamics and active tracers
   USE dom_oce            ! ocean space and time domain
   USE lib_mpp            ! MPP library
   USE timing             ! Timing
   USE lbclnk             ! ocean lateral boundary conditions (or mpp link)
   USE tlu
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC add_extra_halo      ! Indicate which subroutine calls it
   PUBLIC add_extra_halo_2D   ! 
   PUBLIC MPI_collect         ! Indicate which subroutine calls it
   PUBLIC MPI_collect_noi     ! Indicate which subroutine calls it
  !
   INCLUDE 'mpif.h'
   !
   REAL(wp),    PUBLIC    :: mod_var
   !
CONTAINS


   SUBROUTINE add_extra_halo( kt, field_in, field_out, extra_halo, field_in_name )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_newsub  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!-----------------------------------------------------------------------
      INTEGER(i4),                              INTENT(in   ) :: kt
      REAL(wp),         DIMENSION(jpi,jpj,jpk), INTENT(in   ) :: field_in
      REAL(wp),         DIMENSION(jpi + 2*extra_halo,jpj+ 2*extra_halo,jpk), INTENT(  out) :: field_out
      INTEGER(i4),                              INTENT(in   ) :: extra_halo
      CHARACTER(len=2),                         INTENT(in   ) :: field_in_name
      ! MPI variables
      INTEGER :: ierror, status, tag, my_rank
      !
      REAL(wp)  :: lint(        extra_halo,         jpj-2, jpk), rint(        extra_halo,         jpj-2, jpk)
      REAL(wp)  :: lext(        extra_halo,         jpj-2, jpk), rext(        extra_halo,         jpj-2, jpk)
      !
      REAL(wp)  :: bint(jpi + 2*extra_halo, extra_halo +1, jpk), tint(jpi + 2*extra_halo, extra_halo +1, jpk)
      REAL(wp)  :: bext(jpi + 2*extra_halo, extra_halo +1, jpk), text(jpi + 2*extra_halo, extra_halo +1, jpk)
      !
      REAL(wp)  :: field(jpi + 2*extra_halo, jpj + 2*extra_halo, jpk)
      !
      ! Control print (1)
      !
      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'add_extra_halo : Add an extra halo to the fieldi ', TRIM(field_in_name)
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      END IF
      lint = 0._wp
      rint = 0._wp
      tint = 0._wp
      bint = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      lint = field_in(   3                : 3 + extra_halo - 1 , 2:jpj-1, :)
      rint = field_in( jpim1 - extra_halo :          jpim1 - 1 , 2:jpj-1, :)
      lext = 0._wp
      rext = 0._wp

      tag = 65464
      !
      ! Get my rank and do the corresponding job
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
      !
      ! West to East
      !
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_send(rint, size(rint), MPI_double_precision, my_rank+1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_Recv(lext, size(lext), MPI_double_precision, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! East to West
      !
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_send(lint, size(lint), MPI_double_precision, my_rank-1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_Recv(rext, size(rext), MPI_double_precision, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! Safely wait for West2East2West communication
      !
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      tag=1
      !
      ! Ensure transfer of extended borders
      !
      field_out = 0._wp 
      field_out( 1                    :     extra_halo       , 2 + extra_halo : extra_halo + jpj - 1, :) = lext
      field_out( 1 + extra_halo       :     extra_halo + jpi , 2 + extra_halo : extra_halo + jpj - 1, :) = field_in(:, 2 : jpjm1,:)
      field_out( 1 + extra_halo + jpi : 2 * extra_halo + jpi , 2 + extra_halo : extra_halo + jpj - 1, :) = rext
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      tint = field_out(:,      jpjm1     :     extra_halo + jpjm1, :)
      bint = field_out(:, extra_halo + 2 : 2 * extra_halo + 2  , :)
      text = 0._wp
      bext = 0._wp
      !
      ! South to North
      !
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_send(tint, size(tint), MPI_double_precision, my_rank+jpni,         tag, MPI_COMM_WORLD, ierror);
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_Recv(bext, size(bext), MPI_double_precision, my_rank-jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! North to South
      !
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_send(bint, size(bint), MPI_double_precision, my_rank-jpni,         tag, MPI_COMM_WORLD, ierror)
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_Recv(text, size(text), MPI_double_precision, my_rank+jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      !
      field_out(:,                1 :       1 + extra_halo, :) = bext
      field_out(:, extra_halo + jpj : jpj + 2 * extra_halo, :) = text
      !
   END SUBROUTINE add_extra_halo


   SUBROUTINE add_extra_halo_2D( kt, field_in, field_out, extra_halo, field_in_name )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_newsub  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!-----------------------------------------------------------------------
      INTEGER(i4),                              INTENT(in   ) :: kt
      REAL(wp),             DIMENSION(jpi,jpj), INTENT(in   ) :: field_in
      REAL(wp),         DIMENSION(jpi + 2*extra_halo,jpj+ 2*extra_halo), INTENT(  out) :: field_out
      INTEGER(i4),                              INTENT(in   ) :: extra_halo
      CHARACTER(len=6),                         INTENT(in   ) :: field_in_name
      ! MPI variables
      INTEGER :: ierror, status, tag, my_rank
      !
      REAL(wp)  :: lint(        extra_halo,         jpj-2), rint(        extra_halo,         jpj-2)
      REAL(wp)  :: lext(        extra_halo,         jpj-2), rext(        extra_halo,         jpj-2)
      !
      REAL(wp)  :: bint(jpi + 2*extra_halo, extra_halo +1), tint(jpi + 2*extra_halo, extra_halo +1)
      REAL(wp)  :: bext(jpi + 2*extra_halo, extra_halo +1), text(jpi + 2*extra_halo, extra_halo +1)
      !
      REAL(wp)  :: field(jpi + 2*extra_halo, jpj + 2*extra_halo)
      !
      ! Control print (1)
      !
       IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'add_extra_halo : Add an extra halo to the fieldi ', TRIM(field_in_name)
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      END IF
      lint = 0._wp
      rint = 0._wp
      tint = 0._wp
      bint = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      lint = field_in(   3                : 3 + extra_halo - 1 , 2:jpj-1)
      rint = field_in( jpim1 - extra_halo :          jpim1 - 1 , 2:jpj-1)
      lext = 0._wp
      rext = 0._wp

      tag = 65464
      
      !
      ! Get my rank and do the corresponding job
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
      !
      ! West to East
      !
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_send(rint, size(rint), MPI_double_precision, my_rank+1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_Recv(lext, size(lext), MPI_double_precision, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! East to West
      !
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_send(lint, size(lint), MPI_double_precision, my_rank-1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_Recv(rext, size(rext), MPI_double_precision, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! Safely wait for West2East2West communication
      !
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      tag=1
      !
      ! Ensure transfer of extended borders
      !
      field_out = 0._wp 
      field_out( 1                    :     extra_halo       , 2 + extra_halo : extra_halo + jpj - 1) = lext
      field_out( 1 + extra_halo       :     extra_halo + jpi , 2 + extra_halo : extra_halo + jpj - 1) = field_in(:, 2 : jpjm1)
      field_out( 1 + extra_halo + jpi : 2 * extra_halo + jpi , 2 + extra_halo : extra_halo + jpj - 1) = rext
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      tint = field_out(:,      jpjm1     :     extra_halo + jpjm1)
      bint = field_out(:, extra_halo + 2 : 2 * extra_halo + 2    )    
      text = 0._wp
      bext = 0._wp
      !
      ! South to North
      !
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_send(tint, size(tint), MPI_double_precision, my_rank+jpni,         tag, MPI_COMM_WORLD, ierror);
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_Recv(bext, size(bext), MPI_double_precision, my_rank-jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! North to South
      !
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_send(bint, size(bint), MPI_double_precision, my_rank-jpni,         tag, MPI_COMM_WORLD, ierror)
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_Recv(text, size(text), MPI_double_precision, my_rank+jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      !
      field_out(:,                1 :       1 + extra_halo) = bext
      field_out(:, extra_halo + jpj : jpj + 2 * extra_halo) = text
      !
   END SUBROUTINE add_extra_halo_2D

   SUBROUTINE MPI_collect( kt, field_in, field_out, ehalo, field_in_name )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_newsub  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!-----------------------------------------------------------------------
      INTEGER(i4),                              INTENT(in   ) :: kt
      REAL(wp),         DIMENSION(jpi,jpj,jpk), INTENT(in   ) :: field_in
      REAL(wp),         DIMENSION(jpi + 2*ehalo,jpj+ 2*ehalo,jpk), INTENT(  out) :: field_out
      INTEGER(i4),                              INTENT(in   ) :: ehalo
      CHARACTER(len=2),                         INTENT(in   ) :: field_in_name
      INTEGER(i4)                                             :: buf
      ! MPI variables
      INTEGER :: ierror, status, tag, my_rank
      !
      REAL(wp)  :: lint( ehalo + 1,         jpj-2, jpk), rint( ehalo + 1,         jpj-2, jpk)
      REAL(wp)  :: lext( ehalo + 1,         jpj-2, jpk), rext( ehalo + 1,         jpj-2, jpk)
      !
      REAL(wp)  ::  bint(jpi + 2*ehalo, ehalo +1, jpk), tint(jpi + 2*ehalo, ehalo +1, jpk)
      REAL(wp)  ::  bext(jpi + 2*ehalo, ehalo +1, jpk), text(jpi + 2*ehalo, ehalo +1, jpk)
      !
      REAL(wp)  :: field(jpi + 2*ehalo, jpj + 2*ehalo, jpk)
      !
      ! Control print (1)
      !
      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'add_extra_halo : Add an extra halo to the fieldi ', TRIM(field_in_name)
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      END IF
      !
      buf = ehalo + 1
      !
      lint = 0._wp
      rint = 0._wp
      tint = 0._wp
      bint = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      lint = field_in(   2                : 2 + ehalo , 2:jpj-1, :)
      rint = field_in( jpim1 - ehalo :          jpim1 , 2:jpj-1, :)
      lext = 0._wp
      rext = 0._wp

      tag = 65464
      !
      ! Get my rank and do the corresponding job
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
      !
      ! West to East
      !
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_send(rint, size(rint), MPI_double_precision, my_rank+1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_Recv(lext, size(lext), MPI_double_precision, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! East to West
      !
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_send(lint, size(lint), MPI_double_precision, my_rank-1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_Recv(rext, size(rext), MPI_double_precision, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! Safely wait for West2East2West communication
      !
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      tag=1
      !
      ! Ensure transfer of extended borders
      !
      field_out = 0._wp
      field_out( 1               :     buf             , 1 + buf : buf + jpjm1 - 1, :) = lext
      field_out( 1 + buf         :     buf + jpim1 - 1 , 1 + buf : buf + jpjm1 - 1, :) = field_in( 2 : jpim1, 2 : jpjm1, :)
      field_out(     buf + jpim1 : 2 * ehalo + jpi     , 1 + buf : buf + jpjm1 - 1, :) = rext
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      tint = field_out(:, jpjm1     :     ehalo + jpjm1, :)
      bint = field_out(:, ehalo + 2 : 2 * ehalo + 2  , :)
      text = 0._wp
      bext = 0._wp
      !
      ! South to North
      !
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_send(tint, size(tint), MPI_double_precision, my_rank+jpni,         tag, MPI_COMM_WORLD, ierror);
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_Recv(bext, size(bext), MPI_double_precision, my_rank-jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! North to South
      !
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_send(bint, size(bint), MPI_double_precision, my_rank-jpni,         tag, MPI_COMM_WORLD, ierror)
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_Recv(text, size(text), MPI_double_precision, my_rank+jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      !
      field_out(:,           1 :       1 + ehalo, :) = bext
      field_out(:, ehalo + jpj : jpj + 2 * ehalo, :) = text
      !
   END SUBROUTINE MPI_collect


   SUBROUTINE MPI_collect_noi( kt, field_in, field_out, ehalo, field_in_name )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_newsub  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   
      !!
      !!-----------------------------------------------------------------------
      INTEGER(i4),                          INTENT(in   ) :: kt
      REAL(wp),         DIMENSION(jpi,jpj), INTENT(in   ) :: field_in
      REAL(wp),         DIMENSION(jpi + 2*ehalo,jpj+ 2*ehalo), INTENT(inout) :: field_out
      INTEGER(i4),                          INTENT(in   ) :: ehalo
      CHARACTER(len=2),                     INTENT(in   ) :: field_in_name
      INTEGER(i4)                                         :: buf
      ! MPI variables
      INTEGER :: ierror, status, tag, my_rank
      !
      REAL(wp)  :: lint( ehalo + 1,         jpj-2), rint( ehalo + 1,         jpj-2)
      REAL(wp)  :: lext( ehalo + 1,         jpj-2), rext( ehalo + 1,         jpj-2)
      !
      REAL(wp)  ::  bint(jpi + 2*ehalo, ehalo +1), tint(jpi + 2*ehalo, ehalo +1)
      REAL(wp)  ::  bext(jpi + 2*ehalo, ehalo +1), text(jpi + 2*ehalo, ehalo +1)
      !
      REAL(wp)  :: field(jpi + 2*ehalo, jpj + 2*ehalo)
      !
      ! Control print (1)
      !
      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'add_extra_halo : Add an extra halo to the fieldi ', TRIM(field_in_name)
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      END IF
      !
      buf = ehalo + 1
      !
   !   lint = 0._wp
   !   rint = 0._wp
   !   tint = 0._wp
   !   bint = 0._wp
      ! 
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      lint = field_in(   2                : 2 + ehalo , 2:jpj-1)
      rint = field_in( jpim1 - ehalo :          jpim1 , 2:jpj-1)
      lext = field_out( 1               :     buf             , 1 + buf : buf + jpjm1 - 1) ! 0._wp
      rext = field_out(     buf + jpim1 : 2 * ehalo + jpi     , 1 + buf : buf + jpjm1 - 1) ! 0._wp

      tag = 65464
      !
      ! Get my rank and do the corresponding job
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
      !
      ! West to East
      !
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_send(rint, size(rint), MPI_double_precision, my_rank+1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_Recv(lext, size(lext), MPI_double_precision, my_rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! East to West
      !
      IF ( MOD(my_rank, jpni) .ne.      0 ) CALL MPI_send(lint, size(lint), MPI_double_precision, my_rank-1,         tag, MPI_COMM_WORLD, ierror)
      IF ( MOD(my_rank, jpni) .ne. jpni-1 ) CALL MPI_Recv(rext, size(rext), MPI_double_precision, my_rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! Safely wait for West2East2West communication
      !
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      tag=1
      !
      ! Ensure transfer of extended borders
      !
   !   field_out = 0._wp
      field_out( 1               :     buf             , 1 + buf : buf + jpjm1 - 1) = lext
      field_out( 1 + buf         :     buf + jpim1 - 1 , 1 + buf : buf + jpjm1 - 1) = field_in( 2 : jpim1, 2 : jpjm1)
      field_out(     buf + jpim1 : 2 * ehalo + jpi     , 1 + buf : buf + jpjm1 - 1) = rext
      !
      ! Copy array into contiguous memory space to avoid MPI problems
      !
      tint = field_out(:, jpjm1     :     ehalo + jpjm1)
      bint = field_out(:, ehalo + 2 : 2 * ehalo + 2  )
      text = field_out(:, ehalo + jpj : jpj + 2 * ehalo) ! 0._wp
      bext = field_out(:,           1 :       1 + ehalo) ! 0._wp
      !
      ! South to North
      !
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_send(tint, size(tint), MPI_double_precision, my_rank+jpni,         tag, MPI_COMM_WORLD, ierror);
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_Recv(bext, size(bext), MPI_double_precision, my_rank-jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      ! North to South
      !
      IF ( INT(my_rank/jpni)    .ne.    0 ) CALL MPI_send(bint, size(bint), MPI_double_precision, my_rank-jpni,         tag, MPI_COMM_WORLD, ierror)
      IF ( INT(my_rank/jpni) +1 .ne. jpnj ) CALL MPI_Recv(text, size(text), MPI_double_precision, my_rank+jpni, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
      !
      !
      field_out(:,           1 :       1 + ehalo) = bext
      field_out(:, ehalo + jpj : jpj + 2 * ehalo) = text
      !
   END SUBROUTINE MPI_collect_noi


END MODULE tluMPI




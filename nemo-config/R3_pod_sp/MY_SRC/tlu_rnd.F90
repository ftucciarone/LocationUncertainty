!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_rnd
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                F. Tucciarone
!!
!> @version
!!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Contains the random number generation routines
!! 
!!------------------------------------------------------------------------------
MODULE tlurnd
   ! [mod_dep]
   USE par_kind           ! data types defined in par_kind modul
   USE in_out_manager     ! I/O manager
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC Seed
   PUBLIC Real_Uniform_Interval
   PUBLIC Integer_Uniform_Interval
   PUBLIC Standard_Gaussian
   PUBLIC Standard_Gaussian_2D
   !
CONTAINS

   SUBROUTINE Seed( n, static )
      INTEGER(4)                                        :: n
      INTEGER(4), DIMENSION(n), OPTIONAL, INTENT(in   ) :: static
      !
      CALL random_seed( size=n )
      IF ( PRESENT(static) ) CALL random_seed( size=n, put=static )
      !
   END SUBROUTINE Seed

   SUBROUTINE Real_Uniform_Interval( a, b, x, n )
      INTEGER(4),               INTENT(in   ) :: n
      INTEGER(4),               INTENT(in   ) :: a, b
      INTEGER(4), DIMENSION(n), INTENT(  out) :: x
      REAL(wp),   DIMENSION(n)                :: u
      !
      CALL random_number(u)
      x = INT( a + FLOOR( (b-a)*u ) , i4)
      !
   END SUBROUTINE Real_Uniform_Interval

   SUBROUTINE Integer_Uniform_Interval( a, b, x, n )
      INTEGER(4),               INTENT(in   ) :: n
      REAL(wp),                 INTENT(in   ) :: a, b
      REAL(wp),   DIMENSION(n), INTENT(  out) :: x
      REAL(wp),   DIMENSION(n)                :: u
      !
      CALL random_number(u)
      x = (b-a)*u + a
      !
   END SUBROUTINE Integer_Uniform_Interval

   SUBROUTINE Standard_Gaussian( n, x, allones )
      INTEGER(4),                          INTENT(in   ) :: n
      LOGICAL,                  OPTIONAL,  INTENT(in   ) :: allones
      REAL(wp),   DIMENSION(n),            INTENT(  out) :: x
      REAL(wp),   DIMENSION(n)                           :: u, v
      REAL(wp),   PARAMETER                              :: pi = 4._wp * atan(1._wp)
      !
      CALL random_number(u)
      CALL random_number(v)
      !
      x = SQRT( -2*LOG(u) ) * COS( 2*pi*v )
      !
      IF ( PRESENT(allones) ) x = 1._wp
      !
   END SUBROUTINE Standard_Gaussian

   SUBROUTINE Standard_Gaussian_2D( n, m, x, allones )
      INTEGER(4),                            INTENT(in   ) :: n, m
      LOGICAL,                    OPTIONAL,  INTENT(in   ) :: allones
      REAL(wp),   DIMENSION(n,m),            INTENT(  out) :: x
      REAL(wp),   DIMENSION(n,m)                           :: u, v
      REAL(wp),   PARAMETER                                :: pi = 4._wp * atan(1._wp)
      !
      CALL random_number(u)
      CALL random_number(v)
      !
      x = SQRT( -2*LOG(u) ) * COS( 2*pi*v )
      !
      IF ( PRESENT(allones) ) x = 1._wp
      !
   END SUBROUTINE Standard_Gaussian_2D

  SUBROUTINE  r8vec_normal_01 ( kt, n, x, allones )
   !!-------------------------------------------------------------------------------
   !!
   !! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
   !!
   !!  Discussion:
   !!
   !!    An R8VEC is an array of double precision real values.
   !!
   !!    The standard normal probability distribution function (PDF) has
   !!    mean 0 and standard deviation 1.
   !!
   !!  Licensing:
   !!
   !!    This code is distributed under the GNU LGPL license.
   !!
   !!  Modified:
   !!
   !!    18 May 2014
   !!
   !!  Author:
   !!
   !!    John Burkardt
   !!
   !!  Parameters:
   !!
   !!    Input, integer N, the number of values desired.
   !!
   !!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
   !!
   !!  Local:
   !!
   !!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
   !!    random values.  Its dimension is N+1, but really it is only needed
   !!    to be the smallest even number greater than or equal to N.
   !!
   !!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range
   !!    of entries of X that we need to compute
   !!
   !!-------------------------------------------------------------------------------
      IMPLICIT NONE
      ! Dummy variables
      INTEGER,  INTENT(in   )                 :: kt         ! ocean time-step index
      LOGICAL,  INTENT(in   )                 :: allones
      INTEGER,  INTENT(in   )                 :: n
      REAL(wp), INTENT(  out), DIMENSION( n ) :: x 
      ! Local variables
      INTEGER                                 :: m
      INTEGER                                 :: x_lo_index
      INTEGER                                 :: x_hi_index
      REAL(wp),                DIMENSION(n+1) :: r
      REAL(wp), PARAMETER                     :: r8_pi = 4._wp * atan(1._wp)

      !
      !  Record the range of X we need to fill in.
      !
      x_lo_index = 1
      x_hi_index = n
      !
      !  If we need just one new value, do that here to avoid null arrays.
      !
      IF ( x_hi_index - x_lo_index + 1 == 1 ) THEN

         CALL random_number ( harvest = r(1:2) )

         x(x_hi_index) = SQRT ( - 2._wp * LOG ( r(1) ) ) * COS ( 2._wp * r8_pi * r(2) )
      !
      !  If we require an even number of values, that's easy.
      !
      ELSE IF ( MOD ( x_hi_index - x_lo_index, 2 ) == 1 ) THEN

         m = ( x_hi_index - x_lo_index + 1 ) / 2

         call random_number ( harvest = r(1:2*m) )

         x(x_lo_index:x_hi_index-1:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-1:2) ) ) &
                                    & *  COS (   2._wp * r8_pi * r(2:2*m:2) )

         x(x_lo_index+1:x_hi_index:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-1:2) ) ) &
                                    & *  SIN (   2._wp * r8_pi * r(2:2*m:2) )
      !
      !  If we require an odd number of values, we generate an even number,
      !  and handle the last pair specially, storing one in X(N), and
      !  saving the other for later.
      !
      ELSE

         x_hi_index = x_hi_index - 1

         m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

         CALL Random_number ( harvest = r(1:2*m) )

         x(x_lo_index:x_hi_index-1:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-3:2) ) ) &
                                    & *  COS (   2._wp * r8_pi * r(2:2*m-2:2) )

         x(x_lo_index+1:x_hi_index:2) = SQRT ( - 2._wp * LOG ( r(1:2*m-3:2) ) )  &
                                    & * SIN  (   2._wp * r8_pi * r(2:2*m-2:2) )

         x(n) = SQRT ( - 2._wp * LOG ( r(2*m-1) ) ) * COS ( 2._wp * r8_pi * r(2*m) )

      END IF

      IF ( allones ) THEN
        IF (lwp .AND. mod(kt,1) .eq. 0) THEN
            print *, '   Gaussian RV not used: modes all mutliplied by ones' 
        END IF
        x(:) = 1._wp
      END IF

      RETURN
   END SUBROUTINE 

END MODULE tlurnd

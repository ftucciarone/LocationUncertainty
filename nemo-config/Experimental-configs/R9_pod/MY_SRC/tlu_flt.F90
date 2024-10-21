MODULE tluflt
   !!===========================================================================
   !!                       ***  MODULE  tlusvd  ***
   !! Transport under Location Uncertainty: stochastic noise SVD
   !!                       noise generation based on pseudo-observations and SVD.
   !!===========================================================================
   !! History : 0.0  !  2017-..  (P. DERIAN)  Original code
   !!
   !! [TODO]    - read namelist from _cfg, not only _ref
   !!           - write error messages to appropriate unit
   !!           - use NEMO's work arrays?
   !!---------------------------------------------------------------------------
   ! [mod_dep]
   USE par_kind           ! data types defined in par_kind module
   USE in_out_manager     ! I/O manager
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PUBLIC :: filter
   !
   TYPE, PUBLIC :: filter
      CHARACTER(LEN=3)                             :: kind
      INTEGER(i4)                                  :: hsiz   !> @public  Gaussian filter half-size
      INTEGER(i4)                                  :: wdth   !> @private Gaussian filter width  
      INTEGER(i4)                                  :: cntr   !> @private Gaussian kernel center
      INTEGER(i4),     ALLOCATABLE, DIMENSION(:)   :: kidx   !> @private Gaussian kernel indexes
      REAL(wp),        ALLOCATABLE, DIMENSION(:,:) :: kern   !> @public  Gaussian kernel weights
      REAL(wp)                                     :: wtot   !> @private Gaussian kernel total weight
      REAL(wp)                                     :: alpha  !> @private Gaussian kernel isd parameters
      REAL(wp)                                     :: beta   !> @private Gaussian kernel isd parameters
   CONTAINS
      PROCEDURE, PASS(self) :: kern_init
      PROCEDURE, PASS(self) :: kern_apply_1fld
      PROCEDURE, PASS(self) :: kern_apply_2fld
   END TYPE
   !
CONTAINS

   SUBROUTINE kern_init( self, kind, hsiz )
      CLASS(filter),     INTENT(inout) :: self
      CHARACTER(LEN=3),  INTENT(in   ) :: kind
      INTEGER(i4),       INTENT(in   ) :: hsiz   !> @public  Gaussian filter half-size

      INTEGER                      :: i

      self % kind = kind
      self % hsiz = hsiz
      self % wdth = ( self % hsiz * 2 ) + 1
      self % cntr =   self % hsiz       + 1
      self % kidx = [( i, i = - self % hsiz, self % hsiz )]

      ALLOCATE ( self % kern ( self % wdth, self % wdth ) )

      SELECT CASE( self % kind )
         CASE("box")
            CALL box_kern ( self % wdth, self % kern, self % wtot )
         CASE("gau")
            CALL gau_kern ( self % cntr, self % hsiz, self % wdth, self % kern, self % alpha, self % beta, self % wtot)
      END SELECT

      IF (lwp) THEN 
         WRITE(numout,*) '   kern_init : Initialize filter kernel '
         WRITE(numout,*) '   ~~~~~~~~~ '
         WRITE(numout,*) '         Type of kernel (gaussian or box)     kind = ', self % kind
         WRITE(numout,*) '                      half-size of kernel     hsiz = ', self % hsiz
         WRITE(numout,*) '                          width of kernel     wdth = ', self % wdth
         WRITE(numout,*) '                         centre of kernel     cntr = ', self % cntr
         WRITE(numout,*) '                   total weight of kernel     wtot = ', self % wtot
      END IF

   END SUBROUTINE kern_init

   SUBROUTINE kern_apply_2fld( self, fldin_u, fldin_v, fldout_u, fldout_v )
      USE tluMPI
      CLASS(filter),                   INTENT(in   ) :: self
      REAL(wp),      DIMENSION(:,:,:), INTENT(in   ) ::  fldin_u
      REAL(wp),      DIMENSION(:,:,:), INTENT(in   ) ::  fldin_v
      REAL(wp),      DIMENSION(:,:,:), INTENT(  out) :: fldout_u
      REAL(wp),      DIMENSION(:,:,:), INTENT(  out) :: fldout_v
      !
      INTEGER                                        :: i, j
      INTEGER                                        :: isiz, jsiz, ksiz
      !
      isiz = SIZE(fldout_u, 1)
      jsiz = SIZE(fldout_u, 2)
      ksiz = SIZE(fldout_u, 3)
      !
      ! Initialize 
      !
      fldout_u = 0._wp
      fldout_v = 0._wp
      !
      DO i = 1, self % wdth 
         !
         DO j = 1, self % wdth
            !
            fldout_u(:,:,:) = fldout_u(:,:,:)  +  self % kern( self % cntr + self % kidx (i),        self % cntr + self % kidx (j) ) &
            &                                  *      fldin_u( self % cntr + self % kidx (i): isiz + self % hsiz + self % kidx (i) , &
            &                                                  self % cntr + self % kidx (j): jsiz + self % hsiz + self % kidx (j) , 1:ksiz ) 
            !
            fldout_v(:,:,:) = fldout_v(:,:,:)  +  self % kern( self % cntr + self % kidx (i),        self % cntr + self % kidx (j) ) &
            &                                  *      fldin_v( self % cntr + self % kidx (i): isiz + self % hsiz + self % kidx (i) , &
            &                                                  self % cntr + self % kidx (j): jsiz + self % hsiz + self % kidx (j) , 1:ksiz )            !
         END DO
         !
      END DO
      !
   END SUBROUTINE kern_apply_2fld
   ! [kern_apply_2fld]

   SUBROUTINE kern_apply_1fld( self, fldin, fldout )
      USE tluMPI
      CLASS(filter),                   INTENT(in   ) :: self
      REAL(wp),      DIMENSION(:,:,:), INTENT(in   ) :: fldin
      REAL(wp),      DIMENSION(:,:,:), INTENT(  out) :: fldout
      !
      INTEGER                                        :: ki, kj ! index in kernel
      INTEGER                                        :: fi, li ! first and last index in i
      INTEGER                                        :: fj, lj ! first and last index in j
      INTEGER                                        :: i, j, ofs
      INTEGER                                        :: isiz, jsiz, ksiz
      !
      isiz = SIZE(fldout, 1)
      jsiz = SIZE(fldout, 2)
      ksiz = SIZE(fldout, 3)
      !
      ofs = 1 + ( SIZE(fldin, 1) - isiz ) / 2
      !
      ! Initialize 
      !
      fldout = 0._wp
      !
      DO i = 1, self % wdth 
         !
         ki = self % cntr + self % kidx (i)
         fi =         ofs + self % kidx (i)
         li =         ofs + self % kidx (i) + isiz - 1 
         !
         DO j = 1, self % wdth
            !
            kj = self % cntr + self % kidx (j)
            fj =         ofs + self % kidx (j)
            lj =         ofs + self % kidx (j) + jsiz - 1
            !
            fldout(:,:,:) = fldout(:,:,:) + self % kern( ki, kj ) * fldin( fi:li, fj:lj, 1:ksiz )
            !
         END DO
         !
      END DO
      !
   END SUBROUTINE kern_apply_1fld
   ! [kern_apply_1fld]


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
   SUBROUTINE gau_kern (cntr, hsiz, wdth, kerl, alpha, beta, wt)
      ! I/O arguments
      INTEGER,         INTENT(in   ) :: cntr
      INTEGER,         INTENT(in   ) :: hsiz
      INTEGER,         INTENT(in   ) :: wdth
      REAL(KIND = wp), INTENT(  out) :: kerl(wdth,wdth)
      REAL(KIND = wp), INTENT(  out) :: alpha
      REAL(KIND = wp), INTENT(  out) :: beta
      REAL(KIND = wp), INTENT(  out) :: wt
      ! Internal
      INTEGER                        :: wd
      INTEGER                        :: i,j
      REAL(KIND = wp)                :: pi = 4._wp * atan(1._wp) 

      wd = 2 * hsiz
      wt = 0._wp
      alpha = 0._wp
      beta = 0._wp

      IF ( wd .gt. 0 ) THEN
         DO i = 1, wd+1  
            DO j = 1, wd+1
!              Gaussian convolution kernel        
               kerl(i,j) =      (   6._wp / ( pi * (wd**2) ) ) *      &
               &             exp( - 6._wp * ( (i-cntr)**2 +           &
               &                              (j-cntr)**2 ) / (wd**2) ) 
!              Total weight of kernel          
               wt = wt + kerl(i,j)
               
               alpha = alpha + kerl(i,j)**2
               beta = beta + ( kerl(i,j)**2 ) * ( 144._wp * ( (i-cntr)**2 +           &
               &                                              (j-cntr)**2 ) / (wd**4) ) 
            END DO
         END DO
      ELSE
         kerl(1,1) = 1._wp
         alpha = 1._wp
         beta = 1._wp         
         wt = 1._wp
      END IF

   END SUBROUTINE gau_kern

   SUBROUTINE box_kern ( wdth, kerl, wt)
      ! I/O arguments
      INTEGER,         INTENT(in   ) :: wdth
      REAL(KIND = wp), INTENT(  out) :: kerl(wdth,wdth)
      REAL(KIND = wp), INTENT(  out) :: wt

      wt = 0._wp
      
      IF ( wdth .gt. 1 ) THEN
         kerl = 1._wp / ( wdth )**2
         wt = SUM(kerl)
      ELSE
         kerl(1,1) = 1._wp
         wt = 1._wp
      END IF

   END SUBROUTINE box_kern

END MODULE tluflt

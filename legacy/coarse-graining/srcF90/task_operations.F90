!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: coarse_subs
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
!> @details       This module includes all subroutines of filtering for performing
!!                the program coarse_grain.
!!
!!------------------------------------------------------------------------------
MODULE task_operations

!    Modules
   USE class_precision

   IMPLICIT NONE
   PRIVATE


!  Array dimensions    
!  ----------------
   INTEGER, PUBLIC, SAVE :: nxto        ! numbers of oceanic inner gridcells in x (without boundary points)
   INTEGER, PUBLIC, SAVE :: nyto        ! numbers of oceanic inner gridcells in y (without boundary points)
   INTEGER, PUBLIC, SAVE :: nso         ! factor of both nxto and nyto (Ex. nso = 9 for data from R27 to R3) 
   INTEGER, PUBLIC, SAVE :: co1, co2    !
   
   REAL(KIND = wp), PUBLIC, SAVE :: km1, km2    ! 

!  Derived grid parameters 
!  -------------------------------------- 
   INTEGER, PUBLIC, SAVE :: nxpo  ! = nxto + 2    oceanic complete grid with boundaries in x 
   INTEGER, PUBLIC, SAVE :: nypo  ! = nyto + 2    oceanic complete grid with boundaries in y
   INTEGER, PUBLIC, SAVE :: nxtoc ! = nxto/nso    coarse inner gridcells in x
   INTEGER, PUBLIC, SAVE :: nytoc ! = nyto/nso    coarse inner gridcells in y
   INTEGER, PUBLIC, SAVE :: nxpoc ! = nxtoc + 2   complete coarse grid points in x
   INTEGER, PUBLIC, SAVE :: nypoc ! = nytoc + 2   complete coarse grid points in y

   LOGICAL, PARAMETER :: rescale = .TRUE.

!  Subroutines      
   PUBLIC :: gaussf  
   PUBLIC :: smooth    
   PUBLIC :: filt      
   PUBLIC :: DifferenceOfGaussians       
   PUBLIC :: subsample    
   PUBLIC :: uv2w      
   PUBLIC :: bcn2btc      
   PUBLIC :: eke
!   PUBLIC :: tau2v

CONTAINS

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE gaussf (kerl, wt, wd)

!     Construct 2D Gaussian filter kernel kerl of width wd. 
!     The total weight wt is also calculated.
!     [ See https://en.wikipedia.org/wiki/Filter_(large_eddy_simulation) ]

      IMPLICIT NONE    
      
!     I/O arguments
      INTEGER,         INTENT(IN ) :: wd
      REAL(KIND = wp), INTENT(OUT) :: kerl(wd+1,wd+1),wt

!     Local variables      
      INTEGER                      :: i,j,c
      REAL(KIND = wp)              :: pi

      wt = 0._wp
      pi = 4._wp * atan(1._wp) !3.14159265D0
      c  = wd/2 + 1 

      IF ( wd .gt. 0 ) THEN
         DO i = 1, wd+1  
            DO j = 1, wd+1
!              Gaussian convolution kernel        
               kerl(i,j) = ( 6._wp/(pi*(wd**2)) ) *  &
         &                exp( - 6._wp*( (i-c)**2 + (j-c)**2 ) / (wd**2) ) 
!              Total weight of kernel          
               wt = wt + kerl(i,j)
            END DO
         END DO
      ELSE
         kerl(1,1) = 1._wp
         wt = 1._wp
      END IF

   END SUBROUTINE gaussf

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE smooth (pf, p, fk, wt, wd, idendx, idendy)

!     Coarse-graining of eddy-resolving snapshot 'p' using the convolution
!     kernel fk of total weight 'wt' and of width 'wd'.

      IMPLICIT NONE

!     I/O arguments
      INTEGER,         INTENT(IN)                          :: wd
      INTEGER,         INTENT(IN)                          :: idendx,idendy 
      REAL(KIND = wp), INTENT(IN)                          :: wt
      REAL(KIND = wp), INTENT(IN)                          :: fk(wd+1,wd+1)
      REAL(KIND = wp), INTENT(IN)                          :: p(nxpo,nypo)
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpoc,nypoc) :: pf
!
!     fk is filter kernel of total weight wt and width wd
!     p is eddy-resolving (high-resolution) velocity component of no-slip 
!     boundary condition
!     pf and pr are corase-grained and residual data (only computed for
!     coarse gird points)

!     Local variables
      INTEGER nh,i,j,ii,jj,il,jl,ig,jg
      REAL(KIND = wp) psum,wsum


      nh = wd/2 + 1 !! size of half-kernel
      pf = 0._wp !! null boundary (no-slip)

      !$OMP PARALLEL PRIVATE (i,j,ii,jj,il,jl,ig,jg,psum,wsum)
      !$OMP&         SHARED  (p,pf,pr,fk,wt,wd,nh,idendx,idendy)

!     Filtering inner points     
!     ----------------------
      !$OMP DO SCHEDULE(STATIC)
      DO i=2,idendx !! coarse-grid index (without boundaries)
        ii = 1 + (i-1)*nso !! corresp. index on fine-grid
        DO j=2,idendy
          jj = 1 + (j-1)*nso
          psum = 0._wp !! convoluted state
          wsum = 0._wp !! accounted weight
          DO il=1,wd+1 !! local index of kernel points
            ig = ii + il - nh !! global index corresp. to fine-grid
            IF ( (ig.lt.1) .or. (ig.gt.nxpo) ) THEN
                CYCLE !! skip iter. if outside of fine-grid domain
            ENDIF
            DO jl=1,wd+1
              jg = jj + jl - nh
              IF ( (jg.lt.1) .or. (jg.gt.nypo) ) THEN
                CYCLE 
              ENDIF  
!             Convolution of p by fk
              psum = psum + fk(il,jl)*p(ig,jg)
              wsum = wsum + fk(il,jl)            
            ENDDO
          ENDDO
!         Rescale pf by wt          
          IF (rescale) THEN          
            pf(i,j) = psum*wt/wsum
          ELSE
            pf(i,j) = psum
          ENDIF
        ENDDO
      ENDDO
      !$OMP END DO 
      !$OMP END PARALLEL
   !
   END SUBROUTINE smooth


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE filt (pf, pr, p, fk, wt, wd, idendx, idendy)

!     Coarse-graining of eddy-resolving snapshot 'p' using the convolution
!     kernel fk of total weight 'wt' and of width 'wd'.

      IMPLICIT NONE

!     I/O arguments
      INTEGER,         INTENT(IN)                          :: wd
      INTEGER,         INTENT(IN)                          :: idendx,idendy 
      REAL(KIND = wp), INTENT(IN)                          :: wt
      REAL(KIND = wp), INTENT(IN)                          :: fk(wd+1,wd+1)
      REAL(KIND = wp), INTENT(IN)                          :: p(nxpo,nypo)
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpoc,nypoc) :: pf,pr
!
!     fk is filter kernel of total weight wt and width wd
!     p is eddy-resolving (high-resolution) velocity component of no-slip 
!     boundary condition
!     pf and pr are corase-grained and residual data (only computed for
!     coarse gird points)

!     Local variables
      INTEGER nh,i,j,ii,jj,il,jl,ig,jg
      REAL(KIND = wp) psum,wsum


      nh = wd/2 + 1 !! size of half-kernel
      pf = 0._wp !! null boundary (no-slip)

      !$OMP PARALLEL PRIVATE (i,j,ii,jj,il,jl,ig,jg,psum,wsum)
      !$OMP&         SHARED  (p,pf,pr,fk,wt,wd,nh,idendx,idendy)

!     Filtering inner points     
!     ----------------------
      !$OMP DO SCHEDULE(STATIC)
      DO i=2,idendx !! coarse-grid index (without boundaries)
        ii = 1 + (i-1)*nso !! corresp. index on fine-grid
        DO j=2,idendy
          jj = 1 + (j-1)*nso
          psum = 0._wp !! convoluted state
          wsum = 0._wp !! accounted weight
          DO il=1,wd+1 !! local index of kernel points
            ig = ii + il - nh !! global index corresp. to fine-grid
            IF ( (ig.lt.1) .or. (ig.gt.nxpo) ) THEN
                CYCLE !! skip iter. if outside of fine-grid domain
            ENDIF
            DO jl=1,wd+1
              jg = jj + jl - nh
              IF ( (jg.lt.1) .or. (jg.gt.nypo) ) THEN
                CYCLE 
              ENDIF  
!             Convolution of p by fk
              psum = psum + fk(il,jl)*p(ig,jg)
              wsum = wsum + fk(il,jl)            
            ENDDO
          ENDDO
!         Rescale pf by wt          
          IF (rescale) THEN          
            pf(i,j) = psum*wt/wsum
          ELSE
            pf(i,j) = psum
          ENDIF
!         Substract pr
          pr(i,j) = p(ii,jj) - pf(i,j)
        ENDDO
      ENDDO
      !$OMP END DO 
      !$OMP END PARALLEL
   !
   END SUBROUTINE filt

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE subsample (ps, p, idendx, idendy)

!     Coarse-graining of eddy-resolving snapshot 'p' using the convolution
!     kernel fk of total weight 'wt' and of width 'wd'.

      IMPLICIT NONE

!     I/O arguments
      INTEGER,         INTENT(IN)                          :: idendx,idendy 
      REAL(KIND = wp), INTENT(IN)                          :: p(nxpo,nypo)
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpoc,nypoc) :: ps
!
!     fk is filter kernel of total weight wt and width wd
!     p is eddy-resolving (high-resolution) velocity component of no-slip 
!     boundary condition
!     pf and pr are corase-grained and residual data (only computed for
!     coarse gird points)

!     Local variables
      INTEGER nh,i,j,ii,jj,il,jl,ig,jg
      REAL(KIND = wp) psum,wsum


      ps = 0._wp !! null boundary (no-slip)

      !$OMP PARALLEL PRIVATE (i,j,ii,jj,il,jl,ig,jg,psum,wsum)
      !$OMP&         SHARED  (p,pf,pr,fk,wt,wd,nh,idendx,idendy)

!     Filtering inner points     
!     ----------------------
      !$OMP DO SCHEDULE(STATIC)
      DO i=2,idendx !! coarse-grid index (without boundaries)
        ii = 1 + (i-1)*nso !! corresp. index on fine-grid
        DO j=2,idendy
          jj = 1 + (j-1)*nso
!         Substract pr
          ps(i,j) = p(ii,jj)
        ENDDO
      ENDDO
      !$OMP END DO 
      !$OMP END PARALLEL
   !
   END SUBROUTINE subsample


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE DifferenceOfGaussians(band, pf1, pf2, p, fk1, wt1, wd1, &
                                   &                   fk2, wt2, wd2, idendx, idendy)

!     Coarse-graining of eddy-resolving snapshot 'p' using the convolution
!     kernel fk of total weight 'wt' and of width 'wd'.

      IMPLICIT NONE

!     I/O arguments
      INTEGER,         INTENT(IN)                          :: wd1, wd2
      INTEGER,         INTENT(IN)                          :: idendx,idendy 
      REAL(KIND = wp), INTENT(IN)                          :: wt1, wt2
      REAL(KIND = wp), INTENT(IN)                          :: fk1(wd1+1,wd1+1)
      REAL(KIND = wp), INTENT(IN)                          :: fk2(wd2+1,wd2+1)
      REAL(KIND = wp), INTENT(IN)                          :: p(nxpo,nypo)
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpoc,nypoc) :: pf1,pf2
      REAL(KIND = wp), INTENT(OUT), DIMENSION(nxpoc,nypoc) :: band

!
!     fk is filter kernel of total weight wt and width wd
!     p is eddy-resolving (high-resolution) velocity component of no-slip 
!     boundary condition
!     pf and pr are corase-grained and residual data (only computed for
!     coarse gird points)

!     Local variables
      INTEGER nh1,nh2,i,j,ii,jj,il,jl,ig,jg
      REAL(KIND = wp) psum1, wsum1, psum2, wsum2


      nh1 = wd1/2 + 1 !! size of half-kernel
      nh2 = wd2/2 + 1 !! size of half-kernel

      pf1 = 0._wp  !! null boundary (no-slip)
      pf2 = 0._wp  !! null boundary (no-slip)
      band = 0._wp !! null boundary (no-slip)

      !$OMP PARALLEL PRIVATE (i,j,ii,jj,il,jl,ig,jg,psum1,psum2,wsum1,wsum2)
      !$OMP&         SHARED  (p,pf1,pf2,band,fk1,fk2,wt1,wt2,wd1,wd2,nh1,nh2,idendx,idendy)

!     Filtering inner points     
!     ----------------------
      !$OMP DO SCHEDULE(STATIC)
      DO i=2,idendx !! coarse-grid index (without boundaries)
        ii = 1 + (i-1)*nso !! corresp. index on fine-grid
        DO j=2,idendy
          jj = 1 + (j-1)*nso
          !
          ! First Gaussian
          psum1 = 0._wp !! convoluted state
          wsum1 = 0._wp !! accounted weight
          !
          DO il=1,wd1+1 !! local index of kernel points
            ig = ii + il - nh1 !! global index corresp. to fine-grid
            IF ( (ig.lt.1) .or. (ig.gt.nxpo) ) THEN
                CYCLE !! skip iter. if outside of fine-grid domain
            ENDIF
            DO jl=1,wd1+1
              jg = jj + jl - nh1
              IF ( (jg.lt.1) .or. (jg.gt.nypo) ) THEN
                CYCLE 
              ENDIF  
!             Convolution of p by fk
              psum1 = psum1 + fk1(il,jl)*p(ig,jg)
              wsum1 = wsum1 + fk1(il,jl)            
            ENDDO
          ENDDO
          !
          ! Second Gaussian
          psum2 = 0._wp
          wsum2 = 0._wp
          DO il=1,wd2+1 !! local index of kernel points
            ig = ii + il - nh2 !! global index corresp. to fine-grid
            IF ( (ig.lt.1) .or. (ig.gt.nxpo) ) THEN
                CYCLE !! skip iter. if outside of fine-grid domain
            ENDIF
            DO jl=1,wd2+1
              jg = jj + jl - nh2
              IF ( (jg.lt.1) .or. (jg.gt.nypo) ) THEN
                CYCLE 
              ENDIF  
!             Convolution of p by fk
              psum2 = psum2 + fk2(il,jl)*p(ig,jg)
              wsum2 = wsum2 + fk2(il,jl)            
            ENDDO
          ENDDO
!         Rescale pf by wt          
          IF (rescale) THEN          
            pf1(i,j) = psum1*wt1/wsum1
            pf2(i,j) = psum2*wt2/wsum2
          ELSE
            pf1(i,j) = psum1
            pf2(i,j) = psum2
          ENDIF
!         Substract pr
          band(i,j) = pf1(i,j) - pf2(i,j)
        ENDDO
      ENDDO
      !$OMP END DO 
      !$OMP END PARALLEL
   !
   END SUBROUTINE DifferenceOfGaussians


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE uv2w (w, u, v, e1t, e2t, e1v, e2u, e3u, e3v, e3t, nx, ny, nlo)

      IMPLICIT NONE    
      
!     I/O arguments
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(nx,ny,nlo) :: u,v
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(nx,ny)     :: e1t,e2t
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(nx,ny)     :: e2u,e1v
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(nx,ny,nlo) :: e3u,e3v
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(nx,ny,nlo) :: e3t
      INTEGER,         INTENT (IN   )                       :: nx, ny,nlo

      REAL(KIND = wp), INTENT (  OUT), DIMENSION(nx,ny,nlo) :: w
!     Local variables      
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)       :: hdivn
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)       :: tmask
      INTEGER                                              :: ierr
      INTEGER                                              :: ji,jj,jk


      ALLOCATE( hdivn(nx,ny,nlo) ,  &
              & tmask(nx,ny,nlo) ,  STAT=ierr )
      IF( ierr .ne. 0 )  print *, 'Allocation error in uv2w'

      w = 0._wp !! null boundaries
      hdivn = 0._wp
      tmask = 0._wp


      DO jk = 1, nlo-1                                      !==  Horizontal divergence  ==!
         DO jj = 2, ny-1
            DO ji = 2, nx-1   ! vector opt.
               hdivn(ji,jj,jk) =  ( e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk) * u(ji  ,jj  ,jk)       &
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


   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE bcn2btc (vel, vel_b, e3_vel, nx, ny, nlo)

!     Given horizontal velocities u and v, compute the vertical
!     component w using divergence-free condition.


      IMPLICIT NONE    
      
!     I/O arguments
      INTEGER,         INTENT (IN   )                       :: nx, ny, nlo
      REAL(KIND = wp), INTENT(in   ), DIMENSION(nx,ny,nlo) :: vel
      REAL(KIND = wp), INTENT(in   ), DIMENSION(nx,ny,nlo) :: e3_vel
      REAL(KIND = wp), INTENT(  out), DIMENSION(nx,ny)     :: vel_b
!     Local variables     
      REAL(KIND = wp),                DIMENSION(nx,ny)     :: sum_e3 
      INTEGER                                              :: ierr
      INTEGER                                              :: ji,jj,jk

      vel_b = 0._wp !! null boundaries
      sum_e3 = 0._wp
      sum_e3 = SUM(e3_vel, DIM=3)
      vel_b = SUM(vel*e3_vel, DIM=3) / sum_e3

   END SUBROUTINE bcn2btc






   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F.L. Tucciarone 
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
   SUBROUTINE eke (KE, uin, vin, win, e1u, e1v, e1t, e2u, e2v, e2t, e3u, e3v, e3w)

      IMPLICIT NONE    
      
!     I/O arguments
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(:,:,:) :: uin, vin, win
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(:,:)   :: e1u, e1v, e1t
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(:,:)   :: e2u, e2v, e2t
      REAL(KIND = wp), INTENT (IN   ), DIMENSION(:,:,:) :: e3u, e3v, e3w

      REAL(KIND = wp), INTENT (  OUT)                   :: KE

      KE = 0._wp !! null boundaries

      KE = SUM( ( uin**2 ) * (spread(e1u,3,SIZE(e3u,3))*spread(e2u,3,SIZE(e3u,3))*e3u) + &
              & ( vin**2 ) * (spread(e1v,3,SIZE(e3v,3))*spread(e2v,3,SIZE(e3v,3))*e3v) + &
              & ( win**2 ) * (spread(e1t,3,SIZE(e3w,3))*spread(e2t,3,SIZE(e3w,3))*e3w) ) 

      KE = SQRT( KE ) 

   END SUBROUTINE eke







END MODULE task_operations    




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
MODULE coarse_subs

!    Modules
   USE class_precision

   IMPLICIT NONE
   PRIVATE


!  Subroutines      
   PUBLIC :: gaussf, filt, uv2w, tau2v

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

      DO i = 1, wd+1  
         DO j = 1, wd+1
!           Gaussian convolution kernel        
            kerl(i,j) = ( 6._wp/(pi*(wd**2)) ) *  &
     &                exp( - 6._wp*( (i-c)**2 + (j-c)**2 ) / (wd**2) ) 
!           Total weight of kernel          
            wt = wt + kerl(i,j)
         END DO
      END DO
   END SUBROUTINE gaussf

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
   SUBROUTINE filt (pf, pr, p, fk, wt, wd, idendx, idendy)

!     Coarse-graining of eddy-resolving snapshot 'p' using the convolution
!     kernel fk of total weight 'wt' and of width 'wd'.

!     Modules          
      USE param

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
   SUBROUTINE uv2w (w, u, v)

!     Given horizontal velocities u and v, compute the vertical
!     component w using divergence-free condition.

!     Modules
      USE param
      USE state

      IMPLICIT NONE    
      
!     I/O arguments
      REAL(KIND = wp), INTENT(  in ), DIMENSION(nxpoc,nypoc,nlo) :: u,v
      REAL(KIND = wp), INTENT( out ), DIMENSION(nxpoc,nypoc,nlo) :: w
!     Local variables      
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)        :: hdivn
      REAL(KIND = wp), ALLOCATABLE, DIMENSION(:,:,:)        :: tmask
      INTEGER                                               :: ierr
      INTEGER                                               :: ji,jj,jk


      ALLOCATE( hdivn(nxpoc,nypoc,nlo) ,  &
              & tmask(nxpoc,nypoc,nlo) ,  STAT=ierr )
      IF( ierr .ne. 0 )  print *, 'Allocation error in uv2w'

      w = 0._wp !! null boundaries
      hdivn = 0._wp
      tmask = 1._wp
    
      DO jk = 1, nlo-1                                      !==  Horizontal divergence  ==!
         DO jj = 2, nypoc-1
            DO ji = 2, nxpoc-1   ! vector opt.
               hdivn(ji,jj,jk) = (  e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk) * u(ji  ,jj  ,jk)      &
                  &               - e2u(ji-1,jj  ) * e3u(ji-1,jj  ,jk) * u(ji-1,jj  ,jk)      &
                  &               + e1v(ji  ,jj  ) * e3v(ji  ,jj  ,jk) * v(ji  ,jj  ,jk)      &
                  &               - e1v(ji  ,jj-1) * e3v(ji  ,jj-1,jk) * v(ji  ,jj-1,jk)  )   & 
                  &             / ( e3t(ji  ,jj  ,jk) * e1t(ji,jj) * e2t(ji,jj) )
!                  &            * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)

!               IF ( ( b_lev(ji,jj) <= jk ) .or. ( jk <= t_lev(ji,jj) ) ) THEN
!                  tmask(ji,jj,jk) = 1._wp
!               END IF 

            END DO  
         END DO  
      END DO

      DO jk = nlo-1, 1, -1                       ! integrate from the bottom the hor. divergence
         ! computation of w
         w(:,:,jk) = w(:,:,jk+1) - ( e3t(:,:,jk) * hdivn(:,:,jk) ) * tmask(:,:,jk)
      END DO

   END SUBROUTINE uv2w

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> L.Li, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief    Computes the shear velocity of the wind
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
   SUBROUTINE tau2v (tau, v)

!     Modules
      USE param
      USE state
     
      IMPLICIT NONE     

!     I/O arguments
      REAL(KIND = wp), INTENT(  in ), DIMENSION(nxpo,nypo) :: tau
      REAL(KIND = wp), INTENT( out ), DIMENSION(nxpo,nypo) :: v
!     Local variables      
      REAL(KIND = wp), PARAMETER                             :: rhoair = 1.225_wp ! kg/m^3
      INTEGER :: ji, jj

      v(:,:) =  tau(:,:) / rhoair

      DO jj = 1, nypo
         DO ji = 1, nxpo
            v(ji,jj) = ( v(ji,jj) / ABS( v(ji,jj) ) ) * & ! Sign 
                     &   SQRT( ABS(v(ji,jj)) )   
         END DO
      END DO

   END SUBROUTINE tau2v

END MODULE coarse_subs    


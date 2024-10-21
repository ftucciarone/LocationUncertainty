MODULE tlugss
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
   USE iom                ! I/O manager library
   USE par_oce            ! ocean dynamics parameters
   USE oce                ! ocean dynamics and active tracers
   USE dom_oce            ! ocean space and time domain
   USE lib_mpp            ! MPP library
   USE timing             ! Timing
   USE lbclnk             ! ocean lateral boundary conditions (or mpp link)
   USE tlu
   USE tluflt
   !   USE tludgns
   ! [mod_dep]
   !
   IMPLICIT NONE          ! turn off implicit variable declaration
   PRIVATE                ! and make stuff private by default
   !
   PUBLIC tlu_init_gss    ! called by tlu_init in tlu.f90
   PUBLIC gen_white_noise ! called by tlu_fields in tlu_noi.F90
   !
   INCLUDE 'mpif.h'
   !
   ! [ ]
   INTEGER(i4), PUBLIC                                      :: nn_tlu_hsiz = 2 !> @public size of halo
   INTEGER(i4), PUBLIC                                      :: nn_tlu_hbox = 2 !> @public 
   ! [ ]   
   !
   ! [tlu_spnoise]
   REAL(wp)                                                 :: amp = 0.125_wp    !> @private amplitude in m/s
   !
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: u_wht           !> @public U-velocity noise
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: v_wht           !> @public V-velocity noise
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: w_wht           !> @public W-velocity noise
   !
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: u_gsv           !> @public U-velocity Girsenov drift
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: v_gsv           !> @public V-velocity Girsenov drift
   REAL(wp),    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: w_gsv           !> @public W-velocity Girsenov drift
   ! [tlu_spnoise]
   !
   ! [tlu_filt]
   INTEGER(i4)                                              :: ctr             !> @private kernel center
   INTEGER(i4),         ALLOCATABLE,       DIMENSION(:)     :: kidx            !> @private kernel indexes
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: kerf            !> @private kernel weights
   REAL(wp)                                                 :: wtot            !> @private kernel total weight
   ! [tlu_filt] 
   !
   INTEGER(i4)                                              :: jpi_ext, jpj_ext  !> @private extra halo max dimension
   INTEGER(i4)                                              :: fih, fjh          !> @private first index after extra halo
   INTEGER(i4)                                              :: lih, ljh          !> @private last index begore extra halo
   !
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: urnd, vrnd      !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wext            !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk1, wrk2      !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:,:) :: wrk3, wrk4      !> @private working vectors
   REAL(wp),            ALLOCATABLE,       DIMENSION(:,:)   :: wrk5, wrk6      !> @private working vectors
   !
   TYPE(filter)                                             :: gaussian
   TYPE(filter)                                             :: box_3by3
   !
CONTAINS

   SUBROUTINE tlu_init_gss
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_init_svd  ***
      !!
      !! ** Purpose :   Initialize module variables.
      !!
      !! ** Method  :   Read namelist, compute other parameters and allocate arrays.
      !!
      !! ** Action  :   Arrays allocated.
      !!
      !!------------------------------------------------------------------------
      INTEGER(i4) :: ios, ierr(9), chk   ! namelist output, allocation statistics
      INTEGER     :: i
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_gss_init : TLU, Gaussian noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '     half-size of Gaussian kernel      nn_tlu_hsiz = ', nn_tlu_hsiz
      END IF
      !
      ! Allocate Internal variables
      !
      ierr = 0
      !
      ! Initialize filter
      !
      CALL gaussian%kern_init( "gau", nn_tlu_hsiz )
      CALL box_3by3%kern_init( "box", 1 )
      !
      ! Define extended dimensions
      !
      jpi_ext = jpi + 2*nn_tlu_hsiz
      jpj_ext = jpj + 2*nn_tlu_hsiz
      fih = nn_tlu_hsiz + 1
      fjh = nn_tlu_hsiz + 1
      lih = fih + jpi - 1
      ljh = fjh + jpj - 1
      !
      ! Allocate white noise field and girsanov drift
      !
      ALLOCATE( u_wht( jpi, jpj, jpk ), &
      &         v_wht( jpi, jpj, jpk ), &
      &         w_wht( jpi, jpj, jpk ), &
      &         u_gsv( jpi, jpj, jpk ), &
      &         v_gsv( jpi, jpj, jpk ), &
      &         w_gsv( jpi, jpj, jpk ), stat=ierr(2) )
      !
      ! Allocate Internal variables 
      !
      ALLOCATE( urnd(             jpi_ext,             jpj_ext, jpk ), &
      &         vrnd(             jpi_ext,             jpj_ext, jpk ), &
      &         wext( jpi + 2*nn_tlu_hbox, jpj + 2*nn_tlu_hbox, jpk ), stat=ierr(3) )
            !
      ALLOCATE( wrk1( jpi, jpj, jpk ), &
      &         wrk2( jpi, jpj, jpk ), &
      &         wrk3( jpi, jpj, jpk ), &
      &         wrk4( jpi, jpj, jpk ), &
      &         wrk5( jpi, jpj      ), &
      &         wrk6( jpi, jpj      ), stat=ierr(4) )   ! the singular values
      !
      IF (SUM(ierr)/=0) THEN   ! check allocation success (0=OK)
         IF(lwp) WRITE(numout,*) ' tlu_pod_init(): allocation failed = ', ierr  
         !CALL ctl_stop( ctmp1 )   ! [TODO] enable
         STOP
      END IF
      !     
   END SUBROUTINE tlu_init_gss


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE gen_white_noise ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2023 -  ( F. TUCCIARONE   )  Original code and documentation
   !! 
   !> @brief         
   !!                
   !!                
   !!
   !> @details
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!                
   !!             + dt_delay 
   !! 
   !> @param[in]     
   !> @param[in]     
   !> @param[in]     
   !> @param[inout]  
   !! 
   !! @result        
   !!
   !! @warning      
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this gen_white_noise
   ! [gen_white_noise]
   SUBROUTINE gen_white_noise( kt )
      USE tluMPI
      USE tlurnd
      INTEGER,                                  INTENT(in   ) :: kt                     ! ocean time-step index
      !
      INTEGER                                                 :: i, j, k, jj, jk
      !
      IF ( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'gen_white_noise : TLU, filtering of white noise field '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      END IF
      !
      ! Gaussian field generation
      !
      wrk1 = 0._wp
      wrk2 = 0._wp
      wrk3 = 0._wp
      wrk4 = 0._wp
      wrk5 = 1._wp
      wrk6 = 0._wp
      !
      DO jk = 1, jpk
         !
         CALL Standard_Gaussian_2D( jpi, jpj, wrk1(:,:,jk) ) !, allones = .False. )
         CALL Standard_Gaussian_2D( jpi, jpj, wrk2(:,:,jk) ) !, allones = .False. )
         !
      END DO
      !
      ! Collect extra information to compute spatial filtering
      !
      CALL MPI_collect( kt, wrk1, urnd, nn_tlu_hsiz, 'urnd' )
      CALL MPI_collect( kt, wrk2, vrnd, nn_tlu_hsiz, 'vrnd' )
      CALL MPI_collect( kt,   wn, wext, nn_tlu_hbox, 'wext' )
      !
      ! Generate the random modes by filtering
      !
      CALL gaussian%kern_apply_2fld(  urnd, vrnd, u_wht, v_wht )
      !
      ! Compute vertical velocity profile
      !
      CALL box_3by3%kern_apply_1fld( wext, wrk3 )
      wrk3 = wrk3**2
      !
      !                                ! ================
      DO jk = 1, jpk                   ! Horizontal slab
         !                             ! ================
         !
         wrk5(:,:) = wrk5(:,:) + e3w_n(:,:,jk) * wrk3(:,:,jk)
         !
         wrk6(:,:) = wrk6(:,:) + e3w_n(:,:,jk)
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      ! Compute the vertical profile to prescribe
      !
      wrk5(:,:) = e3t_n(:,:,jpk)
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         wrk5(:,:) = wrk5(:,:) + e3t_n(:,:,jk)
         wrk3(:,:,jk) = wrk3(:,:,jk+1) + e3t_n(:,:,jk)
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      wrk3 = wrk3 / spread(wrk5, 3,jpk)
!     if (lwp) print *, wrk3(2,2,:)
      u_wht = u_wht * amp * wrk3
      v_wht = v_wht * amp * wrk3

      CALL iom_put( "wbia_n", wrk3 )   
      !
      w_wht(:,:,:) = 0._wp
      CALL build_wnoi(w_wht, u_wht, v_wht)
      !
      ! Compute Girsenov Drift
      !
      u_gsv(:,:,:) = 0._wp
      v_gsv(:,:,:) = 0._wp
      w_gsv(:,:,:) = 0._wp
      !
      ! Transfer vertical fields across nodes
      !
      CALL lbc_lnk_multi( 'gen_white_noise', w_wht , 'T', -1., w_gsv , 'T', -1.   )
      !
   END SUBROUTINE gen_white_noise
   ! [gen_white_noise]



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE build_wnoi ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2021 -  ( F. TUCCIARONE   )  Original code <br> 
   !! 
   !> @brief         
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     kcmp: Define velocity type for calculation
   !> @param[in]     cdtype: Define velocity type
   !> @param[inout]  pcn, pca:  modified advection term for this component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this build_wnoi
   ! [build_wnoi]
   SUBROUTINE build_wnoi(wout, uin, vin)
      REAL(wp), DIMENSION(:,:,:), INTENT(in   )  :: uin
      REAL(wp), DIMENSION(:,:,:), INTENT(in   )  :: vin
      REAL(wp), DIMENSION(:,:,:), INTENT(  out)  :: wout
      !
      INTEGER                                    :: jk
      INTEGER                                    :: jei, jej
      !
      !!----------------------------------------------------------------------
      !
      jei = SIZE(uin, 1)
      jej = SIZE(uin, 2)
      !                                ! ================
      DO jk = jpkm1, 1, -1             ! Horizontal slab
         !                             ! ================
         !
         wout( 2:jei-1, 2:jej-1, jk) =     wout( 2:jei-1, 2:jej-1, jk + 1) -   &
         !
         !  dx( dy*dz*u_in )
         &                           (      e2u( 2:jei-1, 2:jej-1)             &
         &                           *    e3u_n( 2:jei-1, 2:jej-1, jk )        &
         &                           *      uin( 2:jei-1, 2:jej-1, jk )        &
         &                           -      e2u( 1:jei-2, 2:jej-1)             &
         &                           *    e3u_n( 1:jei-2, 2:jej-1, jk )        &
         &                           *      uin( 1:jei-2, 2:jej-1, jk )        &
         !  dy ( dx*dz*v_isd )
         &                           +      e1v( 2:jei-1, 2:jej-1)             &
         &                           *    e3v_n( 2:jei-1, 2:jej-1, jk )        &
         &                           *      vin( 2:jei-1, 2:jej-1, jk )        &
         &                           -      e1v( 2:jei-1, 1:jej-2)             &
         &                           *    e3v_n( 2:jei-1, 1:jej-2, jk )        &
         &                           *      vin( 2:jei-1, 1:jej-2, jk ) )      &
         !  .../dV
         &                           * r1_e1e2t( 2:jei-1, 2:jej-1 )
         !
         !                             ! ================
      END DO                           !   End of slab
      !                                ! ================
      !
      !
   END SUBROUTINE build_wnoi
   ! [build_wnoi]

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

END MODULE tlugss




!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tludgns
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!
!> @authors
!>                P. Derian, P. Chandramouli, F. Tucciarone
!!
!> @version
!!                2017 -  ( P. DERIAN       )  Original code <br>
!!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
!!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         
!! 
!! 
!! @par           Procedure specifics      
!> @details       
!!
!!------------------------------------------------------------------------------
MODULE tludgns
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
   USE iom                ! I/O manager library
   USE par_oce
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE lib_mpp          ! MPP library
   USE timing           ! Timing
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE tlu              ! Initialization of stochastic structures
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_dgns      ! Called by nemogcm.f90
   PUBLIC tlu_tvor      ! Called by step.f90
   PUBLIC int2d_ene3c   ! Called by step.f90
   PUBLIC int2d_3c      ! Called by step.f90
   PUBLIC eigensum      ! Called by step.F90
   PUBLIC cml_int2d_ene3c
   PUBLIC sumvars       ! Called by step.F90
   ! [public_sub]
   !
   INCLUDE 'mpif.h'
   !
   ! [tlu_StokesDrift]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   udva          !< @public Stokes drift: x component 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vdva          !< @public Stokes drift: y component
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   wdva          !< @public Stokes drift: z component
   ! [tlu_StokesDrift]
   !
   ! [tlu_Vorticity]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   tvor          !< @public Vorticity
   ! [tlu_Vorticity]
   !
   ! [statistics]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   usum,    usum2      !< @public 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vsum,    vsum2      !< @public 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   wsum,    wsum2      !< @public
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   tsum,    tsum2      !< @public 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   ssum,    ssum2      !< @public 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   vortsum, vortsum2   !< @public  
   INTEGER, PUBLIC :: stat_counter
   INTEGER, PUBLIC :: it0_counter
   ! [statistics]


   !
   !! * Substitutions [TODO] enable substitutions for NEMO
! #  include "domzgr_substitute.h90"
   !!---------------------------------------------------------------------------

CONTAINS

   !!---------------------------------------------------------------------------  
   !!> @author 
   !!> P. Derian, P. Chandramouli, F. Tucciarone 
   !!
   !!  DESCRIPTION: 
   !!  
   !!> @brief   Initialize module variables
   !!
   !!  @details
   !!  Reads namelist, compute other parameters and allocate arrays.
   !! 
   !!
   !! 
   !! 
   !> @param[in]     ln_tlu: logical switch for location uncertainty
   !! @param[in]     jpi: the first dimension of the spatial arrays
   !! @param[in]     jpj: the first dimension of the spatial arrays
   !! @param[in]     jpk: the first dimension of the spatial arrays
   !! @param[in]     numout: writing unit 
   !! @result        Arrays allocated
   !! @warning  
   !!  
   !! @note   
   !! @todo    
   !!               
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_dgns 
   ! [tlu_dgns]
   SUBROUTINE tlu_dgns
      INTEGER(i4) :: ierr   ! namelist output, allocation statistics
      !
      !!------------------------------------------------------------------------
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_dgns : Transport under Location Uncertainty, Diagnostic routine '
         WRITE(numout,*) '~~~~~~~~'

      END IF
      !
      !
      ! Allocate
      !------------------------------
      ALLOCATE( udva(jpi,jpj,jpk),      &
          &     vdva(jpi,jpj,jpk),      &
          &     wdva(jpi,jpj,jpk),      &
          &     tlu_wcorr(jpi,jpj,jpk), stat=ierr )     
           
      IF (ierr/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_dgns(): allocation failed = ', ierr   ! [TODO] should be ctmp1? instead of numout
         !CALL ctl_stop( ctmp1 )   !< @todo enable
         STOP
      END IF
      !
      !
   END SUBROUTINE tlu_dgns
   ! [tlu_dgns]



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tvor ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     kvor: Define typology of vorticity
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
   ! @snippet this tlu_tvor
   ! [tlu_tvor]
   SUBROUTINE tlu_tvor( kt, kvor, pun, pvn )
      USE dynvor
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pun, pvn    ! now velocities
      !
      INTEGER                                         ::   ji, jj, jk  ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj)                    ::   zwz         ! 2D workspace

      INTEGER, PARAMETER                              ::   vor_CV  = 1 
      INTEGER, PARAMETER                              ::   vor_RV  = 2
      INTEGER, PARAMETER                              ::   vor_CRV = 3
      INTEGER, PARAMETER                              ::   vor_PV  = 4

      REAL(wp), PARAMETER                             ::   r1_4 = 0.25_wp

      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dgns:tlu_tvor : vorticity term: energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'

         ALLOCATE( tvor(jpi,jpj,jpk) )   
      ENDIF

 

      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( vor_CV )                           !* Coriolis (planetary vorticity)
            zwz(:,:) = ff_f(:,:) 
         CASE ( vor_RV )                           !* Relative vorticity 
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) =  (  e2v(ji+1,jj) * pvn(ji+1,jj,jk) - e2v(ji,jj) * pvn(ji,jj,jk)      &
                     &           - e1u(ji,jj+1) * pun(ji,jj+1,jk) + e1u(ji,jj) * pun(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
               END DO
            END DO
         CASE ( vor_CRV )                           !* Coriolis + relative vorticity
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ff_f(ji,jj) + (  e2v(ji+1,jj) * pvn(ji+1,jj,jk) - e2v(ji,jj) * pvn(ji,jj,jk)      &
                     &                        - e1u(ji,jj+1) * pun(ji,jj+1,jk) + e1u(ji,jj) * pun(ji,jj,jk)  ) * r1_e1e2f(ji,jj)
               END DO
            END DO
         CASE ( vor_PV )                           !* (Coriolis + relative vorticity ) / h
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ff_f(ji,jj) + (  e2v(ji+1,jj) * pvn(ji+1,jj,jk) - e2v(ji,jj) * pvn(ji,jj,jk)      &
                     &                        - e1u(ji,jj+1) * pun(ji,jj+1,jk) + e1u(ji,jj) * pun(ji,jj,jk)  )   &
                     &                        * r1_e1e2f(ji,jj) / gdept_n(ji,jj,jk)
               END DO
            END DO
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','tlu_tvor: wrong value for kvor'  )
         END SELECT
         !
         IF( ln_dynvor_msk ) THEN          !==  mask/unmask vorticity ==!
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END DO
            END DO
         ENDIF

         IF( ln_sco ) THEN   ! Is this really necessary? FT
            zwz(:,:) = zwz(:,:) / e3f_n(:,:,jk)
         ENDIF
         !                                   !==  Interpolate the vorticity on the T-grid  =!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               tvor(ji,jj,jk) = r1_4 * (zwz(ji-1,jj  ) + zwz(ji  ,jj  )      &
                              &       + zwz(ji  ,jj-1) + zwz(ji  ,jj  ) )
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      CALL iom_put( "vort", tvor )  

   END SUBROUTINE tlu_tvor
   ! [tlu_tvor]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE int2d_ene3c ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes a 2-dimensional "energy-like" integral of a 
   !!                3-components field, i.e.
   !!
   !!                    \int\int u_{x}^{2} + u_{y}^{2} + u_{z}^{2} dxdy
   !!                
   !!                
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     cx: x-component 
   !> @param[in]     cy: y-component
   !> @param[in]     cz: z-component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this int2d_ene3c
   ! [int2d_ene3c]
   SUBROUTINE int2d_ene3c( cx, cy, cz, fieldname, norm_const )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::  cx, cy, cz    ! components 
      CHARACTER(len = * )         , INTENT(in   ) ::  fieldname     ! string to print
      REAL(wp)                                    ::  norm_const    ! components 
      !
      INTEGER                                     ::  ji, jj        ! dummy loop indices
      REAL(wp)                                    ::  sub_int       ! 2D workspace
      REAL(wp)                                    ::  dom_int       ! 2D workspace
      REAL(wp), DIMENSION(mppsize)                ::  zw1
      !
      INTEGER                                     ::  ierr, myID

      !!----------------------------------------------------------------------
      !
      sub_int = 0._wp
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            sub_int = sub_int + ( e1u(ji,jj) * e2u(ji,jj) * cx(ji,jj) * cx(ji,jj) &
                            &   + e1v(ji,jj) * e2v(ji,jj) * cy(ji,jj) * cy(ji,jj) &
                            &   + e1t(ji,jj) * e2t(ji,jj) * cz(ji,jj) * cz(ji,jj) ) * norm_const

         END DO
      END DO

      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLGATHER( sub_int, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)
      dom_int = SUM(zw1)

      call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
      IF ( myID == 0 ) THEN
         print '(A36, E16.7)', TRIM(fieldname)//' 2D domain integral value:', dom_int 
      END IF

   END SUBROUTINE int2d_ene3c
   ! [int2d_ene3c]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE cml_int2d_ene3c ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes a 2-dimensional "energy-like" integral of a 
   !!                3-components field, i.e.
   !!
   !!                    \int\int u_{x}^{2} + u_{y}^{2} + u_{z}^{2} dxdy
   !!                
   !!                
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     cx: x-component 
   !> @param[in]     cy: y-component
   !> @param[in]     cz: z-component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this cml_int2d_ene3c
   ! [cml_int2d_ene3c]
   SUBROUTINE cml_int2d_ene3c( cx, cy, cz, norm_const, outval)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::  cx, cy, cz    ! components 
      REAL(wp)                                    ::  norm_const    ! components 
      REAL(wp)                    , INTENT(inout) ::  outval        ! components 

      !
      INTEGER                                     ::  ji, jj        ! dummy loop indices
      REAL(wp)                                    ::  sub_int       ! 2D workspace
      REAL(wp)                                    ::  dom_int       ! 2D workspace
      REAL(wp), DIMENSION(mppsize)                ::  zw1
      !
      INTEGER                                     ::  ierr

      !!----------------------------------------------------------------------
      !
      sub_int = 0._wp
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            sub_int = sub_int + ( e1u(ji,jj) * e2u(ji,jj) * cx(ji,jj) * cx(ji,jj) &
                            &   + e1v(ji,jj) * e2v(ji,jj) * cy(ji,jj) * cy(ji,jj) &
                            &   + e1t(ji,jj) * e2t(ji,jj) * cz(ji,jj) * cz(ji,jj) ) * norm_const

         END DO
      END DO

      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLGATHER( sub_int, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)
      dom_int = SUM(zw1)

      outval = outval + dom_int


   END SUBROUTINE cml_int2d_ene3c
   ! [cml_int2d_ene3c]



   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE int2d_3c ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes a 2-dimensional integral of a 3-components field, i.e.
   !!
   !!                    \int\int u_{x} + u_{y} + u_{z} dxdy
   !!                
   !!                
   !!
   !> @details
   !!                
   !! 
   !> @param[in]     cx: x-component 
   !> @param[in]     cy: y-component
   !> @param[in]     cz: z-component
   !! 
   !! @result  Update (ua,va) with the now noise based coriolis term trend
   !!
   !! @warning       This code has some strange index in variance tensor. 
   !!  
   !! @note          
   !! @todo              
   !!
   !!---------------------------------------------------------------------------
   ! @snippet this int2d_3c
   ! [int2d_3c]
   SUBROUTINE int2d_3c( cx, cy, cz, fieldname, norm_const )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::  cx, cy, cz    ! components 
      CHARACTER(len = * )         , INTENT(in   ) ::  fieldname     ! string to print
      REAL(wp)                                    ::  norm_const    ! components 
      !
      INTEGER                                     ::  ji, jj        ! dummy loop indices
      REAL(wp)                                    ::  sub_int       ! 2D workspace
      REAL(wp)                                    ::  dom_int       ! 2D workspace
      REAL(wp), DIMENSION(mppsize)                ::  zw1
      !
      INTEGER                                     ::  ierr, myID

      !!----------------------------------------------------------------------
      !
      sub_int = 0._wp
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            sub_int = sub_int + ( e1u(ji,jj) * e2u(ji,jj) * cx(ji,jj) &
                            &   + e1v(ji,jj) * e2v(ji,jj) * cy(ji,jj) &
                            &   + e1t(ji,jj) * e2t(ji,jj) * cz(ji,jj) ) * norm_const
         END DO
      END DO

      call MPI_BARRIER(mpi_comm_oce, ierr)
      call MPI_ALLGATHER( sub_int, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)
      dom_int = SUM(zw1)

      call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
      IF ( myID == 0 ) THEN
         print '(A36, E16.7)', fieldname//'      2D domain integral value:', dom_int 
      END IF

   END SUBROUTINE int2d_3c
   ! [int2d_3c]


   SUBROUTINE eigensum(nn_tco_dat_in, nm_mod_nc_in, tlu_lam_spc_in  )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE   ***
      !!
      !! ** Purpose :  
      !!
      !! ** Method  :   
      !!                
      !!                
      !!------------------------------------------------------------------------
      INTEGER                             :: ji, jj
      INTEGER :: nn_tco_dat_IN
      CHARACTER(len = * )         , INTENT(in   ) ::  nm_mod_nc_in     ! string to print
      REAL(wp)                    , INTENT(in   ) ::  tlu_lam_spc_in     ! string to print

      integer, allocatable                ::   supp(:)
      integer, allocatable                ::   nmz2(:)
      integer, allocatable                ::   midz2(:)

      real(wp), allocatable :: eigens1(:,:)
      real(wp), allocatable :: eigens2(:,:)

      INTEGER(i4) ::  ierr(2) 
      !
      INTEGER                                     ::   myID
 
      ALLOCATE( supp(nn_tco_dat_in),       &
             &  nmz2(jpk-1), midz2(jpk) , &
             &  eigens1(nn_tco_dat_in, jpk-1),       &
             &  eigens2(nn_tco_dat_in, jpk-1), stat=ierr(1) )
      OPEN(131, FILE=TRIM(nm_mod_nc_in),FORM='FORMATTED', STATUS = 'OLD')
      !
      ! Read Eigenvalues and cumulate:
      !      eigens(k) = \sum_{j=1}^{k} \sqrt(\lambda_{j})
      !
      READ(131,*) (eigens1(1,jj) , jj = 1, jpk-1) 
      DO ji=2, nn_tco_dat_in
         READ(131,*)(eigens1(ji,jj) , jj = 1, jpk-1) 
         eigens1(ji,:) = eigens1(ji,:) + eigens1(ji-1,:)
      ENDDO
      eigens2 = eigens1
      !
      ! Column by column (i.e. layer by layer) normalize by \sum_{k}^{N}\lambda_{k}
      !
      DO ji = 1, jpk-1
         supp(:) = 0
         eigens1(:,ji) = eigens1(:,ji)/eigens1(nn_tco_dat_in,ji)
         !
         ! Set the treshold "tlu_lam_spec" 
         !
         WHERE (eigens1(:,ji) <= tlu_lam_spc_in)  supp = 1
         nmz2(ji) = sum( supp ) 
      END DO
      !
      ! Shifts the array down one position, to ensure the right amount of energy in
      ! those components requiring a vertical interpolation between modes
      !
      nmz2 = [nmz2(1), nmz2(1:jpk-1)]
      !
      ! Compute index vector
      !
      midz2(1) = 1
      DO ji = 2, jpk
         midz2(ji) = midz2(ji - 1) + nmz2(ji)
      END DO

      call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
      IF ( myID == 0 ) THEN
         print '(A20, I4, A12, E16.7)', ' Sum of first ',nmz2(1), ' eigenvalues', eigens2(nmz2(1),1)
      END IF
      CLOSE(131)
      DEALLOCATE( eigens1, supp, stat=ierr(2) )


   END SUBROUTINE eigensum


   SUBROUTINE sumvars( kt )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE   ***
      !!
      !! ** Purpose :  
      !!
      !! ** Method  :   
      !!                
      !!                
      !!------------------------------------------------------------------------
      INTEGER, INTENT(in   )                          :: kt
      INTEGER                                         ::   ji, jj, jk  ! dummy loop indices


      IF ( kt == it0_counter) stat_counter = 0
  
      IF (.not. ALLOCATED( usum )) THEN
         ALLOCATE( usum (jpi,jpj,jpk),      &
             &     usum2(jpi,jpj,jpk), stat=ierr )
         usum = 0._wp
         usum2 = 0._wp
      END IF 
 
      IF (.not. ALLOCATED( vsum )) THEN
         ALLOCATE( vsum (jpi,jpj,jpk),      &
             &     vsum2(jpi,jpj,jpk), stat=ierr )
         vsum = 0._wp
         vsum2 = 0._wp
      END IF 
 
      IF (.not. ALLOCATED( tsum )) THEN
         ALLOCATE( tsum (jpi,jpj,jpk),      &
             &     tsum2(jpi,jpj,jpk), stat=ierr )
         tsum = 0._wp
         tsum2 = 0._wp
      END IF 
 
      IF (.not. ALLOCATED( ssum )) THEN
         ALLOCATE( ssum (jpi,jpj,jpk),      &
             &     ssum2(jpi,jpj,jpk), stat=ierr )
         ssum = 0._wp
         ssum2 = 0._wp
      END IF 
 
      IF (.not. ALLOCATED( vortsum )) THEN
         ALLOCATE( vortsum (jpi,jpj,jpk),      &
             &     vortsum2(jpi,jpj,jpk), stat=ierr )
         vortsum = 0._wp
         vortsum2 = 0._wp
      END IF

      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
                  usum (ji,jj) = usum (ji,jj,jk) + un(ji,jj,jk)
                  usum2(ji,jj) = usum2(ji,jj,jk) + un(ji,jj,jk)**2
                  vsum (ji,jj) = vsum (ji,jj,jk) + vn(ji,jj,jk)
                  vsum2(ji,jj) = vsum2(ji,jj,jk) + vn(ji,jj,jk)**2
                  wsum (ji,jj) = wsum (ji,jj,jk) + wn(ji,jj,jk)
                  wsum2(ji,jj) = wsum2(ji,jj,jk) + wn(ji,jj,jk)**2
                  tsum (ji,jj) = tsum (ji,jj,jk) + tsn(ji,jj,jk,1)
                  tsum2(ji,jj) = tsum2(ji,jj,jk) + tsn(ji,jj,jk,1)**2
                  ssum (ji,jj) = ssum (ji,jj,jk) + tsn(ji,jj,jk,2)
                  ssum2(ji,jj) = ssum2(ji,jj,jk) + tsn(ji,jj,jk,2)**2
                  vortsum (ji,jj) = vortsum (ji,jj,jk) + tvor(ji,jj,jk)
                  vortsum2(ji,jj) = vortsum2(ji,jj,jk) + tvor(ji,jj,jk)**2
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============


      

   END SUBROUTINE sumvars





END MODULE tludgns




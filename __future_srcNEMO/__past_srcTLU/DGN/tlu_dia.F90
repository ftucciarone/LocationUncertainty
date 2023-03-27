!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tludia
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
MODULE tludia
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
   USE tlunoi           ! needed to compute the noise
   ! [mod_dep]
   !
   IMPLICIT NONE
   PRIVATE              ! Make stuff private by default
   !
   ! [public_sub]
   PUBLIC tlu_dia
   PUBLIC noise_energy
   
   
   
   
   
   ! [public_sub]
   !
   INCLUDE 'mpif.h'
   !
   ! [tlu_EnergyBalance]
   INTEGER,  PARAMETER                          ::   n_samples = 10     !< @private number of samples in the ensemble
   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      ::   int_TEnoi_hh          !< @public Energy of the noise 
   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      ::   int_TEdif_hh          !< @public Energy of the diffusion 
   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      ::   int_TEnoi_zz          !< @public Energy of the noise 
   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      ::   int_TEdif_zz          !< @public Energy of the diffusion 


   INTEGER,  PARAMETER                          ::   sim_time = 200    !< @private number of samples in the ensemble

   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      :: Rt_NoD_hh = 1._wp
   REAL(wp), PUBLIC, SAVE, DIMENSION(jpts)      :: Rt_NoD_zz = 1._wp



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
   !> @snippet this tlu_dia 
   ! [tlu_dia]
   SUBROUTINE tlu_dia
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


      !
      !
      ! Initialise
      !------------------------------




      !
      !
      IF (ierr/=0) THEN   ! check allocation success (0=OK)
         WRITE(numout,*) ' tlu_dgns(): allocation failed = ', ierr   ! [TODO] should be ctmp1? instead of numout
         !CALL ctl_stop( ctmp1 )   !< @todo enable
         STOP
      END IF
      !
      !
   END SUBROUTINE tlu_dia
   ! [tlu_dia]

   SUBROUTINE div_cpt( kt, pun, pvn, pwn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cpt  ***
      !!                    
      !! ** Purpose :   compute the horizontal divergence at now time-step
      !!
      !! ** Method  :   the now divergence is computed as :
      !!         hdivn = 1/(e1e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      and correct with runoff inflow (div_rnf) and cross land flow (div_cla) 
      !!
      !! ** Action  : - update hdivn, the now horizontal divergence
      !!----------------------------------------------------------------------
      INTEGER,                          INTENT(in   ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pun, pvn, pwn
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)         ::   tdivn
      !         
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zraur, zdep   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('div_cpt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cpt : horizontal velocity divergence '
         IF(lwp) WRITE(numout,*) '~~~~~~~   '
      ENDIF
      !
      ALLOCATE( tdivn(jpi,jpj,jpk) )

      DO jk = 1, jpkm1                                      !==  Horizontal divergence  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               tdivn(ji,jj,jk) = (  e2u(ji  ,jj  ) * e3u_n(ji  ,jj  ,jk) * pun(ji  ,jj  ,jk  )      &
                  &               - e2u(ji-1,jj  ) * e3u_n(ji-1,jj  ,jk) * pun(ji-1,jj  ,jk  )      &
                  &               + e1v(ji  ,jj  ) * e3v_n(ji  ,jj  ,jk) * pvn(ji  ,jj  ,jk  )      &
                  &               - e1v(ji  ,jj-1) * e3v_n(ji  ,jj-1,jk) * pvn(ji  ,jj-1,jk  )      &
                  &               + e1t(ji  ,jj  ) * e2t  (ji  ,jj     ) * pwn(ji  ,jj  ,jk  )      &
                  &               - e1t(ji  ,jj  ) * e2t  (ji  ,jj-1   ) * pwn(ji  ,jj-1,jk+1)  )   &
                  &            * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
            END DO  
         END DO  
      END DO
      !
      CALL lbc_lnk( 'div_cpt', tdivn, 'T', 1. )   !   (no sign change)
      !
   END SUBROUTINE div_cpt

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE noise_energy ***
   !!
   !> @authors
   !>                F. Tucciarone
   !!
   !> @version
   !!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         
   !!
   !!               
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
   ! @snippet this noise_energy
   ! [noise_energy]
   SUBROUTINE noise_energy( kt, r_Scale )
   USE tluhdifTRA
   USE tlunoi

      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpts)            , INTENT(in   ) ::   r_Scale               ! Schmidt number
      !
      INTEGER                                     ::  js            ! sample index
      INTEGER                                     ::  jn            ! tracer index
      INTEGER                                     ::  ji, jj, jk    ! dummy loop indices
      REAL(wp), ALLOCATABLE, DIMENSION(:)         ::  sub_int_hh    ! 2D workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:)         ::  sub_int_zz    ! 2D workspace
      REAL(wp)                                    ::  sum_noi       ! 2D workspace
      REAL(wp), DIMENSION(mppsize)                ::  gather_hh, gather_zz
      !
      INTEGER                                     ::  ierr, myID
      REAL(wp), DIMENSION(jpi,jpj,jpk,2)          ::  zwa0                 ! Average in i of pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk,1)          ::  fakeazz                ! Average in i of pcn



      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwa1, zwa2, zwa3
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwpu, zwpv, zwpw   ! workspace
      !!----------------------------------------------------------------------
      int_TEnoi_hh = 0._wp
      int_TEnoi_zz = 0._wp
      int_TEdif_hh = 0._wp
      int_TEdif_zz = 0._wp
      Rt_NoD_hh    = 1._wp
      Rt_NoD_zz    = 1._wp


      ALLOCATE( sub_int_hh(jpts), sub_int_zz(jpts) )
      ALLOCATE( zwpu(jpi,jpj,jpk),   & 
              & zwpv(jpi,jpj,jpk),   &
              & zwpw(jpi,jpj,jpk),   &
              & zwa1(jpi,jpj,jpk),   &
              & zwa2(jpi,jpj,jpk),   &
              & zwa3(jpi,jpj,jpk)   )

      DO js = 1, n_samples
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         CALL tlu_noi( kt ) 
         ! 
         sub_int_hh = 0._wp
         sub_int_zz = 0._wp
         !
         gather_hh = 0._wp
         gather_zz = 0._wp
         !
         zwa1 = 0._wp
         zwa2 = 0._wp
         zwa3 = 0._wp
         !
         zwa1(:,:,1:jpkm1) = unoi(:,:,1:jpkm1) / spread(e1u,3,jpkm1)        !adding scale factor
         zwa2(:,:,1:jpkm1) = vnoi(:,:,1:jpkm1) / spread(e2v,3,jpkm1)        !adding scale factor
         zwa3(:,:,2:jpkm1) = wnoi(:,:,2:jpkm1) / e3w_b(:,:,2:jpkm1)         !adding scale factor
         !
         !                             ! =========== !
         DO jn = 1, jpts               ! tracer loop !
            !                          ! =========== ! 
            zwpu = 0._wp
            zwpv = 0._wp 
            zwpw = 0._wp
            !
            ! x-Directed gradient with interpolation in x
            !
            zwpu(2:jpim1,2:jpjm1,1:jpkm1) = ( zwa1(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( tsn(3:jpi  ,2:jpjm1,1:jpkm1,jn) - tsn(2:jpim1  ,2:jpjm1,1:jpkm1,jn) )   &
            &                               + zwa1(1:jpim1-1,2:jpjm1  ,1:jpkm1) * ( tsn(2:jpim1,2:jpjm1,1:jpkm1,jn) - tsn(1:jpim1-1,2:jpjm1,1:jpkm1,jn) ) ) &
            &                              * tmask(2:jpim1,2:jpjm1,1:jpkm1)
            !
            ! y-Directed gradient with interpolation in y
            !
            zwpv(2:jpim1,2:jpjm1,1:jpkm1) = ( zwa2(2:jpim1  ,2:jpjm1  ,1:jpkm1) * ( ptn(2:jpim1,3:jpj  ,1:jpkm1,jn) - ptn(2:jpim1,2:jpjm1  ,1:jpkm1,jn) ) &
            &                               + zwa2(2:jpim1  ,1:jpjm1-1,1:jpkm1) * ( ptn(2:jpim1,2:jpjm1,1:jpkm1,jn) - ptn(2:jpim1,1:jpjm1-1,1:jpkm1,jn) ) ) &
            &                              * tmask(2:jpim1,2:jpjm1,1:jpkm1)
            !
            ! z-Directed gradient
            !
            zwpw(2:jpim1,2:jpjm1,2:jpkm1) = ( zwa3(2:jpim1,2:jpjm1,2:jpkm1) * ( tsn(2:jpim1,2:jpjm1,1:jpkm1-1,jn) - tsn(2:jpim1,2:jpjm1,2:jpkm1,jn) ) &
            &                               + zwa3(2:jpim1,2:jpjm1,3:jpk  ) * ( tsn(2:jpim1,2:jpjm1,2:jpkm1  ,jn) - tsn(2:jpim1,2:jpjm1,3:jpk  ,jn) ) ) &
            &                              * tmask(2:jpim1,2:jpjm1,1:jpkm1) 
            !
            ! Ensure bottom and surface boundary conditions
            zwpw(:,:,jpk) = 0._wp
            zwpw(:,:, 1 ) = 0._wp
            sub_int_hh(jn) = SUM( ( ( 0.5_wp * ( zwpu + zwpv + zwpw ) / r_Scale(jn) )**2 - zwpw(ji,jj,jk)**2 ) * spread(e1e2t,3,jpk) * e3t_n ) 
            !
            ! Energy of vertical component
            !
            sub_int_zz(jn) =  SUM( ( ( zwpw(2:jpim1,2:jpjm1,2:jpkm1) / r_Scale(jn) )**2 ) &
                           & * spread(e1e2t(2:jpim1,2:jpjm1),3,jpk) * e3t_n(2:jpim1,2:jpjm1,2:jpkm1)  )
            !
            call MPI_BARRIER(mpi_comm_oce, ierr)
            call MPI_ALLGATHER( sub_int_zz(jn), 1, mpi_double_precision, gather_hh, 1, mpi_double_precision, mpi_comm_oce, ierr)
            call MPI_ALLGATHER( sub_int_hh(jn), 1, mpi_double_precision, gather_zz, 1, mpi_double_precision, mpi_comm_oce, ierr)
            call MPI_BARRIER(mpi_comm_oce, ierr)
            int_TEnoi_hh(jn) = int_TEnoi_hh (jn) + SUM(gather_hh)
            int_TEnoi_zz(jn) = int_TEnoi_zz (jn) + SUM(gather_zz)
         END DO
         !
      END DO
      !
      DEALLOCATE(zwa1, zwa2, zwa3, zwpw)
      !
      int_TEnoi_hh = int_TEnoi_hh / (n_samples - 1)
      int_TEnoi_zz = int_TEnoi_zz / (n_samples - 1)
      !
      !
      zwa0 = 0._wp
      sub_int_hh = 0._wp
      sub_int_zz = 0._wp
      CALL tlu_trahhdiff( kt, tsn, zwa0, r_Scale)
      CALL tlu_trahzdiff( kt, tsn, zwa0, r_Scale)
      !
      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !  
         sub_int_hh(jn) = SUM( 2._wp * zwa0(2:jpim1,2:jpjm1,2:jpkm1,jn)   &
                        &             * tsn(2:jpim1,2:jpjm1,2:jpkm1,jn)   &
                        &    * spread(e1e2t(2:jpim1,2:jpjm1,2:jpkm1),3,jpk) * e3t_n(2:jpim1,2:jpjm1,2:jpkm1) ) / rn_rdt
         !
         ! Energy of vertical component
         !
         sub_int_zz(jn) = SUM( ( ( ( var_ten(:,:,1:jpkm1-1, 3) + var_ten(:,:,2:jpkm1  , 3) )                        &
                        &        * (     tsn(:,:,1:jpkm1-1,jn) -     tsn(:,:,2:jpkm1  ,jn) ) / e3w_n(:,:,2:jpkm1) ) &
                        &      - ( ( var_ten(:,:,2:jpkm1  , 3) + var_ten(:,:,3:jpk    , 3) )                        &
                        &        * (     tsn(:,:,2:jpkm1  ,jn) -     tsn(:,:,3:jpk    ,jn) ) / e3w_n(:,:,3:jpk  ) ) &
                        &        / ( 2._wp * e3t_n(:,:,2:jpkm1) * r_Scale(jn)**2 ) )                                &
                        &        * spread(e1e2t(2:jpim1,2:jpjm1),3,jpk) * e3t_n(2:jpim1,2:jpjm1,2:jpkm1)            &
                        &        * ( tsn(2:jpim1,2:jpjm1,2:jpkm1) ) ) / rn_rdt 
         !
         call MPI_BARRIER(mpi_comm_oce, ierr)
         call MPI_ALLGATHER( sub_int_hh(jn), 1, mpi_double_precision, gather_hh, 1, mpi_double_precision, mpi_comm_oce, ierr)
         call MPI_ALLGATHER( sub_int_zz(jn), 1, mpi_double_precision, gather_zz, 1, mpi_double_precision, mpi_comm_oce, ierr)
         call MPI_BARRIER(mpi_comm_oce, ierr)
         !
         int_TEdif_hh(jn) = -1._wp * SUM(gather)
         int_TEdif_zz(jn) = -1._wp * SUM(gather)
         !
      END DO

      IF (int_TEnoi_hh(1) == 0) int_TEnoi_hh(1) = 1._wp
      IF (int_TEnoi_hh(2) == 0) int_TEnoi_hh(2) = 1._wp
      IF (int_TEnoi_zz(1) == 0) int_TEnoi_zz(1) = 1._wp
      IF (int_TEnoi_zz(2) == 0) int_TEnoi_zz(2) = 1._wp

      Rt_NoD_hh(1) = int_TEdif_hh(1)/int_TEnoi_hh(1)
      Rt_NoD_hh(2) = int_TEdif_hh(2)/int_TEnoi_hh(2)
      Rt_NoD_zz(1) = int_TEdif_zz(1)/int_TEnoi_zz(1)
      Rt_NoD_zz(2) = int_TEdif_zz(2)/int_TEnoi_zz(2) 

      IF (Rt_NoD_hh(1) == 0) Rt_NoD_hh(1) = 1._wp
      IF (Rt_NoD_hh(2) == 0) Rt_NoD_hh(2) = 1._wp
      IF (Rt_NoD_zz(1) == 0) Rt_NoD_zz(1) = 1._wp
      IF (Rt_NoD_zz(2) == 0) Rt_NoD_zz(2) = 1._wp


      call MPI_Comm_rank(MPI_COMM_WORLD, myID, ierr)
      IF ( myID == 0 ) THEN

         Rt_NoD_hh(kt - nn_it000 + 1, 1) = int_TEnoi_hh(1)/int_TEdif_hh(1)
         Rt_NoD_hh(kt - nn_it000 + 1, 2) = int_TEnoi_hh(2)/int_TEdif_hh(2)
         Rt_NoD_zz(kt - nn_it000 + 1, 1) = int_TEnoi_zz(1)/int_TEdif_zz(1)
         Rt_NoD_zz(kt - nn_it000 + 1, 2) = int_TEnoi_zz(2)/int_TEdif_zz(2)

         print '(A37, E16.7, E16.7, E16.7, E16.7)', 'Empirical expectation noise (3D int):', &
                & int_TEnoi_hh(1), int_TEnoi_hh(2), int_TEnoi_zz(1), int_TEnoi_zz(2)
         print '(A37, E16.7, E16.7, E16.7, E16.7)', ' Variance expectation noise (3D int):', &
                & int_TEdif_hh(1), int_TEdif_hh(2), int_TEdif_zz(1), int_TEdif_zz(2)
         print '(A37, E16.7, E16.7, E16.7, E16.7)', '  ----> Noise vs Variance ratio <----', &
                & int_TEnoi_hh(1)/int_TEdif_hh(1), int_TEnoi_hh(2)/int_TEdif_hh(2), &
                & int_TEnoi_zz(1)/int_TEdif_zz(1), int_TEnoi_zz(2)/int_TEdif_zz(2)
      END IF

   END SUBROUTINE noise_energy
   ! [noise_energy]





END MODULE tludia




!! ------------------------------------------------------------------------------
!! 
!!                               MODULE: tlu_noi
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
!> @brief         Contains the initialization of the noise and the procedure-independent
!!                routines for noise generation
!! 
!! 
!! @par           Procedure specifics      
!> @details       Each routine in this modue is either independent of the noise 
!!                chosen (viz. the Ito-Stokes drift) or discriminates between the 
!!                noise models and refers to the relative subroutines.
!!
!> @snippet       this tlu_noise
!!                \f{align*}{
!!                  \texttt{unoi} &\leftarrow d\boldsymbol{\eta}_{t,x} \\
!!                  \texttt{vnoi} &\leftarrow d\boldsymbol{\eta}_{t,y} \\
!!                  \texttt{wnoi} &\leftarrow d\boldsymbol{\eta}_{t,z} 
!!               \f}
!!
!> @snippet      this tlu_bias
!!               \f{align*}{
!!                  \texttt{ubias} &\leftarrow d\boldsymbol{\eta}_{t,x}^{b} \\
!!                  \texttt{vbias} &\leftarrow d\boldsymbol{\eta}_{t,y}^{b} \\
!!                  \texttt{wbias} &\leftarrow d\boldsymbol{\eta}_{t,z}^{b} 
!!               \f}
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
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
!> @todo          read namelist from _cfg, not only _ref
!! @todo          write error messages to appropriate unit
!!
!!------------------------------------------------------------------------------
MODULE tlunoi
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
   PUBLIC tlu_noi_init    ! called by tlu_init in tlu.f90
   PUBLIC tlu_fields      ! called by step.f90
   PUBLIC tlu_noi         ! called by tlu_fields and noise_energy in tlu_dgns.F90
   PUBLIC r8vec_normal_01
   !
   INCLUDE 'mpif.h'
   !
   REAL(wp),    PUBLIC    :: cof_tke     
   REAL(wp)               :: KE_ratio_b, KE_ratio_n 
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_noi_init ***
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
   !!                
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
   ! @snippet this tlu_noi_init
   ! [tlu_noi_init]
   SUBROUTINE tlu_noi_init
      USE tlupod
      USE tludmd
      USE tlusvd
      !
      INTEGER(i4) :: ios, ierr(4), chk   ! namelist output, allocation statistics
      INTEGER     :: ji
      LOGICAL :: file_exists
      !
      ! Control print (1)
      !
      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_noi_init : TLU, Proper orthogonal decomposition  noise '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtlu_noi'
         WRITE(numout,*) '      data-based noise model           ln_tlu_pod  = ', ln_tlu_pod
         WRITE(numout,*) '      data-based noise model           ln_tlu_dmd  = ', ln_tlu_dmd
         WRITE(numout,*) '       data-free noise model           ln_tlu_svd  = ', ln_tlu_svd
      END IF

      KE_ratio_b = 1._wp
      KE_ratio_n = 1._wp
      !
      ! SVD-Based noise
      !
      IF ( ln_tlu_pod ) THEN
         !
         ! Allocate index vectors
         !
          CALL tlu_init_pod
         !
      END IF
      !
      ! DMD-Based noise
      !
      IF ( ln_tlu_dmd ) THEN
         !
         ! Allocate index vectors
         !
          CALL tlu_init_dmd
         !
      END IF
      !
      ! Pseudo-Observations model
      !
      IF ( ln_tlu_svd ) THEN
         !
         ! Allocate index vectors
         !
          CALL tlu_init_svd
         !
      END IF      
      !
   END SUBROUTINE tlu_noi_init
   ! [tlu_noi_init]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_fields ***
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
   !!                
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
   ! @snippet this tlu_fields
   ! [tlu_fields]
   SUBROUTINE tlu_fields( kt )
      USE tlupod
      USE tludmd
      USE tlusvd
      USE tluprj
      !
      INTEGER, INTENT(in   ) :: kt         ! ocean time-step index
      REAL(wp)                                  :: tic, toc      
      !
      LOGICAL                :: initialize = .TRUE.
      INTEGER                :: ierr
      !!------------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_noi') ![NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_fields : TLU fields construction'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF
 !     IF ( .NOT. initialize ) THEN
         unoi = 0._wp
         vnoi = 0._wp
         wnoi = 0._wp
 !     END iF
      !
      ! POD-Based noise
      !
      IF ( ln_tlu_pod ) THEN
         ! 
         IF( (kt == nit000) ) THEN      ! Stationary assumption here
            !
            ! Compute the variance tensor a = sum_{k} \phi_{k}^{0}\phi_{k}^{0} * dt and Ito-Stokes drift
            !
            !    tlu_var(   x_mod    y_mod,   z_mod,     nn_nmod,     dt, initialize,  lbc_comm )
            CALL tlu_var( spx_mod, spy_mod, spz_mod, nn_tlu_nmod, rn_rdt,     .TRUE.,    .TRUE. )
            CALL tlu_isd( uisd_n, visd_n, wisd_n, initialize )
            !           
            IF ( ln_pyp ) CALL tlu_tnsf_isd
            IF ( ln_tlu_bia ) CALL tlu_tnsf_bia	    
            !
         END IF
         !
         ! Draw k random (standard) gaussian distribution
         !
         CALL r8vec_normal_01( kt,   nn_tlu_nmod,  brwn_rv, .FALSE. )  
         CALL MPI_BCAST(  brwn_rv,   nn_tlu_nmod, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr )
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         !    tlu_noi(   x_mod,   y_mod,   z_mod,     nn_nmod, time_coeff, initialize, lbc_comm )
         CALL tlu_noi( spx_mod, spy_mod, spz_mod, nn_tlu_nmod,   rbrwn_rv,     .TRUE.,   .TRUE. )
         !         
      END IF
      !
      ! DMD-Based noise
      !
      IF ( ln_tlu_dmd ) THEN
         ! 
         IF( (kt == nit000) ) THEN      ! Stationary assumption here
            !
            ! Compute the variance tensor a = sum_{k} \phi_{k}^{0}\phi_{k}^{0} * dt and Ito-Stokes drift
            !
            !    tlu_var(     x_mod,     y_mod,     z_mod,       nn_nmod,     dt, initialize,  lbc_comm )      
            CALL tlu_var( rpx_mod_r, rpy_mod_r, rpz_mod_r, nn_tlu_nmod_r, rn_rdt,     .TRUE.,  .FALSE. )
            CALL tlu_var( ipx_mod_r, ipy_mod_r, ipz_mod_r, nn_tlu_nmod_r, rn_rdt,    .FALSE.,   .TRUE. )
            CALL tlu_isd( uisd_n, visd_n, wisd_n, initialize )
            !          
            IF ( ln_pyp ) CALL tlu_tnsf_isd
            !
         END IF
         !
         ! Draw k random (standard) gaussian distribution
         !
         CALL r8vec_normal_01( kt, nn_tlu_nmod_r, rbrwn_rv, .FALSE. ) 
         CALL r8vec_normal_01( kt, nn_tlu_nmod_r, ibrwn_rv, .FALSE. ) 
         CALL MPI_BCAST( rbrwn_rv, nn_tlu_nmod_r, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
         CALL MPI_BCAST( ibrwn_rv, nn_tlu_nmod_r, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
         !
         ! Compute the time coefficients
         !
         !    tlu_tcoef_dmd(            time,    freq,    rgaus,    igaus, rBtcos_iBtsin, rBtsin_iBtcos )
         CALL tlu_tcoef_dmd(  REAL( kt, wp ), omega_r, rbrwn_rv, ibrwn_rv,    rtime_coef,    itime_coef )
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         !    tlu_noi(     x_mod,     y_mod,     z_mod,       nn_nmod, time_coeff, initialize, lbc_comm )
         CALL tlu_noi( rpx_mod_r, rpy_mod_r, rpz_mod_r, nn_tlu_nmod_r, rtime_coef,     .TRUE.,  .FALSE. )
         CALL tlu_noi( ipx_mod_r, ipy_mod_r, ipz_mod_r, nn_tlu_nmod_r, itime_coef,    .FALSE.,   .TRUE. )
         !
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         IF (ln_tlu_bia) CALL tlu_bia_dmd ( REAL( kt, wp ), omega_c )
         !   
      END IF
      !
      ! SVD-Based noise
      !
      IF ( ln_tlu_svd ) THEN
         !
         ! Generate the random modes
         !
         CALL gen_vel_mod( kt )
         !
         ! Draw k random (standard) gaussian distribution
         !
         CALL r8vec_normal_01( kt,   nn_tlu_pobs,  gaus_rv, .FALSE. )  
         CALL MPI_BCAST(  gaus_rv,   nn_tlu_pobs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
         ! 
         ! Create one realisation of sigma dBt = sum_{k} \phi_{k} \xi_{k}
         !
         !    tlu_noi(   x_mod,   y_mod,   z_mod,     nn_nmod, time_coeff, initialize, lbc_comm )
         CALL tlu_noi( psx_mod, psy_mod, psz_mod, nn_tlu_pobs,    gaus_rv,     .TRUE.,   .TRUE. )
         !        

         !    tlu_var(   x_mod    y_mod,   z_mod,     nn_nmod,     dt, initialize,  lbc_comm )
         CALL tlu_var( psx_mod, psy_mod, psz_mod, nn_tlu_pobs, rn_rdt,     .TRUE.,    .TRUE. )  
         CALL tlu_isd( uisd_n, visd_n, wisd_n, initialize )
         !          
         IF ( ln_pyp ) CALL tlu_tnsf_isd
         !
      END IF      
      !
      ! Projection on isopycnals
      !
      IF ( ln_pyp ) THEN
         !CALL tlu_var_proj( var_ten, var_ten_n )
         CALL  tlu_vel_proj( uisd_0, visd_0, wisd_0, uisd_n, visd_n, wisd_n )
         CALL  tlu_vel_proj(   unoi,   vnoi,   wnoi,   unoi,   vnoi,   wnoi )
         !    
         IF (ln_tlu_bia) CALL tlu_vel_proj( ubia_n, vbia_n, wbia_n, ubia_n, vbia_n, wbia_n )
         IF (ln_tlu_bia) wbia_n(:,:,1) = 0._wp
         !
      END IF
      ! 
      ! Rescale the noise by some manipulation
      ! dBt = C * (      sum_{k} \phi_{k} \xi_{k} )
      !   a = C * ( dt * sum_{k} \phi_{k}\phi_{k} )
      !
      call tlu_norm_mod ( kt, .TRUE. )
      !
      IF(  (kt == nit000) .AND. (.NOT. ln_pyp) ) THEN 
         !
         ! Print variance tensor (a) at time t0 i.e. stationary
         !
         CALL iom_put( 'var_uu_once', var_ten(:,:,:,1) )  
         CALL iom_put( 'var_vv_once', var_ten(:,:,:,2) )  
         CALL iom_put( 'var_ww_once', var_ten(:,:,:,3) )  
         CALL iom_put( 'var_uv_once', var_ten(:,:,:,4) )  
         CALL iom_put( 'var_uw_once', var_ten(:,:,:,5) )  
         CALL iom_put( 'var_vw_once', var_ten(:,:,:,6) ) 
         !
         ! Print Ito-Stokes Drift
         !
         CALL iom_put( 'U_isd_once', uisd_n )  
         CALL iom_put( 'V_isd_once', visd_n )  
         CALL iom_put( 'W_isd_once', wisd_n )   
         !
      ELSE IF ( ln_pyp ) THEN
         !
         ! Print variance tensor (a) at time t0 i.e. stationary
         !
         CALL iom_put( 'var_uu', var_ten(:,:,:,1) )  
         CALL iom_put( 'var_vv', var_ten(:,:,:,2) )  
         CALL iom_put( 'var_ww', var_ten(:,:,:,3) )  
         CALL iom_put( 'var_uv', var_ten(:,:,:,4) )  
         CALL iom_put( 'var_uw', var_ten(:,:,:,5) )  
         CALL iom_put( 'var_vw', var_ten(:,:,:,6) ) 
         !
         ! Print Ito-Stokes Drift
         !
         CALL iom_put( 'U_isd', uisd_n )  
         CALL iom_put( 'V_isd', visd_n )  
         CALL iom_put( 'W_isd', wisd_n )   
         !
      END IF
      !
      ! Print noise fields
      !
      CALL iom_put( "unoi", unoi )
      CALL iom_put( "vnoi", vnoi )
      CALL iom_put( "wnoi", wnoi )
      !
      IF (ln_tlu_bia) THEN
         !
         ! Print bias fields
         !
         CALL iom_put( "ubia_n", ubia_n )
         CALL iom_put( "vbia_n", vbia_n )
         CALL iom_put( "wbia_n", wbia_n )
	 !
      END IF
      !
      IF( ln_timing )  CALL timing_stop('tlu_noi')   ! [NEMO] check
      !
   END SUBROUTINE tlu_fields
   ! [tlu_fields]   


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tnsf_isd ***
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
   !!                
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
   ! @snippet this tlu_tnsf_isd
   ! [tlu_tnsf_isd]
   SUBROUTINE tlu_tnsf_isd
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_tnsf_isd  ***
      !!
      !! ** Purpose :  Copy the initial modes into the now modes
      !!
      !! ** Action  :  Simple assignation and deallocation 
      !!
      !!------------------------------------------------------------------------
      !
      uisd_0 = uisd_n
      visd_0 = visd_n
      wisd_0 = wisd_n
      !
   END SUBROUTINE tlu_tnsf_isd
   ! [tlu_tnsf_isd]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_tnsf_bia ***
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
   !!                
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
   ! @snippet this tlu_tnsf_bia
   ! [tlu_tnsf_bia]
   SUBROUTINE tlu_tnsf_bia
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_tnsf_bia  ***
      !!
      !! ** Purpose :  Copy the initial modes into the now modes
      !!
      !! ** Action  :  Simple assignation and deallocation 
      !!
      !!------------------------------------------------------------------------
      !
      ubia_n = ubia_0
      vbia_n = vbia_0
      wbia_n = wbia_0
      !
      DEALLOCATE( ubia_0, &
      &           vbia_0, &
      &           wbia_0  )
      !
   END SUBROUTINE tlu_tnsf_bia
   ! [tlu_tnsf_bia]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_noi ***
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
   !!                
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
   ! @snippet this tlu_noi
   ! [tlu_noi]
   SUBROUTINE tlu_noi( x_mod, y_mod, z_mod, nn_nmod, time_coeff, initialize, lbc_comm )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE Create one realisation of noise  ***
      !!
      !! ** Purpose :   Formulate noise from spatial modes and time coeff/random numbers
      !!
      !! ** Method  :   Project spatial modes on corresponding temporal coefficient or
      !!                on gaussian random number as provided by input.
      !!                The noise is divided by the time step because later on in the 
      !!                NEMO deterministic kernel it will multiplied by rn_dt
      !!------------------------------------------------------------------------
      REAL(wp)                                  :: tic, toc
      
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: x_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: y_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: z_mod
      REAL(wp), DIMENSION(:),     INTENT(in   ) :: time_coeff
      INTEGER,                    INTENT(in   ) :: nn_nmod
      LOGICAL,                    INTENT(in   ) :: initialize
      LOGICAL,                    INTENT(in   ) :: lbc_comm
      INTEGER                                   :: jm, m_idx      
      !
      IF( ln_timing ) CALL timing_start('tlu_noi') ![NEMO] check
      !      
      IF ( initialize ) THEN
         unoi = 0._wp
         vnoi = 0._wp
         wnoi = 0._wp
      END IF
      !
      DO jm = 1, nn_nmod
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
         !
         unoi(:,:,:) = unoi(:,:,:) + x_mod(:,:,m_idx + 1 : m_idx + jpk ) * umask * time_coeff(jm)
         vnoi(:,:,:) = vnoi(:,:,:) + y_mod(:,:,m_idx + 1 : m_idx + jpk ) * vmask * time_coeff(jm)
         wnoi(:,:,:) = wnoi(:,:,:) + z_mod(:,:,m_idx + 1 : m_idx + jpk ) * wmask * time_coeff(jm)
         !
      ENDDO
      !
      ! Top and bottom boundary condition
      !
      wnoi(:,:,  1) = 0._wp
      wnoi(:,:,jpk) = 0._wp
      !
      ! Lateral boundary condition transfer across nodes
      !
      IF ( lbc_comm ) THEN
         !
         CALL lbc_lnk_multi( 'tlu_noi', unoi , 'U', -1., vnoi , 'V', -1., wnoi, 'W', -1.)
         !
      END IF
      !
   END SUBROUTINE tlu_noi
   ! [tlu_noi]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_var ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2023 -  ( F. TUCCIARONE   )  Ongoing work and documentation
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
   !!                
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
   ! @snippet this tlu_var
   ! [tlu_var]
   SUBROUTINE tlu_var( x_mod, y_mod, z_mod, nn_nmod, dt, initialize, lbc_comm )
      !!----------------------------------------------------------------------
      !!            ***  ROUTINE Create one realisation of noise  ***
      !!
      !! ** Purpose :   Formulate variance tensor from spatial modes and variance of time_coeff
      !!
      !! ** Method  :   Average spatial modes onto corresponding mesh weighted by the variance
      !!                calculated from each corresponding time-coefficient value (MM method)
      !!                or a constant (POD method)
      !!------------------------------------------------------------------------
      REAL(wp),                   INTENT(in  ) :: dt        ! Time interval
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: x_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: y_mod
      REAL(wp), DIMENSION(:,:,:), INTENT(in  ) :: z_mod

      INTEGER,                    INTENT(in  ) :: nn_nmod
      LOGICAL,                    INTENT(in  ) :: initialize
      LOGICAL,                    INTENT(in  ) :: lbc_comm

      INTEGER                                  :: jm, m_idx

      INTEGER, PARAMETER                       :: ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                       :: ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      ! Initialize variance tensor
      !
      IF ( initialize ) THEN
         var_ten(:,:,:,:) = 0._wp
      END IF
      !
      DO jm = 1, nn_nmod
         !
         ! Define zero-th indexed mode
         !
         m_idx = ( jm - 1 ) * jpk 
         !
         ! Create realization of sigma dBt = \sum_{k} \phi_{k}\xi_{k}
         !
         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia11) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia11)                              &
         !
         &                                               + ( x_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2   & 
         &                                                 + x_mod( 1:jpi-2, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia22) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia22)                              &
         !
         &                                               + ( y_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )**2   & 
         &                                                 + y_mod( 2:jpi-1, 1:jpj-2, m_idx + 1 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia33) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1,ia33)                              &
         !
         &                                               + ( z_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 1 )**2   & 
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk     )**2 ) &
         &                                               * 0.5_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia12) = var_ten(2:jpi-1, 2:jpj-1, 1:jpk-1, ia12)                             &
         !
         &                                               + ( x_mod( 2:jpi-1, 3:jpj-1, m_idx + 1 : m_idx + jpk     )      &
         &                                                 + x_mod( 2:jpi-1, 2:jpj-2, m_idx + 1 : m_idx + jpk     ) )    &
         &                                               * ( y_mod( 3:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk     )      &
         &                                                 + y_mod( 2:jpi-2, 2:jpj-1, m_idx + 1 : m_idx + jpk     ) )    &
         &                                               * 0.25_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia13) = var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia13)                             &
         !
         &                                               + ( x_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 2 )      &
         &                                                 + x_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * ( z_mod( 3:jpi  , 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 )      &
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * 0.25_wp * dt

         var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia23) = var_ten(2:jpi-1, 2:jpj-1, 2:jpk-1, ia23)                             &
         !
         &                                               + ( y_mod( 2:jpi-1, 2:jpj-1, m_idx + 1 : m_idx + jpk - 2 )      &
         &                                                 + y_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * ( z_mod( 2:jpi-1, 3:jpj  , m_idx + 2 : m_idx + jpk - 1 )      &
         &                                                 + z_mod( 2:jpi-1, 2:jpj-1, m_idx + 2 : m_idx + jpk - 1 ) )    &
         &                                               * 0.25_wp * dt
         !
      END DO
      !
      ! Mask the components
      !
      var_ten(:,:,:,ia11) = var_ten(:,:,:,ia11) * tmask
      var_ten(:,:,:,ia22) = var_ten(:,:,:,ia22) * tmask
      var_ten(:,:,:,ia33) = var_ten(:,:,:,ia33) * tmask
      var_ten(:,:,:,ia12) = var_ten(:,:,:,ia12) * fmask
      var_ten(:,:,:,ia13) = var_ten(:,:,:,ia13) * wumask
      var_ten(:,:,:,ia23) = var_ten(:,:,:,ia23) * wvmask
      !
      ! Lateral boundary condition transfer across nodes
      !
      IF ( lbc_comm ) THEN
         !
         CALL lbc_lnk_multi( 'tlu_var_pod', var_ten(:,:,:,ia11) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia22) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia33) , 'T', -1.,   &
         &                                  var_ten(:,:,:,ia12) , 'F', -1.,   &
         &                                  var_ten(:,:,:,ia13) , 'U', -1.,   &
         &                                  var_ten(:,:,:,ia23) , 'V', -1.    )
         !
      END IF
      !
   END SUBROUTINE tlu_var
   ! [tlu_var]
   

   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_isd ***
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
   ! @snippet this tlu_isd
   ! [tlu_isd]
   SUBROUTINE tlu_isd( pcu, pcv, pcw, initialize )
      LOGICAL,                          INTENT(in   )   :: initialize
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcu                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcv                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out)   ::   pcw                   ! before field (so its euler)
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)         ::   flux_var 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3u   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3v   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)           ::   r1_e1e2e3w        
      INTEGER                                           ::   ierr                
      !
      INTEGER                                           ::   jk                    ! dummy loop indices
      INTEGER, PARAMETER                                ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia21 = ndiffidx(2,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia31 = ndiffidx(3,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia32 = ndiffidx(3,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                ::   ia33 = ndiffidx(3,3)  ! Mapping the indeces of a
      !
      !!----------------------------------------------------------------------
      !
      !
      ! Initialize variance tensor
      !
      IF ( initialize ) THEN
         pcu = 0._wp
         pcv = 0._wp
         pcw = 0._wp
      END IF

      ALLOCATE( flux_var(jpi,jpj,jpk,9), stat=ierr )

      ALLOCATE( r1_e1e2e3u(jpi,jpj,jpk), &
      &         r1_e1e2e3v(jpi,jpj,jpk), &
      &         r1_e1e2e3w(jpi,jpj,jpk), stat=ierr )

      flux_var(:,:,:,:) = 0._wp
      flux_var(:,:,:,3) = var_ten(:,:,:,ia31) * 0.5_wp
      flux_var(:,:,:,6) = var_ten(:,:,:,ia32) * 0.5_wp
      flux_var(:,:,:,9) = var_ten(:,:,:,ia33) * 0.5_wp
      !                     
      ! Divergence is by column
      !
      DO jk = 1, jpkm1
         !
         r1_e1e2e3u(:,:,jk) = r1_e1e2u(:,:) * umask(:,:,jk) / e3u_n(:,:,jk) 
         r1_e1e2e3v(:,:,jk) = r1_e1e2v(:,:) * vmask(:,:,jk) / e3v_n(:,:,jk) 
         r1_e1e2e3w(:,:,jk) = r1_e1e2t(:,:) * wmask(:,:,jk) / e3w_n(:,:,jk)
         !
         flux_var(:,:,jk,1) = var_ten(:,:,jk,ia11) * e2t(:,:) * e3t_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,2) = var_ten(:,:,jk,ia21) * e1f(:,:) * e3f_n(:,:,jk) * 0.5_wp
         !
         flux_var(:,:,jk,4) = var_ten(:,:,jk,ia12) * e2f(:,:) * e3f_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,5) = var_ten(:,:,jk,ia22) * e1t(:,:) * e3t_n(:,:,jk) * 0.5_wp
         !
         flux_var(:,:,jk,7) = var_ten(:,:,jk,ia13) * e2u(:,:) * e3uw_n(:,:,jk) * 0.5_wp
         flux_var(:,:,jk,8) = var_ten(:,:,jk,ia23) * e1v(:,:) * e3vw_n(:,:,jk) * 0.5_wp
         !
         !
         !
         pcu(2:jpi -1, 2:jpj -1, jk   ) =        pcu( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3u( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_11)
         &                              * ( flux_var( 3:jpi   , 2:jpj -1, jk   , 1)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk   , 1)    &
         !  d_y(a_21)
         &                                + flux_var( 2:jpi -1, 2:jpj -1, jk   , 2)    &
         &                                - flux_var( 2:jpi -1, 1:jpj -2, jk   , 2) )  &
         !  d_z(a_31)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 3)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 3) )  &
         &                               /      e3u_n( 2:jpi -1, 2:jpj -1, jk   )
         !
         !
         !
         pcv(2:jpi -1, 2:jpj -1, jk   ) =        pcv( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3v( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_12)
         &                              * ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 4)    &
         &                                - flux_var( 1:jpi -2, 2:jpj -1, jk   , 4)    &
         !  d_y(a_22)
         &                                + flux_var( 2:jpi -1, 3:jpj -1, jk   , 5)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -2, jk   , 5) )  &
         !  d_z(a_32)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 6)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 6) )  &
         &                              /      e3v_n( 2:jpi -1, 2:jpj -1, jk   )
         !
         !
         !
         pcw(2:jpi -1, 2:jpj -1, jk   ) =        pcw( 2:jpi -1, 2:jpj -1, jk   )       &
         &                              + r1_e1e2e3w( 2:jpi -1, 2:jpj -1, jk   )       &
         !  d_x(a_13)
         &                              * ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 7)    &
         &                                - flux_var( 1:jpi -2, 2:jpj -1, jk   , 7)    &
         !  d_y(a_23)
         &                                + flux_var( 2:jpi -1, 2:jpj -1, jk   , 8)    &
         &                                - flux_var( 2:jpi -1, 1:jpj -2, jk   , 8) )

         pcw(2:jpi -1, 2:jpj -1, jk +1) =        pcw( 2:jpi -1, 2:jpj -1, jk +1)       &
         !  d_z(a_33)
         &                              + ( flux_var( 2:jpi -1, 2:jpj -1, jk   , 9)    &
         &                                - flux_var( 2:jpi -1, 2:jpj -1, jk +1, 9) )  &
         &                              *      wmask( 2:jpi -1, 2:jpj -1, jk +1)       &
         &                              /      e3w_n( 2:jpi -1, 2:jpj -1, jk +1)
         !
      END DO
      !
      ! Top boundary condition
      !
      pcw(2:jpi -1, 2:jpj -1, 1) =      pcw( 2:jpi -1, 2:jpj -1, 1)                &
      !  d_z(a_33)|_1
      &                          - flux_var( 2:jpi -1, 2:jpj -1, 1, 9)             &
      &                          *    wmask( 2:jpi -1, 2:jpj -1, 1)                &
      &                          /    e3w_n( 2:jpi -1, 2:jpj -1, 1)
      !
      ! Bottom boundary condition
      !
      pcw(:,:,jpk) = 0._wp
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_isd', pcu , 'U', -1., pcv , 'V', -1., pcw , 'T', -1.  )
      !
   END SUBROUTINE tlu_isd
   ! [tlu_isd]


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


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_norm_mod ***
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
   !!                
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
   ! @snippet this tlu_norm_mod
   ! [tlu_norm_mod]
   SUBROUTINE tlu_norm_mod ( kt , prt_noi )
      USE daymod
      !!------------------------------------------------------------------------
      !!                    ***  ROUTINE tlu_norm_mod  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method :    
      !!                
      !!                
      !!
      !! ** Action :    
      !!                
      !!                
      !!                
      !!                
      !!                
      !!                
      !!
      !!------------------------------------------------------------------------
      INTEGER                           , INTENT(in   ) ::   kt         ! ocean time-step index
      LOGICAL                           , INTENT(in   ) ::   prt_noi    ! ocean time-step index

      REAL(wp), DIMENSION(jpi,jpj,jpk)                  :: zwke_noi, zwke_vel
      REAL(wp), DIMENSION(mppsize)                      :: zw1, zw2
      REAL(wp), DIMENSION(12,mppsize)                   :: zw3
      REAL(wp)                                          :: tke_vel_sub
      REAL(wp)                                          :: tke_noi_sub
      REAL(wp)                                          :: tke_noi
      REAL(wp)                                          :: tke_vel
      REAL(wp), DIMENSION(12) ::    zmaxl,    zmaxg
      REAL(wp), DIMENSION(12) ::   zmeanl,   zmeang

      INTEGER                                           :: i, ierr
      !!------------------------------------------------------------------------
      !

!      zwke_noi = 0._wp
!      zwke_noi = ( unoi**2 ) * (spread(e1e2u,3,jpk)*e3u_n) + &
!               & ( vnoi**2 ) * (spread(e1e2v,3,jpk)*e3v_n) + &
!               & ( wnoi**2 ) * (spread(e1e2t,3,jpk)*e3w_n) 

!      tke_noi_sub = SUM( zwke_noi(2:jpim1, 2:jpjm1, 1:jpkm1) )

!      zwke_vel = 0._wp
!      zwke_vel = ( ubia_n**2 ) * (spread(e1e2u,3,jpk)*e3u_n) + &
!               & ( vbia_n**2 ) * (spread(e1e2v,3,jpk)*e3v_n) + &
!               & ( wbia_n**2 ) * (spread(e1e2t,3,jpk)*e3w_n)  

!      tke_vel_sub = SUM( zwke_vel(2:jpim1, 2:jpjm1, 1:jpkm1) )

      zmeanl(1) = SUM( unoi(2:jpim1, 2:jpjm1, 1:jpkm1) )
      zmeanl(2) = SUM( vnoi(2:jpim1, 2:jpjm1, 1:jpkm1) )
      zmeanl(3) = SUM( wnoi(2:jpim1, 2:jpjm1, 1:jpkm1) )

      zmeanl(4) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,1) )
      zmeanl(5) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,2) )
      zmeanl(6) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,3) )

      zmeanl(7) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,4) )
      zmeanl(8) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,5) )
      zmeanl(9) = SUM( var_ten(2:jpim1, 2:jpjm1, 1:jpkm1,6) )

      IF (ln_tlu_bia) THEN 
         zmeanl(10) = SUM( ubia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
         zmeanl(11) = SUM( vbia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
         zmeanl(12) = SUM( wbia_n(2:jpim1, 2:jpjm1, 1:jpkm1) )
      ELSE
         zmeanl(10) = 0._wp 
         zmeanl(11) = 0._wp 
         zmeanl(12) = 0._wp
      END IF

      call MPI_BARRIER(mpi_comm_oce, ierr)
!      call MPI_ALLGATHER( tke_noi_sub, 1, mpi_double_precision, zw1, 1, mpi_double_precision, mpi_comm_oce, ierr)
!      call MPI_ALLGATHER( tke_vel_sub, 1, mpi_double_precision, zw2, 1, mpi_double_precision, mpi_comm_oce, ierr)
      
      call MPI_ALLREDUCE(zmeanl, zmeang, 12, mpi_double_precision, MPI_SUM, mpi_comm_oce, ierr)
      call MPI_BARRIER(mpi_comm_oce, ierr)

      zmeang = zmeang / ( mppsize * jpim1 * jpjm1 * jpkm1 )
!      tke_noi = SUM(zw1)
!      tke_vel = SUM(zw2)

      KE_ratio_n = 0.2 !* SQRT( tke_vel / tke_noi )

      unoi = unoi! * KE_ratio_n
      vnoi = vnoi! * KE_ratio_n
      wnoi = wnoi! * KE_ratio_n

      var_ten(:,:,:,:) = var_ten(:,:,:,:)! * (KE_ratio_n**2) / (KE_ratio_b**2)
!      KE_ratio_b = KE_ratio_n
      IF (prt_noi) THEN
        zmaxl(1) = maxval(unoi); zmaxl(2) = maxval(vnoi); zmaxl(3) = maxval(wnoi);
        zmaxl(4) = maxval(var_ten(:,:,:,1)); zmaxl(5) = maxval(var_ten(:,:,:,2)); zmaxl(6) = maxval(var_ten(:,:,:,3));
        zmaxl(7) = maxval(var_ten(:,:,:,4)); zmaxl(8) = maxval(var_ten(:,:,:,5)); zmaxl(9) = maxval(var_ten(:,:,:,6));
        IF (ln_tlu_bia) THEN 
	   zmaxl(10) = maxval(ubia_n); zmaxl(11) = maxval(vbia_n); zmaxl(12) = maxval(wbia_n);
	END IF
    
        call MPI_BARRIER(mpi_comm_oce, ierr)
        call MPI_REDUCE(zmaxl, zmaxg, 12, mpi_double_precision, MPI_MAX, 0, mpi_comm_oce, ierr)

        IF (lwp .AND. mod(kt,500) .eq. 0) THEN
            print *, ""
            print '(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i3.3)', &
                  &    '======================>> time-step =', kt,'      DATE Y/M/D = ', nyear, '/', nmonth, '/', nday, '      nday_year = ', nday_year
            print *, '--------------------------------------------------------------------------------------------- Maximum values'
            print *, '     Max noise: ',  zmaxg( 1), zmaxg( 2), zmaxg( 3)
            print *, ' Max trace var: ',  zmaxg( 4), zmaxg( 5), zmaxg( 6)
            print *, ' Max cross var: ',  zmaxg( 7), zmaxg( 8), zmaxg( 9)
            print *, '      Max bias: ',  zmaxg(10), zmaxg(11), zmaxg(12)
            print *, '------------------------------------------------------------------------------------------------ Mean values'
            print *, '    Mean noise: ', zmeang( 1),zmeang( 2),zmeang( 3)
            print *, 'Mean trace var: ', zmeang( 4),zmeang( 5),zmeang( 6)
            print *, 'Mean cross var: ', zmeang( 7),zmeang( 8),zmeang( 9)
            print *, '     Mean bias: ', zmeang(10),zmeang(11),zmeang(12)
            print *, '------------------------------------------------------------------------------------------ Ratio of energies'
            print *, 'KEbias/KEnoise: '!, SQRT( tke_vel / tke_noi )
            print *, ""
        ENDIF
      ENDIF

   END SUBROUTINE tlu_norm_mod
   ! [tlu_norm_mod]

END MODULE tlunoi

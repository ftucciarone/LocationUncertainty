MODULE tluadvDYN
   !!==============================================================================
   !!                       ***  MODULE  tludynadv  ***
   !! Ocean active tracers:  advection scheme control
   !!==============================================================================
   !! History :  1.0  !  2006-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option 
   !!            4.0  !  2017-07  (G. Madec)  add a linear dynamics option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv      : compute the momentum advection trend 
   !!   dyn_adv_init : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce
   USE dom_oce         ! ocean space and time domain
   USE dynkeg          ! kinetic energy gradient          (dyn_keg      routine)
   USE dynzad          ! vertical advection               (dyn_zad      routine)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE tlu
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS
   REAL(wp), PARAMETER :: gamma2 = 1._wp/32._wp ! =0   2nd order  ; =1/32 4th order centred

   PUBLIC tlu_dyn_adv       ! routine called by step module
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynadv.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tlu_dyn_adv( kt )
      USE dynadv, ONLY : n_dynadv, nn_dynkeg, np_VEC_c2, np_FLX_c2, np_FLX_ubs
      USE oce,    ONLY : un, vn, wn                                ! ocean  
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv  ***
      !!                
      !! ** Purpose :   compute the ocean momentum advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following n_dynadv
      !!
      !!      NB: in flux form advection (ln_dynadv_cen2 or ln_dynadv_ubs=T) 
      !!      a metric term is add to the coriolis term while in vector form 
      !!      it is the relative vorticity which is added to coriolis term
      !!      (see dynvor module).
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  uadv, vadv, wadv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  udif, vdif
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'tlu_dyn_adv' )
      !
      SELECT CASE( n_dynadv )    !==  compute advection trend and add it to general trend  ==!
      CASE( np_VEC_c2  )     
         CALL dyn_keg     ( kt, nn_dynkeg )    ! vector form : horizontal gradient of kinetic energy
         CALL dyn_zad     ( kt )               ! vector form : vertical advection
      CASE( np_FLX_c2  )
         uadv(:,:,:) = un(:,:,:) - uisd_n(:,:,:) + unoi(:,:,:)
         vadv(:,:,:) = vn(:,:,:) - visd_n(:,:,:) + vnoi(:,:,:)
         wadv(:,:,:) = wn(:,:,:) - wisd_n(:,:,:) + wnoi(:,:,:)
         IF ( ln_tlu_bia ) THEN
            uadv = uadv(:,:,:) + biaSign * ubia_n(:,:,:)
            vadv = vadv(:,:,:) + biaSign * vbia_n(:,:,:)
            wadv = wadv(:,:,:) + biaSign * wbia_n(:,:,:)
         END IF
         CALL tlu_adv_cen2( kt, uadv, vadv, wadv, un, vn )               ! 2nd order centered scheme
      CASE( np_FLX_ubs )  
         uadv(:,:,:) = un(:,:,:) - uisd_n(:,:,:) + unoi(:,:,:)
         vadv(:,:,:) = vn(:,:,:) - visd_n(:,:,:) + vnoi(:,:,:)
         wadv(:,:,:) = wn(:,:,:) - wisd_n(:,:,:) + wnoi(:,:,:)
         IF ( ln_tlu_bia ) THEN
            uadv = uadv(:,:,:) + biaSign * ubia_n(:,:,:)
            vadv = vadv(:,:,:) + biaSign * vbia_n(:,:,:)
            wadv = wadv(:,:,:) + biaSign * wbia_n(:,:,:)
         END IF
         CALL tlu_adv_ubs ( kt, uadv, vadv, wadv )               ! 3rd order UBS      scheme (UP3)
      END SELECT
      !
      udif = 0._wp
      vdif = 0._wp
      !
      CALL tlu_dyn_nonCons( kt, 'U',            un, udif,  unoi, vnoi, wnoi, 1)
      CALL tlu_dyn_nonCons( kt, 'U', 0.5_wp * udif,   ua,  unoi, vnoi, wnoi, 1)
      !
      CALL tlu_dyn_nonCons( kt, 'V',            vn, vdif,  unoi, vnoi, wnoi, 2)
      CALL tlu_dyn_nonCons( kt, 'V', 0.5_wp * vdif,   va,  unoi, vnoi, wnoi, 2)
      !
      IF( ln_timing )   CALL timing_stop( 'tlu_dyn_adv' )
      !
   END SUBROUTINE tlu_dyn_adv


   SUBROUTINE tlu_dyn_nonCons( kt, cdtype, pcn, pca,  uadv, vadv, wadv, kcmp)
      USE tlusbcDYN
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tlu_adv_sto  ***
      !!
      !! ** Purpose :   Compute the stochastic advection velocity -0.5*div(a)'
      !!                for the given component.
      !!
      !! ** Method  :   Set pc to -0.5*divc
      !!                where divc = 1/(e1c e2c)*(di[e2c a1] + dj[e1c a2])
      !!                           + 1/(e3c)*dk[a3]
      !!                and a1, a2, a3 are the appropriate components of the
      !!                stochastic diffusion tensor var_ten.
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER                         , INTENT(in   ) ::   kt             ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kcmp           ! Define velocity type for calculation
      CHARACTER(len=1)                , INTENT(in   ) ::   cdtype         ! =U, V or W (component indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   uadv, vadv, wadv !The components advecting the velocity field
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pcn
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pca             ! modified advection term for this component
      !
      INTEGER  ::   ji, jj, jk, jn                               ! dummy loop indices
      REAL(wp) ::   zdiv                                         ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwe1, zwe2               ! workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwa1, zwa2, zwa3, zwpw   ! workspace
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_dyn_nonCons')   ! [NEMO] check
      !
      IF( kt == nit000 + dt_delay )  THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'tlu_dyn_nonCons : modified adv. for ', cdtype
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( zwa1(jpi,jpj,jpk) , zwa2(jpi,jpj,jpk) , zwa3(jpi,jpj,jpk))
      ALLOCATE( zwe1(jpi,jpj,jpk) , zwe2(jpi,jpj,jpk) )
      !
      zwa1 = 0._wp
      zwa2 = 0._wp
      zwa3 = 0._wp

      SELECT CASE( kcmp ) 
      CASE ( np_ucmp )
         !
         zwe1 = spread(r1_e1t,3,jpkm1)
         zwe2 = spread(r1_e2f,3,jpkm1)       
         !
         zwa1( 2:jpi  , 2:jpjm1, 1:jpkm1) = ( uadv(2:jpi  , 2:jpjm1, 1:jpkm1 )   &
         &                                  + uadv(1:jpi-1, 2:jpjm1, 1:jpkm1 ) ) &
         &                                * (  pcn(2:jpi  , 2:jpjm1, 1:jpkm1 )   &
         &                                  -  pcn(1:jpi-1, 2:jpjm1, 1:jpkm1 ) ) &
         &                                *   zwe1(2:jpi  , 2:jpjm1, 1:jpkm1 ) * 0.5_wp
         !
         zwa2( 2:jpim1, 1:jpjm1, 1:jpkm1) = ( vadv(2:jpi  , 1:jpjm1, 1:jpkm1 )   &
         &                                  + vadv(1:jpi-1, 1:jpjm1, 1:jpkm1 ) ) &
         &                                * (  pcn(2:jpim1, 2:jpj  , 1:jpkm1 )   &
         &                                  -  pcn(2:jpim1, 1:jpj-1, 1:jpkm1 ) ) &
         &                                *   zwe2(2:jpim1, 1:jpjm1, 1:jpkm1 ) * 0.5_wp
         !
         zwa3( 2:jpim1, 2:jpjm1, 2:jpkm1) = ( wadv(2:jpi-1, 1:jpjm1, 2:jpkm1 )   &
         &                                  + wadv(3:jpi  , 1:jpjm1, 2:jpkm1 ) ) &
         &                                * (  pcn(2:jpim1, 2:jpjm1, 1:jpkm1 )   &
         &                                  -  pcn(2:jpim1, 2:jpjm1, 2:jpk   ) ) &
         &                                / e3uw_n(2:jpim1, 1:jpjm1, 2:jpkm1 )
         !
         !
         zwa3( :, :, jpk) = 0._wp
         zwa3( :, :, : ) = zwa3( :, :, :) * 0.5_wp
         !
         ! Surface boundary conditions
         CALL tlu_dynadv_sbc( kt, cdtype, kcmp, zwa3(:,:,1) ) 
         !
         pca( 2:jpim1, 2:jpjm1, 1:jpkm1 ) =   pca( 2:jpim1, 2:jpjm1, 1:jpkm1 )  &
         !
         &                                - ( zwa1( 3:jpi  , 2:jpjm1, 1:jpkm1 ) &
         &                                  + zwa1( 2:jpim1, 2:jpjm1, 1:jpkm1 ) &
         !
         &                                  + zwa2( 2:jpim1, 2:jpj-1, 1:jpkm1 ) &
         &                                  + zwa2( 2:jpim1, 1:jpj-2, 1:jpkm1 ) &
         !
         &                                  + zwa3( 2:jpim1, 2:jpjm1, 1:jpkm1 ) &
         &                                  + zwa3( 2:jpim1, 2:jpjm1, 2:jpk   ) )  * 0.5_wp
         !
      CASE ( np_vcmp )
         !
         zwe1 = spread(r1_e1f,3,jpkm1)
         zwe2 = spread(r1_e2t,3,jpkm1)    
         !
         zwa1( 2:jpi  , 2:jpjm1, 1:jpkm1) = ( uadv(1:jpim1, 2:jpj  , 1:jpkm1 )   &
         &                                  + uadv(1:jpim1, 1:jpj-1, 1:jpkm1 ) ) &
         &                                * (  pcn(3:jpi  , 2:jpjm1, 1:jpkm1 )   &
         &                                  -  pcn(2:jpi-1, 2:jpjm1, 1:jpkm1 ) ) &
         &                                *   zwe1(2:jpi  , 2:jpjm1, 1:jpkm1 ) * 0.5_wp
         !
         zwa2( 2:jpim1, 1:jpjm1, 1:jpkm1) = ( vadv(2:jpim1, 2:jpj-1, 1:jpkm1 )   &
         &                                  + vadv(2:jpim1, 1:jpj-2, 1:jpkm1 ) ) &
         &                                * (  pcn(2:jpim1, 2:jpj-1, 1:jpkm1 )   &
         &                                  -  pcn(2:jpim1, 1:jpj-2, 1:jpkm1 ) ) &
         &                                *   zwe2(2:jpim1, 1:jpjm1, 1:jpkm1 ) * 0.5_wp
         !
         zwa3( 2:jpim1, 2:jpjm1, 2:jpkm1) = ( wadv(2:jpim1, 3:jpj  , 2:jpkm1 )   &
         &                                  + wadv(2:jpim1, 2:jpjm1, 2:jpkm1 ) ) &
         &                                * (  pcn(2:jpim1, 2:jpjm1, 1:jpkm1 )   &
         &                                  -  pcn(2:jpim1, 2:jpjm1, 2:jpk   ) ) &
         &                                / e3vw_n(2:jpim1, 1:jpjm1, 2:jpkm1 )
         !
         zwa3( :, :, jpk) = 0._wp
         zwa3( :, :, : ) = zwa3( :, :, :) * 0.5_wp
         !
         ! Surface boundary conditions
         CALL tlu_dynadv_sbc( kt, cdtype, kcmp, zwa3(:,:,1) ) 
         !
         pca( 2:jpim1, 2:jpjm1, 1:jpkm1 ) =   pca( 2:jpim1, 2:jpjm1, 1:jpkm1 )  &
         !
         &                                - ( zwa1( 2:jpi-1, 2:jpjm1, 1:jpkm1 ) &
         &                                  + zwa1( 1:jpi-2, 2:jpjm1, 1:jpkm1 ) &
         !
         &                                  + zwa2( 2:jpim1, 3:jpj  , 1:jpkm1 ) &
         &                                  + zwa2( 2:jpim1, 2:jpj-1, 1:jpkm1 ) &
         !
         &                                  + zwa3( 2:jpim1, 2:jpjm1, 1:jpkm1 ) &
         &                                  + zwa3( 2:jpim1, 2:jpjm1, 2:jpk   ) )  * 0.5_wp
         !
      CASE DEFAULT                                             ! error
         CALL ctl_stop('STOP','tlu_dyn_nonCons: wrong value for kcmp'  )
      END SELECT
      !
      DEALLOCATE(zwa1, zwa2, zwa3)
      DEALLOCATE(zwe1, zwe2)
      !
      IF( ln_timing ) CALL timing_stop('tlu_dyn_nonCons')   ! [NEMO] check
      !
   END SUBROUTINE tlu_dyn_nonCons

   SUBROUTINE tlu_adv_cen2( kt, uadv, vadv, wadv, u, v)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_cen2  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!
      !! ** Action  :   (ua,va) updated with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   )  ::  uadv, vadv, wadv
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   )  ::     u,    v
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfu_t, zfu_f, zfu_uw, zfu
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfv_t, zfv_f, zfv_vw, zfv, zfw
      !!----------------------------------------------------------------------
      !
      IF( (kt == nit000 + dt_delay) .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tlu_adv_cen2 : 2nd order flux form momentum advection'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF
      !
      !                             !==  Horizontal advection  ==!
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u_n(:,:,jk) * ( uadv(:,:,jk) ) 
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v_n(:,:,jk) * ( vadv(:,:,jk) )
         DO jj = 1, jpjm1                 ! horizontal momentum fluxes (at T- and F-point)
            DO ji = 1, fs_jpim1   ! vector opt.
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( u(ji,jj,jk) + u(ji+1,jj  ,jk) )
               zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( u(ji,jj,jk) + u(ji  ,jj+1,jk) )
               zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( v(ji,jj,jk) + v(ji+1,jj  ,jk) )
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( v(ji,jj,jk) + v(ji  ,jj+1,jk) )
            END DO
         END DO
         DO jj = 2, jpjm1                 ! divergence of horizontal momentum fluxes
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
                  &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
                  &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN           ! trends: send trend to trddyn for diagnostic
         zfu_uw(:,:,:) = ua(:,:,:) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt )
         zfu_t(:,:,:) = ua(:,:,:)
         zfv_t(:,:,:) = va(:,:,:)
      ENDIF
      !
      !                             !==  Vertical advection  ==!
      !
      DO jj = 2, jpjm1                    ! surface/bottom advective fluxes set to zero
         DO ji = fs_2, fs_jpim1
            zfu_uw(ji,jj,jpk) = 0._wp   ;   zfv_vw(ji,jj,jpk) = 0._wp
            zfu_uw(ji,jj, 1 ) = 0._wp   ;   zfv_vw(ji,jj, 1 ) = 0._wp
         END DO
      END DO
      IF( ln_linssh ) THEN                ! linear free surface: advection through the surface
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wadv(ji,jj,1) + e1e2t(ji+1,jj) * wadv(ji+1,jj,1) ) * u(ji,jj,1)
               zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wadv(ji,jj,1) + e1e2t(ji,jj+1) * wadv(ji,jj+1,1) ) * v(ji,jj,1)
            END DO
         END DO
      ENDIF
      DO jk = 2, jpkm1                    ! interior advective fluxes
         DO jj = 2, jpj                       ! 1/4 * Vertical transport
            DO ji = 2, jpi
               zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * wadv(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji+1,jj  ,jk) ) * ( u(ji,jj,jk) + u(ji,jj,jk-1) )
               zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji  ,jj+1,jk) ) * ( v(ji,jj,jk) + v(ji,jj,jk-1) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1                    ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                 ! trends: send trend to trddyn for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                   ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' cen2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE tlu_adv_cen2

   SUBROUTINE tlu_adv_ubs( kt, uadv, vadv, wadv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_ubs  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   The scheme is the one implemeted in ROMS. It depends 
      !!      on two parameter gamma1 and gamma2. The former control the 
      !!      upstream baised part of the scheme and the later the centred 
      !!      part:     gamma1 = 0    pure centered  (no diffusive part)
      !!                       = 1/4  Quick scheme
      !!                       = 1/3  3rd order Upstream biased scheme
      !!                gamma2 = 0    2nd order finite differencing 
      !!                       = 1/32 4th order finite differencing
      !!      For stability reasons, the first term of the fluxes which cor-
      !!      responds to a second order centered scheme is evaluated using  
      !!      the now velocity (centered in time) while the second term which  
      !!      is the diffusive part of the scheme, is evaluated using the 
      !!      before velocity (forward in time). 
      !!      Default value (hard coded in the begining of the module) are 
      !!      gamma1=1/3 and gamma2=1/32.
      !!
      !! ** Action : - (ua,va) updated with the 3D advective momentum trends
      !!
      !! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling. 
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   )  ::  uadv, vadv, wadv
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfu_t, zfu_f, zfu_uw, zfu
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfv_t, zfv_f, zfv_vw, zfv, zfw
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlu_uu, zlu_uv
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlv_vv, zlv_vu
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 + dt_delay ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_adv_ubs : UBS flux form momentum advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      zfu_t(:,:,:) = 0._wp
      zfv_t(:,:,:) = 0._wp
      zfu_f(:,:,:) = 0._wp
      zfv_f(:,:,:) = 0._wp
      !
      zlu_uu(:,:,:,:) = 0._wp
      zlv_vv(:,:,:,:) = 0._wp 
      zlu_uv(:,:,:,:) = 0._wp 
      zlv_vu(:,:,:,:) = 0._wp 
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF
      !                                      ! =========================== !
      DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
         !                                   ! =========================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = e2u(:,:) * e3u_n(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = e1v(:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
         !            
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zlu_uu(ji,jj,jk,1) = ( ub (ji+1,jj  ,jk) - 2.*ub (ji,jj,jk) + ub (ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,1) = ( vb (ji  ,jj+1,jk) - 2.*vb (ji,jj,jk) + vb (ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,1) = ( ub (ji  ,jj+1,jk) - ub (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( ub (ji  ,jj  ,jk) - ub (ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,1) = ( vb (ji+1,jj  ,jk) - vb (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( vb (ji  ,jj  ,jk) - vb (ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
               !
               zlu_uu(ji,jj,jk,2) = ( zfu(ji+1,jj  ,jk) - 2.*zfu(ji,jj,jk) + zfu(ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,2) = ( zfv(ji  ,jj+1,jk) - 2.*zfv(ji,jj,jk) + zfv(ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,2) = ( zfu(ji  ,jj+1,jk) - zfu(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfu(ji  ,jj  ,jk) - zfu(ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,2) = ( zfv(ji+1,jj  ,jk) - zfv(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfv(ji  ,jj  ,jk) - zfv(ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( 'tlu_adv_ubs', zlu_uu(:,:,:,1), 'U', 1. , zlu_uv(:,:,:,1), 'U', 1.,  &
                      &   zlu_uu(:,:,:,2), 'U', 1. , zlu_uv(:,:,:,2), 'U', 1.,  & 
                      &   zlv_vv(:,:,:,1), 'V', 1. , zlv_vu(:,:,:,1), 'V', 1.,  &
                      &   zlv_vv(:,:,:,2), 'V', 1. , zlv_vu(:,:,:,2), 'V', 1.   )
      !
      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u_n(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
         !
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, fs_jpim1   ! vector opt.
               zui = ( un(ji,jj,jk) + un(ji+1,jj  ,jk) )
               zvj = ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) )
               !
               IF( zui > 0 ) THEN   ;   zl_u = zlu_uu(ji  ,jj,jk,1)
               ELSE                 ;   zl_u = zlu_uu(ji+1,jj,jk,1)
               ENDIF
               IF( zvj > 0 ) THEN   ;   zl_v = zlv_vv(ji,jj  ,jk,1)
               ELSE                 ;   zl_v = zlv_vv(ji,jj+1,jk,1)
               ENDIF
               !
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,jk,2) + zlu_uu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,jk,2) + zlv_vv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
               !
               zfuj = ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) )
               zfvi = ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) )
               IF( zfuj > 0 ) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,jk,1)
               ELSE                  ;    zl_v = zlv_vu( ji+1,jj,jk,1)
               ENDIF
               IF( zfvi > 0 ) THEN   ;    zl_u = zlu_uv( ji,jj  ,jk,1)
               ELSE                  ;    zl_u = zlu_uv( ji,jj+1,jk,1)
               ENDIF
               !
               zfv_f(ji  ,jj  ,jk) = ( zfvi - gamma2 * ( zlv_vu(ji,jj,jk,2) + zlv_vu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) - gamma1 * zl_u )
               zfu_f(ji  ,jj  ,jk) = ( zfuj - gamma2 * ( zlu_uv(ji,jj,jk,2) + zlu_uv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) - gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
                  &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
                  &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      IF( l_trddyn ) THEN                          ! trends: send trends to trddyn for diagnostic
         zfu_uw(:,:,:) = ua(:,:,:) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt )
         zfu_t(:,:,:) = ua(:,:,:)
         zfv_t(:,:,:) = va(:,:,:)
      ENDIF
      !                                      ! ==================== !
      !                                      !  Vertical advection  !
      !                                      ! ==================== !
      DO jj = 2, jpjm1                             ! surface/bottom advective fluxes set to zero                  
         DO ji = fs_2, fs_jpim1
            zfu_uw(ji,jj,jpk) = 0._wp
            zfv_vw(ji,jj,jpk) = 0._wp
            zfu_uw(ji,jj, 1 ) = 0._wp
            zfv_vw(ji,jj, 1 ) = 0._wp
         END DO
      END DO
      IF( ln_linssh ) THEN                         ! constant volume : advection through the surface
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji+1,jj) * wn(ji+1,jj,1) ) * un(ji,jj,1)
               zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji,jj+1) * wn(ji,jj+1,1) ) * vn(ji,jj,1)
            END DO
         END DO
      ENDIF
      DO jk = 2, jpkm1                          ! interior fluxes
         DO jj = 2, jpj
            DO ji = 2, jpi
               zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * wn(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji+1,jj,jk) ) * ( un(ji,jj,jk) + un(ji,jj,jk-1) )
               zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1                          ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) =  ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) =  va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                       ! save the vertical advection trend for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                         ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ubs2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE tlu_adv_ubs

END MODULE tluadvDYN

!!-----------------------------------------------------------------------------------------------------------------------------------
!! 
!!                               MODULE: tluhdifTRA
!!   Transport under Location Uncertainty: stochastic parametrization of N-S
!!   equations.
!!
!> @authors
!>                P. Derian, P. Chandramouli, F. Tucciarone
!!
!> @version
!!                2017 -  ( P. DERIAN       )  Original code <br>
!!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
!!                2022 -  ( F. TUCCIARONE   )  Ongoing work and documentation
!! 
!> @brief         Transport under Location Uncertainty: TRACER diffusion.
!! 
!! @par           Procedure specifics  
!!                Computes all the terms of the Stochastic Diffusion of the Tracer, with only exception
!!                of the purely vertical diffusive term, computed with an implicit (in time) procedure.    
!>                \f{align*}{
!!                {\color{gray} \mathrm{d}_{t} \theta } + 
!!                {\color{gray} \nabla\cdot\left\lbrace \left[\boldsymbol{u}-\boldsymbol{u}_{s}+\boldsymbol{\mu}_{t}\right)\theta\mathrm{d}t
!!                                 +  \boldsymbol{\sigma}\mathrm{d}\mathbf{B}_{t}\theta \right\rbrace } - 
!!                {\color{black} \dfrac{1}{2}\nabla\cdot\left( \boldsymbol{a}\nabla\theta \right)\mathrm{d}t } =
!!                {\color{gray} \mathrm{D}^{\theta} } +
!!                {\color{gray} \mathrm{F}^{\theta} } 
!!                       
!!                \f}
!!                The Stochastic Diffusion Term is computed with a finite volume approach, defining 
!!                diffusive fluxes entering and exiting the volume cell and computing the divergence
!!                using the balance of fluxes provided by the divergence theorem:
!!                \f{align*}{
!!                \int_{V}\nabla\cdot \boldsymbol{F}\,\mathrm{d}V = 
!!                       \int_{\partial V}\boldsymbol{F}\cdot\boldsymbol{n}\,\mathrm{d}s.
!!                \f}  
!!                The diffusive flux is defined as
!!                \f{align*}{
!!                \boldsymbol{F}_{T} = \left(
!!                \begin{array}{c}
!!                      F_{T,x} \\
!!                      F_{T,y} \\
!!                      F_{T,z}+a_{33}\partial_{z}\theta
!!                      \end{array}
!!                      \right) = \left(
!!                \begin{array}{c}
!!                      \left( a_{11}\partial_{x}\theta +a_{12}\partial_{y}\theta +a_{13}\partial_{z}\theta \right)\mathrm{d}y\mathrm{d}z \\
!!                      \left( a_{21}\partial_{x}\theta +a_{22}\partial_{y}\theta +a_{23}\partial_{z}\theta \right)\mathrm{d}x\mathrm{d}z \\
!!                      \left( a_{31}\partial_{x}\theta +a_{32}\partial_{y}\theta +a_{33}\partial_{z}\theta \right)\mathrm{d}x\mathrm{d}y
!!                      \end{array}
!!                      \right)
!!                \f}
!!                of which the terms @f$ F_{x} @f$ and @f$ F_{y} @f$ are computed in \ref tlu_trahhdiff,
!!                while the term @f$ F_{z} @f$ (which does not contain the purely vertical diffusion term) 
!!                is computed in \ref tlu_trahzdiff,

!!                \f{align*}{
!!                \begin{array}{c}
!!                      a_{11}\partial_{x}\theta +a_{12}\partial_{y}\theta +a_{13}\partial_{z}\theta \\
!!                      a_{21}\partial_{x}\theta +a_{22}\partial_{y}\theta +a_{23}\partial_{z}\theta \\
!!                      a_{31}\partial_{x}\theta +a_{32}\partial_{y}\theta +a_{33}\partial_{z}\theta
!!                      \end{array}
!!                      \right)
!!                \f}
!> @par           Code specifics
!!                This module is called by tlu_bcdyn in tlu.f90. The terms @f$ F_{x} @f$ and @f$ F_{y} @f$ are computed in \ref tlu_trahhdiff,
!!                while the term @f$ F_{z} @f$ (which does not contain the purely vertical diffusion term) is computed in \ref tlu_trahzdiff,
!!                so that in terms of NEMO operators it becomes
!!                \f{align*}{
!!                \nabla\cdot\boldsymbol{F}_{T} = \underbrace{\dfrac{1}{e_{1t}e_{2t}e_{3t}}\left\lbrace 
!!                       \delta_{i}\left[ F_{T,x} \right]
!!                     + \delta_{j}\left[ F_{T,y} \right] 
!!                       \right\rbrace}_{\texttt{tlu\_trahhdiff}}
!!                     + \underbrace{\dfrac{1}{e_{3t}}\delta_{k}\left[ F_{T,z} \right]}_{\texttt{tlu\_trahzdiff}}
!!                \f}
 
!!                
!!
!> @param[in]     ln_tlu: logical switch for location uncertainty
!! @param[in]     jpi: the first dimension of the spatial arrays
!! @param[in]     jpj: the first dimension of the spatial arrays
!! @param[in]     jpk: the first dimension of the spatial arrays
!!
!> @par           Other modules dependencies
!> @snippet       this mod_dep
!!
!> @par           Public routines
!> @snippet       this public_sub
!!
!> @todo          Write error messages to appropriate unit
!!
!!-----------------------------------------------------------------------------------------------------------------------------------
MODULE tluhdifTRA
   !
   ! [mod_dep]
   USE par_kind         ! data types defined in par_kind module
   USE in_out_manager   ! I/O manager
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
   PRIVATE         
   !
   ! [public_sub]
   PUBLIC tlu_trahhdiff      ! Called by tlu_bcdyn in tlu.f90
   PUBLIC tlu_trahzdiff      ! Called by tlu_bcdyn in tlu.f90
   ! [public_sub]
   !
   !
   INCLUDE 'mpif.h'
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !
CONTAINS


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trahhdiff ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes the Horizontal components of the Stochastic Diffusion
   !!                of the Tracer 
   !!
   !> @par           Procedure specifics
   !>                Computes the horizontal terms of the stochastic diffusion, using a finite volumes approach
   !!                for the diffusive fluxes @f$ F_{T,x} @f$ and @f$ F_{T,y} @f$ centred on a @f$ T- @f$ point
   !!                \f{align*}{
   !!                \nabla\cdot \boldsymbol{F}_{T} &= \dfrac{1}{e_{1t}e_{2t}e_{3t}}\delta_{i}\left[ F_{T,x} \right] 
   !!                                                              + \dfrac{1}{e_{1t}e_{2t}e_{3t}}\delta_{j}\left[ F_{T,y} \right]  
   !!                \f} 
   !!                \f{align*}{
   !!                F_{T,x} &= \overline{a_{11}^{T}}^{\,i}\partial_{x}T \,e_{2u}e_{3u}^{n} 
   !!                         + \overline{a_{12}^{f}\overline{\partial_{y}T \,e_{2v}e_{3v}^{n}}^{\,i}}^{ \,j} 
   !!                         + \overline{a_{13}^{uw}\overline{\partial_{z}T \,e_{2t}e_{3w}^{n}}^{\,i}}^{ \,k} \\
   !!                F_{T,y} &= \overline{a_{21}^{f}\overline{\partial_{x}T \,e_{1u}e_{3u}^{n}}^{\,j}}^{ \,i} 
   !!                         + \overline{a_{22}^{T}}^{ \,j}\partial_{y}T \,e_{1v}e_{3v}^{n} 
   !!                         + \overline{a_{23}^{vw}	\overline{\partial_{z}T \,e_{1u}e_{3uw}^{n}}^{ \,j}}^{ \,k} 
   !!                \f}   
   !! 
   !> @par           Code specifics               
   !!                The following naming conventions are employed:
   !!                \f{align*}{
   !!                \texttt{ztuu} & = \overline{a_{11}^{T}}^{\,i}\partial_{x}T \,e_{2u}e_{3u}^{n} \qquad 
   !!                \texttt{ztuv}   = \overline{a_{12}^{f}\overline{\partial_{y}T \,e_{2v}e_{3v}^{n}}^{\,i}}^{ \,j} \qquad 
   !!                \texttt{ztuw}   = \overline{a_{13}^{uw}\overline{\partial_{z}T \,e_{2t}e_{3w}^{n}}^{\,i}}^{ \,k}   \\
   !!                \texttt{ztvu} & = \overline{a_{21}^{f}\overline{\partial_{x}T \,e_{1u}e_{3u}^{n}}^{\,j}}^{ \,i} \qquad  
   !!                \texttt{ztvv}   = \overline{a_{22}^{T}}^{ \,j}\partial_{y}T \,e_{1v}e_{3v}^{n} \qquad  
   !!                \texttt{ztvw}   = \overline{a_{23}^{vw}	\overline{\partial_{z}T \,e_{1u}e_{3uw}^{n}}^{ \,j}}^{ \,k} 
   !!                \f} 
   !!
   !!
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     ptb: Tracer before field (Euler scheme)
   !> @param[inout]  pta: Tracer trend to be updated
   !> @param[in]     r_Scale: Scaling factor for diffusion
   !! 
   !! @result        Update (pta) with the diffusion term
   !!
   !> @par           Diagnostic
   !!                With Homogeneous Neumann boundary conditions the integrated divergence of the diffusive fluxes
   !!                is zero.        
   !!                \f{align*}{
   !!                \int_{V}\nabla\cdot \boldsymbol{F}_{T}\,\mathrm{d}V &=
   !!                \int_{\partial V}\boldsymbol{F}\cdot\boldsymbol{n}\,\mathrm{d}s= 0
   !!                \f} 
   !!                The computational expression of this equation is        
   !!                \f{align*}{
   !!                \int_{V}\nabla\cdot \boldsymbol{F}_{T}\,\mathrm{d}V &=
   !!                \sum_{i,j,k} \nabla\cdot \boldsymbol{F}_{T} \mathrm{d} x \mathrm{d} y \mathrm{d} z = 
   !!                \sum_{i,j,k} \left\lbrace \dfrac{1}{e_{1t}e_{2t}e_{3t}}\delta_{i}\left[ F_{T,x} \right] 
   !!                             + \dfrac{1}{e_{2t}e_{2t}e_{3t}}\delta_{j}\left[ F_{T,y} \right] 
   !!                  \right\rbrace e_{1t}e_{2t}e_{3t} = 0 
   !!                \f} 
   !  @note          
   !  @todo              
   !!
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_trahhdiff
   ! [tlu_trahhdiff]
   SUBROUTINE tlu_trahhdiff( kt, ptb, pta, r_Scale)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      REAL(wp), DIMENSION(jpts)            , INTENT(in   ) ::   r_Scale               ! Schmidt number
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw0                   ! accumulation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw1,zw2               ! accumulation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw3,zw4,zw5,zw6       ! accumulation array
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var11             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var12             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var21             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var22             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztuu                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztuv                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztuw                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztvu                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztvv                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztvw                  ! workspace array
      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia11 = ndiffidx(1,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia12 = ndiffidx(1,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia13 = ndiffidx(1,3)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia21 = ndiffidx(2,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia22 = ndiffidx(2,2)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia23 = ndiffidx(2,3)  ! Mapping the indeces of a

      REAL(wp), DIMENSION(4) :: zmaxl, zmaxg
      REAL(wp), DIMENSION(mppsize) :: ztau
      REAL(wp), DIMENSION(2)             ::   swl, swg
      REAL(wp), DIMENSION(2)             ::   mxl, mxg
      INTEGER  ::   ierr, norm_typ                   ! dummy loop argument
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_trahhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trahhdiff : Horizontal diffusion of TRACER '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF


      ALLOCATE( int_var11(jpi,jpj,jpk), &
              & int_var22(jpi,jpj,jpk)  )

      ALLOCATE( ztuu(jpi,jpj,jpk), &
              & ztuv(jpi,jpj,jpk), &
              & ztuw(jpi,jpj,jpk), &
              & ztvu(jpi,jpj,jpk), &
              & ztvv(jpi,jpj,jpk), &
              & ztvw(jpi,jpj,jpk)  )

      ALLOCATE( zw0(jpi,jpj,jpk,2), &
              & zw1(jpi,jpj,jpk,2), &
              & zw2(jpi,jpj,jpk,2), &
              & zw3(jpi,jpj,jpk,2), &
              & zw4(jpi,jpj,jpk,2), &
              & zw5(jpi,jpj,jpk,2), &
              & zw6(jpi,jpj,jpk,2)   )
      !
      ! Inizialization of accumulation arrays
      !
      zw0 = 0._wp
      zw1 = 0._wp
      zw2 = 0._wp
      zw3 = 0._wp
      zw4 = 0._wp
      zw5 = 0._wp
      zw6 = 0._wp
      !
      ! Inizialization of local and global Diagnostic arrays
      !  
      swl = 0._wp
      swg = 0._wp
      mxl = 0._wp
      mxg = 0._wp
      !
      DO jk = 1, jpkm1              !==  Interpolation of variance tensor  ==!
         DO jj = 2, jpjm1
            DO ji =  2, jpim1
               int_var11(ji,jj,jk) = ( var_ten(ji+1,jj  ,jk,ia11) + var_ten(ji  ,jj  ,jk,ia11) ) * 0.5_wp
               int_var22(ji,jj,jk) = ( var_ten(ji  ,jj+1,jk,ia22) + var_ten(ji  ,jj  ,jk,ia22) ) * 0.5_wp
            END DO
         END DO
      END DO 
      !
      ! Lateral boundary condition transfer across nodes
      !
      CALL lbc_lnk_multi( 'tlu_trahhdiff', int_var11 , 'U', 1.,   &
                        &                  int_var22 , 'V', 1.    )
      !
      !
      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !   
         ztuu = 0._wp 
         ztuv = 0._wp 
         ztuw = 0._wp
         ztvu = 0._wp 
         ztvv = 0._wp
         ztvw = 0._wp
         !
         !                               
         DO jk = 1, jpkm1              !==  First derivative (gradient)  ==!
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ! Derivative of the tracer
                  ztuu(ji,jj,jk) = r1_e1u(ji,jj) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztuv(ji,jj,jk) = r1_e2v(ji,jj) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztvu(ji,jj,jk) = r1_e1u(ji,jj) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztvv(ji,jj,jk) = r1_e2v(ji,jj) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
                  ! Diffusive flux
                  ztuu(ji,jj,jk) = ztuu(ji,jj,jk) * e2u(ji,jj) * e3u_n(ji,jj,jk) * int_var11(ji,jj,jk) 
                  ztuv(ji,jj,jk) = ztuv(ji,jj,jk) * e2f(ji,jj) * e3f_n(ji,jj,jk) *   var_ten(ji,jj,jk,ia12)
                  ztvu(ji,jj,jk) = ztvu(ji,jj,jk) * e1f(ji,jj) * e3f_n(ji,jj,jk) *   var_ten(ji,jj,jk,ia21)
                  ztvv(ji,jj,jk) = ztvv(ji,jj,jk) * e1v(ji,jj) * e3v_n(ji,jj,jk) * int_var22(ji,jj,jk)
                  ! Masking
                  ztuu(ji,jj,jk) = ztuu(ji,jj,jk) * umask(ji,jj,jk) 
                  ztuv(ji,jj,jk) = ztuv(ji,jj,jk) * fmask(ji,jj,jk) 
                  ztvu(ji,jj,jk) = ztvu(ji,jj,jk) * fmask(ji,jj,jk) 
                  ztvv(ji,jj,jk) = ztvv(ji,jj,jk) * vmask(ji,jj,jk) 
               END DO
            END DO
         END DO
         !
         ! Surface boundary condition
         ! ztuw(:,:,1) = 0._wp
         ! ztvw(:,:,1) = 0._wp
         !
         DO jk = 2, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ! Derivative of the tracer and its interpolation
                  ztuw(ji,jj,jk) = 0.5_wp * ( ( ptb(ji  ,jj  ,jk  ,jn) - ptb(ji  ,jj  ,jk+1,jn) ) /  e3w_n(ji  ,jj,jk)  &
                                 &          + ( ptb(ji+1,jj  ,jk  ,jn) - ptb(ji+1,jj  ,jk+1,jn) ) /  e3w_n(ji+1,jj,jk)  )
                  ztvw(ji,jj,jk) = 0.5_wp * ( ( ptb(ji  ,jj  ,jk  ,jn) - ptb(ji  ,jj  ,jk+1,jn) ) /  e3w_n(ji,jj  ,jk)  &
                                 &          + ( ptb(ji  ,jj+1,jk  ,jn) - ptb(ji  ,jj+1,jk+1,jn) ) /  e3w_n(ji,jj+1,jk)  )
                  ! Diffusive flux 
                  ztuw(ji,jj,jk) = e2u(ji,jj) * e3uw_n(ji,jj,jk) * var_ten(ji,jj,jk,ia13) * ztuw(ji,jj,jk)
                  ztvw(ji,jj,jk) = e1v(ji,jj) * e3vw_n(ji,jj,jk) * var_ten(ji,jj,jk,ia23) * ztvw(ji,jj,jk)
                  ! Masking
                  ztuw(ji,jj,jk) = ztuw(ji,jj,jk) * wumask(ji,jj,jk) 
                  ztvw(ji,jj,jk) = ztvw(ji,jj,jk) * wvmask(ji,jj,jk) 
               END DO
            END DO
         END DO

         !
         CALL lbc_lnk_multi( 'tlu_trahhdiff', ztuv, 'F', 1.,   &
         &                                    ztuw, 'U', 1.,   &
         &                                    ztvu, 'F', 1.,   &
         &                                    ztvw, 'V', 1.     )
         !
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !
                  ! Interpolation of extradiagonal term a_12\partial_y\theta
                  !
                  zw1(ji,jj,jk,jn) = 0.5_wp * (  ztuv(ji,jj,jk) + ztuv(ji  ,jj-1,jk) )   &
                  &                           * umask(ji,jj,jk)
                  !
                  ! Interpolation of extradiagonal term a_21\partial_x\theta
                  !
                  zw2(ji,jj,jk,jn) = 0.5_wp * (  ztvu(ji,jj,jk) + ztvu(ji-1,jj,jk) )   &
                  &                           * vmask(ji,jj,jk)
                  !
               END DO
            END DO
         END DO  
         !
         CALL lbc_lnk_multi( 'tlu_trahhdiff', zw1(:,:,:,jn), 'U', 1.,   &
                        &                     zw2(:,:,:,jn), 'V', 1.    )
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  !
                  ! Interpolation of extradiagonal term a_13\partial_z\theta
                  !
                  zw3(ji,jj,jk,jn) = 0.5_wp * (  ztuw(ji,jj,jk) + ztuw(ji,jj,jk+1) )   &
                  &                            * umask(ji,jj,jk)
                  !
                  ! Interpolation of extradiagonal term a_23\partial_z\theta
                  !
                  zw4(ji,jj,jk,jn) = 0.5_wp * (  ztvw(ji,jj,jk) + ztvw(ji,jj,jk+1) )   &
                  &                           * vmask(ji,jj,jk)                                                
               END DO
            END DO
         END DO
         !
         !
         !
         DO jk = 1, jpkm1              !==  Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !
                  ! Diagonal terms of the diffusion (0.5 factor because of 0.5\nabla\cdot(a\nabla\theta)
                  !
                  zw0(ji,jj,jk,jn) =  0.5_wp * (  ztuu(ji,jj,jk) - ztuu(ji-1,jj  ,jk)     &
                  &                            +  ztvv(ji,jj,jk) - ztvv(ji  ,jj-1,jk) )   &
                  &                            * tmask(ji,jj,jk) / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
                  !
                  ! Extra-diagonal terms (0.5 factor because of 0.5\nabla\cdot(a\nabla\theta)
                  !
                  zw5(ji,jj,jk,jn) = 0.5_wp * (  zw1(ji,jj,jk,jn) - zw1(ji-1,jj  ,jk,jn)     &
                  &                           +  zw2(ji,jj,jk,jn) - zw2(ji  ,jj-1,jk,jn) )   &
                  &                          * tmask(ji,jj,jk) / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
                  !
                  ! Vertical terms (0.5 factor because of 0.5\nabla\cdot(a\nabla\theta)
                  !
                  zw6(ji,jj,jk,jn) = 0.5_wp * (  zw3(ji,jj,jk,jn) - zw3(ji-1,jj  ,jk,jn)     &
                  &                           +  zw4(ji,jj,jk,jn) - zw4(ji  ,jj-1,jk,jn) )   &
                  &                          * tmask(ji,jj,jk) / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
                  !
                  ! Pointwise diffusion 
                  !
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ( zw0(ji,jj,jk,jn) + zw5(ji,jj,jk,jn) + zw6(ji,jj,jk,jn) ) / r_Scale(jn)
                  !
                  ! Diagnostic: swl =  \sum_{i,j,k} \nabla\cdot \boldsymbol{F}_{T} \mathrm{d} x \mathrm{d} y \mathrm{d} z
                  !
                  swl(jn) = swl(jn) + ( zw0(ji,jj,jk,jn) + zw5(ji,jj,jk,jn) + zw6(ji,jj,jk,jn) ) * e1t(ji,jj) * e2t(ji,jj) * e3t_n(ji,jj,jk)  
               END DO
            END DO
         END DO 
         !
         ! Diagnostic: Maximum value of diffusion
         !
         mxl(jn) = MAXVAL( ( zw0(:,:,:,jn) + zw5(:,:,:,jn) + zw6(:,:,:,jn) ) / r_Scale(jn) )

         !                          ! ==================
      END DO                        ! end of tracer loop
      !                             ! ==================

      call MPI_BARRIER(mpi_comm_oce, ierr)    
      call MPI_REDUCE(mxl, mxg, 2, mpi_double_precision, MPI_MAX, 0, mpi_comm_oce, ierr)
      call MPI_ALLREDUCE(swl, swg, 2, mpi_double_precision, MPI_SUM,  mpi_comm_oce, ierr)
      IF (lwp .AND. mod(kt,100) .eq. 0) THEN
         print *, '  MAXVAL: Divergence of a\nabla\theta: ', mxg(1),mxg(2), 'DIF'
         print *, '  INTEGR: Divergence of a\nabla\theta: ', swg(1),swg(2), 'DIF'
      ENDIF
      call MPI_BARRIER(mpi_comm_oce, ierr)

      DEALLOCATE( int_var11, int_var22 )

      DEALLOCATE( ztuu, ztuv, ztuw, ztvu, ztvv, ztvw )

      DEALLOCATE( zw0, zw1, zw2, zw3, zw4, zw5, zw6 )

   END SUBROUTINE tlu_trahhdiff
   ! [tlu_trahhdiff]


   !!---------------------------------------------------------------------------  
   !!            ***  ROUTINE tlu_trahzdiff ***
   !!
   !> @authors
   !>                P. Derian, P. Chandramouli, F. Tucciarone
   !!
   !> @version
   !!                2017 -  ( P. DERIAN       )  Original code <br>
   !!                2019 -  ( P. CHANDRAMOULI )  Modification of Original code <br>
   !!                2021 -  ( F. TUCCIARONE   )  Ongoing work and documentation
   !! 
   !> @brief         Computes the Horizontal components of the Stochastic Diffusion
   !!                of the Tracer 
   !!
   !> @par           Procedure specifics
   !>                Computes the horizontal terms of the stochastic diffusion, using a finite volumes approach
   !!                for the diffusive flux @f$ F_{T,z} @f$ centred on a @f$ T- @f$ point
   !!                \f{align*}{
   !!                \nabla\cdot \boldsymbol{F}_{T} &= \dfrac{1}{e_{3t}}\delta_{k}\left[ F_{T,z} \right]  
   !!                \f} 
   !!                \f{align*}{
   !!                F_{T,z} &= \overline{a_{31}^{uw}\overline{\partial_{x}T \,e_{1t}e_{2t}}^{\,k}}^{ \,i}
   !!                         + \overline{a_{32}^{vw}\overline{\partial_{y}T \,e_{1f}e_{2f}}^{ \,k}}^{\,j} 
   !!                \f}   
   !! 
   !> @par           Code specifics               
   !!                The following naming conventions are employed:
   !!                \f{align*}{
   !!                \texttt{ztwu} & = \overline{a_{31}^{uw}\overline{\partial_{x}T \,e_{1t}e_{2t}}^{\,k}}^{ \,i}\qquad 
   !!                \texttt{ztwv}   = \overline{a_{32}^{vw}\overline{\partial_{y}T \,e_{1f}e_{2f}}^{ \,k}}^{\,j}
   !!                \f} 
   !!
   !!
   !> @param[in]     kt: ocean time-step integer index
   !> @param[in]     ptb: Tracer before field (Euler scheme)
   !> @param[inout]  pta: Tracer trend to be updated
   !> @param[in]     r_Scale: Scaling factor for diffusion
   !! 
   !! @result        Update (pta) with the diffusion term
   !!
   !> @par           Diagnostic
   !!                With Homogeneous Neumann boundary conditions the integrated divergence of the diffusive fluxes
   !!                is zero.        
   !!                \f{align*}{
   !!                \int_{V}\nabla\cdot \boldsymbol{F}_{T}\,\mathrm{d}V &=
   !!                \int_{\partial V}\boldsymbol{F}\cdot\boldsymbol{n}\,\mathrm{d}s= 0
   !!                \f} 
   !!                The computational expression of this equation is        
   !!                \f{align*}{
   !!                \int_{V}\nabla\cdot \boldsymbol{F}_{T}\,\mathrm{d}V &=
   !!                \sum_{i,j,k} \nabla\cdot \boldsymbol{F}_{T} \mathrm{d} x \mathrm{d} y \mathrm{d} z = 
   !!                \sum_{i,j,k} \left\lbrace \dfrac{1}{e_{3t}}\delta_{k}\left[ F_{T,z} \right]  
   !!                  \right\rbrace e_{3t} = 0 
   !!                \f} 
   !  @note          
   !  @todo              
   !!
   !!---------------------------------------------------------------------------
   !> @snippet this tlu_trahzdiff
   ! [tlu_trahzdiff]
   SUBROUTINE tlu_trahzdiff( kt, ptb, pta, r_Scale)
      INTEGER                              , INTENT(in   ) ::   kt                    ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) ::   ptb                   ! before field (so its euler)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) ::   pta                   ! field trend
      REAL(wp), DIMENSION(jpts)            , INTENT(in   ) ::   r_Scale               ! Schmidt number
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw0                   ! accumulation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::   zw1,zw2               ! accumulation array
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var11             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var12             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var21             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   int_var22             ! interpolation array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztwu                  ! workspace array
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)              ::   ztwv                  ! workspace array
      !
      INTEGER                                              ::   ji, jj, jk            ! dummy loop indices
      INTEGER                                              ::   jn                    ! tracer index
      !
      INTEGER, PARAMETER                                   ::   ia31 = ndiffidx(3,1)  ! Mapping the indeces of a
      INTEGER, PARAMETER                                   ::   ia32 = ndiffidx(3,2)  ! Mapping the indeces of a

      REAL(wp), DIMENSION(4) :: zmaxl, zmaxg
      REAL(wp), DIMENSION(mppsize) :: ztau
      REAL(wp), DIMENSION(2)             ::   swl, swg
      REAL(wp), DIMENSION(2)             ::   mxl, mxg
      INTEGER  ::   ierr, norm_typ                   ! dummy loop argument
      !
      !!----------------------------------------------------------------------
      !
      IF( ln_timing ) CALL timing_start('tlu_trahhdiff')   ! [NEMO] check
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tlu_trahzdiff : Vertical/Horizontal diffusion of TRACER '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF

      ALLOCATE( ztwu(jpi,jpj,jpk), &
              & ztwv(jpi,jpj,jpk)  )

      ALLOCATE( zw0(jpi,jpj,jpk,2), &
              & zw1(jpi,jpj,jpk,2), &
              & zw2(jpi,jpj,jpk,2)  )
      !
      ! Inizialization of accumulation arrays
      !
      zw0 = 0._wp
      zw1 = 0._wp
      zw2 = 0._wp
      !
      ! Inizialization of local and global Diagnostic arrays
      !  
      swl = 0._wp
      swg = 0._wp
      mxl = 0._wp
      mxg = 0._wp
      !
      !
      !                             ! =========== !
      DO jn = 1, jpts               ! tracer loop !
         !                          ! =========== !   
         ztwu = 0._wp 
         ztwv = 0._wp 
         !
         ! Surface boundary condition
         ! ztuw(:,:,1) = 0._wp
         ! ztvw(:,:,1) = 0._wp
         !
         DO jk = 2, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ! Derivative of the tracer and its interpolation
                  ztwu(ji,jj,jk) = 0.5_wp * ( ( ptb(ji+1,jj  ,jk  ,jn) - ptb(ji  ,jj  ,jk  ,jn) ) /  e1u(ji,jj)  &
                                 &          + ( ptb(ji+1,jj  ,jk+1,jn) - ptb(ji  ,jj  ,jk+1,jn) ) /  e1u(ji,jj)  )
                  ztwv(ji,jj,jk) = 0.5_wp * ( ( ptb(ji  ,jj+1,jk  ,jn) - ptb(ji  ,jj  ,jk  ,jn) ) /  e2v(ji,jj)  &
                                 &          + ( ptb(ji  ,jj+1,jk+1,jn) - ptb(ji  ,jj  ,jk+1,jn) ) /  e2v(ji,jj)  )
                  ! Diffusive flux 
                  ztwu(ji,jj,jk) = e1u(ji,jj) * e2u(ji,jj) * var_ten(ji,jj,jk,ia31) * ztwu(ji,jj,jk)
                  ztwv(ji,jj,jk) = e1v(ji,jj) * e2v(ji,jj) * var_ten(ji,jj,jk,ia32) * ztwv(ji,jj,jk)
                  ! Masking
                  ztwu(ji,jj,jk) = ztwu(ji,jj,jk) * umask(ji,jj,jk) 
                  ztwv(ji,jj,jk) = ztwv(ji,jj,jk) * vmask(ji,jj,jk) 
               END DO
            END DO
         END DO
         !
         CALL lbc_lnk_multi( 'tlu_trahhdiff', ztwu, 'U', 1.,   &
                        &                     ztwv, 'V', 1.     )
         !
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !
                  ! Interpolation of extradiagonal term a_12\partial_y\theta
                  !
                  zw1(ji,jj,jk,jn) = 0.5_wp * (  ztwu(ji,jj,jk) + ztwu(ji-1,jj,jk) )   &
                  &                           * umask(ji,jj,jk)
                  !
                  ! Interpolation of extradiagonal term a_21\partial_x\theta
                  !
                  zw2(ji,jj,jk,jn) = 0.5_wp * (  ztwv(ji,jj,jk) + ztwv(ji,jj-1,jk) )   &
                  &                           * vmask(ji,jj,jk)
                  !
               END DO
            END DO
         END DO  
         !
         CALL lbc_lnk_multi( 'tlu_trahhdiff', zw1(:,:,:,jn), 'U', 1.,   &
                        &                     zw2(:,:,:,jn), 'V', 1.    )
         !
         ! Surface and Bottom boundary conditions
         !
         zw1(:,:,  1,:) = 0._wp      ! Avoid boudary fluxes
         zw1(:,:,jpk,:) = 0._wp      ! Avoid boundary fluxes
         zw2(:,:,  1,:) = 0._wp      ! Avoid boundary fluxes
         zw2(:,:,jpk,:) = 0._wp      ! Avoid boundary fluxes

         DO jk = 1, jpkm1           !==  Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  ! Vertical term of the diffusion
                  !
                  zw0(ji,jj,jk,jn) =  0.5_wp * ( zw1(ji,jj,jk,jn) - zw1(ji,jj,jk+1,jn) + &
                  &                              zw2(ji,jj,jk,jn) - zw2(ji,jj,jk+1,jn) )  &
                  &                         *  tmask(ji,jj,jk) / (  e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
                  !
                  !
                  ! Vertical term of the diffusion
                  !
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zw0(ji,jj,jk,jn) / r_Scale(jn)
                  !
                  ! Diagnostic: swl =  \sum_{i,j,k} \nabla\cdot \boldsymbol{F}_{T} \mathrm{d} x \mathrm{d} y \mathrm{d} z
                  !
                  swl(jn) = swl(jn) + zw0(ji,jj,jk,jn) * e1t(ji,jj) * e2t(ji,jj) * e3t_n(ji,jj,jk)
                  !
               END DO
            END DO
         END DO
         !
         ! Diagnostic: Maximum value of diffusion
         !
         mxl(jn) = MAXVAL( ( zw0(:,:,:,jn) ) / r_Scale(jn) )

         !                          ! ==================
      END DO                        ! end of tracer loop
      !                             ! ==================

      call MPI_BARRIER(mpi_comm_oce, ierr)    
      call MPI_REDUCE(mxl, mxg, 2, mpi_double_precision, MPI_MAX, 0, mpi_comm_oce, ierr)
      call MPI_ALLREDUCE(swl, swg, 2, mpi_double_precision, MPI_SUM,  mpi_comm_oce, ierr)
      IF (lwp .AND. mod(kt,100) .eq. 0) THEN
         print *, '  MAXVAL: Divergence of a\nabla\theta: ', mxg(1),mxg(2), 'DIF'
         print *, '  INTEGR: Divergence of a\nabla\theta: ', swg(1),swg(2), 'DIF'
      ENDIF
      call MPI_BARRIER(mpi_comm_oce, ierr)

      DEALLOCATE( ztwu, ztwv )

      DEALLOCATE( zw0, zw1, zw2 )

   END SUBROUTINE tlu_trahzdiff
   ! [tlu_trahzdiff]

END MODULE tluhdifTRA




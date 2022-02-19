MODULE dynldf_lap
   !!======================================================================
   !!                       ***  MODULE  dynldf_lap  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  ! 1990-09 (G. Madec) Original code
   !!            4.0  ! 1991-11 (G. Madec)
   !!            6.0  ! 1996-01 (G. Madec) statement function for e3 and ahm
   !!   NEMO     1.0  ! 2002-06 (G. Madec)  F90: Free form and module
   !!             -   ! 2004-08 (C. Talandier) New trends organization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap  : update the momentum trend with the lateral diffusion
   !!                  using an iso-level harmonic operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE ldftra_oce, ONLY : lk_ldfeke   ! GEOMETRIC param. activation
   USE zdf_oce         ! ocean vertical physics
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !
   USE in_out_manager  ! I/O manager
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by step.F90
   
!! jm (28 Oct 18): to tidy up, definitions should go in some other file   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   eke_keS   !: Source term of EKE equation used in ldfeke module

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is an
      !!      harmonic operator (laplacian type) which separates the divergent
      !!      and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!         difu = 1/e1u di[ahmt hdivb] - 1/(e2u*e3u) dj-1[e3f ahmf rotb]
      !!         difv = 1/e2v dj[ahmt hdivb] + 1/(e1v*e3v) di-1[e3f ahmf rotb]
      !!      in the rotational part of the diffusion.
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!
      !! ** Action : - Update (ua,va) with the iso-level harmonic mixing trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   zua, zva, ze2u, ze1v   ! local scalars
!!
!! jm (28 Oct 18): to tidy up, definitions should go in some other file   
      REAL(wp), DIMENSION(jpi,jpj) ::  zah_div2, zah_cur2
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_lap')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze2u = rotb (ji,jj,jk) * fsahmf(ji,jj,jk) * fse3f(ji,jj,jk) ! jm: this is "zcur"
               ze1v = hdivb(ji,jj,jk) * fsahmt(ji,jj,jk)                   ! jm: this is "zdiv"
               ! horizontal diffusive trends
               zua = - ( ze2u - rotb (ji,jj-1,jk)*fsahmf(ji,jj-1,jk)*fse3f(ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     + ( hdivb(ji+1,jj,jk)*fsahmt(ji+1,jj,jk) - ze1v                   ) / e1u(ji,jj)

               zva = + ( ze2u - rotb (ji-1,jj,jk)*fsahmf(ji-1,jj,jk)*fse3f(ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     + ( hdivb(ji,jj+1,jk)*fsahmt(ji,jj+1,jk) - ze1v                   ) / e2v(ji,jj)

               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
               !
               IF( lk_ldfeke ) THEN        ! GEOMETRIC source term
                  zah_cur2(ji,jj) = zah_cur2(ji,jj) +                     ze2u**2 / MAX( 1._wp , fsahmf(ji,jj,jk) ) * fmask(ji,jj,jk)
                  zah_div2(ji,jj) = zah_div2(ji,jj) + fse3t_b(ji,jj,jk) * ze1v**2 / MAX( 1._wp , fsahmt(ji,jj,jk) ) * tmask(ji,jj,jk)
               ENDIF
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      
      IF( lk_ldfeke ) THEN        ! GEOMETRIC source term        
         CALL lbc_lnk( zah_cur2, 'F', 1. )
         zah_cur2(:,:) = zah_cur2(:,:) * e1f(:,:) * e2f(:,:)
         DO jj = 2, jpjm1
            DO ji = fs_2, jpim1   ! vector opt.
               eke_keS(ji,jj) = zah_div2(ji,jj) + (  zah_cur2(ji-1,jj  )   + zah_cur2(ji,jj  )        &
                  &                                + zah_cur2(ji-1,jj-1)   + zah_cur2(ji,jj-1)    )   &
                  &                 / MAX(  1._wp ,  fmask   (ji-1,jj  ,1) + fmask   (ji,jj  ,1)      &
                  &                                + fmask   (ji-1,jj-1,1) + fmask   (ji,jj-1,1)  )   * r1_e12t(ji,jj)
            END DO  
         END DO  
         CALL lbc_lnk( eke_keS, 'T', 1. )
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_lap')
      !
   END SUBROUTINE dyn_ldf_lap

   !!======================================================================
END MODULE dynldf_lap

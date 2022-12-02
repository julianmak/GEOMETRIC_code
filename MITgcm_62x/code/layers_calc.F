C $Header: /u/gcmpack/MITgcm/pkg/layers/layers_calc.F,v 1.8 2011/05/12 15:09:54 dfer Exp $
C $Name: checkpoint62x $

#include "LAYERS_OPTIONS.h"
#ifdef ALLOW_GMREDI
#include "GMREDI_OPTIONS.h"
#endif

CBOP 0
C !ROUTINE: LAYERS_CALC

C !INTERFACE:
      SUBROUTINE LAYERS_CALC(
     I                  myTime, myIter, myThid )

C !DESCRIPTION:
C ===================================================================
C     Calculate the transport in isopycnal layers.
C     This is the meat of the LAYERS package.
C ===================================================================

C !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "LAYERS_SIZE.h"
#include "LAYERS.h"
#ifdef ALLOW_GMREDI
# include "GMREDI.h"
#endif

C !INPUT PARAMETERS:
C     myTime :: Current time in simulation
C     myIter :: Current iteration number
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_LAYERS

C !LOCAL VARIABLES:
C     bi, bj   :: tile indices
C     i,j      :: horizontal indices
C     k        :: vertical index for model grid
C     kci      :: index from CellIndex
C     kg       :: index for looping though layers_G
C     kk       :: vertical index for ZZ (fine) grid
C     kgu,kgv  :: vertical index for isopycnal grid
C     level    :: counter for vertical level in prho calculation
C     prho     :: potential density referenced to layers_kref pressure
C     TatV     :: temperature at U point
C     TatV     :: temperature at V point

      INTEGER bi, bj
      INTEGER i,j,k,kk,kg,kci
      INTEGER level
      INTEGER kgu(sNx+1,sNy+1), kgv(sNx+1,sNy+1)
      _RL TatU, TatV
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#if (defined ALLOW_GMREDI) && (defined GM_BOLUS_ADVEC)
      INTEGER kcip1
      _RL delPsi, maskp1
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C --- The tile loops
      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)

C     Initialize the search indices
      DO j = 1,sNy+1
        DO i = 1,sNx+1
C       The temperature index (layer_G) goes from cold to warm.
C       The water column goes from warm (k=1) to cold (k=Nr).
C       So initialize the search with the warmest value.
          kgu(i,j) = Nlayers
          kgv(i,j) = Nlayers
        ENDDO
      ENDDO

C     Reset the arrays
      DO kg=1,Nlayers
       DO j = 1,sNy+1
        DO i = 1,sNx+1
#ifdef LAYERS_UFLUX
         layers_UFlux(i,j,kg,bi,bj) = 0. _d 0
#ifdef ALLOW_GMREDI
         layers_GMU(i,j,kg,bi,bj) = 0. _d 0
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         layers_HU(i,j,kg,bi,bj) = 0. _d 0
         layers_UU(i,j,kg,bi,bj) = 0. _d 0
         layers_ZU(i,j,kg,bi,bj) = 0. _d 0
         layers_ZZU(i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
         layers_VFlux(i,j,kg,bi,bj) = 0. _d 0
#ifdef ALLOW_GMREDI
         layers_GMV(i,j,kg,bi,bj) = 0. _d 0
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         layers_HV(i,j,kg,bi,bj) = 0. _d 0
         layers_VV(i,j,kg,bi,bj) = 0. _d 0
         layers_ZV(i,j,kg,bi,bj) = 0. _d 0
         layers_ZZV(i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */
        ENDDO
       ENDDO
      ENDDO

#ifdef LAYERS_PRHO_REF
C     For LAYER_nb = 3, calculate the potential density referenced to
C     the model level given by layers_kref.
      IF (LAYER_nb .EQ. 3) THEN
       DO level = 1,Nr
        CALL FIND_RHO_2D( 1-OLx, sNx+OLx, 1-OLy, sNy+OLy,
     &     layers_kref,
     &     theta(1-OLx,1-OLy,level,bi,bj),
     &     salt(1-OLx,1-OLy,level,bi,bj),
     &     prho(1-OLx,1-OLy,level,bi,bj),
     &     level, bi, bj, myThid )
        prho(1-OLx:sNx+OLy,1-OLy:sNy+OLy,level,bi,bj) = rhoConst
     &   + prho(1-OLx:sNx+OLy,1-OLy:sNy+OLy,level,bi,bj)
       ENDDO
      ELSE
      ENDIF
#endif

C      _RL  theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C     Sometimes it is done this way
C      DO j=1-Oly+1,sNy+Oly-1
C       DO i=1-Olx+1,sNx+Olx-1
      DO kk=1,NZZ
       k = MapIndex(kk)
       kci = CellIndex(kk)
       DO j = 1,sNy+1
        DO i = 1,sNx+1

#ifdef LAYERS_UFLUX
C ------ Find theta at the U point (west) on the fine Z grid
         IF (LAYER_nb .EQ. 1) THEN
C        Use theta to determine layers.
	  IF (hFacW(i,j,k+1,bi,bj) .EQ. 0) THEN
	   TatU = 0.5 _d 0 * (theta(i-1,j,k,bi,bj)+theta(i,j,k,bi,bj))
	  ELSE
           TatU = MapFact(kk) *
     &      0.5 _d 0 * (theta(i-1,j,k,bi,bj)+theta(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (theta(i-1,j,k+1,bi,bj)+theta(i,j,k+1,bi,bj))
	  ENDIF
         ELSEIF (LAYER_nb .EQ. 2) THEN
C        Use salinity to determine layers.
	  IF (hFacW(i,j,k+1,bi,bj) .EQ. 0) THEN
	   TatU = 0.5 _d 0 * (salt(i-1,j,k,bi,bj)+salt(i,j,k,bi,bj))
	  ELSE
           TatU = MapFact(kk) *
     &      0.5 _d 0 * (salt(i-1,j,k,bi,bj)+salt(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (salt(i-1,j,k+1,bi,bj)+salt(i,j,k+1,bi,bj))
	  ENDIF
#ifdef LAYERS_PRHO_REF
         ELSEIF (LAYER_nb .EQ. 3)THEN
C        Use potential density to determine layers.
	  IF (hFacW(i,j,k+1,bi,bj) .EQ. 0) THEN
	   TatU = 0.5 _d 0 * (prho(i-1,j,k,bi,bj)+prho(i,j,k,bi,bj))
	  ELSE
           TatU = MapFact(kk) *
     &      0.5 _d 0 * (prho(i-1,j,k,bi,bj)+prho(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (prho(i-1,j,k+1,bi,bj)+prho(i,j,k+1,bi,bj))
	  ENDIF
#endif
         ENDIF

C ------ Now that we know T everywhere, determine the binning.

         IF (TatU .GE. layers_G(Nlayers)) THEN
C        the point is in the hottest bin or hotter
          kgu(i,j) = Nlayers
         ELSE IF (TatU .LT. layers_G(2)) THEN
C        the point is in the coldest bin or colder
          kgu(i,j) = 1
         ELSE IF ( (TatU .GE. layers_G(kgu(i,j)))
     &    .AND. (TatU .LT. layers_G(kgu(i,j)+1)) ) THEN
C        already on the right bin -- do nothing
         ELSE IF (TatU .GE. layers_G(kgu(i,j))) THEN
C        have to hunt for the right bin by getting hotter
          DO WHILE (TatU .GE. layers_G(kgu(i,j)+1))
           kgu(i,j) = kgu(i,j) + 1
          ENDDO
C         now layers_G(kgu(i,j)+1) < TatU <= layers_G(kgu(i,j)+1)
         ELSE IF (TatU .LT. layers_G(kgu(i,j)+1)) THEN
C        have to hunt for the right bin by getting colder
          DO WHILE (TatU .LT. layers_G(kgu(i,j)))
           kgu(i,j) = kgu(i,j) - 1
          ENDDO
C         now layers_G(kgu(i,j)+1) <= TatU < layers_G(kgu(i,j)+1)
         ELSE
C        that should have covered all the options
          WRITE(msgBuf,'(A,1E14.6)')
     &     'S/R LAYERS_CALC: Couldnt find a bin in layers_G for TatU=',
     &     TatU
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R LAYERS_INIT_FIXED'
         END IF

C ------ Augment the bin values
         layers_UFlux(i,j,kgu(i,j),bi,bj) =
     &    layers_UFlux(i,j,kgu(i,j),bi,bj) +
     &    dZZf(kk) * uVel(i,j,kci,bi,bj) * hFacW(i,j,kci,bi,bj)

#ifdef ALLOW_GMREDI
         IF ( GM_AdvForm .AND. useBOLUS ) THEN
           kcip1 = MIN(kci+1,Nr)
           maskp1 = 1.
           IF (kci.GE.Nr) maskp1 = 0.
           delPsi = GM_PsiX(i,j,kcip1,bi,bj)*maskp1
     &            - GM_PsiX(i,j, kci, bi,bj)
           layers_GMU(i,j,kgu(i,j),bi,bj) =
     &      layers_GMU(i,j,kgu(i,j),bi,bj)
     &      + delPsi*recip_drF(kci)*_recip_hFacW(i,j,kci,bi,bj)
     &      * dZZf(kk)*hFacW(i,j,kci,bi,bj)
         ELSEIF ( useBOLUS .AND. .NOT. GM_AdvForm ) THEN
           kcip1 = MIN(kci+1,Nr)
           maskp1 = 1.
           IF (kci.GE.Nr) maskp1 = 0.
             delPsi = 0.5 _d 0 *(
     &         0.5 _d 0 * (rA(i-1,j,bi,bj)*Kwx(i-1,j,kcip1,bi,bj)
     &                    +rA( i ,j,bi,bj)*Kwx( i ,j,kcip1,bi,bj))
     &                  * maskW(i,j,kcip1,bi,bj) * maskp1
     &        -0.5 _d 0 * (rA(i-1,j,bi,bj)*Kwx(i-1,j, kci ,bi,bj)
     &                    +rA( i ,j,bi,bj)*Kwx( i ,j, kci ,bi,bj))
     &                  * maskW(i,j, kci ,bi,bj)
     &                               ) / rAw(i,j,bi,bj)
           layers_GMU(i,j,kgu(i,j),bi,bj) =
     &      layers_GMU(i,j,kgu(i,j),bi,bj)
     &      + delPsi*recip_drF(kci)*_recip_hFacW(i,j,kci,bi,bj)
     &      * dZZf(kk)*hFacW(i,j,kci,bi,bj)
         ENDIF
#endif

#ifdef LAYERS_THICKNESS
         layers_HU(i,j,kgu(i,j),bi,bj) = layers_HU(i,j,kgu(i,j),bi,bj)
     &    + dZZf(kk) * hFacW(i,j,kci,bi,bj)
#endif /* LAYERS_THICKNESS */

#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
C ------ Find theta at the V point (south) on the fine Z grid
         IF (LAYER_nb .EQ. 1) THEN
C        Use theta to determine layers.
	  IF (hFacS(i,j,k+1,bi,bj) .EQ. 0) THEN
           TatV = 0.5 _d 0 * (theta(i,j-1,k,bi,bj)+theta(i,j,k,bi,bj))
	  ELSE
           TatV = MapFact(kk) *
     &      0.5 _d 0 * (theta(i,j-1,k,bi,bj)+theta(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (theta(i,j-1,k+1,bi,bj)+theta(i,j,k+1,bi,bj))
          ENDIF
         ELSEIF (LAYER_nb .EQ. 2) THEN
C        Use salinity to determine layers.
	  IF (hFacS(i,j,k+1,bi,bj) .EQ. 0) THEN
           TatV = 0.5 _d 0 * (salt(i,j-1,k,bi,bj)+salt(i,j,k,bi,bj))
	  ELSE
           TatV = MapFact(kk) *
     &      0.5 _d 0 * (salt(i,j-1,k,bi,bj)+salt(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (salt(i,j-1,k+1,bi,bj)+salt(i,j,k+1,bi,bj))
          ENDIF
#ifdef LAYERS_PRHO_REF
         ELSEIF (LAYER_nb .EQ. 3)THEN
C        Use potential density to determine layers.
	  IF (hFacS(i,j,k+1,bi,bj) .EQ. 0) THEN
           TatV = 0.5 _d 0 * (prho(i,j-1,k,bi,bj)+prho(i,j,k,bi,bj))
	  ELSE
           TatV = MapFact(kk) *
     &      0.5 _d 0 * (prho(i,j-1,k,bi,bj)+prho(i,j,k,bi,bj)) +
     &      (1-MapFact(kk)) *
     &      0.5 _d 0 * (prho(i,j-1,k+1,bi,bj)+prho(i,j,k+1,bi,bj))
          ENDIF
#endif
         ENDIF

C ------ Now that we know T everywhere, determine the binning
         IF (TatV .GE. layers_G(Nlayers)) THEN
C         the point is in the hottest bin or hotter
          kgv(i,j) = Nlayers
         ELSE IF (TatV .LT. layers_G(2)) THEN
C         the point is in the coldest bin or colder
          kgv(i,j) = 1
         ELSE IF ( (TatV .GE. layers_G(kgv(i,j)))
     &    .AND. (TatV .LT. layers_G(kgv(i,j)+1)) ) THEN
C         already on the right bin -- do nothing
         ELSE IF (TatV .GE. layers_G(kgv(i,j))) THEN
C         have to hunt for the right bin by getting hotter
          DO WHILE (TatV .GE. layers_G(kgv(i,j)+1))
           kgv(i,j) = kgv(i,j) + 1
          ENDDO
C         now layers_G(kgv(i,j)+1) < TatV <= layers_G(kgv(i,j)+1)
         ELSE IF (TatV .LT. layers_G(kgv(i,j)+1)) THEN
C         have to hunt for the right bin by getting colder
          DO WHILE (TatV .LT. layers_G(kgv(i,j)))
           kgv(i,j) = kgv(i,j) - 1
          ENDDO
C         now layers_G(kgv(i,j)+1) <= TatV < layers_G(kgv(i,j)+1)
         ELSE
C         that should have covered all the options
          WRITE(msgBuf,'(A,1E14.6)')
     &     'S/R LAYERS_CALC: Couldnt find a bin in layers_G for TatV=',
     &     TatV
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R LAYERS_INIT_FIXED'
         END IF

C ------ Augment the bin values
         layers_VFlux(i,j,kgv(i,j),bi,bj) =
     &    layers_VFlux(i,j,kgv(i,j),bi,bj)
     &    + dZZf(kk) * vVel(i,j,kci,bi,bj) * hFacS(i,j,kci,bi,bj)

#ifdef ALLOW_GMREDI
         IF ( GM_AdvForm .AND. useBOLUS ) THEN
           kcip1 = MIN(kci+1,Nr)
           maskp1 = 1.
           IF (kci.GE.Nr) maskp1 = 0.
           delPsi = GM_PsiY(i,j,kcip1,bi,bj)*maskp1
     &            - GM_PsiY(i,j, kci, bi,bj)
           layers_GMV(i,j,kgv(i,j),bi,bj) =
     &      layers_GMV(i,j,kgv(i,j),bi,bj)
     &      + delPsi*recip_drF(kci)*_recip_hFacS(i,j,kci,bi,bj)
     &      * dZZf(kk)*hFacS(i,j,kci,bi,bj)
         ELSEIF ( useBOLUS .AND. .NOT. GM_AdvForm ) THEN
           kcip1 = MIN(kci+1,Nr)
           maskp1 = 1.
           IF (kci.GE.Nr) maskp1 = 0.
             delPsi = 0.5 _d 0 *(
     &         0.5 _d 0 * (rA(i,j-1,bi,bj)*Kwy(i,j-1,kcip1,bi,bj)
     &                    +rA(i, j ,bi,bj)*Kwy(i, j ,kcip1,bi,bj))
     &                  * maskS(i,j,kcip1,bi,bj) * maskp1
     &        -0.5 _d 0 * (rA(i,j-1,bi,bj)*Kwy(i,j-1, kci ,bi,bj)
     &                    +rA(i, j ,bi,bj)*Kwy(i, j , kci ,bi,bj))
     &                  * maskS(i,j, kci ,bi,bj)
     &                               ) / rAs(i,j,bi,bj)
           layers_GMV(i,j,kgv(i,j),bi,bj) =
     &      layers_GMV(i,j,kgv(i,j),bi,bj)
     &      + delPsi*recip_drF(kci)*_recip_hFacS(i,j,kci,bi,bj)
     &      * dZZf(kk)*hFacS(i,j,kci,bi,bj)
         ENDIF
#endif

#ifdef LAYERS_THICKNESS
         layers_HV(i,j,kgv(i,j),bi,bj) = layers_HV(i,j,kgv(i,j),bi,bj)
     &    + dZZf(kk) * hFacS(i,j,kci,bi,bj)
#endif /* LAYERS_THICKNESS */

#endif /* LAYERS_VFLUX */

        ENDDO
       ENDDO
      ENDDO

#ifdef LAYERS_THICKNESS
      DO kg = 1,Nlayers
       DO j = 1,sNy+1
        DO i = 1,sNx+1
#ifdef LAYERS_UFLUX
C--   Calculate the U velocity within the layer.
         IF ( layers_HU(i,j,kg,bi,bj) .NE. 0.0 ) THEN
          layers_UU(i,j,kg,bi,bj) = layers_UFlux(i,j,kg,bi,bj) 
     &                              / layers_HU(i,j,kg,bi,bj)
         ELSE
          layers_UU(i,j,kg,bi,bj) = 0.0
         ENDIF
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
C--   Calculate the V velocity within the layer.
         IF ( layers_HV(i,j,kg,bi,bj) .NE. 0.0 ) THEN
          layers_VV(i,j,kg,bi,bj) = layers_VFlux(i,j,kg,bi,bj) 
     &                              / layers_HV(i,j,kg,bi,bj)
         ELSE
          layers_VV(i,j,kg,bi,bj) = 0.0
         ENDIF
#endif /* LAYERS_VFLUX */
        ENDDO
       ENDDO
      ENDDO 
#endif /* LAYERS_THICKNESS */

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Julian's dirty hacks
C
C--   Compute Z(rho) and Z^2(rho) on the fly
C--   Depending on use of potential density or temperature,
C--     cumulative sum may need to be done in the opposite directions
#ifdef LAYERS_THICKNESS
      do j = 1,sNy+1
       do i = 1,sNx+1
#ifdef LAYERS_UFLUX
C--   Using potential density as a gridding field, starting from
C--     the *** most *** dense counter
        if (LAYER_nb .eq. 3) then
         layers_ZU(i,j,Nlayers,bi,bj) = layers_HU(i,j,Nlayers,bi,bj)
         do kg = Nlayers-1,1,-1
          layers_ZU(i,j,kg,bi,bj) = layers_HU(i,j,kg,bi,bj) 
     &                            + layers_ZU(i,j,kg+1,bi,bj)
         end do
C--     ...take off the *** lightest *** layer to get depth
C       (!!! assume lightest layer really is at the top though!)
         do kg = Nlayers,1,-1
          layers_ZU(i,j,kg,bi,bj) = layers_ZU(i,j,kg,bi,bj)
     &                            - layers_ZU(i,j,1,bi,bj)
          layers_ZZU(i,j,kg,bi,bj) = layers_ZU(i,j,kg,bi,bj)
     &                             * layers_ZU(i,j,kg,bi,bj)
         end do
C--   Using temperature as a gridding field, starting from coldest
C--     temperature (most dense) so sum as normal...
        else
         layers_ZU(i,j,1,bi,bj) = layers_HU(i,j,1,bi,bj)
         do kg = 2,Nlayers
          layers_ZU(i,j,kg,bi,bj) = layers_HU(i,j,kg,bi,bj)
     &                            + layers_ZU(i,j,kg-1,bi,bj)
         end do
C--     ...take off the *** lightest *** layer to get depth
C       (!!! assume lightest layer really is at the top though!)
         do kg = 1,Nlayers
          layers_ZU(i,j,kg,bi,bj) = layers_ZU(i,j,kg,bi,bj)
     &                            - layers_ZU(i,j,Nlayers,bi,bj)
          layers_ZZU(i,j,kg,bi,bj) = layers_ZU(i,j,kg,bi,bj)
     &                             * layers_ZU(i,j,kg,bi,bj)
         end do
        end if
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
C--   Calculate the V velocity within the layer.
        if (LAYER_nb .eq. 3) then
         layers_ZV(i,j,Nlayers,bi,bj) = layers_HV(i,j,Nlayers,bi,bj)
         do kg = Nlayers-1,1,-1
          layers_ZV(i,j,kg,bi,bj) = layers_HV(i,j,kg,bi,bj) 
     &                            + layers_ZV(i,j,kg+1,bi,bj)
         end do
C--     ...take off the *** lightest *** layer to get depth
C       (!!! assume lightest layer really is at the top though!)
         do kg = Nlayers,1,-1
          layers_ZV(i,j,kg,bi,bj) = layers_ZV(i,j,kg,bi,bj)
     &                            - layers_ZV(i,j,1,bi,bj)
          layers_ZZV(i,j,kg,bi,bj) = layers_ZV(i,j,kg,bi,bj)
     &                             * layers_ZV(i,j,kg,bi,bj)
         end do
C--   Using temperature as a gridding field, starting from coldest
C--     temperature (most dense) so sum as normal...
        else
         layers_ZV(i,j,1,bi,bj) = layers_HV(i,j,1,bi,bj)
         do kg = 2,Nlayers
          layers_ZV(i,j,kg,bi,bj) = layers_HV(i,j,kg,bi,bj)
     &                            + layers_ZV(i,j,kg-1,bi,bj)
         end do
C--     ...take off the *** lightest *** layer to get depth
C       (!!! assume lightest layer really is at the top though!)
         do kg = 1,Nlayers
          layers_ZV(i,j,kg,bi,bj) = layers_ZV(i,j,kg,bi,bj)
     &                            - layers_ZV(i,j,Nlayers,bi,bj)
          layers_ZZV(i,j,kg,bi,bj) = layers_ZV(i,j,kg,bi,bj)
     &                             * layers_ZV(i,j,kg,bi,bj)
         end do
        end if
#endif /* LAYERS_VFLUX */
       end do
      end do
#endif /* LAYERS_THICKNESS */

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef ALLOW_TIMEAVE
C--   Time-average
      IF ( layers_taveFreq.GT.0. ) THEN

#ifdef LAYERS_UFLUX
         CALL TIMEAVE_CUMULATE( layers_UFlux_T, layers_UFlux, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#ifdef ALLOW_GMREDI
         CALL TIMEAVE_CUMULATE( layers_GMU_T, layers_GMU, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         CALL TIMEAVE_CUMULATE( layers_HU_T, layers_HU, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_UU_T, layers_UU, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_ZU_T, layers_ZU, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_ZZU_T, layers_ZZU, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
         CALL TIMEAVE_CUMULATE( layers_VFlux_T, layers_VFlux, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#ifdef ALLOW_GMREDI
         CALL TIMEAVE_CUMULATE( layers_GMV_T, layers_GMV, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         CALL TIMEAVE_CUMULATE( layers_HV_T, layers_HV, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_VV_T, layers_VV, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_ZV_T, layers_ZV, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( layers_ZZV_T, layers_ZZV, Nlayers,
     &                          deltaTclock, bi, bj, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
         CALL TIMEAVE_CUMULATE( prho_tave, prho, Nr,
     &                          deltaTclock, bi, bj, myThid )
#endif /* LAYERS_PRHO_REF */

         layers_TimeAve(bi,bj)=layers_TimeAve(bi,bj)+deltaTclock

      ENDIF
#endif /* ALLOW_TIMEAVE */

C --- End bi,bj loop
      ENDDO
      ENDDO

#endif /* ALLOW_LAYERS */

      RETURN
      END

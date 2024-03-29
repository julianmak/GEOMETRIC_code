C $Header: /u/gcmpack/MITgcm/pkg/layers/layers_output.F,v 1.4 2009/12/28 02:37:52 jmc Exp $
C $Name: checkpoint62x $

#include "LAYERS_OPTIONS.h"

CBOP 0
C     !ROUTINE: LAYERS_OUTPUT

C     !INTERFACE:
      SUBROUTINE LAYERS_OUTPUT( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE LAYERS_OUTPUT
C     | o general routine for LAYERS output
C     *==========================================================*
C     |   write time-average & snap-shot output
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "LAYERS_SIZE.h"
#include "LAYERS.h"

C     !INPUT PARAMETERS:
C     == Routine arguments ==
C     myTime :: Current time of simulation ( s )
C     myIter :: Iteration number
C     myThid :: my Thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_LAYERS

C     !LOCAL VARIABLES:
C     == Local variables ==
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      CHARACTER*(MAX_LEN_MBUF) suff
#ifdef ALLOW_TIMEAVE
      INTEGER bi, bj
#endif
#ifdef ALLOW_MNC
      CHARACTER*(1) pf
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_MNC
      IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
        pf(1:1) = 'D'
      ELSE
        pf(1:1) = 'R'
      ENDIF
#endif /* ALLOW_MNC */

      IF ( DIFFERENT_MULTIPLE(layers_diagFreq,myTime,deltaTClock)
     &  .AND. myIter .GT. 0 ) THEN

        IF ( layers_MDSIO ) THEN
          WRITE(suff,'(I10.10)') myIter
#ifdef LAYERS_UFLUX
            CALL WRITE_FLD_3D_RL( 'layers_UFlux.', suff, Nlayers,
     &                             layers_UFlux, myIter, myThid )
#ifdef ALLOW_GMREDI
            CALL WRITE_FLD_3D_RL( 'layers_GMU.', suff, Nlayers,
     &                             layers_GMU, myIter, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
            CALL WRITE_FLD_3D_RL( 'layers_HU.', suff, Nlayers,
     &                             layers_HU, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_UU.', suff, Nlayers,
     &                             layers_UU, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_ZU.', suff, Nlayers,
     &                             layers_ZU, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_ZZU.', suff, Nlayers,
     &                             layers_ZZU, myIter, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
            CALL WRITE_FLD_3D_RL( 'layers_VFlux.', suff, Nlayers,
     &                             layers_VFlux, myIter, myThid )
#ifdef ALLOW_GMREDI
            CALL WRITE_FLD_3D_RL( 'layers_GMV.', suff, Nlayers,
     &                             layers_GMV, myIter, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
            CALL WRITE_FLD_3D_RL( 'layers_HV.', suff, Nlayers,
     &                             layers_HV, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_VV.', suff, Nlayers,
     &                             layers_VV, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_ZV.', suff, Nlayers,
     &                             layers_ZV, myIter, myThid )
            CALL WRITE_FLD_3D_RL( 'layers_ZZV.', suff, Nlayers,
     &                             layers_ZZV, myIter, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */
#ifdef LAYERS_PRHO_REF
            CALL WRITE_FLD_3D_RL( 'layers_prho.', suff, Nr,
     &                             prho, myIter, myThid )
#endif /* LAYERS_PRHO_REF */
          ENDIF
        ENDIF

#ifdef ALLOW_MNC
#ifdef LAYERS_MNC
        IF ( layers_MNC) THEN
C           Do MNC output...
C           But how?
        ENDIF
#endif /* LAYERS_MNC */
#endif /* ALLOW_MNC */

#ifdef ALLOW_TIMEAVE
      IF ( layers_taveFreq.GT.0. ) THEN

c set arrays to zero if first timestep
       IF ( myIter.EQ.nIter0 ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
#ifdef LAYERS_UFLUX
          CALL TIMEAVE_RESET(layers_UFlux_T,Nlayers,bi,bj,myThid)
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_RESET(layers_GMU_T,Nlayers,bi,bj,myThid)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_RESET(layers_HU_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_UU_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_ZU_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_ZZU_T,Nlayers,bi,bj,myThid)
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
          CALL TIMEAVE_RESET(layers_VFlux_T,Nlayers,bi,bj,myThid)
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_RESET(layers_GMV_T,Nlayers,bi,bj,myThid)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_RESET(layers_HV_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_VV_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_ZV_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_ZZV_T,Nlayers,bi,bj,myThid)
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */
#ifdef LAYERS_PRHO_REF
          CALL TIMEAVE_RESET(prho_tave,Nr,bi,bj,mythid)
#endif /* LAYERS_PRHO_REF */
          layers_TimeAve(bi,bj) = 0.
         ENDDO
        ENDDO

C     Dump files and restart average computation if needed
       ELSEIF (
     &  DIFFERENT_MULTIPLE(layers_taveFreq,myTime,deltaTClock)
     &        ) THEN

C      Normalize by integrated time
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)

#ifdef LAYERS_UFLUX
          CALL TIMEAVE_NORMALIZE( layers_UFlux_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_NORMALIZE( layers_GMU_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_NORMALIZE( layers_HU_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_UU_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_ZU_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_ZZU_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
          CALL TIMEAVE_NORMALIZE( layers_VFlux_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_NORMALIZE( layers_GMV_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_NORMALIZE( layers_HV_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_VV_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_ZV_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
          CALL TIMEAVE_NORMALIZE( layers_ZZV_T, layers_timeave,
     &                            Nlayers, bi, bj, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
          CALL TIMEAVE_NORMALIZE( prho_tave, layers_timeave,
     &                            Nr, bi, bj, myThid )
#endif /* LAYERS_PRHO_REF */

         ENDDO
        ENDDO

        IF ( layers_MDSIO ) THEN
         WRITE(suff,'(I10.10)') myIter
#ifdef LAYERS_UFLUX
         CALL WRITE_FLD_3D_RL( 'layers_UFlux-tave.', suff, Nlayers,
     &                          layers_UFlux_T, myIter, myThid )
#ifdef ALLOW_GMREDI
         CALL WRITE_FLD_3D_RL( 'layers_GMU-tave.', suff, Nlayers,
     &                          layers_GMU_T, myIter, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         CALL WRITE_FLD_3D_RL( 'layers_HU-tave.', suff, Nlayers,
     &                          layers_HU_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_UU-tave.', suff, Nlayers,
     &                          layers_HU_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_ZU-tave.', suff, Nlayers,
     &                          layers_ZU_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_ZZU-tave.', suff, Nlayers,
     &                          layers_ZZU_T, myIter, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */
#ifdef LAYERS_VFLUX
         CALL WRITE_FLD_3D_RL( 'layers_VFlux-tave.', suff, Nlayers,
     &                          layers_VFlux_T, myIter, myThid )
#ifdef ALLOW_GMREDI
         CALL WRITE_FLD_3D_RL( 'layers_GMV-tave.', suff, Nlayers,
     &                          layers_GMV_T, myIter, myThid )
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
         CALL WRITE_FLD_3D_RL( 'layers_HV-tave.', suff, Nlayers,
     &                          layers_HV_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_VV-tave.', suff, Nlayers,
     &                          layers_VV_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_ZV-tave.', suff, Nlayers,
     &                          layers_ZV_T, myIter, myThid )
         CALL WRITE_FLD_3D_RL( 'layers_ZZV-tave.', suff, Nlayers,
     &                          layers_ZZV_T, myIter, myThid )
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
         CALL WRITE_FLD_3D_RL( 'layers_prho-tave.', suff, Nr,
     &                          prho_tave, myIter, myThid )
#endif /* LAYERS_PRHO_REF */

        ENDIF

#ifdef ALLOW_MNC
C     Do MNC output.
#endif

C      Reset averages to zero
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
#ifdef LAYERS_UFLUX
          CALL TIMEAVE_RESET(layers_UFlux_T,Nlayers,bi,bj,myThid)
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_RESET(layers_HU_T,Nlayers,bi,bj,myThid)
#endif /* LAYERS_THICKNESS */
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_RESET(layers_GMU_T,Nlayers,bi,bj,myThid)
#endif /* ALLOW_GMREDI */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
          CALL TIMEAVE_RESET(layers_VFlux_T,Nlayers,bi,bj,myThid)
#ifdef ALLOW_GMREDI
          CALL TIMEAVE_RESET(layers_GMV_T,Nlayers,bi,bj,myThid)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
          CALL TIMEAVE_RESET(layers_HV_T,Nlayers,bi,bj,myThid)
          CALL TIMEAVE_RESET(layers_VV_T,Nlayers,bi,bj,myThid)
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
          CALL TIMEAVE_RESET(prho_tave,Nr,bi,bj,myThid)
#endif /* LAYERS_PRHO_REF */

          layers_TimeAve(bi,bj) = 0.
         ENDDO
        ENDDO

C--   end of bloc: if time is a multiple of layers_taveFreq
       ENDIF

      ENDIF
#endif /* ALLOW_TIMEAVE */

#endif /* ALLOW_LAYERS */

      RETURN
      END

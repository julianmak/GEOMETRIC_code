C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_check.F,v 1.24 2011/02/10 21:24:19 jmc Exp $
C $Name: checkpoint62x $

#include "GMREDI_OPTIONS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: GMREDI_CHECK
C     !INTERFACE:
      SUBROUTINE GMREDI_CHECK( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_CHECK
C     | o Check consistency with model configuration
C     *==========================================================*
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GMREDI.h"
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid :: my Thread Id number
      INTEGER myThid

#ifdef ALLOW_GMREDI
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP
#ifdef ALLOW_PTRACERS
      INTEGER iTr
      LOGICAL redFlag
#endif

      _BEGIN_MASTER(myThid)

       WRITE(msgBuf,'(A)') 'GMREDI_CHECK: #define GMREDI'
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT , myThid )

C- print out some kee parameters :
       CALL WRITE_0D_L( GM_AdvForm, INDEX_NONE,
     &  'GM_AdvForm =', '     /* if FALSE => use SkewFlux Form */')
       CALL WRITE_0D_L( GM_InMomAsStress, INDEX_NONE,
     &  'GM_InMomAsStress =', ' /* if TRUE => apply as Eddy Stress */')
       CALL WRITE_0D_L( GM_AdvSeparate, INDEX_NONE,
     & 'GM_AdvSeparate =',' /* Calc Bolus & Euler Adv. separately */')
       CALL WRITE_0D_L( GM_ExtraDiag, INDEX_NONE,
     &  'GM_ExtraDiag =','   /* Tensor Extra Diag (line 1&2) non 0 */')
       CALL WRITE_0D_RL( GM_isopycK, INDEX_NONE, 'GM_isopycK =',
     &  '    /* Background Isopyc. Diffusivity [m^2/s] */')
       CALL WRITE_0D_RL( GM_background_K*GM_skewflx, INDEX_NONE,
     &  'GM_skewflx*K =',
     &  '  /* Background GM_SkewFlx Diffusivity [m^2/s] */')
       CALL WRITE_0D_RL( GM_background_K*GM_advect, INDEX_NONE,
     &  'GM_advec*K =',
     &  '    /* Backg. GM-Advec(=Bolus) Diffusivity [m^2/s]*/')
       CALL WRITE_0D_RL( GM_Kmin_horiz, INDEX_NONE, 'GM_Kmin_horiz =',
     &  ' /* Minimum Horizontal Diffusivity [m^2/s] */')
       CALL WRITE_0D_RL( GEOM_alpha, INDEX_NONE,
     &  'GEOM_alpha =', ' /* GEOM alpha coeff. [-] */')
       CALL WRITE_0D_RL( GM_Small_Number, INDEX_NONE,
     &  'GM_Small_Number =', '  /* epsilon used in slope calc */')
       CALL WRITE_0D_RL( GM_slopeSqCutoff, INDEX_NONE,
     &  'GM_slopeSqCutoff =', ' /* Slope^2 cut-off value */')
       CALL WRITE_0D_C( GM_taper_scheme, 0, INDEX_NONE,
     &  'GM_taper_scheme =',
     &  '  /* Type of Tapering/Clipping scheme */')
       CALL WRITE_0D_RL( GM_maxSlope, INDEX_NONE,
     &  'GM_maxSlope =', '  /* Maximum Slope (Tapering/Clipping) */')
       CALL WRITE_0D_RL( GM_facTrL2dz, INDEX_NONE,
     &  'GM_facTrL2dz =',
     &  ' /* Minimum Trans.Layer Thick. (factor of dz) */')
       CALL WRITE_0D_RL( GM_facTrL2ML, INDEX_NONE,
     &  'GM_facTrL2ML =',
     &  ' /* Max.Trans.Layer Thick. (factor of MxL Depth)*/')
       CALL WRITE_0D_RL( GM_maxTransLay, INDEX_NONE,
     &  'GM_maxTransLay =',
     &  ' /* Maximum Transition Layer Thickness [m] */')
       CALL WRITE_0D_L( GM_UseBVP, INDEX_NONE,
     &  'GM_UseBVP =',
     &  ' /* if TRUE => use bvp a la Ferrari et al. (2010) */')
       CALL WRITE_0D_I( GM_BVP_ModeNumber, INDEX_NONE,
     &  'GM_BVP_ModeNumber =',
     &  ' /* Vertical mode number for BVP wave speed */')
       CALL WRITE_0D_RL( GM_BVP_cMin, INDEX_NONE,
     &  'GM_BVP_cMin =',
     &  ' /* Minimum wave speed for BVP [m/s] */')

C--  Check parameters:

C-     GM/Redi needs implicit diffusion (will be packaged later)
      IF ( .NOT.implicitDiffusion ) THEN
        WRITE(msgBuf,'(A)') 'GM/Redi needs implicitDiffusion=.true.'
        CALL PRINT_ERROR( msgBuf , myThid )
        STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF

#ifndef GM_GEOM_VARIABLE_K
C     Make sure we are not trying to use something that is unavailable
      IF ( GEOM_alpha.NE.0. ) THEN
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: GEOM variables used in data.gmredi'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: without #define GM_GEOM_VARIABLE_K'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF
#endif

#ifndef GM_BOLUS_ADVEC
C     Make sure we are not trying to use some arrays that are unavailable
      IF ( GM_AdvForm ) THEN
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: GM Advection form used in data.gmredi'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: without #define GM_BOLUS_ADVEC'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF
#endif

#ifndef GM_EXTRA_DIAGONAL
C     Make sure we are not trying to use some arrays that are unavailable
      IF ( GM_ExtraDiag ) THEN
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: GM_skew_Flux_K & GM_isopycK not equal'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: without #define GM_EXTRA_DIAGONAL'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF
#endif

#ifndef GM_NON_UNITY_DIAGONAL
      IF ( GM_iso2dFile .NE. ' ' .OR.
     &     GM_iso1dFile .NE. ' ' ) THEN
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: needs #define GM_NON_UNITY_DIAGONAL'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: to use GM_iso2dFile or GM_iso1dFile'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF
#endif

      IF ( GM_InMomAsStress ) THEN
#ifdef ALLOW_EDDYPSI
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: Bug in ALLOW_EDDYPSI with GM_InMomAsStress'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
#else /* ALLOW_EDDYPSI */
       WRITE(msgBuf,'(2A)')
     &  ' GMREDI_CHECK: need to define ALLOW_EDDYPSI in CPP_OPTIONS.h',
     &  ' to use GM_InMomAsStress'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
#endif /* ALLOW_EDDYPSI */
      ENDIF
      IF ( GM_InMomAsStress .AND. .NOT.GM_AdvForm ) THEN
       WRITE(msgBuf,'(A)')
     &   ' GMREDI_CHECK: need GM_AdvForm=T to use GM_InMomAsStress'
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R GMREDI_CHECK'
      ENDIF

#ifdef ALLOW_PTRACERS
      IF ( GM_AdvForm .AND. .NOT.GM_AdvSeparate
     &       .AND. usePTRACERS ) THEN
        redFlag = .FALSE.
        DO iTr=1,PTRACERS_numInUse
         IF ( .NOT.PTRACERS_useGMRedi(iTr) ) THEN
          redFlag = .TRUE.
          WRITE(msgBuf,'(2A,I3,A,L5)') ' GMREDI_CHECK:',
     &     ' pTracers_useGMRedi(',iTr,' )=', PTRACERS_useGMRedi(iTr)
          CALL PRINT_ERROR( msgBuf, myThid )
         ENDIF
        ENDDO
        IF ( redFlag ) THEN
          WRITE(msgBuf,'(2A)') ' GMREDI_CHECK:',
     &     ' but GM Advective Form applies to all tracers !'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R GMREDI_CHECK'
        ENDIF
      ENDIF
#endif /* ALLOW_PTRACERS */

#ifdef ALLOW_GENERIC_ADVDIFF
C     Check size of overlap region
      IF ( GM_AdvForm .AND. .NOT.GM_AdvSeparate
     &       .AND. GEOM_alpha.NE.0.
     &       .AND. useMultiDimAdvec ) THEN
C       variable K requires 1 more row/column in the overlap:
C       might need to increase Olx,Oly from 2 to 3 if GM advective
C       form & multi-dim advection are used. This happens when:
C       a) using a 5 points stencil advection scheme ; or
C       b) using a 3 points stencil advection scheme on CS-grid
        GAD_OlMinSize(2) = MAX( GAD_OlMinSize(2), 1)
        WRITE(msgBuf,'(A,9I3)')
     &      'GMREDI_CHECK: GAD_OlMinSize=', GAD_OlMinSize
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT , myThid )
      ENDIF
#endif /* ALLOW_GENERIC_ADVDIFF */

      _END_MASTER(myThid)

#endif /* ALLOW_GMREDI */
      RETURN
      END

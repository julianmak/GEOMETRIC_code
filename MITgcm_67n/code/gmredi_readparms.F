#include "GMREDI_OPTIONS.h"

CBOP
C     !ROUTINE: GMREDI_READPARMS
C     !INTERFACE:
      SUBROUTINE GMREDI_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_READPARMS
C     | o Routine to initialize GM/Redi variables and constants.
C     *==========================================================*
C     | Initialize GM/Redi parameters, read in data.gmredi
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GMREDI.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
      INTEGER myThid

#ifdef ALLOW_GMREDI
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
CEOP

C--   GM/Redi parameter
C     GM_Small_Number  :: epsilon used in computing the slope
C     GM_slopeSqCutoff :: slope^2 cut-off value
      NAMELIST /GM_PARM01/
     &          GM_AdvForm, GM_AdvSeparate,
     &          GM_InMomAsStress,
     &          GM_isopycK,
     &          GM_background_K,
     &          GM_iso2dFile, GM_iso1dFile,
     &          GM_bol2dFile, GM_bol1dFile,
     &          GM_background_K3dFile,
     &          GM_isopycK3dFile,
     &          GM_taper_scheme,
     &          GM_maxSlope,
     &          GM_Kmin_horiz,
     &          GM_Small_Number, GM_slopeSqCutoff,
     &          GM_Visbeck_alpha, GM_Visbeck_length,
     &          GM_Visbeck_depth,
     &          GM_Visbeck_minDepth, GM_Visbeck_maxSlope,
     &          GM_Visbeck_minVal_K, GM_Visbeck_maxVal_K,
     &          GM_facTrL2dz, GM_facTrL2ML, GM_maxTransLay,
     &          GM_Scrit, GM_Sd,
     &          GM_MNC,
     &          GM_UseBVP, GM_BVP_cMin, GM_BVP_ModeNumber,
     &          GM_useSubMeso, subMeso_Ceff, subMeso_invTau,
     &          subMeso_LfMin, subMeso_Lmax, GM_K3D_b1,
     &          GM_K3D_gamma, GM_useK3D, GM_K3D_vecFreq,
     &          GM_K3D_EadyMinDepth, GM_K3D_EadyMaxDepth,
     &          GM_K3D_Lambda, GM_K3D_smallK, GM_K3D_surfK,
     &          GM_K3D_maxC, GM_K3D_minCori, GM_K3D_minN2,
     &          GM_K3D_surfMinDepth, GM_maxK3D, GM_K3D_Rmax,
     &          GM_K3D_constK, GM_K3D_use_constK, GM_K3D_smooth,
     &          GM_K3D_beta_eq_0, GM_K3D_Rmin, GM_K3D_ThickSheet,
     &          GM_K3D_constRedi, GM_K3D_minRenorm,
     &          GM_K3D_maxRenorm, GM_useLeithQG,
     &          GM_useGEOM,
     &          GEOM_lmbda, GEOM_alpha, ene_init,
     &          ene_kappa, ene_local,
     &          vert_struc, vert_struc_min, vert_struc_max,
     &          GEOM_minval_K, GEOM_maxval_K,
     &          GEOM_pickup_write_mdsio, GEOM_pickup_write_mnc,
     &          GEOM_pickup_read_mdsio, GEOM_pickup_read_mnc

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useGMRedi ) THEN
C-    pkg GMREDI is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.gmredi is found
         CALL PACKAGES_UNUSED_MSG( 'useGMRedi', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default values GM/Redi
      GM_AdvForm          = .FALSE.
      GM_AdvSeparate      = .FALSE.
      GM_InMomAsStress    = .FALSE.
      GM_isopycK          = -999.
      GM_background_K     = 0. _d 0
      GM_maxSlope         = 1. _d -2
      GM_Kmin_horiz       = 0. _d 0
      GM_Small_Number     = 1. _d -20
      GM_slopeSqCutoff    = 1. _d +48
      GM_taper_scheme     = ' '
      GM_facTrL2dz        = 1.
      GM_facTrL2ML        = 5.
      GM_maxTransLay      = 500.
      GM_Scrit            = 0.004 _d 0
      GM_Sd               = 0.001 _d 0
      GM_MNC              = useMNC
      GM_iso2dFile        = ' '
      GM_iso1dFile        = ' '
      GM_bol2dFile        = ' '
      GM_bol1dFile        = ' '
      GM_background_K3dFile = ' '
      GM_isopycK3dFile    = ' '
      GM_useLeithQG       = .FALSE.

C--   Default values GM/Redi I/O control
c     GM_dumpFreq         = -1.
c     GM_taveFreq         = -1.

C--   Default values Visbeck
      GM_Visbeck_alpha    =    0. _d 0
      GM_Visbeck_length   =  200. _d 3
      GM_Visbeck_depth    = 1000. _d 0
      GM_Visbeck_minDepth =    0. _d 0
      GM_Visbeck_maxSlope = UNSET_RL
      GM_Visbeck_minVal_K =    0. _d 0
      GM_Visbeck_maxVal_K = 2500. _d 0

C--   Default values GEOMETRIC
      GM_useGEOM              = .false.
      GEOM_lmbda              = 2. _d -7
      GEOM_alpha              = 1. _d -1
      ene_init                = 1. _d -6
      ene_kappa               = 1. _d +3
      ene_local               = .false.
      vert_struc              = .false.
      vert_struc_min          = 1. _d -1
      vert_struc_max          = 1. _d 0
      GEOM_minVal_K           = 0. _d 0
      GEOM_maxVal_K           = 2500. _d 0
      GEOM_pickup_write_mdsio = .false.
      GEOM_pickup_write_mnc   = .false.
      GEOM_pickup_read_mdsio  = .false.
      GEOM_pickup_read_mnc    = .false.

C--   Default values BVP
      GM_UseBVP           = .FALSE.
      GM_BVP_ModeNumber   = 1
      GM_BVP_cMin         = 1. _d -1

C--   Default values Sub-Meso (corresponding tau = 5.8 days):
      GM_useSubMeso = .FALSE.
      subMeso_invTau =  2.0 _d -6
      subMeso_LfMin  =  1.0 _d +3
      subMeso_Ceff   =  7.0 _d -2
      subMeso_Lmax   = 110. _d +3

C--   Default values for K3D
      GM_useK3D        = .FALSE.
      GM_K3D_beta_eq_0 = .TRUE.
      GM_K3D_use_constK= .FALSE.
      GM_K3D_ThickSheet= .FALSE.
      GM_K3D_smooth    = .TRUE.
      GM_K3D_surfK     = .FALSE.
      GM_K3D_constRedi = .FALSE.
      GM_K3D_constK    = 1000.0
      GM_K3D_gamma     = 1.0
      GM_K3D_b1        = 4.0
      GM_K3D_EadyMinDepth = 50.0
      GM_K3D_EadyMaxDepth = 1000.0
      GM_K3D_Lambda    = 1.0
      GM_K3D_smallK    = 100.0
      GM_K3D_maxC      = 0.15
      GM_maxK3D        = 20e3
      GM_K3D_Rmin      = 30e3
      GM_K3D_Rmax      = 75e3
      GM_K3D_minCori   = zeroRL
      GM_K3D_minN2     = 1.0e-8
      GM_K3D_surfMinDepth = 100.0
      GM_K3D_vecFreq   = 2592000.0
      GM_K3D_minRenorm = oneRL
      GM_K3D_maxRenorm = 20.0

      WRITE(msgBuf,'(A)') ' GM_READPARMS: opening data.gmredi'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                          'data.gmredi', 'GM_READPARMS',
     O                          iUnit,
     I                          myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=GM_PARM01)
      WRITE(msgBuf,'(A)') ' GM_READPARMS: finished reading data.gmredi'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
C     Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C     Default value for GM_isopycK is equal to GM_background_K :
      IF (GM_isopycK.EQ.-999.) GM_isopycK = GM_background_K

C     Default value for GM_Visbeck_maxSlope is equal to GM_maxSlope :
      IF ( GM_Visbeck_maxSlope .EQ. UNSET_RL )
     &     GM_Visbeck_maxSlope = GM_maxSlope

C     Some constants
      GM_rMaxSlope = 0.
      if (GM_maxSlope.NE.0.) GM_rMaxSlope = 1. _d 0 / GM_maxSlope

      IF (GM_AdvForm) THEN
        GM_skewflx = 0.
        GM_advect  = 1.
        GM_ExtraDiag = GM_Visbeck_alpha.NE.0. .OR. GM_isopycK.NE.0.
     &                 .OR. GEOM_alpha.NE.0.
      ELSE
        GM_skewflx = 1.
        GM_advect  = 0.
        GM_ExtraDiag = GM_isopycK.NE.GM_background_K
      ENDIF
      IF ( GM_iso2dFile .NE. GM_bol2dFile .OR.
     &     GM_iso1dFile .NE. GM_bol1dFile ) THEN
        GM_ExtraDiag = .TRUE.
      ENDIF

      IF (GM_UseBVP) THEN
        GM_BVP_rModeNumber = 1. _d 0/DBLE(GM_BVP_ModeNumber)
        GM_BVP_cHat2Min = GM_BVP_cMin**2*rhoConst/gravity
      ENDIF

C     Make sure that we locally honor the global MNC on/off flag
      GM_MNC = GM_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      GM_MNC = .FALSE.
#endif
      GM_MDSIO = (.NOT. GM_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_GMREDI */

      RETURN
      END

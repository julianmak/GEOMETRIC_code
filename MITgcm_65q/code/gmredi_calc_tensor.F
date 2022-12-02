C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_calc_tensor.F,v 1.45 2015/10/15 23:06:36 gforget Exp $
C $Name:  $

#include "GMREDI_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_KPP
# include "KPP_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: GMREDI_CALC_TENSOR
C     !INTERFACE:
      SUBROUTINE GMREDI_CALC_TENSOR(
     I             iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_CALC_TENSOR
C     | o Calculate tensor elements for GM/Redi tensor.
C     *==========================================================*
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GMREDI.h"
#include "GMREDI_TAVE.h"
#ifdef ALLOW_CTRL
# include "CTRL_FIELDS.h"
#endif
#ifdef ALLOW_KPP
# include "KPP.h"
#endif

#ifdef ALLOW_AUTODIFF_TAMC
#include "tamc.h"
#include "tamc_keys.h"
#endif /* ALLOW_AUTODIFF_TAMC */

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi, bj    :: tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number in simulation
C     myThid    :: My Thread Id. number
C
      INTEGER iMin,iMax,jMin,jMax
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_GMREDI

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,k
      _RL SlopeX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL SlopeY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dSigmaDx(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dSigmaDy(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dSigmaDr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL SlopeSqr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL taperFct(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ldd97_LrhoC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ldd97_LrhoW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ldd97_LrhoS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Cspd, LrhoInf, LrhoSup, fCoriLoc
      _RL Kgm_tmp, isopycK, bolus_K

      INTEGER kLow_W (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER kLow_S (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL locMixLayer(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL baseSlope  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL hTransLay  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL recipLambda(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER  km1
#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
      INTEGER kp1
      _RL maskp1
#endif

#ifdef GM_GEOM_VARIABLE_K
      _RL dSigmaH, dSigmaR
      _RL Sloc, M2loc

      _RL recipMaxSlope
      _RL deltaH, integrDepth
      _RL N2loc, SNloc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- JM + JRM hacks
C     general variables
      integer, save :: energy_initialised = 0
      logical, save :: length_initialised = .false.
      _RL, save :: surface_area, zonal_length, domain_vol,
     &   max_depth
      _RL :: energy_int
      _RL :: Utrans_tot, Utrans_bt, Utrans_bc, energy_avg
      _RL, dimension(:, :), allocatable, save :: integrated_SNloc, 
     &   energy, energy_old, trd_energy_gen, kappa0
      _RL, dimension(:, :), allocatable, save :: depthW, depthS, 
     &   depthC

      integer :: unit
      _RL, save :: t_dum_init = 0.0 _d 0
     
C     variables for local energy budget
      _RL, dimension(:, :), allocatable, save :: 
     & ene_adv_x, ene_adv_y, trd_energy_adv, 
     & ene_lap_x, ene_lap_y, trd_energy_lap, 
     & U_depth_avg, V_depth_avg
      _RL, dimension(:, :, :), allocatable, save :: energy_rhs

C     structure function
      _RL, dimension(:, :, :), allocatable, save :: structure_func
C
C-- end JM + JRM hacks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
#endif /* GM_GEOM_VARIABLE_K */

#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  doDiagRediFlx
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
      _RL dTdz
      _RL tmp1k(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif
#endif /* ALLOW_DIAGNOSTICS */
      
#ifdef GM_GEOM_VARIABLE_K
      if(gm_restart) then
        t_dum_init = myTime
        print *, "JRM+JM mod = starting from a previous energy field"
        open(newunit = unit, file = "gm_checkpoint.dat",
     &    status = "old", access = "stream", action = "read")
        if(.not. allocated(energy_old)) then
          allocate(energy_old(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
          read(unit) energy_old
        end if
        if(energy_local) then
          if(.not. allocated(energy_rhs)) then
            allocate(energy_rhs(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 3))
            read(unit) energy_rhs
          end if
        end if
        read(unit) energy_initialised
        close(unit)
        gm_restart = .false.
      end if
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_AUTODIFF_TAMC
          act1 = bi - myBxLo(myThid)
          max1 = myBxHi(myThid) - myBxLo(myThid) + 1
          act2 = bj - myByLo(myThid)
          max2 = myByHi(myThid) - myByLo(myThid) + 1
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = ikey_dynamics - 1
          igmkey = (act1 + 1) + act2*max1
     &                      + act3*max1*max2
     &                      + act4*max1*max2*max3
#endif /* ALLOW_AUTODIFF_TAMC */

#ifdef ALLOW_DIAGNOSTICS
      doDiagRediFlx = .FALSE.
      IF ( useDiagnostics ) THEN
        doDiagRediFlx = DIAGNOSTICS_IS_ON('GM_KuzTz', myThid )
        doDiagRediFlx = doDiagRediFlx .OR.
     &                  DIAGNOSTICS_IS_ON('GM_KvzTz', myThid )
      ENDIF
#endif

#ifdef GM_GEOM_VARIABLE_K
      recipMaxSlope = 0. _d 0
      if(GM_maxSlope /= 0.0 _d 0) then
        recipMaxSlope = 1. _d 0 / abs(GM_maxSlope)
      end if
      GEOMK(:,:,:,bi,bj) = 0.0 _d 0
#endif

C--   set ldd97_Lrho (for tapering scheme ldd97):
      IF ( GM_taper_scheme.EQ.'ldd97' .OR.
     &     GM_taper_scheme.EQ.'fm07' ) THEN
       Cspd = 2. _d 0
       LrhoInf = 15. _d 3
       LrhoSup = 100. _d 3
C-     Tracer point location (center):
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         IF (fCori(i,j,bi,bj).NE.0.) THEN
           ldd97_LrhoC(i,j) = Cspd/ABS(fCori(i,j,bi,bj))
         ELSE
           ldd97_LrhoC(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoC(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoC(i,j),LrhoSup))
        ENDDO
       ENDDO
C-     U point location (West):
       DO j=1-OLy,sNy+OLy
        kLow_W(1-OLx,j) = 0
        ldd97_LrhoW(1-OLx,j) = LrhoSup
        DO i=1-OLx+1,sNx+OLx
         kLow_W(i,j) = MIN(kLowC(i-1,j,bi,bj),kLowC(i,j,bi,bj))
         fCoriLoc = op5*(fCori(i-1,j,bi,bj)+fCori(i,j,bi,bj))
         IF (fCoriLoc.NE.0.) THEN
           ldd97_LrhoW(i,j) = Cspd/ABS(fCoriLoc)
         ELSE
           ldd97_LrhoW(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoW(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoW(i,j),LrhoSup))
        ENDDO
       ENDDO
C-     V point location (South):
       DO i=1-OLx+1,sNx+OLx
         kLow_S(i,1-OLy) = 0
         ldd97_LrhoS(i,1-OLy) = LrhoSup
       ENDDO
       DO j=1-OLy+1,sNy+OLy
        DO i=1-OLx,sNx+OLx
         kLow_S(i,j) = MIN(kLowC(i,j-1,bi,bj),kLowC(i,j,bi,bj))
         fCoriLoc = op5*(fCori(i,j-1,bi,bj)+fCori(i,j,bi,bj))
         IF (fCoriLoc.NE.0.) THEN
           ldd97_LrhoS(i,j) = Cspd/ABS(fCoriLoc)
         ELSE
           ldd97_LrhoS(i,j) = LrhoSup
         ENDIF
         ldd97_LrhoS(i,j) = MAX(LrhoInf,MIN(ldd97_LrhoS(i,j),LrhoSup))
        ENDDO
       ENDDO
      ELSE
C-    Just initialize to zero (not use anyway)
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          ldd97_LrhoC(i,j) = 0. _d 0
          ldd97_LrhoW(i,j) = 0. _d 0
          ldd97_LrhoS(i,j) = 0. _d 0
        ENDDO
       ENDDO
      ENDIF

#ifdef GM_BOLUS_ADVEC
      DO k=1,Nr
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         GM_PsiX(i,j,k,bi,bj)  = 0. _d 0
         GM_PsiY(i,j,k,bi,bj)  = 0. _d 0
        ENDDO
       ENDDO
      ENDDO
#endif /* GM_BOLUS_ADVEC */
#ifdef ALLOW_AUTODIFF
      DO k=1,Nr
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         Kwx(i,j,k,bi,bj)  = 0. _d 0
         Kwy(i,j,k,bi,bj)  = 0. _d 0
         Kwz(i,j,k,bi,bj)  = 0. _d 0
# ifdef GM_NON_UNITY_DIAGONAL
         Kux(i,j,k,bi,bj)  = 0. _d 0
         Kvy(i,j,k,bi,bj)  = 0. _d 0
# endif
# ifdef GM_EXTRA_DIAGONAL
         Kuz(i,j,k,bi,bj)  = 0. _d 0
         Kvz(i,j,k,bi,bj)  = 0. _d 0
# endif
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_AUTODIFF */

C--   Initialise Mixed Layer related array:
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         hTransLay(i,j) = R_low(i,j,bi,bj)
         baseSlope(i,j) =  0. _d 0
         recipLambda(i,j) = 0. _d 0
         locMixLayer(i,j) = 0. _d 0
       ENDDO
      ENDDO
#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         locMixLayer(i,j) = KPPhbl(i,j,bi,bj)
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         locMixLayer(i,j) = hMixLayer(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- JM + JRM hacks
C      
C Notes: the lat-lon energy field is depth integrated
C        energy generation SN_loc is depth integrated
C        advecting flow (local) is depth AVERAGED
C
#ifdef GM_GEOM_VARIABLE_K
C     allocate and initialise
      if(.not. allocated(integrated_SNloc)) then
       allocate(integrated_SNloc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
      end if
      if(.not. allocated(energy_old)) then
       allocate(energy_old(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
        energy_old = 0.0 _d 0
      end if
      if(.not. allocated(depthW)) then
       allocate(depthW(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
       depthW = 0.0 _d 0
      end if
      if(.not. allocated(depthS)) then
       allocate(depthS(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
       depthS = 0.0 _d 0
      end if
      if(.not. allocated(depthC)) then
       allocate(depthC(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
       depthC = 0.0 _d 0
      end if
C     initialise some lengths
      if(.not. length_initialised) then
C     general surface area
       surface_area = 0.0 _d 0
       do j=1,sNy
        do i=1,sNx
         if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
          surface_area = surface_area + 
     &                   dxG(i,j,bi,bj) * dyG(i,j,bi,bj)
         end if
        end do
       end do
C     general zonal length (for rough idea on transport)
       zonal_length = 0.0 _d 0
       do i=1,sNx 
        if(maskC(i,1,1,bi,bj) /= 0.0 _d 0) then
         zonal_length = zonal_length + dxG(i,1,bi,bj)
        end if
       end do
C     domain volume
       domain_vol = 0.0 _d 0
       max_depth = 0.0 _d 0
       do k= 1,Nr
        do j=1,sNy
         do i=1,sNx
          if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
           domain_vol = domain_vol 
     &                + dxG(i,j,bi,bj) * dyG(i,j,bi,bj)
     &                * drF(k) * hFacC(i,j,k,bi,bj)
           depthC(i,j) = depthC(i,j) + 
     &                   drF(k) * hFacC(i,j,k,bi,bj)
          end if
         end do
        end do
        max_depth = max_depth + drF(k)
       end do
C     depths for doing depth-averaging
       do k=1,Nr
        do j=1-OLy, sNy+OLy
         do i=1-OLx+1, sNx+OLx
          if(maskW(i,j,k,bi,bj) /= 0.0 _d 0) then
           depthW(i,j) = depthW(i,j) +
     &                   drF(k) * hFacW(i,j,k,bi,bj)
          end if
         end do
        end do
        do j=1-OLy+1, sNy+OLy
         do i=1-OLx, sNx+OLx
          if(maskS(i,j,k,bi,bj) /= 0.0 _d 0) then
           depthS(i,j) = depthS(i,j) +
     &                   drF(k) * hFacS(i,j,k,bi,bj)
          end if
         end do
        end do
       end do
       
       length_initialised = .true.
      end if
C     read in domain-avg energy and turn into depth-integrated energy
      if(energy_initialised == 0) then
       if(energy_local) then
        energy_old = energy_init * max_depth
       else
        energy_old = energy_init * domain_vol
       end if
      end if
      integrated_SNloc = 0.0 _d 0

      if(.not. allocated(kappa0)) then
       allocate(kappa0(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
      end if
      kappa0 = 0.0 _d 0

C     allocate the structure function and do a loop on this first
      if(.not. allocated(structure_func)) then
       allocate(structure_func(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr))
      end if
      structure_func = 1.0 _d 0
      
      do k=Nr,2,-1
       do j=1-OLy+1,sNy+OLy-1
        do i=1-OLx+1,sNx+OLx-1
         dSigmaDx(i,j) = op25 * (sigmaX(i+1,j,k-1) + sigmaX(i,j,k-1)
     &                        +  sigmaX(i+1,j,k  ) + sigmaX(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
         dSigmaDy(i,j) = op25 * (sigmaY(i,j+1,k-1) + sigmaY(i,j,k-1)
     &                        +  sigmaY(i,j+1,k  ) + sigmaY(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
        end do
       end do
       do j=1-OLy,sNy+OLy
        do i=1-OLx,sNx+OLx
         dSigmaDr(i,j) = min(sigmaR(i,j,k), 0.0 _d 0)
        end do
       end do
       
C      NOTE: Ignores boundary cells (convenient with the 5-point filter used)
       do j=1-OLy+1, sNy+OLy-1
        do i=1-OLx+1, sNx+OLx-1
         if(maskC(i,j,k,bi,bj) /= 0.0 _d 0) then

C      Compute -N^2 rho_0 / g via averaging, on vertical interfaces
          dSigmaR = (dSigmaDr(i,j) * 4.0 _d 0
     &             + maskC(i-1,j,k,bi,bj)*dSigmaDr(i-1,j)
     &             + maskC(i+1,j,k,bi,bj)*dSigmaDr(i+1,j)
     &             + maskC(i,j-1,k,bi,bj)*dSigmaDr(i,j-1)
     &             + maskC(i,j+1,k,bi,bj)*dSigmaDr(i,j+1)
     &              ) / (4.0 _d 0
     &                 + maskC(i-1,j,k,bi,bj)
     &                 + maskC(i+1,j,k,bi,bj)
     &                 + maskC(i,j-1,k,bi,bj)
     &                 + maskC(i,j+1,k,bi,bj)
     &                  )
C      Compute M^2 rho_0 / g on vertical interfaces
          dSigmaH = sqrt(dSigmaDx(i,j) * dSigmaDx(i,j)
     &                 + dSigmaDy(i,j) * dSigmaDy(i,j))
          if(dSigmaH > 0.0 _d 0) then
           if(abs(dSigmaR) > dSigmaH*recipMaxSlope) then
            Sloc = dSigmaH / abs(dSigmaR)
           else
            Sloc = abs(GM_maxSlope)
           end if
           N2loc = abs(gravity * recip_rhoConst * dSigmaR)
           SNloc = Sloc * sqrt(N2loc)
          else
           SNloc = 0.0 _d 0
          end if
          if(vert_structure) then
           structure_func(i,j,k) = N2loc
          else
           structure_func(i,j,k) = 1.0 _d 0
          end if
         end if
        end do
       end do   
      end do

C     beginning of first k loop
      do k=Nr,2,-1
       do j=1-OLy+1,sNy+OLy-1
        do i=1-OLx+1,sNx+OLx-1
         dSigmaDx(i,j) = op25 * (sigmaX(i+1,j,k-1) + sigmaX(i,j,k-1)
     &                        +  sigmaX(i+1,j,k  ) + sigmaX(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
         dSigmaDy(i,j) = op25 * (sigmaY(i,j+1,k-1) + sigmaY(i,j,k-1)
     &                        +  sigmaY(i,j+1,k  ) + sigmaY(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
        end do
       end do
       do j=1-OLy,sNy+OLy
        do i=1-OLx,sNx+OLx
         dSigmaDr(i,j) = min(sigmaR(i,j,k), 0.0 _d 0)
        end do
       end do
       
C      NOTE: Ignores boundary cells (convenient with the 5-point filter used)
       do j=1-OLy+1, sNy+OLy-1
        do i=1-OLx+1, sNx+OLx-1
         if(maskC(i,j,k,bi,bj) /= 0.0 _d 0) then

C      Compute -N^2 rho_0 / g via averaging, on vertical interfaces
          dSigmaR = (dSigmaDr(i,j) * 4.0 _d 0
     &             + maskC(i-1,j,k,bi,bj)*dSigmaDr(i-1,j)
     &             + maskC(i+1,j,k,bi,bj)*dSigmaDr(i+1,j)
     &             + maskC(i,j-1,k,bi,bj)*dSigmaDr(i,j-1)
     &             + maskC(i,j+1,k,bi,bj)*dSigmaDr(i,j+1)
     &              ) / (4.0 _d 0
     &                 + maskC(i-1,j,k,bi,bj)
     &                 + maskC(i+1,j,k,bi,bj)
     &                 + maskC(i,j-1,k,bi,bj)
     &                 + maskC(i,j+1,k,bi,bj)
     &                  )
C      Compute M^2 rho_0 / g on vertical interfaces
          dSigmaH = sqrt(dSigmaDx(i,j) * dSigmaDx(i,j)
     &                 + dSigmaDy(i,j) * dSigmaDy(i,j))
          if(dSigmaH > 0.0 _d 0) then
           if(abs(dSigmaR) > dSigmaH*recipMaxSlope) then
            Sloc = dSigmaH / abs(dSigmaR)
           else
            Sloc = abs(GM_maxSlope)
           end if
           N2loc = abs(gravity * recip_rhoConst * dSigmaR)
           SNloc = Sloc * sqrt(N2loc)
          else
           SNloc = 0.0 _d 0
          end if
          if(vert_structure) then
C      Cap the (N^2 / N^2_surf) between something (1 and 0.1 default)
           structure_func(i,j,k) = max(min(vert_structure_max, 
     &                structure_func(i,j,k) / structure_func(i,j,2)),
     &                vert_structure_min)
          else
           structure_func(i,j,k) = 1.0 _d 0
          end if
          integrated_SNloc(i,j) = integrated_SNloc(i,j) 
     &                          + SNloc * drC(k) * hFacC(i,j,k,bi,bj)
     &                          * structure_func(i,j,k)
         end if
        end do
       end do   
      end do
      
      if(.not. energy_local) then
C      Horizontal integral
       SNloc = 0.0 _d 0
       do j=1,sNy
        do i=1,sNx
         if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
          SNLoc = SNLoc +
     &     integrated_SNloc(i,j) * dxG(i,j,bi,bj) * dyG(i,j,bi,bj)
         end if
        end do
       end do
       integrated_SNloc = SNloc / surface_area
      end if
#endif /* GM_GEOM_VARIABLE_K */
C
C-- end JM + JRM hacks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef GM_K3D
      IF (GM_useK3D) THEN
C     Calculate the 3D diffusivity as per Bates et al. (2013)
        CALL TIMER_START('GMREDI_K3D      [GMREDI_CALC_TENSOR]',
     &                    myThid)

        CALL GMREDI_K3D(
     I       iMin, iMax, jMin, jMax,
     I       sigmaX, sigmaY, sigmaR,
     I       bi, bj, myTime, myThid)

        CALL TIMER_STOP('GMREDI_K3D      [GMREDI_CALC_TENSOR]',
     &                  myThid)
      ENDIF
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 1rst loop on k : compute Tensor Coeff. at W points.

      DO k=Nr,2,-1

#ifdef ALLOW_AUTODIFF
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         SlopeX(i,j)       = 0. _d 0
         SlopeY(i,j)       = 0. _d 0
         dSigmaDx(i,j)     = 0. _d 0
         dSigmaDy(i,j)     = 0. _d 0
         dSigmaDr(i,j)     = 0. _d 0
         SlopeSqr(i,j)     = 0. _d 0
         taperFct(i,j)     = 0. _d 0
        ENDDO
       ENDDO
#endif /* ALLOW_AUTODIFF */

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
C      Gradient of Sigma at rVel points
         dSigmaDx(i,j)=op25*( sigmaX(i+1,j,k-1)+sigmaX(i,j,k-1)
     &                       +sigmaX(i+1,j, k )+sigmaX(i,j, k )
     &                      )*maskC(i,j,k,bi,bj)
         dSigmaDy(i,j)=op25*( sigmaY(i,j+1,k-1)+sigmaY(i,j,k-1)
     &                       +sigmaY(i,j+1, k )+sigmaY(i,j, k )
     &                      )*maskC(i,j,k,bi,bj)
c        dSigmaDr(i,j)=sigmaR(i,j,k)
        ENDDO
       ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- JM + JRM hacks
C
#ifdef GM_GEOM_VARIABLE_K
       do j=1-OLy,sNy+OLy
        do i=1-OLx,sNx+OLx
         if(maskC(i,j,k,bi,bj) /= 0.0 _d 0) then
          if(GEOM_const_var) then
           GEOMK(i,j,k,bi,bj) = structure_func(i,j,k)
     &                           * GEOM_kappa0_const
           kappa0(i,j) = GEOM_kappa0_const
          else
           GEOMK(i,j,k,bi,bj) = structure_func(i,j,k) * GEOM_alpha
     &      * energy_old(i,j) / 
     &        max(integrated_SNloc(i,j), epsilon(integrated_SNloc(i,j)))
           kappa0(i,j) = GEOM_alpha
     &      * energy_old(i,j) / 
     &        max(integrated_SNloc(i,j), epsilon(integrated_SNloc(i,j)))
           if(.not. energy_local) then
            GEOMK(i,j,k,bi,bj) = GEOMK(i,j,k,bi,bj) / surface_area
            kappa0(i,j) = kappa0(i,j) / surface_area
           end if
          end if
         end if
        end do
       end do
#endif /* GM_GEOM_VARIABLE_K */
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         dSigmaDr(i,j)=sigmaR(i,j,k)
        ENDDO
       ENDDO
C
C-- end JM + JRM hacks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoC, locMixLayer, rF,
     I             kLowC(1-OLx,1-OLy,bi,bj),
     I             k, bi, bj, myTime, myIter, myThid )

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
C      Mask Iso-neutral slopes
         SlopeX(i,j)=SlopeX(i,j)*maskC(i,j,k,bi,bj)
         SlopeY(i,j)=SlopeY(i,j)*maskC(i,j,k,bi,bj)
         SlopeSqr(i,j)=SlopeSqr(i,j)*maskC(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeX(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE SlopeY(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE SlopeSqr(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C      Components of Redi/GM tensor
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
          Kwx(i,j,k,bi,bj)= SlopeX(i,j)*taperFct(i,j)
          Kwy(i,j,k,bi,bj)= SlopeY(i,j)*taperFct(i,j)
          Kwz(i,j,k,bi,bj)= SlopeSqr(i,j)*taperFct(i,j)
        ENDDO
       ENDDO

C-- end 1rst loop on vertical level index k
      ENDDO

#ifdef GM_GEOM_VARIABLE_K
      if(GEOM_alpha /= 0.0 _d 0) then
C   FIXME: k=1 is surface but never used, copy one in for now
       GEOMK(:,:,1,bi,bj) = GEOMK(:,:,2,bi,bj)
       GEOMK(:,:,:,bi,bj) = min(max(GEOMK(:,:,:,bi,bj),
     &                                 GEOM_minVal_K),
     &                              GEOM_maxVal_K)
      end if
#endif /* GM_GEOM_VARIABLE_K */
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- JM + JRM hacks
C   time step the energy
#ifdef GM_GEOM_VARIABLE_K
C     allocate and initialise
      if(.not. allocated(energy)) then
       allocate(energy(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
       energy = 0.0 _d 0
      end if
      if(.not. allocated(trd_energy_gen)) then
       allocate(trd_energy_gen(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
      end if
      trd_energy_gen = 0.0 _d 0

C     loop over k
      do k=Nr,2,-1
       do j=1-OLy+1,sNy+OLy-1
        do i=1-OLx+1,sNx+OLx-1
         dSigmaDx(i,j) = op25 * (sigmaX(i+1,j,k-1) + sigmaX(i,j,k-1)
     &                        +  sigmaX(i+1,j,k  ) + sigmaX(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
         dSigmaDy(i,j) = op25 * (sigmaY(i,j+1,k-1) + sigmaY(i,j,k-1)
     &                        +  sigmaY(i,j+1,k  ) + sigmaY(i,j,k  )
     &                      ) *  maskC(i,j,k,bi,bj)
        end do
       end do
       do j=1-OLy,sNy+OLy
        do i=1-OLx,sNx+OLx
         dSigmaDr(i,j) = min(sigmaR(i,j,k), 0.0 _d 0)
        end do
       end do
       
C      NOTE: Ignores boundary cells (convenient with the 5-point filter used)
       do j=1-OLy+1, sNy+OLy-1
        do i=1-OLx+1, sNx+OLx-1
         if(maskC(i,j,k,bi,bj) /= 0.0 _d 0) then

C      Compute -N^2 rho_0 / g via averaging, on vertical interfaces
          dSigmaR = (dSigmaDr(i,j) * 4.0 _d 0
     &             + maskC(i-1,j,k,bi,bj)*dSigmaDr(i-1,j)
     &             + maskC(i+1,j,k,bi,bj)*dSigmaDr(i+1,j)
     &             + maskC(i,j-1,k,bi,bj)*dSigmaDr(i,j-1)
     &             + maskC(i,j+1,k,bi,bj)*dSigmaDr(i,j+1)
     &              ) / (4.0 _d 0
     &                 + maskC(i-1,j,k,bi,bj)
     &                 + maskC(i+1,j,k,bi,bj)
     &                 + maskC(i,j-1,k,bi,bj)
     &                 + maskC(i,j+1,k,bi,bj)
     &                  )
C      Compute M^2 rho_0 / g on vertical interfaces
          dSigmaH = sqrt(dSigmaDx(i,j) * dSigmaDx(i,j)
     &                 + dSigmaDy(i,j) * dSigmaDy(i,j))
          if(dSigmaH > 0.0 _d 0) then
           if(abs(dSigmaR) > dSigmaH*recipMaxSlope) then
            Sloc = dSigmaH / abs(dSigmaR)
           else
            Sloc = abs(GM_maxSlope)
           end if
           N2loc = abs(gravity * recip_rhoConst * dSigmaR)
           SNloc = Sloc * sqrt(N2loc)
          else
           SNloc = 0.0 _d 0
          end if
          trd_energy_gen(i,j) = trd_energy_gen(i,j)
     &                        + GEOMK(i,j,k,bi,bj) 
     &                        * SNloc * SNloc 
     &                        * drC(k) * hFacC(i,j,k,bi,bj)
         end if
        end do
       end do
      end do
C ---------------------------------------------- 
C ---- start advection of energy ---------------
C ---------------------------------------------- 
      if(energy_local) then
       if(.not. allocated(energy_rhs)) then
         allocate(energy_rhs(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 3))
         energy_rhs = 0.0 _d 0
       end if
       if(.not. allocated(trd_energy_adv)) then
         allocate(trd_energy_adv(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         trd_energy_adv = 0.0 _d 0
       end if
       if(.not. allocated(trd_energy_lap)) then
         allocate(trd_energy_lap(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         trd_energy_lap = 0.0 _d 0
       end if
       if(.not. allocated(ene_adv_x)) then
         allocate(ene_adv_x(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         ene_adv_x = 0.0 _d 0
       end if
       if(.not. allocated(ene_adv_y)) then
         allocate(ene_adv_y(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         ene_adv_y = 0.0 _d 0
       end if
       if(.not. allocated(ene_lap_x)) then
         allocate(ene_lap_x(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         ene_lap_x = 0.0 _d 0
       end if
       if(.not. allocated(ene_lap_y)) then
         allocate(ene_lap_y(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
         ene_lap_y = 0.0 _d 0
       end if
       
       energy_rhs(:, :, 1) = energy_rhs(:, :, 2)
       energy_rhs(:, :, 2) = energy_rhs(:, :, 3)

C ---- form the energy fluxes
C      first form the depth-avg velocity 
       if(.not. allocated(U_depth_avg)) then
        allocate(U_depth_avg(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
        U_depth_avg = 0.0 _d 0
       end if
       if(.not. allocated(V_depth_avg)) then
        allocate(V_depth_avg(1-OLx:sNx+OLx, 1-OLy:sNy+OLy))
        V_depth_avg = 0.0 _d 0
       end if

C      vertically integrate first   
       do k=Nr,2,-1     
        do j=1-OLy, sNy+OLy
         do i=1-OLx+1, sNx+OLx
          if(maskW(i,j,k,bi,bj) /= 0.0 _d 0) then
           U_depth_avg(i, j) = U_depth_avg(i,j) +
     &                         uVel(i,j,k,bi,bj) * drF(k)
     &                         * hFacW(i,j,k,bi,bj)
          end if
         end do
        end do
        do j=1-OLy+1, sNy+OLy
         do i=1-OLx, sNx+OLx
          if(maskS(i,j,k,bi,bj) /= 0.0 _d 0) then
           V_depth_avg(i, j) = V_depth_avg(i,j) +
     &                         vVel(i,j,k,bi,bj) * drF(k)
     &                         * hFacS(i,j,k,bi,bj)
          end if
         end do
        end do
       end do
C      depth average it
       do j=1-OLy, sNy+OLy
        do i=1-OLx+1, sNx+OLx
         if(maskW(i,j,1,bi,bj) /= 0.0 _d 0) then
          U_depth_avg(i,j) = U_depth_avg(i,j) / depthW(i,j)
         end if
        end do
       end do
       do j=1-OLy+1, sNy+OLy
        do i=1-OLx, sNx+OLx
         if(maskS(i,j,1,bi,bj) /= 0.0 _d 0) then
          V_depth_avg(i,j) = V_depth_avg(i,j) / depthS(i,j)
         end if
        end do
       end do
       
       do j=1-OLy, sNy+OLy
        do i=1-OLx+1, sNx+OLx
         if(    maskC(i-1,j,k,bi,bj) /= 0.0 _d 0
     &    .and. maskC(i  ,j,k,bi,bj) /= 0.0 _d 0) then
C ---- second ordered centre difference
          ene_adv_x(i, j) = -U_depth_avg(i,j) *
     &     (dxG(i-1,j,bi,bj)*energy_old(i,j)+
     &      dxG(i,j,bi,bj)*energy_old(i-1,j))/
     &     (dxG(i,j,bi,bj)+dxG(i-1,j,bi,bj)) 
          ene_lap_x(i, j) = energy_kappa *
     &      (2.0 _d 0) * (energy_old(i,j)
     &                   -energy_old(i-1,j))/
     &     (dxG(i,j,bi,bj)+dxG(i-1,j,bi,bj))
         end if
        end do
       end do
       ene_adv_x(-OLx+1:0, :) = ene_adv_x(sNx+1-OLx:sNx, :)
       ene_adv_x(sNx+1:sNx+OLx, :) = ene_adv_x(1:OLx, :)
       ene_lap_x(-OLx+1:0, :) = ene_lap_x(sNx+1-OLx:sNx, :)
       ene_lap_x(sNx+1:sNx+OLx, :) = ene_lap_x(1:OLx, :)
       
       do j=1-OLy+1, sNy+OLy
        do i=1-OLx, sNx+OLx
         if(    maskC(i,j-1,k,bi,bj) /= 0.0 _d 0
     &    .and. maskC(i,j  ,k,bi,bj) /= 0.0 _d 0) then
C ---- second ordered centre difference  
          ene_adv_y(i, j) = -V_depth_avg(i,j) *
     &     (dyG(i,j-1,bi,bj)*energy_old(i,j)+
     &      dyG(i,j,bi,bj)*energy_old(i,j-1))/
     &     (dyG(i,j,bi,bj)+dyG(i,j-1,bi,bj)) 
          ene_lap_y(i, j) = energy_kappa * 
     &      (2.0 _d 0) * (energy_old(i,j)
     &                   -energy_old(i,j-1))/
     &     (dyG(i,j,bi,bj)+dyG(i,j-1,bi,bj))
         end if
        end do
       end do
       ene_adv_y(:, -OLy+1:0) = ene_adv_y(:, sNy+1-OLy:sNy)
       ene_adv_y(:, sNy+1:sNy+OLy) = ene_adv_y(:, 1:OLy)
       ene_lap_y(:, -OLy+1:0) = ene_lap_y(:, sNy+1-OLy:sNy)
       ene_lap_y(:, sNy+1:sNy+OLy) = ene_lap_y(:, 1:OLy)

       do j=1-OLy, sNy+OLy
        do i=1-OLx, sNx+OLx
         trd_energy_adv(i,j) = 
     &         (ene_adv_x(i+1,j)-ene_adv_x(i,j))/dxG(i,j,bi,bj)
     &        +(ene_adv_y(i,j+1)-ene_adv_y(i,j))/dyG(i,j,bi,bj)
C        add in a bit to take into account linear free surface         
c         if(implicitFreeSurface) then
c          if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
c           trd_energy_adv(i,j) = wVel(i,j,1,bi,bj) * energy_old(i,j)
c     &                         / depthC(i,j) 
c          end if
c         end if
         trd_energy_lap(i,j) = 
     &        (ene_lap_x(i+1,j)-ene_lap_x(i,j))/dxG(i,j,bi,bj)
     &       +(ene_lap_y(i,j+1)-ene_lap_y(i,j))/dyG(i,j,bi,bj)
         energy_rhs(i,j,3) = trd_energy_gen(i,j)
     &                     - GEOM_lambda * energy_old(i,j)
     &                     + trd_energy_adv(i,j)
     &                     + trd_energy_lap(i,j)
        end do
       end do
       
       select case(energy_initialised)
        case(0)
         ! AB1 (FE)
         energy = energy_old + deltaT * energy_rhs(:, :, 3)
         energy_initialised = 1
        case(1)
         ! AB2
         energy = energy_old + deltaT *
     &    (+(3.0 _d 0 / 2.0 _d 0) * energy_rhs(:, :, 3)
     &     -(1.0 _d 0 / 2.0 _d 0) * energy_rhs(:, :, 2))
         energy_initialised = 2
        case default
         ! AB3
         energy = energy_old + deltaT *
     &    (+(23.0 _d 0 / 12.0 _d 0) * energy_rhs(:, :, 3)
     &     -( 4.0 _d 0 /  3.0 _d 0) * energy_rhs(:, :, 2)
     &     +( 5.0 _d 0 / 12.0 _d 0) * energy_rhs(:, :, 1))
         energy_initialised = 3
       end select       
C ----------------------------------------------       
C ---- end advection of energy -----------------
C ---------------------------------------------- 
      else  
C      Horizontally integrate and time step the energy

C      form integrand
       energy_int = 0.0 _d 0
       do j=1,sNy
        do i=1,sNx
         if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
          energy_int = energy_int +
     &     trd_energy_gen(i,j) * dxG(i,j,bi,bj) * dyG(i,j,bi,bj)
         end if
        end do
       end do
C      time step and allocate it as a field with same entries
       energy = (deltaT * energy_int + energy_old(1,1))
     &   / (1.0 _d 0 + GEOM_lambda * deltaT)
      
       energy_initialised = 1
       if(minval(energy(1:sNx,1:sNy)) 
     &    /= maxval(energy(1:sNx,1:sNy))) then
         stop "Non constant energy in domain integrated case"
       end if
      end if
      energy_old = energy

#endif /* GM_GEOM_VARIABLE_K */
C
C-- end JM + JRM hacks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-    express the Tensor in term of Diffusivity (= m**2 / s )
      DO k=1,Nr
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
# if (defined (GM_NON_UNITY_DIAGONAL) || \
      defined (GM_GEOM_VARIABLE_K))
CADJ STORE Kwx(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE Kwy(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE Kwz(:,:,k,bi,bj) = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
       km1 = MAX(k-1,1)
       isopycK = GM_isopycK
     &         *(GM_isoFac1d(km1)+GM_isoFac1d(k))*op5
       bolus_K = GM_background_K
     &         *(GM_bolFac1d(km1)+GM_bolFac1d(k))*op5
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
         Kgm_tmp = kapRedi(i,j,k,bi,bj)
#  else
         Kgm_tmp = op5*(kapRedi(i,j,k,bi,bj)+kapRedi(i,j,km1,bi,bj))
#  endif
#else
         Kgm_tmp = isopycK*GM_isoFac2d(i,j,bi,bj)
#endif
#ifdef ALLOW_KAPGM_CONTROL
#  ifdef ALLOW_KAPGM_CONTROL_OLD
     &           + GM_skewflx*kapGM(i,j,k,bi,bj)
#  else
     &     + GM_skewflx*op5*(kapGM(i,j,k,bi,bj)+kapGM(i,j,km1,bi,bj))
#  endif
#else
     &           + GM_skewflx*bolus_K*GM_bolFac2d(i,j,bi,bj)
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &           + GEOMK(i,j,k,bi,bj)*GM_isoFac2d(i,j,bi,bj)
     &           + GM_skewflx*GEOMK(i,j,k,bi,bj)
     &            *GM_bolFac2d(i,j,bi,bj)
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &           + op5*(K3D(i,j,k,bi,bj)+K3D(i,j,km1,bi,bj))
     &             *(1. _d 0 + GM_skewflx)
#endif
         Kwx(i,j,k,bi,bj)= Kgm_tmp*Kwx(i,j,k,bi,bj)
         Kwy(i,j,k,bi,bj)= Kgm_tmp*Kwy(i,j,k,bi,bj)
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
         Kwz(i,j,k,bi,bj)= ( kapRedi(i,j,k,bi,bj)
#  else
         Kwz(i,j,k,bi,bj)= ( op5*(kapRedi(i,j,k,bi,bj)
     &                            +kapRedi(i,j,km1,bi,bj))
#  endif
#else
         Kwz(i,j,k,bi,bj)= ( isopycK*GM_isoFac2d(i,j,bi,bj)
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &                     + GEOMK(i,j,k,bi,bj)
C     &                      *GM_isoFac2d(i,j,bi,bj)
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &                     + op5*(K3D(i,j,k,bi,bj)+K3D(i,j,km1,bi,bj))
#endif
     &                     )*Kwz(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics .AND. GM_taper_scheme.EQ.'fm07' ) THEN
       CALL DIAGNOSTICS_FILL( hTransLay, 'GM_hTrsL', 0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL( baseSlope, 'GM_baseS', 0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(recipLambda,'GM_rLamb', 0,1,2,bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate Stream-Functions used in Advective Form:

#ifdef GM_BOLUS_ADVEC
      IF (GM_AdvForm) THEN
#ifdef GM_BOLUS_BVP
       IF (GM_UseBVP) THEN
        CALL GMREDI_CALC_PSI_BVP(
     I             bi, bj, iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             myThid )
       ELSE
#endif
#ifndef GM_K3D_PASSIVE
         IF (.NOT. GM_useK3D) THEN
#endif
C          If using GM_K3D PsiX and PsiY are calculated in gmredi_k3d
           CALL GMREDI_CALC_PSI_B(
     I              bi, bj, iMin, iMax, jMin, jMax,
     I              sigmaX, sigmaY, sigmaR,
     I              ldd97_LrhoW, ldd97_LrhoS,
     I              myThid )
#ifndef GM_K3D_PASSIVE
         ENDIF
#endif
#ifdef GM_BOLUS_BVP
       ENDIF
#endif
      ENDIF
#endif

#ifndef GM_EXCLUDE_SUBMESO
      IF ( GM_useSubMeso .AND. GM_AdvForm ) THEN
        CALL SUBMESO_CALC_PSI(
     I              bi, bj, iMin, iMax, jMin, jMax,
     I              sigmaX, sigmaY, sigmaR,
     I              locMixLayer,
     I              myIter, myThid )
      ENDIF
#endif /* ndef GM_EXCLUDE_SUBMESO */

#if ( defined (GM_NON_UNITY_DIAGONAL) || defined (GM_EXTRA_DIAGONAL) )
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 2nd  k loop : compute Tensor Coeff. at U point

#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=1-OLy,sNy+OLy
        DO i=2-OLx,sNx+OLx
         locMixLayer(i,j) = ( KPPhbl(i-1,j,bi,bj)
     &                      + KPPhbl( i ,j,bi,bj) )*op5
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=1-OLy,sNy+OLy
        DO i=2-OLx,sNx+OLx
         locMixLayer(i,j) = ( hMixLayer(i-1,j,bi,bj)
     &                      + hMixLayer( i ,j,bi,bj) )*op5
        ENDDO
       ENDDO
      ENDIF
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         hTransLay(i,j) =  0.
         baseSlope(i,j) =  0.
         recipLambda(i,j)= 0.
       ENDDO
       DO i=2-OLx,sNx+OLx
         hTransLay(i,j) = MAX( R_low(i-1,j,bi,bj), R_low(i,j,bi,bj) )
       ENDDO
      ENDDO

      DO k=Nr,1,-1
       kp1 = MIN(Nr,k+1)
       maskp1 = 1. _d 0
       IF (k.GE.Nr) maskp1 = 0. _d 0
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
#endif

C     Gradient of Sigma at U points
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
         dSigmaDx(i,j)=sigmaX(i,j,k)
     &                       *_maskW(i,j,k,bi,bj)
         dSigmaDy(i,j)=op25*( sigmaY(i-1,j+1,k)+sigmaY(i,j+1,k)
     &                       +sigmaY(i-1, j ,k)+sigmaY(i, j ,k)
     &                      )*_maskW(i,j,k,bi,bj)
         dSigmaDr(i,j)=op25*( sigmaR(i-1,j, k )+sigmaR(i,j, k )
     &                      +(sigmaR(i-1,j,kp1)+sigmaR(i,j,kp1))*maskp1
     &                      )*_maskW(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeSqr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE locMixLayer(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoW, locMixLayer, rC,
     I             kLow_W,
     I             k, bi, bj, myTime, myIter, myThid )

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeSqr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

#ifdef GM_NON_UNITY_DIAGONAL
c      IF ( GM_nonUnitDiag ) THEN
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kux(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
     &     ( kapRedi(i,j,k,bi,bj)
#  else
     &     ( op5*(kapRedi(i,j,k,bi,bj)+kapRedi(i-1,j,k,bi,bj))
#  endif
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &       +op5*(GEOMK(i,j,k,bi,bj)+GEOMK(i-1,j,k,bi,bj))
C     &        *GM_isoFac1d(k)
C     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
     &     )*taperFct(i,j)
         ENDDO
        ENDDO
#ifdef ALLOW_AUTODIFF_TAMC
# ifdef GM_EXCLUDE_CLIPPING
CADJ STORE Kux(:,:,k,bi,bj)  = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kux(i,j,k,bi,bj) = MAX( Kux(i,j,k,bi,bj), GM_Kmin_horiz )
         ENDDO
        ENDDO
c      ENDIF
#endif /* GM_NON_UNITY_DIAGONAL */

#ifdef GM_EXTRA_DIAGONAL

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeX(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */
       IF ( GM_ExtraDiag ) THEN
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kuz(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
     &     ( kapRedi(i,j,k,bi,bj)
#  else
     &     ( op5*(kapRedi(i,j,k,bi,bj)+kapRedi(i-1,j,k,bi,bj))
#  endif
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef ALLOW_KAPGM_CONTROL
#  ifdef ALLOW_KAPGM_CONTROL_OLD
     &     - GM_skewflx*kapGM(i,j,k,bi,bj)
#  else
     &     - GM_skewflx*op5*(kapGM(i,j,k,bi,bj)+kapGM(i-1,j,k,bi,bj))
#  endif
#else
     &     - GM_skewflx*GM_background_K*GM_bolFac1d(k)
     &        *op5*(GM_bolFac2d(i-1,j,bi,bj)+GM_bolFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &     +op5*(GEOMK(i,j,k,bi,bj)+GEOMK(i-1,j,k,bi,bj))
C     &         *GM_advect
     &     - GM_skewflx * GEOMK(i,j,k,bi,bj)
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &     +op5*(K3D(i,j,k,bi,bj)+K3D(i-1,j,k,bi,bj))*GM_advect
#endif
     &     )*SlopeX(i,j)*taperFct(i,j)
         ENDDO
        ENDDO
       ENDIF
#endif /* GM_EXTRA_DIAGONAL */

#ifdef ALLOW_DIAGNOSTICS
       IF (doDiagRediFlx) THEN
        km1 = MAX(k-1,1)
        DO j=1,sNy
         DO i=1,sNx+1
C         store in tmp1k Kuz_Redi
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
          tmp1k(i,j) = ( kapRedi(i,j,k,bi,bj)
#  else
          tmp1k(i,j) = ( op5*(kapRedi(i-1,j,k,bi,bj)
     &                       +kapRedi(i,j,k,bi,bj))
#  endif
#else
          tmp1k(i,j) = ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i-1,j,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &     +(GEOMK(i,j,k,bi,bj)+GEOMK(i-1,j,k,bi,bj))*0.5 _d 0
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &     +op5*(K3D(i,j,k,bi,bj)+K3D(i-1,j,k,bi,bj))
#endif
     &                 )*SlopeX(i,j)*taperFct(i,j)
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx+1
C-        Vertical gradients interpolated to U points
          dTdz = (
     &     +recip_drC(k)*
     &       ( maskC(i-1,j,k,bi,bj)*
     &           (theta(i-1,j,km1,bi,bj)-theta(i-1,j,k,bi,bj))
     &        +maskC( i ,j,k,bi,bj)*
     &           (theta( i ,j,km1,bi,bj)-theta( i ,j,k,bi,bj))
     &       )
     &     +recip_drC(kp1)*
     &       ( maskC(i-1,j,kp1,bi,bj)*
     &           (theta(i-1,j,k,bi,bj)-theta(i-1,j,kp1,bi,bj))
     &        +maskC( i ,j,kp1,bi,bj)*
     &           (theta( i ,j,k,bi,bj)-theta( i ,j,kp1,bi,bj))
     &       )      ) * 0.25 _d 0
           tmp1k(i,j) = dyG(i,j,bi,bj)*drF(k)
     &                * _hFacW(i,j,k,bi,bj)
     &                * tmp1k(i,j) * dTdz
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(tmp1k, 'GM_KuzTz', k,1,2,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C-- end 2nd  loop on vertical level index k
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-- 3rd  k loop : compute Tensor Coeff. at V point

#ifdef ALLOW_KPP
      IF ( useKPP ) THEN
       DO j=2-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         locMixLayer(i,j) = ( KPPhbl(i,j-1,bi,bj)
     &                      + KPPhbl(i, j ,bi,bj) )*op5
        ENDDO
       ENDDO
      ELSE
#else
      IF ( .TRUE. ) THEN
#endif
       DO j=2-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         locMixLayer(i,j) = ( hMixLayer(i,j-1,bi,bj)
     &                      + hMixLayer(i, j ,bi,bj) )*op5
        ENDDO
       ENDDO
      ENDIF
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         hTransLay(i,j) =  0.
         baseSlope(i,j) =  0.
         recipLambda(i,j)= 0.
       ENDDO
      ENDDO
      DO j=2-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         hTransLay(i,j) = MAX( R_low(i,j-1,bi,bj), R_low(i,j,bi,bj) )
       ENDDO
      ENDDO

C     Gradient of Sigma at V points
      DO k=Nr,1,-1
       kp1 = MIN(Nr,k+1)
       maskp1 = 1. _d 0
       IF (k.GE.Nr) maskp1 = 0. _d 0
#ifdef ALLOW_AUTODIFF_TAMC
       kkey = (igmkey-1)*Nr + k
#endif

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
         dSigmaDx(i,j)=op25*( sigmaX(i, j ,k) +sigmaX(i+1, j ,k)
     &                       +sigmaX(i,j-1,k) +sigmaX(i+1,j-1,k)
     &                      )*_maskS(i,j,k,bi,bj)
         dSigmaDy(i,j)=sigmaY(i,j,k)
     &                       *_maskS(i,j,k,bi,bj)
         dSigmaDr(i,j)=op25*( sigmaR(i,j-1, k )+sigmaR(i,j, k )
     &                      +(sigmaR(i,j-1,kp1)+sigmaR(i,j,kp1))*maskp1
     &                      )*_maskS(i,j,k,bi,bj)
        ENDDO
       ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE dSigmaDx(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDy(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE dSigmaDr(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE baseSlope(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE hTransLay(:,:)      = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE recipLambda(:,:)    = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

C     Calculate slopes for use in tensor, taper and/or clip
       CALL GMREDI_SLOPE_LIMIT(
     O             SlopeX, SlopeY,
     O             SlopeSqr, taperFct,
     U             hTransLay, baseSlope, recipLambda,
     U             dSigmaDr,
     I             dSigmaDx, dSigmaDy,
     I             ldd97_LrhoS, locMixLayer, rC,
     I             kLow_S,
     I             k, bi, bj, myTime, myIter, myThid )

cph(
#ifdef ALLOW_AUTODIFF_TAMC
cph(
CADJ STORE taperfct(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
cph)
#endif /* ALLOW_AUTODIFF_TAMC */
cph)

#ifdef GM_NON_UNITY_DIAGONAL
c      IF ( GM_nonUnitDiag ) THEN
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kvy(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
     &     ( kapRedi(i,j,k,bi,bj)
#  else
     &     ( op5*(kapRedi(i,j,k,bi,bj)+kapRedi(i,j-1,k,bi,bj))
#  endif
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &      +op5*((GEOMK(i,j,k,bi,bj)+GEOMK(i,j-1,k,bi,bj)))
C     &        *GM_isoFac1d(k)
C     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &     +op5*(K3D(i,j,k,bi,bj)+K3D(i,j-1,k,bi,bj))
#endif
     &     )*taperFct(i,j)
         ENDDO
        ENDDO
#ifdef ALLOW_AUTODIFF_TAMC
# ifdef GM_EXCLUDE_CLIPPING
CADJ STORE Kvy(:,:,k,bi,bj)  = comlev1_bibj_k, key=kkey, byte=isbyte
# endif
#endif
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kvy(i,j,k,bi,bj) = MAX( Kvy(i,j,k,bi,bj), GM_Kmin_horiz )
         ENDDO
        ENDDO
c      ENDIF
#endif /* GM_NON_UNITY_DIAGONAL */

#ifdef GM_EXTRA_DIAGONAL

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE SlopeY(:,:)       = comlev1_bibj_k, key=kkey, byte=isbyte
CADJ STORE taperFct(:,:)     = comlev1_bibj_k, key=kkey, byte=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */
       IF ( GM_ExtraDiag ) THEN
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          Kvz(i,j,k,bi,bj) =
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
     &     ( kapRedi(i,j,k,bi,bj)
#  else
     &     ( op5*(kapRedi(i,j,k,bi,bj)+kapRedi(i,j-1,k,bi,bj))
#  endif
#else
     &     ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef ALLOW_KAPGM_CONTROL
#  ifdef ALLOW_KAPGM_CONTROL_OLD
     &     - GM_skewflx*kapGM(i,j,k,bi,bj)
#  else
     &     - GM_skewflx*op5*(kapGM(i,j,k,bi,bj)+kapGM(i,j-1,k,bi,bj))
#  endif
#else
     &     - GM_skewflx*GM_background_K*GM_bolFac1d(k)
     &        *op5*(GM_bolFac2d(i,j-1,bi,bj)+GM_bolFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &     +op5*(GEOMK(i,j,k,bi,bj)+GEOMK(i,j-1,k,bi,bj))
C     &         *GM_advect
     &     - GM_skewflx * GEOMK(i,j,k,bi,bj)
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &     +op5*(K3D(i,j,k,bi,bj)+K3D(i,j-1,k,bi,bj))*GM_advect
#endif
     &     )*SlopeY(i,j)*taperFct(i,j)
         ENDDO
        ENDDO
       ENDIF
#endif /* GM_EXTRA_DIAGONAL */

#ifdef ALLOW_DIAGNOSTICS
       IF (doDiagRediFlx) THEN
        km1 = MAX(k-1,1)
        DO j=1,sNy+1
         DO i=1,sNx
C         store in tmp1k Kvz_Redi
#ifdef ALLOW_KAPREDI_CONTROL
#  ifdef ALLOW_KAPREDI_CONTROL_OLD
          tmp1k(i,j) = ( kapRedi(i,j,k,bi,bj)
#  else
          tmp1k(i,j) = ( op5*(kapRedi(i,j-1,k,bi,bj)
     &                       +kapRedi(i,j,k,bi,bj))
#  endif
#else
          tmp1k(i,j) = ( GM_isopycK*GM_isoFac1d(k)
     &        *op5*(GM_isoFac2d(i,j-1,bi,bj)+GM_isoFac2d(i,j,bi,bj))
#endif
#ifdef GM_GEOM_VARIABLE_K
C     &     +(GEOMK(i,j,k,bi,bj)+GEOMK(i,j-1,k,bi,bj))*0.5 _d 0
#endif
#if ((defined GM_K3D) && ! (defined GM_K3D_PASSIVE))
     &     +op5*(K3D(i,j,k,bi,bj)+K3D(i,j-1,k,bi,bj))
#endif
     &                 )*SlopeY(i,j)*taperFct(i,j)
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx
C-        Vertical gradients interpolated to U points
          dTdz = (
     &     +recip_drC(k)*
     &       ( maskC(i,j-1,k,bi,bj)*
     &           (theta(i,j-1,km1,bi,bj)-theta(i,j-1,k,bi,bj))
     &        +maskC(i, j ,k,bi,bj)*
     &           (theta(i, j ,km1,bi,bj)-theta(i, j ,k,bi,bj))
     &       )
     &     +recip_drC(kp1)*
     &       ( maskC(i,j-1,kp1,bi,bj)*
     &           (theta(i,j-1,k,bi,bj)-theta(i,j-1,kp1,bi,bj))
     &        +maskC(i, j ,kp1,bi,bj)*
     &           (theta(i, j ,k,bi,bj)-theta(i, j ,kp1,bi,bj))
     &       )      ) * 0.25 _d 0
           tmp1k(i,j) = dxG(i,j,bi,bj)*drF(k)
     &                * _hFacS(i,j,k,bi,bj)
     &                * tmp1k(i,j) * dTdz
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(tmp1k, 'GM_KvzTz', k,1,2,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

C-- end 3rd  loop on vertical level index k
      ENDDO

#endif /* GM_NON_UNITY_DIAGONAL || GM_EXTRA_DIAGONAL */

#ifdef ALLOW_TIMEAVE
C--   Time-average
      IF ( taveFreq.GT.0. ) THEN

         CALL TIMEAVE_CUMULATE( GM_Kwx_T, Kwx, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_Kwy_T, Kwy, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_Kwz_T, Kwz, Nr,
     &                          deltaTclock, bi, bj, myThid )
#ifdef GM_GEOM_VARIABLE_K
       if(GEOM_alpha /= 0.0 _d 0) then
         CALL TIMEAVE_CUMULATE( GEOM_K_T, GEOMK, Nr,
     &                          deltaTclock, bi, bj, myThid )
       ENDIF
#endif
#ifdef GM_BOLUS_ADVEC
       IF ( GM_AdvForm ) THEN
         CALL TIMEAVE_CUMULATE( GM_PsiXtave, GM_PsiX, Nr,
     &                          deltaTclock, bi, bj, myThid )
         CALL TIMEAVE_CUMULATE( GM_PsiYtave, GM_PsiY, Nr,
     &                          deltaTclock, bi, bj, myThid )
       ENDIF
#endif
       GM_timeAve(bi,bj) = GM_timeAve(bi,bj)+deltaTclock

      ENDIF
#endif /* ALLOW_TIMEAVE */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL GMREDI_DIAGNOSTICS_FILL(bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_GMREDI */
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- JM + JRM hacks
C
#ifdef GM_GEOM_VARIABLE_K

C-- checkpointing procedure
C   set the GEOM checkpoint to be the same as something
      
      if(mod(myTime, dumpFreq) == 0. _d 0) then
        if(gm_checkpoint) then
          if(myTime == t_dum_init) then
            print *, "JRM+JM mod: skipping first check point"
          else
            print *, "JRM+JM mod: check point for energy"
           open(newunit = unit, file = "gm_checkpoint.dat",
     &        status = "replace", access = "stream", action = "write")
           write(unit) energy_old
           if(energy_local) then
             write(unit) energy_rhs
            end if
           write(unit) energy_initialised
            close(unit)
          end if
        end if
      end if
C-- diagnostics
      if(mod(myTime, monitorFreq) == 0. _d 0) then
       Utrans_tot = 0. _d 0
       Utrans_bt = 0. _d 0
       Utrans_bc = 0. _d 0
       do k = 1,Nr
        do j = 1,sNy
         do i = 1,sNx
C--   maskW the mask at U points
C--   Calculate the total transport
          Utrans_tot = Utrans_tot + uVel(i,j,k,bi,bj)
     &                  * maskW(i,j,k,bi,bj) * deepFacC(k)
     &                  * dxC(i,j,bi,bj) * dyG(i,j,bi,bj)
     &                  * drF(k) * hFacW(i,j,k,bi,bj)
C--   Calculate the "barotropic" transport 
          Utrans_bt = Utrans_bt + uVel(i,j,Nr,bi,bj)
     &                  * maskW(i,j,k,bi,bj) * deepFacC(k)
     &                  * dxC(i,j,bi,bj) * dyG(i,j,bi,bj)
     &                  * drF(k) * hFacW(i,j,k,bi,bj) 
C--   Calculate the "baroclinic" transport
          Utrans_bc = Utrans_bc
     &                  + (uVel(i,j,k,bi,bj)-uVel(i,j,Nr,bi,bj))
     &                  * maskW(i,j,k,bi,bj) * deepFacC(k)
     &                  * dxC(i,j,bi,bj) * dyG(i,j,bi,bj)
     &                  * drF(k) * hFacW(i,j,k,bi,bj)
         end do
        end do
       end do
      
C--    Zonal average
       Utrans_tot = Utrans_tot / zonal_length
       Utrans_bt = Utrans_bt / zonal_length
       Utrans_bc = Utrans_bc / zonal_length
      
       energy_avg = 0. _d 0
       if(energy_local) then
        do j = 1,sNy
         do i = 1,sNx
          if(maskC(i,j,1,bi,bj) /= 0.0 _d 0) then
           energy_avg = energy_avg + energy(i,j)
     &                * dxG(i,j,bi,bj) * dyG(i,j,bi,bj)
          end if
         end do
        end do
        energy_avg = energy_avg / domain_vol
       else
        energy_avg = energy(1,1) / domain_vol
       end if
C      print *, "JRM", energy_avg, energy0, abs(energy_avg - energy0),
C     &  100.0 _d 0 * abs(energy_avg - energy0) / abs(energy0)

C--    Print out some useful numbers
       write(*,'(a,9(3x,f11.5))'), "JRM+JM mod: = "
     &  , myTime / (3600 * 24 * 360)
     &  , Utrans_tot / 1.0 _d 6
     &  , Utrans_bt / 1.0 _d 6
     &  , Utrans_bc / 1.0 _d 6
     &  , maxval(GEOMK(1:sNx,1:sNy,:,1,1))
     &  , maxval(kappa0(1:sNx,1:sNy))
     &  , energy_avg
      end if

#ifdef ALLOW_DIAGNOSTICS
      if (useDiagnostics) then
C    0 is fill in same order as array
C    Nr is fill Nr levels, 1 is fill top level
C    2 is do the bi-bj loop outside (???)
       call DIAGNOSTICS_FILL(GEOMK, 
     &                       'GEOMkapp', 0,Nr,2,bi,bj,myThid)
       call DIAGNOSTICS_FILL(kappa0,
     &                       'GEOMkap0', 0,1,2,bi,bj,myThid)
       call DIAGNOSTICS_FILL(energy, 
     &                       'GEOMeE  ', 0,1,2,bi,bj,myThid)
       call DIAGNOSTICS_FILL(structure_func, 
     &                       'GEOMstru', 0,Nr,2,bi,bj,myThid)
       if(energy_local) then
        call DIAGNOSTICS_FILL(trd_energy_gen, 
     &                        'trdEgen ', 0,1,2,bi,bj,myThid)
C       minus sign here so "trd_adv" is term on LHS of equation
        call DIAGNOSTICS_FILL(-trd_energy_adv, 
     &                        'trdEadv ', 0,1,2,bi,bj,myThid)
        call DIAGNOSTICS_FILL(trd_energy_lap, 
     &                        'trdElap ', 0,1,2,bi,bj,myThid)
       end if
      end if
#endif /* ALLOW_DIAGNOSTICS */

#endif /* GM_GEOM_VARIABLE_K */
C
C-- end JM + JRM hacks
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C     !ROUTINE: GMREDI_CALC_TENSOR_DUMMY
C     !INTERFACE:
      SUBROUTINE GMREDI_CALC_TENSOR_DUMMY(
     I             iMin, iMax, jMin, jMax,
     I             sigmaX, sigmaY, sigmaR,
     I             bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_CALC_TENSOR_DUMMY
C     | o Calculate tensor elements for GM/Redi tensor.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GMREDI.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER iMin,iMax,jMin,jMax
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_GMREDI
C     !LOCAL VARIABLES:
      INTEGER i, j, k

      DO k=1,Nr
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
         Kwx(i,j,k,bi,bj) = 0.0
         Kwy(i,j,k,bi,bj) = 0.0
         Kwz(i,j,k,bi,bj) = 0.0
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_GMREDI */

      RETURN
      END

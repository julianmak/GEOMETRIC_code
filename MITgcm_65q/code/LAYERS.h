C $Header: /u/gcmpack/MITgcm/pkg/layers/LAYERS.h,v 1.6 2010/12/16 00:56:48 dfer Exp $
C $Name: checkpoint62x $

#ifdef ALLOW_LAYERS

C--   Header for LAYERS package. By Ryan Abernathey.
C--   For computing volume fluxes in isopyncal layers

C --  Parms
      LOGICAL layers_MNC, layers_MDSIO, useBOLUS
      INTEGER LAYER_nb, layers_kref
      _RL layers_taveFreq, layers_diagFreq
      COMMON /LAYERS_PARMS/ layers_MNC, layers_MDSIO,
     &                      layers_taveFreq, layers_diagFreq,
     &                      LAYER_nb, layers_kref, useBOLUS

C     3D Layers fields. The vertical dimension in these fields is nLayers,
C     i.e. the isopycnal coordinate.
C
C      layers_UFlux :: U integrated over layer (m^2/s)
C      layers_GMU   :: GM contribution to U integrated over layer.
C      layers_VFlux :: V integrated over layer (m^2/s)
C      layers_GMV   :: GM contribution to V integrated over layer.
C      layers_HU    :: Layer thickness at the U point (m)
C      layers_UU    :: U velocity in the layer at the U point (m/s)
C      layers_ZU    :: isopycnal layer location at U point (m)
C      layers_ZZU   :: isopycnal layer location squared at U point (m^2)
C      layers_HV    :: Layer thickness at the V point (m)
C      layers_VV    :: V velocity in the layer at the V point (m/s)
C      layers_ZV    :: isopycnal layer location at V point (m)
C      layers_ZZV   :: isopycnal layer location squared at U point (m^2)

#ifdef LAYERS_UFLUX
      _RL layers_UFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#ifdef ALLOW_GMREDI
      _RL layers_GMU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
      _RL layers_HU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_UU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZZU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* LAYERS_THICKNESS */
      COMMON /LAYERS_U/ layers_UFlux
#ifdef ALLOW_GMREDI
     & , layers_GMU
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
     &  , layers_HU, layers_UU, layers_ZU, layers_ZZU
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
      _RL layers_VFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#ifdef ALLOW_GMREDI
      _RL layers_GMV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
      _RL layers_HV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_VV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZZV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* LAYERS_THICKNESS */
      COMMON /LAYERS_V/ layers_VFlux
#ifdef ALLOW_GMREDI
     & , layers_GMV
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
     &  , layers_HV, layers_VV, layers_ZV, layers_ZZV
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
      _RL prho(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /LAYERS_PRHO/ prho
#endif

#ifdef ALLOW_TIMEAVE
C-- The same variables, time-averaged

C     Keep track of time
      _RL layers_TimeAve(nSx,nSy)
      COMMON /LAYERS_TAVE/ layers_TimeAve

#ifdef LAYERS_UFLUX
      _RL layers_UFlux_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#ifdef ALLOW_GMREDI
      _RL layers_GMU_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
      _RL layers_HU_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_UU_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZU_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZZU_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* LAYERS_THICKNESS */
      COMMON /LAYERS_U_TAVE/ layers_UFlux_T
#ifdef ALLOW_GMREDI
     &  , layers_GMU_T
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
     &  , layers_HU_T, layers_UU_T, layers_ZU_T, layers_ZZU_T
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
      _RL layers_VFlux_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#ifdef ALLOW_GMREDI
      _RL layers_GMV_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
      _RL layers_HV_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_VV_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZV_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
      _RL layers_ZZV_T(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nlayers,nSx,nSy)
#endif /* LAYERS_THICKNESS */
      COMMON /LAYERS_V_TAVE/ layers_VFlux_T
#ifdef ALLOW_GMREDI
     &  , layers_GMV_T
#endif /* ALLOW_GMREDI */
#ifdef LAYERS_THICKNESS
     &  , layers_HV_T, layers_VV_T, layers_ZV_T, layers_ZZV_T
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */

#ifdef LAYERS_PRHO_REF
      _RL prho_tave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      COMMON /LAYERS_RPHO_TAVE/ prho_tave
#endif

#endif /* ALLOW_TIMEAVE */

C     Isopycnal grid parameters:
C      layers_G :: boundaries of isopycnal layers
C      dZZf     :: height of fine grid cells
C      NZZ      :: the number of levels to use in the fine grid
C      MapIndex :: indices for mapping ZZ to Z
C      MapFact  :: factors for interpolating T(Z) to T(ZZ)

      _RL layers_G(nLayers+1)
      _RL dZZf(FineGridMax)
      INTEGER MapIndex(FineGridMax), CellIndex(FineGridMax)
      _RL MapFact(FineGridMax)
      INTEGER NZZ
      COMMON /LAYERS_VERT_GRID_I/
     &      NZZ, MapIndex, CellIndex
      COMMON /LAYERS_VERT_GRID_R/
     &      layers_G, MapFact, dZZf


#endif /* ALLOW_LAYERS */

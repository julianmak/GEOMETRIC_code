MODULE ldfdyn_oce
   !!======================================================================
   !!                  ***  MODULE  ldfdyn_oce  ***
   !! Ocean physics:  lateral momentum mixing coefficient defined in memory 
   !!======================================================================
   !! History :  1.0  ! 2002-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PUBLIC

   !                                       !!* Namelist namdyn_ldf : lateral mixing *
   LOGICAL , PUBLIC ::   ln_dynldf_lap      !: laplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_bilap    !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_level    !: iso-level direction
   LOGICAL , PUBLIC ::   ln_dynldf_hor      !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_dynldf_iso      !: iso-neutral direction
   REAL(wp), PUBLIC ::   rn_ahm_0_lap       !: lateral laplacian eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahmb_0          !: lateral laplacian background eddy viscosity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahm_0_blp       !: lateral bilaplacian eddy viscosity (m4/s)
   REAL(wp), PUBLIC ::   ahm0, ahmb0, ahm0_blp         !: OLD namelist names
   REAL(wp), PUBLIC ::   rn_cmsmag_1        !: constant in laplacian Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_cmsmag_2        !: constant in bilaplacian Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_cmsh            !: 1 or 0 , if 0 -use only shear for Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_ahm_m_blp       !: upper limit for bilap  abs(ahm) < min( dx^4/128rdt, rn_ahm_m_blp)
   REAL(wp), PUBLIC ::   rn_ahm_m_lap       !: upper limit for lap  ahm < min(dx^2/16rdt, rn_ahm_m_lap)

   INTEGER , PUBLIC ::   nkahm_smag      =  0          !: 

   !                                                                                  !!! eddy coeff. at U-,V-,W-pts [m2/s]
#if defined key_dynldf_c3d
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahm1, ahm2, ahm3, ahm4   !: ** 3D coefficients **
#elif defined key_dynldf_c2d
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ahm1, ahm2, ahm3, ahm4   !: ** 2D coefficients **
#elif defined key_dynldf_c1d
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     ::   ahm1, ahm2, ahm3, ahm4   !: ** 2D coefficients **
#else
   REAL(wp), PUBLIC                                      ::   ahm1, ahm2, ahm3, ahm4   !: ** 0D coefficients **
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ldfdyn_oce_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION ldfdyn_oce_alloc  ***
      !!----------------------------------------------------------------------
      ldfdyn_oce_alloc = 0
#if defined key_dynldf_c3d
      ALLOCATE( ahm1(jpi,jpj,jpk) , ahm2(jpi,jpj,jpk) , ahm3(jpi,jpj,jpk) , ahm4(jpi,jpj,jpk) , STAT=ldfdyn_oce_alloc )
#elif defined key_dynldf_c2d
      ALLOCATE( ahm1(jpi,jpj    ) , ahm2(jpi,jpj    ) , ahm3(jpi,jpj    ) , ahm4(jpi,jpj    ) , STAT=ldfdyn_oce_alloc )
#elif defined key_dynldf_c1d
      ALLOCATE( ahm1(        jpk) , ahm2(        jpk) , ahm3(        jpk) , ahm4(        jpk) , STAT=ldfdyn_oce_alloc )
#endif
      IF( ldfdyn_oce_alloc /= 0 )   CALL ctl_warn('ldfdyn_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION ldfdyn_oce_alloc

   !!======================================================================
END MODULE ldfdyn_oce

MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficient 
   !!=====================================================================
   !! History :        ! 1997-07  (G. Madec)  from inimix.F split in 2 routines
   !!   NEMO      1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!             2.0  ! 2005-11  (G. Madec)  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE ldfslp          ! ???
   USE in_out_manager  ! I/O manager
   USE ioipsl
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_init   ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !! 
      !! ** Purpose :   initializations of the tracer lateral mixing coeff.
      !!
      !! ** Method  :   the Eddy diffusivity and eddy induced velocity ceoff.
      !!      are defined as follows:
      !!         default option   : constant coef. aht0, aeiv0 (namelist)
      !!        'key_traldf_c1d': depth dependent coef. defined in 
      !!                            in ldf_tra_c1d routine
      !!        'key_traldf_c2d': latitude and longitude dependent coef.
      !!                            defined in ldf_tra_c2d routine
      !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
      !!                            defined in ldf_tra_c3d routine
      !!
      !!      N.B. User defined include files.  By default, 3d and 2d coef.
      !!      are set to a constant value given in the namelist and the 1d
      !!      coefficients are initialized to a hyperbolic tangent vertical
      !!      profile.
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio               ! temporary integer
      INTEGER ::   ios                  ! temporary integer
      LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
      !! 
      NAMELIST/namtra_ldf/ ln_traldf_lap  , ln_traldf_bilap,                  &
         &                 ln_traldf_level, ln_traldf_hor  , ln_traldf_iso,   &
         &                 ln_traldf_grif , ln_traldf_gdia ,                  &
         &                 ln_triad_iso   , ln_botmix_grif ,                  &
         &                 rn_aht_0       , rn_ahtb_0      , rn_aeiv_0,       &
         &                 rn_slpmax      , rn_chsmag      ,    rn_smsh,      &
         &                 rn_aht_m
      !!----------------------------------------------------------------------

      !  Define the lateral tracer physics parameters
      ! =============================================
    

      REWIND( numnam_ref )              ! Namelist namtra_ldf in reference namelist : Lateral physics on tracers
      READ  ( numnam_ref, namtra_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namtra_ldf in configuration namelist : Lateral physics on tracers
      READ  ( numnam_cfg, namtra_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_ldf in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_ldf )

      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_init : lateral tracer physics'
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtra_ldf : lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      laplacian operator            ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '      bilaplacian operator          ln_traldf_bilap = ', ln_traldf_bilap
         WRITE(numout,*) '      iso-level                     ln_traldf_level = ', ln_traldf_level
         WRITE(numout,*) '      horizontal (geopotential)     ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '      iso-neutral                   ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '      iso-neutral (Griffies)        ln_traldf_grif  = ', ln_traldf_grif
         WRITE(numout,*) '      Griffies strmfn diagnostics   ln_traldf_gdia  = ', ln_traldf_gdia
         WRITE(numout,*) '      lateral eddy diffusivity      rn_aht_0        = ', rn_aht_0
         WRITE(numout,*) '      background hor. diffusivity   rn_ahtb_0       = ', rn_ahtb_0
         WRITE(numout,*) '      eddy induced velocity coef.   rn_aeiv_0       = ', rn_aeiv_0
         WRITE(numout,*) '      maximum isoppycnal slope      rn_slpmax       = ', rn_slpmax
         WRITE(numout,*) '      pure lateral mixing in ML     ln_triad_iso    = ', ln_triad_iso
         WRITE(numout,*) '      lateral mixing on bottom      ln_botmix_grif  = ', ln_botmix_grif
         WRITE(numout,*)
      ENDIF

      !                                ! convert DOCTOR namelist names into OLD names
      aht0  = rn_aht_0
      ahtb0 = rn_ahtb_0
      aeiv0 = rn_aeiv_0

      !                                ! Parameter control

      ! ... Check consistency for type and direction :
      !           ==> will be done in traldf module

      ! ... Space variation of eddy coefficients
      ioptio = 0
#if defined key_traldf_c3d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio + 1
#endif
#if defined key_traldf_c2d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude)'
      ioptio = ioptio + 1
#endif
#if defined key_traldf_c1d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( depth )'
      ioptio = ioptio + 1
      IF( .NOT. ln_zco )   CALL ctl_stop( 'key_traldf_c1d can only be used in z-coordinate - full step' )
#endif
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'
        ELSEIF( ioptio > 1 ) THEN
           CALL ctl_stop('          use only one of the following keys:',   &
             &           ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d' )
      ENDIF

      IF( ln_traldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
         IF( aht0 > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal diffusivity coef. aht0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
         IF( aht0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop('The horizontal diffusivity coef. aht0 must be positive' )
      ENDIF


      !  Lateral eddy diffusivity and eddy induced velocity coefficients
      ! ================================================================
#if defined key_traldf_c3d
      IF ( .NOT. lk_ldfeke ) THEN
         CALL ldf_tra_c3d( ll_print )      ! aht = 3D coef. = F( longitude, latitude, depth )
      ENDIF
#elif defined key_traldf_c2d
      CALL ldf_tra_c2d( ll_print )      ! aht = 2D coef. = F( longitude, latitude )
#elif defined key_traldf_c1d
      CALL ldf_tra_c1d( ll_print )      ! aht = 1D coef. = F( depth )
#else
                                        ! Constant coefficients
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '      constant eddy diffusivity coef.   ahtu = ahtv = ahtw = aht0 = ', aht0
      IF( lk_traldf_eiv ) THEN
         IF(lwp)WRITE(numout,*) '      constant eddy induced velocity coef.   aeiu = aeiv = aeiw = aeiv0 = ', aeiv0
      
      ENDIF
#endif

#if defined key_traldf_smag && ! defined key_traldf_c3d
        CALL ctl_stop( 'key_traldf_smag can only be used with key_traldf_c3d' )
#endif
#if defined key_traldf_smag
        IF(lwp) WRITE(numout,*)' SMAGORINSKY DIFFUSION'
        IF(lwp .AND. rn_smsh < 1)  WRITE(numout,*)' only  shear is used '
        IF(lwp.and.ln_traldf_bilap) CALL ctl_stop(' SMAGORINSKY + BILAPLACIAN - UNSTABLE OR NON_CONSERVATIVE' )
#endif

      !
   END SUBROUTINE ldf_tra_init

#if defined key_traldf_c3d
#   include "ldftra_c3d.h90"
#elif defined key_traldf_c2d
#   include "ldftra_c2d.h90"
#elif defined key_traldf_c1d
#   include "ldftra_c1d.h90"
#endif

   !!======================================================================
END MODULE ldftra

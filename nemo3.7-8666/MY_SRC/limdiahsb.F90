MODULE limdiahsb
   !!======================================================================
   !!                       ***  MODULE limdia_hsb   ***
   !!  LIM-3 sea ice model :   diagnostics of ice model 
   !!======================================================================
   !! History :  3.4  ! 2012-10  (C. Rousset)  original code
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                       LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dia_hsb        : computation and output of the time evolution of keys variables
   !!   lim_dia_hsb_init   : initialization and namelist read
   !!----------------------------------------------------------------------
   USE ice             ! LIM-3: sea-ice variable
   USE dom_oce         ! ocean domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: sea-ice fields
   USE daymod          ! model calendar
   USE phycst          ! physical constant
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE iom             ! I/O manager
   USE lib_fortran     ! glob_sum
   USE limrst          ! ice restart

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_diahsb        ! routine called by ice_step.F90
   PUBLIC   lim_diahsb_init   ! routine called in sbcice_lim.F90

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   vol_loc_ini, sal_loc_ini, tem_loc_ini ! initial volume, salt and heat contents
   REAL(wp)                              ::   frc_sal, frc_voltop, frc_volbot, frc_temtop, frc_tembot  ! global forcing trends
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.4 , NEMO Consortium (2012)
   !! $Id: limdiahsb.F90 7646 2017-02-06 09:25:03Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_diahsb( kt )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE lim_diahsb  ***
      !!     
      !! ** Purpose: Compute the ice global heat content, salt content and volume conservation
      !!	
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt    ! number of iteration
      !!
      real(wp)   ::   zbg_ivol, zbg_svol, zbg_area, zbg_isal, zbg_item ,zbg_stem
      REAL(wp)   ::   z_frc_voltop, z_frc_volbot, z_frc_sal, z_frc_temtop, z_frc_tembot  
      REAL(wp)   ::   zdiff_vol, zdiff_sal, zdiff_tem  
      !!---------------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('lim_diahsb')

      ! ----------------------- !
      ! 1 -  Contents !
      ! ----------------------- !
      zbg_ivol = glob_sum( vt_i(:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )                  ! ice volume (km3)
      zbg_svol = glob_sum( vt_s(:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )                  ! snow volume (km3)
      zbg_area = glob_sum( at_i(:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-6 )                  ! area (km2)
      zbg_isal = glob_sum( SUM( smv_i(:,:,:), dim=3 ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 ) ! salt content (pss*km3)
      zbg_item = glob_sum( et_i * e1e2t(:,:) * tmask(:,:,1) * 1.e-20 )                      ! heat content (1.e20 J)
      zbg_stem = glob_sum( et_s * e1e2t(:,:) * tmask(:,:,1) * 1.e-20 )                      ! heat content (1.e20 J)
      
      ! ---------------------------!
      ! 2 - Trends due to forcing  !
      ! ---------------------------!
      z_frc_volbot = r1_rau0 * glob_sum( - ( wfx_ice(:,:) + wfx_snw(:,:) + wfx_err_sub(:,:) ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )  ! freshwater flux ice/snow-ocean 
      z_frc_voltop = r1_rau0 * glob_sum( - ( wfx_sub(:,:) + wfx_spr(:,:) ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )                     ! freshwater flux ice/snow-atm
      z_frc_sal    = r1_rau0 * glob_sum( - sfx(:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )                                            ! salt fluxes ice/snow-ocean
      z_frc_tembot =           glob_sum( hfx_out(:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-20 )                                         ! heat on top of ocean (and below ice)
      z_frc_temtop =           glob_sum( hfx_in (:,:) * e1e2t(:,:) * tmask(:,:,1) * 1.e-20 )                                         ! heat on top of ice-coean
      !
      frc_voltop  = frc_voltop  + z_frc_voltop  * rdt_ice ! km3
      frc_volbot  = frc_volbot  + z_frc_volbot  * rdt_ice ! km3
      frc_sal     = frc_sal     + z_frc_sal     * rdt_ice ! km3*pss
      frc_temtop  = frc_temtop  + z_frc_temtop  * rdt_ice ! 1.e20 J
      frc_tembot  = frc_tembot  + z_frc_tembot  * rdt_ice ! 1.e20 J
            
      ! ----------------------- !
      ! 3 -  Content variations !
      ! ----------------------- !
      zdiff_vol = r1_rau0 * glob_sum( ( rhoic * vt_i(:,:) + rhosn * vt_s(:,:) - vol_loc_ini(:,:)  &  ! freshwater trend (km3) 
         &                            ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 ) 
      zdiff_sal = r1_rau0 * glob_sum( ( rhoic * SUM( smv_i(:,:,:), dim=3 ) - sal_loc_ini(:,:)     &  ! salt content trend (km3*pss)
         &                            ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-9 )
      zdiff_tem =           glob_sum( ( et_i(:,:) + et_s(:,:) - tem_loc_ini(:,:)                  &  ! heat content trend (1.e20 J)
      !  &                            + SUM( qevap_ice * a_i_b, dim=3 ) &     !! clem: I think this line should be commented (but needs a check)
         &                            ) * e1e2t(:,:) * tmask(:,:,1) * 1.e-20 )

      ! ----------------------- !
      ! 4 -  Drifts             !
      ! ----------------------- !
      zdiff_vol = zdiff_vol - ( frc_voltop + frc_volbot )
      zdiff_sal = zdiff_sal - frc_sal
      zdiff_tem = zdiff_tem - ( frc_tembot - frc_temtop )

      ! ----------------------- !
      ! 5 - Diagnostics writing !
      ! ----------------------- !
      !
      IF( iom_use('ibgvolume') )  CALL iom_put( 'ibgvolume' , zdiff_vol        )   ! ice/snow volume  drift            (km3 equivalent ocean water)         
      IF( iom_use('ibgsaltco') )  CALL iom_put( 'ibgsaltco' , zdiff_sal        )   ! ice salt content drift            (psu*km3 equivalent ocean water)
      IF( iom_use('ibgheatco') )  CALL iom_put( 'ibgheatco' , zdiff_tem        )   ! ice/snow heat content drift       (1.e20 J)
      IF( iom_use('ibgheatfx') )  CALL iom_put( 'ibgheatfx' , zdiff_tem /      &   ! ice/snow heat flux drift          (W/m2)
         &                                                    glob_sum( e1e2t(:,:) * tmask(:,:,1) * 1.e-20 * kt*rdt ) )

      IF( iom_use('ibgfrcvoltop') )  CALL iom_put( 'ibgfrcvoltop' , frc_voltop )   ! vol  forcing ice/snw-atm          (km3 equivalent ocean water) 
      IF( iom_use('ibgfrcvolbot') )  CALL iom_put( 'ibgfrcvolbot' , frc_volbot )   ! vol  forcing ice/snw-ocean        (km3 equivalent ocean water) 
      IF( iom_use('ibgfrcsal') )     CALL iom_put( 'ibgfrcsal'    , frc_sal    )   ! sal - forcing                     (psu*km3 equivalent ocean water)   
      IF( iom_use('ibgfrctemtop') )  CALL iom_put( 'ibgfrctemtop' , frc_temtop )   ! heat on top of ice/snw/ocean      (1.e20 J)   
      IF( iom_use('ibgfrctembot') )  CALL iom_put( 'ibgfrctembot' , frc_tembot )   ! heat on top of ocean(below ice)   (1.e20 J)   
      IF( iom_use('ibgfrchfxtop') )  CALL iom_put( 'ibgfrchfxtop' , frc_temtop / & ! heat on top of ice/snw/ocean      (W/m2) 
         &                                                    glob_sum( e1e2t(:,:) * tmask(:,:,1) * 1.e-20 * kt*rdt ) )
      IF( iom_use('ibgfrchfxbot') )  CALL iom_put( 'ibgfrchfxbot' , frc_tembot / & ! heat on top of ocean(below ice)   (W/m2) 
         &                                                    glob_sum( e1e2t(:,:) * tmask(:,:,1) * 1.e-20 * kt*rdt ) )

      IF( iom_use('ibgvol_tot' ) )  CALL iom_put( 'ibgvol_tot'  , zbg_ivol     )   ! ice volume                        (km3)
      IF( iom_use('sbgvol_tot' ) )  CALL iom_put( 'sbgvol_tot'  , zbg_svol     )   ! snow volume                       (km3)
      IF( iom_use('ibgarea_tot') )  CALL iom_put( 'ibgarea_tot' , zbg_area     )   ! ice area                          (km2)
      IF( iom_use('ibgsalt_tot') )  CALL iom_put( 'ibgsalt_tot' , zbg_isal     )   ! ice salinity content              (pss*km3)
      IF( iom_use('ibgheat_tot') )  CALL iom_put( 'ibgheat_tot' , zbg_item     )   ! ice heat content                  (1.e20 J)
      IF( iom_use('sbgheat_tot') )  CALL iom_put( 'sbgheat_tot' , zbg_stem     )   ! snow heat content                 (1.e20 J)
      !
      IF( lrst_ice )   CALL lim_diahsb_rst( numit, 'WRITE' )
      !
      IF( nn_timing == 1 )   CALL timing_stop('lim_diahsb')
      !
   END SUBROUTINE lim_diahsb


   SUBROUTINE lim_diahsb_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE lim_diahsb_init  ***
      !!     
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!	
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      INTEGER            ::   ierror   ! local integer
      !!
      !!NAMELIST/namicehsb/ blabla
      !!----------------------------------------------------------------------
      !
      !!REWIND ( numnam_ice )              ! Read Namelist namicehsb 
      !!READ   ( numnam_ice, namicehsb )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'lim_diahsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !      
      ALLOCATE( vol_loc_ini(jpi,jpj), sal_loc_ini(jpi,jpj), tem_loc_ini(jpi,jpj), STAT=ierror )
      IF( ierror > 0 )  THEN
         CALL ctl_stop( 'lim_diahsb: unable to allocate vol_loc_ini' )
         RETURN
      ENDIF

      CALL lim_diahsb_rst( nstart, 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE lim_diahsb_init

   SUBROUTINE lim_diahsb_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE limdia_rst  ***
     !!                     
     !! ** Purpose :   Read or write DIA file in restart file
     !!
     !! ** Method  :   use of IOM library
     !!----------------------------------------------------------------------
     INTEGER         , INTENT(in) ::   kt     ! ice time-step
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
     INTEGER ::   id1, id2, id3, id4, id5   ! local integer
     !
     !!----------------------------------------------------------------------
     !
     IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
        IF( ln_rstart ) THEN                   !* Read the restart file
           !
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' lim_diahsb_rst read at it= ', kt,' date= ', ndastp
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           id1 = iom_varid( numrir, 'frc_voltop'   , ldstop = .FALSE. )
           id2 = iom_varid( numrir, 'frc_volbot'   , ldstop = .FALSE. )
           id3 = iom_varid( numrir, 'frc_temtop'   , ldstop = .FALSE. )
           id4 = iom_varid( numrir, 'frc_tembot'   , ldstop = .FALSE. )
           id5 = iom_varid( numrir, 'frc_sal'      , ldstop = .FALSE. )
           !
           IF( MIN( id1, id2, id3, id4, id5 ) > 0 ) THEN        ! all required arrays exist
              IF(lwp) write(numout,*) ' all files for lim_diahsb exist, loading ...'
              CALL iom_get( numrir, 'frc_voltop' , frc_voltop  )
              CALL iom_get( numrir, 'frc_volbot' , frc_volbot  )
              CALL iom_get( numrir, 'frc_temtop' , frc_temtop  )
              CALL iom_get( numrir, 'frc_tembot' , frc_tembot  )
              CALL iom_get( numrir, 'frc_sal'    , frc_sal     )
              CALL iom_get( numrir, jpdom_autoglo, 'vol_loc_ini', vol_loc_ini )
              CALL iom_get( numrir, jpdom_autoglo, 'tem_loc_ini', tem_loc_ini )
              CALL iom_get( numrir, jpdom_autoglo, 'sal_loc_ini', sal_loc_ini )
           ELSE                                                 ! one at least array is missing
              IF(lwp) write(numout,*) ' not all files for lim_diahsb exist '
              IF(lwp) write(numout,*) '    --- initialize as for initial state'
              ! set trends to 0
              frc_voltop  = 0._wp
              frc_volbot  = 0._wp
              frc_temtop  = 0._wp
              frc_tembot  = 0._wp
              frc_sal     = 0._wp
              ! record initial ice volume, salt and temp
              vol_loc_ini(:,:) = rhoic * vt_i(:,:) + rhosn * vt_s(:,:)  ! ice/snow volume (kg/m2)
              tem_loc_ini(:,:) = et_i(:,:) + et_s(:,:)                  ! ice/snow heat content (J)
              sal_loc_ini(:,:) = rhoic * SUM( smv_i(:,:,:), dim=3 )     ! ice salt content (pss*kg/m2)
           ENDIF

        ELSE
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           IF(lwp) WRITE(numout,*) ' lim_diahsb at initial state '
           IF(lwp) WRITE(numout,*) '~~~~~~~'
           ! set trends to 0
           frc_voltop  = 0._wp                                          
           frc_volbot  = 0._wp                                          
           frc_temtop  = 0._wp                                                 
           frc_tembot  = 0._wp                                                 
           frc_sal     = 0._wp                                                 
           ! record initial ice volume, salt and temp
           vol_loc_ini(:,:) = rhoic * vt_i(:,:) + rhosn * vt_s(:,:)  ! ice/snow volume (kg/m2)
           tem_loc_ini(:,:) = et_i(:,:) + et_s(:,:)                  ! ice/snow heat content (J)
           sal_loc_ini(:,:) = rhoic * SUM( smv_i(:,:,:), dim=3 )     ! ice salt content (pss*kg/m2)
           
       ENDIF

     ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
        !                                   ! -------------------
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        IF(lwp) WRITE(numout,*) ' lim_diahsb_rst write at it= ', kt,' date= ', ndastp
        IF(lwp) WRITE(numout,*) '~~~~~~~'
        CALL iom_rstput( kt, nitrst, numriw, 'frc_voltop' , frc_voltop  )
        CALL iom_rstput( kt, nitrst, numriw, 'frc_volbot' , frc_volbot  )
        CALL iom_rstput( kt, nitrst, numriw, 'frc_temtop' , frc_temtop  )
        CALL iom_rstput( kt, nitrst, numriw, 'frc_tembot' , frc_tembot  )
        CALL iom_rstput( kt, nitrst, numriw, 'frc_sal'    , frc_sal     )
        CALL iom_rstput( kt, nitrst, numriw, 'vol_loc_ini', vol_loc_ini )
        CALL iom_rstput( kt, nitrst, numriw, 'tem_loc_ini', tem_loc_ini )
        CALL iom_rstput( kt, nitrst, numriw, 'sal_loc_ini', sal_loc_ini )
        !
     ENDIF
     !
   END SUBROUTINE lim_diahsb_rst
 
#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_diahsb          ! Empty routine
   END SUBROUTINE lim_diahsb
#endif
   !!======================================================================
END MODULE limdiahsb

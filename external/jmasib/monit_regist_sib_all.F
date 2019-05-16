      SUBROUTINE MONIT_REGIST_SIB0109_ALL 
!
      USE COM_RUNCONF_SIB0109 , ONLY :
!    1   JCN_SIB_MONIT_DETAIL , 
     1   JCN_SOIL_HEAT_EXPLICIT ,
     1   JCN_SIB_ALBEDO  , 
     1   JCN_SIB_RAD2SIB ,
     1   JCN_SIB_PRESIB  , 
     1   JCN_SIB_SKIN    ,
     1   JCN_SIB_SNOW    ,
     1   JCN_SIB_SOIL 
!
      use canopy_integration_overimp , only :
     1    MR_CANOPY_INTEGRATION_OVERIMP,
     1    MR_CANOPY_INTEGRATION_SS1
!
      use canopy_main , only : mr_canopy_monit_meanvalues
!
      USE SIB_MONIT , ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
!            123456789012345678901234567890
!       CALL MONIT_REGIST_TMP_SS1
!
! ê·Ç†ÇËñ≥ÇµïΩãœìôÇÃì¡ï èoóÕílÇÃìoò^ ( ÉÇÉfÉãàÀë∂ ) 
!
        CALL MR_SIB_AREAMEAN
!
! ÇªÇÍà»äOÇÃìoò^
!
        IF ( JCN_SIB_ALBEDO  .EQ. 1 ) THEN
          CALL MR_ALBEDO_MAIN                 ! çœÇ›
        ENDIF
        IF ( JCN_SIB_RAD2SIB .EQ. 1 ) THEN
          CALL MR_RAD_SHRT_ADD                ! çœÇ›
          CALL MR_RAD_SUM_ALL                 ! çœÇ›
        ENDIF
!
        IF ( JCN_SIB_PRESIB  .EQ. 1 ) THEN
          CALL MR_PRESIB_CONDUC_ALL           ! çœÇ›
        ENDIF


        IF ( JCN_SIB_SKIN    .EQ. 1 ) THEN 
          CALL MR_CANOPY_INTEGRATION_SS1      ! çœÇ›
          CALL MR_CANOPY_INTERC               ! 9 Ç… .. é˚éxä÷åWÇ†ÇË
          CALL MR_CANOPY_STRES                ! çœÇ›
          CALL MR_CANOPY_SENS                 ! çœÇ›
          CALL MR_CANOPY_EVAP                 ! çœÇ›
          CALL MR_CANOPY_HDF                  ! çœÇ›
          CALL MR_CANOPY_TRNS_DIV             ! 9 Ç… .. é˚éxä÷åWÇ†ÇË  
          CALL MR_CANOPY_AVR_NOSSNW           ! çœÇ› 
          CALL MR_CANOPY_HEAT_FLUX 
          CALL MR_CANOPY_FLUX_EQSET
          CALL MR_CANOPY_WATER_FLUX_CORR 
          CALL MR_CANOPY_WATER_FLUX_FINAL
          CALL MR_CANOPY_WATER_DROP 
          CALL MR_CANOPY_INTEGRATION_OVERIMP
          CALL MR_CANOPY_GET_GMUV
          CALL MR_CANOPY_MONIT_MEANVALUES
        ENDIF
!
        CALL MR_MOVE_ICE_INTO_SNOW 
!
        IF ( JCN_SIB_SNOW .EQ. 1 ) THEN 
!         call COM_MONIT_SET_IFLG_SS ( 1 )
          CALL MR_SNOW_PRE
          CALL MR_SNOW_HFLUX_PRE
          CALL MR_SNOW_WFLUX_PRE
          CALL MR_SNOW_AGE_INTG
          CALL MR_SNOW_WFLUX_POST
          CALL MR_SNOW_SLEV_MOD
          CALL MR_SNOW_LITTLE_SNOW
          CALL MR_SNOW_BERG 
          CALL MR_SNOW_BOTTOM_MELT
          CALL MR_SNOW_MONIT_FIXED_DEPTH
        ENDIF
!
        CALL MR_SOIL_PRE_2
        IF ( JCN_SOIL_HEAT_EXPLICIT .EQ. 0 ) THEN
          CALL MR_SOIL_HFLUX_IMPLICIT
        ELSE
          CALL MR_SOIL_HFLUX_EXPLICIT
        ENDIF
        CALL MR_SOIL_WFLUX
        CALL MR_SOIL_FLUX_CONV
        CALL MR_SOIL_MOD
        CALL MR_SOIL_AREA_MEAN
        CALL MR_SOIL_SNOWCVR_MOD

!
        CALL MR_SNOWMOD_MAIN
        CALL MR_MONIT_MOD_MON_LAST 
!!       call COM_MONIT_SET_IFLG_SS ( 0 )
        CALL MR_MONIT_SNAP_CANOPY
!!       call COM_MONIT_SET_IFLG_SS ( 1 )
        CALL MR_MONIT_SNAP_SNOW
        CALL MR_MONIT_SNAP_SOIL
!
        CALL MR_MONIT_TOTAL_H2O_ENG_OF_SKIN
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'TEST1'
      CTITLE = 'test1'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST2'
      CTITLE = 'test2'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST3'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST4'
      CTITLE = 'test1'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST5'
      CTITLE = 'test2'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST6'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST7'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST8'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TEST9'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTA'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTB'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTC'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTD'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTE'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTF'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTG'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'TESTH'
      CTITLE = 'test3'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      RETURN
      END

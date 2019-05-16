      SUBROUTINE MR_SIB_AREAMEAN
!
      USE SIBPRM , ONLY : IDP
!
      USE SIB_MONIT , ONLY : 
     1     IMONIT_LEVEL        , 
     1     IMONIT_MEAN_MRI     ,
     1     IMONIT_FLUX_MRI     ,
     1     IMONIT_SNAP_NPD     ,
     1     IMONIT_FLUX_COMMON
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
      INTEGER        :: K
!
! 雪あり無し平均等の特別出力値のうちで、陸・海共通のモニタはここで登録
!
! ===============================
! >> 共通 : 格子平均フラックス <<
! ===============================
!

      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!  顕熱  -(HSNNF+HSNCF+HSNGF+HSNWF)
        CMARK  = 'FLSH'
        CTITLE = 'Sensible heat flux'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  潜熱  
!    キャノピーの葉の蒸発/昇華・結露/結霜、
!    下草の葉の      蒸発/昇華・結露/結霜、
!    裸地・雪表面の  蒸発/昇華・結露/結霜、
!    キャノピーの葉の蒸散
!    下草の葉の      蒸散
!       - ( WEVNF+WEVSF+WEVGF+WEVBNLD+WEVBNLU+WEVBSWD+WEVBSWU ) * 2507000
!       - ( GSLNF+GSLSF+GSLGF+GSLBNLD+GSLBNLU+GSLBSWD+GSLBSWU ) * 2834500
!       - ( WTRBNLF+WTRBSLF+WTRBGLF) * 2507000 
!
!    下記 3 要素 ( 昇華・遮断損失・蒸散 ) を引くと蒸発である。
!       - (                   WEVBNLD+WEVBNLU+WEVBSWD+WEVBSWU ) * 2507000
!
        CMARK  = 'FLLH'
        CTITLE = 'Latent heat flux'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   蒸発散水フラックス
        CMARK  = 'EVSPS'        !  登録は海と共通
        CTITLE = 'Evapo-transpiration and Sublimation'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  運動量フラックス ( canopy_get_gmuv )
!     CMARK='UMOM' = MFLUG + MFLUW    登録は海と共通
!     CMARK='VMOM' = MFLVG + MFLVW    登録は海と共通
!
        CMARK  = 'UMOM'
        CTITLE = 'Momentum flux (U)'
        CUNIT  = 'N/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'VMOM'
        CTITLE = 'Momentum flux (U)'
        CUNIT  = 'N/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! ==========================
! >> MRI : 格子平均物理量 <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
!
!   粗度長 ( オフラインのみ ... 海洋上でも出すため )
        CMARK  = 'Z0'
        CTITLE = 'Roughness Length'
        CUNIT  = 'm'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   蒸発効率  ( とりあえず枠だけ ) 
        CMARK  = 'BETA'
        CTITLE = 'Evaporation Efficiency'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_MRI ) THEN
!
!   ポテンシャル蒸発量  ( とりあえず枠だけ ) 
        CMARK  = 'EPOT'
        CTITLE = 'Potential Evaporation'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! 本庁・スナップショット
!      
      IF ( IMONIT_LEVEL .GE. IMONIT_SNAP_NPD ) THEN 
!
!   粗度長 ( オフラインのみ ... 海洋上でも出すため )
!
        CMARK  = 'Z0'
        CTITLE = 'Roughness Length'
        CUNIT  = 'm'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! 有効放射温度
!
!!      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
!!        CMARK  = 'TGEF'
!!        CTITLE = 'TEMPERATURE LONG RAD FROM SIB TOP'
!!        CUNIT  = 'K'
!!        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!!      ENDIF
!
      RETURN
      END

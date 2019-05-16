! CANOPY のフラックス諸量のモニタ。
!    evap, subl のモニタ。水・氷移動量のみ
!    corr 後に使う。
!
! 注意。trns は trns_div で行なうため、ここにはない。
!
      SUBROUTINE CANOPY_MNT_WATER_FLUX_CORR (
     I    DELT_CNP , 
     I    INSFLAG  , 
     I    EVAP     , SUBL     , TRNS )
!
! 陰的補正直後の水関連フラックスのモニター。
! いまのところ不要では?
!
      USE SIBPRM , ONLY : IDX , ICG , ICGB 
!
      USE SIBCON , ONLY :
     1   HL_WTR_c , 
     1   HL_W2G_c 
!
      USE SIB_MONIT , ONLY : MONIT_SIB, ICN, ICS,
     1                       IMONIT_LEVEL , IMONIT_ALL
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: INSFLAG 
      REAL(8),INTENT(IN)   :: DELT_CNP
!
      REAL(8),INTENT(IN)   :: EVAP (IDX,ICGB) 
      REAL(8),INTENT(IN)   :: SUBL (IDX,ICGB) 
      REAL(8),INTENT(IN)   :: TRNS (IDX,ICG ) 
!
      CHARACTER(1)      :: CFLAG_NS_C1 
      CHARACTER(1)      :: CFLAG_NS_G1
!
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
      INTEGER           :: IFLG_MON
!
      IF ( INSFLAG .EQ.0 ) THEN
        CFLAG_NS_C1 ='N'
        CFLAG_NS_G1 ='G'
        IFLG_MON    = ICN
      ELSE
        CFLAG_NS_C1 ='S'
        CFLAG_NS_G1 ='W'
        IFLG_MON    = ICS
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ==========
! >> 蒸発 <<
! ==========
!
! キャノピー
!
      CMARK  = 'WEVxC'
      CMARK(4:4) = CFLAG_NS_C1 
      CALL MONIT_SIB ( CMARK, EVAP(1,1) , IFLG_MON, DELT_CNP, ONE)   
!
! 下草
!
      IF ( INSFLAG .EQ. 0 ) THEN
      CMARK  = 'WEVGC'
      CALL MONIT_SIB ( CMARK, EVAP(1,2) , IFLG_MON, DELT_CNP, ONE)   
      ENDIF
!
! 裸地・雪
!
      CMARK  = 'WEVxC'
      CMARK(4:4) = CFLAG_NS_G1 
      CALL MONIT_SIB ( CMARK, EVAP(1,3) , IFLG_MON, DELT_CNP, ONE)   
!
! ==========
! >> 昇華 <<
! ==========
!
! キャノピー
!
      CMARK  = 'GSLxC'
      CMARK(4:4) = CFLAG_NS_C1 
      CALL MONIT_SIB ( CMARK, SUBL(1,1) , IFLG_MON, DELT_CNP, ONE)   
!
! 下草
!
      IF ( INSFLAG .EQ. 0 ) THEN
      CMARK  = 'GSLGC'
      CALL MONIT_SIB ( CMARK, SUBL(1,2) , IFLG_MON, DELT_CNP, ONE)   
      ENDIF
!
! 裸地・雪
!
      CMARK  = 'GSLxC'
      CMARK(4:4) = CFLAG_NS_G1 
      CALL MONIT_SIB ( CMARK, SUBL(1,3) , IFLG_MON , DELT_CNP, ONE)   
!
! ==========
! >> 蒸散 <<
! ==========
!
! キャノピー
!
      CMARK  = 'WTRBxLF'
      CMARK(5:5) = CFLAG_NS_C1 
      CALL MONIT_SIB ( CMARK, TRNS(1,1) , IFLG_MON, DELT_CNP, ONE)   
!
! 下草
!
      IF ( INSFLAG .EQ. 0 ) THEN
      CMARK  = 'WTRBGLF'
      CALL MONIT_SIB ( CMARK, TRNS(1,2) , IFLG_MON, DELT_CNP, ONE)   
      ENDIF
!
      ENDIF
!
      RETURN
      END 
!
!=========================================================
      SUBROUTINE MR_CANOPY_WATER_FLUX_CORR  
!    
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL
 
!
! デバッグ用変数の登録作業
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_CANOPY_MONITOR_WATER_FLUX_CORR' , 
     1             ' 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! --------
! > 蒸発 < 
! --------
!
!  キャノピー
!
      CMARK  = 'WEVNC'
      CTITLE = 'EVAP WTR CNP AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVSC'
      CTITLE = 'EVAP WTR CNP AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'WEVGC'
      CTITLE = 'EVAP WTR GRS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   裸地・雪
!
      CMARK  = 'WEVLC'
      CTITLE = 'EVAP WTR BS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVWC'
      CTITLE = 'EVAP WTR SNW CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! --------
! > 昇華 < 
! --------
!
!  キャノピー
!
      CMARK  = 'GSLNC'
      CTITLE = 'SUBL GLA CNP AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLSC'
      CTITLE = 'SUBL GLA CNP AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'GSLGC'
      CTITLE = 'SUBL GLA GRS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   裸地・雪
!
      CMARK  = 'GSLLC'
      CTITLE = 'SUBL GLA BS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLWC'
      CTITLE = 'SUBL GLA SNW AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! --------
! > 蒸散 < 
! --------
!
!  キャノピー
!
      CMARK  = 'WTRBNLF'
      CTITLE = 'TRANS WTR CNP NOS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBSLF'
      CTITLE = 'TRANS WTR CNP SNW AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'WTRBGLF'
      CTITLE = 'TRANS WTR GRS AFTER CORRECT'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!    
      ENDIF
!
      RETURN
      END

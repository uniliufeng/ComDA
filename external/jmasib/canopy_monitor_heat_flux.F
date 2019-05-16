! CANOPY の熱フラックス諸量のモニタ。
!  (1) canopy_monitor_heat_flux  
!      sens, hdf , rad のモニタ。
!      corr 後と、最終状態の共用
!
! 注意。trns は trns_div で行なうため、ここにはない。
!
! =====================================================================
!
      SUBROUTINE CANOPY_MNT_HEAT_FLUX (
     I    DELT_CNP , 
     I    NSTR     , NEND      ,  
     I    INSFLAG  , IFLAG_CF  , 
     I    SENS     , HDF       , 
     I    RAD_LONG , 
     I    RAD_LONG_2_ATM , 
     O    FLSH     )
!
      USE SIBPRM , ONLY : IDX , ICG 
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN, ICS,
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_FLUX_COMMON ,
     1                      IMONIT_BUDJET ,
     1                      IMONIT_FLUX_MRI
!     1                      MONIT_SIB_6HR
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: INSFLAG 
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: IFLAG_CF        ! 補正後 0 、最終 1  
      REAL(8),INTENT(IN)   :: DELT_CNP
!
      REAL(8),INTENT(IN)   :: SENS (IDX,ICG) 
      REAL(8),INTENT(IN)   :: RAD_LONG  (IDX,ICG) 
      REAL(8),INTENT(IN)   :: RAD_LONG_2_ATM (IDX) 
      REAL(8),INTENT(IN)   :: HDF  (IDX    ) 
      REAL(8),INTENT(OUT)  :: FLSH (IDX    ) 
!
      CHARACTER(1)      :: CFLAG_CF
      CHARACTER(1)      :: CFLAG_NS_C1
      CHARACTER(1)      :: CFLAG_NS_G1
      CHARACTER(2)      :: CFLAG_NS_G2
!
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
!
      INTEGER :: IFLAG_MON     
      INTEGER :: I
!
      DO I=NSTR,NEND
        FLSH(I) = -(SENS(I,1)+SENS(I,2))
      ENDDO
!
#ifdef MONYOS

      IF ( IFLAG_CF .EQ.0 ) THEN
        CFLAG_CF ='C'            ! corr 
      ELSE
        CFLAG_CF ='F'            ! final 
      ENDIF
!
      IF ( INSFLAG .EQ.0 ) THEN
        CFLAG_NS_C1 ='N'
        CFLAG_NS_G1 ='G'
        CFLAG_NS_G2 ='NL'
        IFLAG_MON   = ICN
      ELSE
        CFLAG_NS_C1 ='S'
        CFLAG_NS_G1 ='W'
        CFLAG_NS_G2 ='SW'
        IFLAG_MON   = ICS
      ENDIF
!
! ==================
! >> 全格子平均値 << 
! ==================
!
      IF (       
     1           IMONIT_LEVEL .GE. IMONIT_BUDJET
     1     .AND. IFLAG_CF .EQ. 1                          ! Finel だけ
     1     ) THEN 
!
!  上向き長波放射 ( 下向き正 )
!
        CMARK  = 'ULWBL'
        CALL MONIT_SIB ( CMARK, RAD_LONG_2_ATM, IFLAG_MON, 
     1                   DELT_CNP, ONE ) 
      ENDIF
!
! ============
! >> 領域別 << 
! ============
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! ==========
! >> 放射 <<
! ==========
!
!   キャノピー
!
      CMARK  = 'HLR'
      CMARK(4:4) = CFLAG_NS_C1
      CMARK(5:5) = CFLAG_CF 
      CALL MONIT_SIB ( CMARK , RAD_LONG(1,1), IFLAG_MON, DELT_CNP , ONE)
!
!   地面
!
      CMARK  = 'HLR'
      CMARK(4:4) = CFLAG_NS_G1
      CMARK(5:5) = CFLAG_CF 
      CALL MONIT_SIB ( CMARK , RAD_LONG(1,2), IFLAG_MON, DELT_CNP , ONE)
!      
! ==========
! >> 顕熱 <<
! ==========
!
!   キャノピー
!
      CMARK  = 'HSN'
      CMARK(4:4) = CFLAG_NS_C1 
      CMARK(5:5) = CFLAG_CF 
      CALL MONIT_SIB ( CMARK , SENS(1,1), IFLAG_MON, DELT_CNP , ONE )
!
!   地面、雪
!
      CMARK  = 'HSN'
      CMARK(4:4) = CFLAG_NS_G1
      CMARK(5:5) = CFLAG_CF 
      CALL MONIT_SIB ( CMARK , SENS(1,2), IFLAG_MON, DELT_CNP , ONE )
!
! ==========
! >> 地熱 <<
! ==========
!
      CMARK  = 'HDFB'
      CMARK(5:6) = CFLAG_NS_G2
      CMARK(7:7) = CFLAG_CF 
      CALL MONIT_SIB ( CMARK , HDF  , IFLAG_MON, DELT_CNP , ONE )
!
      ENDIF
!
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_CANOPY_HEAT_FLUX 
!
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_FLUX_COMMON  ,
     1                      IMONIT_BUDJET       ,
     1                      IMONIT_FLUX_MRI
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_MONITOR_HEAT_FLUX 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!     IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!  顕熱  -(HSNNF+HSNCF+HSNGF+HSNWF)
!       CMARK  = 'FLSH'  登録は海と共通
!     ENDIF
!
! ==================
! >> 熱フラックス <<
! ==================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_BUDJET ) THEN
        CMARK  = 'ULWBL'
        CTITLE = 'Upward Long wave rad flux grid mean'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! --------
! > 顕熱 <
! --------
!
!  CORR 
!
      CMARK  = 'HSNNC'
      CTITLE = 'SENS CANOPY NOS AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNSC'
      CTITLE = 'SENS CANOPY SNW AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNGC'
      CTITLE = 'SENS GROUND AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNWC'
      CTITLE = 'SENS SNOW AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  FINAL 
!
      CMARK  = 'HSNNF'
      CTITLE = 'SENS CANOPY NOS FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNSF'
      CTITLE = 'SENS CANOPY SNW FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNGF'
      CTITLE = 'SENS GROUND FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSNWF'
      CTITLE = 'SENS SNOW FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! --------
! > 放射 <
! --------
!
!  CORR 
!
      CMARK  = 'HLRNC'
      CTITLE = 'LONG RAD CANOPY NOS AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRSC'
      CTITLE = 'LONG RAD CANOPY SNW AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRGC'
      CTITLE = 'LONG RAD GROUND AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRWC'
      CTITLE = 'LONG RAD SNOW AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  FINAL 
!
      CMARK  = 'HLRNF'
      CTITLE = 'LONG RAD CANOPY NOS FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRSF'
      CTITLE = 'LONG RAD CANOPY SNW FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRGF'
      CTITLE = 'LONG RAD GROUND FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRWF'
      CTITLE = 'LONG RAD SNOW FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ----------
! > 熱伝導 <
! ----------
!
!  CORR 
!
      CMARK  = 'HDFBNLC'
      CTITLE = 'HEAT DIFFUSE SK-SL1 AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HDFBSWC'
      CTITLE = 'HEAT DIFFUSE SK-SW1 AFTER CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  FINAL 
!
      CMARK  = 'HDFBNLF'
      CTITLE = 'HEAT DIFFUSE SK-SL1 FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HDFBSWF'
      CTITLE = 'HEAT DIFFUSE SK-SW1 FINAL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      RETURN
      END

! CANOPY のフラックス諸量のモニタ。
!  (3) canopy_monitor_water_flux_final
!      evap, subl のモニタ。
!      h2o 移動に伴う熱の移動も。
!
! 注意。蒸散は、格子平均値処理を除いては trns_div で行なうことに注意
!
      SUBROUTINE CANOPY_MNT_WATER_FLUX_FINAL (
     I    NSTR     , NEND      , DELT_CNP , 
     I    INSFLAG  , 
     I    EVAP     , SUBL      , TRNS     , 
     I    TMP_CNP  , TMP_GRS   , TMP_GRD  , 
     O    FLLH , qflux )
!
! ・蒸発・結露/昇華・結霜に伴う H2Oの移動 と それに伴う熱 の最終結果
! ・蒸散は canbopy_trns_div で行なうので、ここでは関係なし。
! ・移動する水の温度は現在次のように考えている。
!     -- 葉  からの蒸発・昇華 = 前ステップでの葉の温度 
!     -- 葉  へのの結露・結霜 = 前ステップでの葉の温度 
!     -- 下面からの蒸発・昇華 = 前ステップでの下面第一層の温度 
!     -- 下面へのの結露・結霜 = 前ステップでの下面第一層の温度 
!
      USE SIBPRM , ONLY : IDX , ICG , ICGB 
!
      USE SIBCON , ONLY :
     1   HL_WTR_c , 
     1   HL_GLA_c , 
     1   HL_W2G_c , 
     1   TMP_FREZ_C , 
     1   CAP_WTR_C , 
     1   CAP_GLA_C  
!
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN, ICS,
     1                      IMONIT_LEVEL, IMONIT_ALL ,
     1                      IMONIT_FLUX_COMMON ,
     1                      IMONIT_FLUX_MRI, IMONIT_ALMA
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: INSFLAG 
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      REAL(8),INTENT(IN)   :: DELT_CNP
!
      REAL(8),INTENT(IN)   :: EVAP     (IDX,ICGB) 
      REAL(8),INTENT(IN)   :: SUBL     (IDX,ICGB) 
      REAL(8),INTENT(IN)   :: TRNS     (IDX,ICG ) 
!
      REAL(8),INTENT(IN)   :: TMP_CNP  (IDX) 
      REAL(8),INTENT(IN)   :: TMP_GRS  (IDX) 
      REAL(8),INTENT(IN)   :: TMP_GRD  (IDX) 
!
      REAL(8),INTENT(OUT)  :: FLLH     (IDX)
      REAL(8),INTENT(OUT)  :: QFLUX   (IDX)
      REAL(8)              :: FTRS     (IDX)
      REAL(8)              :: FSBL     (IDX)
      REAL(8)              :: FINT     (IDX)
!
      CHARACTER(1)      :: CFLAG_NS_C1 
!     CHARACTER(1)      :: CFLAG_NS_G1
      CHARACTER(2)      :: CFLAG_NS_G2
!
      REAL(8)           :: WORK  (IDX)
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE  = 1.D0 
      REAL(8),PARAMETER :: ZERO = 0.D0 
      INTEGER           :: I
      INTEGER           :: IFLG_MON
!
        FLLH(:)=0. 
        FTRS(:)=0. 
        FINT(:)=0. 
        FSBL(:)=0. 
        QFLUX(:)=0. 
!
        DO I=NSTR,NEND
          FLLH(I) = - EVAP(I,3) * HL_WTR_c 
          FTRS(I) = - TRNS(I,1) * HL_WTR_c - TRNS(I,2) * HL_WTR_c 
          FINT(I) = - ( EVAP(I,1) + EVAP(I,2) ) * HL_WTR_c 
     1              - ( SUBL(I,1) + SUBL(I,2) ) * HL_GLA_c  
          FSBL(I) = - SUBL(I,3) * HL_GLA_c  
!
          FLLH(I) = FLLH(I) + FTRS(I) + FINT(I) + FSBL(I) 

!
          QFLUX(I) = - (   EVAP(I,1) + EVAP(I,2) + EVAP(I,3) 
     1                   + SUBL(I,1) + SUBL(I,2) + SUBL(I,3) 
     1                   + TRNS(I,1) + TRNS(I,2) )
        ENDDO
!
#ifdef MONYOS
!
      IF ( INSFLAG .EQ.0 ) THEN
        CFLAG_NS_C1 ='N'
        CFLAG_NS_G2 ='NL'
!       CFLAG_NS_G1 ='L'
        IFLG_MON    = ICN
      ELSE
        CFLAG_NS_C1 ='S'
        CFLAG_NS_G2 ='SW'
!       CFLAG_NS_G1 ='W'
        IFLG_MON    = ICS
      ENDIF
!
! ================
! >> 格子平均量 <<
! ================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!
!!        CMARK = 'FLLH'
!!        CALL MONIT_SIB ( CMARK, FLLH , IFLG_MON, DELT_CNP, ONE)   
!
        CMARK = 'LTRS'
        CALL MONIT_SIB ( CMARK, FTRS , IFLG_MON, DELT_CNP, ONE)   
!
        CMARK = 'LINT'
        CALL MONIT_SIB ( CMARK, FINT , IFLG_MON, DELT_CNP, ONE)   
!
        CMARK = 'LSBL'
        CALL MONIT_SIB ( CMARK, FSBL , IFLG_MON, DELT_CNP, ONE)   
!
      ENDIF
!
!!      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_MRI ) THEN
!!        CMARK = 'EVSPS'
!!        CALL MONIT_SIB ( CMARK, EVSPS, IFLG_MON, DELT_CNP, ONE)   
!!      ENDIF

!
! ==========================
! >> 個々のフラックス出力 <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ==========
! >> 蒸発 <<
! ==========
!
! -----------------------
! > 蒸発/結露した水の量 <
! -----------------------
!
!    キャノピー <---> キャノピー空気
!
        CMARK      = 'WEVxF'
        CMARK(4:4) = CFLAG_NS_C1 
        CALL MONIT_SIB ( CMARK, EVAP(1,1), IFLG_MON, DELT_CNP, ONE)   
!
!    下草  <---> キャノピー空気
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し
          CMARK  = 'WEVGF'
          CALL MONIT_SIB ( CMARK, EVAP(1,2), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
!    裸地 --> キャノピー空気 
! 
        DO I=NSTR,NEND
          WORK(I) = MIN(EVAP(I,3),ZERO)  
        ENDDO 
        CMARK      = 'WEVBxxU'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON , DELT_CNP, ONE)   
!
!    キャノピー空気 --> 裸地・雪 
!
        DO I=NSTR,NEND
          WORK(I) = MAX(EVAP(I,3),ZERO)  
        ENDDO 
        CMARK      = 'WEVBxxD'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
! -----------------------------------
! > 蒸発/結露した水に伴うエネルギー <
! -----------------------------------
!   蒸発/結露/昇華/結霜した水蒸気は、簡単のため、
!   キャノピー等の温度の水蒸気として出入りすると考える。
!   ここでは、水・氷でいる/たときのエネルギーをモニタする。(2002/08/03 変更)

!    キャノピー <---> キャノピー空気 
!
        DO I=NSTR,NEND
! 2002/08/03 保坂変更
!         WORK(I) =   ( CAP_WTR_c * ( TMP_CNP(I) - TMP_FREZ_c ) )
!    1                  + HL_WTR_c )  * EVAP(I,1)               
          WORK(I) =   CAP_WTR_c * ( TMP_CNP(I) - TMP_FREZ_c ) 
     1                          * EVAP(I,1)  
        ENDDO 
!
        CMARK      = 'IEVx'
        CMARK(4:4) = CFLAG_NS_C1 
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
!    下草  <---> キャノピー空気
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し
          DO I=NSTR,NEND
! 2002/08/03 保坂変更
!           WORK(I) =  ( CAP_WTR_c * ( TMP_GRS(I) - TMP_FREZ_c )
!    1                   + HL_WTR_c ) * EVAP(I,2)  
            WORK(I) =     CAP_WTR_c * ( TMP_GRS(I) - TMP_FREZ_c )
     1                              * EVAP(I,2)  
          ENDDO 
!
          CMARK  = 'IEVG'
          CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
!     裸地 --> キャノピー空気 
!
        DO I=NSTR,NEND
          WORK(I) = 
     1              CAP_WTR_c * ( TMP_GRD(I) - TMP_FREZ_c )
     1                * MIN(EVAP(I,3),ZERO)                    ! 結露
        ENDDO 
!
        CMARK  = 'IEVBxxU'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
!     キャノピー空気 --> 裸地・雪 
!
        DO I=NSTR,NEND
          WORK(I) = 
     1              CAP_WTR_c * ( TMP_GRD(I) - TMP_FREZ_c )
     1                * MAX(EVAP(I,3),ZERO)                    ! 結露
        ENDDO 
!
        CMARK  = 'IEVBxxD'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
! ==========
! >> 昇華 <<
! ==========
!
! -----------------------
! > 昇華/結霜した水の量 <
! -----------------------
!
!    キャノピー <---> キャノピー空気
!
        CMARK      = 'GSLxF'
        CMARK(4:4) = CFLAG_NS_C1 
        CALL MONIT_SIB ( CMARK, SUBL(1,1), IFLG_MON, DELT_CNP, ONE)   
!
!    下草  <---> キャノピー空気
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し
          CMARK  = 'GSLGF'
          CALL MONIT_SIB ( CMARK, SUBL(1,2), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
!    裸地 --> キャノピー空気 
! 
        DO I=NSTR,NEND
          WORK(I) = MIN(SUBL(I,3),ZERO)  
        ENDDO 
        CMARK      = 'GSLBxxU'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
!    キャノピー空気 --> 裸地・雪 
!
        DO I=NSTR,NEND
          WORK(I) = MAX(SUBL(I,3),ZERO)  
        ENDDO 
        CMARK  = 'GSLBxxD'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
! -----------------------------------
! > 昇華/結霜した水に伴うエネルギー <
! -----------------------------------
!
!    キャノピー <---> キャノピー空気
!
        DO I=NSTR,NEND
! 2002/08/03 保坂変更
!         WORK(I) =   ( CAP_GLA_c * ( TMP_CNP(I) - TMP_FREZ_c )
!    1                  + HL_WTR_c ) * SUBL(I,1)  
          WORK(I) =   ( CAP_GLA_c * ( TMP_CNP(I) - TMP_FREZ_c )
     1                  - HL_W2G_c ) * SUBL(I,1)  
        ENDDO 
!
        CMARK      = 'JSLx'
        CMARK(4:4) = CFLAG_NS_C1 
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
!    下草  <---> キャノピー空気
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し
          DO I=NSTR,NEND
! 2002/08/03 保坂変更
!            WORK(I) = ( CAP_GLA_c * ( TMP_GRS(I) - TMP_FREZ_c )
!     1                  + HL_WTR_c ) * SUBL(I,2)  
            WORK(I) = ( CAP_GLA_c * ( TMP_GRS(I) - TMP_FREZ_c )
     1                  - HL_W2G_c ) * SUBL(I,2)  
          ENDDO 
!
          CMARK  = 'JSLG'
          CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
!     裸地 --> キャノピー空気 
!
        DO I=NSTR,NEND
          WORK(I) = ( CAP_GLA_c * ( TMP_GRD(I) - TMP_FREZ_c )
     1                - HL_W2G_c ) 
     1              * MIN(SUBL(I,3),ZERO)                    ! 結露
        ENDDO 
!
        CMARK  = 'JSLBxxU'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
!
!     キャノピー空気 --> 裸地・雪 
!
        DO I=NSTR,NEND
          WORK(I) = ( CAP_GLA_c * ( TMP_GRD(I) - TMP_FREZ_c )
     1                - HL_W2G_c ) 
     1                * MAX(SUBL(I,3),ZERO)                    ! 結露
        ENDDO 
!
        CMARK  = 'JSLBxxD'
        CMARK(5:6) = CFLAG_NS_G2 
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, DELT_CNP, ONE)   
      ENDIF
#else
        DO I=NSTR,NEND
          FLLH(I) = - EVAP(I,3) * HL_WTR_c 
          FTRS(I) = - TRNS(I,1) * HL_WTR_c - TRNS(I,2) * HL_WTR_c 
          FINT(I) = - ( EVAP(I,1) + EVAP(I,2) ) * HL_WTR_c 
     1              - ( SUBL(I,1) + SUBL(I,2) ) * HL_GLA_c  
          FSBL(I) = - SUBL(I,3) * HL_GLA_c  
!
          FLLH(I) = FLLH(I) + FTRS(I) + FINT(I) + FSBL(I) 
        ENDDO
!
#endif
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
! Ecanop ( キャノピー・下草 ⇔ 大気 ) の蒸発・結露
!
        CMARK = 'A_EVVEG' 
        CALL MONIT_SIB ( CMARK, EVAP(1,1), IFLG_MON, DELT_CNP, ONE)   
        CALL MONIT_SIB ( CMARK, EVAP(1,2), IFLG_MON, DELT_CNP, ONE)   
!      
! Esoil ( 裸地 ⇔ 大気 ) の蒸発・結露
!
        CMARK = 'A_EVBS' 
        IF ( INSFLAG .EQ. 0 ) THEN
          CALL MONIT_SIB ( CMARK, EVAP(1,3), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
! EvapSnow ( 積雪 ⇔ 大気 ) の蒸発・結露 ( liquid water ) 
!
        CMARK = 'A_EVSS' 
        IF ( INSFLAG .EQ. 1 ) THEN
          CALL MONIT_SIB ( CMARK, EVAP(1,3), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
! SubSnow ( 積雪 ⇔ 大気 ) の昇華・結霜 ( solid water ) 
!
        CMARK = 'A_SLSS' 
        IF ( INSFLAG .EQ. 1 ) THEN
          CALL MONIT_SIB ( CMARK, SUBL(1,3), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
! 
! SubSurf ( 裸地、キャノピー・下草 ⇔ 大気 ) の昇華・結霜 ( solid water ) 
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK = 'A_SLBS' 
          CALL MONIT_SIB ( CMARK, SUBL(1,3), IFLG_MON, DELT_CNP, ONE)   
        ENDIF
!
        CMARK = 'A_SLVEG' 
        CALL MONIT_SIB ( CMARK, SUBL(1,1), IFLG_MON, DELT_CNP, ONE)   
        CALL MONIT_SIB ( CMARK, SUBL(1,2), IFLG_MON, DELT_CNP, ONE)   
!
      ENDIF
#endif
!
      RETURN
      END
!==============================================================
      SUBROUTINE MR_CANOPY_WATER_FLUX_FINAL
!
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL, IMONIT_ALL ,
     1                      IMONIT_FLUX_COMMON  , IMONIT_ALMA
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_MONITOR_WATER_FLUX_FINAL' , 
     1             ' 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
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
!       CMARK = 'FLLH'    登録は海と共通
!
!  昇華 
!       - (                   GSLBNLD+GSLBNLU+GSLBSWD+GSLBSWU ) * 2834500
!
        CMARK  = 'LSBL'
        CTITLE = 'Sublimination'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  遮断損失
!       - ( WEVNF+WEVSF+WEVGF ) * 2507000
!       - ( GSLNF+GSLSF+GSLGF ) * 2834500
!
        CMARK  = 'LINT'
        CTITLE = 'Interception'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  蒸散 
!       - ( WTRBNLF+WTRBSLF+WTRBGLF) * 2507000 
!
        CMARK  = 'LTRS'
        CTITLE = 'Transpiration'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   蒸発散水フラックス
!       CMARK  = 'EVSPS'   登録は海と共通
!
      ENDIF

!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! ================
! >> 水・氷の量 <<
! ================
!
! --------
! > 蒸発 < 
! --------
!
!  キャノピー
!
      CMARK  = 'WEVNF'
      CTITLE = 'EVAP WTR CNP NOS FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVSF'
      CTITLE = 'EVAP WTR CNP SNW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'WEVGF'
      CTITLE = 'EVAP WTR GRS FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   キャノピー空気 -> 下面
!
      CMARK  = 'WEVBNLD'
      CTITLE = 'DOWNWARD EVAP WATER INTO BARE SOIL FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVBSWD'
      CTITLE = 'DOWNWARD EVAP WATER INTO SNOW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下面 -> キャノピー空気 
!
      CMARK  = 'WEVBNLU'
      CTITLE = 'UPWARD EVAP WATER FROM BARE SOIL FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVBSWU'
      CTITLE = 'UPWARD EVAP WATER FROM SNOW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! --------
! > 昇華 < 
! --------
!
!  キャノピー
!
      CMARK  = 'GSLNF'
      CTITLE = 'SUBL GLA CNP NOS FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLSF'
      CTITLE = 'SUBL GLA CNP SNW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'GSLGF'
      CTITLE = 'SUBL GLA GRS FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   キャノピー空気 -> 下面
!
      CMARK  = 'GSLBNLD'
      CTITLE = 'DOWNWARD SUBL GLA INTO BARE SOIL ( SNW ) FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLBSWD'
      CTITLE = 'DOWNWARD SUBL GLA INTO SNOW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下面 -> キャノピー空気 
!
      CMARK  = 'GSLBNLU'
      CTITLE = 'UPWARD SUBL GLA FROM BARE SOIL FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLBSWU'
      CTITLE = 'UPWARD SUBL GLA FROM SNOW FINAL'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ============================
! >> 水・氷に伴うエネルギー <<
! ============================
!
! --------
! > 蒸発 < 
! --------
!
!  キャノピー
!
      CMARK  = 'IEVN'
      CTITLE = 'HEAT OF EVAP WATER IN CANOPY NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IEVS'
      CTITLE = 'HEAT OF EVAP WATER IN CANOPY SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草
!
      CMARK  = 'IEVG'
      CTITLE = 'HEAT OF EVAP WATER IN GRASS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   キャノピー空気 -> 下面
!
      CMARK  = 'IEVBNLD'
      CTITLE = 'HEAT OF DOWNWARD EVAP WATER INTO BARE SOIL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IEVBSWD'
      CTITLE = 'HEAT OF DOWNWARD EVAP WATER INTO SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下面 -> キャノピー空気  
!
      CMARK  = 'IEVBNLU'
      CTITLE = 'HEAT OF UPWARD EVAP WATER FROM BARE SOIL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IEVBSWU'
      CTITLE = 'HEAT OF UPWARD EVAP WATER FROM SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! -------------
! > 昇華/結霜 < 
! -------------
!
! キャノピー
!
      CMARK  = 'JSLN'
      CTITLE = 'HEAT OF SUBL WATER IN CANOPY NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JSLS'
      CTITLE = 'HEAT OF SUBL WATER IN CANOPY SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
      CMARK  = 'JSLG'
      CTITLE = 'HEAT OF SUBL WATER IN GRASS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! キャノピー空気 -> 下面 
!
      CMARK  = 'JSLBNLD'
      CTITLE = 'HEAT OF DOWNWARD SUBL WATER INTO BARE SOIL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JSLBSWD'
      CTITLE = 'HEAT OF DOWNWARD SUBL WATER INTO SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下面  -> キャノピー空気
!
      CMARK  = 'JSLBNLU'
      CTITLE = 'HEAT OF UPWARD SUBL WATER FROM BARE SOIL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JSLBSWU'
      CTITLE = 'HEAT OF UPWARD SUBL WATER FROM SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF 
!
! For ALMA 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
! Ecanop ( キャノピー・下草 ⇔ 大気 ) の蒸発・結露
!
        CMARK  = 'A_EVVEG' 
        CTITLE = 'ALMA - Ecanop evaporation from canopy and grass '
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
!      
! Esoil ( 裸地 ⇔ 大気 ) の蒸発・結露
!
        CMARK = 'A_EVBS' 
        CTITLE = 'ALMA - Esoil evaporation from bare soil'
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
!
! EvapSnow ( 積雪 ⇔ 大気 ) の蒸発・結露 ( liquid water ) 
!
        CMARK = 'A_EVSS' 
        CTITLE = 'ALMA - EvapSnow evaporation from snow'
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
!
! SubSnow ( 積雪 ⇔ 大気 ) の昇華・結霜 ( solid water ) 
!
        CMARK = 'A_SLSS' 
        CTITLE = 'ALMA - SubSnow sublimination from snow'
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
! 
! SubSurf ( 裸地、キャノピー・下草 ⇔ 大気 ) の昇華・結霜 ( solid water ) 
!
        CMARK = 'A_SLBS' 
        CTITLE = 'ALMA - SubSurf sublimination from baresoil'
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
!
        CMARK = 'A_SLVEG' 
        CTITLE = 'ALMA - SubSurf sublimination from canopy and grass'
        CUNIT  = 'kg/m**2/s'
        CALL MONIT_regist_SIB ( CMARK, CTITLE , CUNIT ) 
!
      ENDIF
!
      RETURN
      END

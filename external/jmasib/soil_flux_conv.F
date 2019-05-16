!  土壌中のフラックスの確定
!  99/07/07 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_FLUX_CONV (
     I   NSTR          , NEND          , DELT_SOIL     , 
#ifdef MONYOS
     I   INSFLAG       ,
#endif
     I   HEAT_CONV     , 
     I   WFLUX_H       , WTR_INC       , 
     I   CAP_SOIL_INV  , 
     I   TMP_SOIL      , 
     M   ENG_SOIL      , H2O_SOIL      , TMP_SOIL4  , 
     O   RUNOFF_B     )
!
!  土壌中の積分。
!    温度・水ともフラックスを確定する。
!
!  土壌中の水の移動に伴う熱については何らかの措置が必要。ここでは
!    (1) 各層の取得した熱フラックス収束を比熱の重みをかけて水に分配
!    (2) もとの温度と、もらった熱を持って、水は移動
!    (3) 行った先で温度変化・相変化計算
!  とする。他の方法としては
!    (1) 温度変化・相変化計算
!    (2) 新たな温度を持って水が移動
!    (3) 相変化計算のしなおし
!  が思い付く。

!  前者の方法をとることで、思い付く利点・欠点は以下の通り
!    利点
!      ・相変化計算は、(3) の一度でよい
!    問題点
!      ・0 度のとき、相変化に使われる熱を持って水が動いてしまう
!


      USE SIBPRM , ONLY :
     1   IDX   , 
     1   IDP   
!
      USE SIBCON , ONLY :
     1   CAP_WTR_c  , 
     1   TMP_FREZ_c   
     1  , rho_wtr_c 
     1  , cap_soil4_inv_c
!
      use sibveg , only :
     1  POR_SOIL_co   , DEP_SOIL_co 

!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS,INA, 
     1                       IMONIT_ALL , IMONIT_LEVEL ,
     1                       IMONIT_FLUX_COMMON, IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      REAL(8),INTENT(IN)    :: DELT_SOIL
!
      REAL(8),INTENT(IN)    :: HEAT_CONV(IDX,IDP+1)  ! 熱フラックス収束
!
      REAL(8),INTENT(IN)    :: WFLUX_H  (IDX,IDP)    ! 水フラックス下向き正
      REAL(8),INTENT(IN)    :: WTR_INC  (IDX,IDP)    ! 最終的な水増分
!
      REAL(8),INTENT(IN)    :: TMP_SOIL     (IDX,IDP)   
      REAL(8),INTENT(IN)    :: CAP_SOIL_INV (IDX,IDP) ! 層の比熱合計
!
#ifdef MONYOS
      INTEGER,INTENT(IN)    :: INSFLAG
#endif
!
! OUT 
!
      REAL(8)   :: HEAT_INC  (IDX,IDP+1)    ! 最終的にもらう熱
      REAL(8),INTENT(INOUT)   :: RUNOFF_B (IDX)        ! 下端でのランオフ
!
      REAL(8)   :: R_RUNOFF_WTR(IDX)   ! 最下層から入る水（負）
!
! INOUT 
!
      REAL(8),INTENT(INOUT) :: ENG_SOIL(IDX,IDP)
      REAL(8),INTENT(INOUT) :: H2O_SOIL(IDX,IDP)
      REAL(8),INTENT(INOUT) :: TMP_SOIL4(IDX)
!
! WORK
!
      INTEGER  :: I
      INTEGER  :: K
!
      REAL(8)  :: HFLUX_WH  (IDX,IDP)    ! 水移動に伴う熱フラックス
!
      REAL(8)  :: WFLUX_12                 ! 水フラックス
      REAL(8)  :: WFLUX_21  
      REAL(8)  :: WFLUX_23
      REAL(8)  :: WFLUX_32
      REAL(8)  :: WFLUX_34
      REAL(8)  :: WFLUX_43
!
      REAL(8)  :: FR_HEAT_12         ! 得る熱中の、去る水のもらう割合
      REAL(8)  :: FR_HEAT_21 
      REAL(8)  :: FR_HEAT_23 
      REAL(8)  :: FR_HEAT_32 
      REAL(8)  :: FR_HEAT_34 
!
!     REAL(8)  :: HFLUX_01(IDX)
!     REAL(8)  :: HFLUX_10(IDX)
      REAL(8)  :: HFLUX_12(IDX)
      REAL(8)  :: HFLUX_21(IDX)
      REAL(8)  :: HFLUX_23(IDX)
      REAL(8)  :: HFLUX_32(IDX)
      REAL(8)  :: HFLUX_34(IDX)
      REAL(8)  :: HFLUX_43(IDX)
!
      REAL(8),PARAMETER  :: DZERO = 0.D0
!
      LOGICAL,SAVE :: LFIRST
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
!     REAL(8)           :: WORK(IDX) 
#endif
!
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_FLUX_CONV 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=NSTR,NEND
!
! ==========================================
! >> 各層の取得した熱フラックス収束の分配 << 
! ==========================================
!
!  水の移動量
!
        WFLUX_12  = MAX ( WFLUX_H(I,1) , DZERO ) ! 正:下向きなら値あり
        WFLUX_21  = MIN ( WFLUX_H(I,1) , DZERO ) ! 負:上向きなら値あり
        WFLUX_23  = MAX ( WFLUX_H(I,2) , DZERO ) 
        WFLUX_32  = MIN ( WFLUX_H(I,2) , DZERO ) 
        WFLUX_34  = MAX ( WFLUX_H(I,3) , DZERO ) 
        WFLUX_43  = MIN ( WFLUX_H(I,3) , DZERO ) ! デバッグモニタ用
!
        RUNOFF_B    (I) = WFLUX_34 
        R_RUNOFF_WTR(I) = WFLUX_43                ! 最下層から入る水
!
!  持って行く熱の割合 
!     FR_HEAT_12 は、1 に来た熱のうちで、 2 に行く水が受け取る割合
!     符号は正である!!!
!     CAP_WTR_c * DELT_SOIL * WFLUX_12 等は、移動する水の比熱。
!     CAP_SOIL_INV は層全体の比熱の逆数。
!
        FR_HEAT_12 = CAP_WTR_c * DELT_SOIL * WFLUX_12  
     1                         * CAP_SOIL_INV(I,1)  
        FR_HEAT_21 = CAP_WTR_c * DELT_SOIL * ( - WFLUX_21 ) 
     1                         * CAP_SOIL_INV(I,2)  
        FR_HEAT_23 = CAP_WTR_c * DELT_SOIL * WFLUX_23  
     1                         * CAP_SOIL_INV(I,2)  
        FR_HEAT_32 = CAP_WTR_c * DELT_SOIL * ( - WFLUX_32 ) 
     1                         * CAP_SOIL_INV(I,3)  
        FR_HEAT_34 = CAP_WTR_c * DELT_SOIL * WFLUX_34  
     1                         * CAP_SOIL_INV(I,3)  
!
! 嫁いで行く水につけてやる熱のフラックス? HFLUX_WH
!     時間刻はかかっていないことに注意。
!
        HFLUX_12(I) 
     1      =   WFLUX_12 * CAP_WTR_c * ( TMP_SOIL(I,1) - TMP_FREZ_c )
     1        + FR_HEAT_12 * HEAT_CONV(I,1) 
        HFLUX_21(I) 
     1      =   WFLUX_21 * CAP_WTR_c * ( TMP_SOIL(I,2) - TMP_FREZ_c )
     1        - FR_HEAT_21 * HEAT_CONV(I,2) 
        HFLUX_WH(I,1) = HFLUX_12(I) + HFLUX_21(I)  
!
        HFLUX_23(I) 
     1      =   WFLUX_23 *  CAP_WTR_c * ( TMP_SOIL(I,2) - TMP_FREZ_c )
     1        + FR_HEAT_23 * HEAT_CONV(I,2) 
        HFLUX_32(I) 
     1      =   WFLUX_32 *  CAP_WTR_c * ( TMP_SOIL(I,3) - TMP_FREZ_c )
     1        - FR_HEAT_32 * HEAT_CONV(I,3) 
        HFLUX_WH(I,2) = HFLUX_23(I) + HFLUX_32(I)  
!
        HFLUX_34(I) 
     1      =   WFLUX_34  * CAP_WTR_c * ( TMP_SOIL(I,3) - TMP_FREZ_c )
     1        + FR_HEAT_34 * HEAT_CONV(I,3) 
        HFLUX_43(I) 
     1      =   WFLUX_43  *  CAP_WTR_c * ( TMP_SOIL(I,3) - TMP_FREZ_c )
        HFLUX_WH(I,3) = HFLUX_34(I) + HFLUX_43(I)
!
!       上記について。 
!       もし逆ランオフがあるとすれば、
!       仕方が無いので第 3 層の温度の水が来るとするとした。
!
!       第 4 層には H2O がなく、相変化も考慮していないため。
!
!       もし 4 層目に H2O があるなら下のようになる??
!    1      + WFLUX_43 * CAP_WTR_c * ( TMP_SOIL(I,4) - TMP_FREZ_c )
!    1      - FR_HEAT_43 * HEAT_IN_4 
!
!  もらう熱の合計
!
        HEAT_INC(I,1) = (   HEAT_CONV(I,1) 
     1                                    - HFLUX_WH(I,1) ) * DELT_SOIL
! 
        HEAT_INC(I,2) = (   HEAT_CONV(I,2) 
     1                    + HFLUX_WH(I,1) - HFLUX_WH(I,2) ) * DELT_SOIL
! 
        HEAT_INC(I,3) = (   HEAT_CONV(I,3) 
     1                    + HFLUX_WH(I,2) - HFLUX_WH(I,3) ) * DELT_SOIL
!
        HEAT_INC(I,4) = HEAT_CONV(I,4) * DELT_SOIL
        TMP_SOIL4(I)  = TMP_SOIL4(I) + HEAT_INC(I,4) * CAP_SOIL4_INV_c 
!
!          土壌第四層は水を持たず、第三層との水のやりとりを考えないため。
!
 1000 CONTINUE  
!
! ============================
! >> 各層のエネルギーの変更 << 
! ============================
!
      DO K=1,IDP
      DO I=NSTR,NEND
        ENG_SOIL(I,K) = ENG_SOIL(I,K) + HEAT_INC(I,K) 
        H2O_SOIL(I,K) = H2O_SOIL(I,K) + WTR_INC (I,K) 
      ENDDO 
      ENDDO
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
!
! ----------------------
! > 格子平均フラックス <
! ----------------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN   
        CMARK = 'ROFB' 
        IF ( INSFLAG.EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , RUNOFF_B, ICN, DELT_SOIL, ONE )
        ELSE
          CALL MONIT_SIB ( CMARK , RUNOFF_B, ICS, DELT_SOIL, ONE )
        ENDIF
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN   
        CMARK = 'A_ROB' 
        CALL MONIT_SIB ( CMARK , RUNOFF_B, INA, DELT_SOIL, ONE )
      ENDIF

!
! --------------------
! > 個々のフラックス <
! --------------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN   
!
! ◎水のフラックス。
!   3 層目とその下のみ、ここでは扱う。
!   それ以外は soil_wflux で格納済。
!   上向きと下向きの区別をするため。
!
!   層間の水フラックスに伴う熱
!
        IF ( INSFLAG.EQ.0 ) THEN
!         CMARK = 'IFLN10' 
!         CALL MONIT_SIB ( CMARK ,  HFLUX_10, ICN, DELT_SOIL , ONE)
          CMARK = 'IFLN12' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_12, ICN, DELT_SOIL , ONE)
          CMARK = 'IFLN21' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_21, ICN, DELT_SOIL , ONE)
          CMARK = 'IFLN23' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_23, ICN, DELT_SOIL , ONE)
          CMARK = 'IFLN32' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_32, ICN, DELT_SOIL , ONE)
        ELSE
!         CMARK = 'IFLS10' 
!         CALL MONIT_SIB ( CMARK ,  HFLUX_10, ICS, DELT_SOIL, ONE)      
          CMARK = 'IFLS12' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_12, ICS, DELT_SOIL, ONE)
          CMARK = 'IFLS21' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_21, ICS, DELT_SOIL, ONE)
          CMARK = 'IFLS23' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_23, ICS, DELT_SOIL, ONE)
          CMARK = 'IFLS32' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_32, ICS, DELT_SOIL, ONE)
        ENDIF
!
! 下面ランオフ と、その水のエネルギー
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'WRONB' 
          CALL MONIT_SIB ( CMARK , RUNOFF_B, ICN, DELT_SOIL, ONE )
          CMARK = 'IRONB' 
          CALL MONIT_SIB ( CMARK , HFLUX_34, ICN, DELT_SOIL, ONE )
        ELSE
          CMARK = 'WROSB' 
          CALL MONIT_SIB ( CMARK , RUNOFF_B, ICS, DELT_SOIL, ONE )
          CMARK = 'IROSB' 
          CALL MONIT_SIB ( CMARK , HFLUX_34, ICS, DELT_SOIL, ONE )
        ENDIF
!
!  4 -> 3 層への水 と そのエネルギー
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'WFLN43' 
          CALL MONIT_SIB ( CMARK , R_RUNOFF_WTR,ICN, DELT_SOIL, ONE)
          CMARK = 'IFLN43' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_43  , ICN, DELT_SOIL, ONE)
        ELSE
          CMARK = 'WFLS43' 
          CALL MONIT_SIB ( CMARK , R_RUNOFF_WTR,ICS, DELT_SOIL, ONE)
          CMARK = 'IFLS43' 
          CALL MONIT_SIB ( CMARK ,  HFLUX_43 , ICS, DELT_SOIL, ONE)
        ENDIF
!
      ENDIF
#endif
      RETURN
      END
!===============================================================
      SUBROUTINE MR_SOIL_FLUX_CONV
!
! soil_mod で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , 
     1                       IMONIT_ALL , IMONIT_LEVEL ,
     1                       IMONIT_FLUX_COMMON, IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SOIL_FLUX_CONV 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! ◎層間の水の移動に伴う熱の移動
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!  ランオフ
!     ROFD = WRONB + WROSB 
!
!
        CMARK  = 'ROFB'
        CTITLE = 'Gravitational Bottom Runoff'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN   
!
!  雪無し
!
!     CMARK  = 'IFLN10'
!     CTITLE = 'ENERGY WITH WTR FROM SOIL 1 TO SOIL 0 NOS'
!     CUNIT  = 'W/M**2'
!     CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLN21'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 2 TO SOIL 1 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLN32'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 3 TO SOIL 2 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLN12'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 1 TO SOIL 2 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLN23'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 2 TO SOIL 3 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  雪あり
!
!     CMARK  = 'IFLS10'
!     CTITLE = 'ENERGY WITH WTR FROM SOIL 1 TO SOIL 0 SNW'
!     CUNIT  = 'W/M**2'
!     CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLS21'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 2 TO SOIL 1 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLS32'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 3 TO SOIL 2 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLS12'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 1 TO SOIL 2 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLS23'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 2 TO SOIL 3 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 地下ランオフ と その熱
!
      CMARK  = 'WRONB'
      CTITLE = 'RUNOFF FROM BOTTOM NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WROSB'
      CTITLE = 'RUNOFF FROM BOTTOM SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IRONB'
      CTITLE = 'HEAT OF RUNOFF FROM BOTTOM NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IROSB'
      CTITLE = 'HEAT OF RUNOFF FROM BOTTOM SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  水フラックス 4->3  とその熱
!
      CMARK  = 'WFLN43'
      CTITLE = 'H2O FROM BOTTOM TO SOIL 3 NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS43'
      CTITLE = 'H2O FROM BOTTOM TO SOIL 3 SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )

      CMARK  = 'IFLN43'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 4 TO SOIL 3 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLS43'
      CTITLE = 'ENERGY WITH WTR FROM SOIL 4 TO SOIL 3 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_ROB'
        CTITLE = 'ALMA Qsb : Gravitational Bottom Runoff'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF

#endif
!
      RETURN
      END



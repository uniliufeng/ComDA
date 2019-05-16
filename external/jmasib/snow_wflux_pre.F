! 水まわり準備
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/31 保坂征宏 rho_newsnow_inv の初期値を与えるようにする。
!                   後で参照されることがあるため。参照のされ方チェック要。
!
      SUBROUTINE SNOW_WFLUX_PRE (
     I   NMAX    , DELT_SNOW         , 
     I   NMAX_0  , NMAX_01 , NMAX_2  , NMAX_3  , NMAX_4  , 
     I   GLA_SNOW, WTR_SNOW, TMP_SNOW, 
     I   CVR_SNW , 
     I   HFLUX_CONV     , 
     I   H2O_BUCKET_INTO_SNOW        ,ENH_BUCKET_INTO_SNOW        ,
     I   ENG_SNOW_BUCKET , H2O_SNOW_BUCKET , 
     I   U_REF       , V_REF       , TMP_REF      , 
     I   RHO_SNOW_INV_IN   , 
     O   ENG_PRECIPI , H2O_PRECIPI , RHO_PRECIPI_INV  , 
     O   ENG_SNOW    , H2O_SNOW    , RHO_SNOW_INV , 
     O   IFLG_SLEV_CHANGE      )                      ! 初期化するだけ

!
! 注意
!
!  ・ ここで、格子平均値として扱うことに注意。
!     01 については *CVR_SNW している。
!     降雨・降雪はもうされているとした。
!
!  ・ 雪無し域に雪が降ったものは、旧雪の切り直しと同様に扱う。
!
      USE SIBPRM , ONLY :
     1   IDX   , 
     1   ISN  
!
      USE SIBCON , ONLY : 
     1    CAP_WTR_c ,                             ! 水の比熱
     1    CAP_GLA_c ,                             ! 水の比熱
     1    HL_W2G_c    ,                             ! 氷 -> 水
     1    TMP_FREZ_c                              ! 凍結温度
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, INA, ICS ,
     1                      IMONIT_LEVEL , IMONIT_ALL   
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR 
#endif

!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_0 (3)
      INTEGER,INTENT(IN)   :: NMAX_01(3)
      INTEGER,INTENT(IN)   :: NMAX_2 (3)
      INTEGER,INTENT(IN)   :: NMAX_3 (3)
      INTEGER,INTENT(IN)   :: NMAX_4 (3)
!
      REAL(8),INTENT(IN)   :: CVR_SNW  (IDX)
!
!  バケツ （雪無し域も通した格子平均、TMP_SNOW がかかっている）
! 
      REAL(8),INTENT(IN)   :: ENH_BUCKET_INTO_SNOW (IDX)! 雨水エネルギー
      REAL(8),INTENT(IN)   :: H2O_BUCKET_INTO_SNOW (IDX)    ! 雨の量
!
      REAL(8),INTENT(IN)   :: TMP_SNOW  (IDX,ISN+1)   ! 雪の温度
      REAL(8),INTENT(IN)   :: GLA_SNOW  (IDX,ISN)     ! 雪の氷の量
      REAL(8),INTENT(IN)   :: WTR_SNOW  (IDX,ISN)     ! 雪の水の量
!
      REAL(8),INTENT(IN)   :: HFLUX_CONV (IDX,ISN) 
!                                                !雪層の受け取るエネルギー
!
      REAL(8),INTENT(IN)   :: RHO_SNOW_INV_IN    (IDX,ISN)
      REAL(8),INTENT(IN)   :: DELT_SNOW  
!
      REAL(8),INTENT(IN)   :: ENG_SNOW_BUCKET(IDX)
      REAL(8),INTENT(IN)   :: H2O_SNOW_BUCKET(IDX)
!
      REAL(8),INTENT(IN)   :: U_REF    (IDX)
      REAL(8),INTENT(IN)   :: V_REF    (IDX)
      REAL(8),INTENT(IN)   :: TMP_REF  (IDX)
!
! OUTPUT
!
!   雪積分用
      REAL(8),INTENT(OUT)  :: ENG_SNOW    (IDX,ISN)   ! 旧雪のエネルギー
      REAL(8),INTENT(OUT)  :: H2O_SNOW    (IDX,ISN)   ! 旧雪の持つ水
      REAL(8),INTENT(OUT)  :: ENG_PRECIPI (IDX)       ! 降雪のエネルギー
      REAL(8),INTENT(OUT)  :: H2O_PRECIPI (IDX)       ! 降雪の持つ水
!
      REAL(8),INTENT(OUT)  :: RHO_SNOW_INV   (IDX,ISN)
      REAL(8),INTENT(OUT)  :: RHO_PRECIPI_INV(IDX)
!
      INTEGER,INTENT(OUT)  :: IFLG_SLEV_CHANGE (IDX) ! 雪層の変更フラグ 
!
!
! WORK
!
      REAL(8)      :: ENX_SNOW1    ! 第一層の ENX 
      INTEGER      :: I
      INTEGER      :: K
      LOGICAL,SAVE :: LFIRST 
!
#ifdef MONYOS
      REAL(8),PARAMETER  :: ONE =  1.D0 
      CHARACTER(7)       :: CMARK
#endif
!
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_WFLUX_PRE POCHAKA 01/08/17' 
        LFIRST = .FALSE.
      ENDIF
!
!  初期化
!
      ENG_SNOW(:,:) = 0.D0
      H2O_SNOW(:,:) = 0.D0
      IFLG_SLEV_CHANGE (:) = 0        ! リセット
!
! ======================
! >> 降雪・降水の計算 <<   
! ======================
!  降雪・降水  ... 入力は格子平均値である。
!                  出力もそれで良い
!
      DO I=1,NMAX
        ENG_PRECIPI(I) =   ENH_BUCKET_INTO_SNOW (I) 
!
        H2O_PRECIPI(I) =   H2O_BUCKET_INTO_SNOW (I) 
!
        RHO_PRECIPI_INV(I) 
     1      = 1 / 
     1            ( 3.6 * SQRT(U_REF(I)**2 + V_REF(I)**2 )
     1              - 0.2 * ( TMP_REF(I) - TMP_FREZ_c ) 
     1              + 62 )
      ENDDO
!
!
! ==========================
! >> 雪の密度変化（加重） <<       
! ==========================
!   check  雪の密度積分（）をここでしてしまう？
!
      CALL SNOW_RHO_INTG (
     I    NMAX  , NMAX_01 , NMAX_2  , NMAX_3  , NMAX_4  , 
     I    RHO_SNOW_INV_IN , 
     I    DELT_SNOW       ,
     I    WTR_SNOW        , GLA_SNOW      , TMP_SNOW   , 
     O    RHO_SNOW_INV    )
!
! ============================
! >> 各層のエネルギーの計算 <<   
! ============================
!
!  第一層
!
      DO I=NMAX_01(2),NMAX
        ENG_SNOW (I,1) =  
     1      CAP_GLA_c * ( TMP_SNOW(I,1) - TMP_FREZ_c ) * GLA_SNOW (I,1) 
     1        - HL_W2G_c  * GLA_SNOW (I,1) 
     1        + HFLUX_CONV(I,1) * DELT_SNOW
        H2O_SNOW (I,1) = GLA_SNOW(I,1)    + WTR_SNOW(I,1) 
      ENDDO
!
!   01 は格子平均値にする （ CVR_SNW 倍 ）
!
      DO I=NMAX_01(2),NMAX_01(3)              
        ENG_SNOW (I,1) = ENG_SNOW (I,1) * CVR_SNW(I) 
        H2O_SNOW (I,1) = H2O_SNOW (I,1) * CVR_SNW(I) 
      ENDDO
!
!  第二層
!
      DO I=NMAX_2(2),NMAX
        ENG_SNOW (I,2) =  
     1      CAP_GLA_c * ( TMP_SNOW(I,2) - TMP_FREZ_c ) * GLA_SNOW (I,2) 
     1        - HL_W2G_c  * GLA_SNOW (I,2) 
     1        + HFLUX_CONV(I,2) * DELT_SNOW
        H2O_SNOW (I,2) = GLA_SNOW(I,2)    + WTR_SNOW(I,2) 
      ENDDO
!
!  第三層
!
      DO I=NMAX_3(2),NMAX
        ENG_SNOW (I,3) =  
     1      CAP_GLA_c * ( TMP_SNOW(I,3) - TMP_FREZ_c ) * GLA_SNOW (I,3) 
     1        - HL_W2G_c  * GLA_SNOW (I,3) 
     1        + HFLUX_CONV(I,3) * DELT_SNOW
        H2O_SNOW (I,3) = GLA_SNOW(I,3)    + WTR_SNOW(I,3) 
      ENDDO
!
!  第四層
!
      DO I=NMAX_4(2),NMAX
        ENG_SNOW (I,4) =  
     1      CAP_GLA_c * ( TMP_SNOW(I,4) - TMP_FREZ_c ) * GLA_SNOW (I,4) 
     1        - HL_W2G_c  * GLA_SNOW (I,4) 
     1        + HFLUX_CONV(I,4) * DELT_SNOW
        H2O_SNOW (I,4) = GLA_SNOW(I,4)    + WTR_SNOW(I,4) 
      ENDDO
!
!  ================
!  >> 年齢の補正 <<
!  ================
!  新雪層が熱で融けてきていたら、次のステップは旧雪 
!
! ==============
! > 詳細モニタ <
! ==============
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 水バケツ from canopy 
        CMARK  = 'WFLBWS' 
        CALL MONIT_SIB ( CMARK , H2O_BUCKET_INTO_SNOW, INA , ONE , ONE )
!
! 水バケツ from canopy に伴う熱
        CMARK  = 'IFLBWS' 
        CALL MONIT_SIB ( CMARK , ENH_BUCKET_INTO_SNOW, INA , ONE, ONE)
!
! 小雪バケツ  
!
        CMARK  = 'WLWBGET' 
        CALL MONIT_SIB ( CMARK , H2O_SNOW_BUCKET, INA , ONE   , ONE )  !ver2
!
        CMARK  = 'HLWBGET' 
        CALL MONIT_SIB ( CMARK , ENG_SNOW_BUCKET, INA , ONE   , ONE )  !ver2
!
! 瞬間値
!      
        IF ( ICN_SIB0109_MONITOR.EQ.1 ) THEN
!
        DO K=1,ISN
!
! エネルギー
!
          CMARK  = 'ENGWxY'                           ! スナップ
          WRITE ( CMARK(5:5) , '(I1)' ) K
          CALL MONIT_SIB ( CMARK , ENG_SNOW(1,K), ICS, DELT_SNOW, ONE )  
!
! H2O
!
          CMARK  = 'H2OWxY'                           ! スナップ
          WRITE ( CMARK(5:5) , '(I1)' ) K
          CALL MONIT_SIB ( CMARK , H2O_SNOW(1,K), ICS, DELT_SNOW, ONE )
!
        ENDDO  
        ENDIF
!
      ENDIF
#endif
!
! ==============================
! >> 少雪バケツを降雪に混ぜる <<
! ==============================
!
      DO I = 1,NMAX
        H2O_PRECIPI(I) = H2O_PRECIPI(I) + H2O_SNOW_BUCKET(I) 
        ENG_PRECIPI(I) = ENG_PRECIPI(I) + ENG_SNOW_BUCKET(I) 
      ENDDO
!
!  ====================
!  >> 0 層で雪等あり <<
!  ====================
!  
      DO I=1,NMAX_0(3)
        IF (       H2O_BUCKET_INTO_SNOW(I) .GT. 0. 
     1        .OR. H2O_SNOW_BUCKET     (I) .GT. 0. ) THEN 
          IFLG_SLEV_CHANGE(I) = 1
        ENDIF
      ENDDO
!

      RETURN
      END
!============================================================
      SUBROUTINE MR_SNOW_WFLUX_PRE
!
! SNOW_WFLUX_PRE で出力するデバッグ用変数の登録作業
!
      USE SIBPRM , ONLY :
     1   ISN 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL   
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
      INTEGER :: K 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_WFLUX_PRE 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 降水
!
      CMARK  = 'WFLBWS'
      CTITLE = 'BUCKET H2O INTO SNOW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 降水に伴う熱
!
      CMARK  = 'IFLBWS'
      CTITLE = 'ENG WITH BUCKET H2O INTO SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 少雪バケツ
!
      CMARK  = 'WLWBGET'
      CTITLE = 'H2O GET FROM LITTLE SNOW BUCKET '
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 少雪バケツ
!
      CMARK  = 'HLWBGET'
      CTITLE = 'ENG WITH H2O IN LITTLE SNOW BUCKET '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 中間での H2O snap
!
      CMARK  = 'H2OWxY'
      CTITLE = 'H2O SNOW x - BEFORE SNOW_WATER_CALC SNAP'
      CUNIT  = 'KG/M**2/S'
      DO K=1,ISN
        WRITE ( CMARK (5:5)   , '(I1)' ) K
        WRITE ( CTITLE(10:10) , '(I1)' ) K
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDDO

! 中間での ENG snap
!
      CMARK  = 'ENGWxY'
      CTITLE = 'ENG SNOW x - BEFORE SNOW_WATER_CALC SNAP'
      CUNIT  = 'W/M**2'
      DO K=1,ISN
        WRITE ( CMARK (5:5)   , '(I1)' ) K
        WRITE ( CTITLE(10:10) , '(I1)' ) K
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDDO
!
      ENDIF
#endif
!
      RETURN
      END

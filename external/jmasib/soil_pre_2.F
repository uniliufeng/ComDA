! 土壌の準備
!   ・蒸発・蒸散等により土壌トップから抜ける量を抜く
!   ・上から入る雨・融雪水を土壌にしみこませる
!   ・注意。上記はいずれも、DELT_SOIL のかかった量で計算する。
!     RUNOFF_S を出力するが、これは単位時間当たりの水フラックスにしたいので、
!     最後に単位変換を行なう。
!   
! 99/07/12 保坂征宏、hcap2.F より移動
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/10/04 保坂征宏 大幅変更。soil_pre にする。雪あり無し併用。
! 2000-09-27 豊田英司 不正な値で止める検査を厳しく
!
! 
!
! 古い cpp を使ってしまった場合の対応
#ifndef __LINE__
#  define __LINE__
#endif
!
      SUBROUTINE SOIL_PRE_2 (
     I    NSTR        , NEND        , INSFLAG  , 
     I    NSTR_2      , NEND_2      ,  
     I    CVR_NorS_INV   ,  
     I    NTYP             , 
     I    TMP_SOIL_IN    , TMP_SOIL4_IN , 
     I    FR_WTR_SOIL_IN , 
     I    FR_GLA_SOIL_IN , 
     I    FACT_WGT_H2O , FACT_WGT_H2O_INV , CAP_SOLID  , 
     I    DELT_SOIL   , 
     I    WTR_BUCKET_INTO    , ENH_BUCKET_INTO  , 
     I    WTR_BUCKET_IN_SN2GR, ENH_BUCKET_IN_SN2GR , 
     I    H2O_BUCKET_FROM    , ENH_BUCKET_FROM  , 
     I    TRNS_BUCKET    , 
     O    TMP_SOIL    , TMP_SOIL4   , 
     O    FR_WTR_SOIL , FR_GLA_SOIL , 
     O    CAP_SOIL    , CAP_SOIL_INV, 
     O    ENG_SOIL    , H2O_SOIL    , 
     O    PHASE_SOIL  , RUNOFF_S  )
!
!  注意 
!    ・蒸発・昇華・蒸散分は、はじめに取り除いてしまう。
!      キャノピー計算部ではあるとみこんで評価をしたものが、
!      実は熱計算をした後では無くなってしまっていた、というのでは
!      困るからである。
!      なお、本来、蒸発等で差し引かれる水にともなうエネルギーは分からない
!      はずだが、ここでは土壌の水がその温度で抜かれると仮定しているので、
!      TMP_SOIL からそのまま求めることにする。    
!
!    ・水を抜こうと思ったら水が無くて氷を抜く場合、そのエネルギーは
!      考える必要が無い。 
!             ... って、どういう意味だ? ( 01/09/04 保坂 )
!      蒸発・昇華・蒸散熱バケツと蒸発・昇華・蒸散水バケツに入っている
!      エネルギーと H2O だけ抜けば、削るものが水であれ氷であれ、同じだ
!      からである。（土壌のフラックス評価には効くだろうが、最後の温度や
!      相には直接は効かない。）
!
!    ・H2O そのものが足りない場合は、H2O を足りないだけたしてやる
!
!    ・熱の計算用に必要な変数を出力する
!      
      USE SIBPRM , ONLY : 
     1    IDX     , 
     1    IDP     , 
     1    ITYP 
!
      USE SIBCON , ONLY : 
     1    CAP_WTR_c  ,                            ! CW 
     1    CAP_GLA_c  ,                            ! CICE 
     1    TMP_FREZ_c  , 
     1    HL_W2G_c   ,
     1    RHO_WTR_C  , 
     1    CAP_SOIL4_c   
!
      USE SIBVEG , ONLY : 
     1           POR_SOIL_co            ,  ! 空隙率   (ITYP)  
     1           BEE_co                 ,  ! 土壌水分依存性係数 (ITYP)
     1           DEP_SOIL_co             , ! 土壌層の深さ    (ITYP,IDP)
     1           WCOEF_SAT_SOIL_co         ! 飽和透水係数    (ITYP)

!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN,ICS,INN,INS,INA,
     1                       IMONIT_LEVEL , IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON , IMONIT_ALMA 
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: NSTR_2           ! NMAX_0(2)
      INTEGER,INTENT(IN)    :: NEND_2           ! NMAX_0(3)
      INTEGER,INTENT(IN)    :: INSFLAG 
!
      REAL(8),INTENT(IN)    :: FACT_WGT_H2O     (IDX,IDP)  ! 空隙量x水の比重
      REAL(8),INTENT(IN)    :: FACT_WGT_H2O_INV (IDX,IDP)  ! 空隙量x水の比重
!
      REAL(8),INTENT(IN)    :: CVR_NorS_INV     (IDX)      ! 雪面被覆率
      INTEGER,INTENT(IN)    :: NTYP             (IDX)
!
      REAL(8),INTENT(IN)    :: TRNS_BUCKET      (IDX,IDP)
!
      REAL(8),INTENT(IN)    :: TMP_SOIL_IN      (IDX,IDP)
      REAL(8),INTENT(IN)    :: TMP_SOIL4_IN     (IDX)
!
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL_IN   (IDX,IDP)
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL_IN   (IDX,IDP)
!
      REAL(8),INTENT(IN)    :: H2O_BUCKET_FROM(IDX) 
      REAL(8),INTENT(IN)    :: ENH_BUCKET_FROM(IDX) 
!
      REAL(8),INTENT(IN)    :: CAP_SOLID (IDX,IDP) 
!
!   バケツ関連
!
      REAL(8),INTENT(IN)   ::  DELT_SOIL
      REAL(8),INTENT(IN)   ::  WTR_BUCKET_INTO         (IDX)
      REAL(8),INTENT(IN)   ::  ENH_BUCKET_INTO         (IDX)
!
!        ! 雪なし(0) 域で一旦少雪バケツから融けた水
      REAL(8),INTENT(IN)   ::  WTR_BUCKET_IN_SN2GR     (IDX) 
      REAL(8),INTENT(IN)   ::  ENH_BUCKET_IN_SN2GR     (IDX)
!
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: FR_WTR_SOIL  (IDX,IDP)
      REAL(8),INTENT(OUT)   :: FR_GLA_SOIL  (IDX,IDP)
!
      REAL(8),INTENT(OUT)   :: PHASE_SOIL   (IDX,IDP)  ! 零度の時 0 非 1
      REAL(8),INTENT(OUT)   :: TMP_SOIL     (IDX,IDP) 
      REAL(8),INTENT(OUT)   :: TMP_SOIL4    (IDX) 
      REAL(8),INTENT(OUT)   :: CAP_SOIL     (IDX,IDP) 
      REAL(8),INTENT(OUT)   :: CAP_SOIL_INV (IDX,IDP) 
!
      REAL(8),INTENT(INOUT) :: ENG_SOIL     (IDX,IDP) 
      REAL(8),INTENT(INOUT) :: H2O_SOIL     (IDX,IDP) 
!
      REAL(8),INTENT(INOUT) :: RUNOFF_S     (IDX)      ! 出力時は kg/m**2/s
!
!
! WORK
!
      REAL(8)   :: WTR_SOIL
      REAL(8)   :: GLA_SOIL
!
!     REAL(8)   :: XFR_WTR   (IDX)    ! W 空隙中で液水の占める体積比
!     REAL(8)   :: XPOROS    (IDX)    ! 空隙率 （1- 土 - 氷）
!
! 雪なし域での、少雪の融雪水 を
!
      REAL(8)   :: WTR_MOVE     (IDX) ! 雪無し域で 雪 からの水
      REAL(8)   :: ENH_MOVE     (IDX)   
!
! 蒸発などによる除去分
!
      REAL(8)   :: H2O_ADD_SOIL(IDX,IDP)
      REAL(8)   :: ENH_ADD_SOIL(IDX,IDP)
!
! 水の浸透とランオフ ( ランオフ量そのものはモニタ出力 ) 
!
      REAL(8)   :: WTR_INFILT_1 (IDX)
      REAL(8)   :: WTR_INFILT_2 (IDX)
      REAL(8)   :: ENH_INFILT_1 (IDX)
      REAL(8)   :: ENH_INFILT_2 (IDX)
      REAL(8)   :: ENH_RUNOFF_S (IDX) 
!
      REAL(8)   :: WTR_IN_TRY             ! 雨等
      REAL(8)   :: ENH_IN_TRY
      REAL(8)   :: WTR_INTO               ! 土壌1,2 層へ入る水
      REAL(8)   :: WTR_IN_MAX0            ! 1 層に入ろうとする上限
      REAL(8)   :: WTR_IN_MAX1            ! 1 層のすき間
      REAL(8)   :: WTR_IN_MAX2            ! 2 層のすき間
      REAL(8)   :: FAC_INFILT_1 
      REAL(8)   :: FAC_INFILT_2 
!
! 不足分の追加
!
      REAL(8)   :: ENG_MOD (IDX,IDP)
      REAL(8)   :: H2O_MOD (IDX,IDP)
!
! 温度計算用
!
      REAL(8)   :: FAC
      REAL(8)   :: ENX_SOIL
!
      INTEGER               :: I
      INTEGER               :: J
      INTEGER               :: K
!
      REAL(8),PARAMETER  :: ZERO = 0.D0
!
      REAL(8),PARAMETER :: SMALL  = 1.D-30   
!
#ifdef MONYOS
!     REAL(8),PARAMETER  :: DM1 = -1.D0 
      REAL(8),PARAMETER  :: ONE = 1.D0 
      REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
#ifdef SIB_DEBUG
      REAL(8)  :: H2O_STR(IDX) 
      REAL(8)  :: ENG_STR(IDX) 
      REAL(8)  :: H2O_END(IDX) 
      REAL(8)  :: ENG_END(IDX) 
      REAL(8)  :: H2O_DIFF_MIN
      REAL(8)  :: ENG_DIFF_MIN
#endif
!
      LOGICAL, SAVE:: LFIRST = .TRUE.
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_PRE 99/10/04 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  --------------
!  >> 初期状態 <<
!  --------------
!
      DO 1000 K=1,IDP
      DO 1000 I=NSTR,NEND
        WTR_SOIL         = FR_WTR_SOIL_IN(I,K) * FACT_WGT_H2O(I,K)
        GLA_SOIL         = FR_GLA_SOIL_IN(I,K) * FACT_WGT_H2O(I,K)
        CAP_SOIL(I,K)
     1     =    CAP_SOLID(I,K)
     1       +  CAP_WTR_c  * WTR_SOIL
     1       +  CAP_GLA_c  * GLA_SOIL
!
        H2O_SOIL (I,K) = WTR_SOIL + GLA_SOIL
        ENG_SOIL (I,K) 
     1     = CAP_SOIL(I,K) * ( TMP_SOIL_IN(I,K) - TMP_FREZ_c )
     1       - HL_W2G_c * GLA_SOIL
 1000 CONTINUE           
!
      TMP_SOIL4(:) = 0.
      TMP_SOIL4(NSTR:NEND) = TMP_SOIL4_IN(NSTR:NEND) 
!
#ifdef SIB_DEBUG
      H2O_STR(:) = 0.
      ENG_STR(:) = 0.
      DO K=1,IDP
      DO I=NSTR,NEND
        H2O_STR(I) = H2O_STR(I) + H2O_SOIL(I,K)
        ENG_STR(I) = ENG_STR(I) + ENG_SOIL(I,K)
      ENDDO 
      ENDDO 
!
!      DO I=NSTR,NEND
!        ENG_STR(I) = ENG_STR(I) + CAP_SOIL4_c * TMP_SOIL4(I)
!      ENDDO
!
      IF ( INSFLAG .EQ. 0 ) THEN 
        DO I=NSTR_2,NEND_2
          H2O_STR (I) = H2O_STR (I) + WTR_BUCKET_IN_SN2GR (I)
          ENG_STR (I) = ENG_STR (I) + ENH_BUCKET_IN_SN2GR (I)  
        ENDDO
      ENDIF
!
      DO I=NSTR,NEND
        H2O_STR (I) = H2O_STR(I) + H2O_BUCKET_FROM(I)*CVR_NorS_INV(I) 
        ENG_STR (I) = ENG_STR(I) + ENH_BUCKET_FROM(I)*CVR_NorS_INV(I) 
      ENDDO 
!
      DO K=1,IDP
      DO I=NSTR,NEND
        H2O_STR (I) = H2O_STR(I) 
     1        + TRNS_BUCKET (I,K) * CVR_NorS_INV (I) 
!
        ENG_STR (I) = ENG_STR(I) 
     1        + CAP_WTR_c * TRNS_BUCKET (I,K) * CVR_NorS_INV (I) 
     1           * ( TMP_SOIL_IN(I,K) - TMP_FREZ_c ) 
      ENDDO 
      ENDDO 
!
      DO I=NSTR,NEND
        H2O_STR (I) = H2O_STR(I) 
     1                + WTR_BUCKET_INTO(I) * CVR_NorS_INV(I)
        ENG_STR (I) = ENG_STR(I) 
     1                + ENH_BUCKET_INTO(I) * CVR_NorS_INV(I)
      ENDDO 
#endif        
!
! ======================
! >> バケツ修正の開始 << 
! ======================
!
      CALL RESET ( ENH_ADD_SOIL , IDX*IDP )
      CALL RESET ( H2O_ADD_SOIL , IDX*IDP )
!
      CALL RESET ( WTR_MOVE     , IDX     )
      CALL RESET ( ENH_MOVE     , IDX     )
!
! =============================================
! >> 雪被覆率 = 0 なのに雪層から水がくる場合 <<
! =============================================
!    雪が消えた後で、少雪バケツから水がくる場合など。
!
      IF ( INSFLAG.EQ.0 ) THEN
        DO I=NSTR_2,NEND_2
          WTR_MOVE (I) = WTR_BUCKET_IN_SN2GR (I)
          ENH_MOVE (I) = ENH_BUCKET_IN_SN2GR (I)  
        ENDDO
      ENDIF      
!
! ========================
! >> 強制的除去分の評価 << 
! ========================
!
! --------------------------
! > 蒸発・昇華による除去分 < 
! --------------------------
!    蒸発、昇華の水や氷を引く。
!    氷は昇華だけなので、K=1 の時点で決まる。
!
      DO 2000 I = NSTR , NEND
        H2O_ADD_SOIL (I,1) = H2O_BUCKET_FROM (I) * CVR_NorS_INV (I) 
        ENH_ADD_SOIL (I,1) = ENH_BUCKET_FROM (I) * CVR_NorS_INV (I) 
 2000 CONTINUE
!
! --------------------
! > 蒸散による除去分 < 
! --------------------
!
      DO 2100 K = 1,IDP 
      DO 2100 I = NSTR , NEND
        H2O_ADD_SOIL (I,K) 
     1      = H2O_ADD_SOIL(I,K) 
     1        + TRNS_BUCKET (I,K) * CVR_NorS_INV (I) 
!
        ENH_ADD_SOIL(I,K) 
     1      = ENH_ADD_SOIL(I,K)
     1        + CAP_WTR_c * TRNS_BUCKET (I,K) * CVR_NorS_INV (I) 
     1           * ( TMP_SOIL_IN(I,K) - TMP_FREZ_c ) 
 2100 CONTINUE
!
! ---------------------------------------------------------
! > 除去 < ( エラー処理 : 水の不足があれば加える を含む ) <
! ---------------------------------------------------------
!
      CALL RESET ( ENG_MOD      , IDX*IDP )
      CALL RESET ( H2O_MOD      , IDX*IDP )
!
      DO 2200 K = 1,IDP 
      DO 2200 I = NSTR , NEND
        FAC = 1.D0 
        IF ( H2O_SOIL(I,K) + H2O_ADD_SOIL(I,K) .LT. 0 ) THEN
          FAC = H2O_SOIL(I,K) / MAX ( - H2O_ADD_SOIL(I,K) , SMALL ) 
        ENDIF 
!
        H2O_SOIL (I,K) = H2O_SOIL(I,K) + H2O_ADD_SOIL(I,K) * FAC 
        ENG_SOIL (I,K) = ENG_SOIL(I,K) + ENH_ADD_SOIL(I,K) * FAC 
!
        H2O_MOD  (I,K) = H2O_ADD_SOIL(I,K) * (1-FAC) 
        ENG_MOD  (I,K) = ENH_ADD_SOIL(I,K) * (1-FAC) 
 2200 CONTINUE  
!
! ==================
! >> 中間チェック <<
! ==================
!
#ifdef SIB_DEBUG
      H2O_END(:) = 0.
      ENG_END(:) = 0.
      DO I=NSTR,NEND
        H2O_END(I) = WTR_BUCKET_INTO(I) * CVR_NorS_INV(I)
     1               + WTR_MOVE       (I) 
        ENG_END(I) = ENH_BUCKET_INTO(I) * CVR_NorS_INV(I)
     1               + ENH_MOVE       (I) 
      ENDDO 
!
      DO K=1,IDP
      DO I=NSTR,NEND
        H2O_END(I) = H2O_END(I) + H2O_SOIL(I,K) + H2O_MOD(I,K)
        ENG_END(I) = ENG_END(I) + ENG_SOIL(I,K) + ENG_MOD(I,K)
      ENDDO 
      ENDDO 
!

!
      DO I=NSTR,NEND
        IF ( ABS(H2O_STR(I)-H2O_END(I)) .GT. 1.D-9 ) THEN
          write(6,*) 'soil_pre_2 h2o middle error' , I , 
     1               H2O_STR(I)-H2O_END(I) , H2O_STR(I) , H2O_END(I),
     1               ENG_STR(I)-ENG_END(I) , ENG_STR(I) , ENG_END(I)
          stop 999
        ENDIF
        IF ( ABS(ENG_STR(I)-ENG_END(I)) .GT. 1.D-6 ) THEN
          write(6,*) 'soil_pre_2 eng middle error' , I , 
     1               H2O_STR(I)-H2O_END(I) , H2O_STR(I) , H2O_END(I),
     1               ENG_STR(I)-ENG_END(I) , ENG_STR(I) , ENG_END(I)
          stop 999
        ENDIF
      ENDDO 
#endif        
!
! ======================
! >> 上端からの浸透量 <<
! ======================
!   雪無し域 ～ 降水・結露 ( 結霜は雪面にいくので、ここではない )
!   雪あり域 ～ 融雪水
!
!   土壌第一層が 0 度より高ければ、可能なだけ水は浸透する。
!   ( 0 度以下の場合は空隙量を上限とする )
!   第二層目には、簡単のため、雨又は雪解け水の温度をもって浸透する。
!
      DO 3000 K = 1,2
      DO 3000 I = NSTR , NEND
        ENX_SOIL = ENG_SOIL(I,K) + HL_W2G_c * H2O_SOIL(I,K) 
!
        IF   ( ENG_SOIL(I,K) .GT. 0 ) THEN
          GLA_SOIL = 0. 
          WTR_SOIL = H2O_SOIL(I,K)
!
        ELSEIF ( ENX_SOIL .LT. 0 ) THEN
          GLA_SOIL = H2O_SOIL(I,K)
          WTR_SOIL = 0. 
!
        ELSE  
          FAC = ENX_SOIL / MAX ( HL_W2G_c*H2O_SOIL(I,K), SMALL ) 
                ! 水の割合 
          WTR_SOIL = H2O_SOIL(I,K) * FAC 
          GLA_SOIL = H2O_SOIL(I,K) * (1-FAC)
        ENDIF
!
        FR_WTR_SOIL(I,K) = WTR_SOIL * FACT_WGT_H2O_INV(I,K)
        FR_GLA_SOIL(I,K) = GLA_SOIL * FACT_WGT_H2O_INV(I,K)
 3000 CONTINUE  
!
!#ifdef DEBUG
!      CALL RESET ( WTR_INFILT_1  , IDX )
!      CALL RESET ( WTR_INFILT_2  , IDX )
!      CALL RESET ( ENH_INFILT_1  , IDX )
!      CALL RESET ( ENH_INFILT_2  , IDX )
!#endif
!
      DO 3100 I=NSTR,NEND
!
        J = NTYP(I)
!
! 入ろうとする水
!
        WTR_IN_TRY =   WTR_BUCKET_INTO(I) * CVR_NorS_INV(I)
     1               + WTR_MOVE       (I) 
        ENH_IN_TRY =   ENH_BUCKET_INTO(I) * CVR_NorS_INV(I)
     1               + ENH_MOVE       (I) 
!
! 浸透しうる量 ( 透水係数最大値 x 時間刻 、質量換算 ) 
!
        WTR_IN_MAX0 =  WCOEF_SAT_SOIL_co(J) 
     1                * ( 1 - FR_GLA_SOIL(I,1) ) 
!    1                * XFR_WTR(I)**(2.*BEE_co(J)+3.)   ! w 依存性無し
     1                * DELT_SOIL   
     1                * RHO_WTR_C                       ! 000731 
!
! 実際に土壌表面をつき抜けてくる量
!
        WTR_INTO  = MIN ( WTR_IN_TRY , WTR_IN_MAX0)
!
! ----------------------
! > 土壌第一層への浸透 <
! ----------------------
!
! 第一層のすき間 ( 水の質量にした場合 ) 
!
        WTR_IN_MAX1 = POR_SOIL_co(J) * DEP_SOIL_co(J,1)  ! (1-土-氷)* 体積
     1                * RHO_WTR_C             
     1                * ( 1. - FR_WTR_SOIL(I,1) -FR_GLA_SOIL(I,1) )     
!
! 第一層に入る水
!
        WTR_INFILT_1(I) = MIN ( WTR_INTO , WTR_IN_MAX1 ) 
!
! ----------------------
! > 土壌第二層への浸透 <
! ----------------------
!
! 第一層から第二層に入って来る分
!
        WTR_INTO = MAX ( WTR_INTO - WTR_INFILT_1(I) , ZERO ) 
!
! 第二層のすき間 ( 水の質量にした場合 ) 
!
        WTR_IN_MAX2 = POR_SOIL_co(J) * DEP_SOIL_co(J,2)  ! (1-土-氷)* 体積
     1                * RHO_WTR_C             
     1                * ( 1. - FR_WTR_SOIL(I,2) -FR_GLA_SOIL(I,2) )     
!
        WTR_INFILT_2(I) =  MIN ( WTR_INTO , WTR_IN_MAX2 )
!
! --------------
! > ランオフ等 <
! --------------
!
        FAC_INFILT_1 = WTR_INFILT_1(I) / MAX(WTR_IN_TRY,SMALL) 
        WTR_INFILT_1(I) = WTR_IN_TRY * FAC_INFILT_1
        ENH_INFILT_1(I) = ENH_IN_TRY * FAC_INFILT_1
!
        FAC_INFILT_2 = WTR_INFILT_2(I) / MAX(WTR_IN_TRY,SMALL) 
        WTR_INFILT_2(I) = WTR_IN_TRY * FAC_INFILT_2
        ENH_INFILT_2(I) = ENH_IN_TRY * FAC_INFILT_2
!
!        RUNOFF_S    (I) = WTR_IN_TRY * ( 1-FAC_INFILT_1-FAC_INFILT_2 )
!        ENH_RUNOFF_S(I) = ENH_IN_TRY * ( 1-FAC_INFILT_1-FAC_INFILT_2 )
        RUNOFF_S    (I) = WTR_IN_TRY 
     1         * max ( 1-FAC_INFILT_1-FAC_INFILT_2 , 0.D0 )
        ENH_RUNOFF_S(I) = ENH_IN_TRY 
     1         * max ( 1-FAC_INFILT_1-FAC_INFILT_2 , 0.D0 )
!
        H2O_SOIL(I,1) = H2O_SOIL(I,1) + WTR_INFILT_1(I) 
        ENG_SOIL(I,1) = ENG_SOIL(I,1) + ENH_INFILT_1(I) 
        H2O_SOIL(I,2) = H2O_SOIL(I,2) + WTR_INFILT_2(I) 
        ENG_SOIL(I,2) = ENG_SOIL(I,2) + ENH_INFILT_2(I) 
 3100 CONTINUE  
!
! ================
! >> 相に分ける <<
! ================
!
      DO 4000 K = 1,IDP
      DO 4000 I = NSTR , NEND
        ENX_SOIL  = ENG_SOIL(I,K) + HL_W2G_c * H2O_SOIL(I,K) 
        PHASE_SOIL(I,K) = 1.  
!
        IF   ( ENG_SOIL(I,K) .GT. 0 ) THEN
          GLA_SOIL = 0. 
          WTR_SOIL = H2O_SOIL(I,K)
          CAP_SOIL (I,K)  
     1     =    CAP_SOLID (I,K)
     1       +  CAP_WTR_c  * WTR_SOIL
     1       +  CAP_GLA_c  * GLA_SOIL
          CAP_SOIL_INV(I,K)  = 1./CAP_SOIL(I,K)
          TMP_SOIL(I,K) = TMP_FREZ_c + ENG_SOIL(I,K)*CAP_SOIL_INV(I,K) 
!
        ELSEIF ( ENX_SOIL .LT. 0 ) THEN
          GLA_SOIL = H2O_SOIL(I,K)
          WTR_SOIL = 0. 
          CAP_SOIL (I,K)  
     1     =    CAP_SOLID (I,K)
     1       +  CAP_WTR_c  * WTR_SOIL
     1       +  CAP_GLA_c  * GLA_SOIL
          CAP_SOIL_INV(I,K)  = 1./CAP_SOIL(I,K)
          TMP_SOIL(I,K) = TMP_FREZ_c + ENX_SOIL*CAP_SOIL_INV(I,K) 
!
        ELSE  
          FAC = ENX_SOIL / MAX ( HL_W2G_c*H2O_SOIL(I,K), SMALL ) 
                ! 水の割合 
          WTR_SOIL = H2O_SOIL(I,K) * FAC 
          GLA_SOIL = H2O_SOIL(I,K) * (1-FAC)
          CAP_SOIL (I,K)  
     1     =    CAP_SOLID (I,K)
     1       +  CAP_WTR_c  * WTR_SOIL
     1       +  CAP_GLA_c  * GLA_SOIL
          CAP_SOIL_INV(I,K)  = 1./CAP_SOIL(I,K)
          TMP_SOIL(I,K) = TMP_FREZ_c 
          PHASE_SOIL(I,K) = 0.  
        ENDIF
!
        FR_WTR_SOIL(I,K) = WTR_SOIL * FACT_WGT_H2O_INV(I,K)
        FR_GLA_SOIL(I,K) = GLA_SOIL * FACT_WGT_H2O_INV(I,K)
 4000 CONTINUE  
!
! ==================
! >> 終了チェック <<
! ==================
!
#ifdef SIB_DEBUG
      H2O_END(:) = 0.
      ENG_END(:) = 0.
!
      DO I=NSTR,NEND
        H2O_END(I) = RUNOFF_S     (I)
        ENG_END(I) = ENH_RUNOFF_S (I)
      ENDDO 
!
      DO K=1,IDP
      DO I=NSTR,NEND
        H2O_END(I) = H2O_END(I) + H2O_SOIL(I,K) + H2O_MOD(I,K) 
        ENG_END(I) = ENG_END(I) + ENG_SOIL(I,K) + ENG_MOD(I,K)
      ENDDO 
      ENDDO 
!
      H2O_DIFF_MIN = -1
      ENG_DIFF_MIN = -1
      DO I=NSTR,NEND
        IF ( ABS(H2O_STR(I)-H2O_END(I)) .GT. H2O_DIFF_MIN ) THEN
           H2O_DIFF_MIN = ABS(H2O_STR(I)-H2O_END(I))
        ENDIF 
        IF ( ABS(ENG_STR(I)-ENG_END(I)) .GT. ENG_DIFF_MIN ) THEN
           ENG_DIFF_MIN = ABS(ENG_STR(I)-ENG_END(I))
        ENDIF 
!
        IF ( ABS(H2O_STR(I)-H2O_END(I)) .GT. 1.D-9 ) THEN
          write(6,*) 'soil_pre_2 h2o error' , I , 
     1               H2O_STR(I)-H2O_END(I) , H2O_STR(I) , H2O_END(I),
     1               ENG_STR(I)-ENG_END(I) , ENG_STR(I) , ENG_END(I)
          stop 999
        ENDIF
        IF ( ABS(ENG_STR(I)-ENG_END(I)) .GT. 1.D-6 ) THEN
          write(6,*) 'soil_pre_2 eng error' , I , 
     1               H2O_STR(I)-H2O_END(I) , H2O_STR(I) , H2O_END(I),
     1               ENG_STR(I)-ENG_END(I) , ENG_STR(I) , ENG_END(I)
          stop 999
        ENDIF
      ENDDO 
      write(6,*) 'soil_pre_2 diff_max h2o ' , h2o_diff_min
      write(6,*) 'soil_pre_2 diff_max eng ' , eng_diff_min
#endif        

!
! ================
! >> 詳細モニタ <<   その 1 ver 1
! ================
!
#ifdef MONYOS
!
! 格子平均量 ( ランオフ )
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
        CMARK  = 'ROFS'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , RUNOFF_S, ICN, ONE, ONE )  
        ELSE
          CALL MONIT_SIB ( CMARK , RUNOFF_S, ICS, ONE, ONE )  
        ENDIF
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_ROS'
        CALL MONIT_SIB ( CMARK , RUNOFF_S, INA, ONE, ONE )  
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
        IF ( INSFLAG.EQ.0 ) THEN       
!
! 蒸発・昇華による増加氷量
!
          CMARK  = 'WEVBLN'                                  ! ver2 
          CALL MONIT_SIB ( CMARK , H2O_BUCKET_FROM,INN, ONE, ONE )   
!
! 蒸発・昇華による増加氷に伴うエネルギー
!
          CMARK  = 'IEVBLN' 
          CALL MONIT_SIB ( CMARK , ENH_BUCKET_FROM,INN, ONE, ONE )
!
        ENDIF
!
! 蒸散による増加水量、それに伴うエネルギー量
!
!  一層目
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WTRBL1N' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET,INN, ONE, ONE ) 
        ELSE
          CMARK  = 'WTRBL1S' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET,INS, ONE, ONE ) 
        ENDIF
!
        DO I=NSTR,NEND
          WORK(I) = TRNS_BUCKET(I,1) 
     1              * CAP_WTR_c * ( TMP_SOIL_IN(I,1) - TMP_FREZ_c ) 
        ENDDO
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'ITRBL1N' 
          CALL MONIT_SIB ( CMARK , WORK ,INN, ONE, ONE ) 
        ELSE
          CMARK  = 'ITRBL1S' 
          CALL MONIT_SIB ( CMARK , WORK ,INS, ONE, ONE ) 
        ENDIF
!
! 二層目
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WTRBL2N' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET(1,2),INN, ONE, ONE ) 
        ELSE
          CMARK  = 'WTRBL2S' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET(1,2),INS, ONE, ONE ) 
        ENDIF
!
        DO I=NSTR,NEND
          WORK(I) = TRNS_BUCKET(I,2) 
     1              * CAP_WTR_c * ( TMP_SOIL_IN(I,2) - TMP_FREZ_c ) 
        ENDDO
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'ITRBL2N' 
          CALL MONIT_SIB ( CMARK , WORK ,INN, ONE, ONE ) 
        ELSE
          CMARK  = 'ITRBL2S' 
          CALL MONIT_SIB ( CMARK , WORK ,INS, ONE, ONE ) 
        ENDIF
!
! 三層目
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WTRBL3N' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET(1,3),INN, ONE, ONE ) 
        ELSE
          CMARK  = 'WTRBL3S' 
          CALL MONIT_SIB ( CMARK , TRNS_BUCKET(1,3),INS, ONE, ONE ) 
        ENDIF
!
        DO I=NSTR,NEND
          WORK(I) = TRNS_BUCKET(I,3) 
     1              * CAP_WTR_c * ( TMP_SOIL_IN(I,3) - TMP_FREZ_c ) 
        ENDDO
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'ITRBL3N' 
          CALL MONIT_SIB ( CMARK , WORK ,INN, ONE, ONE ) 
        ELSE
          CMARK  = 'ITRBL3S' 
          CALL MONIT_SIB ( CMARK , WORK ,INS, ONE, ONE ) 
        ENDIF
!
! 上端から一旦入る水  = 雨・結露  or 融雪水
!
        DO I=NSTR,NEND
          WORK(I) = WTR_BUCKET_INTO(I) * CVR_NorS_INV(I) 
        ENDDO
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WFLBLN'
          CALL MONIT_SIB ( CMARK , WORK , ICN ,  ONE, ONE )        
        ELSE
          CMARK  = 'WFLBLW'
          CALL MONIT_SIB ( CMARK , WORK , ICS ,  ONE, ONE )        
        ENDIF
!
!    (雪 -> 雪無し土壌)
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WFLBLNW'
          CALL MONIT_SIB ( CMARK , WTR_MOVE , ICN , ONE, ONE )        
        ENDIF 
!
! 一旦入ろうとする水に伴うエネルギー
!
        DO I=NSTR,NEND
          WORK(I) = ENH_BUCKET_INTO(I) * CVR_NorS_INV(I) 
        ENDDO
        IF ( INSFLAG.EQ.0) THEN
          CMARK = 'IFLBLN' 
          CALL MONIT_SIB ( CMARK , WORK , ICN , ONE, ONE )  
        ELSE
          CMARK = 'IFLBLW' 
          CALL MONIT_SIB ( CMARK , WORK , ICS ,  ONE, ONE )  
        ENDIF
!
!    (雪 -> 雪無し土壌)
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'IFLBLNW'
          CALL MONIT_SIB ( CMARK , ENH_MOVE , ICN, ONE, ONE )        
        ENDIF 
!
! 第一層に入った水
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WINFLN1'
          CALL MONIT_SIB ( CMARK , WTR_INFILT_1 , ICN , ONE, ONE )        
        ELSE
          CMARK  = 'WINFLS1'
          CALL MONIT_SIB ( CMARK , WTR_INFILT_1 , ICS, ONE, ONE )        
        ENDIF
!
! 第一層に入った水に伴うエネルギー
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'IINFLN1'
          CALL MONIT_SIB ( CMARK , ENH_INFILT_1 , ICN, ONE, ONE )        
        ELSE
          CMARK  = 'IINFLS1'
          CALL MONIT_SIB ( CMARK , ENH_INFILT_1 , ICS, ONE, ONE )        
        ENDIF
!
! 第二層に入った水
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WINFLN2'
          CALL MONIT_SIB ( CMARK , WTR_INFILT_2 , ICN, ONE, ONE )        
        ELSE
          CMARK  = 'WINFLS2'
          CALL MONIT_SIB ( CMARK , WTR_INFILT_2 , ICS, ONE, ONE )        
        ENDIF
!
! 第二層に入った水に伴うエネルギー
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'IINFLN2'
          CALL MONIT_SIB ( CMARK , ENH_INFILT_2 , ICN, ONE, ONE )        
        ELSE
          CMARK  = 'IINFLS2'
          CALL MONIT_SIB ( CMARK , ENH_INFILT_2 , ICS, ONE, ONE )        
        ENDIF
!
! 表層ランオフ （正）
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WRONS'
          CALL MONIT_SIB ( CMARK , RUNOFF_S, ICN, ONE, ONE )  
        ELSE
          CMARK  = 'WROSS'
          CALL MONIT_SIB ( CMARK , RUNOFF_S, ICS, ONE, ONE )  
        ENDIF
!
! 表層ランオフした水に伴うエネルギー
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'IRONS'
          CALL MONIT_SIB ( CMARK , ENH_RUNOFF_S, ICN, ONE, ONE )  
        ELSE
          CMARK  = 'IROSS'
          CALL MONIT_SIB ( CMARK , ENH_RUNOFF_S, ICS, ONE, ONE )  
        ENDIF
!
      ENDIF
!
#endif
!
! ================
! >> 詳細モニタ <<   その 2  水・氷エラー処理分
! ================
!
#ifdef MONYOS 
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!
! 本当に与えた水と熱
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'WADLN1B'
          CALL MONIT_SIB ( CMARK , H2O_MOD ,ICN, ONE, ONE )
          CMARK  = 'WADLN2B'
          CALL MONIT_SIB ( CMARK , H2O_MOD(1,2), ICN, ONE, ONE )     
          CMARK  = 'WADLN3B'
          CALL MONIT_SIB ( CMARK , H2O_MOD(1,3), ICN, ONE, ONE )
        ELSE
          CMARK  = 'WADLS1B'
          CALL MONIT_SIB ( CMARK , H2O_MOD     , ICS , ONE, ONE )
          CMARK  = 'WADLS2B'
          CALL MONIT_SIB ( CMARK , H2O_MOD(1,2), ICS , ONE, ONE )     
          CMARK  = 'WADLS3B'
          CALL MONIT_SIB ( CMARK , H2O_MOD(1,3), ICS, ONE, ONE )
        ENDIF
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'HADLN1B'
          CALL MONIT_SIB ( CMARK , ENG_MOD ,     ICN , ONE, ONE )
          CMARK  = 'HADLN2B'
          CALL MONIT_SIB ( CMARK , ENG_MOD(1,2), ICN , ONE, ONE )     
          CMARK  = 'HADLN3B'
          CALL MONIT_SIB ( CMARK , ENG_MOD(1,3), ICN , ONE, ONE )
        ELSE
          CMARK  = 'HADLS1B'
          CALL MONIT_SIB ( CMARK , ENG_MOD     , ICS , ONE, ONE )
          CMARK  = 'HADLS2B'
          CALL MONIT_SIB ( CMARK , ENG_MOD(1,2), ICS , ONE, ONE )     
          CMARK  = 'HADLS3B'
          CALL MONIT_SIB ( CMARK , ENG_MOD(1,3), ICS , ONE, ONE )
        ENDIF
!
      ENDIF
#endif
!
! ==============
! >> 単位変換 <<
! ==============
!
      DO I=NSTR,NEND 
        RUNOFF_S(I) = RUNOFF_S(I) / DELT_SOIL 
      ENDDO
!
      RETURN
      END
!==================================================================
      SUBROUTINE MR_SOIL_PRE_2
!
! soil_pre で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL,
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
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SOIL_PRE 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! --------------------------------------------------------------
!
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!   ランオフ
!     ROFS = WRONS + WROSS 
        CMARK  = 'ROFS'
        CTITLE = 'Surface Runoff'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!  与えた熱
!
      CMARK  = 'HADLN1B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SOIL 1 BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HADLS1B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SNSL   BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HADLN2B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SOIL 2 NOS BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HADLS2B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SOIL 2 SNW BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HADLN3B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SOIL 3 NOS BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HADLS3B'
      CTITLE = 'ADD ENERGY WITH ADD WATER SOIL 3 SNW BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!  与えた水
!
      CMARK  = 'WADLN1B'
      CTITLE = 'ADD WATER SOIL 1 BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADLS1B'
      CTITLE = 'ADD WATER SNSL   BAD'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADLN2B'
      CTITLE = 'ADD WATER SOIL 2 BAD NOS'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADLS2B'
      CTITLE = 'ADD WATER SOIL 2 BAD SNW'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADLN3B'
      CTITLE = 'ADD WATER SOIL 3 BAD NOS'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADLS3B'
      CTITLE = 'ADD WATER SOIL 3 BAD SNW'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!---------------
! ver2
!
!  土壌 1 層からの蒸散水（負）
!
      CMARK  = 'WTRBL1N'
      CTITLE = 'TRANS FROM SOIL 1 NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBL1S'
      CTITLE = 'TRANS FROM SOIL 1 SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  土壌 2 層からの蒸散水（負）
!
      CMARK  = 'WTRBL2N'
      CTITLE = 'TRANS FROM SOIL 2 NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBL2S'
      CTITLE = 'TRANS FROM SOIL 2 SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  土壌 3 層からの蒸散水（負）
!
      CMARK  = 'WTRBL3N'
      CTITLE = 'TRANS FROM SOIL 3 NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBL3S'
      CTITLE = 'TRANS FROM SOIL 3 SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!---------------
! ver2
!
!  土壌 1 層からの蒸散水に伴うエネルギー
!
      CMARK  = 'ITRBL1N'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 1 NOS'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ITRBL1S'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 1 SNW'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  土壌 2 層からの蒸散水に伴うエネルギー
!
      CMARK  = 'ITRBL2N'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 2 NOS'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ITRBL2S'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 2 SNW'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  土壌 3 層からの蒸散水に伴うエネルギー
!
      CMARK  = 'ITRBL3N'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 3 NOS'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ITRBL3S'
      CTITLE = 'ENERGY WITH TRANS FROM SOIL 3 SNW'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 蒸発・昇華バケツ
!
      CMARK  = 'WEVBLN'
      CTITLE = 'EVAP/SUBL WATER FROM SOIL TO CNP_NOS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!
      CMARK  = 'IEVBLN'
      CTITLE = 'ENG WITH EVAP/SUBL WATER FROM SOIL TO CNP_NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  上端からいったんは入る水
!
      CMARK  = 'WFLBLN'
      CTITLE = 'H2O INTO SOIL 1 NOS'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLBLW'
      CTITLE = 'H2O INTO SOIL 1 SNW'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLBLNW'
      CTITLE = 'H2O INTO NOS SOIL FROM SOIL'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  上端からいったんは入る水に伴うエネルギー
!
      CMARK  = 'IFLBLN' 
      CTITLE = 'ENERGY WITH WTR FROM SOIL 0 TO SOIL 1 NOS'
      CUNIT  = 'W/M**2 '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLBLW' 
      CTITLE = 'ENERGY WITH WTR FROM SOIL 0 TO SOIL 1 SNW'
      CUNIT  = 'W/M**2 '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLBLNW'
      CTITLE = 'ENERGY WITH H2O INTO NOS SOIL FROM SOIL'
      CUNIT  = 'W/M**2 '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  第一層への浸透水
!
      CMARK  = 'WINFLN1'
      CTITLE = 'H2O INFILT INTO SOIL 1 NOS'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WINFLS1'
      CTITLE = 'H2O INFILT INTO SOIL 1 SNW'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  第二層への浸透水
!
      CMARK  = 'WINFLN2'
      CTITLE = 'H2O INFILT INTO SOIL 2 NOS'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WINFLS2'
      CTITLE = 'H2O INFILT INTO SOIL 2 SNW'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  第一層への浸透水に伴うエネルギー
!
      CMARK  = 'IINFLN1'
      CTITLE = 'ENG with H2O INFILT INTO SOIL 1 NOS'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IINFLS1'
      CTITLE = 'ENG_with H2O INFILT INTO SOIL 1 SNW'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  第二層への浸透水に伴うエネルギー
!
      CMARK  = 'IINFLN2'
      CTITLE = 'ENG with H2O INFILT INTO SOIL 2 NOS'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IINFLS2'
      CTITLE = 'ENG with H2O INFILT INTO SOIL 2 SNW'
      CUNIT  = 'kg/m**2/s '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )

!
!  表層ランオフ
!
      CMARK  = 'WRONS'
      CTITLE = 'SURFACE RUNOFF NOS'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WROSS'
      CTITLE = 'SURFACE RUNOFF SNW'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  表層ランオフに伴う熱
!
      CMARK  = 'IRONS'
      CTITLE = 'ENERGY OF SURFACE RUNOFF NOS'
      CUNIT  = 'W/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IROSS'
      CTITLE = 'ENERGY OF SURFACE RUNOFF SNW'
      CUNIT  = 'W/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_ROS'
        CTITLE = 'ALMA Qs : Surface Runoff'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END

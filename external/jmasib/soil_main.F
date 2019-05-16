! 土壌の計算 ... 積分まで
! 99/07/12 保坂征宏  第一版
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 2000-09-13 豊田英司 土壌水分多すぎバグの調査
!                     -DDEBUG_PARANOIA をつけると土壌水分が 1 を超えてよい
!                     ところでも検査を行う

      SUBROUTINE SOIL_MAIN  (
     I    NMAX      , NMAX_0      , NMAX_01     , 
     I    NTYP           , 
     I    DELT_SOIL      , 
     I    CVR_SNW        , CVR_SNW_NEXT   , 
     I    CVR_SNW_INV    , CVR_NOS_INV    , 
     I    CVR_SNW_INV_NEXT ,  CVR_NOS_INV_NEXT ,  
     I    HDF_SOIL_H     , HDF_SNSL_H      , 
     I    ENG_BUCKET_NOS , ENG_BUCKET_SNOW_TO_SNSL  , 
     I    WTR_BUCKET_INTO_SOIL     , ENH_BUCKET_INTO_SOIL , 
     I    H2O_BUCKET_FROM_SOIL     , ENH_BUCKET_FROM_SOIL , 
     I    WTR_BUCKET_SNOW_TO_SNSL  , ENH_BUCKET_SNOW_TO_SNSL   , 
     I    TRNS_BUCKET_NOS          , TRNS_BUCKET_SNW        , 
     I    INFO_SNOW_NEXT , 
     M    TMP_SOIL_PROG   , TMP_SNSL_PROG   ,      ! 予報変数 : 温度
     M    FR_WTR_SOIL_PROG, FR_WTR_SNSL_PROG,      ! 予報変数 : 水の割合
     M    FR_GLA_SOIL_PROG, FR_GLA_SNSL_PROG,      ! 予報変数 : 氷の割合
     M    TMP_SOIL4_PROG  ,
     O    RUNOFF )
!
! していないこと
!
!   ・熱伝導係数計算
!        ... 雪ともからむので、その前に
!
! 土壌に関する全体的な注意
!   ・土壌中では水が氷になるときに、体積が増加してほしくない。
!     そこで土壌中体積から質量に換算する際の密度としては 
!     RHO_GLA_c は使わず、RHO_WTR_c を用いる。
!     しかし
!         比熱を使い体積からエネルギー換算する場合
!         潜熱係数を使い体積からエネルギー換算する場合
!         ( 熱伝導係数の計算 ? )
!     は RHO_GLA_c/RHO_WTR_c である RHO_FAC_GW をかけたものを使った方が
!     いいのかも知れない。
!     つまり
!        比熱     CAP_GLA_c のかわりに 
!           CAP_GLA_GW_c ( CAP_WTR_c の誤り?? 保坂 01/09/20 )を使う
!        潜熱係数 HL_W2G_c    のかわりに
!           HL_GLA_GW    ( ??  保坂 01/09/20 )を使う
!     方がいいのかも知れない。
!     と一瞬思ったが、比熱・潜熱係数ではその必要はない気がする。
!     問題があるとすれば熱伝導係数か?
!     このあたりはちょっとペンディングして、CAP_GLA_c を使ってしまう。
!
! ○ 上端からの入力量はいずれも、格子平均値である。
!      * CVR_SNW , または * CVR_NOS された値であることに注意
!
      USE SIBPRM , ONLY :
     1    IDX   , 
     1    IDP   ,
     1    ITYP 
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_SOIL_HEAT_EXPLICIT       ! implicit -> 0 , explicit -> 1 
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NMAX 
      INTEGER,INTENT(IN)    :: NMAX_0 (3) 
      INTEGER,INTENT(IN)    :: NMAX_01(3) 
!
      INTEGER,INTENT(IN)    :: NTYP           (IDX)     ! 植生タイプ
!
      REAL(8),INTENT(IN)    :: DELT_SOIL                ! 土壌用時間刻
!
      REAL(8),INTENT(IN)    :: HDF_SOIL_H   (IDX,IDP) ! 熱伝導係数
      REAL(8),INTENT(IN)    :: HDF_SNSL_H   (IDX,IDP)
!
      REAL(8),INTENT(IN)    :: CVR_SNW       (IDX)          ! 雪面被覆率
      REAL(8),INTENT(IN)    :: CVR_SNW_INV   (IDX)          ! 雪面被覆率
      REAL(8),INTENT(IN)    :: CVR_NOS_INV   (IDX)          ! 雪面被覆率
      REAL(8),INTENT(IN)    :: CVR_SNW_NEXT  (IDX)          ! 雪面被覆率
      REAL(8),INTENT(IN)    :: CVR_NOS_INV_NEXT (IDX)
      REAL(8),INTENT(IN)    :: CVR_SNW_INV_NEXT (IDX)
!
      REAL(8),INTENT(IN)    :: ENG_BUCKET_NOS (IDX)         ! 上端に入る熱量
      REAL(8),INTENT(IN)    :: ENG_BUCKET_SNOW_TO_SNSL(IDX) ! 
!
      REAL(8),INTENT(IN)    :: WTR_BUCKET_INTO_SOIL   (IDX) ! 上端に入る水量
      REAL(8),INTENT(IN)    :: H2O_BUCKET_FROM_SOIL   (IDX) ! 上端から出る水量
      REAL(8),INTENT(IN)    :: WTR_BUCKET_SNOW_TO_SNSL(IDX)
      REAL(8),INTENT(IN)    :: ENH_BUCKET_INTO_SOIL   (IDX) ! 上端から入る水E
      REAL(8),INTENT(IN)    :: ENH_BUCKET_FROM_SOIL   (IDX) ! 
      REAL(8),INTENT(IN)    :: ENH_BUCKET_SNOW_TO_SNSL(IDX)
!
      REAL(8),INTENT(IN)    :: TRNS_BUCKET_NOS       (IDX,IDP) ! 蒸散
      REAL(8),INTENT(IN)    :: TRNS_BUCKET_SNW       (IDX,IDP) ! 
!
!
      REAL(8),INTENT(INOUT) :: TMP_SOIL_PROG   (IDX,IDP)  ! 温度
      REAL(8),INTENT(INOUT) :: TMP_SNSL_PROG   (IDX,IDP)
      REAL(8),INTENT(INOUT) :: FR_WTR_SOIL_PROG(IDX,IDP)  ! 空隙中の水の比率
      REAL(8),INTENT(INOUT) :: FR_WTR_SNSL_PROG(IDX,IDP)
      REAL(8),INTENT(INOUT) :: FR_GLA_SOIL_PROG(IDX,IDP)  ! 氷の比率
      REAL(8),INTENT(INOUT) :: FR_GLA_SNSL_PROG(IDX,IDP)
!
      REAL(8),INTENT(INOUT) :: TMP_SOIL4_PROG  (IDX)
!
      INTEGER,INTENT(IN) :: INFO_SNOW_NEXT      (IDX) 
!
! OUTPUT
!

      REAL(8),INTENT(OUT)   :: RUNOFF         (IDX)      ! ランオフ
!
!
! WORK ( NOS/SNW 別 )
!
      REAL(8)    :: ENG_SOIL_NOS  (IDX,IDP)  ! 層のエネルギー
      REAL(8)    :: ENG_SOIL_SNW  (IDX,IDP)
      REAL(8)    :: H2O_SOIL_NOS  (IDX,IDP)  ! 層の H2O 
      REAL(8)    :: H2O_SOIL_SNW  (IDX,IDP)
!
      REAL(8)    :: ENG_SOIL_NEXT (IDX,IDP)  ! 層のエネルギー
      REAL(8)    :: H2O_SOIL_NEXT (IDX,IDP)
      REAL(8)    :: ENG_SNSL_NEXT (IDX,IDP)  ! 層の H2O 
      REAL(8)    :: H2O_SNSL_NEXT (IDX,IDP)
!
      REAL(8)    :: RUNOFF_S_NOS  (IDX)      ! ランオフ(水フラックス)
      REAL(8)    :: RUNOFF_S_SNW  (IDX)
      REAL(8)    :: RUNOFF_B_NOS  (IDX)
      REAL(8)    :: RUNOFF_B_SNW  (IDX)
!
! WORK ( NOS/SNW 使い回し )
!
      REAL(8)    :: TMP_SOIL       (IDX,IDP)  ! 土壌温度 
      REAL(8)    :: FR_WTR_SOIL    (IDX,IDP)  ! 水の比率 
      REAL(8)    :: FR_GLA_SOIL    (IDX,IDP)  ! 氷の比率 
!
      REAL(8)    :: SIB_MTX_SOIL(IDX,IDP,IDP+1)
      REAL(8)    :: CAP_SOIL    (IDX,IDP)    ! 熱容量    
      REAL(8)    :: CAP_SOIL_INV(IDX,IDP)    
!
      REAL(8)    :: PHASE_SOIL  (IDX,IDP)    ! 零度 0 非 1
      REAL(8)    :: HEAT_CONV   (IDX,IDP+1)  ! 層が得る熱フラックス収束
      REAL(8)    :: HFLUX_H     (IDX,0:IDP)  ! 
!
      REAL(8)    :: WFLUX_H     (IDX,IDP)  ! 水フラックス 
      REAL(8)    :: WTR_INC     (IDX,IDP)  ! 水フラックス収束 * DELT  
!
      REAL(8)    :: TMP_SOIL4_N (IDX)      ! 
      REAL(8)    :: TMP_SOIL4_S (IDX)      ! 
!
! WORK ( 定数的変数 ) 
!
      REAL(8)    :: CAP_SOLID     (IDX,IDP)    ! 固体土の比熱
      REAL(8)    :: FACT_WGT_H2O  (IDX,IDP)    ! 空隙量x水の比重
!                                              ! 水や氷の割合をかけると質量に
      REAL(8)    :: FACT_WGT_H2O_INV (IDX,IDP) ! 空隙量 x 水の比重の逆数
!
      REAL(8)    :: DELT_SOIL_INV 
!
#ifdef SIB_DEBUG
      INTEGER    :: I , K
#endif
!
! その他
!
      REAL(8),PARAMETER :: DUMMY_IDX      (IDX) = 0.
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_MAIN 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DELT_SOIL_INV = 1/DELT_SOIL
!
! ==========
! >> 準備 <<
! ==========
!
! 諸量のリセット
!
      CALL RESET ( ENG_SOIL_NOS , IDX*IDP )
      CALL RESET ( ENG_SOIL_SNW , IDX*IDP )
      CALL RESET ( H2O_SOIL_NOS , IDX*IDP )
      CALL RESET ( H2O_SOIL_SNW , IDX*IDP )
!
      RUNOFF_S_NOS(:) = 0.
      RUNOFF_B_NOS(:) = 0.
      RUNOFF_S_SNW(:) = 0.
      RUNOFF_B_SNW(:) = 0.
!
      TMP_SOIL4_N (:) = 0.
      TMP_SOIL4_S (:) = 0.
!
!   比熱と空隙量の計算
!
      CALL SOIL_PRE_1 (
     I    NMAX            , NTYP       , 
     O    FACT_WGT_H2O     , 
     O    FACT_WGT_H2O_INV ,     
     O    CAP_SOLID      )
!
! ==============
! >> 雪無し域 <<
! ==============
!
      BARESOIL_GRIDS: IF (NMAX_01(3) > 0) THEN
!
      CALL SOIL_PRE_2 (
     I    1              , NMAX_01(3)        , 0       , 
     I    1              , NMAX_0(3)         , 
     I    CVR_NOS_INV    ,  NTYP             , 
     I    TMP_SOIL_PROG  , TMP_SOIL4_PROG   ,
     I    FR_WTR_SOIL_PROG , 
     I    FR_GLA_SOIL_PROG, 
     I    FACT_WGT_H2O   , FACT_WGT_H2O_INV , 
     I    CAP_SOLID      , 
     I    DELT_SOIL      , 
     I    WTR_BUCKET_INTO_SOIL    , ENH_BUCKET_INTO_SOIL    , 
     I    WTR_BUCKET_SNOW_TO_SNSL , ENH_BUCKET_SNOW_TO_SNSL , !雪なし域のみ
     I    H2O_BUCKET_FROM_SOIL    , ENH_BUCKET_FROM_SOIL    , 
     I    TRNS_BUCKET_NOS  , 
     O    TMP_SOIL       , TMP_SOIL4_N   ,
     O    FR_WTR_SOIL    , FR_GLA_SOIL   , 
     O    CAP_SOIL       , CAP_SOIL_INV  , 
     O    ENG_SOIL_NOS   , H2O_SOIL_NOS  , 
     O    PHASE_SOIL     , RUNOFF_S_NOS      ) 
!
! ------------
! > 熱の計算 <  インプリシットの場合
! ------------
!
      BARESOIL_IMPLICIT: IF ( JCN_SOIL_HEAT_EXPLICIT.EQ.0 ) THEN
!
!   熱の方程式を立てる 
!
        CALL SOIL_HFLUX_IMPLICIT_EQSET (
     I  1                , NMAX_01(3)          , 
#ifdef MONYOS
     I  0                ,
#endif
     I  CVR_NOS_INV      , 
     I  TMP_SOIL         , TMP_SOIL4_N         , 
     I  PHASE_SOIL       , 
     I  CAP_SOIL         , 
     I  HDF_SOIL_H       , 
     I  ENG_BUCKET_NOS   , 
     I  DELT_SOIL        , DELT_SOIL_INV       ,
     O  SIB_MTX_SOIL     , HFLUX_H    )
!
!   熱の方程式を解く 
!
        CALL SOIL_HFLUX_IMPLICIT_EQSLV (
     I  1            , NMAX_01(3)     , 
     M  SIB_MTX_SOIL   )
!
!   熱のフラックスを決定 
!     強制は上端から入る、熱バケツに入っていたエネルギー
!      水と氷の量は変わらないとして、熱収束量を計算する。
!
        CALL SOIL_HFLUX_IMPLICIT_CONV (
     I  1               , NMAX_01(3)      , 
#ifdef MONYOS
     I  0               ,
#endif
     I  SIB_MTX_SOIL    , 
     I  HDF_SOIL_H      , 
     I  PHASE_SOIL      , 
     I  DELT_SOIL       ,  
     M  HFLUX_H         ,  
     O  HEAT_CONV    )  
!
! ------------
! > 熱の計算 <  イクスプリシットの場合
! ------------
!
      ELSE BARESOIL_IMPLICIT
!
        CALL SOIL_HFLUX_EXPLICIT (
     I  1                , NMAX_01(3)          , 
#ifdef MONYOS
     I  0                ,
#endif
     I  CVR_NOS_INV      , 
     I  TMP_SOIL         , TMP_SOIL4_N         ,  
     I  HDF_SOIL_H       , 
     I  ENG_BUCKET_NOS   , 
     I  DELT_SOIL        , DELT_SOIL_INV       ,
     O  HEAT_CONV )

      ENDIF BARESOIL_IMPLICIT
!
! ----------------------------
! > 水のフラックスを解く準備 <
! ----------------------------
!    抜くだけの水や氷が土壌に存在しない場合に措置を施す。
!    合計では存在するならば、相変化を起こさせる。
!    （万一どちらも存在しないなら、仕方がないので外から与える。）    
!
! ------------------------
! > 水のフラックスを解く <
! ------------------------
!    強制は上とのやりとり、蒸散。
!    水の移動フラックスを解く。
!
      CALL SOIL_WFLUX (
     I   1               , NMAX_01(3)      , 
     I   0               ,                       ! INSFLAG = 0 非雪, 1 雪
     I   FR_WTR_SOIL     , 
     I   FR_GLA_SOIL     , 
     I   FACT_WGT_H2O    , 
     I   DELT_SOIL       , 
     I   NTYP            , 
     I   TMP_SOIL        , 
     O   WFLUX_H         , WTR_INC ) 
!
! ------------------------
! > フラックス収束の確定 <
! ------------------------
!  ENG と H2O の積分もしている!! 
!
      CALL SOIL_FLUX_CONV (
     I   1             , NMAX_01(3)    , DELT_SOIL     , 
#ifdef MONYOS
     I   0             ,
#endif
     I   HEAT_CONV     , 
     I   WFLUX_H       , WTR_INC      , 
     I   CAP_SOIL_INV  , 
     I   TMP_SOIL      , 
     M   ENG_SOIL_NOS  , H2O_SOIL_NOS , TMP_SOIL4_N  ,  
     O   RUNOFF_B_NOS) 
!
      ENDIF BARESOIL_GRIDS
!
! ====================
! >> 雪あり域の計算 <<
! ====================
!
      SNOW_GRIDS: IF ( NMAX_0(1) .LT. NMAX ) THEN
!
      CALL SOIL_PRE_2 (
     I    NMAX_01(2)    , NMAX         , 1   , 
     I    0              , 0              ,           ! ダミー
     I    CVR_SNW_INV   ,  NTYP , 
     I    TMP_SNSL_PROG    , TMP_SOIL4_PROG , 
     I    FR_WTR_SNSL_PROG , 
     I    FR_GLA_SNSL_PROG , 
     I    FACT_WGT_H2O   , FACT_WGT_H2O_INV   , CAP_SOLID  , 
     I    DELT_SOIL         ,
     I    WTR_BUCKET_SNOW_TO_SNSL , ENH_BUCKET_SNOW_TO_SNSL , ! 部分積雪域のみ
     I    DUMMY_IDX      , DUMMY_IDX         , 
     I    DUMMY_IDX      , DUMMY_IDX         , 
     I    TRNS_BUCKET_SNW  , 
     O    TMP_SOIL         , TMP_SOIL4_S    ,
     O    FR_WTR_SOIL      , FR_GLA_SOIL    , 
     O    CAP_SOIL         , CAP_SOIL_INV   , 
     O    ENG_SOIL_SNW     , H2O_SOIL_SNW   , 
     O    PHASE_SOIL       , RUNOFF_S_SNW ) 

#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_H2O(
     I  CVR_SNW_INV_NEXT(1:NMAX), FACT_WGT_H2O(1:NMAX, 1),
     I  H2O_SOIL_SNW(1:NMAX, 1), 'soil_pre_2 result SNW')
#endif

! ------------
! > 熱の計算 <  インプリシットの場合
! ------------
!
      SNOW_IMPLICIT: IF ( JCN_SOIL_HEAT_EXPLICIT.EQ.0 ) THEN

!  熱の方程式を立てる 
!
        CALL SOIL_HFLUX_IMPLICIT_EQSET (
     I  NMAX_01(2)       , NMAX                , 
#ifdef MONYOS
     I  1                ,
#endif
     I  CVR_SNW_INV      , 
     I  TMP_SOIL         , TMP_SOIL4_S         ,  
     I  PHASE_SOIL       , 
     I  CAP_SOIL         , 
     I  HDF_SNSL_H       , 
     I  ENG_BUCKET_SNOW_TO_SNSL   , 
     I  DELT_SOIL        , DELT_SOIL_INV       , 
     O  SIB_MTX_SOIL     , HFLUX_H        )
!
!  熱の方程式を解く 
!
        CALL SOIL_HFLUX_IMPLICIT_EQSLV (
     I  NMAX_01(2)       , NMAX            , 
     M  SIB_MTX_SOIL   )
!
!  熱のフラックスを決定 
!
        CALL SOIL_HFLUX_IMPLICIT_CONV (
     I  NMAX_01(2)      , NMAX            , 
#ifdef MONYOS
     I  1               ,
#endif
     I  SIB_MTX_SOIL    , 
     I  HDF_SNSL_H      , 
     I  PHASE_SOIL      , 
     I  DELT_SOIL       , 
     M  HFLUX_H         , 
     O  HEAT_CONV    )  
!
! ------------
! > 熱の計算 <  イクスプリシットの場合
! ------------
!
      ELSE SNOW_IMPLICIT

        CALL SOIL_HFLUX_EXPLICIT (
     I  NMAX_01(2)       , NMAX                , 
#ifdef MONYOS
     I  1                ,
#endif
     I  CVR_SNW_INV      , 
     I  TMP_SOIL         , TMP_SOIL4_S         , 
     I  HDF_SNSL_H       ,
     I  ENG_BUCKET_SNOW_TO_SNSL   , 
     I  DELT_SOIL        , DELT_SOIL_INV       , 
     O  HEAT_CONV   )

      ENDIF SNOW_IMPLICIT
!
! ------------------------
! > 水のフラックスを解く <
! ------------------------
!
      CALL SOIL_WFLUX (
     I   NMAX_01(2)      , NMAX            , 
     I   1               ,                       ! INSFLAG = 0 非雪 , 1 雪 
     I   FR_WTR_SOIL     , 
     I   FR_GLA_SOIL     , 
     I   FACT_WGT_H2O    , 
     I   DELT_SOIL       , 
     I   NTYP            , 
     I   TMP_SOIL        , 
     O   WFLUX_H         , WTR_INC   )

#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_H2O(
     I  CVR_SNW_INV_NEXT(1:NMAX), FACT_WGT_H2O(1:NMAX, 1),
     I  H2O_SOIL_SNW(1:NMAX, 1), 'soil_flux_conv entry SNW')
#endif

! ------------------------
! > フラックス収束の確定 <
! ------------------------
!  ENG と H2O の積分もしている!! 
!
      CALL SOIL_FLUX_CONV (
     I   NMAX_01(2)    , NMAX          , DELT_SOIL     , 
#ifdef MONYOS
     I   1             ,
#endif
     I   HEAT_CONV     , 
     I   WFLUX_H       , WTR_INC       , 
     I   CAP_SOIL_INV  , 
     I   TMP_SOIL      , 
     M   ENG_SOIL_SNW  , H2O_SOIL_SNW  ,  TMP_SOIL4_S  ,  
     O   RUNOFF_B_SNW) 
!

#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_H2O(
     I  CVR_SNW_INV_NEXT(1:NMAX), FACT_WGT_H2O(1:NMAX, 1),
     I  H2O_SOIL_SNW(1:NMAX, 1), 'soil_flux_conv result SNW')
#endif

      ENDIF SNOW_GRIDS
!
!  ============
!  >> 後処理 <<
!  ============
!
!  --------------------------------
!  > 雪あり・無し域の単純平均処理 <
!  --------------------------------
!   必要なのかどうか分からないが、とりあえず作っておく。
!   00/08/16 現在はランオフのみだが、
!   いずれはモニタ用データ作成のために、
!   より多くの変数でなされることになろう。
!
!   土壌第 4 層の温度の積分はここでしたことになっている
!
      CALL SOIL_AREA_MEAN (
     I   NMAX          , 
     I   CVR_SNW       , 
     I   RUNOFF_B_NOS  , RUNOFF_B_SNW  , 
     I   RUNOFF_S_NOS  , RUNOFF_S_SNW  , 
     I   TMP_SOIL4_N   , TMP_SOIL4_S   , 
     O   RUNOFF        , TMP_SOIL4_PROG ) 
!
! --------------------------
! > CVR_SNW 変更に伴う補正 <
! --------------------------
!
!   ・雪ありなし域共通である第 4 層についても、
!     詳細モニタでは、雪被覆率に伴う雪無し域と雪あり域との熱のやりとり
!     をモニタ
!
!   ・第 1-3 層については、*_SNW と *_NOS から *_NEXT を作成
!     このとき雪被覆率の変化を考慮する。
!     例えば被覆率増加時は、雪無し域の土壌の熱・水の一部を雪あり域に加算 
!     詳細モニタでも、このやりとりをモニタする
!
      CALL SOIL_SNOWCVR_MOD (
     I   NMAX             , 
     I   CVR_SNW          , CVR_SNW_NEXT     , 
     I   CVR_NOS_INV_NEXT , CVR_SNW_INV_NEXT , 
     I   H2O_SOIL_SNW     , H2O_SOIL_NOS     ,
     I   ENG_SOIL_SNW     , ENG_SOIL_NOS     , 
     I   TMP_SOIL4_PROG   , 
     O   H2O_SOIL_NEXT    , H2O_SNSL_NEXT    ,
     O   ENG_SOIL_NEXT    , ENG_SNSL_NEXT     )

#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_H2O(
     I  CVR_NOS_INV_NEXT(1:NMAX), FACT_WGT_H2O(1:NMAX, 1),
     I  H2O_SOIL_NEXT(1:NMAX, 1), 'soil_snowcvr_mod result NOS')
      CALL TEST_H2O(
     I  CVR_SNW_INV_NEXT(1:NMAX), FACT_WGT_H2O(1:NMAX, 1),
     I  H2O_SNSL_NEXT(1:NMAX), 'soil_snowcvr_mod result SNW')
#endif
!
! ----------------------
! > 土壌予報変数の計算 <
! ----------------------
!   ・熱と H2O の未来値から、土壌温度と土壌液水・氷の割合の未来値を算出 
!
      CALL SOIL_INTG (
     I   NMAX               , 
     I   FACT_WGT_H2O_INV   , CAP_SOLID         ,
     I   INFO_SNOW_NEXT     , 
     I   H2O_SOIL_NEXT      , H2O_SNSL_NEXT     ,
     I   ENG_SOIL_NEXT      , ENG_SNSL_NEXT     , 
     O   TMP_SOIL_PROG      , TMP_SNSL_PROG     , 
     O   FR_WTR_SOIL_PROG   , FR_WTR_SNSL_PROG  , 
     O   FR_GLA_SOIL_PROG   , FR_GLA_SNSL_PROG ) 

#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_FRW(
     I  CVR_NOS_INV_NEXT(1:NMAX), FR_WTR_SOIL_NEXT(1:NMAX, 1),
     I  FR_GLA_SOIL_NEXT(1:NMAX, 1), 'soil_intg result NOS')
      CALL TEST_FRW(
     I  CVR_SNW_INV_NEXT(1:NMAX), FR_WTR_SNSL_NEXT(1:NMAX,1),
     I  FR_GLA_SNSL_NEXT(1:NMAX), 'soil_intg result SNW')
#endif

! --------------
! > 土壌の補正 <
! --------------
!   ・土壌水が少なすぎるとき、氷が多すぎるときに補正する。 
!     00/08/16 現在、系外とやりとりする（やりとり量をモニタ）
!
!   ・温度、氷/水比は変わらない
!     # このようにしたいので、soil_intg を先に行なっている。
!
!   ・CVR_SNW の変更に伴う補正はここでは行わない。( soil_snowcvr_mod でする ) 
!
      CALL SOIL_MOD (
     I   NMAX             , 
     I   FACT_WGT_H2O     , INFO_SNOW_NEXT   , 
     M   TMP_SOIL_PROG    , TMP_SNSL_PROG    , 
     M   FR_WTR_SOIL_PROG , FR_WTR_SNSL_PROG , 
     M   FR_GLA_SOIL_PROG , FR_GLA_SNSL_PROG ) 
!
#ifdef SIB_DEBUG
      DO K=1,IDP
      DO I=1,NMAX
        IF (       TMP_SOIL_PROG(I,K) .LT. 150 
     1       .AND. INFO_SNOW_NEXT(I)  .LE. 0   ) THEN
          write (6,*) 'SOIL_MAIN ERROR 1' , I,K,TMP_SOIL_PROG(I,K), 
     1                 INFO_SNOW_NEXT(I),NMAX_0(3),NMAX_01(3) 
          stop 999
        ENDIF
      ENDDO
!
      DO I=1,NMAX
        IF (       TMP_SNSL_PROG(I,K) .LT. 150 
     1       .AND. INFO_SNOW_NEXT(I)  .GE. 0   ) THEN
          write (6,*) 'SOIL_MAIN ERROR 2' , I,K,TMP_SNSL_PROG(I,K) ,
     1                 INFO_SNOW_NEXT(I),NMAX_0(3),NMAX_01(3) 
          stop 999
        ENDIF
      ENDDO
      ENDDO
#endif


#ifdef DEBUG_PARANOIA
      ! 第1層の水溢れをチェック
      CALL TEST_FRW(
     I  CVR_NOS_INV_NEXT(1:NMAX), FR_WTR_SOIL_PROG(1:NMAX, 1),
     I  FR_GLA_SOIL_PROG(1:NMAX, 1), 'soil_main ENDING NOS')
      CALL TEST_FRW(
     I  CVR_SNW_INV_NEXT(1:NMAX), FR_WTR_SNSL_PROG(1:NMAX),
     I  FR_GLA_SNSL_PROG(1:NMAX), 'soil_main ENDING SNW')
#endif

      RETURN

#ifdef DEBUG_PARANOIA
      CONTAINS

      ! 積雪被覆率の逆数の変数は CVR_SNW が 1 または 0 に近いときは
      ! 無限大に発散しないように上限がつけられており、Inf 評価はしない。
      SUBROUTINE TEST_FRW(CVRSNW_INV, FRW, FRG, MSG)
        REAL(8), INTENT(IN):: CVRSNW_INV(:)
        REAL(8), INTENT(IN):: FRW(:)
        REAL(8), INTENT(IN):: FRG(:)
        CHARACTER(LEN = *), INTENT(IN):: MSG
        REAL(8):: TEST(SIZE(FRW))
        INTEGER:: I
      CONTINUE
        TEST(:) = CVRSNW_INV(:) * (FRW(:) + FRG(:))
        IF (ANY(TEST(:) > 1.0d0)) THEN
          WRITE(6, *) '#' // MSG // ': FRW+FRG > 1 FOR <',
     ,      PACK((/(I, I = 1, SIZE(TEST))/), (TEST(:) > 1.0d0)),
     ,      '> MAX AT', MAXLOC(TEST(:))
          WRITE(6, *) ' MAX(FRW + FRG) - 1 =', MAXVAL(TEST(:)) - 1.0d0
        ENDIF
      END SUBROUTINE TEST_FRW

      SUBROUTINE TEST_H2O(CVRSNW_INV, FACT_WGT_H2O, H2O, MSG)
        REAL(8), INTENT(IN):: CVRSNW_INV(:)
        REAL(8), INTENT(IN):: FACT_WGT_H2O(:)
        REAL(8), INTENT(IN):: H2O(:)
        CHARACTER(LEN = *), INTENT(IN):: MSG
        REAL(8):: TEST(SIZE(H2O))
        INTEGER:: I
      CONTINUE
        TEST(:) = CVRSNW_INV(:) * H2O(:) - FACT_WGT_H2O(:)
        IF (ANY(TEST(:) > 0.0d0)) THEN
          WRITE(6, *) '#' // MSG // ': H2O > SOIL CAPACITY'         
          WRITE(6, *) ' FOR', PACK((/(I, I = 1, SIZE(TEST))/), 
     ,      (TEST(:) > 0.0d0))
          WRITE(6, *) ' MAX EXCESS=', MAXVAL(TEST(:)), ' AT',
     ,      MAXLOC(TEST(:))
        ENDIF
      END SUBROUTINE TEST_H2O
#endif
! 20050517
! /* DEBUG_PARANOIA */

      END SUBROUTINE SOIL_MAIN

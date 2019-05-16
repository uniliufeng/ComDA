! キャノピー計算メイン
! 99/07/xx 保坂征宏      
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/xx 保坂征宏 GMT 等を変更
!                   フラックスや収束・発散の向き決定
! 99/09/02 保坂征宏 HUM と E を区別する変更。 FAC_EVP_2 を導入
!                   CANOPY_IR を CANOPY_RAD に名称変更 
! 99/09/16 保坂征宏 バグ固定。EVAP_TOT と SUBL_TOT の配列
! 99/09/17 保坂征宏 モニタ追加。バケツ関連
! 99/09/28 保坂征宏 canopy_get_ea 分割関連
!
! 99/11/27 保坂征宏 仮積分ができるように変更。
!                   ただし見掛けは イテレーションの如くにした。
      module canopy_main 
!
      contains
!
      SUBROUTINE CANOPY_MAIN__run (
     I    NMAX_0        , NMAX_01       , NMAX_1     , NMAX        , 
     I    NTYP          , MON           , IMSKLT      , 
     I    CVR_SNW       , CVR_NOS       , HDF_GRSK_IN  , HDF_SNSK_IN   , 
     I    DELT_ATM      , DELT_CNP      , 
     I    P_SFC         , U_REF         , V_REF     , Z_REF     , 
     I    Q_REF         , P_DEL     , FAC_CORTMP    , 
     I    GLA_SNOW      , WTR_SNOW      , 
     I    TMP_SNOW      , TMP_SOIL      , TMP_SNSL  , TMP_REF_IN , 
     I    FR_WTR_SOIL   , FR_GLA_SOIL , FR_WTR_SNSL  , FR_GLA_SNSL  , 
     I    CVR_VEG_ms    , CVR_VEG_INV_ms, THERMK_NOS , THERMK_SNW   , 
     M    R_STM_BAR_PAR_NOS , R_STM_BAR_PAR_SNW ,   
     I    DEP_SNOW      , 
     I    RAD_SHRT_NOS  , RAD_SHRT_SNW  , RAD_LONG_NOS  , RAD_LONG_SNW, 
     I    RAD_SHRT_INTO_SNOW_LAYER      ,RAIN_LARGE    , RAIN_CONV    , 
     I    CAP_SOIL_TOP    , CAP_SNOW_TOP     , 
     M    TMP_CNP_NOS_PROG, FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     M    TMP_CNP_SNW_PROG, FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG, 
     M    TMP_GRSK_PROG   , FR_WTR_GRS_PROG    , FR_GLA_GRS_PROG    , 
     M    TMP_SNSK_PROG   , GMT          , GMQ           , GMUV     , 
     O    RAD_LONG_SIBMIX_2_ATM ,RAD_LONG_NOS_2_ATM,RAD_LONG_SNW_2_ATM,
     M    ENH_BUCKET_INTO_SOIL  , ENH_BUCKET_INTO_SNOW , 
     M    ENH_BUCKET_FROM_SOIL  , ENH_BUCKET_FROM_SNOW , 
     M    H2O_BUCKET_INTO_SOIL  , H2O_BUCKET_INTO_SNOW , 
     M    H2O_BUCKET_FROM_SOIL  , H2O_BUCKET_FROM_SNOW , 
     M  TRNS_BUCKET_NOS,TRNS_BUCKET_SNW,ENG_BUCKET_NOS,ENG_BUCKET_SNW, 
     M    SRAD_BUCKET_SNW   ,TMP_SOIL1_TMP       , TMP_SNOW1_TMP    ,
     O    US , VS , TS , QS  , U_FLUX_OUT , V_FLUX_OUT,z0_out ,
     O    flsh_out , fllh_out , qflux_out )
!
!     use com_debug
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IWG    ,
     1    ICG    , 
     1    IEST   , 
     1    ICGB   , 
     1    ICGR   , 
     1    ITYP   , 
     1    IDP    , 
     1    ISN 
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_ITR_MAX  ,
!    1    JCN_WAY_OF_INTEGRATION     ,      ! 1 イテレーション , 2 仮積分
     1    RCN_FACT_SS1  
!
      USE SIBCON , ONLY :     
     1    TMP_FREZ_c ,
     1     EPSFAC_c , delt_fac_cnp_c , delt_fac_atm_c,
     1     stefan_c , hl_gla_c , hl_wtr_c 
     1   , fact_tskin_top_c
!
!     use com_teten , only : 
      use com_teten_sib_0109 , only : 
     1   table , dtable 
!
      use canopy_integration_overimp, only : 
     1   canopy_integration_overimp_itr ,
     1   canopy_integration_overimp_run 
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX_0 (3)         ! 雪無し格子情報
      INTEGER,INTENT(IN)   :: NMAX_01(3)         ! 雪あり無し共存格子情報
      INTEGER,INTENT(IN)   :: NMAX_1 (3)         ! 
      INTEGER,INTENT(IN)   :: NMAX               ! 陸格子個数
!
      REAL(8),INTENT(IN)   :: HDF_GRSK_IN (IDX,2)   ! 地面表面での熱拡散係数
      REAL(8),INTENT(IN)   :: HDF_SNSK_IN (IDX,2)   ! 雪  表面での熱拡散係数
      REAL(8),INTENT(IN)   :: DELT_ATM           ! 大気計算用 DELT 
      REAL(8),INTENT(IN)   :: DELT_CNP           ! キャノピー回り計算用 DELT 
      REAL(8),INTENT(IN)   :: TMP_SNOW (IDX,ISN+1) ! 雪の温度
      REAL(8),INTENT(IN)   :: WTR_SNOW (IDX)     ! 雪の第一層水量
      REAL(8),INTENT(IN)   :: GLA_SNOW (IDX)     ! 雪の第一層氷量
      REAL(8),INTENT(IN)   :: TMP_SOIL (IDX,IDP) ! 地面第一層温度
      REAL(8),INTENT(IN)   :: TMP_SNSL (IDX,IDP) ! 地面第一層温度
!
      REAL(8),INTENT(IN)   :: FR_WTR_SOIL(IDX,IDP)   ! 空隙中の水の割合
      REAL(8),INTENT(IN)   :: FR_GLA_SOIL(IDX,IDP)   ! 空隙中の氷の割合
      REAL(8),INTENT(IN)   :: FR_WTR_SNSL(IDX,IDP)   ! 空隙中の水の割合
      REAL(8),INTENT(IN)   :: FR_GLA_SNSL(IDX,IDP)   ! 空隙中の氷の割合
!
      REAL(8),INTENT(IN)   :: CVR_SNW  (IDX)     
      REAL(8),INTENT(IN)   :: CVR_NOS  (IDX)     
!
      REAL(8),INTENT(IN)   :: CVR_VEG_ms   (IDX,ICG) ! 植生被覆率
      REAL(8),INTENT(IN)   :: CVR_VEG_INV_ms (IDX,ICG) ! 植生被覆率
!
      REAL(8),INTENT(IN)   :: DEP_SNOW(IDX,ISN+1)! 積雪深( m ) 
!
      REAL(8),INTENT(IN)   :: THERMK_NOS   (IDX)    ! 
      REAL(8),INTENT(IN)   :: THERMK_SNW   (IDX)    ! 
!
      REAL(8),INTENT(IN)   :: RAD_SHRT_NOS  (IDX,ICG)  ! 
      REAL(8),INTENT(IN)   :: RAD_SHRT_SNW  (IDX,ICG)  ! 
      REAL(8),INTENT(IN)   :: RAD_LONG_NOS  (IDX,ICG)
      REAL(8),INTENT(IN)   :: RAD_LONG_SNW  (IDX,ICG)
!
      REAL(8),INTENT(IN)   :: RAD_SHRT_INTO_SNOW_LAYER (IDX) 
!
      REAL(8),INTENT(IN)   :: RAIN_LARGE(IDX) 
      REAL(8),INTENT(IN)   :: RAIN_CONV (IDX) 
!
      INTEGER,INTENT(IN)   :: NTYP  (IDX)          ! 植生タイプ
      INTEGER,INTENT(IN)   :: MON 
      INTEGER,INTENT(IN)   :: IMSKLT(IDX)
!
      REAL(8),INTENT(IN)   :: FAC_CORTMP   (IDX)     ! BPS or κ 
!                                                    ! κT_1 = 地面付近でのT 
!
!
      REAL(8),INTENT(INOUT):: R_STM_BAR_PAR_SNW (IDX,ICG)          ! 気孔抵抗
      REAL(8),INTENT(INOUT):: R_STM_BAR_PAR_NOS (IDX,ICG)          ! 気孔抵抗
!
      REAL(8),INTENT(IN)   :: TMP_REF_IN(IDX)   ! 参照温度
      REAL(8),INTENT(IN)   :: P_SFC     (IDX)   ! 
      REAL(8),INTENT(IN)   :: Z_REF     (IDX)  
      REAL(8),INTENT(IN)   :: U_REF     (IDX)  
      REAL(8),INTENT(IN)   :: V_REF     (IDX)  
      REAL(8),INTENT(IN)   :: Q_REF     (IDX)   ! 大気最下層比湿（参照比湿）
      REAL(8),INTENT(IN)   :: P_DEL     (IDX)   ! 大気最下層圧 HPa 
!
!  SiB 予報変数 （入力）
!
      REAL(8),INTENT(INOUT) :: TMP_CNP_NOS_PROG   (IDX)  ! キャノピー温度
      REAL(8),INTENT(INOUT) :: TMP_CNP_SNW_PROG   (IDX)  ! キャノピー温度
      REAL(8),INTENT(INOUT) :: FR_WTR_CNP_SNW_PROG(IDX)  ! キャノピー水割合
      REAL(8),INTENT(INOUT) :: FR_GLA_CNP_SNW_PROG(IDX)  ! キャノピー氷割合
      REAL(8),INTENT(INOUT) :: FR_WTR_CNP_NOS_PROG(IDX)  ! キャノピー水割合
      REAL(8),INTENT(INOUT) :: FR_GLA_CNP_NOS_PROG(IDX)  ! キャノピー氷割合
      REAL(8),INTENT(INOUT) :: TMP_GRSK_PROG      (IDX)  ! 地面・下草温度
      REAL(8),INTENT(INOUT) :: FR_WTR_GRS_PROG    (IDX)  ! 下草水割合
      REAL(8),INTENT(INOUT) :: FR_GLA_GRS_PROG    (IDX)  ! 下草氷割合
      REAL(8),INTENT(INOUT) :: TMP_SNSK_PROG      (IDX)  ! 雪スキン温度
!
      REAL(8),INTENT(INOUT):: GMT      (IDX,2)      
      REAL(8),INTENT(INOUT):: GMQ      (IDX,2)
      REAL(8),INTENT(INOUT):: GMUV     (IDX,3)

!
      REAL(8),INTENT(IN) :: CAP_SOIL_TOP      (IDX)  
      REAL(8),INTENT(IN) :: CAP_SNOW_TOP      (IDX)  
!
! OUTPUT
!
      REAL(8),INTENT(INOUT) :: RAD_LONG_SNW_2_ATM (IDX)
      REAL(8),INTENT(INOUT) :: RAD_LONG_NOS_2_ATM (IDX)
      REAL(8),INTENT(OUT)   :: RAD_LONG_SIBMIX_2_ATM (IDX)  ! 上向き正!!! 
!
!     REAL(8),INTENT(OUT)   :: U_FLUX(IDX)           ! 運動量フラックス
!     REAL(8),INTENT(OUT)   :: V_FLUX(IDX)
!

!
! INOUT  
!
!   バケツ... 基本的に格子平均値とする。( *CVR_SNW などをする ) 
!
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_INTO_SOIL (IDX) 
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_FROM_SOIL (IDX) 
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_INTO_SNOW (IDX) 
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_FROM_SNOW (IDX) 
!
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO_SOIL (IDX) 
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_FROM_SOIL (IDX) 
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO_SNOW (IDX) 
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_FROM_SNOW (IDX) 
!
      REAL(8),INTENT(INOUT) :: ENG_BUCKET_NOS     (IDX) 
      REAL(8),INTENT(INOUT) :: ENG_BUCKET_SNW     (IDX) 
!
      REAL(8),INTENT(INOUT) :: TRNS_BUCKET_NOS    (IDX,IDP) 
      REAL(8),INTENT(INOUT) :: TRNS_BUCKET_SNW    (IDX,IDP) 
      REAL(8),INTENT(INOUT) :: SRAD_BUCKET_SNW    (IDX) 
!
      REAL(8),INTENT(INOUT) :: TMP_SOIL1_TMP      (IDX) 
      REAL(8),INTENT(INOUT) :: TMP_SNOW1_TMP      (IDX) 
!
      real(8),intent(out) :: ts(idx)
      real(8),intent(out) :: qs(idx)
      real(8),intent(out) :: us(idx)
      real(8),intent(out) :: vs(idx)
      real(8),intent(out) :: u_flux_out(idx)
      real(8),intent(out) :: V_flux_out(idx)
!
      real(8),intent(out) :: z0_out   (idx)
      real(8),intent(out) :: flsh_out (idx) 
      real(8),intent(out) :: fllh_out (idx) 
      real(8),intent(out) :: qflux_out(idx)
!
      real(8)   :: ta(IDX) 
      real(8)   :: ea(IDX) 
      real(8)   :: tg(IDX) 
      real(8)   :: fllh (IDX) 
      real(8)   :: flsh (IDX) 
      real(8)   :: qflux(IDX) 
      real(8)   :: u2o  (idx)
!
! WORK  ～ 仮積分・イテレーション用変数    
!
!   イテレーションカウンタ等
!
      INTEGER   :: ITR                  
!
!   インターセプション後の値
!
      REAL(8)   :: TMP_CNP_SAVE    ( IDX )   
      REAL(8)   :: TMP_SKIN_SAVE   ( IDX ) 
      REAL(8)   :: TMP_REF         ( IDX )   
      REAL(8)   :: E_REF_SAVE      ( IDX )  
!     REAL(8)   :: FR_WTR_CNP_SAVE ( IDX )
!     REAL(8)   :: FR_GLA_CNP_SAVE ( IDX )
!     REAL(8)   :: FR_WTR_GRS_SAVE ( IDX )
!     REAL(8)   :: FR_GLA_GRS_SAVE ( IDX )
!
!
! WORK    
!
      integer   :: insflag
      integer   :: nstr 
      integer   :: nend
!
!     INTEGER   :: KLEV_TRNS_MAX(IDX,ICG)
!     REAL(8)   :: SRAD_INTO_SNOW_LAYER_TMP  (IDX)  ! 雪スキンを抜ける短波
!
!
      REAL(8)   :: dta_dtc (IDX) 
      REAL(8)   :: dta_dtg (IDX) 
      REAL(8)   :: dta_dtr (IDX) 
!
      INTEGER   :: I
      INTEGER   :: ITC
      INTEGER   :: ITC_MAX 
!
      REAL(8)  :: RHO_AIR    (IDX)     ! 
!
      REAL(8)  :: TMP_CNP    (IDX)     ! 
      REAL(8)  :: TMP_SKIN   (IDX)     ! 
!
      REAL(8)  :: HDF_SKIN   (IDX,2)     ! 
!
      REAL(8)  :: FLAG_DRY_CNP         (IDX)     ! 
      REAL(8)  :: FLAG_DRY_SKIN        (IDX)     ! 
!
!
      REAL(8)  :: FAC_AC_BAR_INV_ms(IDX)            ! 旧 RBC 
      REAL(8)  :: FAC_AG_INV_ms    (IDX)            ! 旧 RDC 
      REAL(8)  :: FAC_AC_BAR_INV_msx(IDX)           ! 旧 RBC 
      REAL(8)  :: FAC_AG_INV_msx    (IDX)           ! 旧 RDC 
!
      REAL(8)  :: XLAI_VEG_ms    (IDX,ICG)   !  
      REAL(8)  :: XLAI_VEG_msx   (IDX,ICG)   ! 本当は Canopy だけでよい
!
      REAL(8)  :: RHO_CP     (IDX)           ! 旧 RCP 
!
      REAL(8)  :: RIB        (IDX)           ! リチャードソン数
!      
      REAL(8)  :: FUNC_TMP(IDX,ICG)          ! 気孔抵抗の温度依存性
!
      REAL(8)  :: D_ZERO_ms  (IDX)           ! 零面変位
      REAL(8)  :: D_ZERO_msx (IDX)      
!
      REAL(8),PARAMETER :: DUMMY_IDX  (IDX,2) = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX1 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX2 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX3 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX4 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX5 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX6 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX7 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX8 (IDX)   = 0.   ! ダミー用ワーク
      REAL(8)           :: DUMMY_XXX9 (IDX)   = 0.   ! ダミー用ワーク
!
      REAL(8)  :: CVR_VEG_msx     (IDX,ICG)  ! 植生被覆率、雪あり域 G で 0
      REAL(8)  :: CVR_VEG_INV_msx (IDX,ICG)  ! 植生被覆率
!
!  蒸散
      REAL(8)  :: FLAG_TR    (IDX,IDP,ICG) ! 蒸散するしないフラグ
      REAL(8)  :: FLAG_TR_ALL(IDX,    ICG) ! 蒸散するしないフラグ
!                                          ! 整合チェック用にメインまで戻す
      REAL(8)  :: FAC_TRNS_A (IDX,IDP,ICG)   ! 土壌層別割合
      REAL(8)  :: FAC_TRNS_B (IDX,IDP,ICG)   ! 土壌層別割合
!
! フラックス関連諸量
!
!   陽的評価  
!
!     REAL(8)  :: HDF        (IDX)           ! 地面表面での熱拡散
      REAL(8)  :: HDF        (IDX,2)                ! XXX1 
      REAL(8)  :: SENS       (IDX,ICG)   
      REAL(8)  :: EVAP       (IDX,ICGB)             !  990916 fix 
      REAL(8)  :: SUBL       (IDX,ICGB)             !  990916 fix
      REAL(8)  :: TRNS       (IDX,ICG)    
!
!   微分係数（温度依存性）
!
      REAL(8)  :: DSENS_DX   (IDX,ICG,ICGR)   
      REAL(8)  :: DRAD_DX    (IDX,ICGR,ICG) 
      REAL(8)  :: DEVAP_DX   (IDX,IEST,ICGB,ICGR) 
!
!
!   方程式と解
!
!     REAL(8)   :: SIB_MTX    (IDX,4,5) 
      REAL(8)   :: SIB_MTX    (IDX,5,6)                  ! XXX1
!
!   蒸散
!
      REAL(8)   :: TRNS_TMP   (IDX,ICG)
!
!   土壌 または 積雪第一層
!
      REAL(8)   :: TMP_SS1      (IDX) 
      REAL(8)   :: FLAG_PHASE_SS1  (IDX)  ! SS1 の温度変化あり 1 無し 0
!                                         ! TMP_SS1 が TMP_FREZ_C なら変化無し
!
!  --------------------------------------------------------
!
!
      REAL(8)   :: H2O_SAT_CNP  (IDX)     ! キャノピー飽和 H2O  単位 kg/m**2
      REAL(8)   :: H2O_SAT_GRS  (IDX)     ! 地面飽和       H2O 
!
      REAL(8)   :: FAC_EVP      (IDX)     ! ρε/P_s
      REAL(8)   :: FAC_EVP_2    (IDX)     ! P_s/ ε
!
! 蒸発等のフラックス評価用変数 （ 結露/結霜時は 1 ) 
!
      REAL(8)   :: XFR_WTR_CNP(IDX)  ! キャ葉   水 被覆率(結露考慮)
      REAL(8)   :: XFR_GLA_CNP(IDX)  ! キャ葉   氷 被覆率(結露考慮)
      REAL(8)   :: XFR_WTR_GRS(IDX)  ! 下草葉   水 被覆率(結露考慮)
      REAL(8)   :: XFR_GLA_GRS(IDX)  ! 下草葉   氷 被覆率(結露考慮)
      REAL(8)   :: XFR_WTR_BS (IDX)  ! 裸地・雪 水 被覆率(結露考慮)
      REAL(8)   :: XFR_GLA_BS (IDX)  ! 裸地・雪 氷 被覆率(結露考慮)
!
! ループ用変数
!
      REAL(8)   :: FR_H2O_CNP(IDX)  ! 
      REAL(8)   :: FR_H2O_GRS(IDX)  ! 
!
      REAL(8)   :: FR_WBH_CNP(IDX)  ! 
      REAL(8)   :: FR_WBH_GRS(IDX)  ! 
      REAL(8)   :: FR_WBH_SNW(IDX)  ! 
      REAL(8)   :: FR_WBH_BS (IDX)  !                            ! bug 必要!!
!
! 抵抗
!
      REAL(8)   :: R_AIR_INV       (IDX)      ! 1/r_air   (=1/ra)
      REAL(8)   :: R_AIR_CONST_INV (IDX)! 1/r_air のうちの定数部分
      REAL(8)   :: R_AC_BAR        (IDX)   ! 1/r_ac    (=1/rb)
      REAL(8)   :: R_AC_BAR_INV    (IDX)   ! 1/r_ac    (=1/rb)
!     REAL(8)   :: R_CC_INV        (IDX)   ! 1/r_cc    (=1/rc) 
      REAL(8)   :: R_ACX_BAR       (IDX)   ! 1/(r_ac+r_cst) 
      REAL(8)   :: R_ACX_BAR_INV   (IDX)   ! 1/(r_ac+r_cst) 
      REAL(8)   :: R_AG            (IDX)   ! 1/r_d     (=1/rd) 
      REAL(8)   :: R_AG_INV        (IDX)   ! 1/r_d     (=1/rd) 
!     REAL(8)   :: R_BB_INV        (IDX)   ! 1/R_BB    (=1/rg)   
      REAL(8)   :: R_AGX_BAR       (IDX)   ! 1/(r_ag+r_gst)
      REAL(8)   :: R_AGX_BAR_INV   (IDX)   ! 1/(r_ag+r_gst)
      REAL(8)   :: R_AB_INV        (IDX)   ! 1/(r_ag+r_bs)
      REAL(8)   :: R_STM_BAR       (IDX,ICG)  !
!
      REAL(8)   :: R_ROOT_BAR_INV (IDX,IDP,ICG)
      REAL(8)   :: R_BB           (IDX)   ! 雪無土壌表層抵抗 
!                                         ! 雪上は扱わない
      REAL(8)   :: POT_SOIL       (IDX,IDP)    ! 土壌のmatric pot

!
      REAL(8)    :: R_SUM_INV (IDX)     ! 
!
      REAL(8)    :: E_REF (IDX)         ! 大気最下層水蒸気圧
!                                       ! e ～ Ps x q1 / ε の式で評価
!                                       ! Ps を使うので地表付近の値に相当?
!
      REAL(8)    :: HUM_REL_GRSK   (IDX)   ! 雪無し土壌第一層空気相対湿度
      REAL(8)    :: HUM_REL_SNSK   (IDX)   ! 雪の第一層の相対湿度 = 1
!
      REAL(8)    :: E_AIR     (IDX)      ! 
      REAL(8)    :: E_AIR_SAT (IDX)      ! キャ空間温度飽和比湿
      REAL(8)    :: E_CNP_SAT (IDX)      ! キャ温度飽和比湿
      REAL(8)    :: E_GRSK_SAT(IDX)      ! 地面温度飽和比湿
      REAL(8)    :: E_SNSK_SAT(IDX)      ! 地面温度飽和比湿
!
      REAL(8)    :: DEC_SAT_DTC  (IDX)  ! キャ温度飽和比湿 微分
      REAL(8)    :: DEG_SAT_DTG  (IDX)  ! 地面温度飽和比湿 微分
!
      REAL(8)   :: TMP_AIR     (IDX)  
      REAL(8)   :: U2    (IDX) 
      REAL(8)   :: RIB_CONST   (IDX)    ! Richardson 数の定数部分 
      REAL(8)   :: Z0_msx      (IDX)  
      REAL(8)   :: Z0_ms       (IDX)  
      REAL(8)   :: A_COEF      (IDX)  
!
      REAL(8)   :: GMUV_SNW    (IDX,2)  
      REAL(8)   :: GMUV_NOS    (IDX,2)  
!
      REAL(8)   :: GMT_SNW     (IDX)  
      REAL(8)   :: GMT_NOS     (IDX)  
      REAL(8)   :: GMQ_SNW     (IDX)  
      REAL(8)   :: GMQ_NOS     (IDX)  
!
      REAL(8)   :: U_FLUX    (IDX)
      REAL(8)   :: V_FLUX    (IDX)
!
      REAL(8),SAVE :: ZERO_IDX    (IDX)       ! ゼロの入った配列
!
!
      INTEGER   :: ICOUNT
      INTEGER   :: INIT_EA
!
      REAL(8)   :: YI
      REAL(8)   :: X
      INTEGER   :: IY
!
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!     DATA ITC_MAX / 5 /
!     DATA ITC_MAX / 2 /
      DATA ITC_MAX / 1 /
!
      ta (:) = 0.
      ea (:) = 0.
      tg (:) = 0.
      u2o(:) = 0.
      z0_out(:) = 0.
      fllh_out(:) = 0.
      flsh_out(:) = 0.
      qflux_out(:) = 0.
!
!     ISSS = ISSS + 1 
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'CANOPY_MAIN 99/07 POCHAKA'
        WRITE(6,*) 'INTEGRATION ITERATION MAX = ' , JCN_ITR_MAX  
        WRITE(6,*) 'EA TA ITERATION MAX = ' , ITC_MAX , 
     1             ' .... DUMMY ' 
!
        DO 100 I=1,NMAX
          ZERO_IDX(I) = 0.
 100    CONTINUE     
!
        LFIRST = .FALSE. 
      ENDIF
!
      TS(:) = 0.
      QS(:) = 0.
      US(:) = 0.
      VS(:) = 0.
      U_FLUX_OUT(:) = 0.
      V_FLUX_OUT(:) = 0.
!
#ifdef SIB_DEBUG
      do i=1,nmax
        if ( 
     1            rain_large(i) .lt. 0.D0 
     1       .or. rain_conv (i) .lt. 0.D0 ) then
          write(6,*) 'canopy_main rain input error ' , i , 
     1                rain_large(i), rain_conv(i)
          stop 999
        endif
      enddo
#endif
!
      CALL CANOPY_EXTEND (
     I   NMAX      , 
     I   NTYP      , MON        , 
     O   FAC_AG_INV_ms , 
     O   FAC_AC_BAR_INV_ms      , 
     O   XLAI_VEG_ms            , 
     O   D_ZERO_ms              , 
     O   Z0_ms       )
!
!  ------------------------
!  > 参照レベル変数の計算 <
!  ------------------------
!
      CALL CANOPY_ATMOS (
     I   NMAX          , 
     I   P_SFC         , TMP_REF_IN  ,  Q_REF    , 
     O   RHO_AIR       , RHO_CP      ,    
     O   E_REF         , 
     O   FAC_EVP       , FAC_EVP_2 )

!
!====================================================================
!
! ==============
! >> 雪無し域 <<
! ==============
!   方程式の係数・諸項を設定、解く
!
      IF ( NMAX_01(3) .GT. 0 ) THEN
!
        nstr = 1
        nend = nmax_01(3)
        insflag = 0
!
        DO I=NSTR,NEND
          FLAG_PHASE_SS1(I) = 1.D0
          IF (      
!    1              TMP_SOIL1_TMP(I).EQ.TMP_FREZ_C 
!    1         .OR. 
     1              ABS(RCN_FACT_SS1).LT.1.D-5 
     1        ) THEN
            FLAG_PHASE_SS1(I) = 0.D0
          ENDIF
        ENDDO
!
        CALL CANOPY_ROOT (
     I   NSTR            , NEND          , INSFLAG     ,
     I   NTYP            , MON           , 
     I   FR_WTR_SOIL     , 
     I   FR_GLA_SOIL     , 
     I   TMP_SOIL        , 
     O   POT_SOIL        , R_ROOT_BAR_INV, 
     O   HUM_REL_GRSK    , R_BB        )   
!
!  ------------------------------
!  > 雪無しで必要になる量の計算 <
!  ------------------------------
!    裸地での蒸発/昇華比の設定
!
        CALL CANOPY_NOS (
     I    NSTR           , NEND        ,
     I    TMP_SOIL       , FR_WTR_SOIL , FR_GLA_SOIL  , 
     O    FR_WBH_BS       ) 
!
!  ----------------------
!  > インターセプション <
!  ----------------------
!
        CALL CANOPY_INTERC (
     I    NSTR              , NEND              , INSFLAG        , 
     I    NTYP              , 
#ifdef MONYOS
     I    DELT_CNP , 
#endif
     I    RAIN_LARGE        , RAIN_CONV         ,
     I    CVR_VEG_ms        , CVR_NOS           , XLAI_VEG_ms    , 
     I    TMP_REF_IN        , 
     I    TMP_CNP_NOS_PROG  , TMP_GRSK_PROG     , 
     I    FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     I    FR_WTR_GRS_PROG   , FR_GLA_GRS_PROG   , 
     O    TMP_CNP           , TMP_SKIN          , 
     O    FR_H2O_CNP        , FR_WBH_CNP        , 
     O    FR_H2O_GRS        , FR_WBH_GRS        , 
     O    H2O_SAT_CNP       , H2O_SAT_GRS       , 
     M    H2O_BUCKET_INTO_SOIL , ENH_BUCKET_INTO_SOIL  )
!
      CALL CANOPY_COPY ( 
     I    NSTR          , NEND          , INSFLAG   ,  
     I    TMP_CNP       , TMP_SKIN      , TMP_REF_IN, TMP_SOIL1_TMP ,  
     I    E_REF         , 
     O    TMP_CNP_SAVE  , TMP_SKIN_SAVE , TMP_REF   , TMP_SS1 , 
     O    E_REF_SAVE )
!
      CALL CANOPY_CANOPY_INIT_1 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , MON         , 
     I   TMP_REF     , U_REF       , V_REF       , Z_REF      ,   
     I   Z0_ms       , D_ZERO_ms   ,
     I   CVR_VEG_ms  , 
     M   R_STM_BAR_PAR_NOS , 
     O   TMP_AIR     , R_AIR_INV   , R_AIR_CONST_INV  , 
     O   U2          , FAC_CORTMP     , 
     O   A_COEF      , 
     O   R_STM_BAR   , FLAG_DRY_CNP   , FLAG_DRY_SKIN    )
!
!  ------------------------------------
!  > イテレーション（仮積分）ここから <
!  ------------------------------------
!
      ITR = 0 
 8000 CONTINUE
      ITR = ITR + 1  
!
!
!  --------
!  > 準備 <
!  --------
!
! キャノピー空間情報・イテレーション用初期値設定
!
      CALL CANOPY_CANOPY_INIT_2 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , 
     I   TMP_CNP     , TMP_SKIN    , 
     I   U_REF       , V_REF       , Z_REF      ,   
     I   Z0_ms       , D_ZERO_ms   ,
     I   TMP_REF     , FAC_CORTMP , 
     O   RIB_CONST   ,
     O   FUNC_TMP    )  
!
!  ------------------------------------
!  > キャノピー関連イテレーション計算 <
!  ------------------------------------
!
      CALL RESET ( FAC_TRNS_A , IDX*IDP*ICG )
      CALL RESET ( FAC_TRNS_B , IDX*IDP*ICG )
!
      ICOUNT  = 0
      INIT_EA = 1
!
#include "noparallel_fixedform"
      DO ITC = 1,ITC_MAX
!
!  TA , RIB の計算
!
      CALL CANOPY_GET_TA ( 
     I    NSTR        , NEND        , 
     I    NTYP        , MON         , 
     I    TMP_CNP     , TMP_SKIN    , TMP_REF       , 
     I    FAC_AG_INV_ms             , XLAI_VEG_ms   , 
     I    FAC_AC_BAR_INV_ms         , 
     I    U2          , RIB_CONST   , 
     I    FAC_CORTMP  , 
     I    Z_REF       , Z0_ms       , R_AIR_CONST_INV , A_COEF        ,  
     M    ICOUNT      , 
     M    TMP_AIR     , R_AIR_INV   , 
     O    R_AC_BAR    , R_AG        , 
     O    R_AC_BAR_INV, R_AG_INV    , R_SUM_INV     , 
     O    E_AIR_SAT   , E_CNP_SAT   , E_GRSK_SAT    , 
     O    DEC_SAT_DTC , DEG_SAT_DTG , RIB )
!

!
! EA , 気孔抵抗の計算 
!
      CALL CANOPY_GET_EA (
     I    NSTR             , NEND              , INSFLAG ,  
     I    ITR          , JCN_ITR_MAX      , 
     I    NTYP         , 
     I    CVR_VEG_ms   , CVR_VEG_INV_ms, 
     I    R_AC_BAR     , R_AC_BAR_INV , R_AG        , R_AG_INV    , 
     I    R_BB         , R_AIR_INV    , 
     I    R_ROOT_BAR_INV, R_STM_BAR_PAR_NOS, 
     I    E_AIR_SAT    , E_CNP_SAT     , E_REF     , E_GRSK_SAT    ,
     I    HUM_REL_GRSK     , 
     I    TMP_SOIL     , 
     I    FR_H2O_CNP   , FR_H2O_GRS  , 
     I    FUNC_TMP     , 
     I    POT_SOIL     ,
     I    D_ZERO_ms    , 
     I    FAC_EVP      ,
     I    ITC          , ITC_MAX      , 
#ifdef MONYOS
     I   DELT_CNP , 
#endif
     M    ICOUNT       , INIT_EA      ,  
     M    R_STM_BAR    , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN   ,     
     O    E_AIR        , FLAG_TR     ,  FLAG_TR_ALL ,
     O    R_ACX_BAR    , R_AGX_BAR    , R_ACX_BAR_INV , R_AGX_BAR_INV ,
     O    R_AB_INV     , 
!    O    FAC_TRNS_A   , FAC_TRNS_B  , KLEV_TRNS_MAX  )  
     O    FAC_TRNS_A   , FAC_TRNS_B  )
!
      ENDDO
!

!
!  仮積分時の場合、2 度目は値を戻す。
!
      IF ( ITR.EQ.JCN_ITR_MAX) THEN
        CALL CANOPY_COPY ( 
     I    NSTR          , NEND          , INSFLAG   ,  
     I    TMP_CNP_SAVE  , TMP_SKIN_SAVE , TMP_REF_IN, TMP_SOIL1_TMP ,
     I    E_REF_SAVE    ,  
     O    TMP_CNP       , TMP_SKIN      , TMP_REF   , TMP_SS1       ,
     O    E_REF    )
!
        DO  I=NSTR,NEND 
          R_SUM_INV(I) 
     1          = 1 .
     1            / ( R_AG_INV(I) + R_AC_BAR_INV(I)*2 + R_AIR_INV(I) )
!
          TMP_AIR(I) 
     1          =  (   TMP_SKIN(I) * R_AG_INV(I)
     1               + TMP_CNP (I) * R_AC_BAR_INV(I)*2 
     1               + TMP_REF (I) * R_AIR_INV(I)  * FAC_CORTMP(I) )
     1             * R_SUM_INV(I) 
!
          YI = (TMP_AIR(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX ( IY, 1 )
          IY = MIN ( IY, 24999 )
          X = YI - IY
          E_AIR_SAT(I) = ((1.0D0-X)*TABLE(IY)+X*TABLE(IY+1))/EPSFAC_c
! 
          YI = (TMP_CNP(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX( IY, 1 )
          IY = MIN( IY, 24999 )
          X = YI - IY
          E_CNP_SAT   (I)=((1.0D0-X)* TABLE(IY)+X* TABLE(IY+1))/EPSFAC_c
          DEC_SAT_DTC (I)=((1.0D0-X)*DTABLE(IY)+X*DTABLE(IY+1))/EPSFAC_c      
!
          YI = (TMP_SKIN(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX( IY, 1 )
          IY = MIN( IY, 24999 )
          X = YI - IY
          E_GRSK_SAT  (I)=((1.0D0-X)* TABLE(IY)+X* TABLE(IY+1))/EPSFAC_c
          DEG_SAT_DTG (I)=((1.0D0-X)*DTABLE(IY)+X*DTABLE(IY+1))/EPSFAC_c
        ENDDO
!
        CALL CANOPY_AIR_ALL (
     I    NSTR       , NEND       , INSFLAG    , 
     I    CVR_VEG_ms , CVR_VEG_INV_ms, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , FR_H2O_GRS , 
     I    E_CNP_SAT  , E_REF      , E_GRSK_SAT   ,
     I    HUM_REL_GRSK   , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN     , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
!
      ENDIF
!
! イテレーション後の処理
!   結露・蒸発も考慮した、水・氷が覆う割合 XFR_WTR/GLA を計算
!
      CALL CANOPY_SET_AIR  (
     I    NSTR         , NEND         , INSFLAG , 
     I    FR_H2O_CNP   , FR_WBH_CNP   , FR_H2O_GRS  , FR_WBH_GRS  , 
     I    FR_WBH_BS    ,   
     I    FLAG_DRY_CNP , FLAG_DRY_SKIN, 
     O    XFR_WTR_CNP  , XFR_GLA_CNP  , XFR_WTR_GRS , XFR_GLA_GRS , 
     O    XFR_WTR_BS   , XFR_GLA_BS   )
!
!  ------------------
!  > フラックス計算 <
!  ------------------
!
!  長波フラックス
!
      CALL CANOPY_RAD (
     I  NSTR        , NEND        , 
     I  TMP_CNP     , TMP_SKIN    , 
     I  THERMK_NOS  , 
     I  CVR_VEG_ms  , 
!    I  FLAG_TMPMOD_CNP   , FLAG_TMPMOD_SKIN   , 
     O  DRAD_DX )
!
! 顕熱フラックス計算
!
      CALL CANOPY_SENS (
     I    NSTR             , NEND              , 
     I    ITR              , JCN_ITR_MAX           , 
#ifdef MONYOS
     I  DELT_CNP    , INSFLAG       , 
#endif 
     I  R_AIR_INV   , R_AC_BAR_INV  , R_AG_INV  , R_SUM_INV , 
     I  RHO_CP      , FAC_CORTMP    , 
     I  TMP_CNP     , TMP_SKIN      , TMP_AIR     , 
     O  DTA_DTC     , DTA_DTG       , DTA_DTR ,
     O  SENS        , 
     O  DSENS_DX    ) 
!
! 潜熱フラックス設定
!
      CALL CANOPY_EVAP (
     I    NSTR             , NEND              , e_air , 
     I    ITR               , JCN_ITR_MAX           , 
#ifdef MONYOS
     I  DELT_CNP    , INSFLAG        , 
#endif
     I  R_AIR_INV   , 
     I  R_AC_BAR_INV, R_ACX_BAR_INV , 
     I  R_AG_INV    , R_AGX_BAR_INV   , R_AB_INV ,  
     I  E_REF       , E_CNP_SAT    , E_GRSK_SAT    , 
     I  DEC_SAT_DTC , DEG_SAT_DTG , 
     I  CVR_VEG_ms  ,                               ! 被覆率
     I  FAC_EVP     , HUM_REL_GRSK    ,                   ! 雪なら後ろは 1 
     I  FAC_EVP_2   , 
     I  XFR_WTR_CNP , XFR_WTR_GRS  , XFR_WTR_BS   ,   ! 水氷割合（結露考慮）
     I  XFR_GLA_CNP , XFR_GLA_GRS  , XFR_GLA_BS   ,   ! 水氷割合（結露考慮） 
!    I  FLAG_TMPMOD_CNP   , FLAG_TMPMOD_SKIN  ,     
     I  FLAG_DRY_CNP      , FLAG_DRY_SKIN     ,
     I  FLAG_TR_ALL ,
     O  EVAP        , SUBL        , TRNS        , 
     O  DEVAP_DX   )
!
!
! 熱伝導フラックス設定 ... 下向き正
!
      CALL CANOPY_HDF ( 
     I    NSTR      , NEND        , 
     I    ITR       , JCN_ITR_MAX , 
#ifdef MONYOS
     I    INSFLAG   , DELT_CNP    , 
#endif
     I    TMP_SKIN  , TMP_SS1     , TMP_SOIL(1,2)  ,      ! XXX1 
     I    HDF_GRSK_IN  , 
     O    HDF      )                                     ! XXX1 
!
      hdf_skin(1:NMAX_01(3),:) 
     1       = hdf_grsk_in(1:NMAX_01(3),:) * fact_tskin_top_c
!
! ----------
! > 方程式 <
! ----------
!
!  方程式をたてる
!
      CALL CANOPY_FLUX_EQSET (
     I  NSTR         , NEND         , INSFLAG , 
     I  ITR - JCN_ITR_MAX , 
     I  XLAI_VEG_ms  , 
     I  FR_H2O_CNP   , FR_WBH_CNP   , FR_H2O_GRS , FR_WBH_GRS  , 
     I  H2O_SAT_CNP  , H2O_SAT_GRS  , 
     I  SENS         ,                            DSENS_DX  , 
     I  EVAP         , SUBL         , TRNS      , DEVAP_DX  , 
     I  DRAD_DX      , 
     I  RAD_SHRT_NOS , RAD_LONG_NOS , 
     I  HDF          , HDF_skin     , 
     I  DELT_ATM     , DELT_CNP     , 
     I  GMT          , GMQ          , P_DEL     , Q_REF     , 
     I  CAP_SOIL_TOP , FLAG_PHASE_SS1 , 
     O  SIB_MTX   )
!
!! 方程式を解く 
!
      CALL CANOPY_FLUX_EQSLV (
     I    NSTR             , NEND              , 
     M    SIB_MTX    )
!
! --------
! > 積分 <        
! --------
!
      IF ( ITR .LT. JCN_ITR_MAX ) THEN
!
        CALL CANOPY_INTEGRATION_OVERIMP_ITR (
     I      NSTR       , NEND        , 
     I      DELT_CNP   , DELT_ATM    , 
     I      SIB_MTX    , 
     I      FAC_EVP_2  ,
     I      DTA_DTC    , DTA_DTG     , DTA_DTR    , 
     M      TMP_CNP    , TMP_SKIN    , TMP_SS1    , 
     M      TMP_REF    , TMP_AIR     , E_REF  )
!

!
        GOTO 8000
!
      ELSE
!
!
        CALL CANOPY_INTEGRATION_OVERIMP_RUN ( 
     I      NSTR       , NEND        , INSFLAG  , IMSKLT , 
     I      DELT_CNP   , DELT_ATM    , 
     I      SIB_MTX    , 
     I      EVAP       , SUBL        , TRNS     , 
     I      SENS       , 
     I      RAD_LONG_NOS , 
     I      RAD_SHRT_NOS , 
     I      HDF      , 
     I      DEVAP_DX   , DSENS_DX    , DRAD_DX  , HDF_skin , 
     I      H2O_SAT_CNP    , H2O_SAT_GRS    , 
     I      CVR_NOS        , XLAI_VEG_ms    , 
     I      TMP_CNP        , 
     I      FR_H2O_CNP     , FR_WBH_CNP    , 
     I      TMP_SKIN       , 
     I      FR_H2O_GRS     , FR_WBH_GRS    , 
     I      TMP_SOIL       , 
     I      DUMMY_IDX      , 
     I      CAP_SOIL_TOP   , FLAG_PHASE_SS1 , 
     I      GMT            , GMQ            ,  
     O      TMP_CNP_NOS_PROG, FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     O      TMP_GRSK_PROG  , FR_WTR_GRS_PROG, FR_GLA_GRS_PROG, 
     O      TRNS_TMP       , 
     O      GMT_NOS        , GMQ_NOS      , 
     O      FLLH           , FLSH       , qflux , 
     M      RAD_LONG_NOS_2_ATM   , 
     M      ENG_BUCKET_NOS              , 
     M      H2O_BUCKET_INTO_SOIL        , 
     M      H2O_BUCKET_FROM_SOIL        , 
     M      ENH_BUCKET_INTO_SOIL        , 
     M      ENH_BUCKET_FROM_SOIL        , 
     M      DUMMY_XXX1     , 
     M      TMP_SOIL1_TMP  )

      ENDIF
!
! 蒸散の層ごとの割り振り 
!
      CALL CANOPY_TRNS_DIV (
     I   NSTR             , NEND              , INSFLAG  , 
     I   TRNS_TMP    , CVR_NOS    , DELT_CNP  ,     
     I   FAC_TRNS_A  , FAC_TRNS_B , 
!    I   KLEV_TRNS_MAX ,
#ifdef MONYOS
     I   ITR  , JCN_ITR_MAX ,
     I   TMP_SOIL     , 
#endif
     M   TRNS_BUCKET_NOS )
!
! 運動量フラックスの計算
!
      CALL CANOPY_GET_GMUV (
     I    NSTR       , NEND      , INSFLAG   , 
     I    RIB        , DELT_ATM  , DELT_CNP  , 
     I    U_REF      , V_REF     , Z_REF     , P_DEL      , 
     I    R_AIR_CONST_INV        , A_COEF    , RHO_AIR    , 
     I    Z0_ms      , 
     I    GMUV       , 
     O    U_FLUX     , V_FLUX    , 
     O    GMUV_NOS   ) 
!
      do i=1,nmax_01(3)
        ta (i) = ta(i)  +    tmp_air(i)         *(1-cvr_snw(i)) 
        ea (i) = ea(i)  +    e_air(i)           *(1-cvr_snw(i)) 
        tg (i) = tg(i)  +    tmp_grsk_prog(i)   *(1-cvr_snw(i)) 
        u2o(i) = u2o(i) +    u2     (i)         *(1-cvr_snw(i)) 
        z0_out(i) = z0_out(i) + z0_ms(i)        *(1-cvr_snw(i))  
!
        d_zero_msx(i) = d_zero_ms(i) 
!
        u_flux_out(i) = u_flux_out(i) + u_flux(i)*(1-cvr_snw(i))          
        v_flux_out(i) = v_flux_out(i) + v_flux(i)*(1-cvr_snw(i))          
        fllh_out(i)   = fllh_out(i)   + fllh  (i)*(1-cvr_snw(i))          
        flsh_out(i)   = flsh_out(i)   + flsh  (i)*(1-cvr_snw(i))          
        qflux_out(i)  = qflux_out(i)  + qflux (i)*(1-cvr_snw(i)) 
      enddo
!
        CALL CANOPY_CALC_AIR (
     I    NSTR      , NEND       , MON       , NTYP       , CVR_NOS   ,       
#ifdef MONYOS
     I    DELT_CNP  , INSFLAG    ,
#endif
     I    TMP_AIR   , E_AIR      , tmp_grsk_prog          , 
     I    TMP_REF_IN, Q_REF      , Z_REF      , P_SFC     , 
     I    U_REF     , V_REF      , 
     I    Z0_ms     , D_ZERO_ms  , FAC_CORTMP , U2        , 
     M    TS        , QS         , US         , VS   ) 
!
      ENDIF
!
!========================================================
!
! ==============
! >> 雪あり域 <<
! ==============
!
      IF ( NMAX_0(1) .LT. NMAX ) THEN
!
        nstr = nmax_01(2)
        nend = nmax
        insflag = 1
!
        DO I=NSTR,NEND
          FLAG_PHASE_SS1(I) = 1.D0
          IF (      
!     1              TMP_SNOW1_TMP(I).EQ.TMP_FREZ_C 
!     1         .OR. 
     1              ABS(RCN_FACT_SS1).LT.1.D-5 
     1       ) THEN
            FLAG_PHASE_SS1(I) = 0.D0
          ENDIF
        ENDDO
!
        CALL RESET ( FAC_TRNS_A , IDX*IDP*ICG )
        CALL RESET ( FAC_TRNS_B , IDX*IDP*ICG )
!
        CALL CANOPY_ROOT (
     I   NSTR            , NEND          , INSFLAG     ,
     I   NTYP            , MON           , 
     I   FR_WTR_SNSL     , 
     I   FR_GLA_SNSL     , 
     I   TMP_SNSL        , 
     O   POT_SOIL        , R_ROOT_BAR_INV, 
     O   HUM_REL_SNSK    , R_BB        )   ! 雪 INSFLAG=1 では 1 / 0
!
! 雪の存在によるキャノピー空間の属性の変更
!
        CALL CANOPY_MOD_AIR_SNW (
     I    NSTR             , NEND              , 
     I    DEP_SNOW         , NTYP       , MON              , 
     I    FAC_AG_INV_ms    , FAC_AC_BAR_INV_ms             , 
     I    XLAI_VEG_ms      , 
     I    D_ZERO_ms        , Z0_ms        , 
     I    CVR_VEG_ms       , CVR_VEG_INV_ms , 
     I    WTR_SNOW         , GLA_SNOW       , 
     O    Z0_msx           , 
     O    FAC_AG_INV_msx   , FAC_AC_BAR_INV_msx            , 
     O    CVR_VEG_msx      , CVR_VEG_INV_msx, 
     O    XLAI_VEG_msx     , D_ZERO_msx     , 
     O    HUM_REL_SNSK     , FR_WBH_SNW     )
!
! ----------------------
! > インターセプション <
! ----------------------
!
        CALL CANOPY_INTERC (
     I    NSTR             , NEND              , INSFLAG , 
     I    NTYP             , 
#ifdef MONYOS
     I    DELT_CNP , 
#endif
     I    RAIN_LARGE        , RAIN_CONV         ,
     I    CVR_VEG_msx       , CVR_SNW           , XLAI_VEG_msx   , 
     I    TMP_REF_IN        , 
     I    TMP_CNP_SNW_PROG  , DUMMY_IDX         , 
     I    FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG , 
     I    DUMMY_IDX         , DUMMY_IDX         , 
     O    TMP_CNP           , DUMMY_XXX2        , 
     O    FR_H2O_CNP        , FR_WBH_CNP        , 
     O    DUMMY_XXX3        , DUMMY_XXX4        , 
     O    H2O_SAT_CNP       , DUMMY_XXX5        , 
     M    H2O_BUCKET_INTO_SNOW , ENH_BUCKET_INTO_SNOW  )
!
!     CALL COPY ( TMP_SNSK_PROG(NSTR) , TMP_SKIN(NSTR) , NEND-NSTR+1 )
      TMP_SKIN(NSTR:NEND)  = TMP_SNSK_PROG(NSTR:NEND) 
!
      CALL CANOPY_COPY ( 
     I    NSTR         , NEND          , INSFLAG    ,  
     I    TMP_CNP      , TMP_SKIN      , TMP_REF_IN , TMP_SNOW1_TMP ,
     I    E_REF        , 
     O    TMP_CNP_SAVE , TMP_SKIN_SAVE , TMP_REF    , TMP_SS1       ,
     O    E_REF_SAVE ) 
!
      CALL CANOPY_CANOPY_INIT_1 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , MON         , 
     I   TMP_REF     , U_REF       , V_REF       , Z_REF      ,   
     I   Z0_msx      , D_ZERO_msx  ,
     I   CVR_VEG_msx , 
     M   R_STM_BAR_PAR_SNW , 
     O   TMP_AIR     , R_AIR_INV   , R_AIR_CONST_INV  , 
     O   U2          , FAC_CORTMP  , 
     O   A_COEF      , 
     O   R_STM_BAR   , FLAG_DRY_CNP, FLAG_DRY_SKIN    )
!
!  ------------------------------------
!  > イテレーション（仮積分）ここから <
!  ------------------------------------
!
        ITR = 0 
 8100 CONTINUE
        ITR = ITR + 1  
!
! ------------------------------
! > キャノピー空間状況決定準備 <
! ------------------------------
!
! イテレーション準備 
!
      CALL CANOPY_CANOPY_INIT_2 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , 
     I   TMP_CNP     , TMP_SKIN    , 
     I   U_REF       , V_REF       , Z_REF      ,   
     I   Z0_msx      , D_ZERO_msx  ,
     I   TMP_REF     , FAC_CORTMP  ,
     O   RIB_CONST   ,
     O   FUNC_TMP    )   
!
      ICOUNT = 0
      INIT_EA = 1
!
#include "noparallel_fixedform"
      DO ITC = 1,ITC_MAX
!
!  ----------------------------------
!  > TA , RIB のイテレーション計算 <
!  ----------------------------------
!
        CALL CANOPY_GET_TA ( 
     I    NSTR             , NEND              , 
     I    NTYP        , MON         , 
     I    TMP_CNP     , TMP_SKIN    , TMP_REF       , 
     I    FAC_AG_INV_msx            , XLAI_VEG_msx  , 
     I    FAC_AC_BAR_INV_msx        , 
     I    U2          , RIB_CONST   , 
     I    FAC_CORTMP  , 
     I    Z_REF       , Z0_msx      , R_AIR_CONST_INV , A_COEF       ,  
     M    ICOUNT      , 
     M    TMP_AIR     , R_AIR_INV   , 
     O    R_AC_BAR    , R_AG        , 
     O    R_AC_BAR_INV, R_AG_INV    , R_SUM_INV     , 
     O    E_AIR_SAT   , E_CNP_SAT   , E_SNSK_SAT    , 
     O    DEC_SAT_DTC , DEG_SAT_DTG , RIB )
!
!  -------------------------------------
!  > EA , 気孔抵抗のイテレーション計算 <
!  -------------------------------------
!
        CALL CANOPY_GET_EA (
     I    NSTR             , NEND              , INSFLAG , 
     I    ITR          , JCN_ITR_MAX      , 
     I    NTYP         , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR     , R_AC_BAR_INV , R_AG        , R_AG_INV    , 
     I    ZERO_IDX     , R_AIR_INV    ,                             ! R_BB=0
     I    R_ROOT_BAR_INV, R_STM_BAR_PAR_SNW, 
     I    E_AIR_SAT    , E_CNP_SAT     , E_REF     , E_SNSK_SAT   ,
     I    HUM_REL_SNSK     , 
     I    TMP_SNSL     , 
     I    FR_H2O_CNP   , DUMMY_IDX   , 
     I    FUNC_TMP     , 
     I    POT_SOIL     ,
     I    D_ZERO_msx   , 
     I    FAC_EVP      ,
     I    ITC          , ITC_MAX      , 
#ifdef MONYOS
     I   DELT_CNP , 
#endif
     M    ICOUNT       ,  INIT_EA    , 
     M    R_STM_BAR    , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN           ,
     O    E_AIR        , FLAG_TR     , FLAG_TR_ALL ,
     O    R_ACX_BAR    , R_AGX_BAR    , R_ACX_BAR_INV , R_AGX_BAR_INV ,
     O    R_AB_INV     , 
!    O    FAC_TRNS_A   , FAC_TRNS_B  , KLEV_TRNS_MAX  ) 
     O    FAC_TRNS_A   , FAC_TRNS_B  )
!
      ENDDO 
!
!  仮積分時の場合、2 度目は値を戻す。
!
      IF ( ITR.EQ.JCN_ITR_MAX) THEN
        CALL CANOPY_COPY ( 
     I    NSTR              , NEND              , INSFLAG   ,  
     I    TMP_CNP_SAVE , TMP_SKIN_SAVE , TMP_REF_IN , TMP_SNOW1_TMP ,
     I    E_REF_SAVE   ,
     O    TMP_CNP      , TMP_SKIN      , TMP_REF    , TMP_SS1       ,
     O    E_REF ) 
!
        DO  I=NSTR,NEND 
          R_SUM_INV(I) 
     1          = 1 .
     1            / ( R_AG_INV(I) + R_AC_BAR_INV(I)*2 + R_AIR_INV(I) )

          TMP_AIR(I) 
     1          =  (   TMP_SKIN(I) * R_AG_INV(I)
     1               + TMP_CNP (I) * R_AC_BAR_INV(I)*2 
     1               + TMP_REF (I) * R_AIR_INV(I)  * FAC_CORTMP(I) )
     1             * R_SUM_INV(I) 
!
          YI = (TMP_AIR(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX ( IY, 1 )
          IY = MIN ( IY, 24999 )
          X = YI - IY
          E_AIR_SAT(I) = ((1.0D0-X)*TABLE(IY)+X*TABLE(IY+1))/EPSFAC_c
! 
          YI = (TMP_CNP(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX( IY, 1 )
          IY = MIN( IY, 24999 )
          X = YI - IY
          E_CNP_SAT   (I)=((1.0D0-X)* TABLE(IY)+X* TABLE(IY+1))/EPSFAC_c
          DEC_SAT_DTC (I)=((1.0D0-X)*DTABLE(IY)+X*DTABLE(IY+1))/EPSFAC_c      
!
          YI = (TMP_SKIN(I)-123.2D0)*100.0D0
          IY = YI
          IY = MAX( IY, 1 )
          IY = MIN( IY, 24999 )
          X = YI - IY
          E_SNSK_SAT  (I)=((1.0D0-X)* TABLE(IY)+X* TABLE(IY+1))/EPSFAC_c
          DEG_SAT_DTG (I)=((1.0D0-X)*DTABLE(IY)+X*DTABLE(IY+1))/EPSFAC_c
        ENDDO
!
        CALL CANOPY_AIR_ALL (
     I    NSTR       , NEND       , INSFLAG    , 
     I    CVR_VEG_msx, CVR_VEG_INV_msx, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , DUMMY_IDX  , 
     I    E_CNP_SAT  , E_REF      , E_SNSK_SAT   ,
     I    HUM_REL_SNSK   , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN     , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
!
      ENDIF
!
! イテレーション後の処理
!   結露・蒸発も考慮した、水・氷が覆う割合 XFR_WTR/GLA を計算
!
      CALL CANOPY_SET_AIR  (
     I    NSTR         , NEND              , INSFLAG , 
     I    FR_H2O_CNP   , FR_WBH_CNP   , DUMMY_IDX   , DUMMY_IDX , 
     I    FR_WBH_SNW   ,    
     I    FLAG_DRY_CNP , FLAG_DRY_SKIN           , 
     O    XFR_WTR_CNP  , XFR_GLA_CNP   , DUMMY_XXX6  , DUMMY_XXX7, 
     O    XFR_WTR_BS   , XFR_GLA_BS ) 
!
!  ----------------------
!  > 長波フラックス計算 <
!  ----------------------
!
      CALL CANOPY_RAD (
     I   NSTR             , NEND              , 
     I  TMP_CNP     , TMP_SKIN    , 
     I  THERMK_SNW  , 
     I  CVR_VEG_msx , 
!    I  FLAG_TMPMOD_CNP   , FLAG_TMPMOD_SKIN  , 
     O  DRAD_DX    )
!
! ----------------------
! > 顕熱フラックス計算 <
! ----------------------
!
      CALL CANOPY_SENS(
     I    NSTR             , NEND              , 
     I    ITR              , JCN_ITR_MAX           , 
#ifdef MONYOS
     I  DELT_CNP    , INSFLAG        , 
#endif 
     I  R_AIR_INV   , R_AC_BAR_INV, R_AG_INV  , R_SUM_INV , 
     I  RHO_CP      , FAC_CORTMP  , 
     I  TMP_CNP     , TMP_SKIN    , TMP_AIR     , 
!    I  FLAG_TMPMOD_CNP   , FLAG_TMPMOD_SKIN  , 
     o  dta_dtc , dta_dtg , dta_dtr ,
     O  SENS        , 
     O  DSENS_DX    ) 
!
! ----------------------
! > 潜熱フラックス設定 <
! ----------------------
!
      CALL CANOPY_EVAP (
     I    NSTR             , NEND              , e_air ,
     I    ITR               , JCN_ITR_MAX           , 
#ifdef MONYOS
     I  DELT_CNP    , INSFLAG           , 
#endif
     I  R_AIR_INV   , 
     I  R_AC_BAR_INV, R_ACX_BAR_INV , 
     I  R_AG_INV    , R_AGX_BAR_INV , R_AB_INV   ,  
     I  E_REF       , E_CNP_SAT     , E_SNSK_SAT     , 
     I  DEC_SAT_DTC , DEG_SAT_DTG , 
     I  CVR_VEG_msx ,                               ! 被覆率
     I  FAC_EVP     , HUM_REL_SNSK    ,                 ! 雪なら後ろは 1 
     I  FAC_EVP_2   ,
     I  XFR_WTR_CNP  , DUMMY_IDX   , XFR_WTR_BS  ,   ! 水氷割合（結露考慮）
     I  XFR_GLA_CNP  , DUMMY_IDX   , XFR_GLA_BS  ,   ! 水氷割合（結露考慮） 
!    I  FLAG_TMPMOD_CNP   , FLAG_TMPMOD_SKIN  ,   
     I  FLAG_DRY_CNP      , FLAG_DRY_SKIN     ,
     I  FLAG_TR_ALL ,
     O  EVAP        , SUBL        , TRNS        , 
     O  DEVAP_DX   )
!
! 熱伝導フラックス設定 ... 下向き正
!  
      CALL CANOPY_HDF ( 
     I    NSTR             , NEND              , 
     I    ITR     , JCN_ITR_MAX , 
#ifdef MONYOS
     I    INSFLAG , DELT_CNP , 
#endif
!    I    TMP_SKIN  , TMP_SS1     , TMP_SOIL(1,1)  ,   ! snow1/soil1
     I    TMP_SKIN  , TMP_SS1     , TMP_SNOW(1,2)  ,   ! snow1/soil1
     I    HDF_SNSK_IN  , 
     O    HDF  )
!
      hdf_skin(NMAX_01(2):NMAX,:) 
     1      = hdf_snsk_in(NMAX_01(2):NMAX,:) * fact_tskin_top_c

!
! 方程式を立てる
!
      CALL CANOPY_FLUX_EQSET (
     I  NSTR         , NEND         , INSFLAG      , 
     I  ITR - JCN_ITR_MAX , 
     I  XLAI_VEG_msx , 
     I  FR_H2O_CNP   , FR_WBH_CNP   , DUMMY_IDX    , DUMMY_IDX , 
     I  H2O_SAT_CNP  , DUMMY_IDX    , 
     I  SENS         ,                               DSENS_DX  , 
     I  EVAP         , SUBL         , TRNS         , DEVAP_DX  , 
     I  DRAD_DX      , 
     I  RAD_SHRT_SNW , RAD_LONG_SNW , 
     I  HDF          , HDF_skin     , 
     I  DELT_ATM     , DELT_CNP     , 
     I  GMT          , GMQ          ,  P_DEL       , Q_REF     , 
     I  CAP_SNOW_TOP , FLAG_PHASE_SS1 , 
     O  SIB_MTX   )
!
! 方程式を解く
!
      CALL CANOPY_FLUX_EQSLV (
     I  NSTR             , NEND              , 
     M  SIB_MTX    )
!
! 積分
!
      IF ( ITR .LT. JCN_ITR_MAX ) THEN
!
        CALL CANOPY_INTEGRATION_OVERIMP_ITR (
     I      NSTR       , NEND        , 
     I      DELT_CNP   , DELT_ATM    , 
     I      SIB_MTX    , 
     I      FAC_EVP_2  ,
     I      DTA_DTC    , DTA_DTG     , DTA_DTR    , 
     M      TMP_CNP    , TMP_SKIN    , TMP_SS1    , 
     M      TMP_REF    , TMP_AIR     , E_REF     )
!
        GOTO 8100
!
      ELSE 
!
        CALL CANOPY_INTEGRATION_OVERIMP_RUN ( 
     I     NSTR       , NEND        , INSFLAG  ,IMSKLT , 
     I     DELT_CNP   , DELT_ATM    , 
     I     SIB_MTX    , 
     I     EVAP       , SUBL        , TRNS     , 
     I     SENS       , 
     I     RAD_LONG_SNW, 
     I     RAD_SHRT_SNW, 
     I     HDF      , 
     I     DEVAP_DX   , DSENS_DX    , DRAD_DX  , HDF_skin , 
     I     H2O_SAT_CNP    , DUMMY_IDX      , 
     I     CVR_SNW        , XLAI_VEG_msx   , 
     I     TMP_CNP        , 
     I     FR_H2O_CNP     , FR_WBH_CNP     , 
     I     TMP_SKIN       , 
     I     DUMMY_IDX      , DUMMY_IDX      , 
     I     TMP_SNOW       , 
!    I     SRAD_INTO_SNOW_LAYER_TMP   , 
     I     RAD_SHRT_INTO_SNOW_LAYER   , 
     I     CAP_SNOW_TOP   , FLAG_PHASE_SS1 , 
     I     GMT            , GMQ          ,  
     O     TMP_CNP_SNW_PROG, FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG, 
     O     TMP_SNSK_PROG  , DUMMY_XXX8     , DUMMY_XXX9     ,
     O     TRNS_TMP       , 
     O     GMT_SNW        , GMQ_SNW      , 
     O     FLLH           , FLSH       , qflux , 
     M     RAD_LONG_SNW_2_ATM   , 
     M     ENG_BUCKET_SNW              , 
!!     M     WTR_BUCKET_IN_SNW           , GLA_BUCKET_IN_SNW            , 
!!     M     WTR_BUCKET_OUT_SNW          , GLA_BUCKET_OUT_SNW           , 
!!     M     ENG_WTR_BUCKET_IN_SNW       , ENG_GLA_BUCKET_IN_SNW        ,
!!     M     ENG_WTR_BUCKET_OUT_SNW      , ENG_GLA_BUCKET_OUT_SNW       , 
     M      H2O_BUCKET_INTO_SNOW        , 
     M      H2O_BUCKET_FROM_SNOW        , 
     M      ENH_BUCKET_INTO_SNOW        , 
     M      ENH_BUCKET_FROM_SNOW        , 
     M     SRAD_BUCKET_SNW , 
     M     TMP_SNOW1_TMP   ) 
!    M     TMP_SS1       )
!
      ENDIF     
!
! 蒸散の層ごとの割り振り
!
      CALL CANOPY_TRNS_DIV (
     I   NSTR , NEND , INSFLAG ,
     I   TRNS_TMP    , CVR_SNW    , DELT_CNP  ,     
     I   FAC_TRNS_A  , FAC_TRNS_B , 
!    I   KLEV_TRNS_MAX ,
#ifdef MONYOS
     I   ITR  , JCN_ITR_MAX ,
     I   TMP_SNSL     , 
#endif
     M   TRNS_BUCKET_SNW )
!
!
! 運動量フラックスの計算
!
      CALL CANOPY_GET_GMUV (
     I    NSTR       , NEND      , INSFLAG   , 
     I    RIB        , DELT_ATM  , DELT_CNP  , 
     I    U_REF      , V_REF     , Z_REF     , P_DEL  , 
     I    R_AIR_CONST_INV        , A_COEF    , RHO_AIR    , 
     I    Z0_msx     , 
     I    GMUV       , 
     O    U_FLUX     , V_FLUX     , 
     O    GMUV_SNW   )
!
      do i=nmax_01(2),nmax 
        ta (i) = ta(i)  +    tmp_air(i)         *cvr_snw(i) 
        ea (i) = ea(i)  +    e_air(i)           *cvr_snw(i) 
        tg (i) = tg(i)  +    tmp_snsk_prog(i)   *cvr_snw(i) 
        u2o(i) = u2o(i) +    u2     (i)         *cvr_snw(i) 
        z0_out(i) = z0_out(i) + z0_msx(i)       *cvr_snw(i)  
!
        u_flux_out(i) = u_flux_out(i) + u_flux(i)*cvr_snw(i)          
        v_flux_out(i) = v_flux_out(i) + v_flux(i)*cvr_snw(i)          
        fllh_out(i)   = fllh_out(i)   + fllh  (i)*cvr_snw(i)          
        flsh_out(i)   = flsh_out(i)   + flsh  (i)*cvr_snw(i)          
        qflux_out(i)  = qflux_out(i)  + qflux (i)*cvr_snw(i) 
      enddo
!
        CALL CANOPY_CALC_AIR (
     I    NSTR      , NEND       , MON       , NTYP       , CVR_SNW   , 
#ifdef MONYOS
     I    DELT_CNP  , INSFLAG    ,
#endif
     I    TMP_AIR   , E_AIR      , TMP_SNSK_PROG          ,
     I    TMP_REF_IN, Q_REF      , Z_REF      , P_SFC     , 
     I    U_REF     , V_REF      , 
     I    Z0_msx    , D_ZERO_msx , FAC_CORTMP , U2        , 
     M    TS        , QS         , US         , VS  ) 
!
      ENDIF
!
! 大気拡散行列の変更（平均処理）
!
      CALL CANOPY_AVR_NOSSNW  ( 
     I    NMAX       , NMAX_0     , NMAX_01    ,  NMAX_1   , 
#ifdef MONYOS
     I    DELT_CNP   , 
#endif
     I    CVR_NOS    , CVR_SNW    , 
     I    GMUV_NOS   , GMUV_SNW   , 
     I    GMT_NOS    , GMT_SNW    , 
     I    GMQ_NOS    , GMQ_SNW    , 
     I    RAD_LONG_NOS_2_ATM      , RAD_LONG_SNW_2_ATM      , 
     O    GMT        , GMQ        , GMUV       , 
     O    RAD_LONG_SIBMIX_2_ATM   )
! 
      call canopy_monit_meanvalues ( 
     i  nmax         , DELT_CNP   , 
     i  u_flux_out   , v_flux_out , 
     i  fllh_out     , flsh_out   , qflux_out ) 
!
      RETURN      
      END SUBROUTINE CANOPY_MAIN__run 
!==================================================
      SUBROUTINE CANOPY_EXTEND (
     I   NMAX          , 
     I   NTYP          , MON           , 
     O   FAC_AG_INV_ms , 
     O   FAC_AC_BAR_INV_ms  , 
     O   XLAI_VEG_ms   , 
     O   D_ZERO_ms     , 
     O   Z0_ms        ) 
!
! 展開する
!
      USE SIBPRM , ONLY :
     1   IDX    ,
     1   ITYP   , 
     1   ICG  
!
      USE SIBVEG , ONLY :
     1           FAC_AG_INV_mo              , ! R_AG 関連定数   (ITYP,IMON)
     1           FAC_AC_BAR_INV_mo          , ! R_AC 関連定数   (ITYP,IMON)  
     1           XLAI_mo                    , ! LAI (NTYP(I),MON,1)
     1           D_ZERO_mo                  , ! 零面変位 d   (ITYP,IMON)  
     1           Z0_mo                        ! 粗度 z0 (ITYP,IMON)  
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NTYP (IDX)  ! 植生タイプ
      INTEGER,INTENT(IN)   :: MON          ! 月
! 
      REAL(8),INTENT(OUT)  :: FAC_AC_BAR_INV_ms(IDX)   ! 旧 RBC 
      REAL(8),INTENT(OUT)  :: FAC_AG_INV_ms    (IDX)   ! 旧 RDC 
      REAL(8),INTENT(OUT)  :: XLAI_VEG_ms      (IDX,ICG)   !  
      REAL(8),INTENT(OUT)  :: D_ZERO_ms        (IDX)   !  
      REAL(8),INTENT(OUT)  :: Z0_ms            (IDX)   !  
!
! WORK  
!
      INTEGER :: I
!
      LOGICAl,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_EXTEND 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=1,NMAX
        FAC_AC_BAR_INV_ms(I)   = FAC_AC_BAR_INV_mo(NTYP(I),MON) 
        FAC_AG_INV_ms    (I)   = FAC_AG_INV_mo    (NTYP(I),MON) 
        XLAI_VEG_ms      (I,1) = XLAI_mo          (NTYP(I),MON,1)
        XLAI_VEG_ms      (I,2) = XLAI_mo          (NTYP(I),MON,2)
        D_ZERO_ms        (I)   = D_ZERO_mo        (NTYP(I),MON)
        Z0_ms            (I)   = Z0_mo            (NTYP(I),MON)
 1000 CONTINUE
!
      RETURN
      END subroutine canopy_extend
!----------------------------------
      subroutine canopy_monit_meanvalues ( 
     i  nmax         , delt_cnp   , 
     i  uflux        , vflux      , 
     i  fllh         , flsh       , qflux     ) 
!
      USE SIBPRM , ONLY : 
     1    IDX     
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, INA, 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , 
     1     IMONIT_ALMA
#endif
!
      implicit none
!
      integer,intent(in) :: nmax 
      real(8),intent(in) :: delt_cnp
      real(8),intent(in) :: uflux (idx) 
      real(8),intent(in) :: vflux (idx) 
      real(8),intent(in) :: fllh  (idx) 
      real(8),intent(in) :: flsh  (idx) 
      real(8),intent(in) :: qflux (idx)  ! Upward positive 
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
      REAL(8)           :: WORK(IDX) 
#endif

#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
         cmark = 'A_LH'
         CALL MONIT_SIB ( CMARK , FLLH  ,INA, DELT_CNP, ONE ) 
         cmark = 'A_SH'
         CALL MONIT_SIB ( CMARK , FLSH  ,INA, DELT_CNP, ONE ) 
         cmark = 'A_EST'
         CALL MONIT_SIB ( CMARK , QFLUX ,INA, DELT_CNP, ONE ) 
!
         cmark = 'A_QTAU'
         WORK(:) = 0.
         WORK(1:NMAX) = SQRT(UFLUX(1:NMAX)**2 + VFLUX(1:NMAX)**2)
         CALL MONIT_SIB ( CMARK , WORK ,INA, DELT_CNP, ONE ) 
      ENDIF 
#endif
      return
      end subroutine canopy_monit_meanvalues 
!---------------------------------------------
      subroutine mr_canopy_monit_meanvalues 
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , 
     1     IMONIT_ALMA
!
      implicit none
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_CANOPY_MONIT_MEANVALUES 04/09 HOSAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
         cmark = 'A_LH'
         CTITLE = 'ALMA Qle latent heat'
         CUNIT  = 'W/M**2'
         CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
         cmark = 'A_SH'
         CTITLE = 'ALMA Qh sensible heat'
         CUNIT  = 'W/M**2'
         CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
         cmark = 'A_EST'
         CTITLE = 'ALMA Evap : evapotrancipitation + sublim'
         CUNIT  = 'kg/M**2/s'
         CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
         cmark  = 'A_QTAU'
         CTITLE = 'ALMA Qtau momentum flux'
         CUNIT  = 'N/M**2'
         CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif      
      return
      end subroutine mr_canopy_monit_meanvalues 
!----------------------------------
      end module canopy_main


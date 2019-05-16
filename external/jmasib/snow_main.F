! 雪のメインプログラム
! 99/07/11 保坂征宏 第一版
! 99/07/15 保坂征宏 RHO_SNOW を加える
!                   とりあえず rho_snow の圧密や、雪融け補正は保留。
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SNOW_MAIN (
     I    NMAX     , NMAX_0     , NMAX_01    , NMAX_1     , 
     I    NMAX_2   , NMAX_3     , NMAX_4     , imsklt     , 
     I    CVR_SNW  , CVR_SNW_INV, DELT_SNOW  , DELT_CNP   ,
     I    DEP_SNOW , 
     I    HDF_SNOW_H , 
     I    H2O_BUCKET_INTO_SNOW  , H2O_BUCKET_FROM_SNOW    , 
     I    ENH_BUCKET_INTO_SNOW  , ENH_BUCKET_FROM_SNOW    , 
     I    ENG_BUCKET_SNW        , SRAD_BUCKET_SNW   , 
     I    ENG_SNOW_BUCKET       , H2O_SNOW_BUCKET         , 
     I    TMP_REF       , U_REF         , V_REF       ,   
     I    TMP_SNSK_NEXT , 
     I    ENG_SOIL_TOP  , 
     M    WTR_SNOW_PROG , GLA_SNOW_PROG   , TMP_SNOW_PROG   ,      
     M    AGE_SNOW_PROG , RHO_SNOW_INV_PROG     ,  
     M    ENG_BUCKET_INTO_SOIL  , 
     M    H2O_BUCKET_INTO_SOIL  , ENH_BUCKET_INTO_SOIL     , 
     O    CVR_SNW_NEXT          , CVR_SNW_INV_NEXT         ,    
     O    CVR_NOS_NEXT          , CVR_NOS_INV_NEXT         , 
     O    INFO_SNOW_NEXT        , 
     O    DEP_SNOW_NEXT         , 
     O    ENG_BUCKET_SNOW_TO_SNSL      ,  
     O    WTR_BUCKET_SNOW_TO_SNSL  , ENH_BUCKET_SNOW_TO_SNSL  , 
     O    ENG_SNOW_BUCKET_NEXT     , H2O_SNOW_BUCKET_NEXT     , 
     O    BERG                   )
!
!
!  方針
!    ・熱伝導係数は、旧値を使用して、soil_main に入る前に計算しておく。
!
!    ・蒸発・昇華した水分・氷分は soil_pre で差し引く。
!        これにより熱伝導係数は再計算しない。
!        （蒸発、昇華があってもなくても同じとする。）
!        ただし比熱は変化し、解やフラックス補正量は変わる。
!
!    ・soil_h* では、水の移動や相変化はないとして熱伝導方程式を解く。
!      熱伝導等で各層が得るエネルギーフラックス収束量を評価する。
!        ただしもともと零度である相については、温度が変わらないと判断して、
!        陰的に評価することに伴う補正を行なわない。
!    ・この時点で状況を固定する。各層の各相がもつエネルギー量などが
!      フィックスする。
!
!    ・soil_w* では、水の移動を扱う。
!
!  課題
!    ・RHO_SNOW の積分。
!        まだできていない。
!        圧密をどこでするか。
!        それ以外の変化はどこで考慮するか
!
!  注意
!     ・ *_TOP および *SNOW_TO_SNSL は、単位格子あたりであって、
!        雪あり領域あたりではない。
!
!  雪層について、
!    熱の方程式を立てて解く。
!    水の落下を解く。
!    次ステップでの雪の分布を決める。
!
!  なお、ここでは以下のことは行わない 
!    
!    雪のアルベド計算
!      ... これは放射の前。
!
!    雪の熱伝導係数
!      ... これは土壌の熱伝導係数とも関係するので、この以前に計算する。
!
!    雪の上端・内部に入ってくる放射など、上端に入って来るもの
!      ... これはこの以前に計算する。
!
!          短波  ...  SRAD_BUCKET_SNW 
!          熱    ...  ENG_BUCKET_SNW 
!          雪    ...  ENG_GLA_SNW,GLA_SNW,ENG_WTR_SNW,WTR_SNW の各 BUCKET
!                     これら BUCKET 量は、*CVR_SNW をかけた量が入力
!                     雪無し格子に降ってきた雪を扱うため。
!
!    雪面被覆率が変わった場合のアルベドの補正 
!      ... 下草・裸地のアルベドにもかかわるので、この後で行う。
!
!    雪共通配列への詰め替え 雪の時間積分完了
!      ... あとで一気に物理格子に詰め込む。
!
!    CVR_SNW の積分
!
! 注意
!   雪の無かった格子に雪が降ってきた場合は、熱まわりの計算はしない。
!   次のステップでの雪面積を評価するために水まわりの計算はおこなう。
!   その結果の雪の量が少ない場合はバケツに貯めたままとする。
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ISN      
!
      USE COM_STEP_SIB0109 , ONLY :
     1   ICN_SIB0109_SNOW_RESET                  ! 雪層をリセットする 1
!
      USE SIB_MONIT , ONLY :
     1    SIB_MONIT_RUN_SET_SNOWNEXTINFO 
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1   JCN_SNOW_HEAT_EXPLICIT 
!
#ifdef SIB_DEBUG
      use sibcon , only : 
     1   HL_W2G_c , TMP_FREZ_C , CAP_GLA_c , CAP_WTR_c 
#endif
!
      IMPLICIT NONE 
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: IMSKLT(IDX)
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NMAX_0 (3)
      INTEGER,INTENT(IN)    :: NMAX_01(3)
      INTEGER,INTENT(IN)    :: NMAX_1 (3)
      INTEGER,INTENT(IN)    :: NMAX_2 (3)
      INTEGER,INTENT(IN)    :: NMAX_3 (3)
      INTEGER,INTENT(IN)    :: NMAX_4 (3)
!
!  予報変数
!
      REAL(8),INTENT(INOUT) :: TMP_SNOW_PROG    (IDX,ISN+1) ! 雪の温度 
      REAL(8),INTENT(INOUT) :: WTR_SNOW_PROG    (IDX,ISN)
      REAL(8),INTENT(INOUT) :: GLA_SNOW_PROG    (IDX,ISN)
      REAL(8),INTENT(INOUT) :: RHO_SNOW_INV_PROG(IDX,ISN)   ! 雪密度逆数
      REAL(8),INTENT(INOUT) :: AGE_SNOW_PROG   (IDX)        ! 雪の年齢
!
      REAL(8),INTENT(INOUT) :: ENG_BUCKET_INTO_SOIL(IDX)  
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO_SOIL(IDX)  
      REAL(8),INTENT(IN)    :: ENH_BUCKET_INTO_SOIL(IDX)   ! 使わない

!
      REAL(8),INTENT(IN)    :: ENG_SNOW_BUCKET(IDX)
      REAL(8),INTENT(IN)    :: H2O_SNOW_BUCKET(IDX)
!
!   大気からの強制
!
      REAL(8),INTENT(IN)    :: ENH_BUCKET_INTO_SNOW (IDX)  ! 雨水のエネルギー
      REAL(8),INTENT(IN)    :: ENH_BUCKET_FROM_SNOW (IDX)  ! 蒸発水
      REAL(8),INTENT(IN)    :: H2O_BUCKET_INTO_SNOW (IDX)  ! 雨の量
      REAL(8),INTENT(IN)    :: H2O_BUCKET_FROM_SNOW (IDX)  ! 蒸発
      REAL(8),INTENT(IN)    :: SRAD_BUCKET_SNW       (IDX)  ! 短波入射
      REAL(8),INTENT(IN)    :: ENG_BUCKET_SNW        (IDX)  ! 熱フラックス
!
!   参照空気の情報 ( 雪の温度の決定材料 )
!
      REAL(8),INTENT(IN)    :: U_REF    (IDX) 
      REAL(8),INTENT(IN)    :: V_REF    (IDX)
      REAL(8),INTENT(IN)    :: TMP_REF  (IDX)
!
!   その他
!
      REAL(8),INTENT(IN)    :: TMP_SNSK_NEXT (IDX) ! 
!     REAL(8),INTENT(INOUT) :: TMP_GRSK_NEXT (IDX) ! 
!
      REAL(8),INTENT(IN)    :: ENG_SOIL_TOP  (IDX) ! 
!
      REAL(8),INTENT(IN)    :: DEP_SNOW    (IDX,ISN+1) ! 雪の厚さ ISN+1 は合計
      REAL(8),INTENT(IN)    :: CVR_SNW     (IDX)
      REAL(8),INTENT(IN)    :: CVR_SNW_INV (IDX)
!
      REAL(8),INTENT(IN)    :: HDF_SNOW_H  (IDX,ISN)   ! 雪中熱伝導係数
!     REAL(8),INTENT(IN)    :: HDF_SNOW_BH (IDX)       ! 雪B と地面
!
      REAL(8),INTENT(IN)    :: DELT_SNOW               ! 時間刻
      REAL(8),INTENT(IN)    :: DELT_CNP                ! 時間刻
!
! OUTPUT  
!
!
      REAL(8),INTENT(OUT)   :: ENG_SNOW_BUCKET_NEXT(IDX)
      REAL(8),INTENT(OUT)   :: H2O_SNOW_BUCKET_NEXT(IDX)
!
      REAL(8),INTENT(OUT)   :: DEP_SNOW_NEXT  (IDX,ISN+1) ! 雪の厚さ 
!
      REAL(8),INTENT(OUT)   :: BERG (IDX)               ! 雪の流出 (10m以上)
!                 水当量 10m を越えたら、
!                 バッファの層 ( 01/10/29 現在第 2 層 ) から流出させる。
!
!   土壌への強制
!
      REAL(8),INTENT(OUT)   :: ENG_BUCKET_SNOW_TO_SNSL     (IDX) ! 熱伝導
!
      REAL(8),INTENT(OUT)   :: ENH_BUCKET_SNOW_TO_SNSL (IDX) ! 落下水の ENG
      REAL(8),INTENT(OUT)   :: WTR_BUCKET_SNOW_TO_SNSL (IDX) ! 落下水量
!
!   その他
!
      INTEGER,INTENT(OUT)   :: INFO_SNOW_NEXT   (IDX)  ! 層数
      REAL(8),INTENT(OUT)   :: CVR_SNW_NEXT     (IDX)  
      REAL(8),INTENT(OUT)   :: CVR_NOS_NEXT     (IDX)  
      REAL(8),INTENT(OUT)   :: CVR_SNW_INV_NEXT (IDX)  
      REAL(8),INTENT(OUT)   :: CVR_NOS_INV_NEXT (IDX)  
!
!
! WORK 熱関連 ( 雪あり域だけでの平均 ) 
!
      INTEGER   :: I
!     REAL(8)   :: ENG_BUC (IDX)                
!
      REAL(8)   :: TMP_SNOW     (IDX,ISN+1) ! 雪の温度 最後は TMP_SNSL
      REAL(8)   :: WTR_SNOW     (IDX,ISN)
      REAL(8)   :: GLA_SNOW     (IDX,ISN)
!
      REAL(8)   :: HDIFF_FLUX_H    (IDX,0:ISN)     ! 熱伝導フラックス
!
      REAL(8)   :: PHASE_SNOW   (IDX,ISN)        ! 零度なら 0、そうでなければ 1
      REAL(8)   :: SIB_MTX_SNOW (IDX,ISN,ISN+1)  ! 熱用の行列
      REAL(8)   :: CAP_SNOW     (IDX,ISN)        ! 比熱
      REAL(8)   :: HFLUX_CONV(IDX,ISN)       ! 雪の層が得る熱
!
      REAL(8)   :: SRAD_FLUX_CNV (IDX,ISN)       ! 雪の層中の短波フラックス
!
      REAL(8)   :: SRAD_BUCKET_SNOW_TO_SNSL (IDX) ! 短波放射
!
      REAL(8)   :: DELT_SNOW_INV    
!
! WORK 水関連 ( 格子平均 ) 
!
      REAL(8)   :: ENG_SNOW        (IDX,ISN)
      REAL(8)   :: H2O_SNOW        (IDX,ISN)
      REAL(8)   :: RHO_SNOW_INV    (IDX,ISN)
      REAL(8)   :: RHO_PRECIPI_INV (IDX) ! 降雪の雪の密度の逆数
!
      REAL(8)   :: ENG_PRECIPI (IDX)
      REAL(8)   :: H2O_PRECIPI (IDX)
!
      REAL(8)   :: ENG_SNOW_NEXT   (IDX,ISN) 
      REAL(8)   :: H2O_SNOW_NEXT   (IDX,ISN)
!
      REAL(8)   :: ENH_BERG    (IDX)
!
      INTEGER   :: IFLG_SLEV_CHANGE (IDX) ! 雪層のリセットをするかどうか
!
!  不要かも知れない。
      INTEGER   :: LIST_SNOW (IDX,0:4)      ! リストベクトル ( 意味は? )
      INTEGER   :: ICNT_SNOW (    0:4)      ! ? 
!
      INTEGER   :: I_SNOW_LEV      (IDX)    ! 現在の雪の層数 
!
#ifdef SIB_DEBUG
      REAL(8)   :: ENG_SNOW_SUM_SAVE  (IDX)
      REAL(8)   :: H2O_SNOW_SUM_SAVE  (IDX)
      REAL(8)   :: ENG_SNOW_SUM_TMP   (IDX)
      REAL(8)   :: H2O_SNOW_SUM_TMP   (IDX)
      REAL(8)   :: ENG_SNOW_SUM_NEW   (IDX)
      REAL(8)   :: H2O_SNOW_SUM_NEW   (IDX)
      REAL(8)   :: H2O_ERR_MAX
      REAL(8)   :: ENG_ERR_MAX
      REAL(8)   :: TMP_SNOW_SAVE    (IDX,ISN) 
      REAL(8)   :: WTR_SNOW_SAVE    (IDX,ISN)
      REAL(8)   :: GLA_SNOW_SAVE    (IDX,ISN)
#endif
!
      LOGICAL,SAVE :: LFIRST =  .TRUE. 
!
!     integer :: icount = 0
!     icount = icount + 1  
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_MAIN     POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
      DELT_SNOW_INV = 1 / DELT_SNOW
!
!
      ENG_BUCKET_SNOW_TO_SNSL(:) = 0. 
!
      H2O_SNOW     (:,:) = 0.
      ENG_SNOW     (:,:) = 0.
!
#ifdef SIB_DEBUG
      WTR_SNOW_SAVE(:,:) = WTR_SNOW_PROG(:,:)  
      GLA_SNOW_SAVE(:,:) = GLA_SNOW_PROG(:,:)  
      TMP_SNOW_SAVE(:,1:4) = TMP_SNOW_PROG(:,1:4)  
!
      DO I=1,NMAX
        H2O_SNOW_SUM_SAVE(I) = 
     1       H2O_SNOW_BUCKET(I)      
     1     + H2O_BUCKET_INTO_SNOW(I) 
     1     + H2O_BUCKET_INTO_SOIL(I) 
     1     + H2O_BUCKET_FROM_SNOW(I) 
     1     + WTR_SNOW_SAVE    (I,1) 
     1     + WTR_SNOW_SAVE    (I,2) 
     1     + WTR_SNOW_SAVE    (I,3) 
     1     + WTR_SNOW_SAVE    (I,4) 
     1     + GLA_SNOW_SAVE    (I,1) 
     1     + GLA_SNOW_SAVE    (I,2) 
     1     + GLA_SNOW_SAVE    (I,3) 
     1     + GLA_SNOW_SAVE    (I,4) 
!
        ENG_SNOW_SUM_SAVE(I) = 
     1       SRAD_BUCKET_SNW(I)
     1     + ENG_BUCKET_SNW (I)
     1     + ENG_BUCKET_INTO_SOIL (I)
     1     + ENG_SNOW_BUCKET(I) 
     1     + ENH_BUCKET_INTO_SNOW(I) 
     1     + ENH_BUCKET_FROM_SNOW(I) 
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,1) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,1) )
     1           * ( TMP_SNOW_SAVE(I,1) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,1)
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,2) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,2) )
     1           * ( TMP_SNOW_SAVE(I,2) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,2)
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,3) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,3) )
     1           * ( TMP_SNOW_SAVE(I,3) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,3)
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,4) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,4) )
     1           * ( TMP_SNOW_SAVE(I,4) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,4)
      ENDDO
#endif
!
!  ======================
!  >> 蒸発・昇華の除去 <<
!  ======================
!
!
      CALL SNOW_PRE (
     I    NMAX      , NMAX_01   , 
     I    NMAX_2    , NMAX_3    , NMAX_4  , 
     I    CVR_SNW_INV    , 
     I    WTR_SNOW_PROG  , GLA_SNOW_PROG    , 
     I    H2O_BUCKET_FROM_SNOW   , ENH_BUCKET_FROM_SNOW   , 
     I    TMP_SNOW_PROG  , 
     O    TMP_SNOW       , WTR_SNOW       , GLA_SNOW      )   
! 
!  ==================
!  >> 熱まわり計算 <<
!  ==================
!
!  ----------
!  >> 準備 <<
!  ----------
!     ・比熱の計算
!     ・陽に評価される熱収束の計算
!        ・各層の短波吸収量 SRAD_SNOW 計算 
!        ・上端からのエネルギーフラックス ENG_BUC
!        ・熱伝導フラックス HDIFF_FLUX_CNV 計算
!    
      CALL SNOW_HFLUX_PRE (
     I    NMAX      , NMAX_01   , NMAX_1  , 
     I    NMAX_2    , NMAX_3    , NMAX_4  , 
     I    DELT_SNOW , DELT_SNOW_INV         , 
     I    WTR_SNOW  , GLA_SNOW  , 
     I    SRAD_BUCKET_SNW       , ENG_BUCKET_SNW      ,   
     I    DEP_SNOW  , TMP_SNOW  , 
     I    CVR_SNW   , CVR_SNW_INV         , 
     I    HDF_SNOW_H, 
     O    PHASE_SNOW,  
     O    CAP_SNOW  , 
     O    HFLUX_CONV , 
     O    SRAD_BUCKET_SNOW_TO_SNSL , 
     O    HDIFF_FLUX_H  )
!
!  ====================================
!  >> 熱の式を解く ～ インプリシット <<
!  ====================================
!
      IF ( JCN_SNOW_HEAT_EXPLICIT .EQ. 0 ) THEN
!
!  ------------------------
!  >> 熱の方程式を立てる <<
!  ------------------------
!     ・熱伝導方程式をたてる
!     ・陽に評価される熱フラックス収束 HFLUX_CONV 計算
!
      CALL SNOW_HFLUX_IMPLICIT_EQSET (
     I  NMAX       , NMAX_01   , 
     I  NMAX_2     , NMAX_3    , NMAX_4    , 
     I  DELT_SNOW  , 
     I  PHASE_SNOW , 
     I  CAP_SNOW   , HDIFF_FLUX_H ,
     I  HDF_SNOW_H , 
     O  SIB_MTX_SNOW    , 
     M  HFLUX_CONV )
!
!  ----------------------
!  >> 熱の方程式を解く <<
!  ----------------------
!    これで、dT/dt が得られる。
!
      CALL SNOW_HFLUX_IMPLICIT_EQSLV (
     I  NMAX_01   , 
     I  NMAX_1    , NMAX_2    , NMAX_3    , NMAX_4    , 
     M  SIB_MTX_SNOW )
!
!  --------------------------
!  > 熱フラックス収束の補正 <
!  --------------------------
!     ・熱フラックス収束 HFLUX_CONV にインプリシット分を加える。
!     ・SNSL に与える熱バケツ計算 
!    
      CALL SNOW_HFLUX_IMPLICIT_CORR (
     I  NMAX      , NMAX_01   , NMAX_1    , 
     I  NMAX_2    , NMAX_3    , NMAX_4    , 
     I  SIB_MTX_SNOW          , 
     I  HDF_SNOW_H            , 
     I  DELT_SNOW             , PHASE_SNOW            , 
     I  CVR_SNW               ,              
     I  SRAD_BUCKET_SNOW_TO_SNSL     ,  
     I  HDIFF_FLUX_H          , 
     O  ENG_BUCKET_SNOW_TO_SNSL      , 
     M  HFLUX_CONV          )  
! 
!  ======================================
!  >> 熱の式を解く ～ イクスプリシット <<
!  ======================================
!
      ELSE
        CALL SNOW_HFLUX_EXPLICIT (
     I  NMAX      , NMAX_01   , NMAX_1    ,
     I  NMAX_2    , NMAX_3    , NMAX_4    , 
     I  DELT_SNOW , CVR_SNW   ,   
! 20050517
!     I  SRAD_FLUX_CNV  , 
     I  SRAD_BUCKET_SNOW_TO_SNSL     ,  
     I  HDIFF_FLUX_H     , 
     M  HFLUX_CONV    , 
     O  ENG_BUCKET_SNOW_TO_SNSL )


!
      ENDIF
!
!  ==================
!  >> 水まわり計算 <<
!  ==================
!   水の移動に伴う変更
!
!   ここから先は、NMAX_01 の部分は、
!   格子平均値 （*CVR_SNW した値）になる。
!
!  ----------------
!  > 水まわり準備 <
!  ----------------
!     ・格子平均値にする
!     ・新雪・旧雪に分離、各層のエネルギーなど計算
!     ・年齢の加算
!     ・雪層リセットフラグ IFLG_SLEV_CHANGE の初期化等
!
      CALL SNOW_WFLUX_PRE (
     I   NMAX    , DELT_SNOW         , 
     I   NMAX_0  , NMAX_01 , NMAX_2 , NMAX_3  , NMAX_4  , 
     I   GLA_SNOW, WTR_SNOW, TMP_SNOW, 
     I   CVR_SNW , 
     I   HFLUX_CONV     , 
     I   H2O_BUCKET_INTO_SNOW      , ENH_BUCKET_INTO_SNOW        , 
     I   ENG_SNOW_BUCKET , H2O_SNOW_BUCKET , 
     I   U_REF       , V_REF       , TMP_REF         , 
     I   RHO_SNOW_INV_PROG   ,
     O   ENG_PRECIPI , H2O_PRECIPI , RHO_PRECIPI_INV , 
     O   ENG_SNOW    , H2O_SNOW    , RHO_SNOW_INV    ,          
     O   IFLG_SLEV_CHANGE      )                      ! クリア & 雪無しに値
!
      CALL SNOW_AGE_INTG (
     I   NMAX   , NMAX_1 , DELT_SNOW , 
     I   RHO_SNOW_INV , 
     M   AGE_SNOW_PROG ) 
!
      I_SNOW_LEV(1         :NMAX_0(3)) = 0  
      I_SNOW_LEV(NMAX_01(2):NMAX_1(3)) = 1  
!
! SX6 のバグ対策
#ifdef SX6_OPTION
      do i=nmax_2(2),nmax_2(3)
        I_SNOW_LEV(i) = 2  
      ENDDO
#else
      I_SNOW_LEV(NMAX_2 (2):NMAX_2(3)) = 2  
#endif
      I_SNOW_LEV(NMAX_3 (2):NMAX_3(3)) = 3  
!
#ifdef SIB_DEBUG
      DO I=1,NMAX
        H2O_SNOW_SUM_TMP(I) = 
!    1       H2O_SNOW_BUCKET(I)         ! すでに precipi に加算
     1     + H2O_SNOW       (I,1) 
     1     + H2O_SNOW       (I,2) 
     1     + H2O_SNOW       (I,3) 
     1     + H2O_SNOW       (I,4) 
     1     + H2O_PRECIPI    (I) 
     1     + H2O_BUCKET_INTO_SOIL(I)    
        ENG_SNOW_SUM_TMP(I) = 
!    1       ENG_SNOW_BUCKET(I)         ! すでに precipi に加算
     1     + ENG_SNOW       (I,1) 
     1     + ENG_SNOW       (I,2) 
     1     + ENG_SNOW       (I,3) 
     1     + ENG_SNOW       (I,4) 
     1     + ENG_PRECIPI    (I) 
     1     + ENG_BUCKET_SNOW_TO_SNSL (I) 
     1     + ENG_BUCKET_INTO_SOIL(I)    
      ENDDO
#endif
!
!  ---------------------------------------
!  > 雪を解く ( 雪をリセットしない場合 ) < 
!  ---------------------------------------
!
      IF ( ICN_SIB0109_SNOW_RESET.EQ.0 )  THEN 
!
        CALL SNOW_WFLUX (
     I     NMAX      , 
     I     NMAX_01   , NMAX_2    , NMAX_3    , 
     I     H2O_PRECIPI  , ENG_PRECIPI   , RHO_PRECIPI_INV       , 
     I     H2O_SNOW     , ENG_SNOW      , RHO_SNOW_INV          , 
     I     I_SNOW_LEV   , 
     M     IFLG_SLEV_CHANGE      , 
     M     ENG_SNOW_NEXT         , H2O_SNOW_NEXT         , 
     M     RHO_SNOW_INV_PROG     , 
     M     ENH_BUCKET_SNOW_TO_SNSL  , WTR_BUCKET_SNOW_TO_SNSL  , 
     O     INFO_SNOW_NEXT )         ! 初期値 -1 、雪層数 ( 01 なら 1 ) 
!
      ENDIF
!
!  ------------------------------------------------
!  > 雪を解く ～ 雪の層のとり方をリセットする場合 < 
!  ------------------------------------------------
!    雪無しで雪のタネがある場合もここで処理する。
!
      CALL SNOW_SLEV_MOD (
     I    NMAX         , 
     I    ENG_PRECIPI  , H2O_PRECIPI  , RHO_PRECIPI_INV       , 
     I    ENG_SNOW     , H2O_SNOW     , RHO_SNOW_INV          , 
     M    ENG_SNOW_NEXT, H2O_SNOW_NEXT, RHO_SNOW_INV_PROG     ,  
     M    ENH_BUCKET_SNOW_TO_SNSL , WTR_BUCKET_SNOW_TO_SNSL   , 
     M    IFLG_SLEV_CHANGE            , 
     M    INFO_SNOW_NEXT  )                 ! 雪層数を代入。
!
! ======================
! >> 雪後処理ここから <<  
! ======================
!
!  ------------------
!  > 雪下面を融かす <
!  ------------------
!     土壌のエネルギーが正の場合、そのエネルギーをゼロにして、雪を融かす。
!
      CALL SNOW_BOTTOM_MELT (
     I    NMAX          , NMAX_0        , NMAX_01       , 
     I    CVR_SNW       , ENG_SOIL_TOP  , 
     M    INFO_SNOW_NEXT, 
     M    ENG_SNOW_NEXT , H2O_SNOW_NEXT , 
     M    ENG_BUCKET_INTO_SOIL    , ENG_BUCKET_SNOW_TO_SNSL , 
     M    H2O_BUCKET_INTO_SOIL    , WTR_BUCKET_SNOW_TO_SNSL  )  
!
!  --------------------
!  > 少雪バケツに格納 <
!  --------------------
!   ・現状、雪が少なければ少雪バケツへ。
!   ・以下の作業は今後の課題
!     スキンと合わせた温度が高い ( TF 以上 ) 場合、融雪水は土壌へ
!
      CALL SNOW_LITTLE_SNOW (
     I    NMAX                  ,
     M    INFO_SNOW_NEXT        ,                      ! 少雪なら 1 を -2 に
     M    H2O_SNOW_NEXT         , ENG_SNOW_NEXT    , 
     M    AGE_SNOW_PROG         , 
     O    ENG_SNOW_BUCKET_NEXT  , H2O_SNOW_BUCKET_NEXT   ) 
!
! ------------
! > 氷河流出 <
! ------------
!
#ifdef HIGH_SPEED
          BERG    (:) = 0. 
          ENH_BERG(:) = 0. 
#else
      CALL SNOW_BERG (
     I    NMAX                  ,
     M    H2O_SNOW_NEXT         , ENG_SNOW_NEXT    , 
     O    BERG                  , ENH_BERG           )
#endif
!
!  -------------------
!  > 水関係あと始末 <
!  -------------------
!    CVR_SNW_NEXT の計算
!    温度、水・氷の量を時間積分
!    氷河流出
!
      CALL SNOW_WFLUX_POST ( 
     I    NMAX             , 
     I    H2O_SNOW_NEXT    , ENG_SNOW_NEXT    , 
#ifdef MONYOS
     I    DELT_CNP         , DELT_SNOW        ,
#endif
     M    INFO_SNOW_NEXT   ,                         ! 01 なら 0 に
     O    CVR_SNW_NEXT     , CVR_SNW_INV_NEXT , 
     O    CVR_NOS_NEXT     , CVR_NOS_INV_NEXT , 
     O    TMP_SNOW_PROG    , WTR_SNOW_PROG    , GLA_SNOW_PROG   , 
     O    LIST_SNOW        , ICNT_SNOW )
!
! ------------------
! > 雪の深さの計算 <
! ------------------
!   GLA_SNOW_NEXT は雪無し域を含む格子平均なので、深さを求めるには
!   CVR_SNW_INV_NEXT をかける必要がある。
!
      CALL SNOW_POST (
     I    NMAX          , 
     I    GLA_SNOW_PROG , CVR_SNW_INV_NEXT , 
     I    RHO_SNOW_INV_PROG , 
     O    DEP_SNOW_NEXT  )
!
! ----------------------------------------------
! > モニタ用 ～ 次ステップでの雪被覆率等の設定 <
! ----------------------------------------------
!
#ifdef MONYOS
      CALL SIB_MONIT_RUN_SET_SNOWNEXTINFO (
     I   NMAX    , 
     I   CVR_SNW_NEXT , CVR_NOS_NEXT , INFO_SNOW_NEXT )
#endif
!
! ------------------
! > 雪の温度モニタ <
! ------------------
!   固定した深さでの温度モニタ。あくまで目安。温度は線形内挿。
!
      CALL SNOW_MONIT_FIXED_DEPTH (
     I    NMAX   , NMAX_01 , DELT_SNOW ,
     I    CVR_SNW_INV_NEXT  , 
     I    INFO_SNOW_NEXT    , 
     I    TMP_SNOW_PROG     , TMP_SNSK_NEXT     , 
     I    WTR_SNOW_PROG     , GLA_SNOW_PROG     , 
     I    DEP_SNOW_NEXT  )
!
#ifdef SIB_DEBUG
      DO I=1,NMAX
        H2O_SNOW_SUM_NEW(I) = 
     1       H2O_SNOW_BUCKET_NEXT(I)           ! すでに precipi に加算
     1     + H2O_SNOW_NEXT       (I,1) 
     1     + H2O_SNOW_NEXT       (I,2) 
     1     + H2O_SNOW_NEXT       (I,3) 
     1     + H2O_SNOW_NEXT       (I,4) 
     1     + WTR_BUCKET_SNOW_TO_SNSL (I) 
     1     + BERG                    (I) 
     1     + H2O_BUCKET_INTO_SOIL(I)    
!      
        ENG_SNOW_SUM_NEW (I) = 
     1       ENG_SNOW_BUCKET_NEXT(I)           ! すでに precipi に加算
     1     + ENG_SNOW_NEXT       (I,1) 
     1     + ENG_SNOW_NEXT       (I,2) 
     1     + ENG_SNOW_NEXT       (I,3) 
     1     + ENG_SNOW_NEXT       (I,4) 
     1     + ENH_BUCKET_SNOW_TO_SNSL (I)  
     1     + ENG_BUCKET_SNOW_TO_SNSL (I) 
     1     + ENH_BERG                (I) 
     1     + ENG_BUCKET_INTO_SOIL(I)    
      ENDDO
!
      H2O_ERR_MAX = -1
      ENG_ERR_MAX = -1
      DO I=1,NMAX
        IF ( 
     1         ABS(H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I)).GT.1.D-9 
     1     ) THEN
          write(6,*) 'SNOW_MAIN : ERROR HEATCALC H2O' , I , 
     1         H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I) , 
     1         H2O_SNOW_SUM_SAVE(I),H2O_SNOW_SUM_TMP(I)  
        ENDIF 
! 
        IF ( 
     1        ABS(ENG_SNOW_SUM_TMP(I)-ENG_SNOW_SUM_SAVE(I)).GT.1.D-6 
     1     ) THEN
          write(6,*) 'SNOW_MAIN : ERROR HEATCALC ENG' , I , 
     1         ENG_SNOW_SUM_TMP(I)-ENG_SNOW_SUM_SAVE(I) , 
     1         ENG_SNOW_SUM_SAVE(I),ENG_SNOW_SUM_TMP(I)  
!
          write(6,*) 'E_SAVE ' , 
     1       SRAD_BUCKET_SNW(I) ,
     1     + ENG_BUCKET_SNW (I) ,
     1     + ENG_SNOW_BUCKET(I) ,
     1     + ENH_BUCKET_INTO_SNOW(I) ,
     1     + ENH_BUCKET_FROM_SNOW(I) ,
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,1) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,1) )
     1           * ( TMP_SNOW_SAVE(I,1) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,1) ,
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,2) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,2) )
     1           * ( TMP_SNOW_SAVE(I,2) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,2) ,
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,3) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,3) )
     1           * ( TMP_SNOW_SAVE(I,3) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,3) ,
     1     + (   CAP_WTR_c * WTR_SNOW_SAVE(I,4) 
     1         + CAP_GLA_c * GLA_SNOW_SAVE(I,4) )
     1           * ( TMP_SNOW_SAVE(I,4) - TMP_FREZ_c  ) 
     1     - HL_W2G_c * GLA_SNOW_SAVE(I,4) 
!
          write(6,*) 'E_TMP ' , 
     1     + ENG_SNOW       (I,1) ,
     1     + ENG_SNOW       (I,2) ,
     1     + ENG_SNOW       (I,3) ,
     1     + ENG_SNOW       (I,4) ,
     1     + ENG_PRECIPI    (I) ,
     1     + ENG_BUCKET_SNOW_TO_SNSL (I) 
!
        ENDIF 
!
        IF ( ABS(H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I))
     1          .GT.H2O_ERR_MAX) THEN
          H2O_ERR_MAX = ABS(H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I))
        ENDIF
        IF ( ABS(ENG_SNOW_SUM_TMP(I)-ENG_SNOW_SUM_SAVE(I))
     1          .GT.ENG_ERR_MAX) THEN
          ENG_ERR_MAX = ABS(ENG_SNOW_SUM_TMP(I)-ENG_SNOW_SUM_SAVE(I))
        ENDIF
      ENDDO

      write(6,*) 'SNOW_MAIN : water  h2o/eng heat diff ' , 
     1            h2o_err_max , eng_err_max 
! ---------------------------------------------------------------
      H2O_ERR_MAX = -1
      ENG_ERR_MAX = -1
      DO I=1,NMAX
        IF ( 
     1         ABS(H2O_SNOW_SUM_NEW(I)-H2O_SNOW_SUM_TMP(I)).GT.1.D-9 
     1     ) THEN
          write(6,*) 'SNOW_MAIN : ERROR WATERCALC MASS' , I , 
     1               IFLG_SLEV_CHANGE (I)  , 
     1         H2O_SNOW_SUM_NEW(I)-H2O_SNOW_SUM_TMP(I) , 
     1         H2O_SNOW_SUM_TMP(I),H2O_SNOW_SUM_NEW(I)  
!
          write(6,*) ' H2O_TMP' ,
     1     + H2O_SNOW       (I,1) ,
     1     + H2O_SNOW       (I,2) ,
     1     + H2O_SNOW       (I,3) ,
     1     + H2O_SNOW       (I,4) ,
     1     + H2O_PRECIPI    (I) 

          write(6,*) ' H2O_NEW' ,
     1       H2O_SNOW_NEXT       (I,1)     ,
     1       H2O_SNOW_NEXT       (I,2)     ,
     1       H2O_SNOW_NEXT       (I,3)     ,
     1       H2O_SNOW_NEXT       (I,4)     ,
     1       H2O_SNOW_BUCKET_NEXT(I)       ,
     1       WTR_BUCKET_SNOW_TO_SNSL (I)   , 
     1       BERG (I)    
!
!
        ENDIF 
        IF ( 
     1    ABS(ENG_SNOW_SUM_NEW(I)-ENG_SNOW_SUM_TMP(I)).GT.1.D-6 
     1     ) THEN
          write(6,*) 'SNOW_MAIN : ERROR WATERCALC ENG ' , I , 
     1               IFLG_SLEV_CHANGE (I)  , 
     1         ENG_SNOW_SUM_NEW(I)-ENG_SNOW_SUM_TMP(I) , 
     1         ENG_SNOW_SUM_TMP(I),ENG_SNOW_SUM_NEW(I)
!
!
        ENDIF 
!
        IF ( ABS(H2O_SNOW_SUM_NEW(I)-H2O_SNOW_SUM_TMP(I))
     1          .GT.H2O_ERR_MAX) THEN
          H2O_ERR_MAX = ABS(H2O_SNOW_SUM_NEW(I)-H2O_SNOW_SUM_TMP(I))
        ENDIF
        IF ( ABS(ENG_SNOW_SUM_NEW(I)-ENG_SNOW_SUM_TMP(I))
     1          .GT.ENG_ERR_MAX) THEN
          ENG_ERR_MAX = ABS(ENG_SNOW_SUM_NEW(I)-ENG_SNOW_SUM_TMP(I))
        ENDIF
      ENDDO
!
      write(6,*) 'snow_main : water  h2o/eng diff ' , 
     1            h2o_err_max , eng_err_max 
!
      IF ( 
     1            h2o_err_max .gt. 1 
     1      .or.  eng_err_max .gt. 1 
     1     ) THEN
            write(6,*) 'snow_main ERROR too large' 
            stop 999
      ENDIF


#endif
!
      RETURN
      END

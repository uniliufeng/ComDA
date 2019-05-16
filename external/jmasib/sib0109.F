      MODULE SIB0109
!
      use com_debug , only :
     1   ipoint_ij    , 
     1   ipoint_new   ,
     1   ipmask_new
!
      USE PRM , ONLY : 
     1   IJPHY , JLPHY , ISPT , INTRI , INTRJ , KMAX  
!
      USE SIBPRM , ONLY : 
     1    IDX , NUM_BLOCK, IDP , ISN , ICG , IVN , IBD , IWV 
!
      USE SIBVEG , ONLY : SIBVEG_INI 
!
      IMPLICIT NONE
!
      PRIVATE  
!
!                                      このカンマを動かす長いサブルーチンは
!                                      ↓ 規格非適合になります
      PUBLIC 
     1  SIB0109_INI                 , 
     1  SIB0109_RUN_MONITOR_INITIAL ,
     1  SIB0109_RUN_STEPINI         ,
     1  SIB0109_RUN_LOOPINI         ,
     1  SIB0109_RUN_ALBEDO          , 
     1  SIB0109_RUN_SR_CALC         ,
     1  SIB0109_RUN_SR_RESTORE      , 
     1  SIB0109_RUN_LR_SET          , 
     1  SIB0109_RUN_SIBMAIN         ,
     1  SIB0109_RUN_LOOPEND         ,
     1  SIB0109_RUN_STEPEND
!
! 配列について
!   ・mj98 の物理過程では (IJPHY,JLPHY) という切り方をしている。
!     DO JL=1,JLPHY のループで計算され、各 JL のことをバンドと呼ぶ。
!
! ==================
! >> SiB 予報変数 <<
! ==================
!  リスタートファイルに書き込まれる必要がある。
!
      REAL(8),SAVE :: TMP_CNP_NOS_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),SAVE :: FR_WTR_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),SAVE :: FR_GLA_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),SAVE :: TMP_CNP_SNW_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),SAVE :: FR_WTR_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),SAVE :: FR_GLA_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),SAVE :: TMP_GRSK_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),SAVE :: FR_WTR_GRS_ALL     (ISPT*IJPHY,JLPHY)      ! 地面水割合
      REAL(8),SAVE :: FR_GLA_GRS_ALL     (ISPT*IJPHY,JLPHY)      ! 地面氷割合
      REAL(8),SAVE :: TMP_SNSK_ALL       (ISPT*IJPHY,JLPHY)
!
!  土壌
!
      REAL(8),SAVE :: TMP_SOIL_ALL       (ISPT*IJPHY,JLPHY,IDP) ! 
      REAL(8),SAVE :: FR_WTR_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP) ! 空隙中水割合
      REAL(8),SAVE :: FR_GLA_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP) ! 空隙中氷割合
      REAL(8),SAVE :: TMP_SNSL_ALL       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),SAVE :: FR_WTR_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),SAVE :: FR_GLA_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
!
      REAL(8),SAVE :: TMP_SOIL4_ALL       (ISPT*IJPHY,JLPHY) ! 
!
      REAL(8),SAVE :: ALB_SNOW_SNW_ALL       (ISPT*IJPHY,JLPHY,IVN) 
      REAL(8),SAVE :: RAIN_1HOUR_ALL         (ISPT*IJPHY,JLPHY) 

!
!  積雪
!
      REAL(8),SAVE :: TMP_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN) ! 温度
      REAL(8),SAVE :: WTR_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN) ! 水の量kg/m^2
      REAL(8),SAVE :: GLA_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN) ! 氷の量kg/m^2
      REAL(8),SAVE :: RHO_SNOW_INV_ALL   (ISPT*IJPHY,JLPHY,ISN) ! 密度kg/m^3   
!
      REAL(8),SAVE :: AGE_SNOW_ALL       (ISPT*IJPHY,JLPHY)     ! 年齢 sec 
      INTEGER,SAVE :: INFO_SNOW_ALL      (ISPT*IJPHY,JLPHY)     ! 
!
      REAL(8),SAVE :: ENG_SNOW_BUCKET_ALL(ISPT*IJPHY,JLPHY) ! 少雪 ENE バケツ
      REAL(8),SAVE :: H2O_SNOW_BUCKET_ALL(ISPT*IJPHY,JLPHY) ! 少雪 H2O バケツ
!     ... おなじ BUCKET でも、これだけは配列が異なることに注意
!
!  =========================
!  >> SiB SAVE される変数 << 
!  =========================
!    ・短波の計算ごとに変更される。
!    ・陸面の雪の計算よりも時間的に粗いことを想定するので、
!      (ISPT*IJPHY,JLPHY) という切り方をする必要がある。
!
      REAL(8),SAVE :: R_STM_BAR_PAR_NOS_SAVE  (ISPT*IJPHY,JLPHY,ICG)
      REAL(8),SAVE :: R_STM_BAR_PAR_SNW_SAVE  (ISPT*IJPHY,JLPHY,ICG)
      REAL(8),SAVE :: RAD_SHRT_NOS_SAVE       (ISPT*IJPHY,JLPHY,ICG)
      REAL(8),SAVE :: RAD_SHRT_SNW_SAVE       (ISPT*IJPHY,JLPHY,ICG)
!
!  ===============================
!  >> SiB SAVE されるバケツ変数 << 
!  ===============================
!     ( RUN を通じた引渡しの必要はなし。 ) 
!     注意。  
!       ・ここで
!           雪・土壌はおなじステップで計算されること、
!           その計算をしたのちに以下のバケツはクリアされること
!         を前提にしている。
!         # そうでないなら、(ISPT*IJPHY,JLPHY) であり、
!         # 雪層数変更に伴う格子順番の変更に対応する措置が必要。
!         # しかし上の前提により、それは必要でない。
!
      REAL(8),SAVE :: ENH_BUCKET_INTO_SOIL  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: ENH_BUCKET_INTO_SNOW  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: ENH_BUCKET_FROM_SOIL  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: ENH_BUCKET_FROM_SNOW  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: H2O_BUCKET_INTO_SOIL  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: H2O_BUCKET_INTO_SNOW  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: H2O_BUCKET_FROM_SOIL  (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: H2O_BUCKET_FROM_SNOW  (IDX,NUM_BLOCK) 

      REAL(8),SAVE :: ENG_BUCKET_NOS        (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: ENG_BUCKET_SNW        (IDX,NUM_BLOCK) 
!
      REAL(8),SAVE :: TRNS_BUCKET_NOS       (IDX,IDP,NUM_BLOCK)  
      REAL(8),SAVE :: TRNS_BUCKET_SNW       (IDX,IDP,NUM_BLOCK) 
      REAL(8),SAVE :: SRAD_BUCKET_SNW       (IDX,NUM_BLOCK) 
!
!
! ==========================
!  >> SiB SAVE される温度 <<
! ==========================
!   土壌・雪第一層温度
!   各ステップで更新される。
!   土壌・雪計算ステップ後にクリアされ、真の土壌・雪第一層の値になる
!   ジョブ間で引き継いだり、NMAX 並び順が変わることがないので、
!   実格子 (IJPHY*ISPT,NUM_BLOCK) には詰めない
!
      REAL(8),SAVE :: TMP_SOIL1_TMP_ALL (IDX,NUM_BLOCK) 
      REAL(8),SAVE :: TMP_SNOW1_TMP_ALL (IDX,NUM_BLOCK) 
!
! ====================
! >> ジョブ管理変数 <<
! ====================
!
      INTEGER,SAVE  :: IMASK    ( ISPT*IJPHY, JLPHY ) ! 植生・海陸情報
      REAL(8),SAVE  :: RATIO_VEG( ISPT*IJPHY, JLPHY ) ! 植生・海陸情報
!
! ======================== 
! >> ステップによる変数 <<
! ======================== 
!
      REAL(8),SAVE  :: DELT_ATM      ! 大気の時間刻。リープフロッグなら 2 倍?
      REAL(8),SAVE  :: DELT_CNP      ! CANOPY の時間刻
      REAL(8),SAVE  :: DELT_SNOW     
      REAL(8),SAVE  :: DELT_SOIL
!
! ======================
! >> 以下、バンド変数 <<
! ======================
!
! -------------------------------
! > 大気情報 ( SiB バンド変数 ) <
! -------------------------------
!
!   入力 （大気から）
!
      REAL(8),SAVE :: U_REF     (IDX) 
      REAL(8),SAVE :: V_REF     (IDX) 
      REAL(8),SAVE :: Z_REF     (IDX) 
      REAL(8),SAVE :: P_REF     (IDX)                  ! 
      REAL(8),SAVE :: TMP_REF   (IDX)                  ! 
      REAL(8),SAVE :: HUM_REF   (IDX)     ! 大気最下層比湿
      REAL(8),SAVE :: FAC_CORTMP(IDX)     ! 温度補正項。BPS 
      REAL(8),SAVE :: DAYTIME_1HR  (IDX)
      REAL(8),SAVE :: ZMEAN     (IDX)
      REAL(8),SAVE :: ZTEMP     (IDX)
      REAL(8),SAVE :: GMT       (IDX,2)
      REAL(8),SAVE :: GMQ       (IDX,2)
      REAL(8),SAVE :: GMUV      (IDX,3)
      REAL(8),SAVE :: P_DEL     (IDX)      ! 大気最下層圧
      REAL(8),SAVE :: RAIN_LARGE(IDX)
      REAL(8),SAVE :: RAIN_CONV (IDX)
!
!   出力 （大気へ）
!
      REAL(8),SAVE :: GMT_OUT   (IDX,2)
      REAL(8),SAVE :: GMQ_OUT   (IDX,2)
      REAL(8),SAVE :: GMUV_OUT  (IDX,3)
!     REAL(8),SAVE :: U_FLUX    (IDX)      ! 運動量フラックス
!     REAL(8),SAVE :: V_FLUX    (IDX)
!
! ---------------------------------------
! > 予報変数を SiB バンド格子化したもの <
! ---------------------------------------
!
      REAL(8),SAVE :: TMP_CNP_NOS_PROG    (IDX)
      REAL(8),SAVE :: FR_WTR_CNP_NOS_PROG (IDX)    
      REAL(8),SAVE :: FR_GLA_CNP_NOS_PROG (IDX)    
      REAL(8),SAVE :: TMP_CNP_SNW_PROG    (IDX)
      REAL(8),SAVE :: FR_WTR_CNP_SNW_PROG (IDX)    
      REAL(8),SAVE :: FR_GLA_CNP_SNW_PROG (IDX)    
      REAL(8),SAVE :: TMP_GRSK_PROG       (IDX)
      REAL(8),SAVE :: FR_WTR_GRS_PROG     (IDX)      ! 地面水割合
      REAL(8),SAVE :: FR_GLA_GRS_PROG     (IDX)      ! 地面氷割合
      REAL(8),SAVE :: TMP_SNSK_PROG       (IDX)
!
!
!     REAL(8),SAVE :: TMP_SOIL1_TMP       (IDX)  ! 土壌第一層温度 ( canopy 用 )
!     REAL(8),SAVE :: TMP_SNOW1_TMP       (IDX)  ! 積雪第一層温度 ( canopy 用 )
!
      REAL(8),SAVE :: TMP_SOIL4_PROG      (IDX)  ! 土壌第 4 層温度 
      REAL(8),SAVE :: ALB_SNOW_SNW_PROG   (IDX,IVN) 
      REAL(8),SAVE :: RAIN_1HOUR_PROG     (IDX) 
!
      REAL(8),SAVE :: TMP_SOIL_PROG    (IDX,IDP) ! 
      REAL(8),SAVE :: FR_WTR_SOIL_PROG (IDX,IDP) ! 空隙中の水の割合
      REAL(8),SAVE :: FR_GLA_SOIL_PROG (IDX,IDP) ! 空隙中の氷の割合
      REAL(8),SAVE :: TMP_SNSL_PROG    (IDX,IDP) ! 
      REAL(8),SAVE :: FR_WTR_SNSL_PROG (IDX,IDP) ! 空隙中の水の割合
      REAL(8),SAVE :: FR_GLA_SNSL_PROG (IDX,IDP) ! 空隙中の氷の割合
!
      REAL(8),SAVE :: TMP_SNOW_PROG     (IDX,ISN+1)  ! 雪の温度。
!                                                    ! ただし i 層の場合、
!                                                    ! i+1 には TMP_SNSL(1) 
      REAL(8),SAVE :: WTR_SNOW_PROG     (IDX,ISN)
      REAL(8),SAVE :: GLA_SNOW_PROG     (IDX,ISN)
      REAL(8),SAVE :: RHO_SNOW_INV_PROG (IDX,ISN)    
      REAL(8),SAVE :: AGE_SNOW_PROG     (IDX)
      INTEGER,SAVE :: INFO_SNOW      (IDX)          ! 雪層数 0,10,1,2,3,4 

      REAL(8),SAVE :: CVR_SNW        (IDX) 
      REAL(8),SAVE :: CVR_SNW_NEXT   (IDX) 
!
      INTEGER,SAVE :: INFO_SNOW_NEXT      (IDX) 
!
      REAL(8),SAVE  :: ENG_SNOW_BUCKET     (IDX)
      REAL(8),SAVE  :: H2O_SNOW_BUCKET     (IDX)
      REAL(8),SAVE  :: ENG_SNOW_BUCKET_NEXT(IDX)
      REAL(8),SAVE  :: H2O_SNOW_BUCKET_NEXT(IDX)
!
! ----------
! > 雪情報 <
! ----------
!   いろいろなところで参照される
!
      REAL(8),SAVE :: CVR_VEG_ms         (IDX,ICG)    
      REAL(8),SAVE :: CVR_VEG_INV_ms     (IDX,ICG)    
      REAL(8),SAVE :: CVR_NOS            (IDX)               ! 
      REAL(8),SAVE :: CVR_SNW_INV        (IDX)               ! 
      REAL(8),SAVE :: CVR_NOS_INV        (IDX)               ! 
      REAL(8),SAVE :: DEP_SNOW           (IDX,ISN+1)         ! 
!
      REAL(8),SAVE :: DEP_SNOW_NEXT(IDX,ISN+1)

!
      REAL(8),SAVE :: XLAI_VEG_NOS       (IDX,ICG)
      REAL(8),SAVE :: XLAI_VEG_SNW       (IDX,ICG)
!


!
! ------------
! > 放射関連 <
! ------------
!
!  ・ALBEDO 計算時に計算され、RAD2SIB で使用されるもの
!    短波計算時のみ使われるので、(IDX) でよいことに注意。 
!  
      REAL(8),SAVE :: EXTK_NOS        (IDX,ICG,IVN,IBD)  ! 消散係数 
      REAL(8),SAVE :: ABSORB_SHRT_NOS (IDX,ICG,IVN,IBD)  ! 地上下向短波率
      REAL(8),SAVE :: REACH_GRD_NOS   (IDX,IBD)          ! 到達率
!
      REAL(8),SAVE :: EXTK_SNW        (IDX,ICG,IVN,IBD)  ! 消散係数 
      REAL(8),SAVE :: ABSORB_SHRT_SNW (IDX,ICG,IVN,IBD)  ! 地上下向短波率
      REAL(8),SAVE :: REACH_GRD_SNW   (IDX,IBD)          ! 到達率
!
!  ・RAD2SIB で作成され、CANOPY 等で使用されるもの、その 1 
!     *_SAVE にセーブされ、雪計算後に変更される!!! 
!
      REAL(8),SAVE :: R_STM_BAR_PAR_NOS  (IDX,ICG)       ! 気孔抵抗PAR依存性 
      REAL(8),SAVE :: R_STM_BAR_PAR_SNW  (IDX,ICG)       ! 
!
      REAL(8),SAVE :: RAD_SHRT_NOS       (IDX,ICG)       ! SiB に入る短波
      REAL(8),SAVE :: RAD_SHRT_SNW       (IDX,ICG)    
!
      REAL(8),SAVE :: RAD_LONG_NOS       (IDX,ICG)       ! SiB に入る長波
      REAL(8),SAVE :: RAD_LONG_SNW       (IDX,ICG)    
!
      REAL(8),SAVE :: RAD_LONG_SIB_2_ATM(IDX)           ! (下向き正)
      REAL(8),SAVE :: RAD_LONG_NOS_2_ATM(IDX)           ! 
      REAL(8),SAVE :: RAD_LONG_SNW_2_ATM(IDX)           ! 
!
!  ・RAD2SIB で作成され、CANOPY 等で使用されるもの、その 2 
!     毎ステップ新規に作成される。
!
      REAL(8),SAVE :: RAD_LONG_INTO_SIB  (IDX)           ! SiB に入る長波
!
      REAL(8),SAVE :: RAD_SHRT_INTO_SNOW_LAYER  (IDX) 

!
      REAL(8),SAVE :: THERMK_NOS         (IDX)           ! キャノピー射出率
      REAL(8),SAVE :: THERMK_SNW         (IDX)
!
! ------------------
! > バンド格子情報 <
! ------------------
!
      INTEGER,SAVE :: ILIST            (IDX)  !  IP + (IJ-1) * ISPT 
      INTEGER,SAVE :: IMSKLT           (IDX)  !  SiB -> 物理
      INTEGER,SAVE :: IMSKLT_RX        (IDX)  !  SiB -> 放射X
      INTEGER,SAVE :: IMSKLT_RY        (IDX)  !  SiB -> 放射Y
!
      INTEGER,SAVE :: NMAX                     
      INTEGER,SAVE :: NMAX_0 (3)
      INTEGER,SAVE :: NMAX_01(3)
      INTEGER,SAVE :: NMAX_1 (3)
      INTEGER,SAVE :: NMAX_2 (3)
      INTEGER,SAVE :: NMAX_3 (3)
      INTEGER,SAVE :: NMAX_4 (3)
!
      INTEGER,SAVE :: NTYP   (IDX)              ! 現状では IDX=IJPHY 
      INTEGER,SAVE :: MON_NEXT 
!
! 収支チェック
!
      REAL(8),SAVE   :: H2O_CNP_NOS_1         (IDX) 
      REAL(8),SAVE   :: H2O_CNP_SNW_1         (IDX) 
      REAL(8),SAVE   :: H2O_GRSK_1            (IDX) 
      REAL(8),SAVE   :: ENG_CNP_NOS_1         (IDX) 
      REAL(8),SAVE   :: ENG_CNP_SNW_1         (IDX) 
      REAL(8),SAVE   :: ENG_GRSK_1            (IDX) 
      REAL(8),SAVE   :: ENG_SNSK_1            (IDX) 
!
      REAL(8),SAVE   :: H2O_CNP_NOS_2         (IDX) 
      REAL(8),SAVE   :: H2O_CNP_SNW_2         (IDX) 
      REAL(8),SAVE   :: H2O_GRSK_2            (IDX) 
      REAL(8),SAVE   :: ENG_CNP_NOS_2         (IDX) 
      REAL(8),SAVE   :: ENG_CNP_SNW_2         (IDX) 
      REAL(8),SAVE   :: ENG_GRSK_2            (IDX) 
      REAL(8),SAVE   :: ENG_SNSK_2            (IDX) 

      
!
      CONTAINS 
!
! =================================================================
      SUBROUTINE SIB0109_INI  ( 
#ifdef CORRECT
     I    RATIO_VEG_IN ,
#endif
     I    IMASK_IN ,  ID_NEXT ) 
!
!  初期処理。
!    SiB に必要な諸々の情報の設定
!
      USE SIBCON , ONLY :
     1   SIBCON_INI 
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1    CCN_FILENAME_VEG
!
      USE NMAX_INFO , ONLY :
     1    NMAX_INFO_INI  
!
      USE SIB_MONIT , ONLY : SIB_MONIT_INI
!
      IMPLICIT NONE 
!
#ifdef CORRECT
      REAL(8),INTENT(IN) :: RATIO_VEG_IN  ( ISPT*IJPHY,JLPHY ) ! 
#endif
      INTEGER,INTENT(IN) :: IMASK_IN  ( ISPT*IJPHY,JLPHY ) ! 植生・海陸情報
!     REAL(8),INTENT(IN) :: FSECP
!     INTEGER,INTENT(IN) :: KTP 
      INTEGER,INTENT(IN) :: ID_NEXT(5)
!
      INTEGER :: IJ
      INTEGER :: JL
!
! NEWSIB だよ!!!
!
      write(6,*) 'sib0109_ini start' 
!
! -----------------
! > SiB 関連 定数 <
! -----------------
!
!
!  SiB 定数、物理定数等
!
      CALL SIBCON_INI
!
!
!  植生データ読み込み
!
      OPEN ( 90 , FILE=CCN_FILENAME_VEG , form='unformatted' )
      CALL SIBVEG_INI ( 90 )               ! 12 は装置番号
!
      CALL SIB_MONIT_INI
!
! ------------------
! > バケツの初期化 <
! ------------------
!
      call reset ( ENH_BUCKET_INTO_SOIL    , idx*num_block )
      call reset ( ENH_BUCKET_INTO_SNOW    , idx*num_block )
      call reset ( ENH_BUCKET_FROM_SOIL    , idx*num_block )
      call reset ( ENH_BUCKET_FROM_SNOW    , idx*num_block )
!
      call reset ( H2O_BUCKET_INTO_SOIL    , idx*num_block )
      call reset ( H2O_BUCKET_INTO_SNOW    , idx*num_block )
      call reset ( H2O_BUCKET_FROM_SOIL    , idx*num_block )
      call reset ( H2O_BUCKET_FROM_SNOW    , idx*num_block )
!
      call reset ( TRNS_BUCKET_NOS         , idx*num_block*idp )
      call reset ( TRNS_BUCKET_SNW         , idx*num_block*idp )
      call reset ( ENG_BUCKET_NOS          , idx*num_block )
      call reset ( ENG_BUCKET_SNW          , idx*num_block )
      call reset ( SRAD_BUCKET_SNW         , idx*num_block )
!
! ----------------
! > IMASK の設定 <
! ----------------
!
      IMASK(:,:) = IMASK_IN(:,:)
#ifdef CORRECT
        RATIO_VEG(:,:) = RATIO_VEG_IN(:,:)
#else
      RATIO_VEG(:,:) = 0. 
      DO JL=1,JLPHY
      DO IJ=1,IJPHY
        RATIO_VEG(1+(IJ-1)*ISPT,JL)= 1.D0  
      ENDDO
      ENDDO
#endif
!      
! -----------------
! > NMAX 等の設定 <   ( 雪のない場合 ) 
! -----------------
!
      CALL NMAX_INFO_INI ( IMASK )           
!
! ------------------------
! > SiB 全面初期値の設定 <
! ------------------------
!    初期値なら、（現状）適当に作る 
!    リスタートファイルから読む。共有変数に格納
!
      write(6,*) 'sib0109 before SIB0109_FORCAST_INITIAL_SET' , 
!    1           ' KTP FSECP' , KTP , FSECP
     1           ' KTP FSECP' , ID_NEXT 
!
      CALL SIB0109_FCST_INITIAL_SET (
!    I     KTP , FSECP )
     I     ID_NEXT )
!
      RETURN
      END SUBROUTINE SIB0109_INI
! =================================================================
      SUBROUTINE SIB0109_FCST_INITIAL_SET ( ID_NEXT  ) 
!
!  初期値の設定。
!  このサブルーチンは private である。
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_SIB_ISLSCP 

      USE COM_JOBINFO_SIB0109 , ONLY : 
     1    IUNIT_RESTART_IN   ,
     1    IDSTAR             , 
     1    CDIR_RESTART
!
      use restart_sib0109 , only : READ_RESTART_SIB0109

!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: ID_NEXT(5)
      integer :: i,j
      character(100) :: cfile_restart_in 
!
!
!  ---------------------------------------------------
!  > ファイルから読み込まない場合 ( 初期値スタート ) <
!  ---------------------------------------------------
!
      IF     ( IUNIT_RESTART_IN .LE. 0 ) THEN 
        WRITE(6,*) 'INITIAL DATA SETTING FROM ATMOSPHERIC DATA ETC' 
!
        CALL MAKE_SIB_INI (
     I   IMASK              , ID_NEXT            , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
        do j=1,jlphy
        do i=1,ijphy
          if ( info_snow_all(i,j) .le. 0 ) then 
            tmp_soil4_all(i,j) = tmp_soil_all(i,j,3)
          else
            tmp_soil4_all(i,j) = tmp_snsl_all(i,j,3)
          endif
        enddo 
        enddo 
!
!  ----------------------------
!  > ファイルから読み込む場合 <
!  ----------------------------
!
      ELSEIF ( IUNIT_RESTART_IN .GE. 1 ) THEN 
!
        CFILE_RESTART_IN  = TRIM(CDIR_RESTART) // 
     1                      'newsib_restart_YYYY_MM_DD_HH'
        CALL REPLACE_INT(CFILE_RESTART_IN, 'YYYY', IDSTAR(1))
        CALL REPLACE_INT(CFILE_RESTART_IN, 'MM'  , IDSTAR(2))
        CALL REPLACE_INT(CFILE_RESTART_IN, 'DD'  , IDSTAR(3))
        CALL REPLACE_INT(CFILE_RESTART_IN, 'HH'  , IDSTAR(4))
!
        WRITE(6,*) 'Input restart file is ' ,  
     1      CFILE_RESTART_IN , ' unit=' , IUNIT_RESTART_IN
!
        open ( iunit_restart_in , file=cfile_restart_in , 
     1         form='unformatted' )
!
         CALL READ_RESTART_SIB0109 ( 
     I   IUNIT_RESTART_IN   , 
     I   ID_NEXT            , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL,
     O   TMP_SOIL4_ALL      , 
     O   ALB_SNOW_SNW_ALL   , RAIN_1HOUR_ALL    ) 
!
        close ( iunit_restart_in ) 
!
      ENDIF
!
      RETURN
      END SUBROUTINE SIB0109_FCST_INITIAL_SET   
! =================================================================
      SUBROUTINE SIB0109_RUN_MONITOR_INITIAL ( MON , KT0 )   
!
! このサブルーチンは、オンラインではまず呼ばない。
! オフライン専用
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1    JCN_SIB_ISLSCP ,
     1    RCN_DELT_ATM  ,
     1    RCN_DELT_CNP  ,
     1    RCN_DELT_SOIL ,
     1    RCN_DELT_SNOW 
!
      USE SIBCON , ONLY : 
     1   H2O_SAT_CNP_DEP_c , 
     1   H2O_SAT_GRS_DEP_c 

!
      USE COM_STEP , ONLY :      ! io-yoshimura モニタをだますため。
     1   ICNMNTMON               ! 例外的に COM_STEP を使う
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR 
                                 ! 同上 monit_snap_* をだます。
!
      USE SIB_MONIT , ONLY :
     1   SIB_MONIT_RUN_SET_SNOWINFO          ,
     1   SIB_MONIT_RUN_SET_SNOWNEXTINFO    ! モニタ用。
!
      USE NMAX_INFO , ONLY : 
     1   NMAX_INFO_RUN_SNOW 
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN) :: MON      
      INTEGER,INTENT(IN) :: KT0 
      INTEGER            :: NB       ! ブロックラベル
      INTEGER            :: JL       ! ブロックラベル
      INTEGER            :: I
!
! ======================================
! >> 初期値スナップショット作成時のみ <<
! ======================================
!
        DO 6000 NB=1,NUM_BLOCK
!
          WRITE(6,*) 'NB =' , NB 
          JL = NB 
!
! テーブルなど作成
!
         CALL NMAX_INFO_RUN_SNOW (
     I   IMASK(1,NB), INFO_SNOW_ALL(1,NB) , NB         , 
     O   NMAX       , NTYP       ,     
     O   NMAX_0     , NMAX_01    , NMAX_1     , 
     O   NMAX_2     , NMAX_3     , NMAX_4     ,     
     O   IMSKLT     , IMSKLT_RX  , IMSKLT_RY  , ILIST )
!
#ifdef MONYOS
      CALL SIB_MONIT_RUN_SET_SNOWINFO (
     I   NMAX    , 
     I   CVR_SNW , CVR_NOS , IMSKLT , NMAX_01 , NB )
#endif
!
! 全面 SiB 格子 -> SiB 格子
!             NB はしかるべく使用
!
          CALL SET_SIBGRD_SIB (
     I   JL      , ILIST    , 
     I   NMAX    , NMAX_01  , NMAX_1   , NMAX_2  , NMAX_3   , NMAX_4   ,     
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL, 
     I   TMP_SOIL4_ALL      , 
     I   ALB_SNOW_SNW_ALL   , RAIN_1HOUR_ALL     , 
     O   TMP_CNP_NOS_PROG   , FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     O   TMP_CNP_SNW_PROG   , FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG, 
     O   TMP_GRSK_PROG      , FR_WTR_GRS_PROG    , FR_GLA_GRS_PROG    , 
     O   TMP_SNSK_PROG      , 
     O   INFO_SNOW          , 
     O   TMP_SOIL_PROG      , FR_WTR_SOIL_PROG   , FR_GLA_SOIL_PROG   , 
     O   TMP_SNSL_PROG      , FR_WTR_SNSL_PROG   , FR_GLA_SNSL_PROG   ,  
     O   TMP_SNOW_PROG      , WTR_SNOW_PROG      , GLA_SNOW_PROG      , 
     O   RHO_SNOW_INV_PROG  , AGE_SNOW_PROG      , 
     O   ENG_SNOW_BUCKET    , H2O_SNOW_BUCKET    ,
     O   TMP_SOIL4_PROG     ,   
     O   ALB_SNOW_SNW_PROG  , RAIN_1HOUR_PROG )
!
! SiB 格子データ作成 （SiB 格子 -> SiB 格子）
!
          CALL MAKE_SIBGRD_SIB (  
     I  NMAX     , 
     I  NMAX_0   , NMAX_01 , NMAX_1  , NMAX_2   , NMAX_3  , NMAX_4   , 
     I  NTYP     , MON     ,
     I  WTR_SNOW_PROG      , GLA_SNOW_PROG    , 
     I  RHO_SNOW_INV_PROG  ,  
     O  CVR_VEG_ms      , CVR_VEG_INV_ms ,   
     O  CVR_SNW         , CVR_SNW_INV    , 
     O  CVR_NOS         , CVR_NOS_INV    , 
     O  DEP_SNOW      )  
!
          IF ( KT0.EQ.-1 .AND. JCN_SIB_ISLSCP.EQ.0 ) THEN
            CALL MOD_INI_DATA (
     I        NMAX_01 , NMAX , MON  ,  NTYP ,           
     M        FR_WTR_CNP_NOS_PROG , FR_WTR_GRS_PROG       , 
     M        FR_WTR_CNP_SNW_PROG , 
     M        FR_GLA_CNP_NOS_PROG , FR_GLA_GRS_PROG       , 
     M        FR_GLA_CNP_SNW_PROG  ) 
          ENDIF
!
          ICNMNTMON    = 1     ! io-yoshimura 中のサブルーチンをだます
          ICN_SIB0109_MONITOR = 1 
!
#ifdef MONYOS
          CALL SIB_MONIT_RUN_SET_SNOWNEXTINFO  (
     I   NMAX , 
     I   CVR_SNW , CVR_NOS , INFO_SNOW )
#endif
!
!
!  注意。000304 変更。
!    最初の瞬間値を出すとき。
!    平均値も出すが、
!    その際重みとしては、delt_soil ではなく delt_atm を使う。
!    monit_out の時に /delt_atm をするため。
!    そこで info_mean として -1 を代入。 
!
          CALL MONIT_SNAP_SNOW (
     I    NMAX        , 
     I    1           , -1              ,         ! imean , info_mean
     I    RCN_DELT_SNOW   , RCN_DELT_ATM        , 
     I    INFO_SNOW   , 
     I    CVR_SNW     , 
     I    WTR_SNOW_PROG   , GLA_SNOW_PROG  , 
     I    TMP_SNOW_PROG   , DEP_SNOW       , 
     I    H2O_SNOW_BUCKET , 
     I    ENG_SNOW_BUCKET , 
     I    RHO_SNOW_INV_PROG  
     1   )  
!
          CALL MONIT_SNAP_SOIL  (
     I    NMAX        ,  NTYP           , 
     I    1           ,  -1             ,           ! imean , info_mean
     I    RCN_DELT_SOIL   ,  RCN_DELT_ATM       ,
     I    CVR_SNW     , 
     I    TMP_SOIL_PROG    , TMP_SNSL_PROG        , TMP_SOIL4_PROG , 
     I    FR_WTR_SOIL_PROG , FR_WTR_SNSL_PROG     , 
     I    FR_GLA_SOIL_PROG , FR_GLA_SNSL_PROG     
     1   )  
!
          CALL MONIT_SNAP_CANOPY_2 (
     I    NMAX    , MON    , NTYP   , 1   , 
     I    RCN_DELT_CNP   , 
     I    CVR_SNW        , 
     I    DEP_SNOW       ,
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG ,  
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG      
     1     )  
!
          CALL MONIT_SNAP_CANOPY (
     I    NMAX    , MON    , NTYP   , 
     I    RCN_DELT_CNP   , 
     I    DEP_SNOW       ,
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG      
     1     )  
!
 6000   CONTINUE  
!
      RETURN
      END SUBROUTINE SIB0109_RUN_MONITOR_INITIAL
! ===========================================================
      SUBROUTINE SIB0109_RUN_STEPINI ( FSECP )
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     I    RCN_DELT_ATM      ,    
     I    RCN_DELT_CNP      ,  
     I    RCN_DELT_SOIL     ,  
     I    RCN_DELT_SNOW     , 
     I    JCN_SIB_UPPER_MODEL   ,
     1    JCN_IWL_SKIP
!
      USE COM_STEP_SIB0109 , ONLY : 
     I   COM_STEP_SIB0109_RUN ,
     I   SCN_SIB0109_DELT_ATM , 
     I   SCN_SIB0109_DELT
!
      use com_debug
!
      IMPLICIT NONE
!
!     INTEGER,INTENT(IN) :: ID_NOW(5)
!     INTEGER,INTENT(IN) :: ID_PRE(5)
      REAL(8),INTENT(IN) :: FSECP
!

!
! ステップ毎の定数設定
!
      CALL  COM_STEP_SIB0109_RUN (
     I    FSECP )
!
! 時間情報の取得
!
      DELT_ATM = SCN_SIB0109_DELT_ATM 
      DELT_CNP = SCN_SIB0109_DELT   
!
!  時間刻の決め方! 注意!!!!!!!!!! 
!     
      IF ( JCN_SIB_UPPER_MODEL .LT. 0 ) THEN   ! オフラインモデル
        IF ( ABS(RCN_DELT_SOIL-RCN_DELT_CNP).LT.1.D-2 ) THEN
          DELT_SOIL = SCN_SIB0109_DELT
          DELT_SNOW = SCN_SIB0109_DELT
        ELSE 
          DELT_SOIL = RCN_DELT_SOIL
          DELT_SNOW = RCN_DELT_SNOW
        ENDIF
      ELSE                                     ! オンライン 000927 現在
        IF ( ABS(JCN_IWL_SKIP).LE.2.5 ) THEN
          DELT_SOIL = SCN_SIB0109_DELT   
          DELT_SNOW = SCN_SIB0109_DELT   
        ELSE
          write(6,*) 'SIB0109_RUN_STEPINI error' 
          stop 999
        ENDIF
      ENDIF 
!
      RETURN
      END SUBROUTINE SIB0109_RUN_STEPINI 
! ===========================================================
      SUBROUTINE SIB0109_RUN_LOOPINI (  
     I    NB                 , MON                , 
     I    KT0                , 
     I    GPDEL_PHY          , GPHAF_PHY          , GPFUL_PHY       , 
     I    GT_PHY             , GQ_PHY             , 
     I    GU_PHY             , GV_PHY             , 
     I    ZMEAN_PHY          , ZTEMP_PHY          ,
     I    PPLI_PHY           , PPCI_PHY           , 
     I    DAYTIME_1HR_PHY    , 
     I    TMTX_PHY           , QMTX_PHY           , UVMTX_PHY        ) 
!
!  LOOP のはじめにすべきことを行なう。
!    ・IMSKLT, NMAX 等の設定。
!    ・大気情報を 物理格子から SiB 格子に
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1    JCN_SIB_ISLSCP    , 
     I    JCN_SIB_UPPER_MODEL 
!
      USE COM_STEP_SIB0109 , ONLY :
     1    ICN_SIB0109_AFTER_CALC_SOIL_SNW ,
     1    ICN_SIB0109_FIRST_STEP 
!
      USE SIBVEG , ONLY :  
     1   XLAI_mo ,
     1   Z1_mo , Z2_mo 
!
      USE NMAX_INFO , ONLY : 
     1   NMAX_INFO_RUN_SNOW 
!
      USE SIB_MONIT, ONLY : SIB_MONIT_RUN_SET_SNOWINFO
!
      use com_debug
!
      IMPLICIT NONE
!
! 入力
!
      INTEGER,INTENT(IN)     :: NB               ! バンドラベル
      INTEGER,INTENT(IN)     :: MON              ! つき情報
      INTEGER,INTENT(IN)     :: KT0 
!
!   大気情報 ( 物理格子 )
!
      REAL(8),INTENT(IN) ::  ZMEAN_PHY (IJPHY)
      REAL(8),INTENT(IN) ::  ZTEMP_PHY (IJPHY)
      REAL(8),INTENT(IN) ::  PPLI_PHY  (IJPHY)
      REAL(8),INTENT(IN) ::  PPCI_PHY  (IJPHY)
!
      REAL(8),INTENT(IN) ::  GT_PHY    (IJPHY)
      REAL(8),INTENT(IN) ::  GQ_PHY    (IJPHY)
      REAL(8),INTENT(IN) ::  GU_PHY    (IJPHY)
      REAL(8),INTENT(IN) ::  GV_PHY    (IJPHY)
      REAL(8),INTENT(IN) ::  GPFUL_PHY (IJPHY)
      REAL(8),INTENT(IN) ::  GPHAF_PHY (IJPHY)
      REAL(8),INTENT(IN) ::  GPDEL_PHY (IJPHY)
!
      REAL(8),INTENT(IN) ::  DAYTIME_1HR_PHY (IJPHY)
!
      REAL(8),INTENT(IN) ::  TMTX_PHY   (IJPHY,KMAX,2)
      REAL(8),INTENT(IN) ::  QMTX_PHY   (IJPHY,KMAX,2)
      REAL(8),INTENT(IN) ::  UVMTX_PHY  (IJPHY,KMAX,3) 
!
      INTEGER            :: I
      INTEGER            :: JL
!     INTEGER            :: J
!
      REAL(8)            :: XZ 
      REAL(8)            :: SDEP
      REAL(8)            :: FAC
      REAL(8),SAVE       :: ZERO = 0.D0
!
      JL = NB 
!
!  -------------------------
!  > IMSKLT, NMAX 等の設定 <
!  -------------------------
!    本当は、雪の計算直後に行なって、その結果を SAVE しておき、ここでは
!    それをコピーするだけでよいのだが、現状、毎ステップ計算している。
!
         CALL NMAX_INFO_RUN_SNOW (
     I   IMASK(1,NB), INFO_SNOW_ALL(1,NB) , NB         , 
     O   NMAX       , NTYP       ,     
     O   NMAX_0     , NMAX_01    , NMAX_1     , 
     O   NMAX_2     , NMAX_3     , NMAX_4     ,     
     O   IMSKLT     , IMSKLT_RX  , IMSKLT_RY  , ILIST )
!
      ipoint_new = 0
      ipmask_new = 0
!
#ifdef SIB_DEBUG
      do i=1,nmax
        if ( ipoint_ij .eq. imsklt(i) + IJPHY*(NB-1) ) then 
          ipoint_new = i
          ipmask_new = imask(ipoint_ij- IJPHY*(NB-1),NB) 
        endif
      enddo
      write(6,*) 'sib0109_ini_loop ipoint_new infomation ' 
     1           , ipoint_new , ipmask_new 
#endif
!
!  ----------------------------------------------------
!  > 大気状況 物理格子バンドデータ -> バンド SiB 格子 <
!  ----------------------------------------------------
!
      CALL SET_SIBGRD_ATMOS (
     I    NMAX               , 
     I    IMSKLT             ,
     I    GPDEL_PHY          , GPHAF_PHY          , GPFUL_PHY       , 
     I    GT_PHY             , GQ_PHY             , 
     I    GU_PHY             , GV_PHY             , 
     I    ZMEAN_PHY          , ZTEMP_PHY          ,
     I    PPLI_PHY           , PPCI_PHY           , 
     I    TMTX_PHY           , QMTX_PHY           , UVMTX_PHY       , 
     I    DAYTIME_1HR_PHY    , 
     O    ZMEAN              , ZTEMP              , 
     O    RAIN_CONV          , RAIN_LARGE         ,
     O    TMP_REF            , HUM_REF            , 
     O    U_REF              , V_REF              ,
     O    P_DEL              , P_REF              , 
     O    Z_REF              , FAC_CORTMP         ,  
     O    GMT                , GMQ                , GMUV            , 
     O    DAYTIME_1HR    )
!
!  GTM の変換
!    mj98_010901 バージョンに接続する場合。
!    大気側の拡散は s = CpT+gz で行なっている。
!    この ds/dt の式として扱われてきた項を、dT/dt の式の項に変換する。
!
      IF ( JCN_SIB_UPPER_MODEL .GT. 0 ) THEN   ! オンライン 
        CALL YM_TMTX_CONVERT_1 (
     I    HUM_REF   , 
     M    GMT         ) 
      ENDIF
!
!  -----------------------------------------------
!  > SiB 全面格子バンドデータ -> バンド SiB 格子 <
!  -----------------------------------------------
!
      CALL SET_SIBGRD_SIB (
     I   JL      , ILIST    , 
     I   NMAX    , NMAX_01  , NMAX_1   , NMAX_2  , NMAX_3   , NMAX_4  ,     
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL, 
     I   TMP_SOIL4_ALL      , 
     I   ALB_SNOW_SNW_ALL   , RAIN_1HOUR_ALL     , 
     O   TMP_CNP_NOS_PROG   , FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     O   TMP_CNP_SNW_PROG   , FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG, 
     O   TMP_GRSK_PROG      , FR_WTR_GRS_PROG    , FR_GLA_GRS_PROG    , 
     O   TMP_SNSK_PROG      , 
     O   INFO_SNOW          , 
     O   TMP_SOIL_PROG      , FR_WTR_SOIL_PROG   , FR_GLA_SOIL_PROG   , 
     O   TMP_SNSL_PROG      , FR_WTR_SNSL_PROG   , FR_GLA_SNSL_PROG   ,  
     O   TMP_SNOW_PROG      , WTR_SNOW_PROG      , GLA_SNOW_PROG      , 
     O   RHO_SNOW_INV_PROG  , AGE_SNOW_PROG      , 
     O   ENG_SNOW_BUCKET    , H2O_SNOW_BUCKET    , 
     O   TMP_SOIL4_PROG     , 
     O   ALB_SNOW_SNW_PROG  , RAIN_1HOUR_PROG ) 
!
         RAIN_1HOUR_PROG(1:NMAX) = RAIN_1HOUR_PROG(1:NMAX) 
     1           + RAIN_CONV(1:NMAX)  + RAIN_LARGE(1:NMAX)
!
! キャノピー評価用土壌・積雪第一層温度
!   注意。バケツ同様、土壌・雪の計算が終ると clear 
!
      IF ( ICN_SIB0109_AFTER_CALC_SOIL_SNW .EQ. 1 ) THEN  ! 土壌・雪計算直後
        CALL MAKE_TMP_SS1 ( 
     I   NMAX_01          , NMAX     ,  
     I   TMP_SOIL_PROG    , TMP_SNOW_PROG , 
     O   TMP_SOIL1_TMP_ALL(1,NB)     , TMP_SNOW1_TMP_ALL(1,NB) )
      ENDIF
!
#ifdef SIB_DEBUG
      DO I=1,NMAX_01(3)
        IF ( ipoint_new.eq.i ) then
          write(6,*) 'mj_loopini tmpsoil1 ' , 
     1               tmp_soil1_tmp_all(i,NB) , tmp_grsk_prog(i)
        endif
      enddo
      DO I=NMAX_01(2),nmax
        IF ( ipoint_new.eq.i ) then
          write(6,*) 'mj_loopini tmpsnow1 ' , 
     1               tmp_snow1_tmp_all(i,NB) , tmp_snsk_prog(i)
        endif
      ENDDO
#endif
!
! -----------------------------------------------
! > SiB 格子データ作成 （SiB 格子 -> SiB 格子） <
! -----------------------------------------------
!
      CALL MAKE_SIBGRD_SIB (  
     I  NMAX     , 
     I  NMAX_0   , NMAX_01 , NMAX_1  , NMAX_2   , NMAX_3  , NMAX_4   , 
     I  NTYP     , MON     ,
     I  WTR_SNOW_PROG       , GLA_SNOW_PROG   , 
     I  RHO_SNOW_INV_PROG   ,  
     O  CVR_VEG_ms      , CVR_VEG_INV_ms ,   
     O  CVR_SNW         , CVR_SNW_INV    , 
     O  CVR_NOS         , CVR_NOS_INV    , 
     O  DEP_SNOW      )  
!
#ifdef MONYOS
      CALL SIB_MONIT_RUN_SET_SNOWINFO (
     I   NMAX    , 
     I   CVR_SNW , CVR_NOS , IMSKLT , NMAX_01 , NB )
#endif
!
! LAI の変更 ( 雪あり域 )
!
      DO I = 1,NMAX 
        XLAI_VEG_NOS (I,1)    = XLAI_mo (NTYP(I),MON,1) 
        XLAI_VEG_NOS (I,2)    = XLAI_mo (NTYP(I),MON,2) 
        XLAI_VEG_SNW (I,1)    = XLAI_mo (NTYP(I),MON,1) 
        XLAI_VEG_SNW (I,2)    = XLAI_mo (NTYP(I),MON,2) 
      ENDDO
!
      DO I=NMAX_01(2),NMAX 
        XZ   = Z2_mo  ( NTYP(I) , MON ) - Z1_mo  ( NTYP(I) , MON )  
        SDEP = MIN ( DEP_SNOW(I,5) - Z1_mo  ( NTYP(I) , MON )  
     1                  , XZ * 0.95D0 )
        SDEP   = MAX ( SDEP , ZERO ) 
        FAC    = 1.0D0 - SDEP / XZ            
! 
        XLAI_VEG_SNW (I,1)    = XLAI_mo (NTYP(I),MON,1) * FAC
        XLAI_VEG_SNW (I,2)    = 0. 
      ENDDO
!
! -----------------------------------------------------
! > TEST1DIM 、初期値の場合、予報変数変更している措置 < 
! -----------------------------------------------------
!  ここにあるのは正しいか???
!
      IF (       ICN_SIB0109_FIRST_STEP .EQ. 1 
     1     .AND. KT0.EQ.0 
     1     .AND. JCN_SIB_ISLSCP.EQ.0 ) THEN
        CALL MOD_INI_DATA (
     I      NMAX_01 , NMAX , MON  ,  NTYP ,           
     M      FR_WTR_CNP_NOS_PROG , FR_WTR_GRS_PROG, 
     M      FR_WTR_CNP_SNW_PROG , 
     M      FR_GLA_CNP_NOS_PROG , FR_GLA_GRS_PROG, 
     M      FR_GLA_CNP_SNW_PROG  ) 
      ENDIF
!
! ------------------
! > 収支チェック用 <
! ------------------
!
      CALL TOTAL_H2O_ENG_OF_SKIN (
     I    NMAX    , MON       , NTYP   , 
     I    CVR_SNW        , 
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG     , 
     I    DEP_SNOW            , 
     O    H2O_CNP_NOS_1       , H2O_CNP_SNW_1       ,
     O    H2O_GRSK_1      , 
     O    ENG_CNP_NOS_1       , ENG_CNP_SNW_1       ,
     O    ENG_GRSK_1          , ENG_SNSK_1          )  
!
      RETURN
      END SUBROUTINE SIB0109_RUN_LOOPINI
! ===========================================================
      SUBROUTINE YM_TMTX_CONVERT_1 (
     I    HUM_REF   , 
     M    GMT         ) 
!
!     mj98 では、s = CpT+gz として、ds/dt の式を書く
!     その結果与えられた
!          GMT(1) ds_1/dt = GMT(2)         : GMT(1) は無次元
!     の式を
!          GMT(1) dT_1/dt = GMT_mod(2)     : GMT(1) は無次元
!     の形にする必要がある。明らかに
!        GMT_mod(2) = GMT(2) / Cp 
!     である。
!   
      USE SIBPRM , ONLY : IDX 
      USE SIBCON , ONLY : CP_AIR_c
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN)    :: HUM_REF(IDX)
      REAL(8),INTENT(INOUT) :: GMT    (IDX,2)
!
! WORK
      INTEGER               :: I
!
!     注 : どうでもいいと思っているので(^^;) 0.846 という数をモロにだす。
!          真面目にやるなら他の物理定数から算出すべし。
!
      DO I=1,NMAX
!       GMT(I,2) = GMT(I,2) / CP_AIR_c / ( 1 + 0.846 * HUM_REF(I) ) 
        GMT(I,2) = GMT(I,2) 
      ENDDO
!
      RETURN
      END SUBROUTINE YM_TMTX_CONVERT_1 
! ===========================================================
      SUBROUTINE SIB0109_RUN_ALBEDO ( 
     I   MON   , JL    , 
     O   AVISB , ANIRB , AVISD , ANIRD )
!
!  陸のアルベドの計算
!  短波放射計算ステップのみ計算
!  現在は出力 AVISB などは使われていない
!
      USE COM_STEP_SIB0109 , ONLY :
     I   ICN_SIB0109_SW
!
      USE SIBCON, ONLY : TMP_FREZ_C ,
     1                   H2O_SNOW_01_1_INV_c
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)  :: MON
      INTEGER,INTENT(IN)  :: JL 
!
      REAL(8),INTENT(OUT) :: AVISB(IJPHY)
      REAL(8),INTENT(OUT) :: AVISD(IJPHY)
      REAL(8),INTENT(OUT) :: ANIRB(IJPHY)
      REAL(8),INTENT(OUT) :: ANIRD(IJPHY)
!
      INTEGER :: I
      INTEGER :: J
      REAL(8) :: ALB_SIB_SIBMEAN (IDX,IVN,IBD)      ! SiB 格子平均 SALB 
      REAL(8) :: CVR_SNW_FOR_ALBEDO (IDX)   
!
      INTEGER :: IJ
      INTEGER :: IS
!
      IF ( ICN_SIB0109_SW .EQ. 0 ) RETURN
!
!   部分積雪なしの場合 ( 少雪バケツに入っている場合 ) : 
!     雪なし域において、実効的な雪の被覆率を求めておく
!     とりあえず水当量で 4mm のときを 1 にする。
!
      CVR_SNW_FOR_ALBEDO (:) = 0.  
      DO I=1,NMAX_0(3)
        IF ( TMP_GRSK_PROG(I).LE.TMP_FREZ_c ) THEN
          CVR_SNW_FOR_ALBEDO (I) = H2O_SNOW_BUCKET(I) 
     1                             * H2O_SNOW_01_1_INV_c
        ENDIF 
      ENDDO
!
! ALBEDO 等の計算 ( SiB 格子 )
!
      CALL ALBEDO_MAIN  (
     I   NMAX        , NMAX_01    , 
     I   NTYP        , MON        , 
     I   CVR_SNW_FOR_ALBEDO , 
     I   FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I   TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , TMP_SNSK_PROG     , 
     I   DEP_SNOW    , AGE_SNOW_PROG   , ZMEAN          ,        
     I   CVR_SNW     , CVR_VEG_ms , 
     I   XLAI_VEG_NOS,  XLAI_VEG_SNW    ,  
     M   ALB_SNOW_SNW_PROG , RAIN_1HOUR_PROG ,    
     O   EXTK_NOS    , ABSORB_SHRT_NOS  , 
     O   REACH_GRD_NOS            , 
     O   EXTK_SNW    , ABSORB_SHRT_SNW  , 
     O   REACH_GRD_SNW            , 
     O   ALB_SIB_SIBMEAN ) 
!
      DO I=1,NMAX
        IJ = ILIST(I)
        ALB_SNOW_SNW_ALL (IJ,JL,1)= ALB_SNOW_SNW_PROG(I,1)
        ALB_SNOW_SNW_ALL (IJ,JL,2)= ALB_SNOW_SNW_PROG(I,2)
        RAIN_1HOUR_ALL   (IJ,JL)  = RAIN_1HOUR_PROG  (I)
      ENDDO
!
! 大気への出力の計算  ( ISPT 対応 )
!
      IF ( ISPT.EQ.1 ) THEN
        DO I=1,NMAX        
          IJ = ILIST(I)
          AVISB (IJ)= ALB_SIB_SIBMEAN (I,1,1) 
          AVISD (IJ)= ALB_SIB_SIBMEAN (I,1,2) 
          ANIRB (IJ)= ALB_SIB_SIBMEAN (I,2,1) 
          ANIRD (IJ)= ALB_SIB_SIBMEAN (I,2,2) 
        ENDDO
      ELSE   
        AVISB (:)= 0.
        AVISD (:)= 0.
        ANIRB (:)= 0.
        ANIRD (:)= 0. 
        DO IS=1,ISPT 
        DO I=1,NMAX        
          IJ = (ILIST(I)-1) / ISPT+1 
          AVISB (IJ)= AVISB(IJ) 
     1              + ALB_SIB_SIBMEAN (I,1,1) * RATIO_VEG(ILIST(I),JL) 
          AVISD (IJ)= AVISD(IJ) 
     1              + ALB_SIB_SIBMEAN (I,1,2) * RATIO_VEG(ILIST(I),JL) 
          ANIRB (IJ)= ANIRB(IJ) 
     1              + ALB_SIB_SIBMEAN (I,2,1) * RATIO_VEG(ILIST(I),JL) 
          ANIRD (IJ)= ANIRD(IJ) 
     1              + ALB_SIB_SIBMEAN (I,2,2) * RATIO_VEG(ILIST(I),JL) 
        ENDDO
        ENDDO
      ENDIF
!
!     REAL(8) :: ALB_SIB_SIBMEAN (IDX,IVN,IBD)      ! SiB 格子平均 SALB 
!
      RETURN
      END SUBROUTINE SIB0109_RUN_ALBEDO
!=============================================================
      SUBROUTINE SIB0109_RUN_SR_CALC (  
     I   MON         , NB             ,
     I   RVISB       , RVISD          , RNIRB         , RNIRD       ) 
! 
!    短波計算ステップだけ計算
!
!    短波放射計算結果をもとに、
!      気孔抵抗 PAR 依存性などを計算、
!      必要な全面情報を module 変数に格納
!
      USE COM_STEP_SIB0109 , ONLY :
     I   ICN_SIB0109_SW
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: MON 
      INTEGER,INTENT(IN)  :: NB 
      REAL(8),INTENT(IN)  :: RNIRB ( IJPHY/INTRI/INTRJ ) ! 放射、直達可視
      REAL(8),INTENT(IN)  :: RNIRD ( IJPHY/INTRI/INTRJ ) ! 放射、散乱可視
      REAL(8),INTENT(IN)  :: RVISB ( IJPHY/INTRI/INTRJ ) ! 放射、直達近赤
      REAL(8),INTENT(IN)  :: RVISD ( IJPHY/INTRI/INTRJ ) ! 放射、散乱近赤
!
! WORK        
!
      REAL(8)             :: RAD_SHRT_INTO_SIB (IDX,IVN,IBD)    
!
#ifdef SIB_DEBUG
      write(6,*) 'sib0109_RUN_SHORTRAD_CALC try ' , ICN_SIB0109_SW
#endif
!
      IF ( ICN_SIB0109_SW .EQ. 0 ) RETURN
!
#ifdef SIB_DEBUG
      write(6,*) 'sib0109_RUN_SHORTRAD_CALC go ' 
#endif
!
! SiB 格子に詰める
! 
      CALL RAD_SHRT_RGRD_2_SIBGRD (
     I    IMSKLT_RX , IMSKLT_RY , NMAX    ,
     I    RVISB     , RVISD     , RNIRB   , RNIRD   , 
     O    RAD_SHRT_INTO_SIB  )
!
! 短波関係計算
!
      CALL RAD2SIB_CALC_SHRT_MAIN ( 
     I   NMAX        , NMAX_01        , 
     I   NTYP        , MON            , 
     I   RAD_SHRT_INTO_SIB , 
     I   CVR_VEG_ms  , 
     I   EXTK_NOS    , ABSORB_SHRT_NOS, REACH_GRD_NOS , 
     I   EXTK_SNW    , ABSORB_SHRT_SNW, REACH_GRD_SNW , 
     I   ZMEAN       , DAYTIME_1HR    , 
     I   XLAI_VEG_NOS  , XLAI_VEG_SNW  ,  
     O   R_STM_BAR_PAR_NOS  , R_STM_BAR_PAR_SNW      , 
     O   RAD_SHRT_NOS, RAD_SHRT_SNW )
!
! 全面 SiB 格子にセーブ
!
      CALL RAD_SAVE ( 
     I    NMAX               , NMAX_01            ,
     I    NB                 , ILIST              , 
     I    RAD_SHRT_NOS       , RAD_SHRT_SNW       , 
     I    R_STM_BAR_PAR_NOS  , R_STM_BAR_PAR_SNW  , 
     O    RAD_SHRT_NOS_SAVE  , RAD_SHRT_SNW_SAVE  , 
     O    R_STM_BAR_PAR_NOS_SAVE , R_STM_BAR_PAR_SNW_SAVE    )
!
      RETURN
      END SUBROUTINE SIB0109_RUN_SR_CALC
!=====================================================================
      SUBROUTINE SIB0109_RUN_SR_RESTORE ( NB )  
! 
!  短波計算無しのステップだけ計算
!  module にセーブされた必要な全面情報を SiB 格子に取り出す
!
      USE COM_STEP_SIB0109 , ONLY :
     I   ICN_SIB0109_SW
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NB
!
#     ifdef SIB_DEBUG
        write(6,*) 'sib0109_run_shortrad_restore try ',ICN_SIB0109_SW
#     endif
!
      IF ( ICN_SIB0109_SW .EQ. 1 ) RETURN
!
#     ifdef SIB_DEBUG
        write(6,*) 'sib0109_run_shortrad_restore go ' 
#     endif
!
      CALL RAD_RESTORE ( 
     I  NMAX               , NMAX_01            ,
     I  NB                 , ILIST              , 
     I  RAD_SHRT_NOS_SAVE        , RAD_SHRT_SNW_SAVE         , 
     I  R_STM_BAR_PAR_NOS_SAVE   , R_STM_BAR_PAR_SNW_SAVE    , 
     O  RAD_SHRT_NOS             , RAD_SHRT_SNW              , 
     O  R_STM_BAR_PAR_NOS        , R_STM_BAR_PAR_SNW    )
!
      RETURN
      END SUBROUTINE SIB0109_RUN_SR_RESTORE
!=====================================================================
      SUBROUTINE SIB0109_RUN_LR_SET ( 
     I   DLWB  )
! 
!  全ステップ計算
!    現状では長波放射ありのステップだけ計算、にはできない。
!    そうすることは可能だが、
!    雪計算で格子順の並べ替えがあったときの対応をしておく必要あり。
!
!     USE COM_STEP_SIB0109 , ONLY :
!    I   ICN_SIB0109_LW
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN)  :: DLWB  ( IJPHY/INTRI/INTRJ )  
!
!     IF ( ICN_SIB0109_LW .EQ. 0 ) RETURN 
!
      CALL RAD_LONG_RGRD_2_SIBGRD (
     I   IMSKLT_RX , IMSKLT_RY , NMAX      ,
     I   DLWB      , 
     O   RAD_LONG_INTO_SIB   )
!
      RETURN
      END SUBROUTINE SIB0109_RUN_LR_SET
!====================================================================
      SUBROUTINE SIB0109_RUN_SIBMAIN (
     I   MON     , NB   , ID_NEXT     , ID_NOW ) 
!
      USE COM_STEP_SIB0109    , ONLY : 
     1   ICN_SIB0109_CALC_SOIL_SNOW ,
     1   ICN_SIB0109_AFTER_CALC_SOIL_SNW 
!
      USE COM_STEP_SIB0109 , ONLY : 
     1    ICN_SIB0109_MONITOR
!
      use canopy_main, only : canopy_main__run
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)  :: ID_NEXT (5)
      INTEGER,INTENT(IN)  :: ID_NOW  (5)
      INTEGER,INTENT(IN)  :: MON
      INTEGER,INTENT(IN)  :: NB
!
!
! SAVE ( 土壌計算直後のみ計算、セーブ )
!
      REAL(8),SAVE  :: HDF_GRSK_SAVE     (IDX,2  ,NUM_BLOCK)  ! 土   0-1
      REAL(8),SAVE  :: HDF_SNSK_SAVE     (IDX,2  ,NUM_BLOCK)  ! 雪下 0-1
      REAL(8),SAVE  :: HDF_SOIL_H_SAVE   (IDX,IDP,NUM_BLOCK)  ! 土 1-2,2-3,3-
      REAL(8),SAVE  :: HDF_SNSL_H_SAVE   (IDX,IDP,NUM_BLOCK)  ! 雪下土   1-2
!
      REAL(8),SAVE  :: HDF_SNOW_H_SAVE   (IDX,ISN+1,NUM_BLOCK)  
!                                  ! 雪中 1-2,2-3,3-4 間の熱伝導係数
!                                  ! ただし雪層数 i として
!                                  ! i には土壌との熱伝導係数
!
      REAL(8),SAVE  :: CAP_SOIL_TOP_SAVE       (IDX,NUM_BLOCK) 
      REAL(8),SAVE  :: CAP_SNOW_TOP_SAVE       (IDX,NUM_BLOCK) 
      REAL(8),SAVE  :: ENG_SOIL_TOP_SAVE       (IDX,NUM_BLOCK) 
!                      ! 雪無しなら雪無し土壌第一層エネルギー
!                      ! 部分積雪でもあれば、雪あり土壌第一層エネルギー
!                      !   ( 雪あり部分だけでのエネルギー )
!
! WORK
!
! 雪被覆率
!
!     REAL(8)   :: CVR_SNW_NEXT     (IDX)
      REAL(8)   :: CVR_SNW_INV_NEXT (IDX)
      REAL(8)   :: CVR_NOS_NEXT     (IDX)
      REAL(8)   :: CVR_NOS_INV_NEXT (IDX)
!
!  雪 -> 土壌
!
      REAL(8)   :: WTR_BUCKET_SNOW_TO_SNSL       (IDX)
      REAL(8)   :: ENH_BUCKET_SNOW_TO_SNSL       (IDX)
      REAL(8)   :: ENG_BUCKET_SNOW_TO_SNSL       (IDX)
!
!
! とりあえず計算・出力しているだけ。
!   ランオフ ... 大気海洋結合モデルでは必要になろう 
!
      REAL(8)  :: RUNOFF (IDX)               ! ランオフ
      REAL(8)  :: US     (IDX)               ! ダミー
      REAL(8)  :: VS     (IDX)               ! ダミー
      REAL(8)  :: TS     (IDX)               ! ダミー
      REAL(8)  :: QS     (IDX)               ! ダミー
      REAL(8)  :: U_FLUX (IDX)               ! ダミー
      REAL(8)  :: V_FLUX (IDX)               ! ダミー
      REAL(8)  :: Z0     (IDX)               ! ダミー
      REAL(8)  :: FLSH   (IDX)               ! ダミー
      REAL(8)  :: FLLH   (IDX)               ! ダミー
      REAL(8)  :: QFLUX  (IDX)               ! ダミー
!
      REAL(8)  :: BERG   (IDX)               ! ダミー ( 氷河流出 ) 

!
      INTEGER  :: I
!
! ----------------------
! > 放射情報の補正など <
! ----------------------
!   RAD_SHRT_* はここで毎ステップ補正される。
!   特に RAD_SHRT_SNW については、雪スキンを抜けて雪層にはいる分は、
!   rad_shrt_snw(I,2) から除かれ、rad_shrt_into_snow_layer に入る。
!
      CALL RAD2SIB_SUM_MAIN ( 
     I   NMAX        , NMAX_01        , 
#ifdef MONYOS
     I   DELT_CNP    ,
#endif
     I   NTYP        , MON            , 
     I   CVR_VEG_ms  , CVR_VEG_INV_ms , 
     I   XLAI_VEG_NOS , XLAI_VEG_SNW ,
     I   ZMEAN       , ZTEMP          , 
     I   TMP_GRSK_PROG    , TMP_SNSK_PROG       , 
     I   TMP_CNP_NOS_PROG , TMP_CNP_SNW_PROG    , 
     I   RAD_LONG_INTO_SIB  , 
     M   RAD_SHRT_NOS , RAD_SHRT_SNW  ,          
     O   RAD_SHRT_INTO_SNOW_LAYER     , 
     O   THERMK_NOS    ,  
     O   THERMK_SNW    , 
     O   RAD_LONG_NOS_2_ATM , RAD_LONG_SNW_2_ATM , 
     O   RAD_LONG_NOS , RAD_LONG_SNW )   
!
!  ----------------
!  > 陸面計算準備 <
!  ----------------
!
      IF ( ICN_SIB0109_AFTER_CALC_SOIL_SNW .EQ. 1 ) THEN
        CALL PRESIB_MAIN  ( 
     I    NMAX_01     , NMAX_1      , NMAX_2      , NMAX_3      ,
     I    NMAX_4      , NMAX        ,
#ifdef MONYOS
     I    DELT_SOIL   ,
#endif
     I    NTYP        , 
     I    CVR_SNW_INV , 
     I    DEP_SNOW    , RHO_SNOW_INV_PROG   ,
     I                       WTR_SNOW_PROG    , GLA_SNOW_PROG    , 
     I    TMP_SOIL_PROG    , FR_WTR_SOIL_PROG , FR_GLA_SOIL_PROG , 
     I    TMP_SNSL_PROG    , FR_WTR_SNSL_PROG , FR_GLA_SNSL_PROG ,
     O    HDF_GRSK_SAVE  (1,1,NB)  , HDF_SNSK_SAVE  (1,1,NB)  ,
     O    HDF_SOIL_H_SAVE(1,1,NB)  , HDF_SNSL_H_SAVE(1,1,NB)  ,
     O    HDF_SNOW_H_SAVE(1,1,NB)  , 
     O    CAP_SOIL_TOP_SAVE(1,NB)  , CAP_SNOW_TOP_SAVE(1,NB)  ,
     O    ENG_SOIL_TOP_SAVE(1,NB)      )
      ENDIF
!
!  ----------------------
!  > キャノピー周辺計算 <
!  ----------------------
!
#ifdef SIB_DEBUG
      write(6,*) 'sib0109_run_sibmain' , 
     1    delt_cnp , delt_atm , delt_soil , delt_snow 
#endif
!
      CALL CANOPY_MAIN__run (
     I    NMAX_0        , NMAX_01       , NMAX_1     , NMAX        , 
     I    NTYP          , MON           , IMSKLT     ,
     I    CVR_SNW       , CVR_NOS       , 
     I    HDF_GRSK_SAVE(1,1,NB)         , HDF_SNSK_SAVE(1,1,NB)    , 
     I    DELT_ATM      , DELT_CNP      , 
     I    P_REF         , U_REF         , V_REF     , Z_REF     , 
     I    HUM_REF       , P_DEL         , 
     I    FAC_CORTMP    , 
     I    GLA_SNOW_PROG , WTR_SNOW_PROG , 
     I    TMP_SNOW_PROG , TMP_SOIL_PROG , TMP_SNSL_PROG ,  TMP_REF   , 
     I    FR_WTR_SOIL_PROG , FR_GLA_SOIL_PROG   , 
     I    FR_WTR_SNSL_PROG , FR_GLA_SNSL_PROG , 
     I    CVR_VEG_ms    , CVR_VEG_INV_ms, 
     I    THERMK_NOS    , THERMK_SNW    , 
     I    R_STM_BAR_PAR_NOS , R_STM_BAR_PAR_SNW ,  
     I    DEP_SNOW      , 
     I    RAD_SHRT_NOS  , RAD_SHRT_SNW  ,  
     I    RAD_LONG_NOS  , RAD_LONG_SNW  , 
     I    RAD_SHRT_INTO_SNOW_LAYER      ,
     I    RAIN_LARGE    , RAIN_CONV     , 
     I    CAP_SOIL_TOP_SAVE(1,NB)       , CAP_SNOW_TOP_SAVE(1,NB)     , 
     M    TMP_CNP_NOS_PROG , FR_WTR_CNP_NOS_PROG , FR_GLA_CNP_NOS_PROG , 
     M    TMP_CNP_SNW_PROG , FR_WTR_CNP_SNW_PROG , FR_GLA_CNP_SNW_PROG , 
     M    TMP_GRSK_PROG    , FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG    , 
     M    TMP_SNSK_PROG    , 
     M    GMT           , GMQ          , GMUV       , 
     O    RAD_LONG_SIB_2_ATM           , 
     M    RAD_LONG_NOS_2_ATM           , RAD_LONG_SNW_2_ATM           ,
     M    ENH_BUCKET_INTO_SOIL  (1,NB) , ENH_BUCKET_INTO_SNOW  (1,NB) , 
     M    ENH_BUCKET_FROM_SOIL  (1,NB) , ENH_BUCKET_FROM_SNOW  (1,NB) , 
     M    H2O_BUCKET_INTO_SOIL  (1,NB) , H2O_BUCKET_INTO_SNOW  (1,NB) , 
     M    H2O_BUCKET_FROM_SOIL  (1,NB) , H2O_BUCKET_FROM_SNOW  (1,NB) , 
     M    TRNS_BUCKET_NOS     (1,1,NB) , TRNS_BUCKET_SNW     (1,1,NB) , 
     M    ENG_BUCKET_NOS        (1,NB) , ENG_BUCKET_SNW        (1,NB) , 
     M    SRAD_BUCKET_SNW       (1,NB) ,
     M    TMP_SOIL1_TMP_ALL(1,NB)      , TMP_SNOW1_TMP_ALL(1,NB)   ,
     O    TS , QS, US,  VS, U_FLUX, V_FLUX, Z0   , 
     O    FLSH , FLLH  ,qflux  )
!
      call sib0109_monit_fluxes ( 
     I      DELT_CNP , 
     I      flsh , fllh , qflux , u_flux , v_flux ) 
!
      CALL MONIT_SNAP_CANOPY_2 (
     I    NMAX    , MON    , NTYP   , 1 , 
     I    DELT_CNP, 
     I    CVR_SNW        , 
     I    DEP_SNOW,
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG      
     1   )  
!
! ====================
! >> 積雪・土壌計算 <<
! ====================
!
!  降水の積算量の計算 ( 雪無し土壌表面での透水計算時に使用する ) 
!
! ----------------------------------
! > 積雪・土壌計算ステップここから < 
! ----------------------------------
!
      IF ( ICN_SIB0109_CALC_SOIL_SNOW  .EQ. 1 ) THEN 
!
!       雪無し域で、雪や霜を少雪バケツに格納する。
!
        CAll MOVE_ICE_INTO_SNOW (
     I      1   , NMAX_01(3) , 
     M      H2O_BUCKET_INTO_SOIL(1,NB) ,
     M      ENH_BUCKET_INTO_SOIL(1,NB) ,
     M      H2O_snow_bucket  , ENG_snow_bucket , 
     M      INFO_SNOW  )
!
! ------------
! > 雪の計算 <
! ------------
!
        CALL SNOW_MAIN (
     I    NMAX     , NMAX_0     , NMAX_01    , NMAX_1     , 
     I    NMAX_2   , NMAX_3     , NMAX_4     , imsklt     , 
     I    CVR_SNW  , CVR_SNW_INV, DELT_SNOW  , DELT_CNP   , 
     I    DEP_SNOW , 
     I    HDF_SNOW_H_SAVE(1,1,NB) , 
     I    H2O_BUCKET_INTO_SNOW (1,NB) , H2O_BUCKET_FROM_SNOW  (1,NB) , 
     I    ENH_BUCKET_INTO_SNOW (1,NB) , ENH_BUCKET_FROM_SNOW  (1,NB) , 
     I    ENG_BUCKET_SNW(1,NB)        , SRAD_BUCKET_SNW(1,NB)        , 
     I    ENG_SNOW_BUCKET       , H2O_SNOW_BUCKET            , 
     I    TMP_REF  , U_REF      , V_REF      , 
     I    TMP_SNSK_PROG ,   
     O    ENG_SOIL_TOP_SAVE(1,NB), 
     M    WTR_SNOW_PROG , GLA_SNOW_PROG   , TMP_SNOW_PROG   , ! 積分済 
     M    AGE_SNOW_PROG         , RHO_SNOW_INV_PROG        ,  
     M    ENG_BUCKET_NOS      (1,NB)  , 
     M    H2O_BUCKET_INTO_SOIL(1,NB)  , ENH_BUCKET_INTO_SOIL  (1,NB) , 
     O    CVR_SNW_NEXT          , CVR_SNW_INV_NEXT         ,    
     O    CVR_NOS_NEXT          , CVR_NOS_INV_NEXT         , 
     O    INFO_SNOW_NEXT        , 
     O    DEP_SNOW_NEXT         , 
     O    ENG_BUCKET_SNOW_TO_SNSL      ,  
     O    WTR_BUCKET_SNOW_TO_SNSL  , ENH_BUCKET_SNOW_TO_SNSL  , 
     O    ENG_SNOW_BUCKET_NEXT     , H2O_SNOW_BUCKET_NEXT     ,
     O    BERG    )
!  
       ENDIF
!
! --------------
! > 土壌の計算 <
! --------------
!
      IF ( ICN_SIB0109_CALC_SOIL_SNOW  .EQ. 1 ) THEN 
!
        runoff(:) = 0.
!
        CALL SOIL_MAIN  (
     I    NMAX      , NMAX_0      , NMAX_01     , 
     I    NTYP      , 
     I    DELT_SOIL        , 
     I    CVR_SNW          , CVR_SNW_NEXT        , 
     I    CVR_SNW_INV      , CVR_NOS_INV         , 
     I    CVR_SNW_INV_NEXT , CVR_NOS_INV_NEXT    , 
     I    HDF_SOIL_H_SAVE(1,1,NB)     , HDF_SNSL_H_SAVE(1,1,NB)      , 
     I    ENG_BUCKET_NOS      (1,NB)  , ENG_BUCKET_SNOW_TO_SNSL , 
     I    H2O_BUCKET_INTO_SOIL(1,NB)  , ENH_BUCKET_INTO_SOIL  (1,NB) , 
     I    H2O_BUCKET_FROM_SOIL(1,NB)  , ENH_BUCKET_FROM_SOIL  (1,NB) , 
     I    WTR_BUCKET_SNOW_TO_SNSL     , ENH_BUCKET_SNOW_TO_SNSL      , 
     I    TRNS_BUCKET_NOS   (1,1,NB)  , TRNS_BUCKET_SNW     (1,1,NB) , 
     I    INFO_SNOW_NEXT        , 
     M    TMP_SOIL_PROG    , TMP_SNSL_PROG   ,      ! 予報変数 : 温度
     M    FR_WTR_SOIL_PROG , FR_WTR_SNSL_PROG,      ! 予報変数 : 水の割合
     M    FR_GLA_SOIL_PROG , FR_GLA_SNSL_PROG,      ! 予報変数 : 氷の割合
     M    TMP_SOIL4_PROG        ,
     O    RUNOFF         
     1     )  
!
! ----------
! > 後始末 <
! ----------
!    CVR_SNW や DEP_SNOW の変更に伴う
!      アルベドの変更 （放射前のステップを除く）
!      canopy 関連量のエネルギー・水・氷の変更
!
!      CVR_SNW の変化について、(雪と)土壌での変更は、(snow_main 及び
!      soil_main で既になされていることに注意。
!
        CALL SNOWMOD_MAIN (
     I     NMAX            , 
     I     NTYP            , MON                , 
     I     CVR_SNW         , CVR_SNW_NEXT       , 
     I     INFO_SNOW       , INFO_SNOW_NEXT     , 
     I     DEP_SNOW        , DEP_SNOW_NEXT      , 
     I     TMP_SOIL_PROG   , TMP_SNOW_PROG      , 
     M     TMP_CNP_NOS_PROG, FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     M     TMP_CNP_SNW_PROG, FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG,      
     M     TMP_GRSK_PROG   , FR_WTR_GRS_PROG    , FR_GLA_GRS_PROG    , 
     M     TMP_SNSK_PROG    )
!
        CALL MONIT_SNAP_SNOW (
     I    NMAX             , 
     I    1                , ICN_SIB0109_MONITOR  , 
     I    DELT_SNOW        , DELT_CNP , 
     I    INFO_SNOW_NEXT   , 
     I    CVR_SNW_NEXT     , 
     I    WTR_SNOW_PROG    , GLA_SNOW_PROG     , 
     I    TMP_SNOW_PROG    , DEP_SNOW_NEXT       , 
     I    H2O_SNOW_BUCKET_NEXT , 
     I    ENG_SNOW_BUCKET_NEXT , 
     I    RHO_SNOW_INV_PROG    
     1   )  
!
        CALL MONIT_SNAP_SOIL  (
     I    NMAX    ,  NTYP   , 
     I    1       ,  ICN_SIB0109_MONITOR  , 
     I    DELT_SOIL           , DELT_CNP            , 
     I    CVR_SNW_NEXT        , 
     I    TMP_SOIL_PROG       , TMP_SNSL_PROG       , TMP_SOIL4_PROG , 
     I    FR_WTR_SOIL_PROG    , FR_WTR_SNSL_PROG    , 
     I    FR_GLA_SOIL_PROG    , FR_GLA_SNSL_PROG     
     1   )  
!
!#ifdef MONYOS
!        CALL MONIT_CHECK_SS1 ( 
!     I    1           , NMAX_01(3)   , 
!     I    NMAX_01(2)  , NMAX         , 
!     I    TMP_SOIL1_TMP_ALL(1,NB)    , TMP_SOIL  ,
!     I    TMP_SNOW1_TMP_ALL(1,NB)    , TMP_SNOW_PROG,INFO_SNOW_NEXT )
!#endif
!
! ----------------------------------
! > 土壌や雪の計算ステップここまで <
! ----------------------------------
!
! --------------------------------------
! > 土壌や雪の計算無しステップここから <
! --------------------------------------
!
      ELSE             
!
! 雪の高さ DEP_SNOW_NEXT 、雪の被覆率 CVR_SNW_NEXT をコピー
!     （monit_snap_canopy 用）
!
        CALL COPY_SNOW_INFO (
     I      NMAX          , NMAX_01       , 
#ifdef MONYOS
     I   NMAX_0           , NMAX_1        , 
#endif 
     I      CVR_SNW       , DEP_SNOW      ,
     O      CVR_SNW_NEXT  , DEP_SNOW_NEXT  )
!
      ENDIF
!
! --------------------------------------
! > 土壌や雪の計算無しステップここまで <
! --------------------------------------
!
! ==============
! >> 月末処理 <<
! ==============
!    LAI , CVR_VEG 変更に伴う水・氷・エネルギー補正量の計算
!    月末に一度行なう。
!    LAI,CVR_VEG が変わると、葉の水・氷・エネルギーが消える/生まれる
!
!  ここで mon_last 処理を行なう!! 
!
      MON_NEXT = ID_NEXT(2) 
!
      IF ( ID_NOW(2) .NE. ID_NEXT(2) ) THEN
        CALL  MONIT_MOD_MON_LAST (
     I   NMAX             , NTYP             , 
     I   MON              , MON_NEXT         , 
     I   CVR_SNW_NEXT     , 
     I   DEP_SNOW_NEXT    , 
     I   TMP_CNP_NOS_PROG , FR_WTR_CNP_NOS_PROG , FR_GLA_CNP_NOS_PROG ,        
     I   TMP_CNP_SNW_PROG , FR_WTR_CNP_SNW_PROG , FR_GLA_CNP_SNW_PROG ,
     I   TMP_GRSK_PROG    , FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG     ) 
      ENDIF
!  
! ==================
! >> 瞬間値モニタ <<
! ==================
!   注意。月末のステップは、次の月の LAI 等で評価する。
!
      CALL MONIT_SNAP_CANOPY (
!    I    NMAX    , MON    , NTYP   , 
     I    NMAX    , MON_NEXT    , NTYP   , 
     I    DELT_CNP, 
     I    DEP_SNOW_NEXT,
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG      
     1   )  
! 
      RETURN     
      END SUBROUTINE SIB0109_RUN_SIBMAIN
!
!=================================================================
      subroutine sib0109_monit_fluxes ( 
     I      DELT_CNP , 
     I      flsh , fllh , qflux , u_flux , v_flux ) 
!
      use sibprm , only : IDX 
      USE SIB_MONIT, ONLY :  MONIT_SIB, INA , 
     1                      IMONIT_LEVEL , 
     1                      IMONIT_FLUX_COMMON  
!     1                      IMONIT_ALL , 
!     1                      IMONIT_BUDJET       ,
!     1                      IMONIT_FLUX_MRI
!
      implicit none
!
      real(8),intent(in) :: DELT_CNP 
      real(8),intent(in) :: flsh  (idx) 
      real(8),intent(in) :: fllh  (idx) 
      real(8),intent(in) :: qflux (idx) 
      real(8),intent(in) :: u_flux(idx) 
      real(8),intent(in) :: v_flux(idx) 
      integer :: iflag_mon 
!
      character(7) :: cmark 
      real(8),parameter :: one = 1.D0 
!
      IFLAG_MON   = INA
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
        CMARK = 'FLSH'
        CALL MONIT_SIB ( CMARK, FLSH, IFLAG_MON, DELT_CNP, ONE ) 
        CMARK = 'FLLH'
        CALL MONIT_SIB ( CMARK, FLLH, IFLAG_MON, DELT_CNP, ONE ) 
        CMARK = 'EVSPS'
        CALL MONIT_SIB ( CMARK, QFLUX , IFLAG_MON, DELT_CNP, ONE ) 
        CMARK = 'UMOM'
        CALL MONIT_SIB ( CMARK, U_FLUX , IFLAG_MON, DELT_CNP, ONE ) 
        CMARK = 'VMOM'
        CALL MONIT_SIB ( CMARK, V_FLUX , IFLAG_MON, DELT_CNP, ONE ) 
      ENDIF
      return
      end subroutine sib0109_monit_fluxes  
!===================================================================      

!====================================================================
      SUBROUTINE SIB0109_RUN_LOOPEND (
     I   NB          , 
     O   TMTX2L_PHY  , QMTX2L_PHY  , UMTX2L_PHY  , VMTX2L_PHY   ,   
     O   RAD_LONG_SIB_2_ATM_PHY  ) 
!
! 陸面計算後の処理 （全面とのやりとり）
!   ・Yamada-Mellor 用の行列を返す
!   ・予報変数をバンド SiB 格子から全面 SiB 格子に格納
!
      USE COM_STEP_SIB0109 , ONLY : 
     1    ICN_SIB0109_CALC_SOIL_SNOW 
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     I    JCN_SIB_UPPER_MODEL
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NB 
!
      REAL(8),INTENT(OUT) :: TMTX2L_PHY  (IJPHY)
      REAL(8),INTENT(OUT) :: QMTX2L_PHY  (IJPHY)
      REAL(8),INTENT(OUT) :: UMTX2L_PHY  (IJPHY)
      REAL(8),INTENT(OUT) :: VMTX2L_PHY  (IJPHY)
      ! intent(INOUT) では日立でコンパイルエラーになる(バグ)ので授受属性無指定
!
      REAL(8),INTENT(OUT)   :: RAD_LONG_SIB_2_ATM_PHY (IJPHY) ! 上向き正
!
!  --------------------------------
!  > Yamada Mellor 用の行列を返す <
!  --------------------------------
!
      IF ( JCN_SIB_UPPER_MODEL .GT. 0 ) THEN   ! オンライン 
        CALL YM_TMTX_CONVERT_2 (
     I    HUM_REF   , 
     M    GMT        ) 
      ENDIF
!
      CALL SETBACK_TO_ATM (
     I   ILIST       , NMAX           , 
     I   GMT         , GMQ            , GMUV         , 
     I   RAD_LONG_SIB_2_ATM           , 
     O   TMTX2L_PHY  , QMTX2L_PHY     , UMTX2L_PHY   , VMTX2L_PHY   , 
     O   RAD_LONG_SIB_2_ATM_PHY )
!
!  -------------------------------------------
!  > キャノピー予報変数を全面 SiB 格子に格納 <
!  -------------------------------------------
!
      CALL  SETBACK_SIBGRD_CANOPY (
     I   NB              , ILIST              , NMAX               , 
     I   TMP_CNP_NOS_PROG, FR_WTR_CNP_NOS_PROG, FR_GLA_CNP_NOS_PROG, 
     I   TMP_CNP_SNW_PROG, FR_WTR_CNP_SNW_PROG, FR_GLA_CNP_SNW_PROG, 
     I   TMP_GRSK_PROG   , FR_WTR_GRS_PROG    , FR_GLA_GRS_PROG    , 
     I   TMP_SNSK_PROG   , 
     O   TMP_CNP_NOS_ALL , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL    , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL    ) 
!
      IF ( ICN_SIB0109_CALC_SOIL_SNOW .EQ. 1 ) THEN 
!
!  ---------------------------------------------------
!  > 放射セーブ値の変更 （雪面被覆率変更に伴うもの） <
!  ---------------------------------------------------
!    本当は、次のステップでアルベド計算するなら不要
!
!                                       IF ( ICNSW_NEXT .EQ. 0 ) THEN 
        CALL RADSHRT_MOD (
     I      NMAX                , ILIST             , NB         , 
     I      CVR_SNW             , CVR_SNW_NEXT      , 
     I      INFO_SNOW           , INFO_SNOW_NEXT    , 
     M      RAD_SHRT_NOS_SAVE   , R_STM_BAR_PAR_NOS_SAVE  , 
     M      RAD_SHRT_SNW_SAVE   , R_STM_BAR_PAR_SNW_SAVE  ) 
!                                       ENDIF
!
!  -----------------------------------
!  > 雪予報変数を全面 SiB 格子に格納 <
!  -----------------------------------
!
        CALL SETBACK_SIBGRD_SNOW (
     I    NB     , ILIST , NMAX    , 
     I    TMP_SNOW_PROG      , 
     I    WTR_SNOW_PROG      , GLA_SNOW_PROG    , 
     I    RHO_SNOW_INV_PROG  , 
     I    AGE_SNOW_PROG      , 
     I    INFO_SNOW_NEXT     , 
     I    ENG_SNOW_BUCKET_NEXT    , H2O_SNOW_BUCKET_NEXT    , 
     O    TMP_SNOW_ALL      , 
     O    WTR_SNOW_ALL      , GLA_SNOW_ALL    , 
     O    RHO_SNOW_INV_ALL  , 
     O    AGE_SNOW_ALL      , 
     O    INFO_SNOW_ALL     , 
     O    ENG_SNOW_BUCKET_ALL    , H2O_SNOW_BUCKET_ALL    ) 
!
!  -------------------------------------
!  > 土壌予報変数を全面 SiB 格子に格納 <
!  -------------------------------------
!
        CALL  SETBACK_SIBGRD_SOIL (
     I    NB     , ILIST , NMAX    , 
     I    TMP_SOIL_PROG     , TMP_SNSL_PROG      , 
     I    FR_WTR_SOIL_PROG  , FR_GLA_SOIL_PROG   , 
     I    FR_WTR_SNSL_PROG  , FR_GLA_SNSL_PROG   , 
     I    TMP_SOIL4_PROG    ,
     O    TMP_SOIL_ALL      , TMP_SNSL_ALL      , 
     O    FR_WTR_SOIL_ALL   , FR_GLA_SOIL_ALL   , 
     O    FR_WTR_SNSL_ALL   , FR_GLA_SNSL_ALL   ,
     O    TMP_SOIL4_ALL     )
!
! ------------------
! > バケツのクリア <
! ------------------
!
        call reset ( H2O_BUCKET_INTO_SOIL   (1,NB), idx )
        call reset ( H2O_BUCKET_INTO_SNOW   (1,NB), idx )
        call reset ( H2O_BUCKET_FROM_SOIL   (1,NB), idx )
        call reset ( H2O_BUCKET_FROM_SNOW   (1,NB), idx )
        call reset ( ENH_BUCKET_INTO_SOIL   (1,NB), idx )
        call reset ( ENH_BUCKET_INTO_SNOW   (1,NB), idx )
        call reset ( ENH_BUCKET_FROM_SOIL   (1,NB), idx )
        call reset ( ENH_BUCKET_FROM_SNOW   (1,NB), idx )
!
        call reset ( TRNS_BUCKET_NOS        (1,1,NB), idx*idp )
        call reset ( TRNS_BUCKET_SNW        (1,1,NB), idx*idp )
!
        call reset ( ENG_BUCKET_NOS         (1,NB), idx )
        call reset ( ENG_BUCKET_SNW         (1,NB), idx )
        call reset ( SRAD_BUCKET_SNW        (1,NB), idx )
!
      ENDIF
!
! ------------------
! > 収支チェック用 <
! ------------------
!
      CALL TOTAL_H2O_ENG_OF_SKIN (
     I    NMAX    , MON_NEXT  , NTYP   , 
     I    CVR_SNW_NEXT        , 
     I    TMP_CNP_NOS_PROG    , TMP_CNP_SNW_PROG    , 
     I    FR_WTR_CNP_NOS_PROG , FR_WTR_CNP_SNW_PROG , 
     I    FR_GLA_CNP_NOS_PROG , FR_GLA_CNP_SNW_PROG , 
     I    TMP_GRSK_PROG       , TMP_SNSK_PROG       , 
     I    FR_WTR_GRS_PROG     , FR_GLA_GRS_PROG     , 
     I    DEP_SNOW_NEXT       , 
     O    H2O_CNP_NOS_2       , H2O_CNP_SNW_2       ,
     O    H2O_GRSK_2          , 
     O    ENG_CNP_NOS_2       , ENG_CNP_SNW_2       ,
     O    ENG_GRSK_2          , ENG_SNSK_2          )  
!
      CALL MONIT_TOTAL_H2O_ENG_OF_SKIN (
     I    NMAX                , 
     I    H2O_CNP_NOS_1       , H2O_CNP_SNW_1       ,
     I    H2O_GRSK_1      , 
     I    ENG_CNP_NOS_1       , ENG_CNP_SNW_1       ,
     I    ENG_GRSK_1          , ENG_SNSK_1          ,
     I    H2O_CNP_NOS_2       , H2O_CNP_SNW_2       ,
     I    H2O_GRSK_2      , 
     I    ENG_CNP_NOS_2       , ENG_CNP_SNW_2       ,
     I    ENG_GRSK_2          , ENG_SNSK_2          )  
!
      RETURN
      END SUBROUTINE SIB0109_RUN_LOOPEND
!====================================================================-
      SUBROUTINE YM_TMTX_CONVERT_2 (
     I    HUM_REF   , 
     M    GMT_OUT      ) 
!
!     mj98 では、s = CpT+gz として、ds/dt の式を書く
!     陸面で求まった
!          dT_1/dt = GMT(2) 
!     の式を
!          ds_1/dt = GMT_mod(2) 
!     の形にする必要がある。明らかに
!        GMT_mod(2) = GMT(2) x Cp 
!     である。
!   
      USE SIBPRM , ONLY : IDX 
      USE SIBCON , ONLY : CP_AIR_c
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN)    :: HUM_REF(IDX)
      REAL(8),INTENT(INOUT) :: GMT_OUT(IDX)
!
! WORK
      INTEGER               :: I
!
!     注 : どうでもいいと思っているので(^^;) 0.846 という数をモロにだす。
!          真面目にやるなら他の物理定数から算出すべし。
!
      DO I=1,NMAX
        GMT_OUT(I) 
     1      = GMT_OUT(I) 
!       GMT_OUT(I) 
!    1      = GMT_OUT(I) * CP_AIR_c * ( 1 + 0.846 * HUM_REF(I) ) 
      ENDDO
!
      RETURN
      END SUBROUTINE YM_TMTX_CONVERT_2 
!=========================================================
      SUBROUTINE SIB0109_RUN_STEPEND (
!    I    ID_NEXT , ID_NOW , FSECP , KTP )
     I    ID_NEXT , ID_NOW )
!
!  各ステップ終了時にする処理
!
      USE COM_STEP_SIB0109 , ONLY :
     1   ICN_SIB0109_LAST_STEP
!
      use restart_sib0109 , only : 
     1      clean_RESTART_SIB0109  ,
     1      write_RESTART_SIB0109  
!
      USE COM_JOBINFO_SIB0109 , ONLY : 
     1    CDIR_RESTART

!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: ID_NOW (5)
      INTEGER,INTENT(IN)    :: ID_NEXT(5)
!     REAL(8),INTENT(IN)    :: FSECP
!     INTEGER,INTENT(IN)    :: KTP

!
      INTEGER,SAVE          :: IUNIT_RESTART_OUT = 11               
      character(100)        :: cfile_restart_out 
!
!  リスタートデータをきれいにする
!
      CALL CLEAN_RESTART_SIB0109 ( 
     M   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     M   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     M   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     M   TMP_SNSK_ALL       , 
     M   INFO_SNOW_ALL      , 
     M   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     M   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     M   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     M   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       ) 
!
!  リスタートファイル出力
!
      IF (      ICN_SIB0109_LAST_STEP.EQ.1 
     1     .OR. ID_NEXT(2).NE.ID_NOW(2) ) THEN  
        CFILE_RESTART_OUT  = TRIM(CDIR_RESTART) // 
     1                       'newsib_restart_YYYY_MM_DD_HH'
        CALL REPLACE_INT(CFILE_RESTART_OUT, 'YYYY', ID_NEXT(1))
        CALL REPLACE_INT(CFILE_RESTART_OUT, 'MM'  , ID_NEXT(2))
        CALL REPLACE_INT(CFILE_RESTART_OUT, 'DD'  , ID_NEXT(3))
        CALL REPLACE_INT(CFILE_RESTART_OUT, 'HH'  , ID_NEXT(4))
!
        WRITE(6,*) 'Output restart file is ' ,  
     1      CFILE_RESTART_OUT , ' ' , CFILE_RESTART_OUT
!
        open ( IUNIT_RESTART_OUT  , file=cfile_restart_out , 
     1         form='unformatted' ) 
!
         CALL WRITE_RESTART_SIB0109 ( 
     I   IUNIT_RESTART_OUT  , 
     I   ID_NEXT            , 
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL,  
     I   TMP_SOIL4_ALL      , 
     I   ALB_SNOW_SNW_ALL   ,  RAIN_1HOUR_ALL ) 
!
        close ( iunit_restart_out ) 
!
      ENDIF
!
      RETURN
      END SUBROUTINE SIB0109_RUN_STEPEND

!====================================================================
      END MODULE SIB0109

! 蒸散関連計算
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/09/28 保坂征宏 バグだし。icheck 初期化忘れ追加。
! 99/09/28 保坂征宏 終了処理を分割。canopy_set_air に
!
      SUBROUTINE CANOPY_GET_EA (
     I    NSTR         , NEND         , INSFLAG     , 
     I    ITR          , ITR_MAX      , 
     I    NTYP         , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR     , R_AC_BAR_INV , R_AG        , R_AG_INV    , 
     I    R_BB         , R_AIR_INV    , 
     I    R_ROOT_BAR_INV   , R_STM_BAR_PAR    , 
     I    E_AIR_SAT    , E_CNP_SAT     , E_REF     , E_SKIN_SAT   ,
     I    HUM_GRSK     , 
     I    TMP_SOIL     , 
     I    FR_H2O_CNP  , FR_H2O_GRS  , 
     I    FUNC_TMP     , 
     I    POT_SOIL     ,
     I    D_ZERO_msx   , 
     I    FAC_EVP      , 
     I    ITC          , ITC_MAX     , 
#ifdef MONYOS
     I   DELT_CNP , 
#endif
     M    ICOUNT       , INIT_EA     , 
     M    R_STM_BAR    , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN , 
     O    E_AIR        , FLAG_TR     , FLAG_TR_ALL     , 
     O    R_ACX_BAR    , R_AGX_BAR  , R_ACX_BAR_INV, R_AGX_BAR_INV ,
     O    R_AB_INV     , 
!    O    FAC_TRNS_A   , FAC_TRNS_B , KLEV_TRNS_MAX )
     O    FAC_TRNS_A   , FAC_TRNS_B )
!
      USE SIBPRM , ONLY :
     1   IDX , 
     1   IDP , 
     1   ICG , 
     1   ITYP 
!
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR 
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: ITR
      INTEGER,INTENT(IN)    :: ITR_MAX
      INTEGER,INTENT(IN)    :: INSFLAG 
#ifdef MONYOS
      REAL(8),INTENT(IN)   :: DELT_CNP
#endif
!
      REAL(8),INTENT(IN)   :: TMP_SOIL(IDX,IDP)
!
      REAL(8),INTENT(IN)   :: HUM_GRSK (IDX)
!
      REAL(8),INTENT(IN)   :: R_AC_BAR    (IDX) ! R_AC(キャ-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_AC_BAR_INV(IDX) ! R_AC(キャ-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_AG    (IDX)     ! R_AG(地面-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_AG_INV(IDX)     ! R_AG(地面-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_BB    (IDX)     ! 地面抵抗 RSOIL
      REAL(8),INTENT(IN)   :: R_AIR_INV (IDX)     
      REAL(8),INTENT(IN)   :: R_ROOT_BAR_INV(IDX,IDP,ICG) 
!                                               ! 各層の(透水＋吸水)抵抗
!
      REAL(8),INTENT(IN)   :: E_REF     (IDX)   ! 参照空気水蒸気圧
      REAL(8),INTENT(IN)   :: E_CNP_SAT (IDX)   ! キャノピー温度飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_AIR_SAT (IDX)   ! キャノピー空気飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_SKIN_SAT(IDX)   ! 地面温度飽和水蒸気圧
!
      REAL(8),INTENT(IN)   :: FUNC_TMP(IDX,ICG)    ! 気孔抵抗の温度依存性
!
      REAL(8),INTENT(IN)   :: POT_SOIL (IDX,IDP)    ! 土壌のmatric pot
      REAL(8),INTENT(IN)   :: D_ZERO_msx (IDX)        ! 零面変位   (旧；Ｄ)
!
      REAL(8),INTENT(IN)   :: CVR_VEG_msx (IDX,ICG)  ! 植生被覆率
      REAL(8),INTENT(IN)   :: CVR_VEG_INV_msx (IDX,ICG)  ! 植生被覆率
      REAL(8),INTENT(IN)   :: FAC_EVP    (IDX)    ! ρε/P_s

      REAL(8),INTENT(IN)   :: R_STM_BAR_PAR (IDX,ICG) ! 
!
      INTEGER,INTENT(IN)   :: NTYP   (IDX)       ! 植生インデックス

      REAL(8),INTENT(IN)   :: FR_H2O_CNP (IDX)   ! キャノピーの濡れ領域の割合
      REAL(8),INTENT(IN)   :: FR_H2O_GRS (IDX)   ! 下草の濡れ領域の割合
!
! INOUT
!
      REAL(8),INTENT(INOUT) :: R_STM_BAR        (IDX,ICG)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: FLAG_TR     (IDX,IDP,ICG)
      REAL(8),INTENT(OUT)   :: FLAG_TR_ALL (IDX    ,ICG)
!
      REAL(8),INTENT(OUT)   :: R_AGX_BAR        (IDX) 
      REAL(8),INTENT(OUT)   :: R_AGX_BAR_INV    (IDX) 
      REAL(8),INTENT(OUT)   :: R_ACX_BAR        (IDX) 
      REAL(8),INTENT(OUT)   :: R_ACX_BAR_INV    (IDX) 
      REAL(8),INTENT(OUT)   :: R_AB_INV         (IDX) 
!
      REAL(8),INTENT(OUT)  :: FAC_TRNS_A (IDX,IDP,ICG)  
      REAL(8),INTENT(OUT)  :: FAC_TRNS_B (IDX,IDP,ICG)  
!
      REAL(8),INTENT(OUT)  :: E_AIR  (IDX)         ! キャノピー水蒸気圧
!     INTEGER,INTENT(OUT)  :: KLEV_TRNS_MAX(IDX,ICG)
!
      REAL(8),INTENT(INOUT):: FLAG_DRY_CNP      (IDX)    
      REAL(8),INTENT(INOUT):: FLAG_DRY_SKIN     (IDX)    
      INTEGER,INTENT(IN)   :: ITC
      INTEGER,INTENT(IN)   :: ITC_MAX
!
! INOUT
!
      INTEGER,INTENT(INOUT) :: ICOUNT   
      INTEGER,INTENT(INOUT) :: INIT_EA
!
!  WORK 
!
      REAL(8)  ::  POT_LEAF     (IDX,ICG)
!
!     INTEGER,PARAMETER ::  ILP_MAX = 3
      INTEGER,PARAMETER ::  ILP_MAX = 4
      INTEGER  ::  ILP
      INTEGER  ::  IFLG 
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
! 
!  キャノピー周辺の水回りで起こっていることを確定する。
!
!  以下のものは決まっていない
!    ・ 気孔抵抗                  R_STM_BAR
!    ・ キャノピー空気の水蒸気圧  E_AIR
!    ・ 蒸散の on/off             FLAG_TR  ( , FR_H2O_CNP, FR_H2O_GRS は既知 )
!  上のいずれも、どれかを変えると他も変わる。---> イテレーション必要
!
!  ちなみに、
!      E_AIR , FLAG_TR --- CANOPY_STRES --> R_STM_BAR ( , FLAG_TR )  
!      R_STM_BAR   , FLAG_TR --- CANOPY_AIR   --> E_AIR 
!  また
!      FLAG_TR_RESET では、FLAG_TR を土壌温度条件のみ反映させてリセットする。
!
!  とりあえず、以下の手順とする。
!    1) 仮の、小さめの R_STM_BAR を初期値として与える。
!    2) 蒸散ありの状況で CANOPY_AIR で E_AIR を解く。
!    3) FLAG_TR をリセット
!    4) CANOPY_STRES で R_STM_BAR   を解き、FLAG_TR を得る
! 
!    5) CANOPY_AIR   で E_AIR を解く。
!       max 4 回解き、整合するまでやる。
!    6) FLAG_TR をリセット
!    7) CANOPY_STRES で R_STM_BAR   を解き、FLAG_TR を得る
!
!    8) CANOPY_AIR   で E_AIR を解く。
!       max 4 回解き、整合するまでやる。
!  収束性が悪いようなら、
!    ○前のステップでの E_AIR をとっておき、これを初期値として 3 から
!      スタート。
!
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_GET_EA 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF

! ==================
! >> まずは手始め <<
! ==================
!
! ================================
! >> 本格的なイテレーション開始 <<  
! ================================
!
! ------------------
! > E_AIR の計算 <
! ------------------
!
      IF ( INIT_EA .EQ. 1 ) THEN
        CALL CANOPY_AIR_ALL (
     I    NSTR       , NEND       , INSFLAG    , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , FR_H2O_GRS , 
     I    E_CNP_SAT  , E_REF      , E_SKIN_SAT   ,
     I    HUM_GRSK   , 
     O    FLAG_DRY_CNP , FLAG_DRY_SKIN         , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
        INIT_EA = 0
      ENDIF
!
      ILP     = 0 
 4999 CONTINUE
      ILP = ILP + 1 
      ICOUNT = ICOUNT + 1 
!
! ----------------------
! > 蒸散の方程式を解く <
! ----------------------
!     まずは、温度が負の層以外は全層蒸散ありとして見積り。
!     層毎の TR の計算もし、FLAG_TR を設定
!
      IF ( ILP.EQ.1 ) THEN
        CALL FLAG_TR_RESET (
     I     NSTR        , NEND        , INSFLAG , 
     I     TMP_SOIL    , 
!    M     FLAG_TR     , FLAG_TR_ALL , KLEV_TRNS_MAX  )
     M     FLAG_TR     , FLAG_TR_ALL )
      ENDIF
!
      IFLG = 0 
      IF ( ITC.EQ.ITC_MAX .AND. ILP.EQ.ILP_MAX ) IFLG = 1
!
      CALL  CANOPY_STRES (
     I   NSTR        , NEND        , INSFLAG      ,
     I   ITR         , ITR_MAX     , 
!    I   IFLG_AFT_MOD_FLG  , IFLG_LAST    , IFLG_MOD_FLG  , 
     I   1                 , IFLG         , 1             , 
     I   NTYP        , 
     I   R_AC_BAR        , R_AG        ,
     I   E_CNP_SAT   , E_SKIN_SAT    , E_AIR      , E_AIR_SAT      , 
     I   FUNC_TMP    , 
     I   POT_SOIL    ,
     I   D_ZERO_msx  , 
     I   R_ROOT_BAR_INV  ,  CVR_VEG_INV_msx   , 
     I   FAC_EVP     , R_STM_BAR_PAR    , 
     I   FR_H2O_CNP  , FR_H2O_GRS  ,
#ifdef MONYOS
     I   DELT_CNP , 
#endif
     O   POT_LEAF    , R_STM_BAR       , 
!    M   FAC_TRNS_A   , FAC_TRNS_B   , KLEV_TRNS_MAX , 
     M   FAC_TRNS_A   , FAC_TRNS_B   , 
     M   FLAG_TR     , FLAG_TR_ALL )
!
! ------------------
! > E_AIR の計算 <
! ------------------
!
      CALL CANOPY_AIR_ALL (
     I    NSTR       , NEND       , INSFLAG    , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , FR_H2O_GRS , 
     I    E_CNP_SAT  , E_REF      , E_SKIN_SAT   ,
     I    HUM_GRSK   , 
     O    FLAG_DRY_CNP , FLAG_DRY_SKIN         , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
!
      IF ( ILP .LT. ILP_MAX ) GOTO 4999
!
      RETURN
      END

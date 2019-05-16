!
! 大気強制力関連モジュール
!
!   仮定
!     ・JOB 実行は正時スタート
!     ・短波・長波計算は一時間に一度は行なう。
!
! 目的 : 
!   一次元メインプログラムから直接呼ばれる。
!   大気全面データを、ファイルから読み、時間内挿し、モデルに与える。
! 
! 必要な大気強制力データ
!   大気最下層東西風      (m/s)    snap
!   大気最下層南北風      (m/s)    snap
!   大気最下層気温        (K)      snap
!   大気最下層比湿        ( )      snap
!   大気最下層気圧        (hPa)    snap
!   地上気圧              (hPa)    snap
!   大規模降水量          (mm/day) 6hourly sum
!   対流性降水量          (mm/day) 6hourly sum
!   地表面下向き短波放射  (W/m**2) 6hourly avg
!   地表面下向き長波放射  (W/m**2) 6hourly avg
!   全雲量 CLD            (0-1)    snap?
!
! --------------- 大気最下層(=第一層)上端 --------------------HAF(2)
! |
! |--(pdel)    * 大気最下層代表値(uref vref tref qref pref)   FUL(1)
! |            |--(href)
! --------------- 地上気圧(psur)------------------------------HAF(1)
!
!
!  -----------------------------------------------------------
!   
      MODULE ATM
!
      USE PRM , ONLY :                    ! モデル
     1    IDIM  , 
     1    JDIM  , 
     1    JLPHY , 
     1    IJPHY , 
     1    IRAD  ,
     1    JRAD
!
      USE COM_RUNCONF_SIB0109 , ONLY :            ! ジョブ情報定数
!    1    JCN_SIB_MONIT_DETAIL     , 
     1    JCN_IWL_SKIP 
!
      USE ATM_FILE_UNIT , ONLY :       ! ファイル装置番号
     1    IFILE_UREF  , 
     1    IFILE_VREF  , 
     1    IFILE_TREF  , 
     1    IFILE_QREF   , 
     1    IFILE_PREF  , 
     1    IFILE_PSUR  , 
     1    IFILE_LPRC  , 
     1    IFILE_CPRC  , 
     1    IFILE_LWDN  , 
     1    IFILE_SWDN  , 
     1    IFILE_CLD
!
      use sib_monit , only : IMONIT_LEVEL , IMONIT_SELECT

	USE TS, ONLY:
     1     TIME_STEP, RSEC_CLEAN


!
      IMPLICIT NONE   
!
!
! SAVE 変数
!
      REAL(4),SAVE   :: U_PRE    (IDIM,JDIM)         
      REAL(4),SAVE   :: V_PRE    (IDIM,JDIM)         
      REAL(4),SAVE   :: P_PRE    (IDIM,JDIM)
      REAL(4),SAVE   :: PS_PRE   (IDIM,JDIM)
      REAL(4),SAVE   :: TMP_PRE  (IDIM,JDIM)
      REAL(4),SAVE   :: Q_PRE    (IDIM,JDIM)
!
      REAL(4),SAVE   :: U_NEXT   (IDIM,JDIM)
      REAL(4),SAVE   :: V_NEXT   (IDIM,JDIM)
      REAL(4),SAVE   :: P_NEXT   (IDIM,JDIM)
      REAL(4),SAVE   :: PS_NEXT  (IDIM,JDIM)
      REAL(4),SAVE   :: TMP_NEXT (IDIM,JDIM)
      REAL(4),SAVE   :: Q_NEXT   (IDIM,JDIM)
!
      REAL(4),SAVE   :: RLWD_PRE (IDIM,JDIM)
      REAL(4),SAVE   :: RLWD_NEXT(IDIM,JDIM)
!
      REAL(4),SAVE   :: RSWD_SUM (IDIM,JDIM)
      REAL(4),SAVE   :: PLRG_SUM (IDIM,JDIM)
      REAL(4),SAVE   :: PCNV_SUM (IDIM,JDIM)
!
      REAL(4),SAVE   :: CLD_PRE  (IDIM,JDIM)
      REAL(4),SAVE   :: CLD_NEXT (IDIM,JDIM)
!
      REAL(8),SAVE   :: SR_FLUX_WORK(IDIM,JDIM)
!
      LOGICAL,SAVE   :: LFIRST
!
      DATA LFIRST / .TRUE. / 
!
      CONTAINS 
!
! =====================================================================
!
      SUBROUTINE ATM_INI (
     I   ID_NOW     , DELT_ATM    , 
     I   RSEC       , 
     I   IMASK      , 
     I   GLON       , GLAT        ) 
!
! ATM の初期処理
!   前の読み込み段階で読んでおくべきデータを読み込む
!   仮定。現在 ID_NOW は正時
! 
      USE CALENDAR , ONLY : 
     1     CALENDAR_RUN_GETKT ,
     1     CALENDAR_RUN_GETID
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN) :: ID_NOW    (5)         ! 本ステップの時刻
      INTEGER,INTENT(IN) :: IMASK     (IDIM,JDIM)
      REAL(8),INTENT(IN) :: RSEC 
      REAL(8),INTENT(IN) :: DELT_ATM             ! モデルの時間刻
!
      REAL(8),INTENT(IN) :: GLON (IDIM,JDIM)     ! 経度（単位、度）
      REAL(8),INTENT(IN) :: GLAT (IDIM,JDIM)     ! 緯度（単位、度）
!
! WORK
!
      INTEGER :: ID_PRE     (5)        ! 前ステップの時刻
      INTEGER :: ID_PRE_00Z (5) 
      INTEGER :: ID_NEXT_00Z(5)  
      INTEGER :: ID_PRE_03Z (5) 
      INTEGER :: ID_NEXT_03Z(5)  
      INTEGER :: ID_WORK_1  (5) 
      INTEGER :: ID_WORK_2  (5) 
      REAL(8) :: RSEC_OLD 
      REAL(8) :: RDAY_OLD 
      REAL(8) :: RDAY_TOT
      INTEGER :: KT_LAG 
!
! ==============
! >> 初期処理 <<
! ==============
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MODULE ATM_MAIN 00/06/12 POCHAKA '
        LFIRST = .FALSE. 
      ELSE
        WRITE(6,*) 'MODULE ATM_INI ERROR : NOT FIRST TIME'
        stop 999
      ENDIF
!
!
!     SEC_6HR_00Z = MOD ( ID_NOW(4)   , 6 ) * 3600.D0 + DELT_ATM / 2 
!     SEC_6HR_LWDN = MOD ( ID_NOW(4)+6 , 3 ) * 3600.D0 + DELT_ATM / 2 
!
! ==================================
! >> ファイルからデータの読み込み <<
! ==================================
!
      KT_LAG = -1 
      CALL CALENDAR_RUN_GETID ( ID_NOW , ID_PRE ,  KT_LAG )
!
!  読み込み日時の設定
!
      CALL ATM_GETID_XX  ( ID_PRE , 0 , ID_PRE_00Z , ID_NEXT_00Z ) 
      CALL ATM_GETID_XX  ( ID_PRE , 3 , ID_PRE_03Z , ID_NEXT_03Z )  
!
!  00 + 6xn データの読み込み
!
      CALL ATM_FILE_OPEN (ID_PRE_00Z, 3)
      CALL ATM_FILE_OPEN (ID_PRE_00Z, 2)
      write(6,*) 'ATM_INI PREDATA  ' , ID_PRE_00Z
      CALL ATM_RDFL (IFILE_UREF,'UREF',ID_PRE_00Z ,IMASK,U_PRE   )
      CALL ATM_RDFL (IFILE_VREF,'VREF',ID_PRE_00Z ,IMASK,V_PRE   )
      CALL ATM_RDFL (IFILE_PSUR,'PSUR',ID_PRE_00Z ,IMASK,PS_PRE  )
      CALL ATM_RDFL (IFILE_PREF,'PREF',ID_PRE_00Z ,IMASK,P_PRE   )
      CALL ATM_RDFL (IFILE_TREF,'TREF',ID_PRE_00Z ,IMASK,TMP_PRE )
      CALL ATM_RDFL (IFILE_QREF,'QREF',ID_PRE_00Z ,IMASK,Q_PRE   )
      CALL ATM_RDFL (IFILE_CLD, 'CLD ',ID_PRE_00Z ,IMASK,CLD_PRE )
!
      CALL ATM_RDFL (IFILE_LPRC,'LPRC',ID_PRE_00Z ,IMASK,PLRG_SUM )     
      CALL ATM_RDFL (IFILE_CPRC,'CPRC',ID_PRE_00Z ,IMASK,PCNV_SUM )     
      CALL ATM_RDFL (IFILE_SWDN,'SWDN',ID_PRE_00Z ,IMASK,RSWD_SUM )     
!
      write(6,*) 'ATM_INI NEXTDATA ' , ID_NEXT_00Z
      CALL ATM_FILE_OPEN (ID_NEXT_00Z, 3)
      CALL ATM_FILE_OPEN (ID_NEXT_00Z, 2)
      CALL ATM_RDFL (IFILE_UREF,'UREF',ID_NEXT_00Z ,IMASK,U_NEXT   )
      CALL ATM_RDFL (IFILE_VREF,'VREF',ID_NEXT_00Z ,IMASK,V_NEXT   )
      CALL ATM_RDFL (IFILE_PSUR,'PSUR',ID_NEXT_00Z ,IMASK,PS_NEXT  )
      CALL ATM_RDFL (IFILE_PREF,'PREF',ID_NEXT_00Z ,IMASK,P_NEXT   )
      CALL ATM_RDFL (IFILE_TREF,'TREF',ID_NEXT_00Z ,IMASK,TMP_NEXT )
      CALL ATM_RDFL (IFILE_QREF,'QREF',ID_NEXT_00Z ,IMASK,Q_NEXT   )
      CALL ATM_RDFL (IFILE_CLD, 'CLD ',ID_NEXT_00Z ,IMASK,CLD_NEXT )
!
!  6HR PRE データの読み込み （長波）
!     例えば今が 01Z (あるいは 21Z-03Z のうちのどこか)であるとする。
!     ID_PRE_03Z 、_NEXT_03Z には値として 21Z と 03Z が入っている。
!
!     長波は、00Z のデータは 00Z-06Z（実効03Z）の
!     データなので、18Z（実効21Z） と 00Z （実効03Z）のデータを読めばよい。
!
!     つまり ID_PRE_03Z , _NEXT_03Z に -3 時間したデータを読む必要がある。
!     この -3 時間の操作は ATMP_RDFL が行なう。
!       
      write(6,*) 'ATM_INI  LWDN PRE  DATA  ' , ID_PRE_03Z
      CALL ATM_FILE_OPEN (ID_PRE_03Z, 1)
      CALL ATM_RDFL ( IFILE_LWDN,'LWDN',ID_PRE_03Z ,IMASK,RLWD_PRE ) 
!
      write(6,*) 'ATM_INI  LWDN NEXT DATA  ' , ID_NEXT_03Z
      CALL ATM_FILE_OPEN (ID_NEXT_03Z, 1)
      CALL ATM_RDFL ( IFILE_LWDN,'LWDN',ID_NEXT_03Z,IMASK,RLWD_NEXT) 
!
! ============================
! >> 短波重みテーブルの作成 <<
! ============================
!  注。正時スタートを前提としているので
!      短波・長波の 1 時間値は、ATM_INI では用意する必要なし。
!
!   一つ前の (00+6xn)Z での RDAY , RSEC の計算
!
!    RDAY_OLD の計算 ( = ID_PRE における RDAY )
!
      KT_LAG = 0 
      CALL CALENDAR_RUN_GETID ( ID_PRE , ID_WORK_1 , KT_LAG )
      ID_WORK_1(2) = 1 
      ID_WORK_1(3) = 1 
      ID_WORK_1(4) = 0 
      KT_LAG = 0 
      CALL CALENDAR_RUN_GETID ( ID_WORK_1 , ID_WORK_2 , KT_LAG )
      ID_WORK_2(1) = ID_WORK_1(1) + 1 
      CALL CALENDAR_RUN_GETKT ( ID_WORK_1 , ID_WORK_2 , 3 , KT_LAG )
      RDAY_TOT = KT_LAG * 1.D0
!
      CALL CALENDAR_RUN_GETKT ( ID_WORK_1 , ID_PRE    , 3 , KT_LAG )
      RDAY_OLD = KT_LAG * 1.D0 / RDAY_TOT
      write(6,*) 'ATM_INI RDAY_OLD =' , RDAY_OLD , RDAY_TOT , 
     1           KT_LAG,
     1           RDAY_OLD*365
!
!   RSEC_OLD の計算
!
      KT_LAG   = MOD ( ID_PRE(4) , 6 ) 
      RSEC_OLD = RSEC 
!
      IF ( RSEC*24 .LT. 5.99 ) THEN 
         RSEC_OLD = RSEC_OLD + 1.D0 
      ENDIF
!
      RSEC_OLD = RSEC_OLD - ( KT_LAG + 1 ) / 24.D0 
      write(6,*) 'ATM_INI RSEC_OLD =' , RSEC_OLD , RSEC 
      CALL RSEC_CLEAN ( RSEC_OLD ) 
      write(6,*) 'ATM_INI RSEC_OLD =' , RSEC_OLD , RSEC ,
     1           RSEC_OLD*86400. , 86400-RSEC_OLD*86400.
!
!   短波の計算
!
      CALL ATM_SWDN_6HR ( 
     I    RDAY_OLD  , RSEC_OLD , DELT_ATM , 
     I    GLON      , GLAT     , IMASK    , 
     M    RSWD_SUM  , 
     O    SR_FLUX_WORK         )
!
      RETURN
      END SUBROUTINE ATM_INI
! =====================================================================
      SUBROUTINE ATM_MAIN ( 
     I   ID_NOW     , ID_PRE     , DELT_ATM   ,
     I   RDAY       , RSEC        , 
     I   IMASK      , 
     I   GLON       , GLAT       , 
     O   U_PHY      , V_PHY       , 
     O   PD_PHY     , PS_PHY      , PF_PHY      , 
     O   TMP_PHY    , Q_PHY       , 
     O   ZMEAN_PHY  , ZTEMP_PHY   , 
     O   PPLI_PHY   , PPCI_PHY    ,
     O   RVISB , RVISD , RNIRB , RNIRD , DLWB  ) 
!
      USE COM_STEP_SIB0109 , ONLY :
     1    ICN_SIB0109_SW
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN) :: ID_NOW    (5)         ! 本ステップの時刻
      INTEGER,INTENT(IN) :: ID_PRE    (5)         ! 本ステップの時刻
      REAL(8),INTENT(IN) :: RSEC 
      REAL(8),INTENT(IN) :: RDAY 
      REAL(8),INTENT(IN) :: DELT_ATM
!
      REAL(8),INTENT(IN) :: GLON (IDIM,JDIM)     ! 経度（単位、度）
      REAL(8),INTENT(IN) :: GLAT (IDIM,JDIM)     ! 緯度（単位、度）
      INTEGER,INTENT(IN) :: IMASK(IDIM,JDIM)     
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: U_PHY     (IDIM,JDIM)  ! U
      REAL(8),INTENT(OUT)  :: V_PHY     (IDIM,JDIM)  ! V 
      REAL(8),INTENT(OUT)  :: PD_PHY    (IDIM,JDIM)  ! (PS-PH)*2  
      REAL(8),INTENT(OUT)  :: PS_PHY    (IDIM,JDIM)  ! ハーフレベル = Ps HPa
      REAL(8),INTENT(OUT)  :: PF_PHY    (IDIM,JDIM)  ! フルレベル HPa 
      REAL(8),INTENT(OUT)  :: TMP_PHY   (IDIM,JDIM)  ! 温度 
      REAL(8),INTENT(OUT)  :: Q_PHY     (IDIM,JDIM)  ! 比湿 kg/kg
      REAL(8),INTENT(OUT)  :: ZMEAN_PHY (IDIM,JDIM)  ! 一時間?平均天頂角
      REAL(8),INTENT(OUT)  :: ZTEMP_PHY (IDIM,JDIM)  ! 各ステップ天頂角
      REAL(8),INTENT(OUT)  :: PPLI_PHY  (IDIM,JDIM)  ! 大規模凝結性降水
      REAL(8),INTENT(OUT)  :: PPCI_PHY  (IDIM,JDIM)  ! 積雲性降水
!
      REAL(8),INTENT(OUT)  ::  RVISB ( IRAD*JRAD )   ! 放射計算時可視直達
      REAL(8),INTENT(OUT)  ::  RVISD ( IRAD*JRAD )   ! 放射計算時可視散乱
      REAL(8),INTENT(OUT)  ::  RNIRB ( IRAD*JRAD )   ! 放射計算時近赤直達
      REAL(8),INTENT(OUT)  ::  RNIRD ( IRAD*JRAD )   ! 放射計算時近赤散乱
      REAL(8),INTENT(OUT)  ::  DLWB  ( IRAD*JRAD )   ! 放射計算時長波
!      
! WORK
!
      INTEGER :: ID_PRE_00Z (5) 
      INTEGER :: ID_NEXT_00Z(5)  
      INTEGER :: ID_PRE_03Z (5) 
      INTEGER :: ID_NEXT_03Z(5)  
      INTEGER :: ID_WORK    (5)  
      INTEGER :: KT_LAG
      INTEGER :: IDUMMY_1
      INTEGER :: IDUMMY_2
      INTEGER :: IDUMMY_3
      INTEGER :: IDUMMY_4
      REAL(8) :: CLD   (IDIM,JDIM)
      REAL(8) :: RSHRT (IDIM,JDIM)
!
#ifdef MONYOS
      INTEGER           :: JL
      INTEGER           :: IPOINT
      CHARACTER( 7)     :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0
#endif
!
      REAL(8) :: WGT_00Z_NEXT 
      REAL(8) :: WGT_LWDN_NEXT 
!
      REAL(8),SAVE   :: SEC_6HR_00Z 
      REAL(8),SAVE   :: SEC_6HR_LWDN
!
      REAL(8) :: SEC_6HR_TOT
      INTEGER :: NSTEPS_RAD
!
      LOGICAL :: L_6HR_00Z      ! フラグ。(00+6xn)Z に true
      LOGICAL :: L_6HR_03Z      ! フラグ。(03+6xn)Z に true
!
      REAL(8),SAVE :: ZERO = 0.D0 
!
! ==============
! >> チェック <<
! ==============
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'ATM_MAIN ERROR : CALL ATM_INI'
        STOP 999
      ENDIF
!
! ====================
! >> フラグを立てる <<
! ====================
!
      L_6HR_00Z = .FALSE. 
      L_6HR_03Z = .FALSE. 
!
!   0,6,12,18 時になった
!
      IF ( ID_NOW(4).NE.ID_PRE(4) .AND. MOD(ID_NOW(4),6).EQ.0 ) THEN
        L_6HR_00Z = .TRUE. 
      ENDIF
!
!   3,9,15,21 時になった
!
      IF ( ID_NOW(4).NE.ID_PRE(4) .AND. MOD(ID_NOW(4),6).EQ.3 ) THEN
        L_6HR_03Z = .TRUE. 
      ENDIF
!
! =================
! >> (00 + 6xn)Z <<
! =================
!
      IF ( L_6HR_00Z ) THEN 
!
!  SAVE データの更新
!
        U_PRE    (:,:) = U_NEXT    (:,:)
        V_PRE    (:,:) = V_NEXT    (:,:)
        PS_PRE   (:,:) = PS_NEXT   (:,:)
        P_PRE    (:,:) = P_NEXT    (:,:)
        TMP_PRE  (:,:) = TMP_NEXT  (:,:)
        Q_PRE    (:,:) = Q_NEXT    (:,:)
        CLD_PRE  (:,:) = CLD_NEXT  (:,:)
!
        CALL ATM_GETID_XX  ( ID_NOW , 0 , ID_PRE_00Z , ID_NEXT_00Z ) 
!
        IF( ID_NOW(1).NE.ID_NEXT_00Z(1) ) THEN 
           CALL ATM_FILE_OPEN (ID_NEXT_00Z, 3)
        ENDIF
        IF( ID_NOW(1).NE.ID_PRE_00Z(1) ) THEN 
           CALL ATM_FILE_OPEN (ID_PRE_00Z, 2)
        ENDIF
        CALL ATM_RDFL (IFILE_UREF,'UREF',ID_NEXT_00Z ,IMASK,U_NEXT   )
        CALL ATM_RDFL (IFILE_VREF,'VREF',ID_NEXT_00Z ,IMASK,V_NEXT   )
        CALL ATM_RDFL (IFILE_PSUR,'PSUR',ID_NEXT_00Z ,IMASK,PS_NEXT  )
        CALL ATM_RDFL (IFILE_PREF,'PREF',ID_NEXT_00Z ,IMASK,P_NEXT   )
        CALL ATM_RDFL (IFILE_TREF,'TREF',ID_NEXT_00Z ,IMASK,TMP_NEXT )
        CALL ATM_RDFL (IFILE_QREF,'QREF',ID_NEXT_00Z ,IMASK,Q_NEXT   )
        CALL ATM_RDFL (IFILE_CLD, 'CLD ',ID_NEXT_00Z ,IMASK,CLD_NEXT )
!
        CALL ATM_RDFL (IFILE_LPRC,'LPRC',ID_PRE_00Z ,IMASK,PLRG_SUM )     
        CALL ATM_RDFL (IFILE_CPRC,'CPRC',ID_PRE_00Z ,IMASK,PCNV_SUM )     
        CALL ATM_RDFL (IFILE_SWDN,'SWDN',ID_PRE_00Z ,IMASK,RSWD_SUM )     
!
      ENDIF
!
! =================
! >> (03 + 6xn)Z <<
! =================
!   注意。
!    -1) 今 03Z であるとする。
!     0) 長波では、06Z とラベルの張ってあるデータは 00Z-06Z の平均である。
!        03Z ～ 09Z の間は 06Z と 12Z のラベルがあるデータを内挿する。
!        すなわち、03Z では 12Z とラベルが張ってあるデータを読み込み
!        これを RLWD_NEXT とする。
!     1) ということで、まず 、
!        06Z とラベルの張ってあるデータを RLWD_PRE にする
!     2) ISLSCP_GETID_XX で ID_NOW の近辺で 03Z,09Z,15Z,21Z の日時を取得。
!        ID_NOW が 03Z であるので、ID_NEXT_03Z は 09Z になる。
!     3) ISLSCP_RDFL で 'LWDN' については +3 時間して、12Z になる。   
!
      IF ( L_6HR_03Z ) THEN
        RLWD_PRE    (:,:) = RLWD_NEXT    (:,:)
!
        CALL ATM_GETID_XX  ( ID_NOW , 3 , ID_PRE_03Z , ID_NEXT_03Z ) 
        IF( ID_NOW(1).NE.ID_NEXT_03Z(1) ) THEN 
           CALL ATM_FILE_OPEN (ID_NEXT_03Z, 1)
        ENDIF
        CALL ATM_RDFL (IFILE_LWDN,'LWDN',ID_NEXT_03Z,IMASK,RLWD_NEXT)
      ENDIF
!
! ====================
! >> 経過時間の整形 <<
! ====================
!   ジョブを接続しても結果が変わらないようにする。
!   hour が繰り上がったら 3600 の倍数にする。
!
! 積算
      SEC_6HR_00Z = SEC_6HR_00Z + DELT_ATM 
! 値を整える
      IF ( ID_NOW(4).NE.ID_PRE(4) ) THEN        
        SEC_6HR_00Z = MOD ( ID_NOW(4)   , 6 ) * 3600.D0 + DELT_ATM / 2 
      ENDIF
!
!  長波放射用時間 
!
      SEC_6HR_LWDN = MOD ( ID_NOW(4)+3 , 6 ) * 3600.D0 + 1800.D0 
!
!
! ==========
! >> 内挿 <<
! ==========
!
!  ------------------
!  > 内挿重みの計算 <
!  ------------------
!
      SEC_6HR_TOT = 6 * 3600.D0 
      WGT_00Z_NEXT  = SEC_6HR_00Z / SEC_6HR_TOT 
      WGT_LWDN_NEXT = SEC_6HR_LWDN / SEC_6HR_TOT 
!
!  チェック
!
      IF ( WGT_00Z_NEXT.LE.0.D0 .OR. WGT_00Z_NEXT.GE.1.D0 ) THEN
        WRITE(6,*) 'ATM_MAIN ERROR WGT_00Z_NEXT = ' , WGT_00Z_NEXT
        STOP 999
      ENDIF
      IF ( WGT_LWDN_NEXT.LE.0.D0 .OR. WGT_LWDN_NEXT.GE.1.D0 ) THEN
        WRITE(6,*) 'ATM_MAIN ERROR WGT_LWDN_NEXT = ' , WGT_LWDN_NEXT
        STOP 999
      ENDIF
!
!  ==================
!  >> データの加工 <<
!  ==================
!
!  ------------
!  > 時間内挿 <
!  ------------
!
!   単純な 00Z 系物理量のみの線形内挿 
!     放射がらみ（短波、長波、雲量）は一時間に一度のみ行なう
!
      CALL ATM_INTERP ( U_PRE   , U_NEXT   , WGT_00Z_NEXT , U_PHY   )
      CALL ATM_INTERP ( V_PRE   , V_NEXT   , WGT_00Z_NEXT , V_PHY   )
      CALL ATM_INTERP ( PS_PRE  , PS_NEXT  , WGT_00Z_NEXT , PS_PHY  )
      CALL ATM_INTERP ( P_PRE   , P_NEXT   , WGT_00Z_NEXT , PF_PHY  )
      CALL ATM_INTERP ( TMP_PRE , TMP_NEXT , WGT_00Z_NEXT , TMP_PHY )
      CALL ATM_INTERP ( Q_PRE   , Q_NEXT   , WGT_00Z_NEXT , Q_PHY   )
      CALL ATM_INTERP ( CLD_PRE , CLD_NEXT , WGT_00Z_NEXT  , CLD     )
!
! ======================
! >> 短波データの作成 <<                              99/11/17 作成ここから
! ======================
!
!  ----------------------
!  > 6 時間に一度の処理 <
!  ----------------------
!
      IF ( L_6HR_00Z ) THEN
!
        CALL ATM_SWDN_6HR (
     I      RDAY      , RSEC     , DELT_ATM , 
     I      GLON      , GLAT     , IMASK    , 
     M      RSWD_SUM  , 
     O      SR_FLUX_WORK         )
!
      ENDIF
!
!  --------------------------
!  > 放射計算ステップの処理 <
!  --------------------------
!     次の短波放射計算までの平均の天頂角・短波放射 4 成分を計算する。
!
      IF ( ICN_SIB0109_SW .EQ. 1 ) THEN
!
        CALL ATM_SWDN_1HR (
     I    RDAY , RSEC , DELT_ATM  , 
     I    SR_FLUX_WORK,
     I    GLON , GLAT , 
     O    NSTEPS_RAD  , 
     O    RSHRT       ,  ZMEAN_PHY   )
!
!
!   雲量 CLD は 0-1 
        CALL ATM_SWDN_1HR_DEVIDE (
     O     RVISB , RVISD , RNIRB , RNIRD , 
     I     RSHRT , CLD   , ZMEAN_PHY   )
!
        CALL ATM_INTERP( RLWD_PRE, RLWD_NEXT, WGT_LWDN_NEXT , DLWB  )
!
      ENDIF
!
!  --------------------
!  > 毎ステップの処理 <
!  --------------------
!    毎ステップの天頂角の計算
!
      CALL ATM_SUNANG ( 
     I    RDAY     , RSEC    , 
     I    GLON     , GLAT    , 
     O    ZTEMP_PHY    )
!
!  ==================
!  >> その他の加工 <<
!  ==================
!
!   ZTEMP 補正 （負の値を消す）
      ZTEMP_PHY(:,:) = MAX ( ZTEMP_PHY(:,:) , ZERO ) 
!
      PPCI_PHY(:,:) =  PCNV_SUM(:,:)           
     1                    / (3600*24) * DELT_ATM 
      PPLI_PHY(:,:) =  PLRG_SUM(:,:)           
     1                    / (3600*24) * DELT_ATM   
!
!   --------
!   > 設定 <
!   --------
!     大気最下層層厚 を 10HPa に 
      PD_PHY(:,:) = 10.D0
!
! ================
! >> 詳細モニタ <<   
! ================
!
#ifdef MONYOS
!     IF ( JCN_SIB_MONIT_DETAIL .EQ. 1 ) THEN
      IF ( IMONIT_LEVEL .GE. IMONIT_SELECT ) THEN
        DO JL=1,JLPHY
          IPOINT = (JL-1)*IJPHY + 1              ! ポインタ。
!
          IF ( ICN_SIB0109_SW .EQ. 1 ) THEN
!
! 全短波
            CMARK = 'FSR' 
            CALL MONIT_ADD_2 ( CMARK , RSHRT (IPOINT,1) , 
     1                         JL, DELT_ATM*NSTEPS_RAD , ONE )     
!
! 長波
            CMARK = 'FLR' 
            CALL MONIT_ADD_2 ( CMARK , DLWB   (IPOINT) , 
     1                         JL, DELT_ATM*NSTEPS_RAD , ONE )
!
          ENDIF
!
! 温度
          CMARK = 'FTMP' 
          CALL MONIT_ADD_2 ( CMARK , TMP_PHY(IPOINT,1) , 
     1                       JL , DELT_ATM , ONE )
!
! 比湿
          CMARK = 'FQ' 
          CALL MONIT_ADD_2 ( CMARK , Q_PHY(IPOINT,1) , 
     1                       JL , DELT_ATM , ONE )
!
! 地表面気圧
          CMARK = 'FPS' 
          CALL MONIT_ADD_2 ( CMARK , PS_PHY(IPOINT,1) , 
     1                       JL , DELT_ATM , ONE )
!
! 積雲性降水
          CMARK = 'FPC' 
          CALL MONIT_ADD_2 ( CMARK , PPCI_PHY(IPOINT,1) , 
     1                       JL , ONE , ONE )
!
! 大規模凝結性降水
          CMARK = 'FPL' 
          CALL MONIT_ADD_2 ( CMARK , PPLI_PHY(IPOINT,1) , 
     1                       JL , ONE , ONE )
!
! 風速
          CMARK = 'FU' 
          CALL MONIT_ADD_2 ( CMARK , U_PHY(IPOINT,1) , 
     1                       JL , DELT_ATM , ONE )

        ENDDO
      ENDIF
#endif
!
! ==============
! >> 初期処理 <<
! ==============
!
      IF ( LFIRST ) THEN
        LFIRST = .FALSE. 
      ENDIF
!
      RETURN
      END SUBROUTINE ATM_MAIN
!
!======================================================
      SUBROUTINE MONIT_REGIST_ATM_MAIN
!
      USE SIB_MONIT, ONLY : IMONIT_LEVEL , IMONIT_SELECT
!
! ATM_MAIN で出力する変数リスト
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
        WRITE(6,*) 'MONIT_REGIST_ATM_MAIN 99/11 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SELECT ) THEN
!
      CMARK  = 'FSR'
      CTITLE = 'SHORT RADIATION (ATMOSPHERIC FORCING)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FLR'
      CTITLE = 'LONG RADIATION (ATMOSPHERIC FORCING)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FTMP'
      CTITLE = 'TEMPERATURE (ATMOSPHERIC FORCING)'
      CUNIT  = 'K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FQ'
      CTITLE = 'HUMIDITY (ATMOSPHERIC FORCING)'
      CUNIT  = ''
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FPS'
      CTITLE = 'SURFACE PRESSURE (ATMOSPHERIC FORCING)'
      CUNIT  = 'HPa'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FPC'
      CTITLE = 'CONVECTIVE PRECIPITATION (ATMOSPHERIC FORCING)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FPL'
      CTITLE = 'LARGE SCALE PRECIPITATION (ATMOSPHERIC FORCING)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FCLD'
      CTITLE = 'CLOUD AMOUNT (ATMOSPHERIC FORCING)'
      CUNIT  = ''
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FU'
      CTITLE = 'WIND VELOCITY (ATMOSPHERIC FORCING)'
      CUNIT  = 'M/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      end  subroutine monit_regist_atm_main
!
      END MODULE ATM








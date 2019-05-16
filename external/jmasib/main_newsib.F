! デバッグ用メインファイル 99/08/06 バージョン 
! 99/07/xx 保坂征宏
! 99/08/06 保坂征宏 大泉さんに渡すバージョンフィックス
!
      PROGRAM RUN_NEW   
!
      USE SIB0109 , ONLY :
     1  SIB0109_INI          , 
     1  SIB0109_RUN_MONITOR_INITIAL ,
     1  SIB0109_RUN_STEPINI  ,
     1  SIB0109_RUN_LOOPINI  ,
     1  SIB0109_RUN_ALBEDO   , 
     1  SIB0109_RUN_SR_CALC  ,  
     1  SIB0109_RUN_SR_RESTORE , 
     1  SIB0109_RUN_LR_SET   , 
     1  SIB0109_RUN_SIBMAIN       ,
     1  SIB0109_RUN_LOOPEND     ,
     1  SIB0109_RUN_STEPEND 
!
! 鉛直一次元モデル ～ 新 SiB のデバッグ用
!
      USE PRM , ONLY : 
     1   IJPHY , 
     1   JLPHY , 
     1   INTRI , 
     1   INTRJ ,
     1   IDIM  , 
     1   JDIM  , 
     1   JMAX  , 
     1   ISPT  ,
     1   KMP1 
!
!
      USE ISLSCP , ONLY :
     1   ISLSCP_INI , 
     1   ISLSCP_MAIN
!
      USE ATM , ONLY :
     1   ATM_INI , 
     1   ATM_MAIN
!
      USE GSWP2 , ONLY :
     1   GSWP2__INI , 
     1   GSWP2__RUN
!
      USE FORCING, ONLY:
     1   FORCING_SETUP,
     2   FORCING_GET_GEOGRAPHY,
     3   FORCING_MAIN
!
      USE COM_STEP , ONLY :
!    1   ICNSW                     ,      ! mj98 com_step 
     1   ICNMNTMON                 ,      ! mj98 com_step
     1   ICNLAST                   ,       ! mj98 com_step 
     1   SCN_DELT

      USE COM_STEP_SIB0109 , only : 
     1   ICN_SIB0109_CALC_SOIL_SNOW 
!
      use com_debug , only :
     1   istep , com_debug_ini
!
      USE COM_TETEN_sib_0109 , ONLY :
     1    COM_TETEN_sib_0109_INI   
!
      USE CALENDAR , ONLY :
     1    CALENDAR_INI   ,
     1    CALENDAR_RUN_GETID ,
     1    CALENDAR_RUN_GETKT
!
      USE COM_FILEUNIT , ONLY : 
     I   IMONIT_2CFL ,IMONIT_ZCFL ,IMONIT_3CFL ,
     I   IMONIT_ZETACFL ,IMONIT_3ETACFL ,
     I   IMONIT_6HRCFL ,IMONIT_DAYCFL ,IMONIT_GLOBCFL,
     1   imonit_2fl , imonit_2fl8 , imonit_sfl
!
      use com_stdinout_unit , only : com_stdinout_unit_ini
     
      USE COM_RUNCONF , ONLY : 
     1     COM_RUNCONF_INI ,
     1     JCNIMNT
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1     COM_RUNCONF_SIB0109_INI ,
     1     JCN_SIB_ISLSCP , 
     1     JCN_IWL_SKIP   ,
     1     RCN_DELT_ATM   ,
     1     RCN_DELT_CNP   ,
     1     RCN_DELT_SOIL  ,
     1     RCN_DELT_SNOW   
!
      USE COM_JOBINFO_SIB0109 , ONLY :
     1       COM_JOBINFO_SIB0109_INI   ,
     1       IDSTAR            , 
     1       IDATE             , 
     1       IDEND             , 
     1       IOUT_8BYTE        ,
     1       INI_MON           , 
     1       CFILE_MONIT       ,  
     1       CDIR_MONIT  
     
      USE TS, ONLY:
     1     TIME_STEP, RSEC_CLEAN
     
!
      USE MESSAGE, ONLY: MESSAGE_INI, MESSAGE_FLUSH
!
      IMPLICIT NONE
!
! ------------
! > 状況設定 <
! ------------
!
!
      INTEGER   :: IMASK            ( ISPT*IJPHY, JLPHY ) ! 植生・海陸情報
      REAL(8)   :: GLON             ( IDIM, JDIM ) 
      REAL(8)   :: GLAT             ( IDIM, JDIM ) 
#ifdef MONYOS
!     吉村モニタは緯度経度方向の直線格子を想定している
      real(8)   :: GLON_AVERAGE(IDIM)
      real(8)   :: GLAT_AVERAGE(JDIM)
#endif
!
      INTEGER   :: JL                            ! 1 
      INTEGER   :: MON
      REAL(8)   :: DELT_MODEL
!     REAL(8)   :: DELT_SOIL
!
      integer   :: id_pre5
      integer   :: id_now5  
      real(8),parameter   :: rsec_tot = 86400.D0
!
! ----------------------------
! > SAVE 変数 : 大気放射関連 <
! ----------------------------
!
      REAL(8)  :: RNIRB ( IJPHY/INTRI/INTRJ , JLPHY )  
      REAL(8)  :: RNIRD ( IJPHY/INTRI/INTRJ , JLPHY )  
      REAL(8)  :: RVISB ( IJPHY/INTRI/INTRJ , JLPHY )  
      REAL(8)  :: RVISD ( IJPHY/INTRI/INTRJ , JLPHY )  
      REAL(8)  :: DLWB  ( IJPHY/INTRI/INTRJ , JLPHY )  
!
! ---------------------
! > ワーク : 大気関連 <
! ---------------------
!
!   入力
!
      REAL(8)  :: GMT_PHY   (IJPHY,2)
      REAL(8)  :: GMQ_PHY   (IJPHY,2)
      REAL(8)  :: GMUV_PHY  (IJPHY,3)
!
!   出力
!      
      REAL(8)  :: TMTX2L_PHY  (IJPHY)
      REAL(8)  :: QMTX2L_PHY  (IJPHY)
      REAL(8)  :: UMTX2L_PHY  (IJPHY)
      REAL(8)  :: VMTX2L_PHY  (IJPHY)
!
      REAL(8)  :: AVISB   (IJPHY)
      REAL(8)  :: ANIRB   (IJPHY)
      REAL(8)  :: AVISD   (IJPHY)
      REAL(8)  :: ANIRD   (IJPHY)
!
!     REAL(8)  :: U_STR   (IJPHY)
!     REAL(8)  :: V_STR   (IJPHY)
!
      REAL(8)  ::   RAD_LONG_SIB_2_ATM_PHY (IJPHY) 


!
! ----------------
! > 大気情報など <
! ----------------
!
      REAL(8)     ::  ZMEAN_PHY (IJPHY,JLPHY)
      REAL(8)     ::  ZTEMP_PHY (IJPHY,JLPHY)
      REAL(8)     ::  DAYTIME_1HR_PHY (IJPHY,JLPHY) = 3600.D0 
!
      REAL(8)     ::  PPLI_PHY  (IJPHY,JLPHY)
      REAL(8)     ::  PPCI_PHY  (IJPHY,JLPHY)
!
      REAL(8)     ::  GT_PHY    (IJPHY,JLPHY)
      REAL(8)     ::  GQ_PHY    (IJPHY,JLPHY)
      REAL(8)     ::  GU_PHY    (IJPHY,JLPHY)
      REAL(8)     ::  GV_PHY    (IJPHY,JLPHY)
      REAL(8)     ::  GPFUL_PHY (IJPHY,JLPHY)
      REAL(8)     ::  GPHAF_PHY (IJPHY,JLPHY)
      REAL(8)     ::  GPDEL_PHY (IJPHY,JLPHY)
!
      INTEGER     :: KTSTAR
      INTEGER     :: KTEND
!
      INTEGER :: ID_PRE (5)
      INTEGER :: ID_NOW (5)
      INTEGER :: ID_NEXT(5)
!
      REAL(8) :: RSEC       ! 1 日での現在の位置。
      REAL(8) :: RDAY       ! 1 年での現在の位置。短波計算で必要。
!
      LOGICAL :: LFIRST  
!
#ifdef MONYOS
! ------------
! > モニタ用 < for monitor
! ------------
!
      CHARACTER(100) ::
     I  CMONIT_2_NAME, CMONIT_Z_NAME, CMONIT_3_NAME,
     I  CMONIT_ZETA_NAME, CMONIT_3ETA_NAME,
     I  CMONIT_6HR_NAME, CMONIT_DAY_NAME, CMONIT_GLOB_NAME
!
      CHARACTER(8)  :: MODEL
      CHARACTER(11) :: RESL    
      CHARACTER(4)  :: EXPR
!     CHARACTER(1)  :: EXPR1
!     CHARACTER(2)  :: EXPR2
      CHARACTER(80) :: CINF(10)
      REAL(8) :: GW          (JMAX)
      REAL(8) :: COSCLT      (JMAX)
      REAL(8) :: PA(KMP1)
      REAL(8) :: PB(KMP1)    
!
      integer  :: k 
      real(8)  :: TOTMON 
      REAL(8)  :: TOTMON_X 
      INTEGER   :: J
#endif
!
      REAL(8),PARAMETER  :: CONS_MTX_1 = 1.D0   ! いらない?? 
      REAL(8),PARAMETER  :: CONS_MTX_2 = 0.D0
      REAL(8)  :: FSECP 
      REAL(8)  :: FSEC0 
      INTEGER  :: KTP 
      INTEGER  :: KT0 
      INTEGER  :: KTM 
!
      REAL(8)   :: DUMMY_SCALAR = 0.D0 

!
      DATA LFIRST / .TRUE. /
!
! ==============
! >> 初期設定 <<
! ==============

! ------------------------      
! > ジョブ条件などの設定 <
! ------------------------
!
!  モデル設定など
!
!     IOUT_8BYTE = 1             ! 8 バイトファイル出力する
!
      CALL COM_STDINOUT_UNIT_INI(5, 6)
      call com_debug_ini
      CALL MESSAGE_INI
      CALL COM_RUNCONF_INI
      CALL COM_RUNCONF_SIB0109_INI
      CALL COM_JOBINFO_SIB0109_INI
!     CALL TETEN
      CALL COM_TETEN_sib_0109_INI
      CALL CALENDAR_INI
!
      CALL GET_KTINFO (                ! ここではインチキくさい方法をとる
     I    IDSTAR , 
     O    IDATE  , KTP    , KT0    , FSECP ) 
!
      CALL CALENDAR_RUN_GETKT ( IDATE  , IDSTAR  , 4   , KTSTAR )
      CALL CALENDAR_RUN_GETKT ( IDATE  , IDEND   , 4   , KTEND  )
!
!
      DELT_MODEL= RCN_DELT_CNP   !  compdt での時間刻計算に対応
!     DELT_SOIL = RCN_DELT_SOIL
!
      SCN_DELT  = RCN_DELT_CNP   !  初期の  monit_out_2  のため
!
      MON       = IDSTAR(2) 
!
      CALL DELT_CHECK 
     I   ( RCN_DELT_ATM , RCN_DELT_CNP , RCN_DELT_SNOW , RCN_DELT_SOIL )

!
! ------------------------
! > 強制ファイルオープン <
! ------------------------
!
      IF      ( JCN_SIB_ISLSCP .EQ. 0 ) THEN
        WRITE(6,*) 'TEST1DIM SO NO FORCING FILE IS OPENED' 
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 1 ) THEN
        WRITE(6,*) 'OPEN ISLSCP_FILE_OPEN' 
        CALL ISLSCP_FILE_OPEN 
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 2 ) THEN
	CONTINUE
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 3 ) THEN
        WRITE(6,*) 'OPEN FORCING_FILE_OPEN' 
	CALL FORCING_SETUP
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 4 ) THEN
        WRITE(6,*) 'GSWP2 ' 
      ELSE 
        WRITE(6,*) 'MAIN2 ERROR JCN_SIB_ISLSCP = ' , JCN_SIB_ISLSCP 
        STOP 999
      ENDIF
!
!
! ------------------------
! > 陸・植生の分布の決定 <
! ------------------------
!
!  TEST1D 用 
!
      IF ( JCN_SIB_ISLSCP .EQ. 0 ) THEN
!
        DO JL=1,JLPHY
          CALL TEST1D_CONFIG_SET (
     I      JL    ,     
     O      GLON  , GLAT  , 
     O      IMASK   ) 
        ENDDO 
!
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 1 ) THEN
!
        CALL ISLSCP_CONFIG_SET ( 
     O    GLON , GLAT , IMASK  ) 
!
!      
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 2 ) THEN
!
        WRITE(6,*) 'OPEN ATM_CONFIG_SET'
        CALL ATM_CNFILE_OPEN
        CALL ATM_CONFIG_SET ( 
     O    GLON , GLAT , IMASK  ) 
!      
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 3 ) THEN
!
        WRITE(6, *) 'OPEN FORCING: GEOGRAPHY'
        CALL FORCING_GET_GEOGRAPHY(
     O    OUT_GLON=GLON, OUT_GLAT=GLAT, OUT_IMASK=IMASK)
!
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 4 ) THEN
!
        WRITE(6, *) 'MAIN_NEWSIB : CALL GSWP2__INI '
        CALL GSWP2__INI ( 
     I    RCN_DELT_ATM , idstar , 
     O    IMASK , GLON , GLAT ) 
!
      ENDIF
!
      call mod_imask ( imask , ispt, idim , jdim )
!
      CALL SIB0109_INI ( IMASK , IDSTAR ) 
!
!
#ifdef MONYOS
!
! --------------------
! > 吉村モニタ初期化 <
! --------------------
!
      MODEL = 'MJSIB907'
      RESL  = '1DIMversion'
!     CHARACTER(1)  :: EXPR1
!     CHARACTER(2)  :: EXPR2
      DO J=1,JMAX
        GW(J)      = 1.
        COSCLT (J) = (JMAX-2*J+1.) / JMAX  
      ENDDO
      EXPR = 'TEST'  
!
      do 2 k=1,10
        cinf(k) = ' '
 2    continue    
!
      write(6,*) 'main : jcn_iwl_skip , jcnimnt original ' , 
     1           jcn_iwl_skip , jcnimnt   
      IF ( JCN_IWL_SKIP.LE.-1 ) THEN     ! 毎ステップ出力
        IF ( JCNIMNT .NE. -1 ) THEN 
          WRITE(6,*) 'MAIN WARNING ' , 
     1               ' JCNIMNT IS MODIFIED ' , JCNIMNT , ' TO ' , -1
        ENDIF 
        jcnimnt = -1        
      ELSE                               ! JCNIMNT 時間おき出力
        IF ( JCNIMNT .LT. 0 ) THEN 
          WRITE(6,*) 'MAIN WARNING ' , 
     1               ' JCN_IWL_SKIP IS MODIFIED ' 
          JCN_IWL_SKIP = - JCN_IWL_SKIP 
        ENDIF 
      ENDIF
      write(6,*) 'main : jcn_iwl_skip , jcnimnt ' , 
     1           jcn_iwl_skip , jcnimnt   
!
      DO 1 K=1,KMP1
        PA(K) = 0.
        PB(K) = 0.
 1    CONTINUE  
!
! 吉村モニタファイル名の設定
!
      IF ( CFILE_MONIT == ' ' ) then
        CMONIT_2_NAME = 
     1       'newsib_monit_YYYY_MM_DD_HH_' // 
     1                    'YYYY_MM_DD_HH'
        CALL REPLACE_INT(CMONIT_2_NAME, 'YYYY', IDSTAR(1))
        CALL REPLACE_INT(CMONIT_2_NAME, 'MM', IDSTAR(2))
        CALL REPLACE_INT(CMONIT_2_NAME, 'DD', IDSTAR(3))
        CALL REPLACE_INT(CMONIT_2_NAME, 'HH', IDSTAR(4))
        CALL REPLACE_INT(CMONIT_2_NAME, 'YYYY', IDEND(1))
        CALL REPLACE_INT(CMONIT_2_NAME, 'MM', IDEND(2))
        CALL REPLACE_INT(CMONIT_2_NAME, 'DD', IDEND(3))
        CALL REPLACE_INT(CMONIT_2_NAME, 'HH', IDEND(4))
      ELSE 
        CMONIT_2_NAME = CFILE_MONIT
      ENDIF
!
! 吉村モニタファイルオープン
!
      open ( imonit_2fl  , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'.dr',
     1        access='direct' , recl=4*IJPHY*JLPHY , 
     1        form  ='unformatted' ) 
      if ( iout_8byte.eq.1 ) then
        open ( imonit_2fl8 , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'_8.dr',
     1        access='direct' , recl=8*IJPHY*JLPHY , 
     1        form='unformatted' ) 
      endif
      open ( imonit_2cfl , 
     1        file  = TRIM(CDIR_MONIT) // TRIM(cmonit_2_name)//'.ctl')
!
! 吉村モニタ初期化 
!
      call monit_ini
!
      glon_average(:) = sum(glon, dim=2) / JDIM
      glat_average(:) = sum(glat, dim=1) / IDIM
      call monit_grads_ctl (
     I   IMONIT_2CFL ,IMONIT_ZCFL ,IMONIT_3CFL ,
     I   IMONIT_ZETACFL ,IMONIT_3ETACFL ,
     I   IMONIT_6HRCFL ,IMONIT_DAYCFL ,IMONIT_GLOBCFL,
     I   CMONIT_2_NAME ,CMONIT_Z_NAME ,CMONIT_3_NAME ,
     I   CMONIT_ZETA_NAME ,CMONIT_3ETA_NAME ,
     I   CMONIT_6HR_NAME ,CMONIT_DAY_NAME ,CMONIT_GLOB_NAME,
     I   cosclt, idstar, idend, ktstar, ktend,         
     I   model, resl, gw, pa, pb , glon_average , glat_average)
      call monit_clear(totmon) 
!
#endif
!
! ----------------
! > SiB 初期設定 <
! ----------------
!    現状、初期値モニタ出力が抜けている
!
!     CALL SIB0109_INI ( IMASK , FSECP , KTP ) 
!
#ifdef MONYOS
      IF ( INI_MON .EQ.1 ) THEN
        write(6,*) 'initial snap shot monitoring ' , delt_model
!
        CALL SIB0109_RUN_MONITOR_INITIAL ( MON , KT0 )   
!    
        TOTMON_X = DELT_MODEL            ! 積算時間 ..  モニタをだます
!
        ICN_SIB0109_CALC_SOIL_SNOW = 1   ! 雪土壌も計算したステップであると、
                                         ! モニタをだます
!
        CALL MONIT_OUT_2
     1  ( IMONIT_2FL       , IMONIT_2FL8   , IMONIT_SFL    ,
     1    IDATE , MODEL    , RESL  , EXPR  , CINF  ,
     1    KTP   , TOTMON_X , FSECP , IMASK , GW )
!
        CALL MONIT_CLEAR  ( TOTMON_X )
!
        write(6,*) 'initial snap shot monitoring end ' 
!
      ENDIF
#endif 
!
! ============================
! >> 時間積分ループここから <<
! ============================
!
! カウンタリセット
!
      ISTEP = 0
!
      id_now5 = 3600.D0 + 0.01 - delt_model
      if ( delt_model .gt. 3600.D0+0.01 ) then
        write(6,*) 'main_newsib : delt_model gt 3600 delt_model = ' , 
     1              delt_model
        stop 999
      endif
!
 8000 CONTINUE
!
!  時計を進める。ここでフラグもたてていることに注意。
!
      CALL TIME_STEP ( 
     I   DELT_MODEL , DUMMY_SCALAR  , IDSTAR   , IDEND    , 
     M                FSEC0    , FSECP    , 
     M   KTM        , KT0      , KTP      , 
     M   ISTEP      ,
     M   IDATE      , ID_NEXT  , ID_NOW   , ID_PRE   , 
     M   RDAY       , RSEC     
#ifdef MONYOS
     M   , TOTMON 
#endif
     1   ) 
!
      write(6,*) 'main_newsib istep' , istep
!
      mon = id_now(2)
!
      id_pre5   = id_now5 
      id_now5   =                      ! 正時からの秒数
     1              rsec * rsec_tot    ! 1 日のはじめからの秒数
     1            - id_now(4)*3600.D0  
!
      CALL SIB0109_RUN_STEPINI (  FSECP ) 
!
!  -----------------------
!  > ISLSCP 大気条件設定 <
!  -----------------------
!
      IF ( JCN_SIB_ISLSCP .EQ. 1 ) THEN
!
        IF ( LFIRST ) THEN
          CALL ISLSCP_INI (
     I     ID_NOW     , DELT_MODEL   , 
     I     RSEC       , 
     I     IMASK      , 
     I     GLON       , GLAT        ) 
        ENDIF
!
        CALL ISLSCP_MAIN (
     I   ID_NOW     , ID_PRE      , DELT_MODEL    , 
     I   RDAY       , RSEC        , 
     I   IMASK      , 
     I   GLON       , GLAT        ,
     O   GU_PHY     , GV_PHY      , 
     O   GPDEL_PHY  , GPHAF_PHY   , GPFUL_PHY   , 
     O   GT_PHY     , GQ_PHY      , 
     O   ZMEAN_PHY  , ZTEMP_PHY   , 
     O   PPLI_PHY   , PPCI_PHY    ,
     O   RVISB      , RVISD       , RNIRB        , RNIRD       , 
     O   DLWB     ) 
!
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 2 ) THEN
!
        IF ( LFIRST ) THEN
          CALL ATM_INI (
     I     ID_NOW     , DELT_MODEL   , 
     I     RSEC       , 
     I     IMASK      , 
     I     GLON       , GLAT        ) 
        ENDIF
!
        CALL ATM_MAIN (
     I   ID_NOW     , ID_PRE      , DELT_MODEL   , 
     I   RDAY       , RSEC        , 
     I   IMASK      ,
     I   GLON       , GLAT        ,       
     O   GU_PHY     , GV_PHY      , 
     O   GPDEL_PHY  , GPHAF_PHY   , GPFUL_PHY   , 
     O   GT_PHY     , GQ_PHY      , 
     O   ZMEAN_PHY  , ZTEMP_PHY   , 
     O   PPLI_PHY   , PPCI_PHY    ,
     O   RVISB      , RVISD       , RNIRB        , RNIRD       , 
     O   DLWB     ) 
!
      ELSE IF (JCN_SIB_ISLSCP == 3) THEN
!
        CALL FORCING_MAIN (
     I   ID_NOW, ID_PRE, id_now5, id_pre5, 
     I   DELT_MODEL, RDAY, RSEC, IMASK, GLON, GLAT,
     O   GU_PHY, GV_PHY, GPDEL_PHY, GPHAF_PHY, GPFUL_PHY,
     O   GT_PHY, GQ_PHY, 
     O   ZMEAN_PHY, ZTEMP_PHY, PPLI_PHY, PPCI_PHY,
     O   RVISB, RVISD, RNIRB, RNIRD, DLWB) 
!
      ELSE IF ( JCN_SIB_ISLSCP .EQ. 4 ) THEN
!
        GV_PHY(:,:) = 0. 
        CALL GSWP2__run  (
     I   ID_NOW     , RDAY       , RSEC        , 
     O   ZMEAN_PHY  , ZTEMP_PHY   , DAYTIME_1HR_PHY , 
     O   RVISB      , RVISD       , RNIRB        , RNIRD       , 
     O   DLWB       ,  
     O   PPLI_PHY   , PPCI_PHY    ,
     O   GT_PHY     , GQ_PHY      , 
     O   GPHAF_PHY  , GPDEL_PHY  , GPFUL_PHY   , 
     O   GU_PHY         ) 
!
      ENDIF
!
!===============================================================
!
! PHYSCS JL ループ開始
!
      DO 7000 JL=1,JLPHY
!
!     WRITE(6,*) 'JL =' , JL 
!
! 大気の値の決定 （物理格子）
!
      IF ( JCN_SIB_ISLSCP .EQ. 0 ) THEN
        CALL TEST1D_ATMOS_SET (
     I    JL            , 
     I    DELT_MODEL    , 
     O    GPDEL_PHY     , GPHAF_PHY     , GPFUL_PHY      , 
     O    GT_PHY        , GQ_PHY        , 
     O    GU_PHY        , GV_PHY        , 
     O    ZMEAN_PHY     , ZTEMP_PHY     , 
     O    PPLI_PHY      , PPCI_PHY      , 
     O    GMT_PHY       , GMQ_PHY       , GMUV_PHY        )
      ELSE                                    ! ここはいいのかな ？
        CALL CON_SET ( GMT_PHY (1,1) ,   IJPHY , CONS_MTX_1 )
        CALL CON_SET ( GMQ_PHY (1,1) ,   IJPHY , CONS_MTX_1 )
        CALL CON_SET ( GMUV_PHY(1,1) ,   IJPHY , CONS_MTX_1 )
        CALL CON_SET ( GMT_PHY (1,2) ,   IJPHY , CONS_MTX_2 )
        CALL CON_SET ( GMQ_PHY (1,2) ,   IJPHY , CONS_MTX_2 )
        CALL CON_SET ( GMUV_PHY(1,2) , 2*IJPHY , CONS_MTX_2 )
      ENDIF
!
!  LOOP はじめに sib0109 がすべきこと
!     ・IMSKLT・NMAX の設定
!     ・大気情報を SiB 格子に詰める
!     ・SiB 情報をバンド SiB 格子に詰める 
!
      CALL SIB0109_RUN_LOOPINI ( 
     I    JL                 , MON                , 
     I    KT0                , 
     I    GPDEL_PHY(1,JL)    , GPHAF_PHY(1,JL)    , GPFUL_PHY(1,JL)   , 
     I    GT_PHY   (1,JL)    , GQ_PHY   (1,JL)    , 
     I    GU_PHY   (1,JL)    , GV_PHY   (1,JL)    , 
     I    ZMEAN_PHY(1,JL)    , ZTEMP_PHY(1,JL)    ,
     I    PPLI_PHY (1,JL)    , PPCI_PHY (1,JL)    , 
     I    DAYTIME_1HR_PHY(1,JL) ,  
     I    GMT_PHY            , GMQ_PHY            , GMUV_PHY        ) 
!
! ------------------------------------------------------------------------
!          ここまでが、mj の PHYSCS での ZENITH, YMPBL0, SETBND 等にあたる
! ------------------------------------------------------------------------
!
! TEST1DIM の場合、放射の結果
!
      IF ( JCN_SIB_ISLSCP .EQ. 0 ) THEN
        CALL TEST1D_RAD_SET (
     I   JL   ,   
     O   RVISB , RVISD , RNIRB , RNIRD , DLWB )
      ENDIF
!
! ====================
! >> SiB + 放射計算 <<
! ====================
!  本来の順番は
!    ・放射計算するなら
!        アルベド計算  
!        放射計算
!        SiB で必要な計算
!        必要な値をセーブ 
!    ・しないなら
!        必要な値をリストア
!    ・放射補正など、共通作業
!    ・陸面本格計算
!  である。
!
!  OFFLINE では『放射計算値』は既にある
!  
! 短波計算の場合
!
!     IF ( ICNSW .EQ. 1 ) THEN
        CALL SIB0109_RUN_ALBEDO
     I        ( MON , JL , 
     O          AVISB , ANIRB , AVISD , ANIRD )
!
! ★★★★ 本来はここに放射計算 ★★★★
!
! 短波放射計算結果を SiB に与え、必要な情報を計算してセーブ
!
        CALL SIB0109_RUN_SR_CALC (  
     I   MON         , JL             , 
     I   RVISB(1,JL) , RVISD(1,JL)    , RNIRB(1,JL)   , RNIRD(1,JL) ) 
!
!     ELSE
!
! 短波放射関係の必要な情報をリストア （短波計算無しの場合）
!
        CALL SIB0109_RUN_SR_RESTORE ( JL )
!     ENDIF
!
! 長波放射の値を SiB に格納
!
      CALL SIB0109_RUN_LR_SET ( DLWB(1,JL) )
!
! 陸面計算
!
      CALL SIB0109_RUN_SIBMAIN (
     I   MON  , JL   , ID_NEXT , ID_NOW     )
!
      CALL SIB0109_RUN_LOOPEND (
     I   JL          , 
     O   TMTX2L_PHY  , QMTX2L_PHY  , UMTX2L_PHY  , VMTX2L_PHY   ,   
     O   RAD_LONG_SIB_2_ATM_PHY ) 
!
 7000 CONTINUE
!
!  ================
!  >> モニタ出力 <<
!  ================
!
#ifdef MONYOS
      IF ( ICNMNTMON.EQ.1 ) THEN
        CALL MONIT_OUT_2
     1    ( IMONIT_2FL     , IMONIT_2FL8   , IMONIT_SFL    ,
     1      IDATE , MODEL  , RESL  , EXPR  , CINF  ,
     1      KTP   , TOTMON , FSECP , IMASK , GW )
        CALL MONIT_CLEAR  ( TOTMON )
      ENDIF
#endif
!
      CALL SIB0109_RUN_STEPEND (
!    I    ID_NEXT , ID_NOW , FSECP , KTP )
     I    ID_NEXT , ID_NOW )
!
      LFIRST = .FALSE.  
!
      IF ( ICNLAST.EQ.0 ) GOTO 8000 
      !
      CALL MESSAGE_FLUSH
!
      write(6,*) 'normally end' 
!#endif
!
      STOP
      END
! =================================================================
      SUBROUTINE DELT_CHECK (
     I    DELT_ATM  , DELT_CNP  , DELT_SNOW , DELT_SOIL )
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1    JCN_IWL_SKIP
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN)    :: DELT_ATM
      REAL(8),INTENT(IN)    :: DELT_CNP
      REAL(8),INTENT(IN)    :: DELT_SOIL
      REAL(8),INTENT(IN)    :: DELT_SNOW
!
      INTEGER  :: IS
!
!  オイラーで組み込むときは、DELT_ATM = DELT_CNP 
!
      IF ( DELT_ATM .NE. DELT_CNP ) THEN
        WRITE(6,*) 'DELT_check : WARNING DELT_ATM NE DELT_CNP ' ,
     1             DELT_ATM , DELT_CNP
      ENDIF
!
!  DELT_CNP は DELT_ATM の約数
!
      IS = ( DELT_ATM + 0.1 ) / DELT_CNP
      IF ( ABS( IS * DELT_CNP - DELT_ATM ) .GT. 1.D-2 ) THEN
        WRITE(6,*) 'DELT_check ERROR DELT_ATM DELT_CNP ' ,
     1             DELT_ATM , DELT_CNP
        STOP 999
      ENDIF
!
!  DELT_SNOW = DELT_SOIL
!
      IF ( DELT_SNOW .NE. DELT_SOIL ) THEN
        WRITE(6,*) 'DELT_check : ERROR DELT_SNOW DELT_SOIL ' ,
     1             DELT_SNOW , DELT_SOIL
        STOP 999
      ENDIF
!
!  DELT_SOIL は 3600 の約数
!
      IS = 3601. / DELT_SOIL
      IF ( ABS( IS * DELT_SOIL - 3600.) .GT. 1.D-2 ) THEN
        WRITE(6,*) 'DELT_check ERROR DELT_SOIL ' , DELT_SOIL
        STOP 999
      ENDIF
!
!  DELT_CNP は DELT_SOIL の約数
!
      IS = ( DELT_SOIL + 0.1 ) / DELT_CNP
      IF ( ABS( IS * DELT_CNP - DELT_SOIL ) .GT. 1.D-2 ) THEN
        WRITE(6,*) 'DELT_check ERROR DELT_CNP DELT_SOIL ' ,
     1             DELT_CNP , DELT_SOIL
        STOP 999
      ENDIF
! 
!  |JCN_IWL_SKIP| =1,2 の時は毎ステップ土壌計算なので、
!  DELT_SOIL = DELT_CNP 
!
      IF (
     1     (      JCN_IWL_SKIP.EQ.1 .OR. JCN_IWL_SKIP.EQ.-1
     1       .OR. JCN_IWL_SKIP.EQ.2 .OR. JCN_IWL_SKIP.EQ.-2 )
     1     .AND. ( IS .NE. 1 )                               ) THEN
        WRITE(6,*) 'DELT_check ERROR DELT_ATM DELT_SOIL SKIP ' ,
     1             DELT_ATM , DELT_SOIL , JCN_IWL_SKIP
        STOP 999
      ENDIF
!
      RETURN
      END SUBROUTINE DELT_CHECK 
!===========================
      subroutine mod_imask ( imask , ispt, idim , jdim )
!
      use sibveg, only : 
     1   ITYP ,         ! 81
     1   ITYO_veg,      ! 12       
     1   ITYO ,         ! 25 
     1   ITYP_veg       ! 40

!
      implicit none
!
      integer,intent(in)    :: ispt
      integer,intent(in)    :: idim
      integer,intent(in)    :: jdim
      integer,intent(inout) :: imask(ispt*idim,jdim)
!
      INTEGER,PARAMETER :: ITYP_GLACIER = ITYP 
      INTEGER,PARAMETER :: ITYO_GLACIER = ITYO 
!
      integer :: i,j
!
      if ( ispt.ne.1 ) then
        write(6,*) 'mod imask error : not support ispt ne 1 ispt=',ispt 
        stop 999
      endif
!
      do j=1,jdim 
      do i=1,idim 
        if     ( imask(i,j).eq.13 ) then
          imask(i,j) = ityp_glacier  
        elseif ( j.gt.jdim/2 .and. imask(i,j).ne.0 ) then
          imask(i,j) = imask(i,j) + ityp_veg 
        endif
      enddo
      enddo
!
      return
      end
      



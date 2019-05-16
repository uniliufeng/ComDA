      MODULE COM_RUNCONF_SIB0109

!! 実行ジョブ設定情報 ( ジョブ内で一定 )

!! 下の設定値はdefaultの値。
!! set_com_runconfで、いくつかの値が置き換わる。

! run configuration information module file
!
!
!--- FLG8       新  SiB 関連
!  
!  個々のパッケージの選択   : 1 newsib 0 l3sib_mod
!
!     MJ98 では選択可能 ( それ以外では 1 のみ? ) 
!
      INTEGER,SAVE :: JCN_SIB_ALBEDO  = 1  
      INTEGER,SAVE :: JCN_SIB_RAD2SIB = 1  
      INTEGER,SAVE :: JCN_SIB_PRESIB  = 1  
      INTEGER,SAVE :: JCN_SIB_SKIN    = 1  
      INTEGER,SAVE :: JCN_SIB_SNOW    = 1  
      INTEGER,SAVE :: JCN_SIB_SOIL    = 1  
!
!
      INTEGER,SAVE ::
!    1    JCN_SIB_MONIT_DETAIL   = 1  ,    ! 詳細モニタ on なら 1  
     1    JCN_WAY_OF_INTEGRATION = 2  ,    ! canopy 積分 2 overimp, 1 balance
     1    JCN_BAD_WTR_ADD = 1         
!
! 親モデルの指定。 インターフェースで与えること。
!
      INTEGER,SAVE :: JCN_SIB_UPPER_MODEL = -1
!        オフライン  -1   大気比熱ゼロ  
!                    -2   大気比熱有限
!        mj98         1 
!        GSM          2 
!        RSM          3  
!
      INTEGER,SAVE ::                      ! on-line             0
     1    JCN_SNOW_HEAT_EXPLICIT = 0   ,   ! 雪層熱を explicit なら 1 
     1    JCN_SOIL_HEAT_EXPLICIT = 0   ,   ! 土壌熱を explicit なら 1 
     1    JCN_SOIL_WATER_EXPLICIT = 0  ,   ! 土壌水を explicit なら 1 
     1    JCN_SIB_ISLSCP = 0           ,   ! 強制データ種別
!                                          ! online -1 
!                                          ! TEST1DIM 0, 
!                                          ! ISLSCP 1, ERA 2, 
!                                          ! Forcing 3 
     1    JCN_SIB_ISLSCP_CYCLE1987 = 0  ,  ! ISLSCP で 1987 年サイクルラン
     1    JCN_SIB_ATM_CYCLE = 0         ,  ! ERA で年を代入するとその年をcycle
!                                          !  900 以上なら月末
     1    JCN_IWL_SKIP    = 3           ,  ! 
     1    JCN_RRUNOFF_B   = 0           ,  ! 土壌三層下からの流入無し 0 あり 1 
     1    JCN_CAP_CNP_SL  = 1           ,  ! フラックス評価に使う比熱が
!                                          ! 0 小さい、1 実際の値
     1    JCN_ITR_MAX     = 2           ,  ! イテレーションmax ... 仮積分なら2
     1    JCN_INT_MNT     = 999            ! モニタ出力間隔
!
      REAL(8),SAVE :: 
     1    RCN_DELT_ATM     =  3600.D0 ,    
     1    RCN_DELT_CNP     =  3600.D0 ,  
!    1    RCN_DELT_SOIL    = 1800.D0 ,  
!    1    RCN_DELT_SNOW    = 1800.D0 , 
     1    RCN_DELT_SOIL    =  3600.D0 ,  
     1    RCN_DELT_SNOW    =  3600.D0 , 
     1    RCN_FACT_SS1     = 1.D0    ,
     1    RCN_DELT_SR      = 3600.D0     
     
!            ! canopy_* で SS1 の積分時に、かけるファクタ
!            ! ゼロならそもそも変化無し。 
!
! 20050803
!      CHARACTER(40),SAVE   :: CCN_FILENAME_VEG = '../input/GSM9603.VEG' 
	CHARACTER(40),SAVE   :: CCN_FILENAME_VEG = 'input/GSM9603.VEG' 
      !
      ! デバッグ用に各種定数を動的に変更するため
      !    LINUX では ccn_trial の初期値は "" ではなく " " でないと警告が出た
      ! 
      CHARACTER(8), SAVE :: CCN_TRIAL(4) = (/" ", " ", " ", " "/)
      REAL(8)     , SAVE :: RCN_TRIAL(4) = (/0., 0., 0., 0./)
!
      NAMELIST / SIB0109_COM_RUNCONF_SIB0109/       
!    1   JCN_SIB_MONIT_DETAIL  , 
!!   1   JCN_SIB_UPPER_MODEL   ,   
     1   JCN_SNOW_HEAT_EXPLICIT     ,   
     1   JCN_SOIL_HEAT_EXPLICIT     ,   
     1   JCN_SOIL_WATER_EXPLICIT     ,   
     1   JCN_SIB_ISLSCP        ,   
     1   JCN_SIB_ISLSCP_CYCLE1987     ,   
     1   JCN_SIB_ATM_CYCLE     ,   
     1   JCN_IWL_SKIP          , 
     1   JCN_RRUNOFF_B         , 
     1   JCN_CAP_CNP_SL        , 
     1   RCN_DELT_ATM          , 
     1   RCN_DELT_CNP          , 
     1   RCN_DELT_SOIL         , 
     1   RCN_DELT_SNOW         , 
     1   RCN_DELT_SR           , 
     1   RCN_FACT_SS1          , 
     1   JCN_ITR_MAX           ,
     1   CCN_FILENAME_VEG, CCN_TRIAL, RCN_TRIAL
!
      CONTAINS
!======================================================================
        SUBROUTINE COM_RUNCONF_SIB0109_INI
!
        use com_stdinout_unit , only : iunit_in , iunit_out
!
        implicit none 

         !--- ネームリスト読み込みによる変更       
          write(6,*) 'com_runconf_sib0007_ini' 
          REWIND(iunit_in) 
          READ  (iunit_in,SIB0109_COM_RUNCONF_SIB0109)  ! read parameter values
													  !	from stdin
!
          IF      ( JCN_WAY_OF_INTEGRATION .EQ. 1 ) THEN   ! イテレーション
           JCN_ITR_MAX     = 5    ! 
           JCN_CAP_CNP_SL  = 0    ! フラックス評価に使う比熱がゼロ
          ELSE IF ( JCN_WAY_OF_INTEGRATION .EQ. 2 ) THEN   ! 仮積分
           JCN_ITR_MAX     = 2    ! 
           JCN_CAP_CNP_SL  = 1    ! フラックス評価に使う比熱が普通の値
          ENDIF
!
          WRITE(iunit_out,SIB0109_COM_RUNCONF_SIB0109)
!
        END SUBROUTINE COM_RUNCONF_SIB0109_INI
      END MODULE COM_RUNCONF_SIB0109





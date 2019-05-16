! SiB 用定数ファイル
! 99/05/xx 保坂征宏    まだまだ変更の必要あり
! 99/07/15 大泉三津夫
! 99/08/04 保坂征宏    _c つきに変更途中
!
      MODULE SIBCON
!
!
      IMPLICIT NONE
!
      PRIVATE  CAP_SKIN_c  
!
      REAL(8),SAVE ::  PIE_c  = 3.141592653589793D0
      REAL(8),SAVE ::  COSZ_MIN_c = 0.01746D0 
      REAL(8),SAVE ::  GRAV_c = 9.80665D0
      REAL(8),SAVE ::  CP_AIR_c = 1004.6D0    ! 乾燥空気の定圧比熱(J/kg/K)
      REAL(8),SAVE ::  AKAPPA_c               !  R_DRY_c / CP_AIR_c
!
      REAL(8),SAVE ::  STEFAN_c = 5.67D-8  ! Stefan-Boltzmann 定数(W/m^2/K^4)
      REAL(8),SAVE ::  STB4_c              ! 4*STEFAN_c
      REAL(8),SAVE ::  STB8_c              ! 8*STEFAN_c
!
      REAL(8),SAVE ::  TMP_FREZ_c= 273.15D0 ! 凍結温度
!
      REAL(8),SAVE ::  EPSFAC_c             ! R_DRY_c / R_VAP_c 
!
      REAL(8),SAVE ::  CPSY_c               ! CP_AIR_c/(HL_WTR_c*EPSFAC_c) 

      REAL(8),SAVE ::  VKRMN_c      = 0.4D0   ! カルマン定数
      REAL(8),SAVE ::  VKRMN_INV_c  = 2.5D0   ! 同逆数
!
! 密度
!
!     REAL(8),SAVE ::  RHO_SOIL_c = 2650.0D0 !土粒子（粘土鉱物）の密度(kg/m^3) 
!                                           CAMPBELLのテキストより
!    積雪における氷の密度のようなもの。土壌のバルク密度（積雪の濡れ雪密度
!    に相当）にはこの密度の他に、土壌水分、空気、土壌氷の密度と存在比がからむ
!
      REAL(8),SAVE ::  RHO_GLA_c = 917.0D0     ! 氷の密度(kg/m^3)
      REAL(8),SAVE ::  RHO_WTR_c = 1000.0D0    ! 水の密度(kg/m^3)
      REAL(8),SAVE ::  RHO_WTR_INV_c 
!
      REAL(8),SAVE ::  RHO_SNW_c = 200.D0      ! ？
      REAL(8),SAVE ::  RHO_SNW_INV_c 
!
!
! 比熱
!
!   JCN_CAP_CNP_SL = 0 の場合の雪面スキンの熱容量（flux 評価時に使用） 
      REAL(8),SAVE ::  CAP_CNPX_c = 0.1 
!
      REAL(8),SAVE ::  CAP_LAI_c = 840.D0    ! 葉単位面積当たりの比熱：
!                                         CLAI (J/K/m^2)・・0.2mmの水の厚さのCW
      REAL(8),SAVE ::  CAP_WTR_c = 4218.0D0  ! 水比熱：CW (J/K/kg)
      REAL(8),SAVE ::  CAP_GLA_c = 2106.0D0  ! 氷比熱：CICE (J/K/kg)
      REAL(8),SAVE ::  CAP_SOIL_BY_VOL_c= 2.39D6
                                        ! 単位深単位断面積当粘土鉱物比熱
!            CSOIL (J/K/m^3)   CAMPBELLのテキストより
!              これは水や空隙を含まない鉱物である。 
!              空気の比熱を無視している。
!              これに (１－空隙率）をかけると乾燥土壌の比熱になる。
!
      REAL(8),SAVE ::  CAP_GLA_c_INV          ! CICEの逆数
!
! 融解・蒸発・昇華潜熱
!
      REAL(8),SAVE ::  SNOMEL_c =  333624.2D3    ! 氷融解熱 (J/kg)*1000
                                                 ! いずれ不要。
!                                                ! HL_GLA_c-HL_WTR_c  
! 
!                            -> つまり (J/融解時に高さ1m*断面1m^2の水になる氷)
!                    -> 高さ1m*断面1m^2の氷を融かすには、SNOMEL*RHO_GLA_c/1000
      REAL(8),SAVE ::  HL_W2G_c            ! 氷融解熱 (J/kg)*1000
      REAL(8),SAVE ::  HL_W2G_INV_c        ! 氷融解熱 (J/kg)*1000
!
      REAL(8),SAVE ::  HL_WTR_c     = 2.507D6    !水の０℃での蒸発潜熱 (J/kg)
      REAL(8),SAVE ::  HL_GLA_c     = 2.8345D6   !氷の０℃での昇華潜熱 (J/kg)
      REAL(8),SAVE ::  HL_WTR_INV_c  
      REAL(8),SAVE ::  HL_GLA_INV_c 
!
      REAL(8),SAVE ::  H2O_SAT_CNP_DEP_c = 1.D-4  ! 0.1mm 
      REAL(8),SAVE ::  H2O_SAT_GRS_DEP_c = 1.D-4  ! 0.1mm 
!     REAL(8),SAVE ::  H2O_SAT_CNP_DEP_c = 2.D-3  ! 2mm       ppp check
!     REAL(8),SAVE ::  H2O_SAT_GRS_DEP_c = 2.D-3  ! 2mm       ppp check 
!
! 熱伝導率
!
!
      REAL(8),SAVE ::  HDF_SOIL_c = 2.92      ! 粘土鉱物熱伝導率 (W/m/K)
!                                          CAMPBELLのテキストより
!
      REAL(8),SAVE ::  HDF_WTR_c  = 0.561      ! 水熱伝導率   (W/m/K) pppqqq
!     REAL(8),SAVE ::  HDF_WTR_c  = 3.         ! 水熱伝導率   (W/m/K) ppp
      REAL(8),SAVE ::  HDF_GLA_c  = 2.24       ! 氷熱伝導率   (W/m/K)
      REAL(8),SAVE ::  HDF_AIR_c  = 0.0241     ! 空気熱伝導率 (W/m/K)
!
!     REAL(8),SAVE ::  RAMSOL = 2.2       ! 以前の土壌熱伝導率。 使わない 
!
! 気体定数
!
      REAL(8),SAVE ::  R_VAP_c    = 461.0D0    ! 水蒸気気体定数   (J/K/kg)
      REAL(8),SAVE ::  R_DRY_c    = 287.04D0    ! 乾燥空気気体定数 (J/K/kg)
      REAL(8),SAVE ::  R_DRY_INV_c             ! 乾燥空気気体定数 (J/K/kg)
      REAL(8),SAVE ::  R_BY_G_c                ! R_DRY/GRAV_c 
!
! 降水遮断
!
      REAL(8),SAVE :: B_PR_COEF_c = 5.       ! 降水の分布係数 (はてな？） 
      REAL(8),SAVE :: A_PC_COEF_c 
      REAL(8),SAVE :: A_PL_COEF_c 
      REAL(8),SAVE :: C_PC_COEF_c 
      REAL(8),SAVE :: C_PL_COEF_c 
!
!     REAL(8),SAVE :: A_PR_COEF_S_c 
!     REAL(8),SAVE :: C_PR_COEF_S_c 
      REAL(8),SAVE :: FR_PREITC_c = 0.5 ! インターセプト前の値を混ぜる割合 0-1
!                                       ! ダミー
!
      REAL(8),PARAMETER  :: CAP_SKIN_c = 2.5D4 
!
! ==========
! >> 積雪 <<
! ==========
!
!   被覆率
!
      REAL(8),SAVE :: CVR_SNW_INV_MAX_c  = 1.D30   ! 雪の被覆率逆数（ダミー）   
      REAL(8),SAVE :: CVR_VEG_INV_MAX_c  = 1.D30   ! 
      REAL(8),SAVE :: CVR_VEG_MIN_c      = 1.D-30  ! 
!
!   雪の物性
!
      REAL(8),SAVE :: FR_SNOWWTR_MAX_c = 0.1 ! GLA 中に保有可能な液水の割合上限
!                                            ! WTR/GLA < FR_SNOWWTR_MAX 
!
      REAL(8),SAVE :: CAP_SNSK_c  = CAP_SKIN_c   ! 雪スキン熱容量=CAP_SKIN_c 
!!      REAL(8),SAVE :: CAP_SNSK_c  = 2.5D5
!
!      JCN_CAP_CNP_SL = 0 の場合の雪面スキンの熱容量（flux 評価時に使用） 
      REAL(8),SAVE :: CAP_SNSKX_c =  0.1    
!
!   年齢 ( 今はほとんどダミー ) 
!
      REAL(8),SAVE :: AGE_SNOW_MAX1_c   = 1*24*3600.D0     ! 新雪層  年齢上限 
      REAL(8),SAVE :: AGE_SNOW_MAX1_M_c = 1*24*3600.D0 - 1 
!
      REAL(8),SAVE :: AGE_SNOW_MAX2_c   = 10*24*3600.D0    ! "新雪層' 年齢上限"
      REAL(8),SAVE :: AGE_SNOW_MAX2_M_c = 10*24*3600.D0 -1 
!
!     雪アルベドのagingの時定数を別冊46号P60に従い変更 8時間→4日
!       ただしこれは ITUNE_ALBEDO_SNOW=1 のときのみ有効
!!#if    ( SIB0109_OFFLINE_JMA == 1 || SIB0109_GSM_JMA  == 1 )
      REAL(8),SAVE :: AGE_SNOW_ETIME_c     = 4*86400.D0  ! 雪アルベド etime
      REAL(8),SAVE :: AGE_SNOW_ETIME_INV_c = 1.0D0/(4.*24.*60.*60.D0)
!!#else
!!      REAL(8),SAVE :: AGE_SNOW_ETIME_c     = 8*3600.D0   ! 雪アルベド etime 
!!      REAL(8),SAVE :: AGE_SNOW_ETIME_INV_c = 3.472*1.D-5    
!!#endif
!
!   年齢クリア ( 今はダミー ) 
!
      REAL(8),SAVE :: SNOWFALL_MAKE_NEW_c  = 0.25 ! 新雪層作成の降雪量閾値  
!                                     ! 0.5kg/m**2/step ～ 5cm/step 
!                                              
      REAL(8),SAVE :: SNOWFALL_RESET_NEW_c = 0.05 
                                      ! 新雪層の年齢リセットの降雪量閾値  
!
!  itune_albedo_snow=3 or 4
      REAL(8),SAVE :: RAIN_1HOUR_FRESH_c = 0.02    ! kg=mm/1hour ~ 0.5mm/day
!
!  01 と 1 のさかいめ
!
      REAL(8),SAVE :: H2O_SNOW_01_1_c  = 4.D0  ! これ以下なら 01 にする。 
      REAL(8),SAVE :: H2O_SNOW_01_1_INV_c   
!
!  少雪バケツ行きかどうかのさかいめ ( 気象研と数値とで異なる ) 
!
! 20050517
!      REAL(8),SAVE :: H2O_SNOW_MIN_c
#if   ( SIB0109_OFFLINE_MRI == 1 ) 
      REAL(8),SAVE :: H2O_SNOW_MIN_c   = 4.D0  ! 混在無しの場合
#elif ( SIB0109_OFFLINE_ALMA == 1 ) 
      REAL(8),SAVE :: H2O_SNOW_MIN_c   = 4.D0  ! 混在無しの場合
#elif ( SIB0109_MJ98_MRI    == 1 ) 
      REAL(8),SAVE :: H2O_SNOW_MIN_c   = 4.D0  ! 混在無しの場合
#elif ( SIB0109_OFFLINE_JMA == 1 ) 
      REAL(8),SAVE :: H2O_SNOW_MIN_c   = 1.D-3 ! 混在型デフォルト
#elif ( SIB0109_GSM_JMA     == 1 ) 
      REAL(8),SAVE :: H2O_SNOW_MIN_c   = 1.D-3 ! 混在型デフォルト
#endif
!
      REAL(8),SAVE :: ENG_MELT_MAX_c 

!
! 層の切り方
!
      REAL(8),SAVE :: H2O_SNOW_1_2a_c  = 40.  ! これ以上なら 2 層。1,2 層は折半
      REAL(8),SAVE :: H2O_SNOW_2a_2b_c = 80.  ! これ以上なら 2 層で 1 層は定厚
      REAL(8),SAVE :: H2O_SNOW_2b_3_c  = 200. ! これ以上なら 3 層。
!
      REAL(8),SAVE :: H2O_SNOW1_STD_c         ! H2O_OLDSNOW_2a_2b /2
      REAL(8),SAVE :: H2O_SNOW3_STD_c         ! H2O_OLDSNOW_2a_2b /2
!   
      REAL(8),SAVE :: H2O_SNOW1_MIN_c  = 10.  ! 2 層以上のときの 1 層の最小値
      REAL(8),SAVE :: H2O_SNOW1_MAX_c  = 50.  ! 
      REAL(8),SAVE :: H2O_SNOW2_MIN_c  = 15.  
      REAL(8),SAVE :: H2O_SNOW2_MAX_c  = 1.D10 
      REAL(8),SAVE :: H2O_SNOW3_MIN_c  = 20.  
      REAL(8),SAVE :: H2O_SNOW3_MAX_c  = 50.  
!
! ==========
! >> 土壌 <<
! ==========
!
!    JCN_CAP_CNP_SL = 1 の場合の土壌スキンの熱容量（flux 評価時に使用） 
      REAL(8),SAVE :: CAP_GRSKX_c =  0.1  
!
      REAL(8),SAVE :: CAP_GRSK_c = CAP_SKIN_c 
                                   ! 土壌スキンの熱容量 = cap_skin_c 
!
      REAL(8),SAVE :: FR_SOIL_H2O_MIN_c = 0.05
      REAL(8),SAVE :: FR_SOIL_GLA_MAX_c = 0.95
!
! 土壌水分の拡散のタイプ ( ダミー )
!   
!     INTEGER   :: JCN_COEF_WFLUX_SOIL_c  = 0   
!               ! 0  MILLY AND EAGLESON (1982)?   デフォルト
!               ! 1  拡散ではよくあるタイプ
!
      REAL(8),SAVE :: DEP_SOIL4_c     = 10.D0    ! 10m 
      REAL(8),SAVE :: DEP_SOIL4_INV_c 
      REAL(8),SAVE :: CAP_SOIL4_c   
      REAL(8),SAVE :: CAP_SOIL4_INV_c   
      integer,save :: itune_soilfluxbottom_c = 1 ! 土壌底面熱フラックスあり 
!
! 氷床裸地アルベド(気候値ファイルではVIS0.1 NIR=0.15と小さすぎ？)
!#if    ( SIB0109_OFFLINE_JMA == 1 || SIB0109_GSM_JMA  == 1 )
      REAL(8),SAVE :: ALBEDO_ICESHEET_SOIL_VIS_c = 0.80D0
      REAL(8),SAVE :: ALBEDO_ICESHEET_SOIL_NIR_c = 0.60D0
!#endif
!
! ========================================
! >> 時間刻ファクタ（フラックス評価時） <<
! ========================================
!     delt_fac_* = 1 でインプリシット、 > 1 でオーバーインプリシット
!
! 20050517
!	REAL(8),SAVE   :: DELT_FAC_ATM_c
!	REAL(8),SAVE   :: DELT_FAC_CNP_c
#if   ( SIB0109_OFFLINE_MRI == 1 ) 
      REAL(8),SAVE   :: DELT_FAC_ATM_c = 1.5D0
      REAL(8),SAVE   :: DELT_FAC_CNP_c = 2.0D0
#elif ( SIB0109_OFFLINE_ALMA == 1 ) 
      REAL(8),SAVE   :: DELT_FAC_ATM_c = 1.5D0
      REAL(8),SAVE   :: DELT_FAC_CNP_c = 2.0D0
#elif ( SIB0109_MJ98_MRI == 1 ) 
      REAL(8),SAVE   :: DELT_FAC_ATM_c = 1.5D0
      REAL(8),SAVE   :: DELT_FAC_CNP_c = 2.0D0
#elif ( SIB0109_OFFLINE_JMA == 1 ) 
      REAL(8),SAVE   :: DELT_FAC_ATM_c = 2.0D0
      REAL(8),SAVE   :: DELT_FAC_CNP_c = 2.0D0
#elif ( SIB0109_GSM_JMA  == 1 ) 
      REAL(8),SAVE   :: DELT_FAC_ATM_c = 2.0D0
      REAL(8),SAVE   :: DELT_FAC_CNP_c = 2.0D0
#endif

      REAL(8),SAVE   :: DELT_FAC_SNOW_c  = 1.   
      REAL(8),SAVE   :: DELT_FAC_SOIL_c  = 1.  
      REAL(8),SAVE   :: DELT_FAC_SOIL_INV_c
!
! ====================
! >> 4 倍精度の変数 <<
! ====================
!
      INTEGER, PARAMETER:: QUADRUPLE = MAX(SELECTED_REAL_KIND(P=20),
     1  KIND(0.0d0))
!
      REAL(QUADRUPLE)::    HL_W2G_c_16       = 0.3275e6_QUADRUPLE 
      REAL(QUADRUPLE)::    HL_W2G_INV_c_16         
      REAL(QUADRUPLE)::    HL_GLA_c_16      = 2.8345e6_QUADRUPLE
      REAL(QUADRUPLE)::    HL_GLA_c_INV_16  
      REAL(QUADRUPLE)::    TMP_FREZ_16       = 273.15e0_QUADRUPLE
!
! ==================
! >> チューニング <<
! ==================
!
! 雪の加齢方法
!   加齢は、Sib0109/Snow/snow_wflux_pre.F で行なう。
!     10 日が上限で足していく。
!   年齢クリア
!     Sib0109/Snow/snow_age_clear.F 参照
!
!       0 ... 加齢のみ(クリアしない)
!       1 ... 密度が RHO_THR1 より小さければ 0 
!       2 ... 
!
      INTEGER,SAVE  ::  ITUNE_SNOW_AGE = 0  
!
!
!   雪年齢クリア条件を第一層の雪密度で決める場合 ( ITUNE_SNOW_AGE=1 )
!   に有効になる、新雪の敷居値
!
      REAL(8),SAVE ::  RHO_SNOW_AGECLEAR_THR_c  = 78.D0 
      REAL(8),SAVE ::  RHO_SNOW_AGECLEAR_THR_INV_c 
!
! 雪のアルベド計算方法 
!   詳しくは Sib0109/Albedo/albedo_snow.F 参照
!     0 ... 温度のみに依存
!     1 ... 年齢に依存
!     2 ... 
!     3 ... Aoki
!     4 ... 3 + Glacier special treatment 
!
      INTEGER,SAVE  ::  ITUNE_ALBEDO_SNOW = 3  
!
!     ITUNE_ALBEDO_SNOW = 0 の時 ( 温度のみに依存の時 ) に意味あり
!
        REAL(8),SAVE :: ALB_SNOW_VIS_c       = 0.8D0
        REAL(8),SAVE :: ALB_SNOW_NIR_c       = 0.4D0 
        REAL(8),SAVE :: ALB_SNOW_FACT_TEMP_c = 0.6D0   
!
!     ITUNE_ALBEDO_SNOW = 1 ( 年齢に依存の時 ), 3, 4 に意味あり
!                                                   3,4は(_FOR_cのみ使用)
!      
        REAL(8),SAVE :: ALB_SNOW_VIS_INI_c = 0.8   ! 可視アルベド初期値
        REAL(8),SAVE :: ALB_SNOW_NIR_INI_c = 0.6   ! 近赤アルベド初期値
!       REAL(8),SAVE :: ALB_SNOW_VIS_FOR_c = 0.5   ! 可視アルベド向かう値
        REAL(8),SAVE :: ALB_SNOW_VIS_FOR_c = 0.6   ! 可視アルベド向かう値
        REAL(8),SAVE :: ALB_SNOW_NIR_FOR_c = 0.4   ! 近赤アルベド向かう値
!
! 雪圧密予報に使う関数形での敷居値。
!
      REAL(8),SAVE ::  RHO_SNOW_INTG_THR_c  = 78.D0 
      REAL(8),SAVE ::  RHO_SNOW_INTG2_THR_c = 400.D0
!
! 雪スキンと雪の間の拡散係数ファクタ ( presib_conduc_all.F )
      REAL(8),SAVE ::  FACT_HDF_SNSK_c = 1.D0  
      REAL(8),SAVE ::  FACT_HDF_SNOW_c = 1.D0  
!
! 雪スキンでの短波透過ファクタ ( 0 透過なし、1 透過 ) ( rad_sum_all.F )
!
!!(NPD)
!!#if    ( SIB0109_OFFLINE_JMA == 1 || SIB0109_GSM_JMA  == 1 )
      REAL(8),SAVE ::  FACT_TRANS_SR_SNOW_SKIN_c = 1.D0  
!(MRI)
!!#else
!!      REAL(8),SAVE ::  FACT_TRANS_SR_SNOW_SKIN_c = 0.D0  
!!#endif
!
! RSiB での雪の比熱評価用積雪水当量
!
      REAL(8),SAVE :: WTR_EQUIV_SNOW_c = 0.05  ! 5cm 
!
! 底面融解あり無し ( snow_bottom_melt.F )
      INTEGER,SAVE  ::  ITUNE_SNOW_BOTTOM_MELT = 1   
!
! TMP_SKIN の積分時に、下の温度の陰的補正の係数。
      REAL(8),SAVE  :: FACT_HDF_GRSK_c = 1.D0     ! grsk と soil1 
      REAL(8),SAVE ::  FACT_HDF_SOIL_c = 1.D0     ! soil* 間 
!
! skin と 第一層の係数の陰的補正の修正( canopy_main.F )
      REAL(8),SAVE  :: FACT_TSKIN_TOP_c = 1.D0    ! 
!
      INTEGER,SAVE :: JCN_RUNOFF_GRAV_c =  1
      REAL(8),SAVE :: COEF_SOIL_WFLUX_MIN_C = 0. 
      REAL(8),SAVE :: POT_SOIL_MIN_c = - 1.D5 
!
      real(8),save :: fact_precipi_c = 1.D0  
!
      CONTAINS 
!======================================================
      SUBROUTINE SIBCON_INI 
!
      NAMELIST / SIB0109_SIBCON /  
     1     H2O_SNOW_MIN_c , fact_precipi_c , DEP_SOIL4_c 
     1   , CAP_GRSK_c      , CAP_SNSK_c     
     1   , FR_SOIL_H2O_MIN_c , COEF_SOIL_WFLUX_MIN_c,POT_SOIL_MIN_c
     1   , DELT_FAC_ATM_c  , DELT_FAC_CNP_c 
     1   , DELT_FAC_SNOW_c , DELT_FAC_SOIL_c 
     1   , RHO_SNOW_INTG_THR_c  
     1   , RHO_SNOW_AGECLEAR_THR_c 
     1   , ITUNE_SNOW_AGE  
     1   , ITUNE_ALBEDO_SNOW 
     1   , ITUNE_SNOW_BOTTOM_MELT 
     1   , ALB_SNOW_VIS_INI_c , ALB_SNOW_NIR_INI_c 
     1   , ALB_SNOW_VIS_FOR_c , ALB_SNOW_NIR_FOR_c 
     1   , ALB_SNOW_VIS_c     , ALB_SNOW_NIR_c 
     1   , ALB_SNOW_FACT_TEMP_c 
     1   , FACT_TRANS_SR_SNOW_SKIN_c
     1   , FACT_TSKIN_TOP_c 
     1   , FACT_HDF_SNSK_c 
     1   , FACT_HDF_SNOW_c 
     1   , FACT_HDF_GRSK_c
     1   , FACT_HDF_SOIL_c
     1   , WTR_EQUIV_SNOW_c 
     1  , ITUNE_SOILFLUXBOTTOM_c, RAIN_1HOUR_FRESH_c, JCN_RUNOFF_GRAV_c
!
      READ (5,SIB0109_SIBCON)      
      WRITE(6,SIB0109_SIBCON)      
!
!
      IF ( H2O_SNOW_MIN_c .GT. H2O_SNOW_01_1_c  ) THEN
        write(6,*) 'sibcon_ini error : snow_min ' , 
     1             ' H2O_SNOW_MIN_c > H2O_SNOW_01_1_c '  , 
     1               H2O_SNOW_MIN_c, H2O_SNOW_01_1_c 
        stop 999
      ENDIF
!
!  ====================
!  >> 逆数などの設定 <<
!  ====================
!
      HL_WTR_INV_c = 1/HL_WTR_c  
      HL_GLA_INV_c = 1/HL_GLA_c  
      HL_W2G_c     = HL_GLA_c - HL_WTR_c 
      HL_W2G_INV_c = 1/HL_W2G_c
!
      CAP_GLA_c_INV  = 1/CAP_GLA_c
!
      STB4_c         = STEFAN_c * 4
      STB8_c         = STEFAN_c * 8
!
      EPSFAC_c     = R_DRY_c / R_VAP_c                   ! 約 0.622 
! 
      CPSY_c       = CP_AIR_c / ( HL_WTR_c * EPSFAC_c ) 
      AKAPPA_c     = R_DRY_c / CP_AIR_c
!
      R_DRY_INV_c    = 1 / R_DRY_c
      RHO_WTR_INV_c  = 1 / RHO_WTR_c
!
      RHO_SNW_INV_c  = 1.D0 / RHO_SNW_c            ! 01/07/01
!
      RHO_SNOW_AGECLEAR_THR_INV_c = 1.D0 / RHO_SNOW_AGECLEAR_THR_c
!
      R_BY_G_c     = R_DRY_c/GRAV_c 
!
      DELT_FAC_SOIL_INV_c = 1 / DELT_FAC_SOIL_c
!
!
      H2O_SNOW_01_1_INV_c = 1/ H2O_SNOW_01_1_c
!
      H2O_SNOW1_STD_c = H2O_SNOW_2a_2b_c/2
      H2O_SNOW3_STD_c = H2O_SNOW_2a_2b_c/2
!
      HL_W2G_INV_c_16  = 1.0_QUADRUPLE /  HL_W2G_c_16   
      HL_GLA_c_INV_16  = 1.0_QUADRUPLE /  HL_GLA_c_16  
!
!
      DEP_SOIL4_INV_c = 1/DEP_SOIL4_c
      CAP_SOIL4_c     = 0.5 * DEP_SOIL4_c * CAP_SOIL_BY_VOL_c
!                       0.5 は、代表的空隙率
      CAP_SOIL4_INV_c = 1./CAP_SOIL4_c
!
! ------------------------------------------------
! > 降水の分布関数の係数 （チューニングパラメタ）<
! ------------------------------------------------
!    注意。佐藤(1989) とは値が異なる
!
        A_PC_COEF_c = 5.0D0                          ! 対流系降水 ピークは 5 
        C_PL_COEF_c = 0.99998D0                      ! 大規模凝結 一様
        A_PL_COEF_c = ( 1.D0 - C_PL_COEF_c ) * B_PR_COEF_c 
     1               / ( 1 - DEXP(-B_PR_COEF_c) )   
        C_PC_COEF_c =   1.D0 - A_PC_COEF_c   / B_PR_COEF_c  
     1                       * ( 1 - DEXP(-B_PR_COEF_c) ) 
!
#ifdef SIB_DEBUG
         write(6,*) 'sibcon_ini : rain distribution info' , 
     1             ' a_pc c_pc ' , A_PC_COEF_c , C_PC_COEF_c ,      
     1             ' a_pl c_pl ' , A_PL_COEF_c , C_PL_COEF_c ,
     1             ' integ c l ' , 
     1               A_PC_COEF_c / B_PR_COEF_c 
     1                      * ( 1.D0 - DEXP(- B_PR_COEF_c) )  
     1                + C_PC_COEF_c ,
     1               A_PL_COEF_c / B_PR_COEF_c 
     1                      * ( 1.D0 - DEXP(- B_PR_COEF_c) ) 
     1                + C_PL_COEF_c  
#endif
!
!  雪の場合
!
!       A_PR_COEF_S_c = 3.D0                         ! ピークは 3 倍
!       C_PR_COEF_S_c = 1.D0 - A_PR_COEF_S_c / B_PR_COEF_c  
!    1                       * ( 1 - EXP(-B_PR_COEF_c) ) 
!
!  ====================
!  >> 種々のチェック <<
!  ====================
!
! --------------------------------------
! > 降水分布のパラメタに関するチェック <
! --------------------------------------
!
      IF (        A_PL_COEF_c   .LT.0. 
     1       .OR. C_PC_COEF_c   .LT.0. 
!    1       .OR. C_PR_COEF_S_c .LT.0  
     1   ) THEN
          write(6,*) 'sibcon_ini error rain distribution A '  
     1                                   , A_PL_COEF_c 
     1                                   , C_PC_COEF_c  
!    1                                   , C_PR_COEF_S_c 
          stop 999
      ENDIF
!
      IF ( 
     1          ( ABS( A_PL_COEF_c / B_PR_COEF_c 
     1                      * ( 1.D0 - DEXP(- B_PR_COEF_c) ) 
     1            + C_PL_COEF_c - 1.D0 ) 
     1           .GT. 1.D-15 )
     1       .OR.  
     1          ( ABS( A_PC_COEF_c / B_PR_COEF_c 
     1                      * ( 1.D0 - DEXP(- B_PR_COEF_c) ) 
     1            + C_PC_COEF_c -1.D0 ) 
     1            .GT. 1.D-15 ) 
!    1       .OR.  
!    1          ( ABS( A_PR_COEF_S_c / B_PR_COEF_c 
!    1                      * ( 1 - EXP(-B_PR_COEF_c)) 
!    1            + C_PR_COEF_S_c -1.D0 ) .GT. 1.D-30 ) 
     1      ) THEN
          WRITE(6,*) 'sibcon_ini ERROR rain distribution B' 
     1     ,  A_PL_COEF_c / B_PR_COEF_c * ( 1 - DEXP(- B_PR_COEF_c) ) 
     1          + C_PL_COEF_c - 1 
     1     ,  A_PC_COEF_c / B_PR_COEF_c * ( 1 - DEXP(- B_PR_COEF_c) ) 
     1          + C_PC_COEF_c - 1 
!    1     ,  A_PR_COEF_S_c / B_PR_COEF_c * ( 1 - EXP(- B_PR_COEF_c) ) 
!    1          + C_PR_COEF_S_c - 1 
          stop 999
      ENDIF
!
! --------------------------------------------------
! > 雪の年齢をクリアしないのに、アルベドは年齢依存 <
! --------------------------------------------------
!
      IF ( ITUNE_ALBEDO_SNOW .EQ.1 .AND. ITUNE_SNOW_AGE.EQ.0 ) THEN
        WRITE(6,*) 'SIBCON_INI ERROR : '
        WRITE(6,*) '  SNOW AGE IS NOT CLEARED, ' 
        WRITE(6,*) '  BUT ALBEDO DEPEND ON THE SNOW AGE.'
        STOP 999
      ENDIF 
!
      ENG_MELT_MAX_c = H2O_SNOW1_MIN_c * 0.1D0 * HL_W2G_c 
!
      END SUBROUTINE SIBCON_INI 
!
      END MODULE SIBCON

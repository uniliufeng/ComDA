! 蒸発・昇華用の係数および、陽に評価した蒸発・昇華。雪無し域。
!
! 99/06/xx 保坂征宏
!          蒸散を蒸発から分離していない。
! 99/07/04 蒸散 TRNS_N は蒸発と別扱いに
! 99/07/17 雪あり・雪無しを統合
!          雪域については未チェック。
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/30 保坂征宏 フラックスを下向き正に変更
! 99/09/02 保坂征宏 
! 99/09/16 保坂征宏 バグ固定。EVAP_TOT と SUBL_TOT の配列
!
      SUBROUTINE CANOPY_EVAP (
     I  NSTR        , NEND        , e_air , 
     I  ITR         , ITR_MAX           , 
#ifdef MONYOS
     I  DELT_CNP    , INSFLAG     , 
#endif
     I  R_AIR_INV   , 
     I  R_AC_BAR_INV, R_ACX_BAR_INV , 
     I  R_AG_INV    , R_AGX_BAR_INV   , R_AB_INV   ,  
     I  E_REF       , E_TC_SAT       , E_TG_SAT     , 
     I  DEC_SAT_DTC , DEG_SAT_DTG , 
     I  CVR_VEG_ms  ,                                          ! 被覆率
     I  FAC_EVP     , HUM_REL_SKIN    ,                            ! 
     I  FAC_EVP_2   , 
     I  XFR_WTR_CNP  , XFR_WTR_GRS  , XFR_WTR_BS   ,   ! 水氷割合
     I  XFR_GLA_CNP  , XFR_GLA_GRS  , XFR_GLA_BS   ,   ! 水氷割合 
     I  FLAG_DRY_CNP      , FLAG_DRY_SKIN     ,
     I  FLAG_TR_ALL ,
     O  EVAP        , SUBL        , TRNS        , 
     O  DEVAP_DX   )
!
! 蒸発・昇華用の係数および、陽に評価した蒸発・昇華。雪無し域。
! 下草からの蒸散は含まない。
! 
! 注意
!
!   ・蒸発量は熱、水の両フラックスの計算に使う。
!     共通に使えるように、潜熱係数はかけない、水蒸気量の変化率などを出力する。
!
!   ・地面と雪面の蒸発等には、(1-Cs) と Cs とをかけていない。
!     それぞれの単位面積あたり、という扱いである。
!
!   ・氷は蒸散しない。
!
!   ・符号の向きは？
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ICG    ,
     1    IWG    ,
     1    IEST   ,
     1    ICGB   ,
     1    ICGR    
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS ,
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
#endif
!
      IMPLICIT NONE 
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
!
#ifdef MONYOS
      REAL(8),INTENT(IN)    :: DELT_CNP 
      INTEGER,INTENT(IN)    :: INSFLAG             ! 雪無し 0 あり  1
#endif     
! 
      REAL(8),INTENT(IN)    :: R_AIR_INV (IDX)   ! 1/r_air   (=1/ra)
      REAL(8),INTENT(IN)    :: R_AC_BAR_INV  (IDX)   ! 1/r_ac    (=1/rb)
!     REAL(8),INTENT(IN)    :: R_CC_INV  (IDX)   ! 1/r_cc    (=1/rc) 
      REAL(8),INTENT(IN)    :: R_ACX_BAR_INV (IDX)   ! 1/(r_ac+r_cst) 
      REAL(8),INTENT(IN)    :: R_AG_INV  (IDX)   ! 1/r_d     (=1/rd) 
!     REAL(8),INTENT(IN)    :: R_BB_INV  (IDX)   ! 1/R_BB    (=1/rg)   
      REAL(8),INTENT(IN)    :: R_AGX_BAR_INV (IDX)   ! 1/(r_ag+r_gst)
      REAL(8),INTENT(IN)    :: R_AB_INV (IDX)   ! 1/(r_ag+r_bs)
!
      REAL(8),INTENT(IN)    :: E_REF   (IDX)   ! 大気最下層比湿（参照比湿）
      REAL(8),INTENT(IN)    :: E_TC_SAT (IDX)   ! キャ温度飽和比湿
      REAL(8),INTENT(IN)    :: E_TG_SAT (IDX)   ! 地面温度飽和比湿
!
      REAL(8),INTENT(IN)    :: DEC_SAT_DTC  (IDX)  ! キャ温度飽和比湿 微分
      REAL(8),INTENT(IN)    :: DEG_SAT_DTG  (IDX)  ! 地面温度飽和比湿 微分
!
      REAL(8),INTENT(IN)    :: CVR_VEG_ms (IDX,ICG)    ! 植生被覆率 雪なら 0
!
      REAL(8),INTENT(IN)    :: FAC_EVP   (IDX)  ! ρε/P_s ... 水蒸気圧 -> Pot
      REAL(8),INTENT(IN)    :: FAC_EVP_2 (IDX)  ! Ps/ε ... 比湿 -> 水蒸気圧
!
      REAL(8),INTENT(IN)    :: XFR_WTR_CNP(IDX)  ! キャ葉   水 被覆率(結露考慮)
      REAL(8),INTENT(IN)    :: XFR_GLA_CNP(IDX)  ! キャ葉   氷 被覆率(結露考慮)
      REAL(8),INTENT(IN)    :: XFR_WTR_GRS(IDX)  ! 下草葉   水 被覆率(結露考慮)
      REAL(8),INTENT(IN)    :: XFR_GLA_GRS(IDX)  ! 下草葉   氷 被覆率(結露考慮)
      REAL(8),INTENT(IN)    :: XFR_WTR_BS (IDX)  ! 裸地or雪 水被覆率 (結露考慮)
      REAL(8),INTENT(IN)    :: XFR_GLA_BS (IDX)  ! 裸地or雪 氷 被覆率(結露考慮)
!                                               ! 注意。雪では XFR_???_SNSK 
!
      REAL(8),INTENT(IN)    :: HUM_REL_SKIN (IDX)   
                                                ! 土壌の比湿相当ファクタ hr
!                                               ! 注意。雪では 1 
      REAL(8)               :: XHUM_REL_SKIN (IDX)   
!
      REAL(8),INTENT(IN)    :: FLAG_DRY_CNP  (IDX)  
      REAL(8),INTENT(IN)    :: FLAG_DRY_SKIN (IDX) 
!
      REAL(8),INTENT(IN)    :: FLAG_TR_ALL (IDX,ICG) ! 蒸散有無フラグ
!
      INTEGER,INTENT(IN)   :: ITR                  
      INTEGER,INTENT(IN)   :: ITR_MAX 

!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: EVAP       (IDX,ICGB) ! 蒸発量
      REAL(8),INTENT(OUT)   :: SUBL       (IDX,ICGB) ! 昇華量
      REAL(8),INTENT(OUT)   :: TRNS       (IDX,ICG) ! 蒸散量
!
      REAL(8),INTENT(OUT)   :: DEVAP_DX (IDX,IEST,ICGB,ICGR)  
!                                            発昇散    キ地参
!                                                 キ下裸
! WORK 
!
      REAL(8),intent(in)    :: E_AIR      (IDX)    ! キャノピー空間空気比湿
#ifdef SIB_DEBUG
      REAL(8)               :: E_AIR_XXX  (IDX)    ! キャノピー空間空気比湿
#endif
      REAL(8)               :: FAC_RGST   (IDX)    
!
      REAL(8)               :: FAC_EV_CNP (IDX)    ! 基本ファクタ
      REAL(8)               :: FAC_TR_CNP (IDX)     
      REAL(8)               :: FAC_SL_CNP (IDX)
      REAL(8)               :: FAC_EV_GRS (IDX)
      REAL(8)               :: FAC_TR_GRS (IDX)
      REAL(8)               :: FAC_SL_GRS (IDX)
      REAL(8)               :: FAC_EV_BS  (IDX)
      REAL(8)               :: FAC_SL_BS  (IDX)
!
      REAL(8)               :: EVAP_CNP_X         ! 蒸発計算用作業 
      REAL(8)               :: EVAP_GRS_X 
      REAL(8)               :: EVAP_BS_X 
!
      REAL(8)               :: DEA_DTC (IDX)     ! DEVAP_DTX 計算用作業 A
      REAL(8)               :: DEA_DTG (IDX)
      REAL(8)               :: DEA_DQR (IDX)
!
      REAL(8)               :: FAC_E_CNP         ! DEVAP_DTX 計算用作業 B
      REAL(8)               :: FAC_S_CNP 
      REAL(8)               :: FAC_T_CNP 
      REAL(8)               :: FAC_E_GRS
      REAL(8)               :: FAC_S_GRS
      REAL(8)               :: FAC_T_GRS
      REAL(8)               :: FAC_E_BS
      REAL(8)               :: FAC_S_BS
!
!     REAL(8)   :: RNSFLAG                   ! = 1-INSFLAG  雪無し1 雪あり0
!
      REAL(8)   :: EVAP_EV_CNP (IDX)         ! 蒸発量、内訳
      REAL(8)   :: EVAP_SL_CNP (IDX)         ! 昇華
      REAL(8)   :: EVAP_TR_CNP (IDX)         ! 蒸散
      REAL(8)   :: EVAP_EV_GRS (IDX)
      REAL(8)   :: EVAP_SL_GRS (IDX)
      REAL(8)   :: EVAP_TR_GRS (IDX)
      REAL(8)   :: EVAP_EV_BS  (IDX)
      REAL(8)   :: EVAP_SL_BS  (IDX)
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      REAL(8)      :: WORK(IDX) 
      CHARACTER(7) :: CMARK 
      INTEGER      :: IFLG
#endif
!
      INTEGER               :: I
!
      LOGICAL,SAVE :: LFIRST
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_EVAP 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!     RNSFLAG = 1. - INSFLAG 
!
!  ===========================================
!  EA （キャノピー空気比湿） の評価 --- 雪なし
!  ===========================================
!
!     大気中への蒸発 = キャノピー・下草・地面からの蒸発散
!
      DO 1000 I=NSTR,NEND 
!
        XHUM_REL_SKIN(I) = HUM_REL_SKIN(I) * FLAG_DRY_SKIN(I) 
     1                     +  (1-FLAG_DRY_SKIN(I))
!
!
!       FAC_AIRN    = R_AIR_INV(I)
!
        FAC_EV_CNP(I) = R_AC_BAR_INV (I) * XFR_WTR_CNP(I)        ! C 蒸発
!
        FAC_SL_CNP(I) = R_AC_BAR_INV (I) * XFR_GLA_CNP(I)        ! C 昇華
!
        FAC_TR_CNP(I) = R_ACX_BAR_INV(I)                         ! C 蒸散
     1                  * (1-XFR_WTR_CNP(I)-XFR_GLA_CNP(I)) 
!!!  1                  * FLAG_DRY_CNP(I)
     1                  * FLAG_TR_ALL(I,1)
!
        FAC_EV_GRS(I) = R_AG_INV (I) * XFR_WTR_GRS(I)            ! GRS 蒸発
     1                               * CVR_VEG_ms(I,2)              ! 雪0

        FAC_SL_GRS(I) = R_AG_INV (I) * XFR_GLA_GRS(I)            ! GRS 昇華
     1                               * CVR_VEG_ms(I,2)              ! 雪0

        FAC_TR_GRS(I) = R_AGX_BAR_INV(I)                         ! 雪0
     1                  * (1-XFR_WTR_GRS(I)-XFR_GLA_GRS(I))      ! GRS 蒸散
!!!  1                  * FLAG_DRY_SKIN(I)
     1                  * FLAG_TR_ALL(I,2)
!
        FAC_EV_BS (I) = R_AB_INV(I) * XFR_WTR_BS(I) 
     1                                 * (1-CVR_VEG_ms(I,2))    ! 裸地or雪 蒸発
        FAC_SL_BS (I) = R_AB_INV(I) * XFR_GLA_BS(I) 
     1                                 * (1-CVR_VEG_ms(I,2))    ! 裸地or雪 昇華
!
        FAC_RGST(I) = 1 / (    
     1                    R_AIR_INV(I)                       ! A->R
     1                  + FAC_EV_CNP (I)
     1                  + FAC_SL_CNP (I)
     1                  + FAC_TR_CNP (I)
     1                  + FAC_EV_GRS (I) 
     1                  + FAC_SL_GRS (I) 
     1                  + FAC_TR_GRS (I) 
     1                  + FAC_EV_BS  (I) 
     1                  + FAC_SL_BS  (I) )
!
!  keisanzumi
!
#ifdef SIB_DEBUG
        E_AIR_XXX (I) = FAC_RGST(I) 
     1    * (
     1         E_REF(I) * R_AIR_INV(I)
     1       + E_TC_SAT(I) 
     1         * ( FAC_EV_CNP(I) + FAC_SL_CNP(I) + FAC_TR_CNP(I) ) 
     1       + E_TG_SAT(I) 
     1         * (   ( FAC_EV_GRS(I) + FAC_SL_GRS(I) + FAC_TR_GRS(I) )      
     1             + ( FAC_EV_BS(I)  + FAC_SL_BS(I) ) 
     1          * XHUM_REL_SKIN(I))      
     1      )
#endif
!
!
#ifdef SIB_DEBUG 
       if ( abs((e_air(i)-e_air_xxx(i))/e_air(i)) .gt. 1.D-6 ) then
        write(6,*) 'CCCXXXX canopy_evap warning e_air diff ',
     1             i,e_air(i),e_air_xxx(i),
     1             abs((e_air(i)-e_air_xxx(i))/e_air(i))
!       stop 999
       endif
#endif
!
 1000 CONTINUE 
!
!  ======================================================
!   蒸発散評価 
!  ======================================================
!
      DO 2000 I=NSTR,NEND 
!
!   蒸発量・昇華量・蒸散量
!
        EVAP_CNP_X = FAC_EVP(I) * ( E_AIR(I) - E_TC_SAT(I) )   
!
        EVAP_GRS_X = FAC_EVP(I) * ( E_AIR(I) - E_TG_SAT(I) )   
        EVAP_BS_X  = FAC_EVP(I) 
     1                   * ( E_AIR(I) - XHUM_REL_SKIN(I)*E_TG_SAT(I) ) 
!       
#ifdef KILLEVAP
          fac_ev_cnp(i) = 0 
          fac_tr_cnp(i) = 0
          fac_sl_cnp(i) = 0
          fac_ev_grs(i) = 0
          fac_tr_grs(i) = 0
          fac_sl_grs(i) = 0
          fac_ev_bs (i) = 0
          fac_sl_bs (i) = 0
#endif
!
        EVAP_EV_CNP (I) = EVAP_CNP_X * FAC_EV_CNP(I) 
        EVAP_SL_CNP (I) = EVAP_CNP_X * FAC_SL_CNP(I) 
        EVAP_TR_CNP (I) = EVAP_CNP_X * FAC_TR_CNP(I) 
!

!       
        EVAP_EV_GRS (I) = EVAP_GRS_X * FAC_EV_GRS(I)              ! 雪0
        EVAP_SL_GRS (I) = EVAP_GRS_X * FAC_SL_GRS(I)              ! 雪0
        EVAP_TR_GRS (I) = EVAP_GRS_X * FAC_TR_GRS(I)              ! 雪0
!       
        EVAP_EV_BS  (I) = EVAP_BS_X  * FAC_EV_BS (I)
        EVAP_SL_BS  (I) = EVAP_BS_X  * FAC_SL_BS (I)
!
!   蒸発散量の評価 （時間重み、HL をかけていないことに注意）
!
!    EVAP    (IDX,ICGB) 
!    SUBL    (IDX,ICGB) 
!    TRNS    (IDX,ICG) 
!
!
        EVAP(I,1) = EVAP_EV_CNP(I) 
        TRNS(I,1) = EVAP_TR_CNP(I)  
        SUBL(I,1) = EVAP_SL_CNP(I) 
!
        EVAP(I,2) = EVAP_EV_GRS(I) 
        TRNS(I,2) = EVAP_TR_GRS(I) 
        SUBL(I,2) = EVAP_SL_GRS(I)
!
        EVAP(I,3) = EVAP_EV_BS (I)       
        SUBL(I,3) = EVAP_SL_BS (I) 
!
 2000 CONTINUE
!
!  ======================================================
!   蒸発散の各物理量依存性 
!  ======================================================
!
      DO 3000 I=NSTR,NEND
!
!  ---------------------------------------------------------
!  作業変数  FAC_x_yyy .. 単位 Ty 変化で、x からの蒸発等増分
!  ---------------------------------------------------------
!
        FAC_E_CNP = FAC_EVP(I) * FAC_EV_CNP(I)       ! キャ蒸発
        FAC_S_CNP = FAC_EVP(I) * FAC_SL_CNP(I)       ! 
        FAC_T_CNP = FAC_EVP(I) * FAC_TR_CNP(I)       ! 
!
        FAC_E_GRS = FAC_EVP(I) * FAC_EV_GRS(I)       !          雪では0
        FAC_S_GRS = FAC_EVP(I) * FAC_SL_GRS(I)       !          雪では0
        FAC_T_GRS = FAC_EVP(I) * FAC_TR_GRS(I)       !          雪では0
!
!       FAC_E_BS  = FAC_EVP(I) * FAC_EV_BS (I) * XHUM_REL_SKIN(I)   
!       FAC_S_BS  = FAC_EVP(I) * FAC_SL_BS (I) * XHUM_REL_SKIN(I)  
        FAC_E_BS  = FAC_EVP(I) * FAC_EV_BS (I) 
        FAC_S_BS  = FAC_EVP(I) * FAC_SL_BS (I) 
!
!  --------------------------------
!  EA の TC,TG,TS,QR 依存性 DEA_Dxx 
!  --------------------------------  
!
        DEA_DTC(I)
     1     = FAC_RGST(I) * DEC_SAT_DTC(I) 
     1          * (   R_AC_BAR_INV(I)  
     1                    * (  XFR_WTR_CNP(I)+XFR_GLA_CNP(I))     ! 蒸発・昇華
     1              + R_ACX_BAR_INV(I) 
     1                    * (1-XFR_WTR_CNP(I)-XFR_GLA_CNP(I)) )  ! 蒸散
!    1          * FLAG_TMPMOD_CNP (I)
!
        DEA_DTG(I)
     1   = FAC_RGST(I) * DEG_SAT_DTG(I) 
     1       * (                                        
     1          CVR_VEG_ms(I,2) * R_AG_INV (I) 
     1               * (  XFR_WTR_GRS(I)+XFR_GLA_GRS(I))        ! 下草蒸発 雪0
     1        + R_AGX_BAR_INV(I) 
     1               * (1-XFR_WTR_GRS(I)-XFR_GLA_GRS(I))        ! 下草蒸散 雪0
     1        +(1-CVR_VEG_ms(I,2))                              ! 裸地蒸発
     1               * (  R_AB_INV (I) * XHUM_REL_SKIN(I) )       
     1         )
!    1         * FLAG_TMPMOD_SKIN (I)
!
        DEA_DQR(I) = FAC_RGST(I) * R_AIR_INV (I) 
     1               * FAC_EVP_2(I)                       
!
!  ------------------------------------------------------
!  蒸発散の各物理量依存性 DEVAP_DX(IDX,IEST,ICGB,ICGR)
!  ------------------------------------------------------  
!
!    TC 依存性
!
        DEVAP_DX(I,1,1,1) = FAC_E_CNP * ( DEA_DTC(I) - DEC_SAT_DTC(I) )
        DEVAP_DX(I,2,1,1) = FAC_S_CNP * ( DEA_DTC(I) - DEC_SAT_DTC(I) )
        DEVAP_DX(I,3,1,1) = FAC_T_CNP * ( DEA_DTC(I) - DEC_SAT_DTC(I) )
!
        DEVAP_DX(I,1,2,1) = FAC_E_GRS * ( DEA_DTC(I) - 0 )  
        DEVAP_DX(I,2,2,1) = FAC_S_GRS * ( DEA_DTC(I) - 0 )  
        DEVAP_DX(I,3,2,1) = FAC_T_GRS * ( DEA_DTC(I) - 0 )  
!                                      ! 以上 3 項、雪では 0 
!
        DEVAP_DX(I,1,3,1) = FAC_E_BS  * ( DEA_DTC(I) - 0 )
        DEVAP_DX(I,2,3,1) = FAC_S_BS  * ( DEA_DTC(I) - 0 )
        DEVAP_DX(I,3,3,1) = 0. 
!
!    TG 依存性
! 
        DEVAP_DX(I,1,1,2) = FAC_E_CNP * ( DEA_DTG(I) - 0 )       
        DEVAP_DX(I,2,1,2) = FAC_S_CNP * ( DEA_DTG(I) - 0 )
        DEVAP_DX(I,3,1,2) = FAC_T_CNP * ( DEA_DTG(I) - 0 )
!
        DEVAP_DX(I,1,2,2) = FAC_E_GRS * ( DEA_DTG(I) - DEG_SAT_DTG(I) )  
        DEVAP_DX(I,2,2,2) = FAC_S_GRS * ( DEA_DTG(I) - DEG_SAT_DTG(I) )  
        DEVAP_DX(I,3,2,2) = FAC_T_GRS * ( DEA_DTG(I) - DEG_SAT_DTG(I) )  
!                                      ! 以上 3 項、雪では 0 
!
!       DEVAP_DX(I,1,3,2) = FAC_E_BS  * ( DEA_DTG(I) - DEG_SAT_DTG(I) )
!       DEVAP_DX(I,2,3,2) = FAC_S_BS  * ( DEA_DTG(I) - DEG_SAT_DTG(I) )
        DEVAP_DX(I,1,3,2) = FAC_E_BS  
     1                      * ( DEA_DTG(I) 
     1                          - DEG_SAT_DTG(I) * XHUM_REL_SKIN(I)   )
        DEVAP_DX(I,2,3,2) = FAC_S_BS  
     1                      * ( DEA_DTG(I) 
     1                          - DEG_SAT_DTG(I) * XHUM_REL_SKIN(I)   )
        DEVAP_DX(I,3,3,2) = 0. 
!
!    QR 依存性
! 
        DEVAP_DX(I,1,1,3) = FAC_E_CNP * ( DEA_DQR(I) - 0 )
        DEVAP_DX(I,2,1,3) = FAC_S_CNP * ( DEA_DQR(I) - 0 )
        DEVAP_DX(I,3,1,3) = FAC_T_CNP * ( DEA_DQR(I) - 0 )
!
        DEVAP_DX(I,1,2,3) = FAC_E_GRS * ( DEA_DQR(I) - 0 )   ! 雪では 0      
        DEVAP_DX(I,2,2,3) = FAC_S_GRS * ( DEA_DQR(I) - 0 )   ! 雪では 0
        DEVAP_DX(I,3,2,3) = FAC_T_GRS * ( DEA_DQR(I) - 0 )   ! 雪では 0
!
        DEVAP_DX(I,1,3,3) = FAC_E_BS  * ( DEA_DQR(I) - 0 )
        DEVAP_DX(I,2,3,3) = FAC_S_BS  * ( DEA_DQR(I) - 0 )
        DEVAP_DX(I,3,3,3) = 0. 
!
 3000 CONTINUE
!
! ================
! >> 詳細モニタ <<
! ================
! 
#ifdef MONYOS
      IF ( ITR.EQ.ITR_MAX ) THEN
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! --------------------------
! > 各種フラックス陽的評価 <
! --------------------------
!
!  CNP 蒸発
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'WEVN0'
          CALL MONIT_SIB ( CMARK , EVAP_EV_CNP ,ICN, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'WEVS0'
          CALL MONIT_SIB ( CMARK , EVAP_EV_CNP ,ICS, DELT_CNP , ONE )      
        ENDIF
!
!  GRS 蒸発
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'WEVG0'
          CALL MONIT_SIB ( CMARK , EVAP_EV_GRS ,ICN, DELT_CNP , ONE )      
        ENDIF
!
!  BS 蒸発
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'WEVL0'
          CALL MONIT_SIB ( CMARK , EVAP_EV_BS ,ICN, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'WEVW0'
          CALL MONIT_SIB ( CMARK , EVAP_EV_BS ,ICS, DELT_CNP , ONE )      
        ENDIF
!
!  CNP 昇華
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'GSLN0'
          CALL MONIT_SIB ( CMARK , EVAP_SL_CNP ,ICN, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'GSLS0'
          CALL MONIT_SIB ( CMARK , EVAP_SL_CNP ,ICS, DELT_CNP , ONE )      
        ENDIF
!
!  GRS 昇華
!
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'GSLG0'
          CALL MONIT_SIB ( CMARK , EVAP_SL_GRS ,ICN, DELT_CNP , ONE )      
        ENDIF
!
!  BS 昇華
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'GSLL0'
          CALL MONIT_SIB ( CMARK , EVAP_SL_BS ,ICN, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'GSLW0'
          CALL MONIT_SIB ( CMARK , EVAP_SL_BS ,ICS, DELT_CNP , ONE )      
        ENDIF
!
! 蒸散
!
!   キャノピー
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'WTRBNL0'
          CALL MONIT_SIB ( CMARK , EVAP_TR_CNP ,ICN, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'WTRBSL0'
          CALL MONIT_SIB ( CMARK , EVAP_TR_CNP ,ICS, DELT_CNP , ONE )      
        ENDIF
!
!    下草
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'WTRBGL0'
          CALL MONIT_SIB ( CMARK , EVAP_TR_GRS ,ICN, DELT_CNP , ONE )      
        ENDIF
!
! -----------------------------------
! > 裸地・積雪上からの蒸発/昇華詳細 <
! -----------------------------------
!
        IF ( INSFLAG .EQ. 0 ) THEN
          IFLG = ICN
        ELSE
          IFLG = ICS
        ENDIF
!
        DO I=NSTR,NEND
          WORK(I) = E_AIR(I) - XHUM_REL_SKIN(I)*E_TG_SAT(I)  
        ENDDO 
!
        CMARK  = 'EAIR'
        CALL MONIT_SIB ( CMARK , E_AIR , IFLG , DELT_CNP , ONE )      
        CMARK  = 'RELSKIN'
        CALL MONIT_SIB ( CMARK , XHUM_REL_SKIN, IFLG, DELT_CNP, ONE )      
        CMARK  = 'ETGS'
        CALL MONIT_SIB ( CMARK , E_TG_SAT, IFLG, DELT_CNP, ONE )      
        CMARK  = 'ELDIF'
        CALL MONIT_SIB ( CMARK , WORK, IFLG, DELT_CNP, ONE )      
!
        DO I=NSTR,NEND
          WORK(I) = R_AB_INV(I) * XFR_WTR_BS(I) 
     1                                 * (1-CVR_VEG_ms(I,2))
        ENDDO 
!
        CMARK  = 'ELFACT'
        CALL MONIT_SIB ( CMARK , WORK , IFLG , DELT_CNP , ONE )      
        CMARK  = 'RABI'
        CALL MONIT_SIB ( CMARK , R_AB_INV, IFLG, DELT_CNP, ONE )      
        CMARK  = 'FRWGZ'
        CALL MONIT_SIB ( CMARK , XFR_WTR_BS, IFLG, DELT_CNP, ONE )      
        CMARK  = 'FRCVRG'
        CALL MONIT_SIB ( CMARK , CVR_VEG_ms(1,2),IFLG,DELT_CNP,ONE )      
!
! ------------------------------
! > キャノピーからの蒸発・蒸散 <
! ------------------------------
!
        DO I=NSTR,NEND
          WORK(I) = E_AIR(I) - E_TC_SAT(I)  
        ENDDO 
!
        CMARK  = 'ETCS'
        CALL MONIT_SIB ( CMARK , E_TC_SAT, IFLG, DELT_CNP, ONE )      
!
        CMARK  = 'ECDIF'
        CALL MONIT_SIB ( CMARK , WORK, IFLG, DELT_CNP, ONE )      
!
!       FAC_TR_CNP(I) = R_ACX_BAR_INV(I) * FLAG_TR_ALL(I,1)
!    1                  * (1-XFR_WTR_CNP(I)-XFR_GLA_CNP(I)) 
!
        DO I=NSTR,NEND
          WORK(I) = R_ACX_BAR_INV(I) * FLAG_TR_ALL(I,1)
        ENDDO 
        CMARK   = 'RACXI'
        CALL MONIT_SIB ( CMARK , WORK, IFLG, DELT_CNP, ONE )      
!
      ENDIF
!
      ENDIF
!
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_CANOPY_EVAP
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
#endif
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_EVAP 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! ----------------------------------------------------------------------
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! --------------------------
! > 各種フラックス陽的評価 <
! --------------------------
!
! CNP での蒸発量 E_CNP 
!
      CMARK  = 'WEVN0'
      CTITLE = 'CANOPY EVAPORATION NOS explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WEVS0'
      CTITLE = 'CANOPY EVAPORATION SNW explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草 での蒸発量 E_GRS
!
      CMARK  = 'WEVG0'
      CTITLE = 'GRASS EVAPORATION NOS explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 裸地 での蒸発量 E_BS
!
      CMARK  = 'WEVL0'
      CTITLE = 'BARE SOIL EVAPORATION explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 雪面での蒸発量 E_BS
!
      CMARK  = 'WEVW0'
      CTITLE = 'SNOW EVAPORATION explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! CNP での昇華量 SUBL_CNP 
!
      CMARK  = 'GSLN0'
      CTITLE = 'CANOPY SUBL NOS explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GSLS0'
      CTITLE = 'CANOPY SUBL SNW explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草 での昇華量 SUBL_GRS
!
      CMARK  = 'GSLG0'
      CTITLE = 'GRASS SUBL explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 裸地 での昇華量 SUBL_BS
!
      CMARK  = 'GSLL0'
      CTITLE = 'BARE SOIL SUBL explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 雪 での昇華量 SUBL_BS
!
      CMARK  = 'GSLW0'
      CTITLE = 'SNOW SUBL explicit'
      CUNIT  = 'Kg/M**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!
! 蒸散
!
      CMARK  = 'WTRBNL0'
      CTITLE = 'TRANS CANOPY NOS explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBSL0'
      CTITLE = 'TRANS CANOPY SNW explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRBGL0'
      CTITLE = 'TRNS GRASS explicit'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! -----------------------------------
! > 裸地・積雪上からの蒸発/昇華詳細 <
! -----------------------------------
!
!  水蒸気圧
!
      CMARK  = 'EAIR'
      CTITLE = 'E_AIR'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ETGS'
      CTITLE = 'Saturated E of TG '
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ETCS'
      CTITLE = 'Saturated E of TC '
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  裸地/雪蒸発
!
      CMARK  = 'ELDIF'
      CTITLE = 'EA-EG'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RELSKIN'
      CTITLE = 'RELATIVE HUMIDITY OF SOIL 1E_AIR'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ELFACT'
      CTITLE = 'Factor from BS/SNOW Evap'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RABI'
      CTITLE = 'Inverse of Resist of AB'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRWGZ'
      CTITLE = 'Fraction of BS/Snow surf (normal=1)'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRCVRG'
      CTITLE = 'Fraction of CVRG'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  キャノピー 蒸発散
!
      CMARK  = 'ECDIF'
      CTITLE = 'EA-ETC'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RACXI'
      CTITLE = 'Inverse of Resist of ACX'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
! 
      ENDIF
#endif
!
      RETURN
      END

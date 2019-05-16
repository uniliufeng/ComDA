      module CANOPY_INTEGRATION_OVERIMP 
!
      contains 
!
! オーバーインプリシット積分用積分ツール
!
! 作戦
!   フラックスと変化量の評価は
!        C dT/dt = F ( T + dT/dt x α Δt ) 、α～1.5  
!   を用いる。この右辺がフラックス収束、左辺/C x Δt が変化量。
!
! 構成
!   canopy_main 
!     
!     canopy_integration_overimp_itr    仮積分時に呼ばれるサブルーチン ok 
!     canopy_integration_overimp        積分時に呼ばれるサブルーチン
!       canopy_integration_overimp_pre  前処理     
!       canopy_integration_overimp_veg  植生   
!       canopy_integration_overimp_snow 雪
!
      SUBROUTINE CANOPY_INTEGRATION_OVERIMP_RUN ( 
     I   NSTR       , NEND        , INSFLAG    , IMSKLT , 
     I   DELT_CNP   , DELT_ATM    , 
     I   SMTX       , 
     I   EVAP       , SUBL        , TRNS       , 
     I   SENS       , 
     I   RAD_LONG   , 
     I   RAD_SHRT   , 
     I   HDF        , 
     I   DEVAP_DX   , DSENS_DX    , DRAD_DX    , HDF_SKIN  , 
     I   H2O_SAT_CNP    , H2O_SAT_GRS    , 
     I   CVR_S_or_N     , XLAI_VEG_msx   , 
     I   TMP_CNP_OLD    , 
     I   FR_H2O_CNP     , FR_WBH_CNP     , 
     I   TMP_SKIN_OLD   , 
     I   FR_H2O_GRS     , FR_WBH_GRS     , 
     I   TMP_BOTTOM     , 
     I   RAD_SHRT_INTO_SNOW_LAYER   , 
     I   CAP_SS1        , FLAG_PHASE_SS1 , 
     I   GMT_IN         , GMQ_IN     , 
     O   TMP_CNP_NEXT   , FR_WTR_CNP_NEXT, FR_GLA_CNP_NEXT, 
     O   TMP_SKIN_NEXT  , FR_WTR_GRS_NEXT, FR_GLA_GRS_NEXT, 
     O   TRNS_TMP       , 
     O   GMT_OUT        , GMQ_OUT    , 
     O   FLLH           , FLSH       , qflux , 
     M   RAD_LONG_NorS_2_ATM  , 
     M   ENG_BUCKET                  , 
     M   H2O_BUCKET_INTO             , 
     M   H2O_BUCKET_FROM             , 
     M   ENH_BUCKET_INTO             , 
     M   ENH_BUCKET_FROM             , 
     M   SRAD_BUCKET_SNW , 
     M   TMP_SS1         ) 

!
      USE SIBPRM , ONLY :
     1    IDX  , ICG , ICGB  , IEST , ICGR  
!
      USE SIBCON , ONLY :
     1   tmp_frez_c , hl_wtr_c , hl_gla_c 
!
      USE SIB_MONIT , ONLY :
     1   IMONIT_LEVEL  ,
     1   IMONIT_SELECT ,
     1   IMONIT_ALL
!
      IMPLICIT NONE 
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: INSFLAG
      INTEGER,INTENT(IN)   :: IMSKLT(IDX)
!
!     REAL(8),INTENT(IN)   :: SMTX       (IDX,4,5)
      REAL(8),INTENT(IN)   :: SMTX       (IDX,5,6)
!
      REAL(8),INTENT(IN)   :: DELT_CNP
      REAL(8),INTENT(IN)   :: DELT_ATM
      REAL(8),INTENT(IN)   :: XLAI_VEG_msx(IDX,ICG) 
      REAL(8),INTENT(IN)   :: CVR_S_OR_N  (IDX) 
!
      REAL(8),INTENT(IN)   :: TMP_BOTTOM  (IDX)
      REAL(8),INTENT(IN)   :: H2O_SAT_GRS (IDX)
      REAL(8),INTENT(IN)   :: H2O_SAT_CNP (IDX)
!
      REAL(8),INTENT(IN)   :: TMP_CNP_OLD   (IDX)
      REAL(8),INTENT(IN)   :: TMP_SKIN_OLD  (IDX)
      REAL(8),INTENT(IN)   :: FR_H2O_CNP (IDX)  ! 
      REAL(8),INTENT(IN)   :: FR_H2O_GRS (IDX)  ! 
      REAL(8),INTENT(IN)   :: FR_WBH_CNP (IDX)  ! 
      REAL(8),INTENT(IN)   :: FR_WBH_GRS (IDX)  ! 
!
      REAL(8),INTENT(IN)   :: EVAP       (IDX,ICGB)
      REAL(8),INTENT(IN)   :: SUBL       (IDX,ICGB)
      REAL(8),INTENT(IN)   :: TRNS       (IDX,ICG)
      REAL(8),INTENT(IN)   :: DEVAP_DX   (IDX,IEST,ICGB,ICGR)
!
      REAL(8),INTENT(IN)   :: SENS       (IDX,ICG)
      REAL(8),INTENT(IN)   :: DSENS_DX   (IDX,ICG,ICGR) 
!
      REAL(8),INTENT(IN)   :: RAD_LONG   (IDX,ICG)     ! 放射
      REAL(8),INTENT(IN)   :: RAD_SHRT   (IDX,ICG)     ! 放射
      REAL(8),INTENT(IN)   :: DRAD_DX    (IDX,ICGR,ICG) ! 放射吸収量変化率
!
      REAL(8),INTENT(IN)   :: HDF        (IDX,2) 
      REAL(8),INTENT(IN)   :: HDF_SKIN   (IDX,2) 
!
      REAL(8),INTENT(IN)   :: RAD_SHRT_INTO_SNOW_LAYER (IDX)
!
!     REAL(8),INTENT(IN)   :: FLAG_TMPMOD_SKIN     (IDX) 
!     REAL(8),INTENT(IN)   :: FLAG_TMPMOD_CNP      (IDX) 
!
      REAL(8),INTENT(IN)   :: CAP_SS1 (IDX) 
      REAL(8),INTENT(IN)   :: FLAG_PHASE_SS1 (IDX) 
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: TMP_CNP_NEXT  (IDX)    ! 
      REAL(8),INTENT(OUT)  :: TMP_SKIN_NEXT (IDX)    ! 
!
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_NEXT (IDX)  ! 
      REAL(8),INTENT(OUT)  :: FR_WTR_GRS_NEXT (IDX)  ! 
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_NEXT (IDX)  ! 
      REAL(8),INTENT(OUT)  :: FR_GLA_GRS_NEXT (IDX)  ! 
!
      REAL(8),INTENT(OUT)  :: TRNS_TMP (IDX,ICG ) 
!
      REAL(8),INTENT(INOUT)  :: RAD_LONG_NorS_2_ATM (IDX)  
      REAL(8),INTENT(IN)     :: GMT_IN (IDX,2) 
      REAL(8),INTENT(IN)     :: GMQ_IN (IDX,2) 
      REAL(8),INTENT(OUT)    :: GMT_OUT(IDX  )  
      REAL(8),INTENT(OUT)    :: GMQ_OUT(IDX  )
!
      REAL(8),INTENT(OUT)    :: FLSH(IDX)
      REAL(8),INTENT(OUT)    :: FLLH(IDX)    ! positive = upward 
      REAL(8),INTENT(INOUT)  :: qflux(IDX)   ! positive = upward 

!
! WORK
!
      REAL(8)  :: EVAP_TMP    (IDX,ICGB) 
      REAL(8)  :: SUBL_TMP    (IDX,ICGB) 
      REAL(8)  :: RAD_LONG_TMP(IDX,ICG ) 
      REAL(8)  :: SENS_TMP    (IDX,ICG ) 
      REAL(8)  :: HDF_TMP     (IDX,2   ) 
      REAL(8)  :: DUMMY_XXX1 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_XXX2 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_XXX3 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_XXX4 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_XXX5 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_XXX6 (IDX)   = 0.D0 
      REAL(8)  :: DUMMY_IDX (IDX)   = 0.D0
      REAL(8)  :: DUMMY_IDX2(IDX,2) = 0.D0
!
      REAL(8)  :: DRAD_DX_CONV (IDX)       

!
! INOUT 各種バケツ
!
      REAL(8),INTENT(INOUT):: ENG_BUCKET            (IDX)
      REAL(8),INTENT(INOUT):: H2O_BUCKET_INTO       (IDX)
      REAL(8),INTENT(INOUT):: H2O_BUCKET_FROM       (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_INTO       (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_FROM       (IDX)
      REAL(8),INTENT(INOUT):: SRAD_BUCKET_SNW       (IDX)
!
      REAL(8),INTENT(INOUT):: TMP_SS1               (IDX)
!                                                
! WORK 
!
      INTEGER   :: I
      REAL(8)   :: WTR_LEAF_DROP(IDX,ICG)
      REAL(8)   :: GLA_LEAF_DROP(IDX,ICG)
!
      REAL(8)   :: XDTMP_CNP (IDX) 
      REAL(8)   :: XDTMP_SKIN(IDX) 
!
      LOGICAL,SAVE :: LFIRST = .TRUE. 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP 00/02/23 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! ==============
! >> 初期処理 <<
! ==============
!   フラックスの見積もり。インプリシット補正分を加える。
!
      CALL CANOPY_INTEGRATION_OVERIMP_PRE (
     I   NSTR       , NEND        , INSFLAG    , 
     I   DELT_CNP   , DELT_ATM    , 
     I   SMTX       , 
     I   EVAP       , SUBL        , TRNS       , 
     I   SENS       , 
     I   RAD_LONG   , 
!    I   RAD_SHRT    , 
     I   HDF        , 
     I   DEVAP_DX   , DSENS_DX    , DRAD_DX    , HDF_SKIN  , 
     O   EVAP_TMP   , SUBL_TMP    , TRNS_TMP   , 
     O   RAD_LONG_TMP , 
     O   SENS_TMP    , HDF_TMP    ,
     O   XDTMP_CNP   , XDTMP_SKIN , 
     M   RAD_LONG_NorS_2_ATM )
!
#ifdef MONYOS
! ========================
! >> フラックスチェック <<    
! ========================
!   収支チェックには不要。
!     この後にフラックス修正作業がされるが、
!     その大きさをチェックするのに使われる。
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SELECT ) THEN
!
!  熱フラックスインプリシット補正後の量
!
        CALL CANOPY_MNT_HEAT_FLUX (
     I    DELT_CNP , 
     I    NSTR     , NEND      ,  
     I    INSFLAG  , 0         , 
     I    SENS_TMP , HDF_TMP   , 
!    I    RAD_SHRT , rad_shrt_into_snow_layer ,  
     I    RAD_LONG_TMP         ,
     I    RAD_LONG_NorS_2_ATM  ,                    ! ダミー
     O    FLSH     )
!
!  水フラックスインプリシット補正後の量
!
        CALL CANOPY_MNT_WATER_FLUX_CORR (
     I    DELT_CNP , 
     I    INSFLAG  , 
     I    EVAP_TMP , SUBL_TMP  , TRNS_TMP )
!
      ENDIF
#endif
!
! ====================
! >> キャノピー積分 << 
! ====================
!
        DO I=NSTR,NEND
          DRAD_DX_CONV(I) = DRAD_DX(I,3,1) - DRAD_DX(I,2,1)
        ENDDO
!
        CALL CANOPY_INTEGRATION_OVERIMP_VEG (
     I   NSTR         , NEND         , 1            , CVR_S_OR_N   , 
     I   INSFLAG      , IMSKLT       , 
     I   TMP_CNP_OLD  , XDTMP_CNP    , 
     I   FR_H2O_CNP   , FR_WBH_CNP   , 
     I   DUMMY_IDX    , H2O_SAT_CNP  , DELT_CNP    , XLAI_VEG_msx(1,1), 
     M   RAD_LONG_TMP(1,1) , 
     I   RAD_SHRT(1,1) , 
     I   SENS_TMP(1,1), DUMMY_XXX1   , 
     I   DRAD_DX_CONV  , 
     I   DSENS_DX(1,1,1), DUMMY_IDX2 , 
     I   TRNS_TMP(1,1) ,  
     M   EVAP_TMP(1,1) , SUBL_TMP(1,1), 
     I   DUMMY_IDX     , DUMMY_IDX    , 
     O   TMP_CNP_NEXT  , FR_WTR_CNP_NEXT, FR_GLA_CNP_NEXT, 
     O   WTR_LEAF_DROP(1,1) , GLA_LEAF_DROP(1,1) , 
     M   RAD_LONG_NorS_2_ATM,   
     M   DUMMY_XXX2   , 
     M   H2O_BUCKET_INTO              , 
     M   DUMMY_XXX3                   , 
     M   ENH_BUCKET_INTO              , 
     M   DUMMY_XXX5                   ) 
!
! ====================
! >> 下草・裸地積分 << 
! ====================
!
      IF ( INSFLAG .EQ. 0 ) THEN 
!
!  スキン
!
        DO I=NSTR,NEND
          DRAD_DX_CONV(I) = DRAD_DX(I,1,2) + DRAD_DX(I,3,2)
        ENDDO
!
        CALL CANOPY_INTEGRATION_OVERIMP_VEG (
     I   NSTR         , NEND         , 2            , CVR_S_OR_N     , 
     I   INSFLAG      , IMSKLT       , 
     I   TMP_SKIN_OLD     , XDTMP_SKIN   , 
     I   FR_H2O_GRS   , FR_WBH_GRS   , 
     I   TMP_BOTTOM   , H2O_SAT_GRS  , DELT_CNP     , XLAI_VEG_msx(1,2), 
     I   RAD_LONG_TMP(1,2) , 
     I   RAD_SHRT(1,2) , 
     M   SENS_TMP(1,2), HDF_TMP      , 
     I   DRAD_DX_CONV, 
     I   DSENS_DX(1,2,2), HDF_SKIN  , 
     I   TRNS_TMP(1,2),  
     M   EVAP_TMP(1,2), SUBL_TMP(1,2), 
     I   EVAP_TMP(1,3), SUBL_TMP(1,3), 
     O   TMP_SKIN_NEXT, FR_WTR_GRS_NEXT, FR_GLA_GRS_NEXT, 
     O   WTR_LEAF_DROP(1,2) , GLA_LEAF_DROP(1,2) , 
     M   RAD_LONG_NorS_2_ATM,   
     M   ENG_BUCKET                  , 
     M   H2O_BUCKET_INTO             , 
     M   H2O_BUCKET_FROM             , 
     M   ENH_BUCKET_INTO             , 
     M   ENH_BUCKET_FROM             ) 
!
!  土壌第一層
!
        CALL CANOPY_INTEGRATION_SS1 (
     I    NSTR     , NEND    , INSFLAG    ,
     I    DELT_CNP , CAP_SS1 , HDF_TMP    , FLAG_PHASE_SS1 , 
     M    TMP_SS1 ) 
!
! ============
! >> 雪積分 << 
! ============
!
      ELSE
!
        DO I=NSTR,NEND
          DRAD_DX_CONV(I) = DRAD_DX(I,1,2) + DRAD_DX(I,3,2)
        ENDDO

!  雪スキン
!
        CALL CANOPY_INTEGRATION_OVERIMP_SNOW (
     I   NSTR       , NEND      , 
     I   DELT_CNP   , CVR_S_OR_N, 
     I   TMP_SKIN_OLD   , XDTMP_SKIN , 
     I   TMP_BOTTOM ,              
     I   EVAP_TMP(1,3)  , SUBL_TMP(1,3), 
     M   SENS_TMP(1,2)  , 
     M   RAD_LONG_TMP (1,2)   , 
     I   RAD_SHRT(1,2) , 
     M   HDF_TMP   , 
     I   DSENS_DX(1,2,2), 
     I   DRAD_DX_CONV       , 
     I   HDF_SKIN  , 
     I   RAD_SHRT_INTO_SNOW_LAYER    ,
     O   TMP_SKIN_NEXT           , 
     M   RAD_LONG_NorS_2_ATM,   
     M   ENG_BUCKET              , 
     M   H2O_BUCKET_INTO         , H2O_BUCKET_FROM        , 
     M   ENH_BUCKET_INTO         , ENH_BUCKET_FROM        ,
     M   SRAD_BUCKET_SNW   )
!
!  雪面第一層
!
        CALL CANOPY_INTEGRATION_SS1 (
     I    NSTR     , NEND    , INSFLAG    , 
     I    DELT_CNP , CAP_SS1 , HDF_TMP    , FLAG_PHASE_SS1 , 
     M    TMP_SS1 ) 
!
      ENDIF
!
! ====================
! >> GMT,GMQ の計算 << 
! ====================
!
      CALL CANOPY_INTEGRATION_GET_GMTQ (
     I     NSTR , NEND , 
     I     SENS_TMP , EVAP_TMP , SUBL_TMP , TRNS_TMP ,  
     I     GMT_IN   , GMQ_IN   , 
     O     GMT_OUT  , GMQ_OUT   )
!
! ========================
! >> フラックスチェック << 
! ========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SELECT ) THEN  ! 必ずとおる
!
! 熱フラックス最終量
!
        CALL CANOPY_MNT_HEAT_FLUX (
     I    DELT_CNP , 
     I    NSTR     , NEND      ,  
     I    INSFLAG  , 1         , 
     I    SENS_TMP , HDF_TMP   , 
     I    RAD_LONG_TMP        ,
     I    RAD_LONG_NorS_2_ATM ,   
     O    FLSH     )
!
! 水フラックス最終量
!
        CALL CANOPY_MNT_WATER_FLUX_FINAL (
     I    NSTR     , NEND      , DELT_CNP , 
     I    INSFLAG  , 
     I    EVAP_TMP , SUBL_TMP  , TRNS_TMP ,  
     I    TMP_CNP_NEXT , TMP_SKIN_NEXT , TMP_BOTTOM ,
     O    FLLH , qflux )
!
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_SELECT ) THEN  ! 必ずとおる
!
! 水フラックス落下量
!
        CALL CANOPY_MNT_WATER_DROP (
     I    NSTR     , NEND      , 
     I    INSFLAG  , 
     I    WTR_LEAF_DROP        , GLA_LEAF_DROP       , 
     I    TMP_CNP_NEXT         , TMP_SKIN_NEXT  )  
!
      ENDIF
#endif
!
      RETURN
      END SUBROUTINE CANOPY_INTEGRATION_OVERIMP_RUN  
!================================================================
!
! 仮積分時に、見積もられた温度等変化分を加算する。
! ただしゼロ度を跨ぐときは、相変化の有無に依らず、ゼロ度にする。
!
      SUBROUTINE CANOPY_INTEGRATION_OVERIMP_ITR (
     I   NSTR       , NEND        , 
     I   DELT_CNP   , DELT_ATM    , 
     I   SMTX       , 
     I   FAC_EVP_2  , 
     I   DTA_DTC    , DTA_DTG     , DTA_DTR    , 
     M   TMP_CNP    , TMP_SKIN    , TMP_SS1    , 
     M   TMP_REF    , TMP_AIR     , E_REF      )
!
      USE SIBPRM , ONLY :
     1    IDX  , ICG , ICGB  , IEST , ICGR  
!
      USE SIBCON , ONLY :
     1   TMP_FREZ_C
!
      IMPLICIT NONE 
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
!
      REAL(8),INTENT(IN)   :: SMTX       (IDX,5,6)
!
      REAL(8),INTENT(IN)   :: DELT_CNP
      REAL(8),INTENT(IN)   :: DELT_ATM
!
      REAL(8),INTENT(IN)   :: DTA_DTC (IDX) 
      REAL(8),INTENT(IN)   :: DTA_DTG (IDX) 
      REAL(8),INTENT(IN)   :: DTA_DTR (IDX) 
!
      REAL(8),INTENT(IN)   :: FAC_EVP_2(IDX) 
!
      REAL(8),INTENT(INOUT):: TMP_CNP  (IDX)
      REAL(8),INTENT(INOUT):: TMP_AIR  (IDX)
      REAL(8),INTENT(INOUT):: TMP_REF  (IDX)
      REAL(8),INTENT(INOUT):: TMP_SKIN (IDX)
      REAL(8),INTENT(INOUT):: TMP_SS1  (IDX)
      REAL(8),INTENT(INOUT):: E_REF    (IDX)
!
! WORK
!
      REAL(8)  :: YDTMP_CNP
      REAL(8)  :: YDTMP_SKIN
      REAL(8)  :: YDTMP_SS1
      REAL(8)  :: YDTMP_AIR 
      REAL(8)  :: YDTMP_REF
      REAL(8)  :: YDQ_REF
!
      REAL(8)  :: TMP_CNP_SAVE(IDX)
      REAL(8)  :: TMP_SKIN_SAVE(IDX)
      REAL(8)  :: TMP_SS1_SAVE(IDX)
!
      REAL(8)  :: DELT_A
      REAL(8)  :: DELT_C
      INTEGER  :: I
!
      LOGICAL,SAVE :: LFIRST = .TRUE. 
!
      integer,SAVE :: icount = 0
!
#ifdef SIB_DEBUG
      integer :: IMSG_NUM 
      real(8) :: XMSG_SQ  
#endif

      icount = icount+1
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP_ITR 00/06 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DELT_A = DELT_ATM 
      DELT_C = DELT_CNP 
!
! ================
! >> 温度の積分 <<
! ================
!
      DO I=NSTR,NEND
!
! ------------------
! > 仮の温度変化量 <
! ------------------
!
        YDTMP_CNP     =   SMTX(I,1,6) * DELT_C  
!
        YDTMP_SKIN    =   SMTX(I,2,6) * DELT_C  
!
        YDTMP_SS1     =   SMTX(I,3,6) * DELT_C  
!
!       YDTMP_REF     =   SMTX(I,5,6) * DELT_A  
        YDTMP_REF     =   SMTX(I,4,6) * DELT_A  
!
        YDQ_REF       =   SMTX(I,5,6) * DELT_A 
!
        YDTMP_AIR     =    DTA_DTC(I) * YDTMP_CNP 
     1                   + DTA_DTG(I) * YDTMP_SKIN
     1                   + DTA_DTR(I) * YDTMP_REF 
!
! ----------------------
! > 推定平衡温度の算出 <
! ----------------------
!    ただし零度を跨ぐときは零度にする。（イテレーション中の場合）
!
!    注意。
!      TMP_AIR はイテレーションループのはじめの方で計算されるが、
!      その初期推定値として一応変更。
!
!      ・水蒸気圧 e ( 単位 HPa ) の評価には、
!             e ～ pq/ε という近似式を用いている。fac_evp_2=Ps/ε
!
        TMP_CNP_SAVE(I)  = TMP_CNP (I)
        TMP_SKIN_SAVE(I) = TMP_SKIN(I)
        TMP_SS1_SAVE(I)  = TMP_SS1 (I)
        TMP_CNP (I)      = TMP_CNP (I) + YDTMP_CNP 
        TMP_SKIN(I)      = TMP_SKIN(I) + YDTMP_SKIN
        TMP_SS1 (I)      = TMP_SS1 (I) + YDTMP_SS1
        TMP_REF (I)      = TMP_REF (I) + YDTMP_REF 
#ifdef CORRECT
        TMP_AIR (I)      = TMP_AIR (I) + YDTMP_AIR 
#else
        TMP_AIR (I)      = TMP_AIR (I) + 0 
#endif
!
        E_REF   (I)      = E_REF(I)    + YDQ_REF * FAC_EVP_2(I)  
      ENDDO
!
      DO I=NSTR,NEND 
!
        IF ( ( TMP_CNP(I)-TMP_FREZ_C )
     1        * ( TMP_CNP_SAVE(I) -TMP_FREZ_C ) .LE. 0. ) THEN
#ifdef CORRECT
          TMP_AIR (I) =  TMP_AIR(I) 
     1                    + DTA_DTC(I) * ( TMP_FREZ_C - TMP_CNP(I) ) 
#endif
          TMP_CNP (I) = TMP_FREZ_C 
        ENDIF
!
        IF ( ( TMP_SKIN(I)-TMP_FREZ_C )
     1        * ( TMP_SKIN_SAVE(I)-TMP_FREZ_C ) .LE. 0. ) THEN
#ifdef CORRECT
          TMP_AIR (I) =  TMP_AIR(I) 
     1                    + DTA_DTG(I) * ( TMP_FREZ_C - TMP_SKIN(I) ) 
#endif
          TMP_SKIN(I) = TMP_FREZ_C 
        ENDIF
!
        IF ( ( TMP_SS1(I)-TMP_FREZ_C )
     1        * ( TMP_SS1_SAVE(I)-TMP_FREZ_C ) .LE. 0. ) THEN
          TMP_SS1(I) = TMP_FREZ_C 
        ENDIF
      ENDDO
!
! ================
! > 値のチェック <
! ================
!
#ifdef SIB_DEBUG
!
#include "noparallel_fixedform"
      DO I=NSTR,NEND 
        IF (      ABS(SMTX(I,1,6)*DELT_C) .GT. 50.
     1       .OR. ABS(SMTX(I,2,6)*DELT_C) .GT. 50. ) THEN
          write(6,*) 'CANOPY_INTEGRATION_OVERIMP_ITR ERROR B1 ' , 
     1             'I ' , I , ' NSTR NEND' , NSTR , NEND , 
     1             'DTC ' , SMTX(I,1,6)*DELT_C,
     1             'DTG ' , SMTX(I,2,6)*DELT_C
          stop 999
        ENDIF
      ENDDO
!
      IMSG_NUM = 0
      XMSG_SQ  = 0
!
      DO I=NSTR,NEND 
        IF 
     1     (      ABS(SMTX(I,1,6)*DELT_C) .GT. 10.
     1       .OR. ABS(SMTX(I,2,6)*DELT_C) .GT. 10. ) THEN
          IMSG_NUM = IMSG_NUM+1 
          XMSG_SQ  = XMSG_SQ  + ( SMTX(I,2,6)*DELT_C -10 )**2 
     1                        + ( SMTX(I,1,6)*DELT_C -10 )**2 
!         write(6,*) 'CANOPY_INTEGRATION_OVERIMP_ITR WARINIG B1 ' , 
!    1             'I ' , I , ' NSTR NEND' , NSTR , NEND 
!         write(6,*)  
!    1             'TCO TCN DTC ' , TMP_CNP_SAVE(I) , 
!    1                              TMP_CNP (I) , SMTX(I,1,6)*DELT_C
!         write(6,*)  
!    1             'TGO TGN DTG ' , TMP_SKIN_SAVE(I) ,  
!    1                              TMP_SKIN(I) , SMTX(I,2,6)*DELT_C
        ENDIF
      ENDDO
!
      if ( imsg_num.gt.0 ) then
        write(6,*) 'canopy_integration_overimp_itr mod more than 10K' , 
     1           'num_point/sum_sq' , IMSG_NUM , XMSG_SQ 
      endif
#endif
!
      RETURN
      END SUBROUTINE CANOPY_INTEGRATION_OVERIMP_ITR 
!=======================================================================
!
! 本積分時（仮積分ではない）に、各種フラックスを計算する。
!    .. インプリシット補正項を加える。
!
      SUBROUTINE CANOPY_INTEGRATION_OVERIMP_PRE (
     I   NSTR       , NEND        , INSFLAG    , 
     I   DELT_CNP   , DELT_ATM    , 
     I   SMTX       , 
     I   EVAP       , SUBL        , TRNS       , 
     I   SENS       , 
     I   RAD_LONG   , 
!    I   RAD_SHRT   , 
     I   HDF        , 
     I   DEVAP_DX   , DSENS_DX    , DRAD_DX    , HDF_SKIN  , 
     O   EVAP_TOT   , SUBL_TOT    , TRNS_TMP   , 
     O   RAD_LONG_TMP  , 
     O   SENS_TMP    , HDF_TMP    ,
     O   DTMP_CNP   , DTMP_SKIN   , 
     M   RAD_LONG_NorS_2_ATM )
!
      USE SIBPRM , ONLY :
     1    IDX  , ICG , ICGB  , IEST , ICGR  
!
      USE SIBCON , ONLY :
     1    DELT_FAC_ATM_c , 
     1    DELT_FAC_CNP_c ,
     1    HL_WTR_c
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, ICN, ICS,
     1    IMONIT_LEVEL ,
     1    IMONIT_ALL   
#endif
      !
!     USE MESSAGE
!
      IMPLICIT NONE 
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: INSFLAG 
!
!     REAL(8),INTENT(IN)   :: SMTX       (IDX,4,5)
      REAL(8),INTENT(IN)   :: SMTX       (IDX,5,6)
!
      REAL(8),INTENT(IN)   :: DELT_CNP
      REAL(8),INTENT(IN)   :: DELT_ATM
!
      REAL(8),INTENT(IN)   :: EVAP       (IDX,ICGB)
      REAL(8),INTENT(IN)   :: SUBL       (IDX,ICGB)
      REAL(8),INTENT(IN)   :: TRNS       (IDX,ICG)
      REAL(8),INTENT(IN)   :: DEVAP_DX   (IDX,IEST,ICGB,ICGR)
!
      REAL(8),INTENT(IN)   :: SENS       (IDX,ICG)
      REAL(8),INTENT(IN)   :: DSENS_DX   (IDX,ICG,ICGR) 
!
      REAL(8),INTENT(IN)   :: RAD_LONG   (IDX,ICG)
!     REAL(8),INTENT(IN)   :: RAD_SHRT   (IDX,ICG)
      REAL(8),INTENT(IN)   :: DRAD_DX    (IDX,ICGR,ICG) ! 放射吸収量変化率
!
      REAL(8),INTENT(IN)   :: HDF        (IDX,2) 
      REAL(8),INTENT(IN)   :: HDF_SKIN   (IDX,2) 
!
      REAL(8),INTENT(OUT)  :: EVAP_TOT (IDX,ICGB) 
      REAL(8),INTENT(OUT)  :: SUBL_TOT (IDX,ICGB) 
      REAL(8),INTENT(OUT)  :: TRNS_TMP (IDX,ICG ) 
      REAL(8),INTENT(OUT)  :: RAD_LONG_TMP  (IDX,ICG ) 
      REAL(8),INTENT(OUT)  :: SENS_TMP (IDX,ICG ) 
      REAL(8),INTENT(OUT)  :: HDF_TMP  (IDX,2) 
!
      REAL(8),INTENT(OUT)  :: DTMP_CNP (IDX)
      REAL(8),INTENT(OUT)  :: DTMP_SKIN(IDX)
!
      REAL(8),INTENT(INOUT):: RAD_LONG_NorS_2_ATM(IDX)
!
      REAL(8)              :: XDTMP_CNP (IDX)
      REAL(8)              :: XDTMP_SKIN(IDX)
      REAL(8)              :: XDTMP_SS1 (IDX)
      REAL(8)              :: XDTMP_REF (IDX)
!
      REAL(8)              :: XDHUM_REF (IDX)
!
      REAL(8)              :: BAD_WTR (IDX,ICG)
!                                                
! WORK 
!
      INTEGER   :: I    
      REAL(8)   :: DELT_A
      REAL(8)   :: DELT_C
!     REAL(8),PARAMETER  :: DZERO = 0.D0 
!
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      CHARACTER(7) :: CMARK 
      real(8)      :: work(idx) 
#endif
!
#ifdef SIB_DEBUG
      integer :: IMSG_NUM 
      real(8) :: XMSG_SQ  
#endif

!
      LOGICAL,SAVE       :: LFIRST = .TRUE.  
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP_PRE 00/06 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DELT_A = DELT_ATM * DELT_FAC_ATM_c
      DELT_C = DELT_CNP * DELT_FAC_CNP_C
!
#ifdef SIB_DEBUG
      write(6,*) 'canopy_integration_overimp_pre ' , 
     1  ' DELT_A  DELT_ATM  DELT_FAC_ATM DELT_C DELT_CNP DELT_FAC_CNP ',
     1  DELT_A , DELT_ATM , DELT_FAC_ATM_c ,
     1  DELT_C , DELT_CNP , DELT_FAC_CNP_c
#endif
!
! ================
! >> 温度変化量 <<
! ================
!   方程式の解から得られる温度・比湿変化見積もり量。
!
      DO I=NSTR,NEND
        DTMP_CNP  (I) =   SMTX(I,1,6) * DELT_CNP  
        DTMP_SKIN (I) =   SMTX(I,2,6) * DELT_CNP  
        XDTMP_CNP (I) =   SMTX(I,1,6) * DELT_C  
        XDTMP_SKIN(I) =   SMTX(I,2,6) * DELT_C  
        XDTMP_SS1 (I) =   SMTX(I,3,6) * DELT_C  
        XDTMP_REF (I) =   SMTX(I,4,6) * DELT_A  
!
        XDHUM_REF (I) =   SMTX(I,5,6) * DELT_A  
      ENDDO
!
! ==============
! >> チェック <<
! ==============
!
! ----------------------------
! > 変化量の大きすぎチェック <
! ----------------------------
!  温度変化量が大きすぎないか。目安として 30K にしているが、意味はない。
!  甘すぎる気がとてもする。
!
#ifdef SIB_DEBUG
!
#include "noparallel_fixedform"
      DO I=NSTR,NEND 
        IF (      ABS(DTMP_CNP(I) ) .GT. 40.              ! 000926 で 30->40 
     1       .OR. ABS(DTMP_SKIN(I)) .GT. 40. ) THEN
          write(6,*) 'CANOPY_INTEGRATION_OVERIMP_PRE ERROR B1 ' , 
     1             'I ' , I , ' INSFLAG ' , INSFLAG , 
     1               DTMP_CNP(I) , DTMP_SKIN(I)  
          write(6,*) 'dtmp ookisugi. suguochiru hazu'
          stop 999
        ENDIF
      ENDDO
!
      IMSG_NUM = 0
      XMSG_SQ  = 0
!
      DO I=NSTR,NEND 
        IF (      ABS(DTMP_CNP(I) ) .GT. 10.       
     1       .OR. ABS(DTMP_SKIN(I)) .GT. 10. ) THEN
          IMSG_NUM = IMSG_NUM + 1  
          XMSG_SQ  = XMSG_SQ  + (DTMP_CNP(I) -10.)**2       
     1                        + (DTMP_SKIN(I)-10.)**2       
!         write(6,*) 'CANOPY_INTEGRATION_OVERIMP_PRE ERROR B1 ' , 
!    1             'I ' , I , ' INSFLAG ' , INSFLAG , 
!    1               DTMP_CNP(I) , DTMP_SKIN(I)  
!         write(6,*) 'dtmp ookisugi. suguochiru hazu'
        ENDIF
      ENDDO
!
      if ( imsg_num .gt. 0 ) then
        write(6,*) 'canopy_integration_overimp_pre mod more than 10K' , 
     1           'num_point/sum_sq' , IMSG_NUM , XMSG_SQ 
      endif
#endif
!
! ========================
! >> 各種フラックス計算 <<
! ========================
!
!   まずエラー補正をせずに見積もる。残差はあるはず。
!
! 蒸発 （下草・裸値、雪はそれぞれ合計。）
!   フラックス、下向き正
!       EVAP_TMP(IDX,ICGB) EVAP(IDX,ICGB), DEVAP_DX(IDX,IEST,ICGB,ICGR)
! 
      DO I=NSTR,NEND
        EVAP_TOT (I,1) 
     1       =         EVAP  (I,1) 
     1              +  DEVAP_DX(I,1,1,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,1,1,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,1,1,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,1,1,3) * XDHUM_REF (I)
!
        EVAP_TOT (I,2) 
     1        =        EVAP  (I,2) 
     1              +  DEVAP_DX(I,1,2,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,1,2,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,1,2,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,1,2,3) * XDHUM_REF (I)
!
        EVAP_TOT (I,3) 
     1        =        EVAP  (I,3) 
     1              +  DEVAP_DX(I,1,3,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,1,3,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,1,3,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,1,3,3) * XDHUM_REF (I)
!
! 昇華 
!       SUBL_TOT(IDX,ICGB) SUBL(IDX,ICGB),
!       DEVAP_DX(IDX,IEST,ICGB,ICGR)
!
        SUBL_TOT (I,1)  
     1        =        SUBL  (I,1) 
     1              +  DEVAP_DX(I,2,1,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,2,1,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,2,1,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,2,1,3) * XDHUM_REF (I)
!
        SUBL_TOT (I,2) 
     1        =        SUBL  (I,2) 
     1              +  DEVAP_DX(I,2,2,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,2,2,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,2,2,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,2,2,3) * XDHUM_REF (I)
!
        SUBL_TOT (I,3) 
     1        =        SUBL  (I,3) 
     1              +  DEVAP_DX(I,2,3,1) * XDTMP_CNP (I)
     1              +  DEVAP_DX(I,2,3,2) * XDTMP_SKIN(I)
!    1              +  DEVAP_DX(I,2,3,3) * XDTMP_REF (I)
     1              +  DEVAP_DX(I,2,3,3) * XDHUM_REF (I)
!
! 蒸散、とりあえず全層からの和 
!       TRNS_TMP(IDX,ICG) TRNS(IDX,ICG), DEVAP_DX(IDX,IEST,ICGB,ICGR)
!       裸地からの蒸散はない
!
        TRNS_TMP (I,1)  
     1        =       TRNS  (I,1) 
     1              + DEVAP_DX(I,3,1,1) * XDTMP_CNP (I)
     1              + DEVAP_DX(I,3,1,2) * XDTMP_SKIN(I)
!    1              + DEVAP_DX(I,3,1,3) * XDTMP_REF (I)
     1              + DEVAP_DX(I,3,1,3) * XDHUM_REF (I)
!
        TRNS_TMP (I,2) 
     1        =       TRNS  (I,2) 
     1              + DEVAP_DX(I,3,2,1) * XDTMP_CNP (I)
     1              + DEVAP_DX(I,3,2,2) * XDTMP_SKIN(I)
!    1              + DEVAP_DX(I,3,2,3) * XDTMP_REF (I)
     1              + DEVAP_DX(I,3,2,3) * XDHUM_REF (I)
!
! 顕熱計算
!       フラックス、下向き正
!       SENS_TMP (IDX,ICG) , SENS (IDX,ICG), DSENS_DX (IDX,ICG,ICGR) 
!
        SENS_TMP (I,1)  
     1        =       SENS (I,1) 
     1              + DSENS_DX(I,1,1) * XDTMP_CNP (I)
     1              + DSENS_DX(I,1,2) * XDTMP_SKIN(I)
     1              + DSENS_DX(I,1,3) * XDTMP_REF (I)
!
        SENS_TMP (I,2)  
     1        =       SENS (I,2) 
     1              + DSENS_DX(I,2,1) * XDTMP_CNP (I)
     1              + DSENS_DX(I,2,2) * XDTMP_SKIN(I)
     1              + DSENS_DX(I,2,3) * XDTMP_REF (I)
!
! 放射収束計算  ... 受け取る向きが正。
!                   RAD, DRAD_DX は受け取る向きが正であることに注意。
!       RAD  (IDX,ICG) , DRAD_DX_TMP (IDX,ICG,ICG) 
!

        RAD_LONG_TMP (I,1)  
     1       =        RAD_LONG (I,1) 
     1              - DRAD_DX(I,2,1) * XDTMP_CNP (I)
     1              + DRAD_DX(I,3,1) * XDTMP_CNP (I)
     1              - DRAD_DX(I,1,2) * XDTMP_SKIN(I)
!
        RAD_LONG_TMP (I,2)  
     1       =        RAD_LONG (I,2)
     1              + DRAD_DX(I,2,1) * XDTMP_CNP (I)
     1              + DRAD_DX(I,1,2) * XDTMP_SKIN(I)
     1              + DRAD_DX(I,3,2) * XDTMP_SKIN(I)

!      大気へ戻す長波放射の補正
!
        RAD_LONG_NorS_2_ATM (I)  
     1       = RAD_LONG_NorS_2_ATM (I)  
     1         + DRAD_DX(I,3,1) * XDTMP_CNP (I)
     1         + DRAD_DX(I,3,2) * XDTMP_SKIN(I)
!
! 下面とのやりとり ... 下向き = 失う向き正
!
        HDF_TMP (I,1) 
     1      =         HDF (I,1) 
     1              + HDF_SKIN(I,1) * XDTMP_SKIN (I)
     1              - HDF_SKIN(I,1) * XDTMP_SS1  (I)
!
        HDF_TMP (I,2) 
     1      =         HDF (I,2) 
     1              + HDF_SKIN(I,2) * XDTMP_SS1  (I)
!
      ENDDO
!
! ========================
! >> 下向きの蒸散の処理 <<
! ========================
!
! 下向きの蒸散はあり得ない。エネルギー的に大罪なので、何かで面倒を見たい。
!
!   顕熱で面倒を見る。
!   ただしモニタにはきちんと取っておく
      BAD_WTR(:,:) = 0.D0
      DO I=NSTR,NEND
        IF ( TRNS_TMP(I,1) .GT. 0.D0 ) THEN
          SENS_TMP (I,1) = SENS_TMP(I,1) + TRNS_TMP (I,1) * HL_WTR_c
          BAD_WTR  (I,1) = TRNS_TMP(I,1)
          TRNS_TMP (I,1) = 0.D0 
        ENDIF
        IF ( TRNS_TMP(I,2) .GT. 0.D0 ) THEN
          SENS_TMP (I,2) = SENS_TMP(I,2) + TRNS_TMP (I,2) * HL_WTR_c
          BAD_WTR  (I,2) = TRNS_TMP(I,2)
          TRNS_TMP (I,2) = 0.D0 
        ENDIF
      ENDDO
!
#ifdef AAAAAA
!   蒸発（結露）で面倒を見る。方法 2 
!   ただしモニタにはきちんと取っておく
      BAD_WTR(:,:) = 0.D0
      DO I=NSTR,NEND
        IF ( TRNS_TMP(I,1) .GT. 0.D0 ) THEN
          EVAP_TOT (I,1) = EVAP_TOT(I,1) + TRNS_TMP (I,1)
          BAD_WTR  (I,1) = TRNS_TMP(I,1)
          TRNS_TMP (I,1) = 0.D0 
        ENDIF
        IF ( TRNS_TMP(I,2) .GT. 0.D0 ) THEN
          EVAP_TOT (I,2) = EVAP_TOT(I,2) + TRNS_TMP (I,2)
          BAD_WTR  (I,2) = TRNS_TMP(I,2)
          TRNS_TMP (I,2) = 0.D0 
        ENDIF
      ENDDO
#endif
!
!  --------------
!  > 詳細モニタ <
!  --------------
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL  .GE. IMONIT_ALL ) THEN
!
!  負の蒸散を結露にした量
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'WMVNT2E' 
          CALL MONIT_SIB ( CMARK, BAD_WTR(1,1), ICN, DELT_CNP, ONE)       
          CMARK = 'WMVGT2E' 
          CALL MONIT_SIB ( CMARK, BAD_WTR(1,2), ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'WMVST2E' 
          CALL MONIT_SIB ( CMARK, BAD_WTR(1,1), ICS, DELT_CNP, ONE)       
        ENDIF  
!
!   放射吸収の修正量
!       向きは、下向き正
!
!
!     キャノピー -> 大気 
        DO I=NSTR,NEND 
!         WORK(I) = DRAD_DX_TMP(I,1,1) * XDTMP_CNP (I) / 2 
!                   （ TC 変化で C が得る分の半分 ）  
          WORK(I) = DRAD_DX    (I,3,1) * XDTMP_CNP (I) 
        ENDDO
        IF ( INSFLAG.EQ.0 ) THEN 
          CMARK = 'HLRNAD' 
          CALL MONIT_SIB ( CMARK, WORK, ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'HLRSAD' 
          CALL MONIT_SIB ( CMARK, WORK, ICS, DELT_CNP, ONE)       
        ENDIF
!
!     キャノピー -> 地面 
        DO I=NSTR,NEND 
!                          （ TC 変化で G が得る分 
!                           ちなみに = - TC 変化でCが得る分の半分）
!         WORK(I) = DRAD_DX_TMP(I,2,1) * XDTMP_CNP (I) 
          WORK(I) = DRAD_DX    (I,2,1) * XDTMP_CNP (I) 
        ENDDO
        IF ( INSFLAG.EQ.0 ) THEN 
          CMARK = 'HLRNGD' 
          CALL MONIT_SIB ( CMARK, WORK, ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'HLRSWD' 
          CALL MONIT_SIB ( CMARK, WORK, ICS, DELT_CNP, ONE)       
        ENDIF
!
!     地面       -> 大気 
        DO I=NSTR,NEND 
!                      （ TG 変化で G が得る分 - TG 変化で C が得る分 ）  
!         WORK(I) =   DRAD_DX_TMP(I,2,2) * XDTMP_SKIN(I)
!    1              - DRAD_DX_TMP(I,1,2) * XDTMP_SKIN(I)
          WORK(I) =   DRAD_DX    (I,3,2) * XDTMP_SKIN(I)
        ENDDO
        IF ( INSFLAG.EQ.0 ) THEN 
          CMARK = 'HLRGAD' 
          CALL MONIT_SIB ( CMARK, WORK, ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'HLRWAD' 
          CALL MONIT_SIB ( CMARK, WORK, ICS, DELT_CNP, ONE)       
        ENDIF

!     地面       -> キャノピー 
        DO I=NSTR,NEND 
!         （ - TG 変化で C が得る分 ）  
!         WORK(I) = - DRAD_DX_TMP(I,1,2) * XDTMP_SKIN(I)
          WORK(I) =   DRAD_DX    (I,1,2) * XDTMP_SKIN(I)
        ENDDO
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'HLRGND' 
          CALL MONIT_SIB ( CMARK, WORK, ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'HLRWSD' 
          CALL MONIT_SIB ( CMARK, WORK, ICS, DELT_CNP, ONE)       
        ENDIF
!
      ENDIF
#endif
!
      RETURN
      END SUBROUTINE  CANOPY_INTEGRATION_OVERIMP_PRE 
!================================================================
!  植生の積分
!    手順概要
!    ・水を移動。水が足りないときは蒸発・昇華を減らして顕熱へ。
!    ・温度積分
!    ・ゼロ度を跨ぐ変化のときを意識した補正。
!
!    ・収支チェック （お粗末なバグ発見用）
!
!   注意
!    ・系にもともとあった H2O は、DTMP+... だけ昇温する。
!      蒸発などする水についても、昇温して、相変化し、系から潜熱を奪っ
!      て出ていくと考える。
!    ・大気内にあった H2O が葉にくる場合、昇温後の温度をもってくると
!      考える。
!      葉・地面に結露する水は、系の昇温後の温度をもった水蒸気が、葉・
!      地面の上で相変化し、潜熱を与えるとする。葉に入る水はそのまま
!      系の水となり、地面への水はバケツに入る。
!    ・土壌中の H2O が蒸発・昇華する場合、土壌第一層の温度をもって地
!      面表面に達し、地面スキンから潜熱を奪って土壌の温度の水蒸気と
!      なって出ていくと考える。
!      （canopy 計算用の仮の土壌第一層温度ではなく、真の土壌第一層温度）
!      土壌中から蒸散する水も同様で、土壌の温度のままで水蒸気になり、
!      潜熱を葉から奪うと考える。
!    ・以上からわかるように、H2O の面倒な温度変化は全て大気に押しつけ
!      ている。かつ、この処理は大気の鉛直拡散では現時点(00/09/13) で考
!      慮されていないので、陸面 + 大気鉛直拡散 の合同系ではエネルギーは
!      保存していない。これは今後の解決しなければならない課題である。
!      ( もし難しいとしても、どの程度の量なのか、見積りが欲しい )
!
      SUBROUTINE CANOPY_INTEGRATION_OVERIMP_VEG (
     I   NSTR         , NEND         , ICG_FLAG     , CVR_S_OR_N     , 
     I   INSFLAG      , IMSKLT       , 
     I   TMP_OLD      , DTMP         , 
     I   FR_H2O_OLD   , FR_WBH_OLD   , 
     I   TMP_SOIL     , H2O_SAT_LEAF , DELT_CNP     , XLAI_VEG_msx   , 
     M   RAD_LONG     , 
     I   RAD_SHRT     , 
     M   SENS         , HDF          , 
     I   DRAD_DX      , 
     I   DSENS_DX     , HDF_GRSK     , 
     I   TRNS         ,  
     M   EVAP_LEAF    , SUBL_LEAF    , 
     I   EVAP_GRSK    , SUBL_GRSK    , 
     O   TMP_NEW      , FR_WTR_NEW   , FR_GLA_NEW   , 
     O   WTR_LEAF_DROP, GLA_LEAF_DROP,  
     M   RAD_LONG_NorS_2_ATM,   
     M   ENG_BUCKET                  , 
     M   H2O_BUCKET_INTO             , 
     M   H2O_BUCKET_FROM             , 
     M   ENH_BUCKET_INTO             , 
     M   ENH_BUCKET_FROM             ) 
! 
      USE SIBPRM , ONLY : IDX 
!
      USE SIBCON , ONLY : 
     1      TMP_FREZ_C ,  
     1      CAP_WTR_C  ,  
     1      CAP_GLA_C  ,  
     1      CAP_GRSK_C ,  
     1      CAP_LAI_C  ,  
     1      HL_W2G_C   , 
     1      HL_WTR_C   , 
     1      HL_GLA_C    
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN, ICS,
     1      IMONIT_LEVEL ,
     1      IMONIT_ALL 
#endif
!     USE MESSAGE
!
      IMPLICIT NONE
!
! IN
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: ICG_FLAG          ! C->1 G->2
      INTEGER,INTENT(IN)    :: INSFLAG           ! N->0 S->1 
      INTEGER,INTENT(IN)    :: IMSKLT(IDX) 
      REAL(8),INTENT(IN)    :: DELT_CNP           
      REAL(8),INTENT(IN)    :: DTMP     (IDX)    ! 評価された昇温量。
                                                 ! 積分前温度に対して。
      REAL(8),INTENT(IN)    :: TMP_OLD  (IDX)    ! 積分前の温度
      REAL(8),INTENT(IN)    :: TMP_SOIL (IDX)    ! 土壌第一層の温度(蒸発水温度)
      REAL(8),INTENT(IN)    :: CVR_S_OR_N  (IDX) ! 積雪被覆率 or 1- 
      REAL(8),INTENT(IN)    :: H2O_SAT_LEAF(IDX) ! 葉の飽和時の水量
      REAL(8),INTENT(IN)    :: XLAI_VEG_msx(IDX) ! LAI 
      REAL(8),INTENT(IN)    :: FR_H2O_OLD (IDX) ! 積分前の葉のもつ水量
      REAL(8),INTENT(IN)    :: FR_WBH_OLD (IDX) ! 積分前の葉のもつ氷量
!
! OUT
!
      REAL(8),INTENT(OUT)   :: TMP_NEW    (IDX)  ! 積分後の温度
      REAL(8),INTENT(OUT)   :: FR_WTR_NEW (IDX)  ! 積分後の水の割合
      REAL(8),INTENT(OUT)   :: FR_GLA_NEW (IDX)  ! 積分後の氷の割合
      REAL(8),INTENT(OUT)   :: WTR_LEAF_DROP (IDX) ! 葉から落ちた水の量
      REAL(8),INTENT(OUT)   :: GLA_LEAF_DROP (IDX) ! 葉から落ちた氷の量
!
! INOUT 
!
      REAL(8),INTENT(INOUT) :: RAD_LONG (IDX)    ! 放射
      REAL(8),INTENT(IN)    :: RAD_SHRT (IDX)    ! 放射
      REAL(8),INTENT(INOUT) :: SENS     (IDX)    ! 顕熱 
      REAL(8),INTENT(INOUT) :: HDF      (IDX)    ! 地熱
      REAL(8),INTENT(IN)    :: TRNS     (IDX)    ! 蒸散
      REAL(8),INTENT(IN)    :: DRAD_DX  (IDX)    ! 
      REAL(8),INTENT(IN)    :: DSENS_DX (IDX)    ! 顕熱の温度依存性
      REAL(8),INTENT(IN)    :: HDF_GRSK (IDX,2)  ! 地熱の温度依存性
      REAL(8),INTENT(INOUT) :: EVAP_LEAF(IDX)    ! 葉からの蒸発
      REAL(8),INTENT(IN)    :: EVAP_GRSK(IDX)    ! 裸地からの蒸発
      REAL(8),INTENT(INOUT) :: SUBL_LEAF(IDX)    ! 葉からの昇華
      REAL(8),INTENT(IN)    :: SUBL_GRSK(IDX)    ! 裸地からの昇華
!
      REAL(8),INTENT(INOUT) :: RAD_LONG_NorS_2_ATM(IDX)   
!
! INOUT 各種バケツ
!
      REAL(8),INTENT(INOUT):: ENG_BUCKET            (IDX)
      REAL(8),INTENT(INOUT):: H2O_BUCKET_INTO       (IDX)
      REAL(8),INTENT(INOUT):: H2O_BUCKET_FROM       (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_INTO       (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_FROM       (IDX)
!
! WORK
!
      INTEGER      :: I
      REAL(8)      :: CAP_SYSTEM_OLD  (IDX)     ! 系の比熱
!
      REAL(8)      :: GLA_LEAF_OLD  (IDX)  
      REAL(8)      :: WTR_LEAF_OLD  (IDX) 
      REAL(8)      :: DGLA_LEAF_MOD (IDX)  
      REAL(8)      :: DWTR_LEAF_MOD (IDX) 
      REAL(8)      :: GLA_LEAF_NEW  (IDX)     ! 積分後の葉の持つ氷
      REAL(8)      :: WTR_LEAF_NEW  (IDX)     ! 積分後の葉の持つ水 

      INTEGER      :: IOVER_TMP_ZERO(IDX) 
!
      REAL(8)      :: ENG_PHASE_ALL (IDX) 
      REAL(8)      :: ENG_PHASE     (IDX) 
      REAL(8)      :: ENG_STORE_1   (IDX)     ! 水の不足で溜る熱
      REAL(8)      :: ENG_STORE_2   (IDX)     ! ゼロ度越えがなくなることで
!                                             ! フラックスがかわって溜る熱
      REAL(8)      :: ENG_STORE_3   (IDX)     ! ゼロ度越えがなくなることで
!                                             ! 内部エネルギーで溜る熱
      REAL(8)      :: TMP_MOD       (IDX)     ! ゼロ度越えがなくなることで
!                                             ! 変わる温度上昇量
      REAL(8)      :: SENS_MOD      (IDX)     
      REAL(8)      :: RAD_LONG_MOD  (IDX)     
      REAL(8)      :: HDF_MOD       (IDX)     
!
!     REAL(8)      :: RAD_LONG_MOD_C2A (IDX)     
!     REAL(8)      :: RAD_LONG_MOD_G2A (IDX)     
!     REAL(8)      :: RAD_LONG_MOD_C2G (IDX)     
!     REAL(8)      :: RAD_LONG_MOD_G2C (IDX)     
!
      REAL(8)      :: X 
      REAL(8)      :: FRAC
      REAL(8)      :: FRAC_DROP
      REAL(8)      :: FRAC_STAY
!
      REAL(8)      :: CAP_SYSTEM_TMP(IDX) 
      REAL(8)      :: ENG_OLD     (IDX)       
      REAL(8)      :: ENG_NEW     (IDX)       
      REAL(8)      :: ENG_ADD     (IDX)       
!     REAL(8)      :: ENG_ADD_2   (IDX)       
      REAL(8)      :: CAP_SYSTEM_AFTER(IDX)  ! 葉の比熱（最後）
!
      REAL(8),SAVE :: TMP_MIN = 1.D-10     ! TMP_FREZ_C を越えた変化をするとき
      REAL(8),SAVE :: ZERO    = 0.D0
!
      INTEGER,PARAMETER :: III_040120 = 0
!
#ifdef MONYOS
!     REAL(8),PARAMETER :: ONE_MINUS = -1.D0 
      REAL(8),PARAMETER :: ONE = 1.D0 
      CHARACTER(7) :: CMARK 
      REAL(8) :: WORK_IDX(IDX)
#endif
!
#ifdef SIB_DEBUG
      integer :: IMSG_NUM 
      real(8) :: XMSG_SQ  
#endif
      real(8) :: eng , enx , cap
!
      LOGICAL,SAVE :: LFIRST = .TRUE. 
!
      integer,SAVE :: icount = 0
      icount = icount+1
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP_VEG 00/06 POCHAKA'
        IF (      TMP_FREZ_C .EQ. TMP_FREZ_C+TMP_MIN  
     1       .OR. TMP_FREZ_C .EQ. TMP_FREZ_C-TMP_MIN ) THEN
          WRITE(6,*) 'ERROR CANOPY_INTEGRATION_OVERIMP_VEG INI' 
          STOP 999
        ENDIF   
        LFIRST = .FALSE. 
      ENDIF
!
      TMP_NEW    (:) = 0.   
      FR_WTR_NEW (:) = 0.
      FR_GLA_NEW (:) = 0.
!
! ----------------
! > 値のチェック <
! ----------------
!  昇温量が 30 度を越えていたら、止める。
!  30 度という基準はいい加減。もっと厳しくてもいいと思う。(00/09/13 pochaka) 
!
#ifdef SIB_DEBUG
!
#include "noparallel_fixedform"
      DO I=NSTR,NEND 
        IF ( ABS(DTMP(I)) .GT. 40. ) THEN   ! 000926 30->40
          write(6,*) 'CANOPY_INTEGRATION_OVERIMP_VEG ERROR B1 ' ,  
     1             'I ' , I , 'IMSKLT ' , IMSKLT(I) , 
     1             ' INSFLAG ' , INSFLAG , 
     1             'ICG_FLAG '    , ICG_FLAG     , 
     1             ' CVR_S_OR_N ' , CVR_S_OR_N(I),
     1             'DTMP TMP_OLD TMP_OLD+DTMP' ,    
     1               DTMP(I) , TMP_OLD(I) , TMP_OLD(I)+DTMP(I)  
          stop 999
        ENDIF
      ENDDO
!
      IMSG_NUM = 0
      XMSG_SQ  = 0
!
      DO I=NSTR,NEND 
        IF ( ABS(DTMP(I)) .GT. 10. ) THEN   ! 000926 30->40
          IMSG_NUM = IMSG_NUM + 1  
          XMSG_SQ  = XMSG_SQ  + (DTMP(I) -10.)**2       
!         write(6,*) 'CANOPY_INTEGRATION_OVERIMP_VEG warning B1 ' ,  
!    1             'I ' , I , 'IMSKLT ' , IMSKLT(I) , 
!    1             ' INSFLAG ' , INSFLAG , 
!    1             'ICG_FLAG '    , ICG_FLAG     , 
!    1             ' CVR_S_OR_N ' , CVR_S_OR_N(I),
!    1             'DTMP TMP_OLD TMP_OLD+DTMP' ,    
!    1               DTMP(I) , TMP_OLD(I) , TMP_OLD(I)+DTMP(I)  
        ENDIF
      ENDDO
!
      if (imsg_num.gt.0) then
        write(6,*) 'canopy_integration_overimp_veg mod more than 10K' , 
     1           'num_point/sum_sq' , IMSG_NUM , XMSG_SQ 
      endif
#endif
!
! ========
! > 準備 <
! ========
!  積分前の状態。
!  水/氷の量、比熱。
!  下草・土壌の場合は、葉の比熱に土壌スキンの比熱を加える。 
!
      DO I=NSTR,NEND 
!
        WTR_LEAF_OLD(I) =  FR_H2O_OLD (I) * FR_WBH_OLD(I)  
     1                                    * H2O_SAT_LEAF(I)
        GLA_LEAF_OLD(I) =  FR_H2O_OLD (I) * ( 1 - FR_WBH_OLD(I) ) 
     1                                    * H2O_SAT_LEAF(I)
!
        CAP_SYSTEM_OLD(I) =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1                      + CAP_WTR_c * WTR_LEAF_OLD(I) 
     1                      + CAP_GLA_c * GLA_LEAF_OLD(I) 
!
        IF ( ICG_FLAG .EQ. 2 ) THEN 
          CAP_SYSTEM_OLD(I) = CAP_SYSTEM_OLD(I) + CAP_GRSK_C         
        ENDIF
!
      ENDDO
!
#ifdef SIB_DEBUG
!
! ----------------------
! > 最初の収支チェック <
! ----------------------
!    もちろん省略可能。
!    系に入ったエネルギーフラックス収束 ( 陰的補正済み ) が
!    系の昇温評価量とつじつまがあうことのチェック。
!    
      DO I=NSTR,NEND 

        X = 
     1     (   
     1         RAD_LONG  (I) 
     1       + RAD_SHRT  (I) 
     1       + SENS (I) 
     1       + EVAP_LEAF (I) * HL_WTR_C 
     1       + SUBL_LEAF (I) * HL_GLA_C                                
     1       + TRNS      (I) * HL_WTR_C                                
     1       + EVAP_GRSK (I) * HL_WTR_C 
     1       + SUBL_GRSK (I) * HL_GLA_C                                
     1       - HDF  (I) 
     1     )   * DELT_CNP 
     1     - DTMP(I) * CAP_SYSTEM_OLD  (I) 
!
        IF ( ABS ( X ) .GT. 1.D-8 ) THEN
          write(6,*) 'canopy_intg_overimp_veg error A1 ',
     1               i,icg_flag,
     1               DTMP(I) , X , 
     1               CAP_SYSTEM_OLD(I),TMP_OLD(I)
          write(6,*) 'each elem' ,
     1         RAD_LONG  (I) ,
     1       + RAD_SHRT  (I) ,
     1       + SENS (I)      ,
     1       + EVAP_LEAF (I) * HL_WTR_C ,
     1       + SUBL_LEAF (I) * HL_GLA_C ,                               
     1       + TRNS      (I) * HL_WTR_C ,                               
     1       + EVAP_GRSK (I) * HL_WTR_C ,
     1       + SUBL_GRSK (I) * HL_GLA_C ,                               
     1       - HDF  (I) , delt_cnp , 
     1       - DTMP(I) * CAP_SYSTEM_OLD  (I) 
          stop 999
        endif
!
!     もとの系のエネルギー 
!       後の収支チェックで使うため
!
        ENG_OLD(I) = 
     1     + CAP_SYSTEM_OLD(I)          * (TMP_OLD(I) - TMP_FREZ_C)
     1     - HL_W2G_C  * GLA_LEAF_OLD(I) 
!
      ENDDO
!
#endif
!
! ====================
! > 水の移動量の確定 <
! ====================
!   
      DO I=NSTR,NEND 
!
!   移動後の量の計算
!       旧値に、蒸発/結露、昇華/結氷を加えた量。
!       TMP_FREZ を越えた変化をしている場合は一般に水氷共存
!
        WTR_LEAF_NEW (I) =   WTR_LEAF_OLD(I) 
     1                     + EVAP_LEAF  (I) * DELT_CNP
!
        GLA_LEAF_NEW (I) =   GLA_LEAF_OLD(I) 
     1                     + SUBL_LEAF  (I) * DELT_CNP
!
        DWTR_LEAF_MOD (I) = 0.               ! 初期化
        DGLA_LEAF_MOD (I) = 0.               ! 初期化
        ENG_STORE_1   (I) = 0.               ! 初期化
        ENG_STORE_2   (I) = 0.               ! 初期化
        ENG_STORE_3   (I) = 0.               ! 初期化
        TMP_MOD       (I) = 0.               ! 初期化
        SENS_MOD      (I) = 0.               ! 初期化
        RAD_LONG_MOD  (I) = 0.               ! 初期化
        HDF_MOD       (I) = 0.               ! 初期化
        ENG_PHASE     (I) = 0.
!
      ENDDO
!
#ifdef SIB_DEBUG
        IMSG_NUM = 0
#endif
!
      DO I=NSTR,NEND 
!
!   ゼロ以下になった場合のエラー処理
!       更新後が零なら出ていくフラックスを減らす。（下向き正なので値は増える）
!       注意。簡単のため氷と水は別々に扱う。
!
!       ENG_STORE_1 : たまる熱量 （ 顕熱で出ていくことにする ）  >= 0
!
        IF ( WTR_LEAF_NEW(I) .LT. 0. ) THEN
          DWTR_LEAF_MOD (I) = WTR_LEAF_NEW(I)             ! < 0
          WTR_LEAF_NEW  (I) = 0. 
          EVAP_LEAF     (I) = EVAP_LEAF(I) - DWTR_LEAF_MOD(I)/DELT_CNP 
          ENG_STORE_1   (I) = ENG_STORE_1(I) 
     1                        - DWTR_LEAF_MOD(I) * HL_WTR_C
#ifdef SIB_DEBUG
          IMSG_NUM = IMSG_NUM+1 
!#else
!          write(6,*) 'CANOPY_INTEGRATION_VEG: INSUFFICIENT WTR' 
!          write(6,*) ' i=', i, ' icg=', icg_flag
#endif
        ENDIF
!
        IF ( GLA_LEAF_NEW(I) .LT. 0. ) THEN
          DGLA_LEAF_MOD (I) = GLA_LEAF_NEW(I)             ! < 0
          GLA_LEAF_NEW  (I) = 0. 
          SUBL_LEAF     (I) = SUBL_LEAF(I) - DGLA_LEAF_MOD(I)/DELT_CNP 
          ENG_STORE_1   (I) = ENG_STORE_1(I) 
     1                        - DGLA_LEAF_MOD(I) * HL_GLA_C
#ifdef SIB_DEBUG
          IMSG_NUM = IMSG_NUM+1 
!#else
!          write(6,*) 'CANOPY_INTEGRATION_VEG: INSUFFICIENT GLA'
!          write(6,*) ' i=', i, ' icg=', icg_flag
#endif
        ENDIF
!
!   確認
!
#ifdef SIB_DEBUG
        IF ( ENG_STORE_1(I).LT.0.D0 ) THEN 
          IMSG_NUM = IMSG_NUM+100000
!#else
!         write(6,*) 'canopy_integration_overimp_veg error A2 ' , 
!    1               I , ENG_STORE_1(I) 
!         stop 999
        ENDIF
#endif
!
      ENDDO
!
#ifdef SIB_DEBUG
      if     ( imsg_num.ge. 100000) then
         write(6,*) 'canopy_overintegration_veg error XXX' 
         stop 999
      elseif ( imsg_num.gt. 0) then
         write(6,*) 'canopy_overintegration_veg insuf wtr or gla ' ,
     1              'num of point/tot =' , imsg_num , nend+1-nstr 
      endif
#endif
!
! -----------------------
! > 収支チェック その 2 <
! -----------------------
!
#ifdef SIB_DEBUG
      DO I=NSTR,NEND 
!
        X  = 
     1     (   
     1         RAD_LONG  (I) 
     1       + RAD_SHRT  (I) 
     1       + SENS (I) 
     1       + EVAP_LEAF (I) * HL_WTR_C 
     1       + SUBL_LEAF (I) * HL_GLA_C                                
     1       + TRNS      (I) * HL_WTR_C                                
     1       + EVAP_GRSK (I) * HL_WTR_C 
     1       + SUBL_GRSK (I) * HL_GLA_C                                
     1       - HDF  (I) 
     1     )  * DELT_CNP 
     1     - CAP_SYSTEM_OLD(I) * DTMP(I) 
     1     - ENG_STORE_1   (I) 
!
        IF ( ABS ( X ) .GT. 1.D-6 ) THEN
          write(6,*) 'canopy_intg_overimp_veg error A3 ',i,icg_flag,
     1       X ,
     1     - CAP_SYSTEM_OLD  (I) * DTMP(I) ,
     1     - ENG_STORE_1   (I) 
          stop 999
        endif
      ENDDO
#endif
!
! -------------------------
! > 収支チェック その 2-2 <
! -------------------------
!   質量増分も加味したチェック
!
#ifdef SIB_DEBUG
!
      DO I=NSTR,NEND 
        ENG_ADD(I) 
     1      = (   
     1            CAP_WTR_C * EVAP_LEAF(I) 
     1                      * (TMP_OLD(I)+DTMP(I)-TMP_FREZ_C) 
     1          + HL_WTR_C  * EVAP_LEAF(I)
     1          + CAP_GLA_C * SUBL_LEAF(I) 
     1                      * (TMP_OLD(I)+DTMP(I)-TMP_FREZ_C)
     1          + HL_WTR_C  * SUBL_LEAF(I)
     1          + RAD_LONG  (I) 
     1          + RAD_SHRT  (I) 
     1          + SENS      (I) 
     1          + TRNS      (I) * HL_WTR_C                                
     1          + EVAP_GRSK (I) * HL_WTR_C 
     1          + SUBL_GRSK (I) * HL_GLA_C                                
     1          - HDF  (I) 
     1        ) * DELT_CNP
     1     - ENG_STORE_1   (I)               ! 顕熱増加分
!
!   注意。
!      ここでは tmp_frez_c を越えても、比熱は変わらないとして扱っている。
!      このループの特例である。
!
        TMP_NEW (I) = TMP_OLD (I) + DTMP (I) 
!
        CAP_SYSTEM_TMP(I) =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1                      + CAP_WTR_c * WTR_LEAF_NEW(I) 
     1                      + CAP_GLA_c * GLA_LEAF_NEW(I) 
        IF ( ICG_FLAG .EQ. 2 ) THEN
          CAP_SYSTEM_TMP(I) = CAP_SYSTEM_TMP(I) + CAP_GRSK_C         
        ENDIF

        ENG_NEW (I) = 
     1     + CAP_SYSTEM_TMP(I)         * (TMP_NEW(I) - TMP_FREZ_C)
     1     - HL_W2G_C  * GLA_LEAF_NEW(I) 
!
        X = ENG_NEW(I) - ENG_OLD(I) - ENG_ADD(I) 
!
        IF ( ABS(X) .GT. 1.D-6 ) THEN
          write(6,*) 'canopy_integration_veg error C1' , I,
     1              X, ENG_OLD(I),
     2              ENG_ADD(I) , ENG_NEW(I) 
          stop 999
        ENDIF 
      ENDDO
#endif
!
! ===============
! >  温度の積分 <
! ===============
!
      DO I=NSTR,NEND
        TMP_NEW (I) = TMP_OLD (I) + DTMP (I) 
      ENDDO
!
! ======================================
! > 零度を跨ぐ変化をするかどうかの判定 <
! ======================================
!
!      符号変化無し OVER_TMP_ZERO=0 
!      負 -> 正     OVER_TMP_ZERO=1 
!      正 -> 負     OVER_TMP_ZERO=2 
!
      DO I=NSTR,NEND
        IOVER_TMP_ZERO (I) = 0 
        IF     (       TMP_NEW(I).GE.TMP_FREZ_C 
     1           .AND. TMP_OLD(I).LE.TMP_FREZ_C ) THEN
          IOVER_TMP_ZERO (I) = 1 
          TMP_MOD  (I) = TMP_FREZ_C - ( TMP_OLD(I) + DTMP(I) ) 
        ELSEIF (       TMP_NEW(I).LE.TMP_FREZ_C 
     1           .AND. TMP_OLD(I).GE.TMP_FREZ_C ) THEN
          IOVER_TMP_ZERO (I) = 2 
          TMP_MOD  (I) = TMP_FREZ_C - ( TMP_OLD(I) + DTMP(I) )  
        ENDIF        
      ENDDO  
!
! ==========================
! > ゼロ度を越えていく場合 <
! ==========================
!
#ifdef SIB_DEBUG
      imsg_num = 0
#endif
!
      DO I=NSTR,NEND
        IF ( IOVER_TMP_ZERO(I) .EQ. 1 ) THEN
!
!   熱のたまり量 （ゼロ度で停止した場合）
!     DRAD_DX,DSENS_DX<0  & HDF_GRSK>0
!
          ENG_STORE_2(I) 
     1           = (        
     1               + DRAD_DX (I) * TMP_MOD(I)    
     1               + DSENS_DX(I) * TMP_MOD(I) 
     1               - HDF_GRSK(I,1) * TMP_MOD(I)  
     1                               * III_040120
     1             ) * DELT_CNP 
     1             -   CAP_SYSTEM_OLD(I) * TMP_MOD(I) 
!
#ifdef SIB_DEBUG
          IF ( ENG_STORE_2(I) .LE. 0. ) THEN 
!           imsg_num = imsg_num + 1 
!#else
             write(6,*) 'canopy_integraton_overimp_veg error A31' , 
     1                  i  , 
     1                  'eng_store2' , ENG_STORE_2(I) ,
     1                  'drad_dx'    , DRAD_DX (I) ,
     1                  'dsens_dx'   , DSENS_DX(I) , 
     1                  'hdf_grsk'   , HDF_GRSK(I,1) , 
     1                  'delt_cnp'   , DELT_CNP ,
     1                  'cap'        , CAP_SYSTEM_OLD(I) , 
     1                  'tmp_mod'    , TMP_MOD(I) 

             stop 999
          ENDIF
#endif
!
!   相変化に要する熱量
!
          ENG_PHASE_ALL(I) =   HL_W2G_C * GLA_LEAF_NEW(I)       
!
#ifdef SIB_DEBUG
!
          IF ( ENG_PHASE_ALL(I) .LT. 0. ) THEN 
!
!           imsg_num = imsg_num + 10 
!#else
            write(6,*) 'canopy_integraton_overimp_veg error A32' , 
     1                  i  , ENG_PHASE_ALL(I) 
            stop 999
          ENDIF
#endif
!
!   比較
!
          IF (    ENG_STORE_2(I) .GT. ENG_PHASE_ALL(I) ) THEN  ! 相変化完了
            ENG_STORE_2(I) = ENG_STORE_2(I) - ENG_PHASE_ALL(I) 
            TMP_NEW(I)     = TMP_FREZ_C + TMP_MIN 
            WTR_LEAF_NEW(I)= WTR_LEAF_NEW(I) + GLA_LEAF_NEW(I)  
            GLA_LEAF_NEW(I)= 0. 
            ENG_PHASE   (I)= ENG_PHASE_ALL(I)
!
          ELSE                                                 ! 相変化で終る
            WTR_LEAF_NEW(I)= WTR_LEAF_NEW(I) 
     1                       + ENG_STORE_2(I) / HL_W2G_C 
            GLA_LEAF_NEW(I)= GLA_LEAF_NEW(I) 
     1                       - ENG_STORE_2(I) / HL_W2G_C 
            TMP_NEW(I)     = TMP_FREZ_C 
            ENG_PHASE  (I) = ENG_STORE_2(I)
            ENG_STORE_2(I) = 0. 
          ENDIF
!
#ifdef SIB_DEBUG
          IF (      WTR_LEAF_NEW(I).LT.0. 
     1         .OR. GLA_LEAF_NEW(I).LT.0. ) THEN
!           imsg_num = imsg_num + 100 
!#else
             write(6,*) 'canopy_integration_overimp_veg error A33'
             stop 999
          ENDIF
#endif
        ENDIF
      ENDDO
!
#ifdef SIB_DEBUG
      if ( imsg_num .gt. 0 ) then
        write(6,*) 'canopy_integration_overimp ERROR A_zero_part',
     1             imsg_num  
        stop 999
      endif
#endif
!
! ==========================
! > ゼロ度を割っていく場合 <
! ==========================
!
#ifdef SIB_DEBUG
      imsg_num = 0 
#endif
!
      DO I=NSTR,NEND
        IF ( IOVER_TMP_ZERO(I) .EQ. 2 ) THEN
!
!   熱のたまり量 （ゼロ度で停止した場合）
!     DRAD_DX,DSENS_DX<0  & HDF_GRSK>0
!
          ENG_STORE_2(I) 
     1              = (                                        ! < 0 
     1                    DRAD_DX     (I) * TMP_MOD(I) 
     1                  + DSENS_DX(I) * TMP_MOD(I) 
     1                  - HDF_GRSK(I,1) * TMP_MOD(I) 
     1                                  * III_040120
     1                ) * DELT_CNP
     1                  - CAP_SYSTEM_OLD(I) * TMP_MOD(I) 
!
#ifdef SIB_DEBUG
          IF ( ENG_STORE_2(I) .GE. 0. ) THEN 
            imsg_num = imsg_num + 1 
!#else
!            write(6,*) 'canopy_integraton_overimp_veg error A31x' , 
!     1                  i  , ENG_STORE_2(I) 
!            stop 999
          ENDIF
#endif
!

!   相変化に要する熱量
!
          ENG_PHASE_ALL(I) = - HL_W2G_C * WTR_LEAF_NEW(I)   ! < 0 
!
#ifdef SIB_DEBUG
          IF ( ENG_PHASE_ALL(I) .GT. 0. ) THEN 
            imsg_num = imsg_num + 10 
!#else
!            write(6,*) 'canopy_integraton_overimp_veg error A32x' , 
!     1                  i  , ENG_PHASE_ALL(I) 
!            stop 999
          ENDIF
#endif
!
!   比較
!
          IF (    ENG_STORE_2(I) .LT. ENG_PHASE_ALL(I) ) THEN  ! 相変化完了
            ENG_STORE_2(I) = ENG_STORE_2(I) - ENG_PHASE_ALL(I) 
            TMP_NEW(I)     = TMP_FREZ_C - TMP_MIN 
            GLA_LEAF_NEW(I)= WTR_LEAF_NEW(I) + GLA_LEAF_NEW(I)  
            WTR_LEAF_NEW(I)= 0. 
            ENG_PHASE  (I) = ENG_PHASE_ALL(I)
!
          ELSE                                                 ! 相変化で終る
            WTR_LEAF_NEW(I)= WTR_LEAF_NEW(I) 
     1                       + ENG_STORE_2(I) / HL_W2G_C 
            GLA_LEAF_NEW(I)= GLA_LEAF_NEW(I) 
     1                       - ENG_STORE_2(I) / HL_W2G_C 
            TMP_NEW(I)     = TMP_FREZ_C 
            ENG_PHASE  (I) = ENG_STORE_2(I)
            ENG_STORE_2(I) = 0. 
          ENDIF
!
#ifdef SIB_DEBUG
          IF (      WTR_LEAF_NEW(I).LT.0. 
     1         .OR. GLA_LEAF_NEW(I).LT.0. ) THEN
            imsg_num = imsg_num + 100 
!#else
!            write(6,*) 'canopy_integration_overimp_veg error A33x'
!            stop 999
!#endif
          endif
#endif
        ENDIF
      ENDDO
!
#ifdef SIB_DEBUG
      if ( imsg_num .gt. 0 ) then
        write(6,*) 'canopy_integration_overimp ERROR B_zero_part',
     1              imsg_num  
        stop 999
      endif
#endif
!
! =========
! > 修正 <
! =========
!
      DO I=NSTR,NEND
        IF ( IOVER_TMP_ZERO(I) .EQ. 0 ) THEN
          CAP =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1          + CAP_WTR_c * WTR_LEAF_NEW(I) 
     1          + CAP_GLA_c * GLA_LEAF_NEW(I) 
          IF ( ICG_FLAG .EQ. 2 ) THEN
            CAP = CAP + CAP_GRSK_C         
          ENDIF
!
          ENG = CAP * ( TMP_NEW(I) - TMP_FREZ_c ) 
     1          - HL_W2G_c * GLA_LEAF_NEW(I) 
          ENX = ENG + HL_W2G_c * ( WTR_LEAF_NEW(I) + GLA_LEAF_NEW(I) ) 
!
          IF ( ENG.GT.0 ) THEN
            WTR_LEAF_NEW(I) = WTR_LEAF_NEW(I) + GLA_LEAF_NEW(I)
            GLA_LEAF_NEW(I) = 0. 
            CAP =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1            + CAP_WTR_c * WTR_LEAF_NEW(I) 
     1            + CAP_GLA_c * GLA_LEAF_NEW(I) 
            IF ( ICG_FLAG .EQ. 2 ) THEN
              CAP = CAP + CAP_GRSK_C         
            ENDIF
            TMP_NEW(I) = TMP_FREZ_C + ENG / CAP
          ELSEIF ( ENX.LT.0 ) THEN
            GLA_LEAF_NEW(I) = WTR_LEAF_NEW(I) + GLA_LEAF_NEW(I)
            WTR_LEAF_NEW(I) = 0. 
            CAP =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1            + CAP_WTR_c * WTR_LEAF_NEW(I) 
     1            + CAP_GLA_c * GLA_LEAF_NEW(I) 
            IF ( ICG_FLAG .EQ. 2 ) THEN
              CAP = CAP + CAP_GRSK_C         
            ENDIF
            TMP_NEW(I) = TMP_FREZ_C + ENX / CAP
          ELSE
            TMP_NEW(I) = TMP_FREZ_C   ! wtr/gla_leaf_new should be modified!! 
!                                     ! Now neglected!!! 
          ENDIF 
        ENDIF
      ENDDO
!
! ================
! > 顕熱等の変更 <
! ================
!
      DO I=NSTR,NEND
        SENS_MOD     (I) = DSENS_DX(I) * TMP_MOD (I) 
        RAD_LONG_MOD (I) = DRAD_DX     (I) * TMP_MOD (I) 
        HDF_MOD      (I) = HDF_GRSK(I,1) * TMP_MOD (I) 
     1                                   * III_040120
!
        SENS_MOD (I) = 
     1         SENS_MOD (I) 
     1         - ( ENG_STORE_1(I) + ENG_STORE_2(I) ) 
     1            / DELT_CNP 
        SENS(I) = SENS(I) + SENS_MOD(I) 

        RAD_LONG (I) = RAD_LONG (I) + RAD_LONG_MOD (I) 
        HDF (I) = HDF (I) + HDF_MOD (I) 
!
        TMP_MOD  (I) = TMP_NEW(I) - ( TMP_OLD(I) + DTMP(I) ) 
!
!      大気へ戻す長波放射の補正
!
        RAD_LONG_NorS_2_ATM (I)  
     1       = RAD_LONG_NorS_2_ATM (I)  
     1         + DRAD_DX(I) * TMP_MOD (I)
!
      ENDDO
!
!
! --------------------------
! > エネルギー収支チェック 3-1 <
! --------------------------
!  ちなみにここで精度が一気に落ちている。( 1.D-8 -> 1.D-3 )
!  何かエラーがあるのかも知れない。
!
#ifdef SIB_DEBUG
      DO I=NSTR,NEND 
!
        X  =  ( 
     1         RAD_LONG  (I) 
     1       + RAD_SHRT  (I) 
     1       + SENS (I) 
     1       + EVAP_LEAF (I) * HL_WTR_C 
     1       + SUBL_LEAF (I) * HL_GLA_C                                
     1       + TRNS      (I) * HL_WTR_C                                
     1       + EVAP_GRSK (I) * HL_WTR_C 
     1       + SUBL_GRSK (I) * HL_GLA_C                                
     1       - HDF  (I) 
     1     )   * DELT_CNP 
     1     - CAP_SYSTEM_OLD  (I) * ( TMP_NEW(I) - TMP_OLD(I) )
     1     - ENG_PHASE (I)

!
        IF ( ABS ( X ) .GT. 1.D-6 ) THEN
          write(6,*) 'CANOPY_INTG_OVERIMP_VEG WARNING A4 ',
     1             i,icg_flag,
     1             X,
     1             CAP_SYSTEM_OLD(I) , TMP_MOD(I) , 
     1             TMP_OLD(I) , TMP_OLD(I)+DTMP(I) , TMP_NEW(I),
     1             IOVER_TMP_ZERO(I) ,
     1             ' E1 ' , ENG_STORE_1(I) , 
     1             ' E2 ' , ENG_STORE_2(I) , 
     1             ' E3 ' , ENG_STORE_3(I),
     1             ' EP ' , ENG_PHASE  (I),
     1             DSENS_DX(I) * TMP_MOD (I) * DELT_CNP ,
     1             DRAD_DX     (I) * TMP_MOD (I) * DELT_CNP ,
     1             HDF_GRSK(I,1) * TMP_MOD (I) * DELT_CNP ,
     1             X / CAP_SYSTEM_OLD(I) ,
     1         RAD_LONG  (I) ,
     1        RAD_SHRT  (I) ,
     1        SENS (I) ,
     1        EVAP_LEAF (I) * HL_WTR_C ,
     1        SUBL_LEAF (I) * HL_GLA_C ,                               
     1        TRNS      (I) * HL_WTR_C ,                               
     1        EVAP_GRSK (I) * HL_WTR_C ,
     1        SUBL_GRSK (I) * HL_GLA_C ,                               
     1       - HDF  (I) ,
     1      DELT_CNP ,
     1     - CAP_SYSTEM_OLD  (I) , TMP_NEW(I) , TMP_OLD(I) 
     1     - ENG_PHASE (I)

          IF ( ABS ( X ) .GT. 1.D-4 ) THEN
            write(6,*) 'canopy_intg_overimp_veg error A4 ' , X
            stop 999
          endif
        ENDIF
!
      ENDDO
#endif
!
! -------------------------
! > 収支チェック その 3-2 <
! -------------------------
!   質量変化分も加味したチェック。
!   はじめのエネルギー + 入ってくるエネルギー + 入ってくる水のエネルギー  
!        = あとのエネルギー 
!
!   ENG_ADD   
!       = 入ってくるエネルギー + 入ってくる水のエネルギー 
!
!   注意。
!     ・出入りする水は、TMP_NEW の温度をもつ。相は TMP_OLD で決まる。
!   検討を要する点（あっているかどうか怪しい点）
!     ・TMP_OLD と TMP_NEW とが TMP_FREZ_C を越えた変化をする場合。
!       TMP_FREZ_C をはさんで比熱は CAP_WTR_C または CAP_GLA_C をとるはず。
!       しかし、これを実現するのは難しい。
!       今(000917)は、基本的に TMP_FREZ_C を越えた変化はほとんどないとして
!       いるため実害はない。
!       つまり、越えた変化をする場合、ごくわずかにしか越えないようにしてい
!       る。事実上、過去の CAP_*_C を使って問題ない。
!
      DO I=NSTR,NEND 
        ENG_ADD(I) 
     1      = (   
!   入ってくる/出ていく水蒸気のエネルギー その1
!   基準は 0 度の液水
!   比熱は、0 度より高い場合は CAP_WTR_C 、低い場合は CAP_GLA_C 
!   水蒸気と液水の相変化に必要な潜熱は HL_WTR_C 
     1            CAP_WTR_C * EVAP_LEAF(I) * (TMP_NEW(I)-TMP_FREZ_C) 
     1          + HL_WTR_C  * EVAP_LEAF(I)
     1          + CAP_GLA_C * SUBL_LEAF(I) * (TMP_NEW(I)-TMP_FREZ_C)
     1          + HL_WTR_C  * SUBL_LEAF(I)
!   その他の熱フラックス
     1          + RAD_LONG  (I) 
     1          + RAD_SHRT  (I) 
     1          + SENS (I) 
     1          + TRNS      (I) * HL_WTR_C                                
     1          + EVAP_GRSK (I) * HL_WTR_C 
     1          + SUBL_GRSK (I) * HL_GLA_C                                
     1          - HDF  (I) 
     1        ) * DELT_CNP
!
!    新たな比熱 
!
        CAP_SYSTEM_AFTER(I) =   CAP_LAI_c * XLAI_VEG_msx(I) 
     1                        + CAP_WTR_c * WTR_LEAF_NEW(I) 
     1                        + CAP_GLA_c * GLA_LEAF_NEW(I) 
!
        IF ( ICG_FLAG .EQ. 2 ) THEN
          CAP_SYSTEM_AFTER(I) = CAP_SYSTEM_AFTER(I) + CAP_GRSK_C  
        ENDIF
!
        ENG_NEW (I) = 
     1     + CAP_SYSTEM_AFTER(I)  * (TMP_NEW(I) - TMP_FREZ_C)
     1     - HL_W2G_C  * GLA_LEAF_NEW(I) 
!
#ifdef SIB_DEBUG
        X = ENG_NEW(I)-ENG_OLD(I)-ENG_ADD(I)
!
        IF ( ABS(X) .GT. 1.D-4 ) THEN
          write(6,*) 'canopy_integration_veg error C2' , 
     1            I , ICG_FLAG , 
!    1            IOVER_TMP_ZERO(I) ,
     1            X , 
     1            ENG_OLD(I) , ENG_ADD(I) , ENG_NEW(I) 
          stop 999
!
        ELSEIF ( ABS(X) .GT. 1.D-6 ) THEN
          write(6,*) 'canopy_integration_veg warning diff large C2' , 
     1            I , ICG_FLAG , 
!    1            IOVER_TMP_ZERO(I) ,
     1            X , 
     1            ENG_OLD(I) , ENG_ADD(I) , ENG_NEW(I) 
        ENDIF 
#endif
      ENDDO
!
! ========================
! > 水が溜ったときの処理 <
! ========================
!   飽和分以上にたまったら、時間刻の最後に落す。
!
#ifdef SIB_DEBUG
      icount = 0
#endif
      DO I=NSTR,NEND
!
        FRAC = ( WTR_LEAF_NEW (I) + GLA_LEAF_NEW(I) ) / H2O_SAT_LEAF(I) 
        FRAC_DROP = 0.
        IF ( FRAC .GT. 1 ) THEN
#ifdef SIB_DEBUG
          icount = icount + 1 
#endif
          FRAC_DROP = (FRAC-1)/FRAC
        ENDIF
        FRAC_STAY = 1 - FRAC_DROP 
!
        WTR_LEAF_DROP (I) = WTR_LEAF_NEW (I) * FRAC_DROP
        GLA_LEAF_DROP (I) = GLA_LEAF_NEW (I) * FRAC_DROP
        WTR_LEAF_NEW  (I) = WTR_LEAF_NEW (I) * FRAC_STAY
        GLA_LEAF_NEW  (I) = GLA_LEAF_NEW (I) * FRAC_STAY
!
        FR_WTR_NEW (I) = WTR_LEAF_NEW  (I) / H2O_SAT_LEAF(I) 
        FR_GLA_NEW (I) = GLA_LEAF_NEW  (I) / H2O_SAT_LEAF(I) 
      ENDDO
#ifdef SIB_DEBUG
      if ( icount .gt. 0 ) then
        write(6,*) 'canopy_integ_ovimp_veg h2o is too much so drop' ,  
     1                                     ' ins= ' , insflag , 
     1                                     ' icg= ' , icg_flag , 
     1                         ' point/tot= ' , icount , nend-nstr+1
      endif
#endif
!
! ----------------
! > 値のチェック <
! ----------------
!
#ifdef SIB_DEBUG
!
#include "noparallel_fixedform"
      DO I=NSTR,NEND 
        IF ( TMP_NEW(I) .LT. 175. .OR. TMP_NEW(I).GT.340. ) THEN
          WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP_VEG WARNING/ERROR B2', 
     1              'I ' , I , 'IMSKLT' , IMSKLT(I) ,  I , 
     1              'ICG_FLAG' ,     ICG_FLAG ,
     1              ' TMP ' , TMP_NEW(I) , TMP_OLD(I) , DTMP (I) ,
     1              ' CAP ' , CAP_SYSTEM_OLD(I) , 
     1              ' INPUT ' , 
     1        (
     1            RAD_SHRT  (I) 
     1          + RAD_LONG  (I) 
     1          + SENS (I) 
     1          + EVAP_LEAF (I) * HL_WTR_C 
     1          + SUBL_LEAF (I) * HL_GLA_C                                
     1          + TRNS      (I) * HL_WTR_C                                
     1          + EVAP_GRSK (I) * HL_WTR_C 
     1          + SUBL_GRSK (I) * HL_GLA_C                                
     1          - HDF  (I) 
     1        )  * DELT_CNP ,
     1            RAD_SHRT  (I) ,
     1            RAD_LONG  (I) ,
     1          + SENS (I) ,
     1          + EVAP_LEAF (I) * HL_WTR_C ,
     1          + SUBL_LEAF (I) * HL_GLA_C ,                               
     1          + TRNS      (I) * HL_WTR_C ,                               
     1          + EVAP_GRSK (I) * HL_WTR_C ,
     1          + SUBL_GRSK (I) * HL_GLA_C ,                               
     1          - HDF  (I) ,
     1          ' PHASE ' , ENG_PHASE  (I)

          IF ( TMP_NEW(I) .LT. 160. .OR. TMP_NEW(I).GT.350. ) THEN
            STOP 999
          ENDIF
        ENDIF
      ENDDO
#endif
! 
! ==================
! > バケツへの格納 <
! ==================
!
      IF ( ICG_FLAG .EQ. 2 ) THEN                 ! 地面
!
        DO I=NSTR,NEND
!
!  熱バケツ
!
          ENG_BUCKET     (I) = ENG_BUCKET     (I) 
     1                        + HDF (I) * DELT_CNP * CVR_S_OR_N(I)
!
!  水/氷バケツ in  （結露・結霜 = 土壌第一層の温度）
!
          H2O_BUCKET_INTO (I) 
     1       = H2O_BUCKET_INTO(I) 
     1         +  MAX(EVAP_GRSK(I),ZERO) * DELT_CNP * CVR_S_OR_N(I)
!
          ENH_BUCKET_INTO(I) 
     1       = ENH_BUCKET_INTO(I) 
     1         + MAX(EVAP_GRSK(I),ZERO) * DELT_CNP 
     1              * CAP_WTR_C * ( TMP_SOIL(I) - TMP_FREZ_C ) 
     1              * CVR_S_OR_N(I)
!
          H2O_BUCKET_INTO(I) 
     1       = H2O_BUCKET_INTO(I)  
     1         +  MAX(SUBL_GRSK(I),ZERO) * DELT_CNP * CVR_S_OR_N(I)
!
          ENH_BUCKET_INTO(I) 
     1     = ENH_BUCKET_INTO(I)  
     1       + MAX(SUBL_GRSK(I),ZERO) * DELT_CNP 
     1          * (  - HL_W2G_C 
     1               + CAP_GLA_C * ( TMP_SOIL (I) - TMP_FREZ_C ) )
     1          * CVR_S_OR_N(I)
!
!  水/氷バケツ out （蒸発・昇華 = 土壌第一層の温度）
!
          H2O_BUCKET_FROM(I) 
     1       = H2O_BUCKET_FROM(I) 
     1         + MIN(EVAP_GRSK(I),ZERO) * DELT_CNP * CVR_S_OR_N(I) 
!
          ENH_BUCKET_FROM(I) 
     1       = ENH_BUCKET_FROM(I) 
     1         + MIN(EVAP_GRSK(I),ZERO) * DELT_CNP 
     1           * CAP_WTR_C * ( TMP_SOIL(I) - TMP_FREZ_C ) 
     1           * CVR_S_OR_N (I)
!
          H2O_BUCKET_FROM(I) 
     1       = H2O_BUCKET_FROM(I)  
     1         + MIN(SUBL_GRSK(I),ZERO) * DELT_CNP * CVR_S_OR_N(I) 
!
          ENH_BUCKET_FROM(I) 
     1       = ENH_BUCKET_FROM(I)  
     1         + MIN(SUBL_GRSK(I),ZERO) * DELT_CNP 
     1            * (  - HL_W2G_C 
     1                 + CAP_GLA_C * ( TMP_SOIL(I) - TMP_FREZ_C ) )
     1            * CVR_S_OR_N(I)
        ENDDO
      ENDIF
!
!  落下する水・氷をバケツに。
!    温度は葉の最終温度
!
      DO I=NSTR,NEND
!
          H2O_BUCKET_INTO(I) 
     1       = H2O_BUCKET_INTO(I) 
     1         + WTR_LEAF_DROP(I) * CVR_S_OR_N(I)
!
          ENH_BUCKET_INTO(I) 
     1       = ENH_BUCKET_INTO(I) 
     1         + WTR_LEAF_DROP(I) 
     1              * CAP_WTR_C * ( TMP_NEW(I) - TMP_FREZ_C ) 
     1              * CVR_S_OR_N(I)
!
          H2O_BUCKET_INTO(I) 
     1       = H2O_BUCKET_INTO(I)  
     1         + GLA_LEAF_DROP(I) * CVR_S_OR_N(I)
!
          ENH_BUCKET_INTO(I) 
     1     = ENH_BUCKET_INTO(I)  
     1         + GLA_LEAF_DROP(I)  
     1              * (  - HL_W2G_C 
     1                   + CAP_GLA_C * ( TMP_NEW (I) - TMP_FREZ_C ) )
     1              * CVR_S_OR_N(I)
      ENDDO      
!
! ==============
! > 詳細モニタ <
! ==============
!   放射補正量
!     
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!       TC/TG 変更で CNP/GRD から大気にいく放射の変更
!       大気が受け取る分 = - 自分が受け取る分
!       下向き正
!       
        IF  ( ICG_FLAG .EQ. 1 ) THEN         ! キャノピー
          IF     ( INSFLAG.EQ.0 ) THEN          ! 雪なし
            CMARK ='HLRNAE' 
            CALL MONIT_SIB ( CMARK, RAD_LONG_MOD, ICN, DELT_CNP, ONE)       
          ELSEIF ( INSFLAG.EQ.1 ) THEN          ! 雪あり
            CMARK ='HLRSAE' 
            CALL MONIT_SIB ( CMARK, RAD_LONG_MOD, ICS, DELT_CNP, ONE)       
          ENDIF
!
        ELSEIF ( ICG_FLAG .EQ. 2 ) THEN      ! 地面
          CMARK ='HLRGAE' 
          CALL MONIT_SIB ( CMARK, RAD_LONG_MOD, ICN, DELT_CNP, ONE)       
!
!       スキンの蓄える熱 (+下草からでていったH2Oの昇温に要したエネルギー）
!
           CMARK ='HGETSK' 
           DO I=NSTR,NEND
             WORK_IDX(I) =   
     1                       CAP_SYSTEM_AFTER(I) * TMP_NEW(I)
     1                     - CAP_SYSTEM_OLD  (I) * TMP_OLD(I)
           ENDDO
           CALL MONIT_SIB ( CMARK, WORK_IDX , ICN , ONE, ONE)       

!
        ENDIF
!
      ENDIF
#endif      
!
      RETURN
      END SUBROUTINE  CANOPY_INTEGRATION_OVERIMP_VEG 
!===================================================================
      SUBROUTINE CANOPY_INTEGRATION_OVERIMP_SNOW (
     I   NSTR       , NEND      , 
     I   DELT_CNP   , CVR_SNW   , 
     I   TMP_SNSK_OLD           , XDTMP_SKIN , 
     I   TMP_SNOW   ,              
     I   EVAP       , SUBL       , 
     M   SENS       , 
     M   RAD_LONG   , 
     I   RAD_SHRT   , 
     M   HDF      , 
     I   DSENS_DX   , 
     I   DRAD_DX    , 
     I   HDF_SNSK ,
     I   RAD_SHRT_INTO_SNOW_LAYER     ,
     O   TMP_SNSK_NEW            ,    
     M   RAD_LONG_NorS_2_ATM,   
     M   ENG_BUCKET_SNW          , 
     M   H2O_BUCKET_INTO_SNOW    , H2O_BUCKET_FROM_SNOW    , 
     M   ENH_BUCKET_INTO_SNOW    , ENH_BUCKET_FROM_SNOW    , 
     M   SRAD_BUCKET_SNW   )
!
!  雪スキンの計算。
!
      USE SIBPRM , ONLY :
     1   IDX , ICG , ICGR , ICGB , IEST
!
      USE SIBCON , ONLY :
     1   HL_GLA_c    , 
     1   HL_W2G_c    , 
     1   HL_WTR_c    , 
     1   TMP_FREZ_c  , 
     1   CAP_GLA_C   , 
     1   CAP_SNSK_c 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICS,
     1   IMONIT_LEVEL ,
     1   IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR              ! 開始格子
      INTEGER,INTENT(IN)   :: NEND              ! 終了格子
      REAL(8),INTENT(IN)   :: DELT_CNP
      REAL(8),INTENT(INOUT):: SENS        (IDX)  ! 
      REAL(8),INTENT(INOUT):: RAD_LONG    (IDX)  ! 
      REAL(8),INTENT(IN)   :: RAD_SHRT    (IDX)  ! 
      REAL(8),INTENT(INOUT):: HDF         (IDX)  ! 
      REAL(8),INTENT(IN)   :: EVAP        (IDX)  ! 
      REAL(8),INTENT(IN)   :: SUBL        (IDX)  ! 
      REAL(8),INTENT(IN)   :: TMP_SNSK_OLD(IDX)    
      REAL(8),INTENT(IN)   :: TMP_SNOW    (IDX)    ! 積雪第一層温度
      REAL(8),INTENT(IN)   :: XDTMP_SKIN  (IDX)    
!
      REAL(8),INTENT(IN)   :: RAD_SHRT_INTO_SNOW_LAYER (IDX) 
!
      REAL(8),INTENT(IN)   :: DSENS_DX   (IDX) 
      REAL(8),INTENT(IN)   :: DRAD_DX    (IDX) 
      REAL(8),INTENT(IN)   :: HDF_SNSK   (IDX,2) 
!
      REAL(8),INTENT(IN)   :: CVR_SNW   (IDX) 
!     REAL(8),INTENT(IN)   :: FLAG_TMPMOD_SKIN  (IDX)
!
      REAL(8),INTENT(INOUT):: RAD_LONG_NorS_2_ATM(IDX)   
!
! INOUT
!
      REAL(8),INTENT(INOUT):: ENG_BUCKET_SNW        (IDX)
      REAL(8),INTENT(INOUT):: SRAD_BUCKET_SNW       (IDX)
!
      REAL(8),INTENT(INOUT):: H2O_BUCKET_INTO_SNOW     (IDX)
      REAL(8),INTENT(INOUT):: H2O_BUCKET_FROM_SNOW    (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_INTO_SNOW     (IDX)
      REAL(8),INTENT(INOUT):: ENH_BUCKET_FROM_SNOW    (IDX)
!
! OUT 
!
      REAL(8),INTENT(OUT)  :: TMP_SNSK_NEW  (IDX)    
!
! WORK 
! 
      INTEGER   :: I
      REAL(8)   :: X
!
      REAL(8)   :: SENS_MOD  (IDX)
      REAL(8)   :: RAD_LONG_MOD   (IDX)
      REAL(8)   :: HDF_MOD   (IDX)
      REAL(8)   :: ENG_STORE (IDX)
      REAL(8)   :: DTMP_SKIN_MOD (IDX)    ! 温度変更量再変更分
!
#ifdef MONYOS
!     REAL(8),PARAMETER :: ONE_MINUS = -1.D0 
      REAL(8),PARAMETER :: ONE       = 1.D0 
      CHARACTER(7) :: CMARK 
      REAL(8)      :: WORK_IDX(IDX)
#endif
!
      REAL(8)   :: ZERO  = 0.D0
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTEGRATION_OVERIMP_SNOW 00/06 POCHAKA '
        LFIRST = .FALSE. 
      ENDIF
!
      TMP_SNSK_NEW(:) = 0.

!
!  ------------------
!  > 収支のチェック <
!  ------------------
!
#ifdef SIB_DEBUG
      DO I=NSTR,NEND
!
        X = 
     1     (   
     1         RAD_LONG(I)
     1       + RAD_SHRT(I)
     1       + SENS (I) 
     1       + EVAP (I) * HL_WTR_C 
     1       + SUBL (I) * HL_GLA_C                                
     1       - HDF  (I)          
     1     ) * DELT_CNP 
     1     - CAP_SNSK_c * XDTMP_SKIN(I) 
!
!
        IF ( ABS(X).GT.1.D-6 ) THEN
          write(6,*) 'canopy_integration_overimp_snow warning A ',i,X,
     1               X / CAP_SNSK_c 
          write(6,*) 
     1       + RAD_LONG(I)        , 
     1       + RAD_SHRT(I)        , 
     1       + SENS (I)             ,
     1       + EVAP (I) * HL_WTR_C  ,
     1       + SUBL (I) * HL_GLA_C  ,                               
     1       - HDF  (I)             ,
     1     - CAP_SNSK_c * XDTMP_SKIN(I) , X / delt_cnp 

!
          IF ( ABS(X).GT.1.D-4) THEN
            write(6,*) 'canopy_integration_overimp_snow error A ' ,i,X,
     1               X / CAP_SNSK_c 
            stop 999
          ENDIF
        ENDIF
!
      ENDDO
#endif
!
! ==============
! > 温度の積分 <
! ==============
!
      DO I=NSTR,NEND
        DTMP_SKIN_MOD(I) = 0. 
        SENS_MOD     (I) = 0.
        RAD_LONG_MOD (I) = 0.
        HDF_MOD      (I) = 0.
        ENG_STORE    (I) = 0.
!
        TMP_SNSK_NEW (I) = TMP_SNSK_OLD(I) + XDTMP_SKIN(I)
!
        IF ( TMP_SNSK_NEW(I) .GE. TMP_FREZ_C ) THEN
          DTMP_SKIN_MOD(I) = TMP_FREZ_C - TMP_SNSK_NEW (I) 
          TMP_SNSK_NEW (I) = TMP_FREZ_C 
!
          SENS_MOD(I)      = DSENS_DX(I) * DTMP_SKIN_MOD(I) 
          RAD_LONG_MOD (I) = DRAD_DX (I) * DTMP_SKIN_MOD(I) 
          HDF_MOD (I)      = HDF_SNSK(I,1) * DTMP_SKIN_MOD(I) 
!!!!                                       * III_040120
!
          SENS    (I)      = SENS(I) + SENS_MOD(I) 
          RAD_LONG(I)      = RAD_LONG (I) + RAD_LONG_MOD (I) 
          HDF     (I)      = HDF (I) + HDF_MOD (I) 
!
          ENG_STORE (I)    = (   SENS_MOD(I)
     1                         + RAD_LONG_MOD (I)
     1                         - HDF_MOD (I) ) * DELT_CNP 
     1                       - CAP_SNSK_c * DTMP_SKIN_MOD(I) 
!
          HDF     (I)      = HDF (I) + ENG_STORE(I) / DELT_CNP
!
!      大気へ戻す長波放射の補正
!
          RAD_LONG_NorS_2_ATM (I)  
     1       = RAD_LONG_NorS_2_ATM (I)  
     1         + DRAD_DX(I) * DTMP_SKIN_MOD(I) 
!
        ENDIF
      ENDDO
!
!  ------------------
!  > 収支のチェック <
!  ------------------
!
#ifdef SIB_DEBUG
      DO I=NSTR,NEND
!
        X = 
     1     (   RAD_LONG(I) 
     1       + RAD_SHRT(I)
     1       + SENS (I) 
     1       + EVAP (I) * HL_WTR_C 
     1       + SUBL (I) * HL_GLA_C                                
     1       - HDF  (I)          
     1     ) * DELT_CNP 
     1     - CAP_SNSK_c * ( TMP_SNSK_NEW(I) - TMP_SNSK_OLD(I) )
!    1     - ENG_STORE(I) 
!
        IF ( ABS(X).GT.1.D-6 ) THEN
          write(6,*) 'canopy_integration_overimp_snow warning A2' ,i,X,
     1               X / CAP_SNSK_c ,TMP_SNSK_NEW(I) 
          write(6,*) 
     1         RAD_LONG  (I)             ,
     1         RAD_SHRT  (I) , 
     1       + SENS (I)             ,
     1       + EVAP (I) * HL_WTR_C  ,
     1       + SUBL (I) * HL_GLA_C  ,                               
     1       - HDF  (I)             ,
     1       - CAP_SNSK_c * ( TMP_SNSK_NEW(I) - TMP_SNSK_OLD(I) ),
     1       X / delt_cnp 

          IF ( ABS(X).GT.1.D-4) THEN
            write(6,*) 'canopy_integration_overimp_snow error A2' ,i,X,
     1               X / CAP_SNSK_c ,TMP_SNSK_NEW(I) 
            stop 999
          ENDIF
        ENDIF
!
      ENDDO
#endif
!
! ==================
! > バケツへの格納 <
! ==================
!
      DO I=NSTR,NEND
!
!  熱バケツ
!
        ENG_BUCKET_SNW(I) = ENG_BUCKET_SNW(I) 
     1                      + HDF (I)  * DELT_CNP * CVR_SNW(I)
!
!  水/氷バケツ in
!
        H2O_BUCKET_INTO_SNOW(I) 
     1     = H2O_BUCKET_INTO_SNOW(I) 
     1       + MAX(EVAP(I),ZERO) * DELT_CNP * CVR_SNW(I) 
!
        ENH_BUCKET_INTO_SNOW(I) 
     1     = ENH_BUCKET_INTO_SNOW(I) 
!           零度の時のみなので、エネルギーは零。
!
        H2O_BUCKET_INTO_SNOW(I) 
     1     = H2O_BUCKET_INTO_SNOW(I)  
     1       + MAX(SUBL(I),ZERO) * DELT_CNP * CVR_SNW(I) 
!
        ENH_BUCKET_INTO_SNOW(I) 
     1     = ENH_BUCKET_INTO_SNOW(I)  
     1       + MAX(SUBL(I),ZERO) 
     1          * (  - HL_W2G_C 
     1               + CAP_GLA_C * ( TMP_SNOW (I) - TMP_FREZ_C ) )
     1          * DELT_CNP * CVR_SNW(I)
!
!  水/氷バケツ out
!
        H2O_BUCKET_FROM_SNOW(I) 
     1     = H2O_BUCKET_FROM_SNOW(I) 
     1       + MIN(EVAP(I),ZERO) * DELT_CNP * CVR_SNW(I) 
!
        ENH_BUCKET_FROM_SNOW(I) 
     1     = ENH_BUCKET_FROM_SNOW(I) 
!           零度の時のみ なので、エネルギーは零。
!
        H2O_BUCKET_FROM_SNOW(I) 
     1     = H2O_BUCKET_FROM_SNOW(I)  
     1       + MIN(SUBL(I),ZERO) * DELT_CNP * CVR_SNW(I) 
!
        ENH_BUCKET_FROM_SNOW(I) 
     1     = ENH_BUCKET_FROM_SNOW(I)  
     1       + MIN(SUBL(I),ZERO) 
     1          * (  - HL_W2G_C 
     1               + CAP_GLA_C * ( TMP_SNOW (I) - TMP_FREZ_C ) )
     1          * DELT_CNP * CVR_SNW(I)
!
!  雪への入射短波
!
          SRAD_BUCKET_SNW(I) 
     1      = SRAD_BUCKET_SNW(I) 
     1        +  RAD_SHRT_INTO_SNOW_LAYER(I) * CVR_SNW(I) 
     1           * DELT_CNP 
!
      ENDDO
!
! ==============
! > 詳細モニタ <
! ==============
!   放射補正量
!     
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!       TG 変更で SNOW から大気にいく放射の変更
!       大気が受け取る分 = - 自分が受け取る分
!       下向き正
!
        CMARK ='HLRWAE' 
        CALL MONIT_SIB ( CMARK, RAD_LONG_MOD , ICS , DELT_CNP, ONE)       
!
!       スキンの蓄える熱
!
        CMARK ='HGETSK' 
        DO I=NSTR,NEND
          WORK_IDX(I) = CAP_SNSK_c *(TMP_SNSK_NEW(I)-TMP_SNSK_OLD(I))
        ENDDO
        CALL MONIT_SIB ( CMARK, WORK_IDX , ICS , ONE, ONE)       
!
      ENDIF
#endif      
!
      RETURN
      END SUBROUTINE   CANOPY_INTEGRATION_OVERIMP_SNOW 
!==========================================================
      SUBROUTINE CANOPY_INTEGRATION_SS1 (
     I    NSTR     , NEND    , INSFLAG   , 
     I    DELT_CNP , CAP_SS1 , HDF_TMP    , FLAG_PHASE_SS1 , 
     M    TMP_SS1 ) 
!
      USE SIBPRM , ONLY :
     1    IDX  
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1    RCN_FACT_SS1 
!
      USE SIBCON, ONLY :
     1    TMP_FREZ_C 
!
#ifdef MONYOS
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS,INN,INS,
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
#endif
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: INSFLAG
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
!
      REAL(8),INTENT(IN)    :: DELT_CNP 
      REAL(8),INTENT(IN)    :: CAP_SS1 (IDX)
      REAL(8),INTENT(IN)    :: HDF_TMP(IDX,2)
      REAL(8),INTENT(IN)    :: FLAG_PHASE_SS1(IDX)
!
      REAL(8),INTENT(INOUT) :: TMP_SS1 (IDX)
      REAL(8)               :: TMP_SAVE(IDX)
!
      INTEGER   :: I
!
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      CHARACTER(7) :: CMARK 
#endif
!
      DO I=NSTR,NEND
        TMP_SAVE(I) = TMP_SS1(I) 
        TMP_SS1(I)  = TMP_SS1(I) 
     1       + ( HDF_TMP(I,1) - HDF_TMP(I,2) ) * DELT_CNP / CAP_SS1(I) 
     1         * FLAG_PHASE_SS1(I) * RCN_FACT_SS1 
!!        IF ( (TMP_SAVE(I)-TMP_FREZ_C)*(TMP_SS1(I)-TMP_FREZ_C)
!!     1       .LE.0.D0 ) THEN
!!          TMP_SS1(I) = TMP_FREZ_C 
!!        ENDIF 
      ENDDO
!
!  --------------
!  > 詳細モニタ <
!  --------------
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'TMPL1XL' 
          CALL MONIT_SIB ( CMARK, TMP_SS1, ICN, DELT_CNP, ONE)       
        ELSE
          CMARK = 'TMPW1XL' 
          CALL MONIT_SIB ( CMARK, TMP_SS1, ICS, DELT_CNP, ONE)       
        ENDIF  
!
        IF (  ICN_SIB0109_MONITOR .EQ. 1 ) THEN
          IF ( INSFLAG.EQ.0 ) THEN
            CMARK = 'TMPL1XV' 
            CALL MONIT_SIB ( CMARK, TMP_SS1, INN, DELT_CNP, ONE)       
          ELSE
            CMARK = 'TMPW1XV' 
            CALL MONIT_SIB ( CMARK, TMP_SS1, INS, DELT_CNP, ONE)       
          ENDIF  
        ENDIF
      ENDIF
#endif
!
      RETURN
      END SUBROUTINE CANOPY_INTEGRATION_SS1 
!========================================================
      SUBROUTINE MR_CANOPY_INTEGRATION_SS1 
!
      USE SIB_MONIT , ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
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
        WRITE(6,*) 'MONIT_REGIST_TMP_SS1 00/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
        CMARK  = 'TMPL1XV'
        CTITLE = 'TEMPERATURE OF SOIL1 (FOR NOS CANOPY)'
        CUNIT  = 'K'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'TMPW1XV'
        CTITLE = 'TEMPERATURE OF SNOW1 (FOR SNW CANOPY)'
        CUNIT  = 'K'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'TMPL1XL' 
        CTITLE = 'TEMPERATURE OF SOIL1 (FOR NOS CANOPY)'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'TMPW1XL' 
        CTITLE = 'TEMPERATURE OF SNOW1 (FOR SNW CANOPY)'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      RETURN
      END SUBROUTINE MR_CANOPY_INTEGRATION_SS1 
!================================================================
      SUBROUTINE CANOPY_INTEGRATION_GET_GMTQ (
     I     NSTR , NEND , 
     I     SENS_TMP , EVAP_TMP , SUBL_TMP , TRNS_TMP ,  
     I     GMT_IN   , GMQ_IN   ,
     O     GMT_OUT  , GMQ_OUT   )
!
      USE SIBPRM, ONLY : IDX , ICG , ICGB
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
!
!
      REAL(8),INTENT(IN)   :: SENS_TMP(IDX,ICG)
      REAL(8),INTENT(IN)   :: EVAP_TMP(IDX,ICGB)
      REAL(8),INTENT(IN)   :: SUBL_TMP(IDX,ICGB)
      REAL(8),INTENT(IN)   :: TRNS_TMP(IDX,ICG)
!
      REAL(8),INTENT(IN)  :: GMT_IN(IDX,2)
      REAL(8),INTENT(IN)  :: GMQ_IN(IDX,2)
      REAL(8),INTENT(OUT)  :: GMT_OUT(IDX)
      REAL(8),INTENT(OUT)  :: GMQ_OUT(IDX)
!
      INTEGER :: I
!
      DO I=NSTR,NEND
        GMT_OUT(I  ) = ( GMT_IN(I,2) 
     1                   - ( SENS_TMP(I,1) + SENS_TMP(I,2) ) 
     1                 )  
     1                  / GMT_IN(I,1)  
        GMQ_OUT(I)   = ( GMQ_IN(I,2) 
     1                   - (   EVAP_TMP(I,1) + EVAP_TMP(I,2) 
     1                       + EVAP_TMP(I,3)  
     1                       + SUBL_TMP(I,1) + SUBL_TMP(I,2) 
     1                       + SUBL_TMP(I,3) 
     1                       + TRNS_TMP(I,1) + TRNS_TMP(I,2) )
     1                   )  / GMQ_IN(I,1)  
!
!!      GMT_OUT(I,2) = - ( SENS_TMP(I,1) + SENS_TMP(I,2) ) 
!!   1                   / GMT_OUT(I,1)  
!!      GMQ_OUT(I,2) = - (   EVAP_TMP(I,1) + EVAP_TMP(I,2) 
!!   1                     + SUBL_TMP(I,1) + SUBL_TMP(I,2) 
!!   1                     + TRNS_TMP(I,1) + TRNS_TMP(I,2) 
!!   1                   )  / GMQ_OUT(I,1)  
      ENDDO
!
      RETURN
      END SUBROUTINE  CANOPY_INTEGRATION_GET_GMTQ 
!================================================================
      SUBROUTINE MR_CANOPY_INTEGRATION_OVERIMP
!
      USE SIB_MONIT, ONLY :
     1    IMONIT_LEVEL ,
     1    IMONIT_ALL
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
        WRITE(6,*) 'MR_CANOPY_INTEGRATION_OVETIMP ' // 
     1             ' 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL.GE.IMONIT_ALL) THEN
!
      CMARK  = 'WMVNT2E' 
      CTITLE = 'WATER MODIFY FROM TRNS TO EVAP AT CANOPY NOS'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WMVGT2E' 
      CTITLE = 'WATER MODIFY FROM TRNS TO EVAP AT CANOPY GRASS'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WMVST2E' 
      CTITLE = 'WATER MODIFY FROM TRNS TO EVAP AT CANOPY SNW'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!
      CMARK  = 'HLRNAD' 
      CTITLE = 'LONG RAD FROM CNP TO ATM NOS implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRGAD' 
      CTITLE = 'LONG RAD FROM GROUND TO ATM implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRNGD' 
      CTITLE = 'LONG RAD FROM CNP NOS TO GROUND implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRGND' 
      CTITLE = 'LONG RAD FROM GROUND TO CNP NOS implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRSAD' 
      CTITLE = 'LONG RAD FROM CNP TO ATM SNW implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRWAD' 
      CTITLE = 'LONG RAD FROM SNOW TO ATM implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRSWD' 
      CTITLE = 'LONG RAD FROM CNP SNW TO SNOW implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRWSD' 
      CTITLE = 'LONG RAD FROM SNOW TO CNP SNW implicit correct' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 最終補正
!
      CMARK  = 'HLRNAE' 
      CTITLE = 'LONG RAD FROM CNP TO ATM NOS last mod' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRGAE' 
      CTITLE = 'LONG RAD FROM GROUND TO ATM last mod' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRSAE' 
      CTITLE = 'LONG RAD FROM CNP TO ATM SNW last mod' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRWAE' 
      CTITLE = 'LONG RAD FROM SNOW TO ATM last mod' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK   ='HGETSK' 
      CTITLE = 'ENERGY STORE AT SKIN' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )

!
      ENDIF
!
      RETURN
      END SUBROUTINE MR_CANOPY_INTEGRATION_OVERIMP

      end  module CANOPY_INTEGRATION_OVERIMP 



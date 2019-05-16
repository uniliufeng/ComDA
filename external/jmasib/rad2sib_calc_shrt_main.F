! 放射-> SiB 関連メインプログラム
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/30 保坂征宏 コメント等変更。
!
!  1 時間に一度だけ呼ばれる
!
      SUBROUTINE RAD2SIB_CALC_SHRT_MAIN ( 
     I   NMAX        , NMAX_01        , 
     I   NTYP        , MON            , 
     I   RAD_SHRT_INTO_SIB , 
     I   CVR_VEG_ms  , 
     I   EXTK_NOS    , ABSORB_SHRT_NOS, REACH_GRD_NOS , 
     I   EXTK_SNW    , ABSORB_SHRT_SNW, REACH_GRD_SNW , 
     I   ZMEAN       , DAYTIME_1HR    , 
     I   ZLT_NOS        , ZLT_SNW       ,
     O   R_STM_BAR_PAR_NOS  , R_STM_BAR_PAR_SNW       , 
     O   RAD_SHRT_NOS, RAD_SHRT_SNW )
!
      USE SIBPRM , ONLY : 
     1    IDX   , 
     1    ICG   , 
     1    IVN   , 
     1    IBD   , 
     1    ITYP 
!
      IMPLICIT NONE
!
! INPUT  
!
      REAL(8),INTENT(IN) :: 
     1  EXTK_SNW       (IDX,ICG,IVN,IBD) , ! 消散係数 
     1  ABSORB_SHRT_SNW(IDX,ICG,IVN,IBD) , ! 下向き短波の植生吸収率 旧RADFAC
     1  REACH_GRD_SNW  (IDX,        IBD)   ! 地面上端に届く割合 旧 P?F 
!
      REAL(8),INTENT(IN) :: 
     1  EXTK_NOS       (IDX,ICG,IVN,IBD) , ! 消散係数 
     1  ABSORB_SHRT_NOS(IDX,ICG,IVN,IBD) , ! 下向き短波の植生吸収率 旧RADFAC
     1  REACH_GRD_NOS  (IDX,        IBD)   ! 地面上端に届く割合 旧 P?F 
!
      REAL(8),INTENT(IN) :: RAD_SHRT_INTO_SIB ( IDX , IVN , IBD ) 
!
      REAL(8),INTENT(IN) ::
     1  CVR_VEG_ms     ( IDX,ICG )        ! 放射、直達可視
!
      REAL(8),INTENT(IN)   :: ZMEAN  (IDX)                
      REAL(8),INTENT(IN)   :: DAYTIME_1HR (IDX)                
      REAL(8),INTENT(IN)   :: ZLT_NOS(IDX,ICG)                
      REAL(8),INTENT(IN)   :: ZLT_SNW(IDX,ICG)                
!
      INTEGER,INTENT(IN)   :: NTYP   (IDX)   
!
      INTEGER,INTENT(IN)   :: MON 
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_01(3)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: R_STM_BAR_PAR_SNW (IDX,ICG)          ! 気孔抵抗
      REAL(8),INTENT(OUT)   :: R_STM_BAR_PAR_NOS (IDX,ICG)          ! 気孔抵抗
!
! OUTPUT 
!
      REAL(8),INTENT(OUT) ::  RAD_SHRT_NOS   (IDX,ICG)   
      REAL(8),INTENT(OUT) ::  RAD_SHRT_SNW   (IDX,ICG)   
!
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD2SIB_CALC_SHRT_MAIN 00/07 POCHAKA ' 
        LFIRST = .FALSE.
      ENDIF 
!
! ====================
! >> 気孔抵抗の計算 <<
! ====================
!
!  気孔抵抗の PAR 依存性 R_STM_BAR_PAR の計算 （ 4 成分はこのために必要 ）
!
      CALL RAD_R_STM (
     I    1              , NMAX_01(3)     , 0        , 
     I    MON            , NTYP           , 
     I    ZMEAN          , CVR_VEG_ms     , 
     I    EXTK_NOS       , REACH_GRD_NOS  , ZLT_NOS  , 
     I    RAD_SHRT_INTO_SIB   , 
     O    R_STM_BAR_PAR_NOS     )
!
      CALL RAD_R_STM (
     I    NMAX_01(2)     , NMAX           , 1        , 
     I    MON            , NTYP           , 
     I    ZMEAN          , CVR_VEG_ms     , 
     I    EXTK_SNW       , REACH_GRD_SNW  , ZLT_SNW  , 
     I    RAD_SHRT_INTO_SIB   , 
     O    R_STM_BAR_PAR_SNW     )
!
! ================================
! >> SiB 入射短波放射の足し込み <<
! ================================
!      セーブすべきは、成分を加算したもの
!
        CALL RAD_SHRT_ADD (
     I    NMAX              , NMAX_01           , 
     I    ABSORB_SHRT_NOS   , ABSORB_SHRT_SNW   , 
     I    RAD_SHRT_INTO_SIB , 
     I    DAYTIME_1HR       , 
     O    RAD_SHRT_NOS      , RAD_SHRT_SNW      ) 
!
      RETURN
      END

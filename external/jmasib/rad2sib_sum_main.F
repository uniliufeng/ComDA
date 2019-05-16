! 放射-> SiB 関連メインプログラム
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/30 保坂征宏 コメント等変更。
!
!
      SUBROUTINE RAD2SIB_SUM_MAIN ( 
     I   NMAX        , NMAX_01        , 
#ifdef MONYOS
     I   DELT_CNP    ,
#endif
     I   NTYP        , MON            , 
     I   CVR_VEG_ms  , CVR_VEG_INV_ms ,   ZLT_NOS , ZLT_SNW , 
     I   ZMEAN       , ZTEMP          , 
     I   TMP_GRSK    , TMP_SNSK       , 
     I   TMP_CNP_NOS , TMP_CNP_SNW    , 
     I   RAD_LONG_INTO_SIB  , 
     M   RAD_SHRT_NOS , RAD_SHRT_SNW  , 
     O   RAD_SHRT_INTO_SNOW_LAYER     , 
     O   THERMK_NOS    ,  
     O   THERMK_SNW    , 
     O   RAD_LONG_NOS_2_ATM , RAD_LONG_SNW_2_ATM ,     
     O   RAD_LONG_NOS , RAD_LONG_SNW )   
!
      USE SIBPRM , ONLY : 
     1    IDX   , 
     1    ICG    
!
      IMPLICIT NONE
!
! INPUT  
!
      REAL(8),INTENT(IN) ::
     1  CVR_VEG_ms     ( IDX,ICG ) ,        ! 放射、直達可視
     1  CVR_VEG_INV_ms ( IDX,ICG )          ! 放射、直達可視
!
      REAL(8),INTENT(IN) ::
     1  ZMEAN  (IDX)               , 
     1  ZTEMP  (IDX) 

      REAL(8),INTENT(IN) ::
     1  ZLT_NOS  (IDX,ICG) ,
     1  ZLT_SNW  (IDX,ICG)

!
      REAL(8),INTENT(IN)  ::
     1  TMP_CNP_NOS(IDX) , 
     1  TMP_CNP_SNW(IDX) , 
     1  TMP_GRSK   (IDX) , 
     1  TMP_SNSK   (IDX)  
!
      INTEGER,INTENT(IN)  ::  NTYP     (IDX)   

!
      INTEGER,INTENT(IN)   :: MON 
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_01(3)
!
#ifdef MONYOS
      REAL(8),INTENT(IN)   :: DELT_CNP
#endif
!
! INPUT 
!
      REAL(8),INTENT(INOUT) ::  RAD_SHRT_NOS  (IDX,ICG)         
      REAL(8),INTENT(INOUT) ::  RAD_SHRT_SNW  (IDX,ICG)             
!
! OUTPUT 
!
!   CG 、雪での単位時間あたり放射吸収量 （吸収が正）
!
      REAL(8),INTENT(OUT)  :: 
     1    RAD_SHRT_INTO_SNOW_LAYER (IDX)       ! 雪層の放射吸収量
!    1    RAD_SUM_SNW       (IDX,ICG)        , ! CG での放射吸収量 
!    1    RAD_SUM_NOS       (IDX,ICG)        
!
      REAL(8),INTENT(OUT)  :: 
     1    RAD_LONG_NOS      (IDX,ICG)        ,
     1    RAD_LONG_SNW      (IDX,ICG)
!
      REAL(8),INTENT(OUT)  :: 
     1    RAD_LONG_NOS_2_ATM      (IDX) 
      REAL(8),INTENT(OUT)  :: 
     1    RAD_LONG_SNW_2_ATM      (IDX)
!
!
!   キャノピーでの  1-射出率
!
      REAL(8),INTENT(OUT)   :: THERMK_SNW  (IDX)       ! 1 - キャノピー射出率
      REAL(8),INTENT(OUT)   :: THERMK_NOS  (IDX)       ! 1 - キャノピー射出率
!
      REAL(8),INTENT(IN)    :: RAD_LONG_INTO_SIB (IDX) 
!
      REAL(8)     :: DUMMY_IDX(IDX)
!
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD2SIB_SUM_MAIN 00/07 POCHAKA ' , MON 
        LFIRST = .FALSE.
      ENDIF 
!
! ==================
! >> 短波放射補正 <<
! ==================
!  一時間平均の天頂角から求めたフラックスに、
!  各ステップでの天頂角を与えることで補正する。
!
      CALL RAD_SHRT_COR (
     I  NMAX    , NMAX_01 , 
     I  ZMEAN   , ZTEMP   , 
     M  RAD_SHRT_NOS      ,  RAD_SHRT_SNW       )  
!
! ========================
! >> 放射の各成分の加算 <<
! ========================
!
!    短波と長波を加えて、キャノピー・地面が受け取るエネルギーの和を求める。
!    長波については、キャノピーでの射出率（吸収率）も計算する。
!
!    注意。
!      ・雪層に入る、RAD_SHRT_SNW(*,2) については RAD_SUM_SNW に加算しない。
!        かわりに RAD_SHRT_intoSNW に格納する。
!      ・長波についてはここで RAD 格子から SiB 格子に変換する。
!
      CALL RAD_SUM_ALL (
     I  1              , NMAX_01(3)     , 0          ,  
#ifdef MONYOS
     I  DELT_CNP   , 
#endif 
     I  NTYP           , 
     I  TMP_CNP_NOS    , TMP_GRSK       ,  
     I  CVR_VEG_ms     , CVR_VEG_INV_ms , 
     I  RAD_LONG_INTO_SIB               , ZLT_NOS , 
     I  RAD_SHRT_NOS   , 
     O  THERMK_NOS     , 
     O  DUMMY_IDX      , 
     O  RAD_LONG_NOS_2_ATM , 
     O  RAD_LONG_NOS )
!
      CALL RAD_SUM_ALL (
     I  NMAX_01(2)     , NMAX           , 1          ,  
#ifdef MONYOS
     I  DELT_CNP   , 
#endif 
     I  NTYP           , 
     I  TMP_CNP_SNW    , TMP_SNSK       ,  
     I  CVR_VEG_ms     , CVR_VEG_INV_ms , 
     I  RAD_LONG_INTO_SIB               , ZLT_SNW , 
     M  RAD_SHRT_SNW   , 
     O  THERMK_SNW     , 
     O  RAD_SHRT_INTO_SNOW_LAYER        , 
     O  RAD_LONG_SNW_2_ATM ,  
     O  RAD_LONG_SNW  )
!
      RETURN
      END
!=======================================================
      SUBROUTINE RAD2SIB_EXTEND (
     I   NMAX      , 
     I   NTYP      , MON        ) 
!
! 展開する
!
      USE SIBPRM , ONLY :
     1   IDX    ,
     1   ITYP   , 
     1   ICG  
!
      USE SIBVEG , ONLY :
     1           CVR_VEG_mo                ! 被覆率 (ITYP,IMON,ICG)
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NMAX         ! 陸格子数
      INTEGER,INTENT(IN)   :: NTYP (IDX)  ! 植生タイプ
      INTEGER,INTENT(IN)   :: MON          ! 月
!
! WORK  
!
      INTEGER :: I
      INTEGER :: IVEG
!
      LOGICAl,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD2SIB EXTEND 99/07 POCHAKA ' 
     1             , ' DUMMY ' , NMAX , MON , NTYP 
        LFIRST = .FALSE. 
      ELSE
        WRITE(6,*) 'RAD2SIB EXTEND HAS ALREADY CALLED'
        STOP 999
      ENDIF
!
      DO 1000 IVEG=1,ICG
      DO 1000 I=1,NMAX
!       CVR_VEG_ms    (I,IVEG) = CVR_VEG_mo(NTYP(I),MON,IVEG)
!       CVR_VEG_INV_ms(I,IVEG) = 1/CVR_VEG_mo(NTYP(I),MON,IVEG)
 1000 CONTINUE
!
      RETURN
      END

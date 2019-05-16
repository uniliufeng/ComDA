! キャノピーの更新後の値を、物理格子に戻す
! 履歴 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SETBACK_SIBGRD_SNOW (
     I    JL     , ILIST , NMAX    , 
     I    TMP_SNOW_SIB      , 
     I    WTR_SNOW_SIB      , GLA_SNOW_SIB    , 
     I    RHO_SNOW_INV_SIB  , 
     I    AGE_SNOW_SIB      , 
     I    INFO_SNOW_SIB     , 
     I    ENG_SNOW_BUCKET_SIB    , H2O_SNOW_BUCKET_SIB    , 
     O    TMP_SNOW_PHY      , 
     O    WTR_SNOW_PHY      , GLA_SNOW_PHY    , 
     O    RHO_SNOW_INV_PHY  , 
     O    AGE_SNOW_PHY      , 
     O    INFO_SNOW_PHY     , 
     O    ENG_SNOW_BUCKET_PHY    , H2O_SNOW_BUCKET_PHY    ) 
!
      USE PRM , ONLY :
     1   IJPHY  ,
     1   JLPHY  ,
     1   ISPT
!
      USE SIBPRM , ONLY : 
     1   IDX ,   
     1   ISN

      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: JL
      INTEGER,INTENT(IN)  :: NMAX
      INTEGER,INTENT(IN)  :: ILIST              (IDX)
!
      REAL(8),INTENT(IN)  :: TMP_SNOW_SIB       (IDX,ISN)
      REAL(8),INTENT(IN)  :: WTR_SNOW_SIB       (IDX,ISN)
      REAL(8),INTENT(IN)  :: GLA_SNOW_SIB       (IDX,ISN)
      REAL(8),INTENT(IN)  :: RHO_SNOW_INV_SIB   (IDX,ISN) 
      REAL(8),INTENT(IN)  :: AGE_SNOW_SIB       (IDX) 
      INTEGER,INTENT(IN)  :: INFO_SNOW_SIB      (IDX) 
      REAL(8),INTENT(IN)  :: ENG_SNOW_BUCKET_SIB(IDX) 
      REAL(8),INTENT(IN)  :: H2O_SNOW_BUCKET_SIB(IDX) 
!
! OUTPUT
!
      REAL(8),INTENT(INOUT) :: TMP_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)      
      REAL(8),INTENT(INOUT) :: WTR_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(INOUT) :: GLA_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(INOUT) :: RHO_SNOW_INV_PHY   (ISPT*IJPHY,JLPHY,ISN) 
      REAL(8),INTENT(INOUT) :: AGE_SNOW_PHY       (ISPT*IJPHY,JLPHY) 
      INTEGER,INTENT(INOUT) :: INFO_SNOW_PHY      (ISPT*IJPHY,JLPHY) 
      REAL(8),INTENT(INOUT) :: ENG_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY) 
      REAL(8),INTENT(INOUT) :: H2O_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY) 
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: J
      INTEGER   :: K
!
      IF ( NMAX.LT.1 ) RETURN
!
      DO 1000 K=1,ISN
      DO 1000 I=1,NMAX
        J = ILIST(I)
        TMP_SNOW_PHY      (J,JL,K) = TMP_SNOW_SIB     (I,K)     
        WTR_SNOW_PHY      (J,JL,K) = WTR_SNOW_SIB     (I,K)     
        GLA_SNOW_PHY      (J,JL,K) = GLA_SNOW_SIB     (I,K)     
        RHO_SNOW_INV_PHY  (J,JL,K) = RHO_SNOW_INV_SIB (I,K)     
 1000 CONTINUE
!
      DO 2000 I=1,NMAX
        J = ILIST(I)
        AGE_SNOW_PHY        (J,JL) = AGE_SNOW_SIB        (I)     
        INFO_SNOW_PHY       (J,JL) = INFO_SNOW_SIB       (I)     
        ENG_SNOW_BUCKET_PHY (J,JL) = ENG_SNOW_BUCKET_SIB (I)     
        H2O_SNOW_BUCKET_PHY (J,JL) = H2O_SNOW_BUCKET_SIB (I)     
 2000 CONTINUE
!
      RETURN
      END

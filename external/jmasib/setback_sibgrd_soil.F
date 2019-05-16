! キャノピーの更新後の値を、物理格子に戻す
! 履歴 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SETBACK_SIBGRD_SOIL (
     I    JL     , ILIST , NMAX    , 
     I    TMP_SOIL_SIB      , TMP_SNSL_SIB       , 
     I    FR_WTR_SOIL_SIB   , FR_GLA_SOIL_SIB    , 
     I    FR_WTR_SNSL_SIB   , FR_GLA_SNSL_SIB    , 
     I    TMP_SOIL4_SIB     ,
     O    TMP_SOIL_PHY      , TMP_SNSL_PHY       , 
     O    FR_WTR_SOIL_PHY   , FR_GLA_SOIL_PHY    , 
     O    FR_WTR_SNSL_PHY   , FR_GLA_SNSL_PHY    ,
     O    TMP_SOIL4_PHY )
!
      USE PRM , ONLY :
     1   IJPHY  ,
     1   JLPHY  ,
     1   ISPT
!
      USE SIBPRM , ONLY : 
     1   IDX ,   
     1   IDP 

      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: JL
      INTEGER,INTENT(IN)  :: NMAX
      INTEGER,INTENT(IN)  :: ILIST              (IDX)
!
      REAL(8),INTENT(IN)  :: TMP_SOIL_SIB       (IDX,IDP)
      REAL(8),INTENT(IN)  :: FR_WTR_SOIL_SIB    (IDX,IDP)
      REAL(8),INTENT(IN)  :: FR_GLA_SOIL_SIB    (IDX,IDP)
      REAL(8),INTENT(IN)  :: TMP_SNSL_SIB       (IDX,IDP)
      REAL(8),INTENT(IN)  :: FR_WTR_SNSL_SIB    (IDX,IDP)
      REAL(8),INTENT(IN)  :: FR_GLA_SNSL_SIB    (IDX,IDP)
!
      REAL(8),INTENT(IN)  :: TMP_SOIL4_SIB      (IDX)
!
! OUTPUT
!
      REAL(8),INTENT(INOUT) :: TMP_SOIL_PHY     (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT) :: FR_WTR_SOIL_PHY  (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT) :: FR_GLA_SOIL_PHY  (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT) :: TMP_SNSL_PHY     (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT) :: FR_WTR_SNSL_PHY  (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT) :: FR_GLA_SNSL_PHY  (ISPT*IJPHY,JLPHY,IDP)
!
      REAL(8),INTENT(INOUT) :: TMP_SOIL4_PHY    (ISPT*IJPHY,JLPHY)
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: J
      INTEGER   :: K
!
      IF ( NMAX.LT.1 ) RETURN
!
      DO 1000 K=1,IDP
      DO 1000 I=1,NMAX
        J = ILIST(I)
        TMP_SOIL_PHY      (J,JL,K) = TMP_SOIL_SIB   (I,K)     
        FR_WTR_SOIL_PHY   (J,JL,K) = FR_WTR_SOIL_SIB(I,K)      
        FR_GLA_SOIL_PHY   (J,JL,K) = FR_GLA_SOIL_SIB(I,K)      
!
        TMP_SNSL_PHY      (J,JL,K) = TMP_SNSL_SIB   (I,K)     
        FR_WTR_SNSL_PHY   (J,JL,K) = FR_WTR_SNSL_SIB(I,K)      
        FR_GLA_SNSL_PHY   (J,JL,K) = FR_GLA_SNSL_SIB(I,K)      
 1000 CONTINUE
!
      DO I=1,NMAX
        J = ILIST(I)
        TMP_SOIL4_PHY   (J,JL) = TMP_SOIL4_SIB(I) 
      ENDDO
!
      RETURN
      END

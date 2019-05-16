! キャノピーの更新後の値を、物理格子に戻す
! 履歴 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SETBACK_SIBGRD_CANOPY (
     I   JL     , ILIST , NMAX    , 
     I   TMP_CNP_NOS_SIB , FR_WTR_CNP_NOS_SIB , FR_GLA_CNP_NOS_SIB , 
     I   TMP_CNP_SNW_SIB , FR_WTR_CNP_SNW_SIB , FR_GLA_CNP_SNW_SIB , 
     I   TMP_GRSK_SIB    , FR_WTR_GRS_SIB     , FR_GLA_GRS_SIB     , 
     I   TMP_SNSK_SIB    , 
     O   TMP_CNP_NOS_PHY , FR_WTR_CNP_NOS_PHY , FR_GLA_CNP_NOS_PHY , 
     O   TMP_CNP_SNW_PHY , FR_WTR_CNP_SNW_PHY , FR_GLA_CNP_SNW_PHY , 
     O   TMP_GRSK_PHY    , FR_WTR_GRS_PHY     , FR_GLA_GRS_PHY     , 
     O   TMP_SNSK_PHY       ) 
!
      USE PRM , ONLY :
     1   IJPHY ,
     1   JLPHY ,
     1   ISPT 
!
      USE SIBPRM , ONLY : 
     1   IDX   
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: JL
      INTEGER,INTENT(IN)  :: NMAX
      INTEGER,INTENT(IN)  :: ILIST              (IDX)
!
      REAL(8),INTENT(IN)  :: TMP_CNP_NOS_SIB    (IDX)
      REAL(8),INTENT(IN)  :: FR_WTR_CNP_NOS_SIB (IDX)    
      REAL(8),INTENT(IN)  :: FR_GLA_CNP_NOS_SIB (IDX)    
      REAL(8),INTENT(IN)  :: TMP_CNP_SNW_SIB    (IDX)
      REAL(8),INTENT(IN)  :: FR_WTR_CNP_SNW_SIB (IDX)    
      REAL(8),INTENT(IN)  :: FR_GLA_CNP_SNW_SIB (IDX)    
      REAL(8),INTENT(IN)  :: TMP_GRSK_SIB       (IDX)
      REAL(8),INTENT(IN)  :: FR_WTR_GRS_SIB     (IDX)   
      REAL(8),INTENT(IN)  :: FR_GLA_GRS_SIB     (IDX)   
      REAL(8),INTENT(IN)  :: TMP_SNSK_SIB       (IDX)
!
! OUTPUT
!
      REAL(8),INTENT(INOUT) :: TMP_CNP_NOS_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT) :: FR_WTR_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT) :: FR_GLA_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT) :: TMP_CNP_SNW_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT) :: FR_WTR_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT) :: FR_GLA_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT) :: TMP_GRSK_PHY       (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT) :: FR_WTR_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(INOUT) :: FR_GLA_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(INOUT) :: TMP_SNSK_PHY       (ISPT*IJPHY,JLPHY)
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: J
!
      IF ( NMAX.LT.1 ) RETURN
!
! ------------
! > 雪なし域 <
! ------------
!
      DO 2000 I=1,NMAX
        J = ILIST(I)
        TMP_CNP_NOS_PHY   (J,JL) = TMP_CNP_NOS_SIB   (I)     
        FR_WTR_CNP_NOS_PHY(J,JL) = FR_WTR_CNP_NOS_SIB(I)      
        FR_GLA_CNP_NOS_PHY(J,JL) = FR_GLA_CNP_NOS_SIB(I)     
        TMP_GRSK_PHY      (J,JL) = TMP_GRSK_SIB      (I)     
        FR_WTR_GRS_PHY    (J,JL) = FR_WTR_GRS_SIB    (I)     
        FR_GLA_GRS_PHY    (J,JL) = FR_GLA_GRS_SIB    (I)     
 2000 CONTINUE
!
! ------------
! > 雪あり域 <
! ------------
!
      DO 3000 I=1,NMAX   
        J = ILIST(I)
        TMP_CNP_SNW_PHY   (J,JL) = TMP_CNP_SNW_SIB   (I)     
        FR_WTR_CNP_SNW_PHY(J,JL) = FR_WTR_CNP_SNW_SIB(I)     
        FR_GLA_CNP_SNW_PHY(J,JL) = FR_GLA_CNP_SNW_SIB(I)     
        TMP_SNSK_PHY      (J,JL) = TMP_SNSK_SIB      (I)     
 3000 CONTINUE
!
      RETURN
      END

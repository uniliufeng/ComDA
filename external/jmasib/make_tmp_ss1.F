      SUBROUTINE MAKE_TMP_SS1 ( 
     I   NMAX_01          , NMAX     ,  
     I   TMP_SOIL         , TMP_SNOW , 
     O   TMP_SOIL1_TMP    , TMP_SNOW1_TMP )
!
      USE SIBPRM , ONLY :
     1   IDX
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NMAX_01(3)
!     INTEGER,INTENT(IN)    :: ILIST   (IDX)   
!
      REAL(8),INTENT(IN)    :: TMP_SOIL(IDX)
      REAL(8),INTENT(IN)    :: TMP_SNOW(IDX)
!
      REAL(8),INTENT(OUT)   :: TMP_SOIL1_TMP(IDX)
      REAL(8),INTENT(OUT)   :: TMP_SNOW1_TMP(IDX)
!
      INTEGER :: I
!
      DO I=1,NMAX_01(3) 
        TMP_SOIL1_TMP(I) = TMP_SOIL(I)
      ENDDO
      DO I=NMAX_01(2),NMAX 
        TMP_SNOW1_TMP(I) = TMP_SNOW(I)
      ENDDO
!
      RETURN
      END


      SUBROUTINE CANOPY_COPY ( 
     I    NSTR        , NEND         , INSFLAG      ,  
     I    TMP_CNP_IN  , TMP_SKIN_IN  , TMP_REF_IN   , TMP_SS1_IN  ,   
     I    E_REF_IN    ,    
     O    TMP_CNP_OUT , TMP_SKIN_OUT , TMP_REF_OUT  , TMP_SS1_OUT ,
     O    E_REF_OUT  ) 
!
      USE SIBPRM , ONLY : 
     1   IDX 
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR 
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: INSFLAG
!
      REAL(8),INTENT(IN)   :: TMP_CNP_IN      ( IDX )   
      REAL(8),INTENT(IN)   :: TMP_SKIN_IN     ( IDX ) 
      REAL(8),INTENT(IN)   :: TMP_REF_IN      ( IDX ) 
      REAL(8),INTENT(IN)   :: TMP_SS1_IN      ( IDX ) 
!
      REAL(8),INTENT(IN)   :: E_REF_IN        ( IDX ) 
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: TMP_CNP_OUT     ( IDX )   
      REAL(8),INTENT(OUT)  :: TMP_SKIN_OUT    ( IDX ) 
      REAL(8),INTENT(OUT)  :: TMP_REF_OUT     ( IDX ) 
      REAL(8),INTENT(OUT)  :: TMP_SS1_OUT     ( IDX ) 
!
      REAL(8),INTENT(INOUT)  :: E_REF_OUT       ( IDX ) 
!
! WORK
!
      INTEGER       :: I
!
! WORK
!
      DO I=NSTR,NEND
        TMP_CNP_OUT    (I) = TMP_CNP_IN    (I)
        TMP_SKIN_OUT   (I) = TMP_SKIN_IN   (I)
        TMP_REF_OUT    (I) = TMP_REF_IN    (I)
        TMP_SS1_OUT    (I) = TMP_SS1_IN    (I)
        E_REF_OUT      (I) = E_REF_IN      (I)
      ENDDO
!
      RETURN
      END
 

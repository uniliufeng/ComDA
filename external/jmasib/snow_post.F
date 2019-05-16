!===========================================================
      SUBROUTINE SNOW_POST (
     I    NMAX          , 
     I    GLA_SNOW_NEXT , CVR_SNW_INV_NEXT , 
     I    RHO_SNOW_INV_NEXT , 
     O    DEP_SNOW_NEXT  )
!
      USE SIBPRM , ONLY :
     1    IDX  , 
     1    ISN  
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: NMAX 
      REAL(8),INTENT(IN)  :: CVR_SNW_INV_NEXT(IDX)
      REAL(8),INTENT(IN)  :: GLA_SNOW_NEXT    (IDX,ISN)
      REAL(8),INTENT(IN)  :: RHO_SNOW_INV_NEXT(IDX,ISN)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT) :: DEP_SNOW_NEXT    (IDX,ISN+1)
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: K
      LOGICAL,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_POST 99/07/xx POCHAKA'
        LFIRST = .FALSE.
      ENDIF
!
! ------------
! > á‚Ì[‚³ <
! ------------
!
      DEP_SNOW_NEXT(:,ISN+1) = 0. 
!
      DO K=1,ISN
      DO I=1,NMAX 
        DEP_SNOW_NEXT(I,K) = GLA_SNOW_NEXT(I,K) * RHO_SNOW_INV_NEXT(I,K) 
     1                       * CVR_SNW_INV_NEXT(I)       
        DEP_SNOW_NEXT(I,ISN+1) = DEP_SNOW_NEXT(I,ISN+1) 
     1                           + DEP_SNOW_NEXT(I,K)   
      ENDDO
      ENDDO
!
      RETURN
      END


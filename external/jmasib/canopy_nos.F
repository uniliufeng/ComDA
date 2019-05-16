! óáínÇ≈ÇÃèˆî≠Ç∆è∏âÿÇÃäÑçáÇÃåàíË
! 99/10/01 ï€ç‚ê™çG êVãK
!
      SUBROUTINE CANOPY_NOS (
     I    NSTR             , NEND        , 
     I    TMP_SOIL         , FR_WTR_SOIL , FR_GLA_SOIL  , 
     O    FR_WBH_BS      ) 
!
      USE SIBPRM , ONLY :
     1    IDX
!
      USE SIBCON , ONLY :
     1    TMP_FREZ_c
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NSTR 
      INTEGER,INTENT(IN)   :: NEND
      REAL(8),INTENT(IN)   :: TMP_SOIL   (IDX)
      REAL(8),INTENT(IN)   :: FR_WTR_SOIL(IDX)
      REAL(8),INTENT(IN)   :: FR_GLA_SOIL(IDX)
!
      REAL(8),INTENT(OUT)  :: FR_WBH_BS   (IDX)
!
      INTEGER    :: I
!
      DO I=NSTR,NEND
        IF      ( TMP_SOIL(I) .GT. TMP_FREZ_c ) THEN
          FR_WBH_BS(I) = 1. 
        ELSEIF  ( TMP_SOIL(I) .LT. TMP_FREZ_c ) THEN
          FR_WBH_BS(I) = 0. 
        ELSE
          FR_WBH_BS(I) = FR_WTR_SOIL(I) 
     1                   / ( FR_WTR_SOIL(I) + FR_GLA_SOIL(I) )      
        ENDIF
      ENDDO
!
      RETURN
      END



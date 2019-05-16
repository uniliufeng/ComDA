      SUBROUTINE PRESIB_CAP_SSTOP (
     I   NSTR_NOS        , NEND_NOS      ,
     I   NSTR_SNW        , NEND_SNW      , 
     I   NTYP            , 
     I   CVR_SNW_INV     , 
     I   TMP_SOIL        , FR_WTR_SOIL     , FR_GLA_SOIL   , 
     I   TMP_SNSL        , FR_WTR_SNSL     , FR_GLA_SNSL   , 
     I   WTR_SNOW        , GLA_SNOW      , 
     O   CAP_SOIL_TOP    , CAP_SNOW_TOP  ,
     O   ENG_SOIL_TOP  )
!
      USE SIBPRM , ONLY : 
     1   IDX 
!
      USE SIBVEG , ONLY : 
     1           POR_SOIL_co            ,  ! ãÛåÑó¶   (ITYP)  
     1           DEP_SOIL_co               ! ìyèÎëwÇÃê[Ç≥    (ITYP,IDP)
!
      USE SIBCON , ONLY : 
     1    CAP_SOIL_BY_VOL_c ,                     ! CSOIL
     1    CAP_WTR_c  ,                            ! CW 
     1    CAP_GLA_c  ,                            ! CICE 
     1    HL_W2G_c   ,                            ! 
     1    TMP_FREZ_c ,                            ! 
     1    RHO_WTR_c 
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NSTR_NOS
      INTEGER,INTENT(IN)    :: NEND_NOS
      INTEGER,INTENT(IN)    :: NSTR_SNW
      INTEGER,INTENT(IN)    :: NEND_SNW
      INTEGER,INTENT(IN)    :: NTYP        (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_SOIL    (IDX)    
      REAL(8),INTENT(IN)    :: TMP_SNSL    (IDX)    
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL (IDX)    
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL (IDX)    
      REAL(8),INTENT(IN)    :: FR_WTR_SNSL (IDX)    
      REAL(8),INTENT(IN)    :: FR_GLA_SNSL (IDX)    
      REAL(8),INTENT(IN)    :: WTR_SNOW    (IDX)    
      REAL(8),INTENT(IN)    :: GLA_SNOW    (IDX)    
      REAL(8),INTENT(IN)    :: CVR_SNW_INV (IDX)    
!
      REAL(8),INTENT(OUT)   :: CAP_SOIL_TOP(IDX)
      REAL(8),INTENT(OUT)   :: CAP_SNOW_TOP(IDX)
      REAL(8),INTENT(OUT)   :: ENG_SOIL_TOP(IDX)
!
      INTEGER   :: I
      INTEGER   :: J
!
! ê·ñ≥ÇµìyèÎ
!
      DO I=NSTR_NOS,NEND_NOS 
        J = NTYP(I) 
        CAP_SOIL_TOP (I) 
     1      = DEP_SOIL_co(J,1) 
     1           * ( 1 - POR_SOIL_co(J) ) 
     1           * CAP_SOIL_BY_VOL_c   
     1      + DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_WTR_SOIL(I)    
     1           * CAP_WTR_c   
     1      + DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_GLA_SOIL(I)    
     1           * CAP_GLA_c   
!
        ENG_SOIL_TOP(I) 
     1      =   CAP_SOIL_TOP(I) * ( TMP_SOIL(I) - TMP_FREZ_c )
     1        - DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_GLA_SOIL(I)    
     1           * HL_W2G_c 
      ENDDO
!
      DO I=NSTR_SNW,NEND_SNW
        J = NTYP(I) 
!
! ê·Ç†ÇËìyèÎ
!
        CAP_SOIL_TOP (I) 
     1      = DEP_SOIL_co(J,1) 
     1           * ( 1 - POR_SOIL_co(J) ) 
     1           * CAP_SOIL_BY_VOL_c   
     1      + DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_WTR_SNSL(I)    
     1           * CAP_WTR_c   
     1      + DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_GLA_SNSL(I)    
     1           * CAP_GLA_c   
!
        ENG_SOIL_TOP(I) 
     1      =   CAP_SOIL_TOP(I) * ( TMP_SNSL(I) - TMP_FREZ_c )
     1        - DEP_SOIL_co(J,1) 
     1           * POR_SOIL_co(J) * RHO_WTR_c * FR_GLA_SNSL(I)    
     1           * HL_W2G_c 
!
! ê·
!
        CAP_SNOW_TOP(I) 
     1     = (   CAP_GLA_c  * GLA_SNOW(I)
     1         + CAP_WTR_c  * WTR_SNOW(I) ) * CVR_SNW_INV(I) 
      ENDDO
!
      RETURN
      END

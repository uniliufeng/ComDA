! 土壌の積分  （水）
!
      SUBROUTINE SOIL_INTG (
     I   NMAX               , 
     I   FACT_WGT_H2O_INV    , CAP_SOLID         ,
     I   INFO_SNOW_NEXT     , 
     I   H2O_SOIL_NEXT      , H2O_SNSL_NEXT     ,
     I   ENG_SOIL_NEXT      , ENG_SNSL_NEXT     , 
     O   TMP_SOIL_NEXT      , TMP_SNSL_NEXT     , 
     O   FR_WTR_SOIL_NEXT   , FR_WTR_SNSL_NEXT  , 
     O   FR_GLA_SOIL_NEXT   , FR_GLA_SNSL_NEXT ) 
!
      USE SIBPRM , ONLY :
     1   IDX   , 
     1   IDP   
!
      USE SIBCON , ONLY :
     1   RHO_WTR_c   , 
     1   CAP_GLA_c   , 
     1   CAP_WTR_c   , 
     1   HL_W2G_c
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX
!
!
      REAL(8),INTENT(IN)    :: FACT_WGT_H2O_INV  (IDX,IDP)    ! 空隙量x水の比重
      REAL(8),INTENT(IN)    :: CAP_SOLID    (IDX,IDP)    ! 固体土の比熱
!
!
      REAL(8),INTENT(IN)    :: ENG_SOIL_NEXT (IDX,IDP)   
      REAL(8),INTENT(IN)    :: H2O_SOIL_NEXT (IDX,IDP)   
      REAL(8),INTENT(IN)    :: ENG_SNSL_NEXT (IDX,IDP) 
      REAL(8),INTENT(IN)    :: H2O_SNSL_NEXT (IDX,IDP) 
!
      INTEGER,INTENT(IN)    :: INFO_SNOW_NEXT(IDX) 
!
! INOUT
!
      REAL(8),INTENT(OUT)   :: TMP_SOIL_NEXT   (IDX,IDP) 
      REAL(8),INTENT(OUT)   :: TMP_SNSL_NEXT   (IDX,IDP) 
      REAL(8),INTENT(OUT)   :: FR_WTR_SOIL_NEXT(IDX,IDP) 
      REAL(8),INTENT(OUT)   :: FR_GLA_SOIL_NEXT(IDX,IDP) 
      REAL(8),INTENT(OUT)   :: FR_WTR_SNSL_NEXT(IDX,IDP) 
      REAL(8),INTENT(OUT)   :: FR_GLA_SNSL_NEXT(IDX,IDP) 
!
! WORK
!
      INTEGER  :: I
      INTEGER  :: K
      INTEGER  :: INFO_NOS(IDX)
      INTEGER  :: INFO_SNW(IDX)
!     INTEGER,PARAMETER  :: IONE_IDX(IDX) = 1
!
      LOGICAL,SAVE :: LFIRST
!
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_INTG 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  ======================
!  >> 層変化・温度変化 <<
!  ======================
!
      DO I=1,NMAX
        IF ( INFO_SNOW_NEXT(I).LE.0  ) THEN    ! 雪無し
          INFO_NOS(I) = 1 
        ELSE
          INFO_NOS(I) = 0 
        ENDIF
        IF ( INFO_SNOW_NEXT(I).GE.0 ) THEN     ! 雪あり
          INFO_SNW(I) = 1 
        ELSE
          INFO_SNW(I) = 0 
        ENDIF
      ENDDO
!
      TMP_SOIL_NEXT   (:,:) = 0.
      FR_WTR_SOIL_NEXT(:,:) = 0.
      FR_GLA_SOIL_NEXT(:,:) = 0.
      TMP_SNSL_NEXT   (:,:) = 0.
      FR_WTR_SNSL_NEXT(:,:) = 0.
      FR_GLA_SNSL_NEXT(:,:) = 0.
!
      DO K=1,IDP 
!
!  -------------------------
!  > 土壌第 k 層、雪無し域 <
!  -------------------------
!
!
      CALL CALC_SOIL_TMP (
     I   NMAX                  , 
     I   INFO_NOS              , 
     I   CAP_SOLID(1,K)        , FACT_WGT_H2O_INV(1,K)  , 
     I   ENG_SOIL_NEXT(1,K)    , H2O_SOIL_NEXT   (1,K)  , 
     O   TMP_SOIL_NEXT(1,K)    , 
     O   FR_WTR_SOIL_NEXT(1,K) , FR_GLA_SOIL_NEXT(1,K)  )
!
!  -------------------------
!  > 土壌第 k 層、雪あり域 <
!  -------------------------
!
!
      CALL CALC_SOIL_TMP (
     I   NMAX                  ,
     I   INFO_SNW              , 
     I   CAP_SOLID(1,K)        , FACT_WGT_H2O_INV(1,K)  , 
     I   ENG_SNSL_NEXT(1,K)    , H2O_SNSL_NEXT   (1,K)  , 
     O   TMP_SNSL_NEXT(1,K)    , 
     O   FR_WTR_SNSL_NEXT(1,K) , FR_GLA_SNSL_NEXT(1,K)  )
!
      ENDDO 
!
      RETURN
      END
!==================================================================
      SUBROUTINE CALC_SOIL_TMP (
     I   NMAX     , 
     I   INFO     , 
     I   CAP_SOLID, FACT_WGT_H2O_INV   , 
     I   ENG      , H2O       , 
     O   TMP      , FR_WTR    , FR_GLA  )
!  
      USE SIBPRM , ONLY :
     1     IDX 
!
      USE SIBCON , ONLY :
     1     HL_W2G_c , 
     1     TMP_FREZ_c ,
     1     CAP_WTR_c ,
     1     CAP_GLA_c 
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: INFO           (IDX) 
      REAL(8),INTENT(IN)    :: CAP_SOLID      (IDX)
      REAL(8),INTENT(IN)    :: FACT_WGT_H2O_INV(IDX)
      REAL(8),INTENT(IN)    :: ENG(IDX)
      REAL(8),INTENT(IN)    :: H2O(IDX)
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT)   :: TMP   (IDX)
      REAL(8),INTENT(INOUT)   :: FR_WTR(IDX)
      REAL(8),INTENT(INOUT)   :: FR_GLA(IDX)
!
! WORK
!
      INTEGER :: I
      real(8) :: CAP
      real(8) :: WTR_1
      real(8) :: GLA_1
      REAL(8) :: ENG_X(IDX) 
      REAL(8) :: FR   (IDX) 
      REAL(8) :: BUNBO 
      REAL(8),PARAMETER  :: DZERO = 0.D0
      REAL(8),PARAMETER  :: SMALL = 1.D-20
!
      DO 1000 I=1,NMAX
        IF ( INFO(I) .EQ. 1 ) THEN
          ENG_X (I) = ENG(I) + HL_W2G_c * H2O(I) 
        ENDIF
 1000 CONTINUE  
!
      DO 2000 I=1,NMAX
        IF ( INFO(I) .EQ. 1 ) THEN
          IF     ( ENG(I) .GE. 0. ) THEN                      
            FR(I) = 0. 
          ELSEIF ( ENG_X(I) .LE. 0. ) THEN              
            FR(I) = 1. 
          ELSE 
            BUNBO = MAX ( ENG_X(I)-ENG(I) , SMALL )
            FR(I) = - ENG(I) / BUNBO                  ! 氷の割合
          ENDIF
        ENDIF
 2000 CONTINUE  
!
      DO 3000 I=1,NMAX 
        IF ( INFO(I) .EQ. 1 ) THEN
          WTR_1 = H2O(I) * (1.-FR(I))
          GLA_1 = H2O(I) * FR(I)
          CAP   =   CAP_SOLID(I)
     1            + CAP_WTR_c  * WTR_1        
     1            + CAP_GLA_c  * GLA_1        
!
          TMP(I) = TMP_FREZ_c 
     1             + (   MAX ( ENG  (I) , DZERO )  
     1                 + MIN ( ENG_X(I) , DZERO ) ) / CAP 
          FR_WTR(I) = WTR_1 * FACT_WGT_H2O_INV(I)
          FR_GLA(I) = GLA_1 * FACT_WGT_H2O_INV(I)
        ENDIF
 3000 CONTINUE  
!
      RETURN
      END

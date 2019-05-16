      SUBROUTINE SOIL_PRE_1 (
     I    NMAX            , NTYP       , 
     O    FACT_WGT_H2O     , 
     O    FACT_WGT_H2O_INV ,     
     O    CAP_SOLID      )
!
!
      USE SIBPRM , ONLY : 
     1    IDX     , 
     1    IDP     , 
     1    ITYP 
!
      USE SIBCON , ONLY : 
     1    CAP_SOIL_BY_VOL_c ,                            ! CSOIL
     1    CAP_WTR_c  ,                            ! CW 
     1    CAP_GLA_c  ,                            ! CICE 
     1    RHO_WTR_c  , 
     1    RHO_WTR_INV_c   
!
      USE SIBVEG , ONLY : 
     1           POR_SOIL_co            ,  ! ‹óŒ„—¦   (ITYP)  
     1           DEP_SOIL_co            ,  ! “yë‘w‚Ì[‚³    (ITYP,IDP)
     1           POR_SOIL_INV_co        ,  ! ‹óŒ„—¦   (ITYP)  
     1           DEP_SOIL_INV_co           ! “yë‘w‚Ì[‚³    (ITYP,IDP)
!
      IMPLICIT NONE 
! INPUT 
      INTEGER,INTENT(IN)    :: NTYP     (IDX)     ! A¶ƒ^ƒCƒv
      INTEGER,INTENT(IN)    :: NMAX  
! OUTPUT 
      REAL(8),INTENT(OUT)   :: FACT_WGT_H2O  (IDX,IDP)    ! ‹óŒ„—Êx…‚Ì”äd
      REAL(8),INTENT(OUT)   :: FACT_WGT_H2O_INV (IDX,IDP) 
      REAL(8),INTENT(OUT)   :: CAP_SOLID    (IDX,IDP)     ! ŒÅ‘Ì“y‚Ì”ä”M
!
      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
!
#ifdef SIB_DEBUG
      REAL(8):: TEST
#endif
!
      DO K=1,3 
        DO I=1,NMAX
          J = NTYP(I)  
          FACT_WGT_H2O(I,K)     = DEP_SOIL_co(J,K) * RHO_WTR_c 
     1                            * POR_SOIL_co(J)    
          FACT_WGT_H2O_INV(I,K) = DEP_SOIL_INV_co(J,K) * RHO_WTR_INV_c 
     1                            * POR_SOIL_INV_co(J)    
!
#ifdef SIB_DEBUG
          TEST = FACT_WGT_H2O(I, K) * FACT_WGT_H2O_INV(I, K) - 1.0d0
!         IF ( ABS(TEST) > EPSILON(TEST)) THEN
          IF ( ABS(TEST) > 1.D-15 ) THEN
            write(6,*) 'soil_pre_1: bad precision I K ' , I , K , test,
     1                  FACT_WGT_H2O(I, K) , FACT_WGT_H2O_INV(I, K) ,
     1                  DEP_SOIL_co(J,K) , RHO_WTR_c ,
     1                             POR_SOIL_co(J)    ,
     1                  DEP_SOIL_INV_co(J,K) , RHO_WTR_INV_c ,
     1                             POR_SOIL_INV_co(J)    
            stop 999
          ENDIF
#endif
!
          CAP_SOLID  (I,K) = DEP_SOIL_co(J,K) 
     1                            * ( 1 - POR_SOIL_co(J) ) 
     1                            * CAP_SOIL_BY_VOL_c   
!
        ENDDO
      ENDDO
!
      RETURN
      END

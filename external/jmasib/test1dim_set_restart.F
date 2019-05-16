!==============================================================
! ó\ïÒïœêîÇÃílÇÃê›íË
      SUBROUTINE TEST1D_SET_SIB_RESTART (
     O   TMP_CNP_NOS_PHY    , FR_WTR_CNP_NOS_PHY , FR_GLA_CNP_NOS_PHY , 
     O   TMP_CNP_SNW_PHY    , FR_WTR_CNP_SNW_PHY , FR_GLA_CNP_SNW_PHY , 
     O   TMP_GRSK_PHY       , FR_WTR_GRS_PHY     , FR_GLA_GRS_PHY     , 
     O   TMP_SNSK_PHY       , 
     O   INFO_SNOW_PHY      , 
     O   TMP_SOIL_PHY       , FR_WTR_SOIL_PHY    , FR_GLA_SOIL_PHY    , 
     O   TMP_SNSL_PHY       , FR_WTR_SNSL_PHY    , FR_GLA_SNSL_PHY    , 
     O   TMP_SNOW_PHY       , WTR_SNOW_PHY       , GLA_SNOW_PHY       , 
     O   RHO_SNOW_INV_PHY   , AGE_SNOW_PHY       , 
     O   ENG_SNOW_BUCKET_PHY, H2O_SNOW_BUCKET_PHY ) 
!
      USE PRM , ONLY :
     1   IJPHY , 
     1   JLPHY , 
     1   IDIM  , 
     1   JPHY2 ,
     1   ISPT   

      USE SIBPRM , ONLY :
     1   IDP   ,
     1   ISN 
!
      IMPLICIT NONE
!
! INPUT
!
      REAL(8),INTENT(OUT)  :: TMP_CNP_NOS_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: TMP_CNP_SNW_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: TMP_GRSK_PHY       (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: FR_WTR_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(OUT)  :: FR_GLA_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(OUT)  :: TMP_SNSK_PHY       (ISPT*IJPHY,JLPHY)
!
      REAL(8),INTENT(OUT)  :: TMP_SOIL_PHY       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT)  :: FR_WTR_SOIL_PHY    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT)  :: FR_GLA_SOIL_PHY    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT)  :: TMP_SNSL_PHY       (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: FR_WTR_SNSL_PHY    (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT)  :: FR_GLA_SNSL_PHY    (ISPT*IJPHY,JLPHY)    
!
      REAL(8),INTENT(OUT)  :: TMP_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT)  :: WTR_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT)  :: GLA_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT)  :: RHO_SNOW_INV_PHY   (ISPT*IJPHY,JLPHY,ISN)    
!
      REAL(8),INTENT(OUT)  :: AGE_SNOW_PHY       (ISPT*IJPHY,JLPHY)
      INTEGER,INTENT(OUT)  :: INFO_SNOW_PHY      (ISPT*IJPHY,JLPHY) 
!
      REAL(8),INTENT(OUT)  :: ENG_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT)  :: H2O_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY)

!
      INTEGER :: I 
!     INTEGER :: J 
      INTEGER :: IJ 
      INTEGER :: JL 
!
      DO JL=1,JLPHY
      DO IJ=1,IJPHY
!
        I = MOD(IJ-1,IDIM)+1
!       J = (IJ-1)/IDIM + 1 + (JL-1)*JPHY2  
!
        TMP_CNP_NOS_PHY   (IJ,JL) = 305.17477306397348
        FR_WTR_CNP_NOS_PHY(IJ,JL) = 7.84475679633962949D-113      ! íPà mm/m**2
!       FR_WTR_CNP_NOS_PHY(IJ,JL) = 1.314394372902081E-08 * 30000
        FR_WTR_CNP_NOS_PHY(IJ,JL) = 1.066000000000000e-04   *0.999999  !  vvv
!                                     ! Ç∆ÇËÇ†Ç¶Ç∏ÅAëOÇÃíl = ó 
        FR_GLA_CNP_NOS_PHY(IJ,JL) = 0.
!
        TMP_CNP_SNW_PHY   (IJ,JL) = 9.e30
        FR_WTR_CNP_SNW_PHY(IJ,JL) = 0.
        FR_GLA_CNP_SNW_PHY(IJ,JL) = 0.
!
        TMP_GRSK_PHY    (IJ,JL)   = 290.0460573597463 + 0.1          !ppp
!       TMP_GRSK_PHY    (IJ,JL)   = 308.25239056244521
        FR_WTR_GRS_PHY  (IJ,JL)   = 0.D0                         ! íPà  mm/m**2
        FR_WTR_GRS_PHY  (IJ,JL)   = 1.600000000000000D-04  *0.999999   !  vvv
!       FR_WTR_GRS_PHY  (IJ,JL)   = 1.190420383496277E-11 * 800
        FR_GLA_GRS_PHY  (IJ,JL)   = 0.
!
        TMP_SNSK_PHY    (IJ,JL)   = 1.D10          
!
!       TMP_SOIL_PHY    (IJ,JL,1) = 308.25239056244521
        TMP_SOIL_PHY    (IJ,JL,1) = 293.25239056244521
        TMP_SOIL_PHY    (IJ,JL,2) = 294.46443130614671
        TMP_SOIL_PHY    (IJ,JL,3) = 292.07075789044626
        FR_WTR_SOIL_PHY (IJ,JL,1) = 0.26793251158871367
        FR_WTR_SOIL_PHY (IJ,JL,2) = 0.32566380833675668
        FR_WTR_SOIL_PHY (IJ,JL,3) = 0.61199846674273839
        FR_GLA_SOIL_PHY (IJ,JL,1) = 0.
        FR_GLA_SOIL_PHY (IJ,JL,2) = 0.
        FR_GLA_SOIL_PHY (IJ,JL,3) = 0.
!
        TMP_SNSL_PHY    (IJ,JL)   = 0. 
        FR_WTR_SNSL_PHY (IJ,JL)   = 0. 
        FR_GLA_SNSL_PHY (IJ,JL)   = 0.
!
        AGE_SNOW_PHY    (IJ,JL)   = 10*24*3600+1 
        INFO_SNOW_PHY   (IJ,JL)   = -1 
!
        TMP_SNOW_PHY    (IJ,JL,1) = 200.
        TMP_SNOW_PHY    (IJ,JL,2) = 200.
        TMP_SNOW_PHY    (IJ,JL,3) = 200.
        TMP_SNOW_PHY    (IJ,JL,4) = 200.
        WTR_SNOW_PHY    (IJ,JL,1) = 0.
        WTR_SNOW_PHY    (IJ,JL,2) = 0.
        WTR_SNOW_PHY    (IJ,JL,3) = 0.
        WTR_SNOW_PHY    (IJ,JL,4) = 0.
        GLA_SNOW_PHY    (IJ,JL,1) = 0.
        GLA_SNOW_PHY    (IJ,JL,2) = 0.
        GLA_SNOW_PHY    (IJ,JL,3) = 0.
        GLA_SNOW_PHY    (IJ,JL,4) = 0.
        RHO_SNOW_INV_PHY(IJ,JL,1) = 1.
        RHO_SNOW_INV_PHY(IJ,JL,2) = 1.
        RHO_SNOW_INV_PHY(IJ,JL,3) = 1.
        RHO_SNOW_INV_PHY(IJ,JL,4) = 1.
!
        ENG_SNOW_BUCKET_PHY(IJ,JL) = 0.  
        H2O_SNOW_BUCKET_PHY(IJ,JL) = 0.  
!
! ê·
!
        IF ( I.EQ.2 .OR. I.EQ.3 ) THEN 
!
          TMP_CNP_NOS_PHY   (IJ,JL) = 263.32618441409886
          FR_WTR_CNP_NOS_PHY(IJ,JL) = 0.
          FR_GLA_CNP_NOS_PHY(IJ,JL) = 0.
!
          TMP_CNP_SNW_PHY   (IJ,JL) = 263.32618441409886
          FR_WTR_CNP_SNW_PHY(IJ,JL) = 0. 
          FR_GLA_CNP_SNW_PHY(IJ,JL) = 0.
!
          TMP_GRSK_PHY    (IJ,JL)   = 263.32618441409886 
          FR_WTR_GRS_PHY  (IJ,JL)   = 0.
          FR_GLA_GRS_PHY  (IJ,JL)   = 0.
!
          TMP_SNSK_PHY    (IJ,JL)   = 262.98401267482558
!
          TMP_SOIL_PHY    (IJ,JL,1) = 262.98401267482558
          TMP_SOIL_PHY    (IJ,JL,2) = 262.59382787514198
          TMP_SOIL_PHY    (IJ,JL,3) = 273.
!
          FR_WTR_SOIL_PHY (IJ,JL,1) = 0.
          FR_WTR_SOIL_PHY (IJ,JL,2) = 0.
          FR_WTR_SOIL_PHY (IJ,JL,3) = 0.
          FR_GLA_SOIL_PHY (IJ,JL,1) = 0.5
          FR_GLA_SOIL_PHY (IJ,JL,2) = 0.4
          FR_GLA_SOIL_PHY (IJ,JL,3) = 0.3
!
          TMP_SNSL_PHY    (IJ,JL)   = 262.98401267482558 - 0.35
          FR_WTR_SNSL_PHY (IJ,JL)   = 0. 
          FR_GLA_SNSL_PHY (IJ,JL)   = 0.5
!
!            INFO_SNOW_PHY   (IJ,JL)   = 7          ! 3 ëwãåê·
             INFO_SNOW_PHY   (IJ,JL)   = 5          ! 2 ëwãåê·
             AGE_SNOW_PHY    (IJ,JL)   = 10*24*3600+1 
             TMP_SNOW_PHY    (IJ,JL,1) = 262.98401267482558 -0.1
             TMP_SNOW_PHY    (IJ,JL,2) = 262.98401267482558 -0.2
             TMP_SNOW_PHY    (IJ,JL,3) = 0. 
             TMP_SNOW_PHY    (IJ,JL,4) = 0.
             WTR_SNOW_PHY    (IJ,JL,1) = 0.
             WTR_SNOW_PHY    (IJ,JL,2) = 0.
             WTR_SNOW_PHY    (IJ,JL,3) = 0.
             WTR_SNOW_PHY    (IJ,JL,4) = 0.
             GLA_SNOW_PHY    (IJ,JL,1) = 30.
             GLA_SNOW_PHY    (IJ,JL,2) = 70.
             GLA_SNOW_PHY    (IJ,JL,3) = 0. 
             GLA_SNOW_PHY    (IJ,JL,4) = 0.
             RHO_SNOW_INV_PHY(IJ,JL,1) = 1/200.
             RHO_SNOW_INV_PHY(IJ,JL,2) = 1/200.
             RHO_SNOW_INV_PHY(IJ,JL,3) = 1/200.
             RHO_SNOW_INV_PHY(IJ,JL,4) = 1.
!
             ENG_SNOW_BUCKET_PHY(IJ,JL) = 0.  
             H2O_SNOW_BUCKET_PHY(IJ,JL) = 0.  
        ENDIF
      ENDDO
      ENDDO
!
      RETURN
      END




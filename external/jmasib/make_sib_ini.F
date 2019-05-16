      SUBROUTINE MAKE_SIB_INI (
     I   IMASK              , ID_NEXT , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_SIB_ISLSCP 
!
      USE PRM    , ONLY : ISPT , IJPHY , JLPHY 
      USE SIBPRM , ONLY : IDP  , ISN 
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)  :: IMASK              (ISPT*IJPHY,JLPHY)
      INTEGER,INTENT(IN)  :: ID_NEXT(5) 
      REAL(8),INTENT(OUT) :: TMP_CNP_NOS_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT) :: FR_WTR_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: FR_GLA_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: TMP_CNP_SNW_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT) :: FR_WTR_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: FR_GLA_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: TMP_GRSK_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(OUT) :: FR_WTR_GRS_ALL     (ISPT*IJPHY,JLPHY)  
      REAL(8),INTENT(OUT) :: FR_GLA_GRS_ALL     (ISPT*IJPHY,JLPHY)  
      REAL(8),INTENT(OUT) :: TMP_SNSK_ALL       (ISPT*IJPHY,JLPHY)
!
!  土壌
!
      REAL(8),INTENT(OUT) :: TMP_SOIL_ALL       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT) :: FR_WTR_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT) :: FR_GLA_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(OUT) :: TMP_SNSL_ALL       (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: FR_WTR_SNSL_ALL    (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: FR_GLA_SNSL_ALL    (ISPT*IJPHY,JLPHY)    
!
!  積雪
!
      REAL(8),INTENT(OUT) :: TMP_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT) :: WTR_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT) :: GLA_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(OUT) :: RHO_SNOW_INV_ALL   (ISPT*IJPHY,JLPHY,ISN)
!
      REAL(8),INTENT(OUT) :: AGE_SNOW_ALL       (ISPT*IJPHY,JLPHY)    
      INTEGER,INTENT(OUT) :: INFO_SNOW_ALL      (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(OUT) :: ENG_SNOW_BUCKET_ALL(ISPT*IJPHY,JLPHY) 
      REAL(8),INTENT(OUT) :: H2O_SNOW_BUCKET_ALL(ISPT*IJPHY,JLPHY) 
!
!    GCM の場合
!
        IF      ( JCN_SIB_ISLSCP.EQ.-1 ) THEN
!
!         さて、どうやって作りますかねえ。
!
          stop 999
!
!    TEST1DIM の場合
!
        ELSE IF ( JCN_SIB_ISLSCP.EQ.0 ) THEN
!
          CALL TEST1D_SET_SIB_RESTART (
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
!   ISLSCP の場合
!
        ELSE IF ( JCN_SIB_ISLSCP.EQ.1 ) THEN
!
          CALL ISLSCP_SET_SIB_INI (
     I   IMASK , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
!   ERA の場合
!
        ELSE IF ( JCN_SIB_ISLSCP.EQ.2 ) THEN
!
          CALL ATM_SET_SIB_INI (
     I   IMASK , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
!   その他の強制データセットの場合
!
        ELSE IF (JCN_SIB_ISLSCP == 3) THEN
!
          CALL FORCING_INIT_SIB(
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
!
!   ISLSCP の場合
!
        ELSE IF ( JCN_SIB_ISLSCP.EQ.4 ) THEN
!
          CALL GSWP2_SET_SIB_INI (
     I   IMASK              , ID_NEXT            , 
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENG_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL ) 
      ELSE
        write(6,*) 'make_sib_ini : error' 
        stop 999
!
      ENDIF
!
      RETURN
      END


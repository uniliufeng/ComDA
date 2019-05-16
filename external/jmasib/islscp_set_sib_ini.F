      SUBROUTINE ISLSCP_SET_SIB_INI (
     I   IMASK              , 
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
     1   IDIM  , 
     1   JDIM  , 
     1   ISPT   
!
      USE ISLSCP_FILE_UNIT , ONLY :
     1    IFILE_TEMP   
!
      USE SIBPRM , ONLY :
     1   IDP   ,
     1   ISN 
!
      USE SIBCON , ONLY : TMP_FREZ_c 
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: IMASK              (IDIM,JDIM)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: TMP_CNP_NOS_PHY    (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_NOS_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_NOS_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: TMP_CNP_SNW_PHY    (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_CNP_SNW_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_CNP_SNW_PHY (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: TMP_GRSK_PHY       (IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: FR_WTR_GRS_PHY     (IDIM,JDIM)   
      REAL(8),INTENT(OUT)  :: FR_GLA_GRS_PHY     (IDIM,JDIM)   
      REAL(8),INTENT(OUT)  :: TMP_SNSK_PHY       (IDIM,JDIM)
!
      REAL(8),INTENT(OUT)  :: TMP_SOIL_PHY       (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: FR_WTR_SOIL_PHY    (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: FR_GLA_SOIL_PHY    (IDIM,JDIM,IDP)
      REAL(8),INTENT(OUT)  :: TMP_SNSL_PHY       (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_WTR_SNSL_PHY    (IDIM,JDIM)    
      REAL(8),INTENT(OUT)  :: FR_GLA_SNSL_PHY    (IDIM,JDIM)    
!
      REAL(8),INTENT(OUT)  :: TMP_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: WTR_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: GLA_SNOW_PHY       (IDIM,JDIM,ISN)
      REAL(8),INTENT(OUT)  :: RHO_SNOW_INV_PHY   (IDIM,JDIM,ISN)    
!
      REAL(8),INTENT(OUT)  :: AGE_SNOW_PHY       (IDIM,JDIM)
      INTEGER,INTENT(OUT)  :: INFO_SNOW_PHY      (IDIM,JDIM) 
!
      REAL(8),INTENT(OUT)  :: ENG_SNOW_BUCKET_PHY(IDIM,JDIM)
      REAL(8),INTENT(OUT)  :: H2O_SNOW_BUCKET_PHY(IDIM,JDIM)
!
      INTEGER :: I 
      INTEGER :: J 
      REAL(8) :: TMP 
      REAL(8) :: TMP_F 
      REAL(4) :: TMP_4(IDIM,JDIM) 
      CHARACTER(4) :: CMARK 
      INTEGER :: IDFILE(5)
!
! ======================
! >> —¤–Ê‰Šú’l‚ÌÝ’è <<  
! ======================
!
      CMARK = 'TEMP' 
      IDFILE(1) = 1987 
      IDFILE(2) = 1
      IDFILE(3) = 1
      IDFILE(4) = 0
      IDFILE(5) = 0
!
      CALL ISLSCP_RDFL ( IFILE_TEMP , CMARK , IDFILE , IMASK , TMP_4 ) 
!
      DO J=1,JDIM
      DO I=1,IDIM
        TMP   = TMP_4(I,J)
        TMP_F = TMP_FREZ_c + 1. 
        TMP_CNP_NOS_PHY (I,J)   = TMP
!
        IF ( TMP .GT. TMP_F ) THEN 
          TMP_GRSK_PHY    (I,J)   = TMP + 3.0
          TMP_SOIL_PHY    (I,J,1) = TMP + 2.0
          TMP_SOIL_PHY    (I,J,2) = TMP + 1.0
          TMP_SOIL_PHY    (I,J,3) = TMP
!
          FR_WTR_SOIL_PHY (I,J,1) = 0.5
          FR_WTR_SOIL_PHY (I,J,2) = 0.5
          FR_WTR_SOIL_PHY (I,J,3) = 0.5
          FR_GLA_SOIL_PHY (I,J,1) = 0.
          FR_GLA_SOIL_PHY (I,J,2) = 0.
          FR_GLA_SOIL_PHY (I,J,3) = 0.
        ELSE 
          TMP_GRSK_PHY    (I,J)   = TMP - 3.0
          TMP_SOIL_PHY    (I,J,1) = TMP - 2.0
          TMP_SOIL_PHY    (I,J,2) = TMP - 1.0
          TMP_SOIL_PHY    (I,J,3) = TMP
!
          FR_WTR_SOIL_PHY (I,J,1) = 0.
          FR_WTR_SOIL_PHY (I,J,2) = 0.
          FR_WTR_SOIL_PHY (I,J,3) = 0.
          FR_GLA_SOIL_PHY (I,J,1) = 0.5
          FR_GLA_SOIL_PHY (I,J,2) = 0.5
          FR_GLA_SOIL_PHY (I,J,3) = 0.5
        ENDIF
!
        FR_WTR_CNP_NOS_PHY(I,J) = 0.
        FR_GLA_CNP_NOS_PHY(I,J) = 0. 
        FR_WTR_GRS_PHY    (I,J) = 0.
        FR_GLA_GRS_PHY    (I,J) = 0.
!
        TMP_CNP_SNW_PHY   (I,J) = 0.
        FR_WTR_CNP_SNW_PHY(I,J) = 0.
        FR_GLA_CNP_SNW_PHY(I,J) = 0.
!
        TMP_SNSK_PHY      (I,J) = 0. 
        TMP_SNSL_PHY      (I,J) = 0. 
        FR_WTR_SNSL_PHY   (I,J) = 0. 
        FR_GLA_SNSL_PHY   (I,J) = 0.
!
        AGE_SNOW_PHY      (I,J)   = 10*24*3600 
        INFO_SNOW_PHY     (I,J)   = -1 
!
        TMP_SNOW_PHY    (I,J,1) = 0.
        TMP_SNOW_PHY    (I,J,2) = 0.
        TMP_SNOW_PHY    (I,J,3) = 0.
        TMP_SNOW_PHY    (I,J,4) = 0.
        WTR_SNOW_PHY    (I,J,1) = 0.
        WTR_SNOW_PHY    (I,J,2) = 0.
        WTR_SNOW_PHY    (I,J,3) = 0.
        WTR_SNOW_PHY    (I,J,4) = 0.
        GLA_SNOW_PHY    (I,J,1) = 0.
        GLA_SNOW_PHY    (I,J,2) = 0.
        GLA_SNOW_PHY    (I,J,3) = 0.
        GLA_SNOW_PHY    (I,J,4) = 0.
        RHO_SNOW_INV_PHY(I,J,1) = 0.
        RHO_SNOW_INV_PHY(I,J,2) = 0.
        RHO_SNOW_INV_PHY(I,J,3) = 0.
        RHO_SNOW_INV_PHY(I,J,4) = 0.
!
        ENG_SNOW_BUCKET_PHY(I,J) = 0.  
        H2O_SNOW_BUCKET_PHY(I,J) = 0.  
!
      ENDDO
      ENDDO
!
      RETURN
      END

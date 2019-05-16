      SUBROUTINE MOVE_ICE_INTO_SNOW (
     I      NSTR   , NEND    , 
     M      H2O_BUCKET_INTO_SOIL       ,
     M      ENH_BUCKET_INTO_SOIL       , 
     M      H2O_SNOW_BUCKET  , ENH_SNOW_BUCKET , 
     M      INFO_SNOW  )
!
      USE SIBPRM , ONLY : IDX 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, INN ,
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      INTEGER,INTENT(IN)    ::  NSTR   
      INTEGER,INTENT(IN)    ::  NEND
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO_SOIL (IDX)
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_INTO_SOIL (IDX)
      REAL(8),INTENT(INOUT) :: H2O_SNOW_BUCKET      (IDX)
      REAL(8),INTENT(INOUT) :: ENH_SNOW_BUCKET      (IDX)
      INTEGER,INTENT(INOUT) :: INFO_SNOW            (IDX)
!
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      CHARACTER(7) :: CMARK 
      real(8)      :: work_1(idx) 
      real(8)      :: work_2(idx) 
#endif
!
      INTEGER  :: I
!       
#ifdef MONYOS
      WORK_1(:) = 0.
      WORK_2(:) = 0.
#endif
!
      DO I=NSTR,NEND
        IF ( ENH_BUCKET_INTO_SOIL(I).LT.0 ) THEN
          h2o_snow_bucket(i) = h2o_snow_bucket(i) 
     1                        + H2O_BUCKET_INTO_SOIL(I) 
          enh_snow_bucket(i) = enh_snow_bucket(i) 
     1                        + ENH_BUCKET_INTO_SOIL(I) 
#ifdef MONYOS
          WORK_1(I) = H2O_BUCKET_INTO_SOIL(I) 
          WORK_2(I) = ENH_BUCKET_INTO_SOIL(I) 
#endif
          H2O_BUCKET_INTO_SOIL(I) = 0 
          ENH_BUCKET_INTO_SOIL(I) = 0 
          IF ( INFO_SNOW(I) .EQ. -1 ) THEN
            INFO_SNOW(I) = -2 
          ENDIF
        ENDIF
      ENDDO
!
#ifdef MONYOS
      IF ( IMONIT_LEVEl .GE. IMONIT_ALL ) THEN
         CMARK = 'WMVN2LS'
         CALL MONIT_SIB ( CMARK, WORK_1, INN , ONE , ONE )       
         CMARK = 'IMVN2LS'
         CALL MONIT_SIB ( CMARK, WORK_2, INN , ONE , ONE )       
      ENDIF
#endif
!
      RETURN
      END
!================================================================
      SUBROUTINE MR_MOVE_ICE_INTO_SNOW 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
#ifdef MONYOS
      IF ( IMONIT_LEVEl .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'WMVN2LS' 
      CTITLE = 'ICE MOVE BEFORE INTOSOIL TO SNOW'  
      CUNIT  = 'KG'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IMVN2LS' 
      CTITLE = 'ICE MOVE BEFORE INTOSOIL TO SNOW'  
      CUNIT  = 'KG'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      END 


      MODULE RESTART_SIB0109
!
      CONTAINS 
!
      SUBROUTINE WRITE_RESTART_SIB0109 ( 
     I   IUNIT_RESTART_OUT  , 
     I   IDNEXT             , 
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENH_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL,
     I   TMP_SOIL4_ALL , 
     I   ALB_SNOW_SNW_ALL   , RAIN_1HOUR_ALL )
!
      USE PRM , ONLY :
     1   ISPT  , 
     1   IJPHY , 
     1   JLPHY 
!
      USE SIBPRM , ONLY :
     1   IDP  , 
     1   ISN  , IVN
!
      IMPLICIT NONE 
!
      INTEGER      ,intent(in) :: IUNIT_RESTART_OUT  
      INTEGER,intent(in) :: IDNEXT (5)
!
      REAL(8),intent(in) :: TMP_CNP_NOS_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),intent(in) :: FR_WTR_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(in) :: FR_GLA_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(in) :: TMP_CNP_SNW_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),intent(in) :: FR_WTR_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(in) :: FR_GLA_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(in) :: TMP_GRSK_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),intent(in) :: FR_WTR_GRS_ALL     (ISPT*IJPHY,JLPHY)   
      REAL(8),intent(in) :: FR_GLA_GRS_ALL     (ISPT*IJPHY,JLPHY)   
      REAL(8),intent(in) :: TMP_SNSK_ALL       (ISPT*IJPHY,JLPHY)
!
      REAL(8),intent(in) :: TMP_SOIL_ALL       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(in) :: FR_WTR_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(in) :: FR_GLA_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(in) :: TMP_SNSL_ALL       (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),intent(in) :: FR_WTR_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),intent(in) :: FR_GLA_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
!
      REAL(8),intent(in) :: TMP_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(in) :: WTR_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(in) :: GLA_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(in) :: RHO_SNOW_INV_ALL   (ISPT*IJPHY,JLPHY,ISN)    
!
      REAL(8),intent(in) :: TMP_SOIL4_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),intent(in) :: AGE_SNOW_ALL       (ISPT*IJPHY,JLPHY)
      INTEGER,intent(in) :: INFO_SNOW_ALL      (ISPT*IJPHY,JLPHY) 
!
      REAL(8),intent(in) :: ENH_SNOW_BUCKET_ALL   (ISPT*IJPHY,JLPHY) 
      REAL(8),intent(in) :: H2O_SNOW_BUCKET_ALL   (ISPT*IJPHY,JLPHY) 
!
      REAL(8),intent(in) :: ALB_SNOW_SNW_ALL   (ISPT*IJPHY,JLPHY,IVN)
      REAL(8),intent(in) :: RAIN_1HOUR_ALL     (ISPT*IJPHY,JLPHY)
!
      WRITE(IUNIT_RESTART_OUT)  'NS_RS_04' , IDNEXT 
      WRITE(IUNIT_RESTART_OUT)  
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENH_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL 
      WRITE(IUNIT_RESTART_OUT)  
     I   TMP_SOIL4_ALL
      WRITE(IUNIT_RESTART_OUT)  
     I   ALB_SNOW_SNW_ALL  
      WRITE(IUNIT_RESTART_OUT)  
     I   RAIN_1HOUR_ALL
!
      WRITE(6,*) 'WRITE_RESTART END'
!
      RETURN
      END SUBROUTINE WRITE_RESTART_SIB0109  
!=================================================================== 
      SUBROUTINE READ_RESTART_SIB0109 ( 
     I   IUNIT_RESTART_IN   , 
     I   IDNOW              ,
     O   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     O   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     O   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     O   TMP_SNSK_ALL       , 
     O   INFO_SNOW_ALL      , 
     O   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     O   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     O   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     O   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     O   ENH_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL,
     O   TMP_SOIL4_ALL      , 
     O   ALB_SNOW_SNW_ALL   , RAIN_1HOUR_ALL ) 
!
      USE PRM , ONLY :
     1   ISPT  , 
     1   IJPHY , 
     1   JLPHY 
!
      USE SIBPRM , ONLY :
     1   IDP  , 
     1   ISN  , IVN
!
      use com_jobinfo_sib0109, only : nocheck_restartfile_date
!
      IMPLICIT NONE 
!
      INTEGER      ,intent(in) :: IUNIT_RESTART_IN

      INTEGER,intent(in)  :: IDNOW(5)
      INTEGER             :: ID_READ(5)
!
      REAL(8),intent(out) :: TMP_CNP_NOS_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),intent(out) :: FR_WTR_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(out) :: FR_GLA_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(out) :: TMP_CNP_SNW_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),intent(out) :: FR_WTR_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(out) :: FR_GLA_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),intent(out) :: TMP_GRSK_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),intent(out) :: FR_WTR_GRS_ALL     (ISPT*IJPHY,JLPHY)   
      REAL(8),intent(out) :: FR_GLA_GRS_ALL     (ISPT*IJPHY,JLPHY)   
      REAL(8),intent(out) :: TMP_SNSK_ALL       (ISPT*IJPHY,JLPHY)
!
      REAL(8),intent(out) :: TMP_SOIL_ALL       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(out) :: FR_WTR_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(out) :: FR_GLA_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),intent(out) :: TMP_SNSL_ALL       (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),intent(out) :: FR_WTR_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),intent(out) :: FR_GLA_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
!
      REAL(8),intent(out) :: TMP_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(out) :: WTR_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(out) :: GLA_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),intent(out) :: RHO_SNOW_INV_ALL   (ISPT*IJPHY,JLPHY,ISN)    
!
      REAL(8),intent(out) :: AGE_SNOW_ALL       (ISPT*IJPHY,JLPHY)
      INTEGER,intent(out) :: INFO_SNOW_ALL      (ISPT*IJPHY,JLPHY) 
!
      REAL(8),intent(out) :: ENH_SNOW_BUCKET_ALL   (ISPT*IJPHY,JLPHY) 
      REAL(8),intent(out) :: H2O_SNOW_BUCKET_ALL   (ISPT*IJPHY,JLPHY) 
!
      REAL(8),intent(out) :: TMP_SOIL4_ALL   (ISPT*IJPHY,JLPHY) 
      REAL(8),intent(out) :: ALB_SNOW_SNW_ALL   (ISPT*IJPHY,JLPHY,IVN)
      REAL(8),intent(out) :: RAIN_1HOUR_ALL     (ISPT*IJPHY,JLPHY)
!
      CHARACTER(8)        :: C8 
      CHARACTER(8),SAVE   :: C8_MARK_01 = 'NS_RS_01' 
      CHARACTER(8),SAVE   :: C8_MARK_02 = 'NS_RS_02' 
      CHARACTER(8),SAVE   :: C8_MARK_03 = 'NS_RS_03' 
      CHARACTER(8),SAVE   :: C8_MARK_04 = 'NS_RS_04' 
!
      READ(IUNIT_RESTART_IN,ERR=999)  C8 , ID_READ 
!
      IF (       C8(1:7) .NE. C8_MARK_01(1:7) ) THEN 
        write(6,*) 'read_restart_sib0109 error ' , C8 
        stop 999
      ENDIF
!
      IF (  NOCHECK_RESTARTFILE_DATE == 0 .and. 
     1      (     IDNOW(1).NE.ID_READ(1)
     1       .OR. IDNOW(2).NE.ID_READ(2)
     1       .OR. IDNOW(3).NE.ID_READ(3)
     1       .OR. IDNOW(4).NE.ID_READ(4)
     1       .OR. IDNOW(5).NE.ID_READ(5) ) 
     1    )  THEN
        write(6,*) 'read_restart_sib0109 error ' , 
     1             IDNOW , ' ' , ID_READ   
        stop 999
      ENDIF

!
      READ(IUNIT_RESTART_IN,ERR=999)  
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       , 
     I   ENH_SNOW_BUCKET_ALL, H2O_SNOW_BUCKET_ALL 
!
      WRITE(6,*) 'READ_RESTART_SIB0109 END ID=' , IDNOW 
!
      IF     ( C8 == C8_MARK_01 ) THEN
        TMP_SOIL4_ALL(:,:) 
     1        = MAX(TMP_SOIL_ALL(:,:,3),TMP_SNSL_ALL(:,:,3))
      ELSEIF (      C8 == C8_MARK_02 
     1         .OR. C8 == C8_MARK_03 
     1         .OR. C8 == C8_MARK_04 ) THEN  
        READ(IUNIT_RESTART_IN,ERR=999)  
     I   TMP_SOIL4_ALL    
      ELSE
         write(6,*) 'read_restart_sib0109 error ' , C8
         stop 999 
      ENDIF
!
! ALB_SNOW_SNW_ALL 
!
      IF     ( C8 == C8_MARK_01 .OR. C8 == C8_MARK_02) THEN
        ALB_SNOW_SNW_ALL(:,:,1) = 0.8
        ALB_SNOW_SNW_ALL(:,:,2) = 0.7
      ELSEIF ( C8 == C8_MARK_03 .OR. C8 == C8_MARK_04 ) THEN 
        READ(IUNIT_RESTART_IN,ERR=999)  
     I            ALB_SNOW_SNW_ALL 
      ELSE
         write(6,*) 'read_restart_sib0109 error ' , C8
         stop 999 
      ENDIF
!
! RAIN_1HOUR_ALL
!
      IF     (      C8 == C8_MARK_01 .OR. C8 == C8_MARK_02 
     1         .OR. C8 == C8_MARK_03)  THEN
        RAIN_1HOUR_ALL(:,:) = 0.
      ELSEIF ( C8 == C8_MARK_04 )  THEN
        READ(IUNIT_RESTART_IN,ERR=999)  
     I            RAIN_1HOUR_ALL
      ELSE
         write(6,*) 'read_restart_sib0109 error ' , C8
         stop 999 
      ENDIF
!
      RETURN
!
 999  CONTINUE
!
      WRITE(6,*) 'READ_RESTART_SIB0109 ERROR' 
      STOP   
!
      END SUBROUTINE READ_RESTART_SIB0109  
!=================================================================== 
      SUBROUTINE CLEAN_RESTART_SIB0109 ( 
     I   TMP_CNP_NOS_ALL    , FR_WTR_CNP_NOS_ALL , FR_GLA_CNP_NOS_ALL , 
     I   TMP_CNP_SNW_ALL    , FR_WTR_CNP_SNW_ALL , FR_GLA_CNP_SNW_ALL , 
     I   TMP_GRSK_ALL       , FR_WTR_GRS_ALL     , FR_GLA_GRS_ALL     , 
     I   TMP_SNSK_ALL       , 
     I   INFO_SNOW_ALL      , 
     I   TMP_SOIL_ALL       , FR_WTR_SOIL_ALL    , FR_GLA_SOIL_ALL    , 
     I   TMP_SNSL_ALL       , FR_WTR_SNSL_ALL    , FR_GLA_SNSL_ALL    , 
     I   TMP_SNOW_ALL       , WTR_SNOW_ALL       , GLA_SNOW_ALL       , 
     I   RHO_SNOW_INV_ALL   , AGE_SNOW_ALL       ) 
!
      USE PRM , ONLY :
     1   ISPT  , 
     1   IJPHY , 
     1   JLPHY 
!
      USE SIBPRM , ONLY :
     1   IDP  , 
     1   ISN  
!
      USE SIBCON , ONLY :
     1   AGE_SNOW_MAX2_c 
!
      IMPLICIT NONE 
!
      REAL(8),INTENT(INOUT):: TMP_CNP_NOS_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT):: FR_WTR_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT):: FR_GLA_CNP_NOS_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT):: TMP_CNP_SNW_ALL    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT):: FR_WTR_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT):: FR_GLA_CNP_SNW_ALL (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(INOUT):: TMP_GRSK_ALL       (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT):: FR_WTR_GRS_ALL     (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT):: FR_GLA_GRS_ALL     (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(INOUT):: TMP_SNSK_ALL       (ISPT*IJPHY,JLPHY)
!
      REAL(8),INTENT(INOUT):: TMP_SOIL_ALL       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT):: FR_WTR_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(INOUT):: FR_GLA_SOIL_ALL    (ISPT*IJPHY,JLPHY,IDP)     
      REAL(8),INTENT(INOUT):: TMP_SNSL_ALL       (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),INTENT(INOUT):: FR_WTR_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),INTENT(INOUT):: FR_GLA_SNSL_ALL    (ISPT*IJPHY,JLPHY,IDP)    
!
      REAL(8),INTENT(INOUT):: TMP_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(INOUT):: WTR_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(INOUT):: GLA_SNOW_ALL       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(INOUT):: RHO_SNOW_INV_ALL   (ISPT*IJPHY,JLPHY,ISN)    
!
      REAL(8),INTENT(INOUT):: AGE_SNOW_ALL       (ISPT*IJPHY,JLPHY)
      INTEGER,INTENT(IN)   :: INFO_SNOW_ALL      (ISPT*IJPHY,JLPHY)
!
      LOGICAL, SAVE :: LFIRST
!
      INTEGER :: I 
      INTEGER :: J 
!
#ifdef SIB_DEBUG
      integer :: imsg_num
      integer :: ibad_ij
      integer :: ibad_jl
#endif

!
      DATA LFIRST /.TRUE./ 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CLEAN_RESTART_NEWSIB POCHAKA 99/11/14'
        LFIRST = .FALSE.
      ENDIF
!
!  ê·ñ≥ÇµàÊÇÃë|èú
!     
#ifdef SIB_DEBUG
       imsg_num = 0 
#endif
!
      DO J=1,JLPHY
      DO I=1,ISPT*IJPHY 
        IF ( INFO_SNOW_ALL(I,J) .GE. 1 ) THEN   ! ê·Ç™ëSñ àÍëwà»è„
          TMP_CNP_NOS_ALL    (I,J)    = 0. 
          FR_WTR_CNP_NOS_ALL (I,J)    = 0.  
          FR_GLA_CNP_NOS_ALL (I,J)    = 0. 
          TMP_GRSK_ALL       (I,J)    = 0. 
          FR_WTR_GRS_ALL     (I,J)    = 0. 
          FR_GLA_GRS_ALL     (I,J)    = 0. 
          TMP_SOIL_ALL       (I,J,:)  = 0. 
          FR_WTR_SOIL_ALL    (I,J,:)  = 0. 
          FR_GLA_SOIL_ALL    (I,J,:)  = 0. 
        ELSEIF ( INFO_SNOW_ALL(I,J) .LE. -1 ) THEN   ! ê·ñ≥Çµ -1 ( -2 è≠ê· ) 
          TMP_CNP_SNW_ALL    (I,J)    = 0. 
          FR_WTR_CNP_SNW_ALL (I,J)    = 0.  
          FR_GLA_CNP_SNW_ALL (I,J)    = 0. 
          TMP_SNSK_ALL       (I,J)    = 0. 
          TMP_SNSL_ALL       (I,J,:)  = 0. 
          FR_WTR_SNSL_ALL    (I,J,:)  = 0. 
          FR_GLA_SNSL_ALL    (I,J,:)  = 0. 
!         TMP_SNOW_ALL       (I,J,:)  = 0. 
          TMP_SNOW_ALL       (I,J,:)  = 0. 
          WTR_SNOW_ALL       (I,J,:)  = 0. 
          GLA_SNOW_ALL       (I,J,:)  = 0. 
          RHO_SNOW_INV_ALL   (I,J,:)  = 0. 
          AGE_SNOW_ALL       (I,J)    = AGE_SNOW_MAX2_c 
        ELSEIF ( INFO_SNOW_ALL(I,J) .EQ. 0 ) THEN    ! ã§ë∂
!
        ELSE
#ifdef SIB_DEBUG
          imsg_num = 1 
          ibad_ij  = i
          ibad_jl  = j
#endif
        ENDIF
      ENDDO
      ENDDO
!
#ifdef SIB_DEBUG
      if ( imsg_num.gt.0 ) then
          WRITE(6,*) 'CLEAN_RESTART ERROR INFO_SNOW ' , 
     1                ' TOTAL NUMBER = ' , IMSG_NUM ,     
     1                ' FOR EXAMPLE IJ JL ' , IBAD_IJ , IBAD_JL 
          STOP 999
      endif
#endif
!
      RETURN
      END SUBROUTINE CLEAN_RESTART_SIB0109 
! 
      END MODULE RESTART_SIB0109

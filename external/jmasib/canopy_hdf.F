!  
      SUBROUTINE CANOPY_HDF ( 
     I    NSTR  , NEND , 
     I    ITR     , ITR_MAX , 
#ifdef MONYOS
     I    INSFLAG , DELT_CNP , 
#endif
     I    TMP_SKIN , TMP_SS1 , TMP_SS2  , 
     I    HDF_SKIN , 
     O    HDF  )
!
      USE SIBPRM , ONLY :
     1    IDX 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN , ICS ,
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
#endif
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: ITR
      INTEGER,INTENT(IN)   :: ITR_MAX
!
      REAL(8),INTENT(IN)   :: TMP_SS1 (IDX) ! 第一層の温度（時間変化）
      REAL(8),INTENT(IN)   :: TMP_SS2 (IDX) ! 第二層の温度（しばし固定）
      REAL(8),INTENT(IN)   :: TMP_SKIN(IDX) ! スキン温度
      REAL(8),INTENT(IN)   :: HDF_SKIN(IDX,2) ! スキン <-> 1層 熱伝導係数
!
#ifdef MONYOS
      INTEGER,INTENT(IN)   :: INSFLAG  
      REAL(8),INTENT(IN)   :: DELT_CNP
#endif
!
      REAL(8),INTENT(OUT)  :: HDF     (IDX,2)
!
! WORK
!
      INTEGER :: I  
!
#ifdef MONYOS
      REAL(8),PARAMETER    :: ONE = 1. 
      CHARACTER(7)         :: CMARK 
#endif
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_HDF 00/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
! 
      DO I=NSTR,NEND
         HDF (I,1) = HDF_SKIN(I,1)  * ( TMP_SKIN(I) - TMP_SS1(I) )
#ifdef KILLHDF
         HDF (I,1) = 0                                             
#endif
!
         HDF (I,2) = HDF_SKIN(I,2)  * ( TMP_SS1(I) - TMP_SS2(I) )
!
#ifdef SIB_DEBUG
         if ( abs(hdf(i,1)).gt. 200. ) then
           write(6,*) 'canopy_hdf too large sk_s1 ' , i , 
     1                ' flux ' , hdf(i,1) , 
     1                ' coef ' , HDF_SKIN(I,1)  , 
     1                ' tmp  ' , TMP_SKIN(I) , TMP_SS1(I)
         endif
         if ( abs(hdf(i,2)).gt. 100. ) then
           write(6,*) 'canopy_hdf too large s1_s2' , i , hdf(i,2) , 
     1                HDF_SKIN(I,2)  ,  TMP_SS1(I) , TMP_SS2(I)
         endif
#endif

      ENDDO
!
! ================
! >> 詳細モニタ <<
! ================
! 
#ifdef MONYOS
      IF ( ITR.EQ.ITR_MAX) THEN
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( INSFLAG.EQ.0) THEN
          CMARK  = 'HDFBNL0'
          CALL MONIT_SIB ( CMARK , HDF(1,1) ,ICN, DELT_CNP , ONE )  
        ELSE
          CMARK  = 'HDFBSW0'
          CALL MONIT_SIB ( CMARK , HDF(1,1) ,ICS, DELT_CNP , ONE )  
        ENDIF
      ENDIF
      ENDIF
#endif
!
      RETURN
      END
!=====================================================================
      SUBROUTINE MR_CANOPY_HDF
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_CANOPY_HDF 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! HDF 
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'HDFBNL0'
      CTITLE = 'HEAT DIFFUSE GSKIN TO D1 NOS explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HDFBSW0'
      CTITLE = 'HEAT DIFFUSE GSKIN TO D1 SNW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END




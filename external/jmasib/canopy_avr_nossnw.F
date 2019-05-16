      SUBROUTINE CANOPY_AVR_NOSSNW  ( 
     I    NMAX       , NMAX_0     , NMAX_01    ,  NMAX_1   , 
#ifdef MONYOS
     I    DELT_CNP   , 
#endif
     I    CVR_NOS    , CVR_SNW    , 
     I    GMUV_NOS   , GMUV_SNW   , 
     I    GMT_NOS    , GMT_SNW    , 
     I    GMQ_NOS    , GMQ_SNW    , 
     I    RAD_LONG_NOS_2_ATM      , RAD_LONG_SNW_2_ATM     , 
     O    GMT_OUT    , GMQ_OUT    , GMUV_OUT   , 
     O    RAD_LONG_SIB_2_ATM )
!
      USE SIBPRM , ONLY :
     1    IDX 
!
#ifdef MONYOS
      USE SIBCON , ONLY : STEFAN_c 
!
      USE SIB_MONIT , ONLY : MONIT_SIB, INA, 
     1     IMONIT_LEVEL    , 
     1     IMONIT_MEAN_MRI , 
     1     IMONIT_ALL      , IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX 
      INTEGER,INTENT(IN)    :: NMAX_0 (3) 
      INTEGER,INTENT(IN)    :: NMAX_01(3) 
      INTEGER,INTENT(IN)    :: NMAX_1 (3) 
!
#ifdef MONYOS
      REAL(8),INTENT(IN)    :: DELT_CNP 
#endif
      REAL(8),INTENT(IN)    :: CVR_NOS           (IDX)
      REAL(8),INTENT(IN)    :: CVR_SNW           (IDX)
!     REAL(8),INTENT(IN)    :: U_FLUX_NOS        (IDX)
!     REAL(8),INTENT(IN)    :: V_FLUX_NOS        (IDX)
!     REAL(8),INTENT(IN)    :: U_FLUX_SNW        (IDX)
!     REAL(8),INTENT(IN)    :: V_FLUX_SNW        (IDX)
      REAL(8),INTENT(IN)    :: GMUV_NOS          (IDX,2)
      REAL(8),INTENT(IN)    :: GMUV_SNW          (IDX,2)
      REAL(8),INTENT(IN)    :: GMT_NOS           (IDX)
      REAL(8),INTENT(IN)    :: GMT_SNW           (IDX)
      REAL(8),INTENT(IN)    :: GMQ_NOS           (IDX)
      REAL(8),INTENT(IN)    :: GMQ_SNW           (IDX)
!
      REAL(8),INTENT(IN)   :: RAD_LONG_NOS_2_ATM(IDX) 
      REAL(8),INTENT(IN)   :: RAD_LONG_SNW_2_ATM(IDX) 
!
! OUTPUT
!
!     REAL(8),INTENT(OUT)   :: U_FLUX   (IDX)   ! 運動量フラックス
!     REAL(8),INTENT(OUT)   :: V_FLUX   (IDX)   ! 同上
!
!                       M.Hosaka bug fix from here 2003/01/14
!!      REAL(8),INTENT(OUT)   :: GMUV_OUT (IDX,3)
!!      REAL(8),INTENT(OUT)   :: GMT_OUT  (IDX,2)
!!      REAL(8),INTENT(OUT)   :: GMQ_OUT  (IDX,2)
      REAL(8),INTENT(INOUT)   :: GMUV_OUT (IDX,3)
      REAL(8),INTENT(INOUT)   :: GMT_OUT  (IDX,2)
      REAL(8),INTENT(INOUT)   :: GMQ_OUT  (IDX,2)
!                       M.Hosaka bug fix to   here 2003/01/14
      REAL(8),INTENT(OUT)   :: RAD_LONG_SIB_2_ATM(IDX) 
!
! WORK
!
      INTEGER :: I
!
#ifdef MONYOS
      CHARACTER( 7)  :: CMARK 
      REAL(8)        :: WORK(IDX)
      REAL(8),PARAMETER :: ONE = 1.D0 
#endif
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST  / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_AVR_NOSSNW POCAHKA 99/08/12'  
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=1,NMAX_0(3) 
        GMUV_OUT (I,2)  = GMUV_NOS  (I,1) 
        GMUV_OUT (I,3)  = GMUV_NOS  (I,2) 
        GMT_OUT  (I,2)  = GMT_NOS   (I) 
        GMQ_OUT  (I,2)  = GMQ_NOS   (I) 
!       U_FLUX   (I)    = U_FLUX_NOS(I)
!       V_FLUX   (I)    = V_FLUX_NOS(I)
        RAD_LONG_SIB_2_ATM(I) = RAD_LONG_NOS_2_ATM(I)  
 1000 CONTINUE
!
      DO 2000 I=NMAX_01(2),NMAX_01(3) 
        GMUV_OUT (I,2)  =   GMUV_NOS  (I,1) * CVR_NOS(I) 
     1                    + GMUV_SNW  (I,1) * CVR_SNW(I) 
        GMUV_OUT (I,3)  =   GMUV_NOS  (I,2) * CVR_NOS(I) 
     1                    + GMUV_SNW  (I,2) * CVR_SNW(I) 
        GMT_OUT  (I,2)  =   GMT_NOS   (I)   * CVR_NOS(I) 
     1                    + GMT_SNW   (I)   * CVR_SNW(I) 
        GMQ_OUT  (I,2)   =  GMQ_NOS   (I)   * CVR_NOS(I) 
     1                    + GMQ_SNW   (I)   * CVR_SNW(I) 
        RAD_LONG_SIB_2_ATM(I) =   RAD_LONG_NOS_2_ATM(I) * CVR_NOS(I) 
     1                          + RAD_LONG_SNW_2_ATM(I) * CVR_SNW(I) 
 2000 CONTINUE
!
      DO 3000 I=NMAX_1(2),NMAX 
        GMUV_OUT (I,2)   = GMUV_SNW  (I,1) 
        GMUV_OUT (I,3)   = GMUV_SNW  (I,2) 
        GMT_OUT  (I,2)   = GMT_SNW   (I) 
        GMQ_OUT  (I,2)   = GMQ_SNW   (I) 
        RAD_LONG_SIB_2_ATM(I) = RAD_LONG_SNW_2_ATM(I)  
 3000 CONTINUE
!
! 有効放射温度のモニタ ( 登録は monit_regist_sib_all にて ) 
!
#ifdef MONYOS
      IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALMA
     1     ) THEN
        DO I=1,NMAX
#ifdef SIB_DEBUG
          IF (       RAD_LONG_SIB_2_ATM(I) .GT. - STEFAN_c *150.D0 **4 
     1         .OR.  RAD_LONG_SIB_2_ATM(I) .LT. - STEFAN_c *400.D0 **4 
     1        ) THEN
            write(6,*) 'canopy_avr_nossnw error' , i , 
     1                  RAD_LONG_SIB_2_ATM(I) ,
     1                  STEFAN_C * 100.D0 **4
            stop 999
          ENDIF
#endif
          WORK(I) = SQRT(SQRT( (-RAD_LONG_SIB_2_ATM(I))/ STEFAN_c ))
!         WORK(I) = - RAD_LONG_SIB_2_ATM(I)
        ENDDO
!
        IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
          CMARK='TGEF'
          CALL MONIT_SIB ( CMARK , WORK , INA , DELT_CNP , ONE )      
        ENDIF 
        IF (      IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          CMARK='A_TRAD'
          CALL MONIT_SIB ( CMARK , WORK , INA , DELT_CNP , ONE )      
          CMARK='A_LWUP'
          CALL MONIT_SIB ( CMARK , RAD_LONG_SIB_2_ATM , INA , 
     1                     DELT_CNP , - ONE )      
        ENDIF 
      ENDIF
#endif 
!
      RETURN
      END
!============================================================
      SUBROUTINE MR_CANOPY_AVR_NOSSNW
!
      USE SIB_MONIT , ONLY : IMONIT_LEVEL , IMONIT_MEAN_MRI, IMONIT_ALMA
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_AVR_NOSSNW 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  有効放射温度 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
!
        CMARK  = 'TGEF'           
        CTITLE = 'TEMP EFFECTIVE'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
!  ALMA 有効放射温度 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_TRAD'           
        CTITLE = 'ALMA RadT'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'A_LWUP'           
        CTITLE = 'ALMA - LWnet Longwave upward'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      RETURN
      END 

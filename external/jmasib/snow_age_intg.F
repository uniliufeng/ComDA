!
      SUBROUTINE SNOW_AGE_INTG (
     I   NMAX   , NMAX_1 , DELT_SNOW , 
     I   RHO_SNOW_INV , 
     M   AGE_SNOW_PROG ) 
!
      USE SIBPRM , ONLY :
     1   IDX   , 
     1   ISN  
!
      USE SIBCON , ONLY : 
     1    AGE_SNOW_MAX2_c    ,             ! /* 新雪層' の年齢上限 */
     1    RHO_SNOW_AGECLEAR_THR_INV_c ,
     1    ITUNE_SNOW_AGE              , itune_albedo_snow  
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, INA, ICS ,
     1                      IMONIT_LEVEL , IMONIT_ALL   
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR 
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_1(3) 
      REAL(8),INTENT(IN)   :: DELT_SNOW  
!
      REAL(8),INTENT(IN)   :: RHO_SNOW_INV       (IDX,ISN)
!
      REAL(8),INTENT(INOUT):: AGE_SNOW_PROG      (IDX)         ! 雪の年齢
!
      INTEGER      :: I
      LOGICAL,SAVE :: LFIRST 
!
#ifdef MONYOS
      REAL(8),PARAMETER  :: ONE =  1.D0 
      CHARACTER(7)       :: CMARK
#endif
!
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_AGE_INTG POCHAKA 01/08/17' 
        LFIRST = .FALSE.
      ENDIF
!
      IF ( ITUNE_ALBEDO_SNOW .LE. 2 ) THEN
!
!  ==========
!  >> 加齢 <<
!  ==========
!
      DO I=1,NMAX
        AGE_SNOW_PROG(I)  = MIN ( AGE_SNOW_PROG(I) + DELT_SNOW , 
     1                            AGE_SNOW_MAX2_c  )      
      ENDDO
!
! ========================
! >> 雪年齢クリア方法 1 <<
! ========================
!
! 1 層以上あって、第一層の密度が RHO_SNOW_AGECLEAR_THR_c より
! 小さければクリアする。
!
      IF ( ITUNE_SNOW_AGE .EQ. 1 ) THEN 
        DO I=NMAX_1(2)+1,NMAX
          IF ( RHO_SNOW_INV(I,1) 
     1         .GT. RHO_SNOW_AGECLEAR_THR_INV_c ) THEN
            AGE_SNOW_PROG(I) = 0. 
          ENDIF
        ENDDO
      ENDIF 
!
      ENDIF
!
! ==========
! > モニタ <
! ==========
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
! 雪の年齢
        CMARK = 'AGESNWM'
        CALL MONIT_SIB ( CMARK , AGE_SNOW_PROG , INA, DELT_SNOW, ONE)
      ENDIF
#endif
!
!
      RETURN
      END
!============================================================
      SUBROUTINE MR_SNOW_AGE_INTG 
!
      USE SIBPRM , ONLY :
     1   ISN 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL   
#endif
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
!
      INTEGER :: K 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_AGE_INTG 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'AGESNWM'
      CTITLE = 'SNOW AGE'
      CUNIT  = 's'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

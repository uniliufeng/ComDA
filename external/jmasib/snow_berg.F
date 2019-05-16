       SUBROUTINE SNOW_BERG ( 
     I    NMAX                  ,
     M    H2O_SNOW_NEXT         , ENG_SNOW_NEXT    , 
     O    BERG                  , ENH_BERG  )
!
! 氷河流出量の見積り
!   バッファ層が水当量で 10m を越えた分について、バッファ層から流出させる。
!
      USE SIBPRM , ONLY : 
     1     IDX  ,
     1     ISN  
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, INA ,
     1                       IMONIT_LEVEL, IMONIT_ALL, IMONIT_FLUX_MRI
#endif

!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN) :: NMAX 
!
      REAL(8),INTENT(INOUT) :: H2O_SNOW_NEXT(IDX,ISN)
      REAL(8),INTENT(INOUT) :: ENG_SNOW_NEXT(IDX,ISN)
!
      REAL(8),INTENT(OUT)   :: BERG      (IDX)
      REAL(8),INTENT(OUT)   :: ENH_BERG  (IDX)
!
      INTEGER           :: I
      INTEGER,PARAMETER :: KLEV_BUF = 2
      REAL(8),PARAMETER :: H2O_MAX = 1.D4   ! 1.D4kg = 水当量 10m 
      REAL(8)           :: FACT 
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK
      REAL(8),PARAMETER :: ONE = 1.D0 
#endif
!
      BERG    (:) = 0. 
      ENH_BERG(:) = 0. 
!
      DO I=1,NMAX
        IF ( H2O_SNOW_NEXT(I,KLEV_BUF) .GT. H2O_MAX ) THEN
          FACT = H2O_MAX / H2O_SNOW_NEXT(I,KLEV_BUF)
          BERG    (I) = H2O_SNOW_NEXT(I,KLEV_BUF) * ( 1.D0 -FACT ) 
          ENH_BERG(I) = ENG_SNOW_NEXT(I,KLEV_BUF) * ( 1.D0 -FACT ) 
          H2O_SNOW_NEXT(I,KLEV_BUF) 
     1                = H2O_SNOW_NEXT(I,KLEV_BUF) * FACT  
          ENG_SNOW_NEXT(I,KLEV_BUF) 
     1                = ENG_SNOW_NEXT(I,KLEV_BUF) * FACT  
        ENDIF  
      ENDDO
!
#ifdef MONYOS
      IF (      IMONIT_LEVEL .GE. IMONIT_FLUX_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL      ) THEN
        CMARK = 'BERG' 
        CALL MONIT_SIB ( CMARK , BERG     , INA , ONE , ONE )
      ENDIF

      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        CMARK = 'JFLBERG' 
        CALL MONIT_SIB ( CMARK , ENH_BERG , INA , ONE , ONE )
      ENDIF 
#endif
!
      RETURN
      END 
!=================================================================
      SUBROUTINE MR_SNOW_BERG 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_FLUX_MRI   
#endif
!
      IMPLICIT NONE
!
#ifdef MONYOS
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
#endif
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_BERG 01/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
!
      IF (      IMONIT_LEVEL .GE. IMONIT_FLUX_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL      ) THEN
!
! 氷河流出
!
        CMARK  = 'BERG'
        CTITLE = 'Runoff of Ice Sheet'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 流出氷河の氷に伴うエネルギー
!
        CMARK  = 'JFLBERG'
        CTITLE = 'Energy of Runoff Ice'
        CUNIT  = 'W/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END

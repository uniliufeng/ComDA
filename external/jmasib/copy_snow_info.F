      SUBROUTINE COPY_SNOW_INFO (
     I   NMAX          , NMAX_01       , 
#ifdef MONYOS
     I   NMAX_0        , NMAX_1        , 
#endif 
     I   CVR_SNW       , DEP_SNOW      ,
     O   CVR_SNW_NEXT  , DEP_SNOW_NEXT  )
!
      USE SIBPRM , ONLY :
     1   IDX  ,
     1   ISN   
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY :
     1   CNOS , 
     1   CSNW , 
     1   CNOS_NEXT , 
     1   CSNW_NEXT , 
     1   INF_SNW_NEXT
#endif
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NMAX_01(3)
#ifdef MONYOS
      INTEGER,INTENT(IN)    :: NMAX_0(3)
      INTEGER,INTENT(IN)    :: NMAX_1(3)
#endif
!
      REAL(8),INTENT(IN)    :: CVR_SNW (IDX)
      REAL(8),INTENT(IN)    :: DEP_SNOW(IDX,ISN+1)
!
      REAL(8),INTENT(OUT)   :: CVR_SNW_NEXT (IDX)
      REAL(8),INTENT(OUT)   :: DEP_SNOW_NEXT(IDX,ISN+1)
!
      INTEGER :: I
!
      CALL RESET ( CVR_SNW_NEXT  , IDX )
      CALL RESET ( DEP_SNOW_NEXT , IDX*(ISN+1) )
!
      DO I=NMAX_01(2),NMAX
        CVR_SNW_NEXT  (I)    = CVR_SNW  (I)
        DEP_SNOW_NEXT (I,1)  = DEP_SNOW (I,1)
        DEP_SNOW_NEXT (I,2)  = DEP_SNOW (I,2)
        DEP_SNOW_NEXT (I,3)  = DEP_SNOW (I,3)
        DEP_SNOW_NEXT (I,4)  = DEP_SNOW (I,4)
        DEP_SNOW_NEXT (I,5)  = DEP_SNOW (I,5)
      ENDDO
!
#ifdef MONYOS
!     write(6,*) 'copy_snow_info ' , csnw
!     write(6,*) 'copy_snow_info ' , csnw_next 
      DO I=1,NMAX
        CSNW_NEXT  (I) = CSNW(I)  
        CNOS_NEXT  (I) = CNOS(I)  
      ENDDO
!
      DO I=1,NMAX_0 (3)
        INF_SNW_NEXT (I) = -1
      ENDDO
      DO I=NMAX_01(2),NMAX_01(3)
        INF_SNW_NEXT (I) = 0
      ENDDO
      DO I=NMAX_1(2),NMAX
        INF_SNW_NEXT (I) = 1
      ENDDO
#endif
!
      RETURN
      END


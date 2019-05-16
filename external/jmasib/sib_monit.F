      MODULE SIB_MONIT  
!
      USE SIBPRM , ONLY :
     1   IDX 
!
      IMPLICIT NONE
!
      REAL(8),SAVE  ::  CNOS (IDX) 
      REAL(8),SAVE  ::  CSNW (IDX) 
      INTEGER,SAVE  ::  IMSK (IDX)      ! SIB -> 物理
      INTEGER,SAVE  ::  N01  (3) 
      INTEGER,SAVE  ::  NMX  
      INTEGER,SAVE  ::  JLABEL  
!
! NEXT 
!
      REAL(8),SAVE  :: CNOS_NEXT   (IDX)   
      REAL(8),SAVE  :: CSNW_NEXT   (IDX)   
      INTEGER,SAVE  :: INF_SNW_NEXT(IDX) ! -1 雪無し 0 共存 1 雪あり 
!
      INTEGER,SAVE  :: ICN =  0       ! 雪無し域のみ加算、被覆率かける        
      INTEGER,SAVE  :: ICS =  1       ! 雪あり域のみ加算、被覆率かける   
      INTEGER,SAVE  :: INN =  2       ! 雪無し域のみ加算、被覆率かけない
      INTEGER,SAVE  :: INS =  3       ! 雪あり域のみ加算、被覆率かけない   
      INTEGER,SAVE  :: INA = -1       ! 全域加算、被覆率かけない   
!
! 20050517
!      INTEGER,SAVE  :: IMONIT_SELECT
!      INTEGER,SAVE  :: IMONIT_FLUX_COMMON 
!      INTEGER,SAVE  :: IMONIT_FLUX_MRI 
!      INTEGER,SAVE  :: IMONIT_MEAN_MRI
!      INTEGER,SAVE  :: IMONIT_SNAP_NPD
!      INTEGER,SAVE  :: IMONIT_BUDJET
!      INTEGER,SAVE  :: IMONIT_ALL
!      INTEGER,SAVE  :: IMONIT_ALMA
!	INTEGER,SAVE  :: IMONIT_LEVEL
!
#if   ( SIB0109_OFFLINE_MRI == 1 ) 
      INTEGER,SAVE  :: IMONIT_SELECT      = 0    ! 今注目しているもの
      INTEGER,SAVE  :: IMONIT_FLUX_COMMON = 1    ! 
      INTEGER,SAVE  :: IMONIT_FLUX_MRI    = 1    ! 
      INTEGER,SAVE  :: IMONIT_MEAN_MRI    = 1    ! 
      INTEGER,SAVE  :: IMONIT_SNAP_NPD    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_BUDJET      = 3  
      INTEGER,SAVE  :: IMONIT_ALL         = 9    ! 全て
      INTEGER,SAVE  :: IMONIT_ALMA        = 1000 ! 
!
      INTEGER,SAVE  :: IMONIT_LEVEL       = 1    ! デフォルト
!
#elif ( SIB0109_MJ98_MRI == 1 ) 
      INTEGER,SAVE  :: IMONIT_SELECT      = 0    ! 今注目しているもの
      INTEGER,SAVE  :: IMONIT_FLUX_COMMON = 1    ! 
      INTEGER,SAVE  :: IMONIT_FLUX_MRI    = 1    ! 
      INTEGER,SAVE  :: IMONIT_MEAN_MRI    = 1    ! 
      INTEGER,SAVE  :: IMONIT_SNAP_NPD    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_BUDJET      = 3  
      INTEGER,SAVE  :: IMONIT_ALMA        = 1000 ! 
      INTEGER,SAVE  :: IMONIT_ALL         = 9    ! 全て
!
      INTEGER,SAVE  :: IMONIT_LEVEL       = 1    ! デフォルト
!
#elif ( SIB0109_OFFLINE_JMA == 1 ) 
      INTEGER,SAVE  :: IMONIT_SELECT      = 0    ! 今注目しているもの
      INTEGER,SAVE  :: IMONIT_FLUX_COMMON = 1    ! 
      INTEGER,SAVE  :: IMONIT_FLUX_MRI    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_MEAN_MRI    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_SNAP_NPD    = 1    ! 
      INTEGER,SAVE  :: IMONIT_BUDJET      = 3  
      INTEGER,SAVE  :: IMONIT_ALMA        = 1000 ! 
      INTEGER,SAVE  :: IMONIT_ALL         = 9    ! 全て
!
      INTEGER,SAVE  :: IMONIT_LEVEL       = 1    ! デフォルト
!
#elif ( SIB0109_GSM_JMA == 1 ) 
      INTEGER,SAVE  :: IMONIT_SELECT      = 0    ! 今注目しているもの
      INTEGER,SAVE  :: IMONIT_FLUX_COMMON = 1    ! 
      INTEGER,SAVE  :: IMONIT_FLUX_MRI    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_MEAN_MRI    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_SNAP_NPD    = 1    ! 
      INTEGER,SAVE  :: IMONIT_BUDJET      = 3  
      INTEGER,SAVE  :: IMONIT_ALMA        = 1000 ! 
      INTEGER,SAVE  :: IMONIT_ALL         = 9    ! 全て
!
      INTEGER,SAVE  :: IMONIT_LEVEL       = 1    ! デフォルト
!
#elif ( SIB0109_OFFLINE_ALMA == 1 ) 
      INTEGER,SAVE  :: IMONIT_SELECT      = 0    ! 今注目しているもの
      INTEGER,SAVE  :: IMONIT_FLUX_COMMON = 2    ! 
      INTEGER,SAVE  :: IMONIT_FLUX_MRI    = 2    ! 
      INTEGER,SAVE  :: IMONIT_MEAN_MRI    = 2    ! 
      INTEGER,SAVE  :: IMONIT_SNAP_NPD    = 1000 ! 
      INTEGER,SAVE  :: IMONIT_BUDJET      = 3    !
      INTEGER,SAVE  :: IMONIT_ALMA        = 1    ! 
      INTEGER,SAVE  :: IMONIT_ALL         = 9    ! 全て
!
      INTEGER,SAVE  :: IMONIT_LEVEL       = 1    ! デフォルト

#endif
!
      CONTAINS 
!
! =====================================================================
!
      SUBROUTINE SIB_MONIT_INI 
!
      IMPLICIT NONE
!
      NAMELIST /SIB0109_SIB_MONIT_INI/ IMONIT_LEVEL
!
      NAMELIST /SIB0109_SIB_MONIT_INI_TMP/ 
     1   IMONIT_SELECT      
     1 , IMONIT_FLUX_COMMON 
     1 , IMONIT_FLUX_MRI    
     1 , IMONIT_MEAN_MRI    
     1 , IMONIT_SNAP_NPD    
     1 , IMONIT_BUDJET      
     1 , IMONIT_ALL         
!
      READ (5,SIB0109_SIB_MONIT_INI) 
      WRITE(6,SIB0109_SIB_MONIT_INI) 
!
      WRITE(6,SIB0109_SIB_MONIT_INI_TMP) 
!
      RETURN
      END SUBROUTINE SIB_MONIT_INI 
!
! =====================================================================
      SUBROUTINE SIB_MONIT_RUN_SET_SNOWINFO (
     I   NMAX    , 
     I   CVR_SNW , CVR_NOS , IMSKLT , NMAX_01 , NB )
!   
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NMAX 
      INTEGER,INTENT(IN)   :: NB 
      INTEGER,INTENT(IN)   :: IMSKLT ( IDX ) 
      INTEGER,INTENT(IN)   :: NMAX_01(3)
      REAL(8),INTENT(IN)   :: CVR_SNW(IDX)
      REAL(8),INTENT(IN)   :: CVR_NOS(IDX)
!
      INTEGER              :: I 
!
      DO I=1,NMAX
        CSNW(I) = CVR_SNW(I)
        CNOS(I) = CVR_NOS(I)
        IMSK(I) = IMSKLT (I)
      ENDDO
      N01(1) = NMAX_01(1)
      N01(2) = NMAX_01(2)
      N01(3) = NMAX_01(3)
      NMX    = NMAX
      JLABEL = NB
!
      END SUBROUTINE SIB_MONIT_RUN_SET_SNOWINFO
! =====================================================================
      SUBROUTINE SIB_MONIT_RUN_SET_SNOWNEXTINFO (
     I   NMAX    , 
     I   CVR_SNW , CVR_NOS , INFO_SNOW )
!
      USE SIBPRM , ONLY : ISN
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: INFO_SNOW(IDX)
      REAL(8),INTENT(IN)   :: CVR_SNW  (IDX) 
      REAL(8),INTENT(IN)   :: CVR_NOS  (IDX)
!
      INTEGER :: I
! 
      DO I=1,NMAX
        CSNW_NEXT    (I) = CVR_SNW(I)  
        CNOS_NEXT    (I) = CVR_NOS(I)  
        IF     ( INFO_SNOW(I).LE. -1 ) THEN
          INF_SNW_NEXT (I) = -1 
        ELSEIF ( INFO_SNOW(I).EQ.  0 ) THEN
          INF_SNW_NEXT (I) =  0
        ELSEIF (        INFO_SNOW(I) .GE. 1
     1           .AND.  INFO_SNOW(I) .LE. ISN ) THEN
          INF_SNW_NEXT (I) =  1
        ELSE 
          WRITE(6,*) 'SIB_MONIT_RUN_SET_SNOWNEXTINFO ERROR ' ,
     1               ' INFO_SNOW_NEXT ' , 
     1               I , INFO_SNOW(I) 
          STOP 999
        ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE SIB_MONIT_RUN_SET_SNOWNEXTINFO 
!
! =====================================================================
!
      SUBROUTINE MONIT_SIB ( CMARK , DAT , INSFLAG , DELT , FACT )
!
      IMPLICIT NONE 
!
      CHARACTER(7),INTENT(IN) :: CMARK 
      INTEGER,INTENT(IN) :: INSFLAG 
      REAL(8),INTENT(IN) :: DAT  ( IDX )
      REAL(8),INTENT(IN) :: DELT 
      REAL(8),INTENT(IN) :: FACT 
!
      REAL(8)            :: WORK ( IDX )
      REAL(8),PARAMETER  :: ONE = 1.D0 
!
      INTEGER :: I  
!
      WORK(:) = 0.D0 
!
      IF     ( INSFLAG.EQ.ICN ) THEN           ! 雪無し , 面積比をかける
#include "indep_fixedform"
        DO I=1,N01(3)
          WORK(IMSK(I)) = DAT(I) * CNOS(I) * FACT 
        ENDDO
!
      ELSEIF ( INSFLAG.EQ.ICS ) THEN           ! 雪あり , 面積比をかける
#include "indep_fixedform"
        DO I=N01(2),NMX 
          WORK(IMSK(I)) = DAT(I) * CSNW(I) * FACT 
        ENDDO
!
      ELSEIF  ( INSFLAG.EQ.INN ) THEN          ! 雪無し、面積比はかけない
#include "indep_fixedform"
        DO I=1,N01(3)
          WORK(IMSK(I)) = DAT(I) * FACT 
        ENDDO
!
      ELSEIF ( INSFLAG.EQ.INS ) THEN           ! 雪あり 、面積比はかけない
#include "indep_fixedform"
        DO I=N01(2),NMX 
          WORK(IMSK(I)) = DAT(I) * FACT 
        ENDDO
!
      ELSEIF ( INSFLAG.EQ.INA ) THEN           ! どちらも。面積比かけない
#include "indep_fixedform"
        DO I=1,NMX
          WORK(IMSK(I)) = DAT(I) * FACT 
        ENDDO
#ifdef SIB_DEBUG
      ELSE
        write(6,*) 'monit_sib flag error' , insflag 
        stop 999
#endif
      ENDIF
!     
! 20050804  
!      CALL MONIT_ADD_2 ( CMARK , WORK , JLABEL , DELT , ONE )
!
      RETURN
      END SUBROUTINE MONIT_SIB
!===========================================================
      SUBROUTINE MONIT_SIB_NEXT ( CMARK , DAT , INSFLAG , DELT , FACT )
!
      IMPLICIT NONE 
!
      CHARACTER(7),INTENT(IN) :: CMARK 
      INTEGER,INTENT(IN) :: INSFLAG 
      REAL(8),INTENT(IN) :: DAT  ( IDX )
      REAL(8),INTENT(IN) :: DELT 
      REAL(8),INTENT(IN) :: FACT 
!
      REAL(8)            :: WORK ( IDX )
      REAL(8),PARAMETER  :: ONE = 1.D0 
!
      INTEGER :: I  
!
      CALL RESET ( WORK , IDX )
!
      IF     ( INSFLAG.EQ.ICN ) THEN           ! 雪無し
#include "indep_fixedform"
        DO I=1,NMX
          IF ( INF_SNW_NEXT(I).LE.0) THEN
            WORK(IMSK(I)) = DAT(I) * CNOS_NEXT(I) * FACT 
          ENDIF
        ENDDO 
!
      ELSEIF ( INSFLAG.EQ.ICS ) THEN           ! 雪あり
#include "indep_fixedform"
        DO I=1,NMX 
          IF ( INF_SNW_NEXT(I).GE.0) THEN
            WORK(IMSK(I)) = DAT(I) * CSNW_NEXT(I) * FACT 
          ENDIF
        ENDDO
!
      ELSEIF ( INSFLAG.EQ.INN ) THEN           ! 雪無し
#include "indep_fixedform"
        DO I=1,NMX
          IF ( INF_SNW_NEXT(I).LE.0) THEN
            WORK(IMSK(I)) = DAT(I) * FACT 
          ENDIF
        ENDDO 
!
      ELSEIF ( INSFLAG.EQ.INS ) THEN           ! 雪あり
#include "indep_fixedform"
        DO I=1,NMX 
          IF ( INF_SNW_NEXT(I).GE.0) THEN
            WORK(IMSK(I)) = DAT(I) * FACT 
          ENDIF
        ENDDO
!
      ELSEIF ( INSFLAG.EQ.INA ) THEN           ! どちらも
#include "indep_fixedform"
        DO I=1,NMX
          WORK(IMSK(I)) = DAT(I) * FACT 
        ENDDO
!
#ifdef SIB_DEBUG
      ELSE
        write(6,*) 'monit_sib_next flag error' , insflag 
        stop 999
#endif
!
      ENDIF
!
      CALL MONIT_ADD_2 ( CMARK , WORK , JLABEL , DELT , ONE )
!
      RETURN
      END SUBROUTINE MONIT_SIB_NEXT
!
      END MODULE SIB_MONIT

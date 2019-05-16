! 土壌の平均値計算 （主にモニタ用）
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_AREA_MEAN (
     I   NMAX         , 
     I   CVR_SNW     , 
     I   RUNOFF_B_NOS , RUNOFF_B_SNW , 
     I   RUNOFF_S_NOS , RUNOFF_S_SNW , 
     I   TMP_SOIL4_N   , TMP_SOIL4_S   , 
     O   RUNOFF        , TMP_SOIL4 ) 
!
      USE SIBPRM , ONLY : 
     1   IDX     , 
     1   IDP
!
      USE SIBCON, ONLY : CAP_SOIL4_c 
#ifdef MONYOS
      use sib_monit, only : monit_sib, imonit_level, imonit_all, INA
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX 
      REAL(8),INTENT(IN)    :: CVR_SNW     (IDX)
      REAL(8),INTENT(IN)    :: TMP_SOIL4_N  (IDX)
      REAL(8),INTENT(IN)    :: TMP_SOIL4_S  (IDX)
      REAL(8),INTENT(IN)    :: RUNOFF_B_NOS (IDX) !ランオフ(水フラックス)
      REAL(8),INTENT(IN)    :: RUNOFF_B_SNW (IDX) !ランオフ(水フラックス)
      REAL(8),INTENT(IN)    :: RUNOFF_S_NOS (IDX) !ランオフ(水フラックス)
      REAL(8),INTENT(IN)    :: RUNOFF_S_SNW (IDX) !ランオフ(水フラックス)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: RUNOFF   (IDX)     !ランオフ(水フラックス)
      REAL(8),INTENT(OUT)   :: TMP_SOIL4 (IDX) 
!
! WORK
!
      INTEGER   ::  I
!
#ifdef MONYOS
      character(7) :: cmark 
      real(8) :: work(idx)
      real(8),parameter :: one = 1.D0 
#endif

!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_AREA_MEAN 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=1,NMAX
        RUNOFF  (I) =   RUNOFF_S_NOS(I) * ( 1 - CVR_SNW(I) )
     1                + RUNOFF_S_SNW(I) *  CVR_SNW(I) 
     1                + RUNOFF_B_NOS(I) * ( 1 - CVR_SNW(I) )
     1                + RUNOFF_B_SNW(I) *  CVR_SNW(I) 
!
        TMP_SOIL4(I) =  TMP_SOIL4_N (I) * ( 1 - CVR_SNW(I) )
     1                + TMP_SOIL4_S (I) *  CVR_SNW(I) 
 1000 CONTINUE  
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! エネルギーの移動 （雪なしからありへ）
!
        CMARK = 'HL4N2S1' 
        CALL RESET ( WORK , IDX )  
        DO I=1,NMAX  
          WORK (I) = CAP_SOIL4_c * ( TMP_SOIL4_N(I) - TMP_SOIL4(I) ) 
     1                           * ( 1 - CVR_SNW(I) ) 
        ENDDO
        CALL MONIT_SIB ( CMARK , WORK , INA , ONE , ONE )
      ENDIF
#endif
!
      RETURN
      END
!============================================================
      subroutine MR_soil_area_mean 
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
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SOIL_AREA_MEAN 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
      CMARK  = 'HL4N2S1'
      CTITLE = 'ENERGY MOVE SOIL 1 FROM NOS TO SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      end subroutine MR_soil_area_mean 


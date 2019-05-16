! 雪面被覆率変更
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_SNOWCVR_MOD (
     I   NMAX             , 
     I   CVR_SNW          , CVR_SNW_NEXT     , 
     I   CVR_NOS_INV_NEXT , CVR_SNW_INV_NEXT , 
     I   H2O_SOIL_SNW     , H2O_SOIL_NOS     ,
     I   ENG_SOIL_SNW     , ENG_SOIL_NOS     , 
     I   TMP_SOIL4        , 
     O   H2O_SOIL_NEXT    , H2O_SNSL_NEXT     ,
     O   ENG_SOIL_NEXT    , ENG_SNSL_NEXT       )
!
      USE SIBPRM , ONLY :
     1   IDX   , 
     1   IDP
!
      use sibcon, only : cap_soil4_c,tmp_frez_c
!
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, INA , 
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE

!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX   
      REAL(8),INTENT(IN)   :: CVR_SNW      (IDX)
      REAL(8),INTENT(IN)   :: CVR_SNW_NEXT (IDX)
      REAL(8),INTENT(IN)   :: CVR_NOS_INV_NEXT (IDX)
      REAL(8),INTENT(IN)   :: CVR_SNW_INV_NEXT (IDX)
!
      REAL(8),INTENT(IN)   :: H2O_SOIL_NOS  (IDX,IDP)
      REAL(8),INTENT(IN)   :: ENG_SOIL_NOS  (IDX,IDP)
      REAL(8),INTENT(IN)   :: H2O_SOIL_SNW  (IDX,IDP)
      REAL(8),INTENT(IN)   :: ENG_SOIL_SNW  (IDX,IDP)
!
      REAL(8),INTENT(IN)   :: TMP_SOIL4     (IDX)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: H2O_SOIL_NEXT (IDX,IDP)
      REAL(8),INTENT(OUT)  :: ENG_SOIL_NEXT (IDX,IDP)
      REAL(8),INTENT(OUT)  :: H2O_SNSL_NEXT (IDX,IDP)
      REAL(8),INTENT(OUT)  :: ENG_SNSL_NEXT (IDX,IDP)
!
! WORK
!
      INTEGER :: I
      INTEGER :: K
      REAL(8) :: FAC_NOS_NOS(IDX)
      REAL(8) :: FAC_SNW_SNW(IDX)
!
      LOGICAL,SAVE :: LFIRST
!
#ifdef MONYOS
      REAL(8) ::  WORK                (IDX)
      CHARACTER(7) :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
#endif

      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_SNOWCVR_MOD 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!   雪の被覆率が増える場合、
!     雪がない場所の ENG , H2O はそのまま
!     雪がある場所には、ない場所の寄与あり
!   雪の被覆率が減る場合、
!     雪がある場所の ENG , H2O はそのまま
!     雪がない場所には、ある場所の寄与あり
!
!  ・以下の操作で、エネルギーと H2O は不変であることに注意せよ。
!      
      DO I=1,NMAX
        IF     ( CVR_SNW_NEXT(I) .EQ. CVR_SNW(I) ) THEN
          FAC_NOS_NOS(I) = 1.                                ! 無し -> 無し
          FAC_SNW_SNW(I) = 1.                                ! あり -> あり
        ELSEIF ( CVR_SNW_NEXT(I) .GT. CVR_SNW(I) ) THEN  ! 雪が増える
          FAC_NOS_NOS(I) = 1.                                ! 無し -> 無し
          FAC_SNW_SNW(I) = CVR_SNW(I) * CVR_SNW_INV_NEXT(I)  ! あり -> あり
        ELSEIF ( CVR_SNW_NEXT(I) .LT. CVR_SNW(I) ) THEN  ! 雪が減る
          FAC_NOS_NOS(I) = ( 1 - CVR_SNW(I) ) 
     1                       * CVR_NOS_INV_NEXT(I)           ! 無し -> 無し
          FAC_SNW_SNW(I) = 1.                                ! あり -> あり
        ENDIF
      ENDDO
!
#include "noparallel_fixedform"
      DO K=1,IDP
      DO I=1,NMAX
        ENG_SOIL_NEXT(I,K)                                   
     1       =   ENG_SOIL_NOS(I,K) * FAC_NOS_NOS(I)          
     1         + ENG_SOIL_SNW(I,K) * ( 1 - FAC_NOS_NOS(I) )  
        ENG_SNSL_NEXT(I,K)                                     
     1       =   ENG_SOIL_NOS(I,K) * ( 1 - FAC_SNW_SNW(I) )  
     1         + ENG_SOIL_SNW(I,K) * FAC_SNW_SNW(I)          
!
        H2O_SOIL_NEXT(I,K) 
     1       =   H2O_SOIL_NOS(I,K) * FAC_NOS_NOS(I) 
     1         + H2O_SOIL_SNW(I,K) * ( 1 - FAC_NOS_NOS(I) ) 
        H2O_SNSL_NEXT(I,K) 
     1       =   H2O_SOIL_NOS(I,K) * ( 1 - FAC_SNW_SNW(I) )
     1         + H2O_SOIL_SNW(I,K) * FAC_SNW_SNW(I) 
      ENDDO
      ENDDO
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! エネルギーの移動 （雪なしからありへ）
!
        CMARK = 'HLxN2S' 
#include "noparallel_fixedform"
        DO K=1,3
          CALL RESET ( WORK , IDX )  
          DO I=1,NMAX  
            IF     ( CVR_SNW_NEXT(I) .GT. CVR_SNW(I) ) THEN  ! 雪が増える
              WORK (I) = ENG_SOIL_NOS(I,K) 
     1                      * ( CVR_SNW_NEXT(I) - CVR_SNW(I) ) 
            ELSEIF ( CVR_SNW_NEXT(I) .LT. CVR_SNW(I) ) THEN  ! 雪が減る
              WORK (I) = ENG_SOIL_SNW(I,K) 
     1                      * ( CVR_SNW_NEXT(I) - CVR_SNW(I) ) 
            ENDIF
          ENDDO
!
          WRITE(CMARK(3:3),'(I1)' ) K 
          CALL MONIT_SIB ( CMARK , WORK , INA , ONE , ONE )
        ENDDO
!
        CALL RESET ( WORK , IDX )  
        DO I=1,NMAX  
          WORK (I) = (TMP_SOIL4(I)-TMP_FREZ_c) * CAP_SOIL4_c 
     1                      * ( CVR_SNW_NEXT(I) - CVR_SNW(I) ) 
        ENDDO
        CMARK = 'HL4N2S2' 
        CALL MONIT_SIB ( CMARK , WORK , INA , ONE , ONE )
!
! H2O の移動 （雪なしからありへ）
!
        CMARK = 'WLxN2S' 
#include "noparallel_fixedform"
        DO K=1,3
          CALL RESET ( WORK , IDX )  
          DO I=1,NMAX  
            IF     ( CVR_SNW_NEXT(I) .GT. CVR_SNW(I) ) THEN  ! 雪が増える
              WORK (I) = H2O_SOIL_NOS(I,K) 
     1                      * ( CVR_SNW_NEXT(I) - CVR_SNW(I) ) 
            ELSEIF ( CVR_SNW_NEXT(I) .LT. CVR_SNW(I) ) THEN  ! 雪が減る
              WORK (I) = H2O_SOIL_SNW(I,K) 
     1                      * ( CVR_SNW_NEXT(I) - CVR_SNW(I) ) 
            ENDIF
          ENDDO
!
          WRITE(CMARK(3:3),'(I1)' ) K 
          CALL MONIT_SIB ( CMARK , WORK , INA , ONE , ONE )
        ENDDO
!
      ENDIF
#endif
!
      RETURN
      END
!================================================================
      SUBROUTINE MR_SOIL_SNOWCVR_MOD
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
        WRITE(6,*) 'MONIT_REGIST_SOIL_SNOWCVR_MOD 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
      CMARK  = 'HL1N2S'
      CTITLE = 'ENERGY MOVE SOIL 1 FROM NOS TO SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HL2N2S'
      CTITLE = 'ENERGY MOVE SOIL 2 FROM NOS TO SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HL3N2S'
      CTITLE = 'ENERGY MOVE SOIL 3 FROM NOS TO SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HL4N2S2'
      CTITLE = 'ENERGY MOVE SOIL 4 FROM NOS TO SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WL1N2S'
      CTITLE = 'H2O MOVE SOIL 1 FROM NOS TO SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WL2N2S'
      CTITLE = 'H2O MOVE SOIL 2 FROM NOS TO SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WL3N2S'
      CTITLE = 'H2O MOVE SOIL 3 FROM NOS TO SNW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif 
!
      RETURN
      END

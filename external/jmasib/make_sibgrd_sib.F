! SiB の全面配列から、SiB で必要なデータの計算・作成
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス

      SUBROUTINE MAKE_SIBGRD_SIB (  
     I  NMAX     , 
     I  NMAX_0   , NMAX_01 , NMAX_1  , NMAX_2   , NMAX_3  , NMAX_4   , 
     I  NTYP     , MON     ,
     I  WTR_SNOW        , GLA_SNOW       , 
     I  RHO_SNOW_INV    ,  
     O  CVR_VEG_ms      , CVR_VEG_INV_ms ,   
     O  CVR_SNW         , CVR_SNW_INV    , 
     O  CVR_NOS         , CVR_NOS_INV    , 
     O  DEP_SNOW      )  
!
      USE SIBPRM , ONLY : 
     1   IDX    , 
     1   ICG    , 
     1   IDP    , 
     1   ISN 
!
      USE SIBCON , ONLY :
     1   H2O_SNOW_01_1_INV_c    ,      ! これ以下なら雪 01 
     1   CVR_SNW_INV_MAX_c          ,
     1   CVR_VEG_MIN_c 
!
      USE SIBVEG , ONLY : 
     1           CVR_VEG_mo               ! 被覆率 (ITYP,IMON,ICG)
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: MON 
      INTEGER,INTENT(IN)   :: NTYP (IDX) 
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_0 (3)  
      INTEGER,INTENT(IN)   :: NMAX_01(3)  
      INTEGER,INTENT(IN)   :: NMAX_1(3)  
      INTEGER,INTENT(IN)   :: NMAX_2(3)  
      INTEGER,INTENT(IN)   :: NMAX_3(3)  
      INTEGER,INTENT(IN)   :: NMAX_4(3)  
!
      REAL(8),INTENT(IN)   ::  RHO_SNOW_INV      (IDX,ISN)
      REAL(8),INTENT(IN)   ::  WTR_SNOW          (IDX,ISN)
      REAL(8),INTENT(IN)   ::  GLA_SNOW          (IDX,ISN)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  ::  CVR_SNW       (IDX)
      REAL(8),INTENT(OUT)  ::  CVR_SNW_INV   (IDX)
      REAL(8),INTENT(OUT)  ::  CVR_NOS       (IDX)
      REAL(8),INTENT(OUT)  ::  CVR_NOS_INV   (IDX)
      REAL(8),INTENT(OUT)  ::  DEP_SNOW      (IDX,ISN+1)
!
      REAL(8),INTENT(OUT)  ::  CVR_VEG_ms    (IDX,ICG)
      REAL(8),INTENT(OUT)  ::  CVR_VEG_INV_ms(IDX,ICG)
!
! WORK
!
      INTEGER  ::  I 
      INTEGER  ::  J 
      REAL(8)  ::  TMP 
!
!     DEP_SNOW 
!     CVR_SNW 
!
      LOGICAL LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MAKE_SIBGRID_SIB 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF
!
      CALL RESET ( DEP_SNOW , IDX *(ISN+1) )
!
! ---------------
! > 雪の被覆率 <
! ---------------
!
      DO 3000 I=1,NMAX_0(3)
        CVR_SNW(I)     = 0.
        CVR_SNW_INV(I) = CVR_SNW_INV_MAX_c 
        CVR_NOS(I)     = 1.
        CVR_NOS_INV(I) = 1. 
 3000 CONTINUE  
!
      DO 3100 I=NMAX_01(2),NMAX_01(3)
        CVR_SNW(I) = ( WTR_SNOW(I,1) + GLA_SNOW(I,1) ) 
     1                    * H2O_SNOW_01_1_INV_c
        CVR_SNW_INV(I) 
     1      = MIN ( CVR_SNW_INV_MAX_c , 1/CVR_SNW(I) )
        CVR_NOS(I) = 1. - CVR_SNW(I) 
        CVR_NOS_INV(I)  
     1      = MIN ( CVR_SNW_INV_MAX_c , 1/CVR_NOS(I) )
 3100 CONTINUE  
!
      DO 3200 I=NMAX_1(2),NMAX
        CVR_SNW    (I) = 1.
        CVR_SNW_INV(I) = 1.
        CVR_NOS(I)     = 0.
        CVR_NOS_INV(I) = CVR_SNW_INV_MAX_c 
 3200 CONTINUE  
!
! ------------
! > 雪の深さ <
! ------------
!
      DO 2100 I=NMAX_01(2),NMAX   
        DEP_SNOW(I,1) 
     1       = GLA_SNOW(I,1) * RHO_SNOW_INV(I,1)
     1                       * CVR_SNW_INV(I) 
        DEP_SNOW(I,5) = DEP_SNOW(I,1)
 2100 CONTINUE  
!
      DO 2200 I=NMAX_2(2),NMAX   
        DEP_SNOW(I,2) 
     1       = GLA_SNOW(I,2) * RHO_SNOW_INV(I,2)
     1                       * CVR_SNW_INV(I) 
        DEP_SNOW(I,5) = DEP_SNOW(I,5) + DEP_SNOW(I,2)
 2200 CONTINUE  
!
      DO 2300 I=NMAX_3(2),NMAX   
        DEP_SNOW(I,3) 
     1       = GLA_SNOW(I,3) * RHO_SNOW_INV(I,3)
     1                       * CVR_SNW_INV(I) 
        DEP_SNOW(I,5) = DEP_SNOW(I,5) + DEP_SNOW(I,3)
 2300 CONTINUE  
!
      DO 2400 I=NMAX_4(2),NMAX   
        DEP_SNOW(I,4) 
     1       = GLA_SNOW(I,4) * RHO_SNOW_INV(I,4)
     1                       * CVR_SNW_INV(I) 
        DEP_SNOW(I,5) = DEP_SNOW(I,5) + DEP_SNOW(I,4)
 2400 CONTINUE  
!
! --------------
! > 植生被覆率 < 
! --------------
!
      DO 4000 I=1,NMAX 
        J = NTYP(I) 
        CVR_VEG_ms(I,1) = CVR_VEG_mo(J,MON,1)
        CVR_VEG_ms(I,2) = CVR_VEG_mo(J,MON,2)
!
        TMP = MAX ( CVR_VEG_ms(I,1) , CVR_VEG_MIN_c ) 
        CVR_VEG_INV_ms(I,1) = 1 / TMP 
        TMP = MAX ( CVR_VEG_ms(I,2) , CVR_VEG_MIN_c ) 
        CVR_VEG_INV_ms(I,2) = 1 / TMP 
 4000 CONTINUE  
!
!
      RETURN
      END

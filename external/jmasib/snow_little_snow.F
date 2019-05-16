! 新雪と旧雪の量のバランスのチェック
! 雪の量が少ないときは、少雪バケツにいれる。
!      
      SUBROUTINE SNOW_LITTLE_SNOW (
     I    NMAX                  ,
     M    INFO_SNOW_NEXT        , 
     M    H2O_SNOW_NEXT         , ENG_SNOW_NEXT    , 
     M    AGE_SNOW_NEXT         , 
     O    ENG_SNOW_BUCKET_NEXT  , H2O_SNOW_BUCKET_NEXT   ) 
!
      USE SIBPRM , ONLY :
     1    IDX  , 
     1    ISN
!  
!     USE MESSAGE
!
      USE SIBCON , ONLY : 
     1    AGE_SNOW_MAX2_c     ,
     1    H2O_SNOW_MIN_c       
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICS, INS , INA,
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    ::  NMAX
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)    :: H2O_SNOW_BUCKET_NEXT(IDX)
      REAL(8),INTENT(OUT)    :: ENG_SNOW_BUCKET_NEXT(IDX)

!
! INOUT
!
      REAL(8),INTENT(INOUT) ::  H2O_SNOW_NEXT   (IDX,ISN)
      REAL(8),INTENT(INOUT) ::  ENG_SNOW_NEXT   (IDX,ISN)
      REAL(8),INTENT(INOUT) ::  AGE_SNOW_NEXT   (IDX)
      INTEGER,INTENT(INOUT) ::  INFO_SNOW_NEXT  (IDX)
!
! WORK
      INTEGER  :: I 
!
#ifdef MONYOS
      REAL(8),PARAMETER  :: ONE = 1.D0 
!     REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_LITTLE_SNOW POCHAKA 01/08/21' 
        LFIRST = .FALSE.
      ENDIF
!
      CALL RESET ( ENG_SNOW_BUCKET_NEXT , IDX )
      CALL RESET ( H2O_SNOW_BUCKET_NEXT , IDX )
!
!  新雪と旧雪の和がとても小さい場合、バケツに格納する
!       
! VOPTION = 条件を満たす場合のみの実行
      DO 1000 I=1,NMAX
!
        IF ( INFO_SNOW_NEXT (I) .EQ. 1 ) THEN
        IF ( H2O_SNOW_NEXT(I,1) .LT. H2O_SNOW_MIN_c ) THEN
!
!  総和が旧雪の下限 1 よりも少なければ、バケツにいれる
! 
          AGE_SNOW_NEXT(I)         = AGE_SNOW_MAX2_c 
          H2O_SNOW_BUCKET_NEXT (I) = H2O_SNOW_NEXT(I,1) 
          ENG_SNOW_BUCKET_NEXT (I) = ENG_SNOW_NEXT(I,1) 
          H2O_SNOW_NEXT(I,1) = 0.        
          ENG_SNOW_NEXT(I,1) = 0.        
!
          INFO_SNOW_NEXT(I)  = -2 
        ENDIF
        ENDIF
 1000 CONTINUE 
!
! ==============
! > 詳細モニタ <
! ==============
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 小雪バケツ  
!
        CMARK  = 'WLWBPUT' 
        CALL MONIT_SIB ( CMARK , H2O_SNOW_BUCKET_NEXT ,INA, ONE , ONE )
!
        CMARK  = 'HLWBPUT' 
        CALL MONIT_SIB ( CMARK , ENG_SNOW_BUCKET_NEXT ,INA, ONE , ONE )
!
      ENDIF
#endif
!
      RETURN
      END
!============================================================
      SUBROUTINE MR_SNOW_LITTLE_SNOW
!
! SNOW_WFLUX_PRE で出力するデバッグ用変数の登録作業
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
!
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_BALANCE_CHECK 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! 少雪バケツ
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'WLWBPUT'
      CTITLE = 'H2O PUT INTO LITTLE SNOW BUCKET '
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLWBPUT'
      CTITLE = 'ENG WITH H2O OUT : LITTLE SNOW BUCKET '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

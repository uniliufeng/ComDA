      SUBROUTINE SNOW_PRE (
     I    NMAX      , NMAX_01   , 
     I    NMAX_2    , NMAX_3    , NMAX_4  , 
     I    CVR_SNW_INV    , 
     I    WTR_SNOW_IN    , GLA_SNOW_IN    , 
     I    H2O_BUCKET_FROM_SNOW ,  ENH_BUCKET_FROM_SNOW , 
     I    TMP_SNOW_IN     , 
     O    TMP_SNOW       , WTR_SNOW       , GLA_SNOW      )   
!
      USE SIBPRM , ONLY : 
     I    IDX     , 
     I    ISN
!
      USE SIBCON , ONLY : 
     1    CAP_WTR_c  ,                            ! CW 
     1    CAP_GLA_c  ,                            ! CICE 
     1    TMP_FREZ_c , 
     1    HL_W2G_c   
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , INS ,
     1                IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX  
      INTEGER,INTENT(IN)    :: NMAX_01(3)  
      INTEGER,INTENT(IN)    :: NMAX_2(3)  
      INTEGER,INTENT(IN)    :: NMAX_3(3)    
      INTEGER,INTENT(IN)    :: NMAX_4(3)  
!
!     REAL(8),INTENT(IN)    :: DELT_SNOW_INV 
!
      REAL(8),INTENT(IN)    :: WTR_SNOW_IN  (IDX,ISN)
      REAL(8),INTENT(IN)    :: GLA_SNOW_IN  (IDX,ISN)
!
      REAL(8),INTENT(IN)    :: CVR_SNW_INV  (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_SNOW_IN  (IDX,ISN+1)
!
      REAL(8),INTENT(IN)    :: H2O_BUCKET_FROM_SNOW (IDX)
      REAL(8),INTENT(IN)    :: ENH_BUCKET_FROM_SNOW (IDX) 
!
!      REAL(8),INTENT(IN)    :: DELT_SNOW   
!
! OUTPUT
!
      REAL(8),INTENT(OUT)   :: TMP_SNOW     (IDX,ISN+1)
      REAL(8),INTENT(OUT)   :: WTR_SNOW     (IDX,ISN)
      REAL(8),INTENT(OUT)   :: GLA_SNOW     (IDX,ISN)
!
! WORK
! 
      REAL(8) :: ENG_NEW(IDX)
      REAL(8) :: ENX_NEW(IDX)
      REAL(8) :: H2O_NEW(IDX)
      REAL(8) :: FR_W
      REAL(8) :: CAP
!
      INTEGER :: I
!
#ifdef MONYOS
!     REAL(8),PARAMETER  :: DM1 = -1.D0 
      REAL(8),PARAMETER  :: ONE = 1.D0 
!     REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SNOW_PRE 99/10/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  --------------
!  > 積雪第一層 <
!  --------------
!    雪あり域平均にする
!
      DO 1000 I=NMAX_01(2),NMAX
        TMP_SNOW(I,2) = TMP_SNOW_IN(I,2) 
!
        ENG_NEW(I) = ( 
     1                   (   CAP_WTR_c * WTR_SNOW_IN(I,1) 
     1                     + CAP_GLA_c * GLA_SNOW_IN(I,1) )
     1                   * ( TMP_SNOW_IN(I,1) - TMP_FREZ_c  ) 
     1                   - HL_W2G_c * GLA_SNOW_IN(I,1)
     1                   + ENH_BUCKET_FROM_SNOW (I) 
     1               ) * CVR_SNW_INV(I) 
!
        H2O_NEW(I) = (
     1                   WTR_SNOW_IN(I,1) 
     1                 + GLA_SNOW_IN(I,1) 
     1                 + H2O_BUCKET_FROM_SNOW(I) 
     1               ) * CVR_SNW_INV(I) 

        ENX_NEW(I) = ENG_NEW(I) + HL_W2G_c * H2O_NEW(I)
!
#ifdef SIB_DEBUG
        if ( H2O_NEW(I) .LT. 0. ) THEN
          write(6,*) 'snow_pre error : h2o_new < 0 ' , i , 
     1                   'sum' , h2o_new(i),
     1                   'W ' , WTR_SNOW_IN       (I,1)      ,
     1                   'WB' , H2O_BUCKET_FROM_SNOW(I)      ,
     1                   'G ' , GLA_SNOW_IN       (I,1) 
          stop 999
        endif
#endif

 1000 CONTINUE  
!
      DO 1100 I=NMAX_01(2),NMAX
        IF     ( ENX_NEW(I) .LT. 0. ) THEN 
          CAP = CAP_GLA_c * H2O_NEW(I)
          TMP_SNOW(I,1) = TMP_FREZ_c + ENX_NEW(I) / CAP 
          GLA_SNOW(I,1) = H2O_NEW(I)
          WTR_SNOW(I,1) = 0. 
        ELSE
          TMP_SNOW(I,1) = TMP_FREZ_c 
          FR_W = ENX_NEW(I) / ( ENX_NEW(I) - ENG_NEW(I) )  ! 水の割合
          GLA_SNOW(I,1) = H2O_NEW(I) * ( 1 - FR_W )
          WTR_SNOW(I,1) = H2O_NEW(I) *  FR_W 
        ENDIF
!
#ifdef SIB_DEBUG
        if ( GLA_SNOW(I,1).LT.0. .OR. WTR_SNOW(I,1).LT.0. ) THEN
          write(6,*) 'snow_pre error : 1 ', i , 
     1                   'W ' , WTR_SNOW  (I,1)      ,
     1                   'G ' , GLA_SNOW  (I,1)     
          stop 999
        endif
#endif

 1100 CONTINUE     
!
!  --------------
!  > 積雪第二層 <
!  --------------
!
      DO 2000 I=NMAX_2(2),NMAX
        TMP_SNOW(I,3) = TMP_SNOW_IN(I,3) 
!
        WTR_SNOW(I,2) = WTR_SNOW_IN(I,2) 
        GLA_SNOW(I,2) = GLA_SNOW_IN(I,2) 
!
#ifdef SIB_DEBUG
        if ( GLA_SNOW(I,2).LT.0. .OR. WTR_SNOW(I,2).LT.0. ) THEN
          write(6,*) 'snow_pre error : 2 ', i , 
     1                   'W ' , WTR_SNOW  (I,2)      ,
     1                   'G ' , GLA_SNOW  (I,2)     
          stop 999
        endif
#endif

 2000 CONTINUE  
!
!  --------------
!  > 積雪第三層 <
!  --------------
!
      DO 3000 I=NMAX_3(2),NMAX
        TMP_SNOW(I,4) = TMP_SNOW_IN(I,4) 
!
        WTR_SNOW(I,3) = WTR_SNOW_IN(I,3) 
        GLA_SNOW(I,3) = GLA_SNOW_IN(I,3) 
!
#ifdef SIB_DEBUG
        if ( GLA_SNOW(I,3).LT.0. .OR. WTR_SNOW(I,3).LT.0. ) THEN
          write(6,*) 'snow_pre error : 3 ', i , 
     1                   'W ' , WTR_SNOW  (I,3)      ,
     1                   'G ' , GLA_SNOW  (I,3)     
          stop 999
        endif
#endif

 3000 CONTINUE  
!
!  --------------
!  > 積雪第四層 <
!  --------------
!
      DO 4000 I=NMAX_4(2),NMAX
        TMP_SNOW(I,5) = TMP_SNOW_IN(I,5) 
!
        WTR_SNOW(I,4) = WTR_SNOW_IN(I,4) 
        GLA_SNOW(I,4) = GLA_SNOW_IN(I,4) 
!
#ifdef SIB_DEBUG
        if ( GLA_SNOW(I,4).LT.0. .OR. WTR_SNOW(I,4).LT.0. ) THEN
          write(6,*) 'snow_pre error : 4 ', i , 
     1                   'W ' , WTR_SNOW  (I,4)      ,
     1                   'G ' , GLA_SNOW  (I,4)     
          stop 999
        endif
#endif
!
 4000 CONTINUE  
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 蒸発・昇華する水の量
!
        CMARK  = 'WEVBWS'       ! 蒸発 From snoW Mean
        CALL MONIT_SIB ( CMARK, H2O_BUCKET_FROM_SNOW ,INS, ONE , ONE ) 
!
! 蒸発・昇華する水に伴うエネルギー量
!
        CMARK  = 'IEVBWS'
        CALL MONIT_SIB ( CMARK,ENH_BUCKET_FROM_SNOW,INS,ONE,ONE)  ! ver2 
!
      ENDIF
#endif 
!
      RETURN
      END
!==================================================================
      SUBROUTINE MR_SNOW_PRE
!
! SNOW_PRE で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                IMONIT_LEVEL , IMONIT_ALL
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
        WRITE(6,*) 'MONIT_REGIST_SNOW_PRE 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 蒸発する水の量
!
      CMARK  = 'WEVBWS'
      CTITLE = 'EVAP WATER FROM SNOW'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 蒸発する水に伴う熱
!
      CMARK  = 'IEVBWS'
      CTITLE = 'ENERGY WITH EVAP WATER FROM SNOW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

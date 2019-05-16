      SUBROUTINE MONIT_SNAP_CANOPY (
     I    NMAX    , MON    , NTYP   , 
     I    DELT_CNP, 
     I    DEP_SNOW_NEXT,
     I    TMP_CNP_NOS_NEXT    , TMP_CNP_SNW_NEXT    , 
     I    FR_WTR_CNP_NOS_NEXT , FR_WTR_CNP_SNW_NEXT , 
     I    FR_GLA_CNP_NOS_NEXT , FR_GLA_CNP_SNW_NEXT , 
     I    TMP_GRSK_NEXT       , TMP_SNSK_NEXT       , 
     I    FR_WTR_GRS_NEXT     , FR_GLA_GRS_NEXT      
     1   )  
! 
      USE SIBPRM , ONLY : 
     1    IDX , 
     1    ICG ,
     1    ISN
!
      USE SIBCON , ONLY :
     1    CAP_GRSK_c  ,
     1    CAP_SNSK_c  ,
     1    CAP_LAI_c  ,
     1    CAP_WTR_c  ,
     1    CAP_GLA_c  , 
     1    TMP_FREZ_c , 
     1    HL_W2G_c   , 
     1    RHO_WTR_c  , 
     1    RHO_GLA_c  , 
     1    H2O_SAT_CNP_DEP_c ,
     1    H2O_SAT_GRS_DEP_c
!
      USE SIBVEG , ONLY :
     1           Z2_mo                  , ! キャノピー上端高さ (ITYP,IMON)
     1           Z1_mo                  , ! キャノピー下端高さ (ITYP,IMON)
     1           XLAI_mo                  ! LAI, 格子平均だ!! (ITYP,IMON,ICG)
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB_NEXT, ICN, ICS,INA ,
     1                      IMONIT_LEVEL , IMONIT_ALL , IMONIT_ALMA
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif
!
      IMPLICIT NONE
!
! INPUT
!
!     INTEGER,INTENT(IN)    :: IMEAN            ! 1 平均値も出す
!
      INTEGER,INTENT(IN)    :: NMAX       
      INTEGER,INTENT(IN)    :: MON
      INTEGER,INTENT(IN)    :: NTYP(IDX)
!
!     INTEGER,INTENT(IN)    :: NMAX_01(3)
      REAL(8),INTENT(IN)    :: DELT_CNP
!
      REAL(8),INTENT(IN)    :: TMP_CNP_NOS_NEXT     (IDX) 
      REAL(8),INTENT(IN)    :: TMP_CNP_SNW_NEXT     (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_CNP_NOS_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_CNP_SNW_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_CNP_NOS_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_CNP_SNW_NEXT  (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_GRSK_NEXT        (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_GRS_NEXT      (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_GRS_NEXT      (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_SNSK_NEXT        (IDX) 
!
      REAL(8),INTENT(IN)    :: DEP_SNOW_NEXT        (IDX,ISN+1) 

!
! WORK
!     
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE  = 1. 
      CHARACTER(7) :: CMARK 
!     REAL(8)      :: WORK  (IDX)
      REAL(8)      :: WORK_N(IDX)
      REAL(8)      :: WORK_S(IDX)
      REAL(8)      :: WTR_N(IDX)
      REAL(8)      :: WTR_S(IDX)
      REAL(8)      :: GLA_N(IDX)
      REAL(8)      :: GLA_S(IDX)
#endif
!
      INTEGER   ::  I
      REAL(8)   ::  XZ 
      REAL(8)   ::  SDEP 
      REAL(8),PARAMETER  ::  ZERO = 0.D0
      REAL(8)   ::  FAC 
!
      REAL(8)   :: XLAI_VEG_ms  (IDX,ICG)
      REAL(8)   :: XLAI_VEG_msx (IDX,ICG)
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_SNAP_CANOPY 99/10 POCHAKA'

        LFIRST = .FALSE. 
      ENDIF
!
! ================
! >> モニタ準備 <<
! ================
!
      DO I = 1,NMAX 
        XLAI_VEG_ms  (I,1)    = XLAI_mo (NTYP(I),MON,1) 
        XLAI_VEG_ms  (I,2)    = XLAI_mo (NTYP(I),MON,2) 
!
!#ifdef CORRECT 
        XZ     = Z2_mo  ( NTYP(I) , MON ) - Z1_mo  ( NTYP(I) , MON )  
        SDEP   = MIN ( DEP_SNOW_NEXT(I,5) - Z1_mo  ( NTYP(I) , MON )  
     1                  , XZ * 0.95D0 )
!#else
!        XZ     = Z2_mo  ( NTYP(I) , MON )
!        SDEP   = MIN ( DEP_SNOW_NEXT(I,5)  , XZ * 0.95D0 )  
!#endif
        SDEP   = MAX ( SDEP , ZERO ) 
        FAC    = 1.0D0 - SDEP / XZ            
!
        XLAI_VEG_msx (I,1)    = XLAI_mo (NTYP(I),MON,1) * FAC
        XLAI_VEG_msx (I,2)    = 0. 
      ENDDO
!
! ==============================
! >> キャノピー収支チェック用 <<
! ==============================
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1 
     1   ) THEN
!
! 液水量
!
      DO I=1,NMAX
        WTR_N(I) = FR_WTR_CNP_NOS_NEXT(I)
     1             * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
      ENDDO 
      CMARK  = 'WTRCNX'
      CALL MONIT_SIB_NEXT ( CMARK , WTR_N ,ICN, DELT_CNP, ONE)      
!
      DO I=1,NMAX
        WTR_S(I) = FR_WTR_CNP_SNW_NEXT(I)
     1             * XLAI_VEG_msx(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
      ENDDO 
      CMARK  = 'WTRCSX'
      CALL MONIT_SIB_NEXT ( CMARK , WTR_S ,ICS, DELT_CNP, ONE)      
!
! 氷量
!
      DO I=1,NMAX
        GLA_N(I) = FR_GLA_CNP_NOS_NEXT(I)
     1             * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
      ENDDO 
      CMARK  = 'GLACNX'
      CALL MONIT_SIB_NEXT ( CMARK , GLA_N ,ICN, DELT_CNP, ONE)      
!
      DO I=1,NMAX
        GLA_S(I) = FR_GLA_CNP_SNW_NEXT(I)
     1             * XLAI_VEG_msx(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
      ENDDO 
      CMARK  = 'GLACSX'
      CALL MONIT_SIB_NEXT ( CMARK , GLA_S ,ICS, DELT_CNP, ONE)      
!
! エネルギー
!
      DO I=1,NMAX
        WORK_N(I) =
     1             (   CAP_LAI_c * XLAI_VEG_ms(I,1) 
     1               + CAP_WTR_c * WTR_N(I)
     1               + CAP_GLA_c * GLA_N(I) ) 
     1             * ( TMP_CNP_NOS_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_N(I) 
      ENDDO
      CMARK  = 'ENGCNX'
      CALL MONIT_SIB_NEXT ( CMARK , WORK_N ,ICN , DELT_CNP, ONE)      
!
      DO I=1,NMAX
        WORK_S(I) =
     1             (   CAP_LAI_c * XLAI_VEG_msx(I,1) 
     1               + CAP_WTR_c * WTR_S(I)
     1               + CAP_GLA_c * GLA_S(I) ) 
     1             * ( TMP_CNP_SNW_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_S(I) 
!
      ENDDO 
      CMARK  = 'ENGCSX'
      CALL MONIT_SIB_NEXT ( CMARK , WORK_S ,ICS, DELT_CNP, ONE)      
!
! ==============================
! >> 下草・地面収支チェック用 <<
! ==============================
!
! 下草液水量
!
      DO I=1,NMAX
        WTR_N(I) = FR_WTR_GRS_NEXT(I)
     1             * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
      ENDDO 
!
      CMARK  = 'WTRGRX'
      CALL MONIT_SIB_NEXT ( CMARK , WTR_N ,ICN, DELT_CNP, ONE)      
!
! 下草氷量
!
      DO I=1,NMAX
        GLA_N(I) = FR_GLA_GRS_NEXT(I)
     1             * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
      ENDDO 
!
      CMARK  = 'GLAGRX'
      CALL MONIT_SIB_NEXT ( CMARK , GLA_N ,ICN, DELT_CNP, ONE)      
!
! 下草・土壌スキンエネルギー
!
      DO I=1,NMAX
        WORK_N(I) =
     1             (   CAP_GRSK_c      
     1               + CAP_LAI_c * XLAI_VEG_ms(I,2) 
     1               + CAP_WTR_c * WTR_N(I)
     1               + CAP_GLA_c * GLA_N(I) ) 
     1             * ( TMP_GRSK_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_N(I) 
!    1             - HL_W2G_c * CAP_GLA_c * GLA_N(I) 
      ENDDO 
!
      CMARK  = 'ENGGAX'
      CALL MONIT_SIB_NEXT ( CMARK , WORK_N ,ICN, DELT_CNP, ONE)      
!
! ========================
! >> 雪面収支チェック用 <<
! ========================
!
! 雪面エネルギー
!
      DO I=1,NMAX
        WORK_N(I) =
     1                CAP_SNSK_c 
     1             * ( TMP_SNSK_NEXT(I) - TMP_FREZ_c )
      ENDDO 
!
      CMARK  = 'ENGW0X'
      CALL MONIT_SIB_NEXT ( CMARK , WORK_N ,ICS, DELT_CNP, ONE)      
!
      ENDIF
#endif
!
      RETURN
      END
!=======================================================================
      SUBROUTINE MONIT_SNAP_CANOPY_2 (
     I    NMAX    , MON    , NTYP   , IMEAN , 
     I    DELT_CNP, 
     I    CVR_SNW ,  
     I    DEP_SNOW, 
     I    TMP_CNP_NOS_NEXT    , TMP_CNP_SNW_NEXT    , 
     I    FR_WTR_CNP_NOS_NEXT , FR_WTR_CNP_SNW_NEXT , 
     I    FR_GLA_CNP_NOS_NEXT , FR_GLA_CNP_SNW_NEXT , 
     I    TMP_GRSK_NEXT       , TMP_SNSK_NEXT       , 
     I    FR_WTR_GRS_NEXT     , FR_GLA_GRS_NEXT      
     1   )  
!
      USE SIBPRM , ONLY : 
     1    IDX , 
     1    ICG ,
     1    ISN
!
      USE SIBCON , ONLY :
     1    CAP_GRSK_c  ,
     1    CAP_SNSK_c  ,
     1    CAP_LAI_c  ,
     1    CAP_WTR_c  ,
     1    CAP_GLA_c  , 
     1    TMP_FREZ_c , 
     1    HL_W2G_c   , 
     1    RHO_WTR_c  , 
     1    RHO_GLA_c  , 
     1    H2O_SAT_CNP_DEP_c ,
     1    H2O_SAT_GRS_DEP_c
!
      USE SIBVEG , ONLY :
     1           Z2_mo                  , ! キャノピー上端高さ (ITYP,IMON)
     1           Z1_mo                  , ! キャノピー下端高さ (ITYP,IMON)
     1           XLAI_mo                  ! LAI, 格子平均だ!! (ITYP,IMON,ICG)
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB ,INA , ICN,ICS,
     1                      IMONIT_LEVEL , IMONIT_ALL ,
     1                      IMONIT_MEAN_MRI  ,
     1                      IMONIT_SNAP_NPD , IMONIT_ALMA
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: IMEAN             ! 1 平均値も出す
!
      INTEGER,INTENT(IN)    :: NMAX       
      INTEGER,INTENT(IN)    :: MON
      INTEGER,INTENT(IN)    :: NTYP(IDX)
!
      REAL(8),INTENT(IN)    :: CVR_SNW    (IDX) 
      REAL(8),INTENT(IN)    :: DELT_CNP
!
      REAL(8),INTENT(IN)    :: TMP_CNP_NOS_NEXT     (IDX) 
      REAL(8),INTENT(IN)    :: TMP_CNP_SNW_NEXT     (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_CNP_NOS_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_CNP_SNW_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_CNP_NOS_NEXT  (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_CNP_SNW_NEXT  (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_GRSK_NEXT        (IDX) 
      REAL(8),INTENT(IN)    :: FR_WTR_GRS_NEXT      (IDX) 
      REAL(8),INTENT(IN)    :: FR_GLA_GRS_NEXT      (IDX) 
!
      REAL(8),INTENT(IN)    :: TMP_SNSK_NEXT        (IDX) 
      REAL(8),INTENT(IN)    :: DEP_SNOW             (IDX,ISN+1) 
!
! WORK
!     
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      REAL(8),PARAMETER :: D100 = 100.D0 
      CHARACTER(7) :: CMARK 
      REAL(8)      :: WORK  (IDX)
      REAL(8)      :: WORK_N(IDX)
      REAL(8)      :: WORK_S(IDX)
      REAL(8)      :: WTR_N(IDX)
      REAL(8)      :: WTR_S(IDX)
      REAL(8)      :: GLA_N(IDX)
      REAL(8)      :: GLA_S(IDX)
#endif
!
      INTEGER   ::  I
      REAL(8)   ::  XZ 
      REAL(8),PARAMETER  ::  ZERO = 0.D0
      REAL(8)   ::  FAC 
      REAL(8)   ::  SDEP 
!
      REAL(8)   :: XLAI_VEG_ms  (IDX,ICG)
      REAL(8)   :: XLAI_VEG_msx (IDX,ICG)
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_SNAP_CANOPY_2 99/10 POCHAKA'

        LFIRST = .FALSE. 
      ENDIF
!
! ================
! >> モニタ準備 <<
! ================
!
      DO I = 1,NMAX 
        XLAI_VEG_ms  (I,1)    = XLAI_mo (NTYP(I),MON,1) 
        XLAI_VEG_ms  (I,2)    = XLAI_mo (NTYP(I),MON,2) 
!
        XZ     = Z2_mo  ( NTYP(I) , MON ) - Z1_mo  ( NTYP(I) , MON )  
        SDEP   = MIN ( DEP_SNOW(I,5) - Z1_mo  ( NTYP(I) , MON )  
     1                  , XZ * 0.95D0 )
        SDEP   = MAX ( SDEP , ZERO ) 
        FAC    = 1.0D0 - SDEP / XZ            
!
        XLAI_VEG_msx (I,1)    = XLAI_mo (NTYP(I),MON,1) * FAC
        XLAI_VEG_msx (I,2)    = 0. 
      ENDDO
!
#ifdef MONYOS
!
! ==========================
! >> キャノピー・格子平均 <<
! ==========================
!
! --------------------
! > 気象研・時間平均 <
! --------------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
      IF ( IMEAN .EQ. 1 ) THEN
!
! 温度
!
        CMARK  = 'TSC'
        CALL MONIT_SIB ( CMARK , TMP_CNP_NOS_NEXT,ICN, DELT_CNP, ONE)      
        CALL MONIT_SIB ( CMARK , TMP_CNP_SNW_NEXT,ICS, DELT_CNP, ONE)      
!
! 水率
!
        CMARK  = 'MSC'
        CALL MONIT_SIB (CMARK,FR_WTR_CNP_NOS_NEXT,ICN,DELT_CNP,ONE)      
        CALL MONIT_SIB (CMARK,FR_WTR_CNP_SNW_NEXT,ICS,DELT_CNP,ONE)      
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_NOS_NEXT,ICN,DELT_CNP,ONE)      
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_SNW_NEXT,ICS,DELT_CNP,ONE)      
!
! 積雪被覆率
!
        CMARK  = 'CVRSNWC'
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_NOS_NEXT,ICN,DELT_CNP,ONE)      
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_SNW_NEXT,ICS,DELT_CNP,ONE)      
!
      ENDIF
      ENDIF
!
! ----------------
! > 本庁・瞬間値 <
! ----------------
!
      IF ( 
     1           IMONIT_LEVEL .GE. IMONIT_SNAP_NPD
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1 
     1   ) THEN
!
! 温度
!
        CMARK  = 'TSC'
        CALL MONIT_SIB ( CMARK , TMP_CNP_NOS_NEXT,ICN, DELT_CNP, ONE)      
        CALL MONIT_SIB ( CMARK , TMP_CNP_SNW_NEXT,ICS, DELT_CNP, ONE)      
!
! 水率
!
        CMARK  = 'MSC'
        CALL MONIT_SIB (CMARK,FR_WTR_CNP_NOS_NEXT,ICN,DELT_CNP,D100)      
        CALL MONIT_SIB (CMARK,FR_WTR_CNP_SNW_NEXT,ICS,DELT_CNP,D100)      
!
! 氷率
!
        CMARK  = 'ISC'
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_NOS_NEXT,ICN,DELT_CNP,D100)      
        CALL MONIT_SIB (CMARK,FR_GLA_CNP_SNW_NEXT,ICS,DELT_CNP,D100)      
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ========================
! >> キャノピー収支以外 <<
! ========================
!
! 温度（スナップ）
!
        IF ( ICN_SIB0109_MONITOR .EQ. 1 ) THEN
          CMARK  = 'TMPCNV'
          CALL MONIT_SIB ( CMARK , TMP_CNP_NOS_NEXT,ICN, DELT_CNP, ONE)      
!
          CMARK  = 'TMPCSV'
          CALL MONIT_SIB ( CMARK , TMP_CNP_SNW_NEXT,ICS, DELT_CNP, ONE)      
        ENDIF
!
! ----------
! > 平均値 <
! ----------
!
        IF ( IMEAN .EQ. 1 ) THEN
!
! 温度
!
          CMARK  = 'TMPCNL'
          CALL MONIT_SIB ( CMARK , TMP_CNP_NOS_NEXT,ICN, DELT_CNP, ONE)      
          CMARK  = 'TMPCSL'
          CALL MONIT_SIB ( CMARK , TMP_CNP_SNW_NEXT,ICS, DELT_CNP, ONE)      
!
! 液水量
!
          DO I=1,NMAX
            WTR_N(I) = FR_WTR_CNP_NOS_NEXT(I)
     1                * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c      
            WTR_S(I) = FR_WTR_CNP_SNW_NEXT(I)
     1                * XLAI_VEG_msx(I,1)* H2O_SAT_CNP_DEP_c * RHO_WTR_c 
          ENDDO 

          CMARK  = 'WTRCNL'
          CALL MONIT_SIB (CMARK,WTR_N,ICN,DELT_CNP,ONE)      
          CMARK  = 'WTRCSL'
          CALL MONIT_SIB (CMARK,WTR_S,ICS,DELT_CNP,ONE)      
!
! 液水率
!
          CMARK  = 'FRWCNL'
          CALL MONIT_SIB (CMARK,FR_WTR_CNP_NOS_NEXT,ICN,DELT_CNP,ONE)      
          CMARK  = 'FRWCSL'
          CALL MONIT_SIB (CMARK,FR_WTR_CNP_SNW_NEXT,ICS,DELT_CNP,ONE)      
!
! 氷量
!
          DO I=1,NMAX
            GLA_N(I) = FR_GLA_CNP_NOS_NEXT(I)
     1                * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
            GLA_S(I) = FR_GLA_CNP_SNW_NEXT(I)
     1                * XLAI_VEG_msx(I,1)* H2O_SAT_CNP_DEP_c * RHO_WTR_c 
          ENDDO 
!
          CMARK  = 'GLACNL'
          CALL MONIT_SIB (CMARK,GLA_N,ICN,DELT_CNP,ONE)      
          CMARK  = 'GLACSL'
          CALL MONIT_SIB (CMARK,GLA_S,ICS,DELT_CNP,ONE)      
!
! 氷率
!
          CMARK  = 'FRGCNL'
          CALL MONIT_SIB (CMARK,FR_GLA_CNP_NOS_NEXT,ICN,DELT_CNP,ONE)      
          CMARK  = 'FRGCSL'
          CALL MONIT_SIB (CMARK,FR_GLA_CNP_SNW_NEXT,ICS,DELT_CNP,ONE)      
!
! エネルギー
!
          DO I=1,NMAX
            WORK_N(I) =
     1                 (   CAP_LAI_c * XLAI_VEG_ms(I,1) 
     1                   + CAP_WTR_c * WTR_N(I)
     1                   + CAP_GLA_c * GLA_N(I) ) 
     1                 * ( TMP_CNP_NOS_NEXT(I) - TMP_FREZ_c )
     1                 - HL_W2G_c * GLA_N(I) 
!
            WORK_S(I) =
     1                 (   CAP_LAI_c * XLAI_VEG_msx(I,1) 
     1                   + CAP_WTR_c * WTR_S(I)
     1                   + CAP_GLA_c * GLA_S(I) ) 
     1                 * ( TMP_CNP_SNW_NEXT(I) - TMP_FREZ_c )
     1                 - HL_W2G_c * GLA_S(I) 
          ENDDO 
!
          CMARK  = 'ENGCNL'
          CALL MONIT_SIB (CMARK,WORK_N,ICN,DELT_CNP,ONE)      
          CMARK  = 'ENGCSL'
          CALL MONIT_SIB (CMARK,WORK_S,ICS,DELT_CNP,ONE)      
!
        ENDIF
      ENDIF
!
! ==========================
! >> 下草・地面の格子平均 <<
! ==========================
!
! ----------------
! > 気象研平均値 <
! ----------------
!
      IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( IMEAN .EQ. 1 ) THEN
!
! 積雪被覆率
!
          CMARK  = 'CVRSNWL'
          CALL MONIT_SIB ( CMARK , CVR_SNW ,INA, DELT_CNP, ONE)      
        ENDIF
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
        IF ( IMEAN .EQ. 1 ) THEN
!
! 温度
!
          CMARK  = 'TSG'
          CALL MONIT_SIB ( CMARK , TMP_GRSK_NEXT,ICN, DELT_CNP, ONE)      
          CALL MONIT_SIB ( CMARK , TMP_SNSK_NEXT,ICS, DELT_CNP, ONE)      
!
! 下草・水率
!
          CMARK  = 'MSG'
          CALL MONIT_SIB (CMARK,FR_WTR_GRS_NEXT,ICN,DELT_CNP,ONE)      
          CALL MONIT_SIB (CMARK,FR_GLA_GRS_NEXT,ICN,DELT_CNP,ONE)      
        ENDIF
      ENDIF
!
! --------------
! > 本庁瞬間値 <
! --------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SNAP_NPD 
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1 
     1     ) THEN
!
! 温度
!
        CMARK  = 'TSG'
        CALL MONIT_SIB ( CMARK , TMP_GRSK_NEXT,ICN, DELT_CNP, ONE)      
        CALL MONIT_SIB ( CMARK , TMP_SNSK_NEXT,ICS, DELT_CNP, ONE)      
!
! 下草・水率
!
        CMARK  = 'MSG'
        CALL MONIT_SIB (CMARK,FR_WTR_GRS_NEXT,ICN,DELT_CNP,D100)      
!
! 下草・水率
!
        CMARK  = 'ISG'
        CALL MONIT_SIB (CMARK,FR_GLA_GRS_NEXT,ICN,DELT_CNP,D100)      
      ENDIF
!
! ================================
! >> 下草・地面収支チェック以外 <<
! ================================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( IMEAN .EQ. 1 ) THEN
!
          DO I=1,NMAX
!
! 下草液水量
!
          WTR_N(I) = FR_WTR_GRS_NEXT(I)
     1              * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
!
! 下草氷量
!
          GLA_N(I) = FR_GLA_GRS_NEXT(I)
     1              * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
!
! 下草・土壌スキンエネルギー
!
          WORK_N(I) =
     1             (   CAP_GRSK_c      
     1               + CAP_LAI_c * XLAI_VEG_ms(I,2) 
     1               + CAP_WTR_c * WTR_N(I)
     1               + CAP_GLA_c * GLA_N(I) ) 
     1             * ( TMP_GRSK_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_N(I) 
!    1             - HL_W2G_c * CAP_GLA_c * GLA_N(I) 
          ENDDO 
!
!
! 下草・裸地温度
!
          CMARK  = 'TMPGAL'
          CALL MONIT_SIB ( CMARK , TMP_GRSK_NEXT ,ICN, DELT_CNP, ONE)      
!
! 下草液水量
!
          CMARK  = 'WTRGRL'
          CALL MONIT_SIB ( CMARK , WTR_N ,ICN, DELT_CNP, ONE)      
!
! 下草液水率
!
          CMARK  = 'FRWGRL'
          CALL MONIT_SIB ( CMARK , FR_WTR_GRS_NEXT ,ICN, DELT_CNP, ONE)      
!
! 下草氷量
!
          CMARK  = 'GLAGRL'
          CALL MONIT_SIB ( CMARK , GLA_N ,ICN, DELT_CNP, ONE)      
!
! 下草液水率
!
          CMARK  = 'FRGGRL'
          CALL MONIT_SIB ( CMARK , FR_GLA_GRS_NEXT ,ICN, DELT_CNP, ONE)      
!
! 下草・裸地エネルギー
!
          CMARK  = 'ENGGAL'
          CALL MONIT_SIB ( CMARK , WORK_N ,ICN, DELT_CNP, ONE)      
!
        ENDIF
      ENDIF
!
! ==========
! >> 雪面 <<
! ==========
!
!  -------------------
!  > 格子平均?モニタ <
!  -------------------
!     積雪面温度・時間平均
!      
      IF ( 
     1      IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN
        IF ( IMEAN .EQ. 1 ) THEN
          CMARK  = 'TMPW0L'
          CALL MONIT_SIB (CMARK,TMP_SNSK_NEXT,ICS,DELT_CNP,ONE)      
        ENDIF
      ENDIF
!
! 温度瞬間値
!
      IF ( 
     1       (     IMONIT_LEVEL .GE. IMONIT_SNAP_NPD
     1        .OR. IMONIT_LEVEL .GE. IMONIT_ALL )
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1 
     1   ) THEN
        CMARK  = 'TMPW0V'
        CALL MONIT_SIB ( CMARK , TMP_SNSK_NEXT,ICS, DELT_CNP, ONE)      
      ENDIF
!
!  ----------------
!  > 各領域別出力 <
!  ----------------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!  雪被覆率瞬間値
!
        IF ( ICN_SIB0109_MONITOR .EQ. 1 ) THEN 
          CMARK  = 'CVRSNWV'
          CALL MONIT_SIB ( CMARK , CVR_SNW, INA, DELT_CNP, ONE)      
        ENDIF
!
        IF ( IMEAN .EQ. 1 ) THEN
!
! エネルギー
!
          DO I=1,NMAX
            WORK_N(I) =
     1                    CAP_SNSK_c 
     1                 * ( TMP_SNSK_NEXT(I) - TMP_FREZ_c )
          ENDDO 
!
          CMARK  = 'ENGW0L'
          CALL MONIT_SIB ( CMARK , WORK_N ,ICS, DELT_CNP, ONE)      
!
!
        ENDIF
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          CMARK  = 'A_FACSL'
          CALL MONIT_SIB ( CMARK , CVR_SNW ,INA, DELT_CNP, ONE)      
!
          CMARK  = 'A_TSS' 
          CALL MONIT_SIB ( CMARK , TMP_SNSK_NEXT,ICS, DELT_CNP, ONE)      
!
          CMARK  = 'A_TBS' 
          CALL MONIT_SIB ( CMARK , TMP_GRSK_NEXT,ICN, DELT_CNP, ONE)      
!
          CMARK  = 'A_TCP'
          CALL MONIT_SIB ( CMARK, TMP_CNP_NOS_NEXT,ICN, DELT_CNP, ONE)      
          CALL MONIT_SIB ( CMARK, TMP_CNP_SNW_NEXT,ICS, DELT_CNP, ONE)      
!
          WORK(:) = 0. 
          DO I=1,NMAX
            WORK(I) = FR_GLA_CNP_NOS_NEXT(I)
     1             * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
     1             * (1-CVR_SNW(I)) 
     1              + FR_GLA_CNP_SNW_NEXT(I)
     1             * XLAI_VEG_msx(I,1)* H2O_SAT_CNP_DEP_c * RHO_WTR_c 
     1             * CVR_SNW(I) 
          ENDDO 
          CMARK  = 'A_SWEV'
          CALL MONIT_SIB (CMARK,WORK,INA,DELT_CNP,ONE)      
!
          WORK(:) = 0. 
          DO I=1,NMAX
            WORK(I) = 
     1               + (FR_GLA_CNP_NOS_NEXT(I) + FR_WTR_CNP_NOS_NEXT(I))
     1             * XLAI_VEG_ms(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
     1             * (1-CVR_SNW(I)) 
     1               + (FR_GLA_GRS_NEXT(I) + FR_WTR_GRS_NEXT(I))
     1             * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
     1             * (1-CVR_SNW(I)) 
     1              + (FR_GLA_CNP_SNW_NEXT(I) + FR_WTR_CNP_SNW_NEXT(I))      
     1             * XLAI_VEG_msx(I,1)* H2O_SAT_CNP_DEP_c * RHO_WTR_c 
     1             * CVR_SNW(I) 
          ENDDO 
          CMARK  = 'A_MVEG'
          CALL MONIT_SIB (CMARK,WORK,INA,DELT_CNP,ONE)      

      ENDIF
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_MONIT_SNAP_CANOPY
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_MEAN_MRI,
     1                      IMONIT_SNAP_NPD, IMONIT_ALMA
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
        WRITE(6,*) 'MONIT_REGIST_MONIT_SNAP_CANOPY 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!
#ifdef MONYOS
!
! ==========================
! >> 時間平均・格子平均量 <<
! ==========================
!
!   雪表面温度 ( TSG で代用可能 )
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN
        CMARK  = 'TMPW0L'
        CTITLE = 'SNOW SKIN TEMP (TMEAN)'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN
!
!   雪被覆率 
!
        CMARK  = 'CVRSNWL'
        CTITLE = 'SNOW COVER (TMEAN)'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! ==========================
! >> MRI : 格子平均物理量 <<
! ==========================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN
!
!   キャノピー温度 TSC = TMPCNL + TMPCSL  
        CMARK  = 'TSC'
        CTITLE = 'Canopy Temperature'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   キャノピー水率 MSC = FR_WTR/GLA_CNP_NOS/SNW
        CMARK  = 'MSC'
        CTITLE = 'Canopy Water'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   キャノピー雪被覆率 FR_GLA_CNP_NOS/SNW 
        CMARK  = 'CVRSNWC'
        CTITLE = 'Cover Ratio of Snow on Canopy'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   土壌・下草・雪面温度 TGC = TMPGAL + TMPW0L  
        CMARK  = 'TSG'
        CTITLE = 'Ground Skin Temperature'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!   下草 MSG = FR_WTR/GLA_GRS
        CMARK  = 'MSG'
        CTITLE = 'Grass Water'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! ==========================
! >> スナップ・格子平均量 <<
! ==========================
!
      IF (      IMONIT_LEVEL .GE. IMONIT_SNAP_NPD 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN
!
!  雪表面温度
!
        CMARK  = 'TMPW0V'
        CTITLE = 'SNOW SKIN TEMP (TSNAP)'
        CUNIT  = 'K'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_SNAP_NPD ) THEN
!
!   キャノピー温度 TSC = TMPCNL + TMPCSL  
        CMARK  = 'TSC'
        CTITLE = 'Canopy Temperature'
        CUNIT  = 'K'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   キャノピー水率 MSC = FR_WTR_CNP_NOS/SNW * 100
        CMARK  = 'MSC'
        CTITLE = 'Canopy Water'
        CUNIT  = '%'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   キャノピー氷率 ISC = FR_GLA_CNP_NOS/SNW * 100
        CMARK  = 'ISC'
        CTITLE = 'Canopy Ice'
        CUNIT  = '%'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   下草・地面・雪温度 TSG = TMPGAL + TMPW0L  
        CMARK  = 'TSG'
        CTITLE = 'Ground Temperature'
        CUNIT  = 'K'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   下草水率 MSG = FR_WTR_GRS * 100
        CMARK  = 'MSG'
        CTITLE = 'Grass Water'
        CUNIT  = '%'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!   下草氷率 ISC = FR_GLA_GRS * 100
        CMARK  = 'ISG'
        CTITLE = 'Grass Ice'
        CUNIT  = '%'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDIF
!
! ======================
! >> 領域別モニタ登録 <<
! ======================
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN    
!
! ================================
! >> キャノピー、収支チェック用 <<
! ================================
!
! 水
!
      CMARK  = 'WTRCNX'
      CTITLE = 'CANOPY WATER NOS (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRCSX'
      CTITLE = 'CANOPY WATER SNW (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! 氷
!
      CMARK  = 'GLACNX'
      CTITLE = 'CANOPY GLA NOS (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GLACSX'
      CTITLE = 'CANOPY GLA SNW (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! エネルギー
!
      CMARK  = 'ENGCNX'
      CTITLE = 'CANOPY ENERGY NOS (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ENGCSX'
      CTITLE = 'CANOPY ENERGY SNW (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! ==========================
! >> キャノピー、その他用 <<
! ==========================
!
! 温度
!
      CMARK  = 'TMPCNL'
      CTITLE = 'CANOPY TEMP NOS (TMEAN)'
      CUNIT  = 'K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPCSL'
      CTITLE = 'CANOPY TEMP SNW (TMEAN)'
      CUNIT  = 'K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPCNV'
      CTITLE = 'CANOPY TEMP NOS (TSNAP)'
      CUNIT  = 'K'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPCSV'
      CTITLE = 'CANOPY TEMP SNW (TSNAP)'
      CUNIT  = 'K'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! 液水量
!
      CMARK  = 'WTRCNL'
      CTITLE = 'CANOPY WATER NOS (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WTRCSL'
      CTITLE = 'CANOPY WATER SNW (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 液水率
!
      CMARK  = 'FRWCNL'
      CTITLE = 'FRACTION OF CANOPY FRAC WTR NOS (TMEAN)'
      CUNIT  = ' '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRWCSL'
      CTITLE = 'FRACTION OF CANOPY FRAC WTR SNW (TMEAN)'
      CUNIT  = ' '
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 氷量
!
      CMARK  = 'GLACNL'
      CTITLE = 'CANOPY GLA NOS (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GLACSL'
      CTITLE = 'CANOPY GLA SNW (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 氷率
!
      CMARK  = 'FRGCNL'
      CTITLE = 'FRACTION OF CANOPY FRAC GLA NOS (TMEAN)'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'FRGCSL'
      CTITLE = 'FRACTION OF CANOPY FRAC GLA SNW (TMEAN)'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! エネルギー
!
      CMARK  = 'ENGCNL'
      CTITLE = 'ENERGY OF CANOPY NOS (TMEAN)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ENGCSL'
      CTITLE = 'ENERGY OF CANOPY SNW (TMEAN)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ==========================
! >> 下草、収支チェック用 <<
! ==========================
!
! 水
!
      CMARK  = 'WTRGRX'
      CTITLE = 'GRASS WATER NOS = AREA MEAN (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! 氷
!
      CMARK  = 'GLAGRX'
      CTITLE = 'GRASS GLA NOS = AREA MEAN (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! エネルギー
!
      CMARK  = 'ENGGAX'
      CTITLE = 'GRASS ENERGY NOS = AREA MEAN (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! ==================
! >> 下草、その他 <<
! ==================
!
! 温度
!
      CMARK  = 'TMPGAL'
      CTITLE = 'GRASS TEMP NOS (TMEAN)'
      CUNIT  = 'K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'TMPGAV'
      CTITLE = 'GRASS TEMP NOS (TSNAP)'
      CUNIT  = 'K'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! 水量
!
      CMARK  = 'WTRGRL'
      CTITLE = 'GRASS WATER NOS (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 水率
!
      CMARK  = 'FRWGRL'
      CTITLE = 'FRACTION OF GRASS WATER NOS (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 氷量
!
      CMARK  = 'GLAGRL'
      CTITLE = 'GRASS GLA NOS (TMEAN)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 氷率
!
      CMARK  = 'FRGGRL'
      CTITLE = 'FRACTION OF GRASS GLA NOS (TMEAN)'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! エネルギー
!
      CMARK  = 'ENGGAL'
      CTITLE = 'GRASS BARE ENERGY NOS (TMEAN)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ==========================
! >> 雪面収支関連スナップ <<
! ==========================
!
! 雪面のもつエネルギー 
!
      CMARK  = 'ENGW0X'
      CTITLE = 'SNSK ENERGY (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! ==============
! >> 雪その他 <<
! ==============
!
! 被覆率
!
!
      CMARK  = 'CVRSNWV'
      CTITLE = 'SNOW COVER (TSNAP)'
      CUNIT  = '0-1'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
!
! 雪面のもつエネルギー 
!
      CMARK  = 'ENGW0L'
      CTITLE = 'SNOW SKIN ENERGY'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! ALMA 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
          CMARK  = 'A_FACSL'
          CTITLE = 'ALMA (Not sent) cover snow for surface/canopy calc '
          CUNIT  = '0-1'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
          CMARK  = 'A_TSS' 
          CTITLE = 'ALMA SnowT.     devided by a_facsl'
          CUNIT  = 'K'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
          CMARK  = 'A_TBS' 
          CTITLE = 'ALMA BaresoilT. devided by (1-a_facsl)'
          CUNIT  = 'K'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
          CMARK  = 'A_TCP'
          CTITLE = 'ALMA VegT'
          CUNIT  = 'K'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
          CMARK  = 'A_SWEV'
          CTITLE = 'ALMA SWEVeg (Ice on Canopy)'
          CUNIT  = 'kg/m**2'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
          CMARK  = 'A_MVEG'
          CTITLE = 'ALMA CanopInt (Liquid and Ice on Canopy and Grass)'
          CUNIT  = 'kg/m**2'
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END


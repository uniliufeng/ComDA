      SUBROUTINE TOTAL_H2O_ENG_OF_SKIN (
     I    NMAX    , MON_NEXT  , NTYP   , 
     I    CVR_SNW_NEXT        , 
     I    TMP_CNP_NOS_NEXT    , TMP_CNP_SNW_NEXT    , 
     I    FR_WTR_CNP_NOS_NEXT , FR_WTR_CNP_SNW_NEXT , 
     I    FR_GLA_CNP_NOS_NEXT , FR_GLA_CNP_SNW_NEXT , 
     I    TMP_GRSK_NEXT       , TMP_SNSK_NEXT       , 
     I    FR_WTR_GRS_NEXT     , FR_GLA_GRS_NEXT     , 
     I    DEP_SNOW_NEXT       , 
     O    H2O_CNP_NOS         , H2O_CNP_SNW         ,
     O    H2O_GRSK            , 
     O    ENG_CNP_NOS         , ENG_CNP_SNW         ,
     O    ENG_GRSK            , ENG_SNSK            )  
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_MEAN_MRI,
     1                      IMONIT_SNAP_NPD
#endif
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
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: NMAX       
      INTEGER,INTENT(IN)    :: MON_NEXT
      INTEGER,INTENT(IN)    :: NTYP(IDX)
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
      REAL(8),INTENT(IN)    :: CVR_SNW_NEXT        (IDX) 
!
      REAL(8),INTENT(IN)    :: DEP_SNOW_NEXT        (IDX,ISN+1) 
!
      REAL(8),INTENT(OUT)   :: H2O_CNP_NOS          (IDX) 
      REAL(8),INTENT(OUT)   :: H2O_CNP_SNW          (IDX) 
      REAL(8),INTENT(OUT)   :: H2O_GRSK             (IDX) 
      REAL(8),INTENT(OUT)   :: ENG_CNP_NOS          (IDX) 
      REAL(8),INTENT(OUT)   :: ENG_CNP_SNW          (IDX) 
      REAL(8),INTENT(OUT)   :: ENG_GRSK             (IDX) 
      REAL(8),INTENT(OUT)   :: ENG_SNSK             (IDX) 
!
! WORK
!     
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE  = 1. 
      CHARACTER(7) :: CMARK 
      REAL(8)      :: WORK_N(IDX)
      REAL(8)      :: WORK_S(IDX)

#endif
!
! 20050517
      REAL(8)      :: WTR_N(IDX)
      REAL(8)      :: WTR_S(IDX)
      REAL(8)      :: GLA_N(IDX)
      REAL(8)      :: GLA_S(IDX)
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
        WRITE(6,*) 'TOTAL_H2O_ENG_OF_SKIN 02/02 M.Hosaka'
        LFIRST = .FALSE. 
      ENDIF
!
! ================
! >> モニタ準備 <<
! ================
!
      DO I = 1,NMAX 
        XLAI_VEG_ms  (I,1)    = XLAI_mo (NTYP(I),MON_NEXT,1) 
        XLAI_VEG_ms  (I,2)    = XLAI_mo (NTYP(I),MON_NEXT,2) 
!
!#ifdef CORRECT 
        XZ     =   Z2_mo  ( NTYP(I) , MON_NEXT ) 
     1           - Z1_mo  ( NTYP(I) , MON_NEXT )  
        SDEP   = MIN ( DEP_SNOW_NEXT(I,5) 
     1                   - Z1_mo  ( NTYP(I) , MON_NEXT )  
     1                  , XZ * 0.95D0 )
!#else
!        XZ     = Z2_mo  ( NTYP(I) , MON_NEXT )
!        SDEP   = MIN ( DEP_SNOW_NEXT(I,5)  , XZ * 0.95D0 )  
!#endif
        SDEP   = MAX ( SDEP , ZERO ) 
        FAC    = 1.0D0 - SDEP / XZ            
!
        XLAI_VEG_msx (I,1)    = XLAI_mo (NTYP(I),MON_NEXT,1) * FAC
        XLAI_VEG_msx (I,2)    = 0. 
      ENDDO
!
! ================
! >> キャノピー <<
! ================
!
      DO I=1,NMAX
        WTR_N(I) = FR_WTR_CNP_NOS_NEXT(I)
     1             * XLAI_VEG_ms (I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
        WTR_S(I) = FR_WTR_CNP_SNW_NEXT(I)
     1             * XLAI_VEG_msx(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
        GLA_N(I) = FR_GLA_CNP_NOS_NEXT(I)
     1             * XLAI_VEG_ms (I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
        GLA_S(I) = FR_GLA_CNP_SNW_NEXT(I)
     1             * XLAI_VEG_msx(I,1) * H2O_SAT_CNP_DEP_c * RHO_WTR_c 
!
        H2O_CNP_NOS(I) = ( WTR_N(I)+GLA_N(I) ) * ( 1-CVR_SNW_NEXT(I) )     
        H2O_CNP_SNW(I) = ( WTR_S(I)+GLA_S(I) ) * CVR_SNW_NEXT(I)      
!
        ENG_CNP_NOS(I) =
     1             (   CAP_LAI_c * XLAI_VEG_ms(I,1) 
     1               + CAP_WTR_c * WTR_N(I)
     1               + CAP_GLA_c * GLA_N(I) ) 
     1             * ( TMP_CNP_NOS_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_N(I) 
        ENG_CNP_NOS(I) = ENG_CNP_NOS(I) * ( 1-CVR_SNW_NEXT(I) )     
!
        ENG_CNP_SNW(I) =
     1             (   CAP_LAI_c * XLAI_VEG_msx(I,1) 
     1               + CAP_WTR_c * WTR_S(I)
     1               + CAP_GLA_c * GLA_S(I) ) 
     1             * ( TMP_CNP_SNW_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_S(I) 
        ENG_CNP_SNW(I) = ENG_CNP_SNW(I) * CVR_SNW_NEXT(I)      
!
      ENDDO 
!
! ================
! >> 下草・地面 <<
! ================
!
      DO I=1,NMAX
        WTR_N(I) = FR_WTR_GRS_NEXT(I)
     1             * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
        GLA_N(I) = FR_GLA_GRS_NEXT(I)
     1             * XLAI_VEG_ms(I,2) * H2O_SAT_GRS_DEP_c * RHO_WTR_c 
!
        H2O_GRSK   (I) = ( WTR_N(I)+GLA_N(I) ) * ( 1-CVR_SNW_NEXT(I) )     
!
        ENG_GRSK   (I) =
     1             (   CAP_GRSK_c      
     1               + CAP_LAI_c * XLAI_VEG_ms(I,2) 
     1               + CAP_WTR_c * WTR_N(I)
     1               + CAP_GLA_c * GLA_N(I) ) 
     1             * ( TMP_GRSK_NEXT(I) - TMP_FREZ_c )
     1             - HL_W2G_c * GLA_N(I) 
        ENG_GRSK   (I) = ENG_GRSK   (I) * ( 1-CVR_SNW_NEXT(I) )     
!
        ENG_SNSK   (I) =
     1                CAP_SNSK_c 
     1             * ( TMP_SNSK_NEXT(I) - TMP_FREZ_c )
     1             * CVR_SNW_NEXT(I) 
      ENDDO 
!
      RETURN
      END
!=========================================================
      SUBROUTINE MONIT_TOTAL_H2O_ENG_OF_SKIN (
     I    NMAX                , 
     I    H2O_CNP_NOS_1       , H2O_CNP_SNW_1       ,
     I    H2O_GRSK_1      , 
     I    ENG_CNP_NOS_1       , ENG_CNP_SNW_1       ,
     I    ENG_GRSK_1          , ENG_SNSK_1          ,
     I    H2O_CNP_NOS_2       , H2O_CNP_SNW_2       ,
     I    H2O_GRSK_2      , 
     I    ENG_CNP_NOS_2       , ENG_CNP_SNW_2       ,
     I    ENG_GRSK_2          , ENG_SNSK_2          )  
!
      USE SIBPRM , ONLY : 
     1    IDX  
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL  , INA ,
     1                      MONIT_SIB_NEXT
#endif
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN)   :: H2O_CNP_NOS_1         (IDX) 
      REAL(8),INTENT(IN)   :: H2O_CNP_SNW_1         (IDX) 
      REAL(8),INTENT(IN)   :: H2O_GRSK_1            (IDX) 
      REAL(8),INTENT(IN)   :: ENG_CNP_NOS_1         (IDX) 
      REAL(8),INTENT(IN)   :: ENG_CNP_SNW_1         (IDX) 
      REAL(8),INTENT(IN)   :: ENG_GRSK_1            (IDX) 
      REAL(8),INTENT(IN)   :: ENG_SNSK_1            (IDX) 
!
      REAL(8),INTENT(IN)   :: H2O_CNP_NOS_2         (IDX) 
      REAL(8),INTENT(IN)   :: H2O_CNP_SNW_2         (IDX) 
      REAL(8),INTENT(IN)   :: H2O_GRSK_2            (IDX) 
      REAL(8),INTENT(IN)   :: ENG_CNP_NOS_2         (IDX) 
      REAL(8),INTENT(IN)   :: ENG_CNP_SNW_2         (IDX) 
      REAL(8),INTENT(IN)   :: ENG_GRSK_2            (IDX) 
      REAL(8),INTENT(IN)   :: ENG_SNSK_2            (IDX) 
!
      INTEGER,INTENT(IN)   :: NMAX
!
      REAL(8)  :: WORK(IDX)
!      
      INTEGER :: I      
!
      REAL(8),PARAMETER :: ONE  = 1. 
!
      CHARACTER( 7)  :: CMARK 
!
!
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! キャノピー水
!
      DO I=1,NMAX
        WORK(I) =    ( H2O_CNP_NOS_2 (I) + H2O_CNP_SNW_2 (I) )
     1            -  ( H2O_CNP_NOS_1 (I) + H2O_CNP_SNW_1 (I) )
      ENDDO 
      CMARK  = 'H2OCDIF'
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,INA , ONE , ONE )      
!
! キャノピーエネルギー
!
      DO I=1,NMAX
        WORK(I) =    ( ENG_CNP_NOS_2 (I) + ENG_CNP_SNW_2 (I) )
     1            -  ( ENG_CNP_NOS_1 (I) + ENG_CNP_SNW_1 (I) )
      ENDDO 
      CMARK  = 'ENGCDIF'
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,INA , ONE , ONE )      
!
! スキン水
!
      DO I=1,NMAX
        WORK(I) =    H2O_GRSK_2 (I) - H2O_GRSK_1 (I) 
      ENDDO 
      CMARK  = 'H2OGDIF'
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,INA , ONE , ONE )      
!
! スキンエネルギー
!
      DO I=1,NMAX
        WORK(I) =    ( ENG_GRSK_2 (I) + ENG_SNSK_2 (I) )
     1            -  ( ENG_GRSK_1 (I) + ENG_SNSK_1 (I) )
      ENDDO 
      CMARK  = 'ENGGDIF'
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,INA , ONE , ONE )      
!
      ENDIF
#endif
!
      RETURN
      END SUBROUTINE MONIT_TOTAL_H2O_ENG_OF_SKIN 
!=========================================================
      SUBROUTINE MR_MONIT_TOTAL_H2O_ENG_OF_SKIN
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
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MR_MONIT_TOTAL_H2O_ENG_OF_SKIN 02/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN
        CMARK  = 'H2OCDIF'
        CTITLE = 'H2O CNP INCRESE'
        CUNIT  = 'kg/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'ENGCDIF'
        CTITLE = 'ENG CNP INCRESE'
        CUNIT  = 'J/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'H2OGDIF'
        CTITLE = 'H2O SKIN INCRESE'
        CUNIT  = 'kg/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'ENGGDIF'
        CTITLE = 'ENG SKIN INCRESE'
        CUNIT  = 'J/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END


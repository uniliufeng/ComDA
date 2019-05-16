! —¤–Ê‰ß’öŒvŽZ‚Ì‚½‚ß‚Ì”€”õ
! 99/07/xx •ÛâªG
! 99/08/04 •ÛâªG ‘åò‚³‚ñ‚É“n‚·ƒo[ƒWƒ‡ƒ“ƒtƒBƒbƒNƒX
!
      SUBROUTINE PRESIB_MAIN  ( 
     I    NMAX_01     , NMAX_1      , NMAX_2      , NMAX_3      ,
     I    NMAX_4      , NMAX        ,
#ifdef MONYOS
     I    DELT_SOIL   ,
#endif
     I    NTYP        , 
     I    CVR_SNW_INV , 
     I    DEP_SNOW    , RHO_SNOW_INV   ,
     I    WTR_SNOW    , GLA_SNOW    , 
     I    TMP_SOIL    , FR_WTR_SOIL , FR_GLA_SOIL , 
     I    TMP_SNSL    , FR_WTR_SNSL , FR_GLA_SNSL ,
     O    HDF_GRSK    , HDF_SNSK    ,
     O    HDF_SOIL_H  , HDF_SNSL_H  ,
     O    HDF_SNOW_H  , 
     O    CAP_SOIL_TOP, CAP_SNOW_TOP, 
     O    ENG_SOIL_TOP   )
!
      USE SIBPRM , ONLY :
     1    IDX , 
     1    ISN , 
     1    IDP , 
     1    ICG 
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NMAX_01 (3)
      INTEGER,INTENT(IN)    :: NMAX_1  (3)
      INTEGER,INTENT(IN)    :: NMAX_2  (3)
      INTEGER,INTENT(IN)    :: NMAX_3  (3)
      INTEGER,INTENT(IN)    :: NMAX_4  (3)
      INTEGER,INTENT(IN)    :: NMAX
!
      INTEGER,INTENT(IN)    :: NTYP(IDX)
!
      REAL(8),INTENT(IN)    :: DEP_SNOW  (IDX,ISN+1)    ! Ïá‘w‚ÌŒú‚³(m)
      REAL(8),INTENT(IN)    :: RHO_SNOW_INV (IDX,ISN)    ! Ïá‘w‚Ì–§“x(kg/m3)
      REAL(8),INTENT(IN)    :: WTR_SNOW  (IDX)          ! 
      REAL(8),INTENT(IN)    :: GLA_SNOW  (IDX)          ! 
      REAL(8),INTENT(IN)    :: CVR_SNW_INV  (IDX)          ! 
!
#ifdef MONYOS
      REAL(8),INTENT(IN)    :: DELT_SOIL 
#endif
!
      REAL(8),INTENT(IN)    :: TMP_SOIL   (IDX,IDP)  
      REAL(8),INTENT(IN)    :: TMP_SNSL   (IDX,IDP)  
!
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL(IDX,IDP)     ! ‹óŒ„’†‚Ì…‚ÌŠ„‡
!                                                        (‘ÌÏ”ä)
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL(IDX,IDP)     ! ‹óŒ„’†‚Ì•X‚ÌŠ„‡
!                                                        (‘ÌÏ”ä)
      REAL(8),INTENT(IN)    :: FR_WTR_SNSL(IDX,IDP)     ! ‹óŒ„’†‚Ì…‚ÌŠ„‡
!                                                        (Ïá‰º“yë‚P‘w)
      REAL(8),INTENT(IN)    :: FR_GLA_SNSL(IDX,IDP)     ! ‹óŒ„’†‚Ì•X‚ÌŠ„‡
!                                                        (Ïá‰º“yë‚P‘w)!
! OUTPUT
!
! ”M“`“±ŒW”
!
      REAL(8),INTENT(OUT)   :: HDF_GRSK   (IDX,2)       ! “y   SK-L1,L1-L2
      REAL(8),INTENT(OUT)   :: HDF_SNSK   (IDX,2)       ! á‰º SK-S1,S1-X 
      REAL(8),INTENT(OUT)   :: HDF_SNOW_H (IDX,ISN+1)   
!                                     ! á   S1-S2,S2-S3,S3-S4,SX-L1,L1-L2
      REAL(8),INTENT(OUT)   :: HDF_SOIL_H (IDX,IDP)     ! “y 1-2,2-3,3-
      REAL(8),INTENT(OUT)   :: HDF_SNSL_H (IDX,IDP)     ! á‰º“y   1-2
!
! ”ä”M
!
      REAL(8),INTENT(OUT)   :: CAP_SOIL_TOP (IDX) 
      REAL(8),INTENT(OUT)   :: CAP_SNOW_TOP (IDX) 
!
! “yë‘æˆê‘wƒGƒlƒ‹ƒM[ ( á–³‚µ‚È‚çá–³‚µAá—L‚è‚È‚çá—L‚èˆæ‚¾‚¯‚Å‚Ì•½‹Ï )
!
      REAL(8),INTENT(OUT)   :: ENG_SOIL_TOP (IDX) 
!
      LOGICAL, SAVE :: LFIRST 
      DATA LFIRST  / .TRUE. / 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'PRESIB_MAIN 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF 
!
! --------------------                                                   
! > ”M“`“±ŒW”‚ÌŒvŽZ <
! --------------------                                                   
!
        CALL PRESIB_CONDUC_ALL (
     I    NMAX_01     , NMAX_1      , NMAX_2      , NMAX_3      ,
     I    NMAX_4      , NMAX        , 
#ifdef MONYOS
     I    DELT_SOIL   ,
#endif
     I    NTYP        , 
     I    DEP_SNOW    , RHO_SNOW_INV,
     I    FR_WTR_SOIL , FR_GLA_SOIL , FR_WTR_SNSL , FR_GLA_SNSL ,
     O    HDF_GRSK    , HDF_SNSK    ,
     O    HDF_SOIL_H  , HDF_SNSL_H  ,
     O    HDF_SNOW_H  )
!
! -------------------------
! > “yë/Ïá‘æˆê‘w‚Ì”ä”M <
! -------------------------
!
      CALL PRESIB_CAP_SSTOP (
     I   1               , NMAX_01(3)    ,
     I   NMAX_01(2)      , NMAX          , 
     I   NTYP            , 
     I   CVR_SNW_INV     , 
     I   TMP_SOIL        , FR_WTR_SOIL     , FR_GLA_SOIL   , 
     I   TMP_SNSL        , FR_WTR_SNSL     , FR_GLA_SNSL   , 
     I   WTR_SNOW        , GLA_SNOW      , 
     O   CAP_SOIL_TOP    , CAP_SNOW_TOP  , 
     O   ENG_SOIL_TOP )  
!
      RETURN
      END

! “yëEÏáŠe‘w‚Ì”M“`“±ŒW”(=”M“`“±—¦/Œú‚³) `áE“yëŒvŽZŽž‚Ì‚ÝŒvŽZ‚·‚é‚à‚Ì
! 99/06/1x •ÛâªG
! 99/06/2x ‘åò      class “™‚ÅÌ—p‚µ‚Ä‚¢‚éŒ`Ž®‚É•ÏXB
! 99/07/13 •ÛâªG  “yë’†‚Å‚Í•X‚Æ…‚Ì–§“x‚ª‚¨‚È‚¶‚Æ‚·‚é‚±‚Æ‚ÉŠÖ‚·‚éC³
!                    soilcon.F ‚ð‚­‚Á‚Â‚¯‚é
! 99/08/04 •ÛâªG ‘åò‚³‚ñ‚É“n‚·ƒo[ƒWƒ‡ƒ“ƒtƒBƒbƒNƒX
!
      SUBROUTINE PRESIB_CONDUC_ALL (
     I    NMAX_01     , NMAX_1      , NMAX_2      , NMAX_3      ,
     I    NMAX_4      , NMAX        ,
#ifdef MONYOS
     I    DELT_SOIL   ,
#endif
     I    NTYP        , 
     I    DEP_SNOW    , RHO_SNOW_INV,
     I    FR_WTR_SOIL , FR_GLA_SOIL , FR_WTR_SNSL  , FR_GLA_SNSL ,
     O    HDF_GRSK    , HDF_SNSK    ,
     O    HDF_SOIL_H  , HDF_SNSL_H  ,
     O    HDF_SNOW_H  )
!
      USE SIBPRM , ONLY :
     1      IDX ,
     1      IDP ,
     1      ISN
!
      USE SIBCON , ONLY :
     1      HDF_SOIL_c  ,
     1      HDF_WTR_c   ,
     1      HDF_GLA_c   ,
     1      FACT_HDF_GRSK_c ,   
     1      FACT_HDF_SOIL_c ,   
     1      FACT_HDF_SNSK_c ,  
     1      FACT_HDF_SNOW_c ,
     1      CAP_SOIL4_c     ,
     1      DEP_SOIL4_c    
!
      USE SIBVEG , ONLY :
     1           DEP_SOIL_co          , ! “yë‘w‚Ì[‚³    (ITYP,IDP)
     1           DEP_SOIL_INV_co        ! “yë‘w‚Ì[‚³    (ITYP,IDP)
!
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN, ICS,INA,INS,
     1                      IMONIT_LEVEL , IMONIT_ALL  
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: NMAX_01 (3)
      INTEGER,INTENT(IN)    :: NMAX_1  (3)
      INTEGER,INTENT(IN)    :: NMAX_2  (3)
      INTEGER,INTENT(IN)    :: NMAX_3  (3)
      INTEGER,INTENT(IN)    :: NMAX_4  (3)
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NTYP(IDX)
!
#ifdef MONYOS
      REAL(8),INTENT(IN)    :: DELT_SOIL ! ŽŸ‚Ì“yëŒvŽZ‚Ü‚Å‚ÌŽžŠÔ‚Ý
#endif

!
      REAL(8),INTENT(IN)    :: DEP_SNOW  (IDX,ISN)       ! Ïá‘w‚ÌŒú‚³(m)
      REAL(8),INTENT(IN)    :: RHO_SNOW_INV (IDX,ISN)    ! Ïá‘w‚Ì–§“x(kg/m3)
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
      REAL(8),INTENT(OUT)   :: HDF_GRSK    (IDX,2) 
      REAL(8),INTENT(OUT)   :: HDF_SNSK    (IDX,2) 
      REAL(8),INTENT(OUT)   :: HDF_SOIL_H  (IDX,IDP)     ! “y 1-2,2-3,3-
      REAL(8),INTENT(OUT)   :: HDF_SNSL_H  (IDX,IDP)     ! á‰º“y   1-2
      REAL(8),INTENT(OUT)   :: HDF_SNOW_H  (IDX,ISN+1)   
!                                   ! S1-S2,S2-S3,S3-S4,SX-L1,L1-L2 
!
! WORK
!
      REAL(8)      :: RAM_SOIL_INV  (IDX,IDP)    ! “yë”M“`“±—¦(‹t”)
      REAL(8)      :: RAM_SNSL_INV  (IDX,IDP)    ! “yë”M“`“±—¦(‹t”)
      REAL(8)      :: RAM_SOIL      (IDX,IDP)    ! 
      REAL(8)      :: RAM_SNSL      (IDX,IDP)    ! 
!
      REAL(8)      :: RAM_SNOW_INV  (IDX,ISN)    ! áŠe‘w‚ÌƒÉ^-1
      REAL(8)      :: RAM_SNOW_B_INV(IDX)        ! á‰º’[‚ÌƒÉ^-1
      REAL(8)      :: RAM_SNOW      (IDX,ISN)    ! Ïá‘æˆê‘w
!
      INTEGER      :: I
!
!!    REAL(8),SAVE          :: HDF_SKIN_MAX = 100.D0 
      REAL(8),SAVE          :: HDF_SKIN_MAX = 1.D10
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
#endif
      REAL(8),PARAMETER :: ONE = 1.D0 

!
      LOGICAL, SAVE :: LFIRST 
      DATA LFIRST  / .TRUE. / 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'PRESIB_CONDUC_ALL 99/07 POCHAKA'
        write(6,*) '  FACT_HDF_SNSK_c ' , FACT_HDF_SNSK_c 
        write(6,*) '  FACT_HDF_SNOW_c ' , FACT_HDF_SNOW_c 
        LFIRST = .FALSE.
      ENDIF 
!
! ==============================
! >> á‚ÌŠe‘w‚Ì”M“`“±—¦‚Ì‹t” <<
! ==============================
!
!
      ram_snow    (:,:) = 0.D0 
      ram_snow_inv(:,:) = 0.D0 
!
!  check check check !!! á‚Ì”M“`“±—¦‚ÌŠÖ”Œ`
!
      DO 100 I=NMAX_01(2),NMAX
        RAM_SNOW_INV(I,1) 
     1      = 10.D0**(-2.650D0 * 1.0D-3/RHO_SNOW_INV(I,1) + 1.652D0)
        RAM_SNOW      (I,1) = 1 / RAM_SNOW_INV(I,1) 
        RAM_SNOW_B_INV(I)   = RAM_SNOW_INV(I,1)
  100 CONTINUE

      DO 120 I=NMAX_2(2),NMAX
        RAM_SNOW_INV(I,2) 
     1      = 10.D0**(-2.650D0*1.0D-3/RHO_SNOW_INV(I,2) + 1.652D0)
#ifdef MONYOS
        RAM_SNOW      (I,2) = 1 / RAM_SNOW_INV(I,2)
#endif
        RAM_SNOW_B_INV(I) = RAM_SNOW_INV(I,2)
  120 CONTINUE
!
      DO 130 I=NMAX_3(2),NMAX
        RAM_SNOW_INV(I,3) 
     1      = 10.D0**( - 2.650D0*1.0D-3/RHO_SNOW_INV(I,3) + 1.652D0)
#ifdef MONYOS
        RAM_SNOW      (I,3) = 1 / RAM_SNOW_INV(I,3)
#endif
        RAM_SNOW_B_INV(I) = RAM_SNOW_INV(I,3)
  130 CONTINUE
!
      DO 140 I=NMAX_4(2),NMAX
        RAM_SNOW_INV(I,4) 
     1      = 10.D0**( - 2.650D0*1.0D-3/RHO_SNOW_INV(I,4) + 1.652D0)
#ifdef MONYOS
        RAM_SNOW      (I,4) = 1 / RAM_SNOW_INV(I,4)
#endif
        RAM_SNOW_B_INV(I) = RAM_SNOW_INV(I,4)
  140 CONTINUE
!
! ==============================
! >> “yëŠe‘w‚Ì”M“`“±—¦‚Ì‹t” <<
! ==============================
!
!  á‚È‚µˆæ‚Ì“yë
!
      CALL PRESIB_CONDUC_SOIL ( 
     I  1            , NMAX_01(3)  , NTYP    , 
     I  FR_WTR_SOIL  , FR_GLA_SOIL , 
     O  RAM_SOIL_INV , RAM_SOIL )
!
!  á‚ ‚èˆæ‚Ì“yë
!
      CALL PRESIB_CONDUC_SOIL ( 
     I  NMAX_01(2)   , NMAX , NTYP , 
     I  FR_WTR_SNSL  , FR_GLA_SNSL , 
     O  RAM_SNSL_INV , RAM_SNSL )
!
! =======================
! >>  ”M“`“±ŒW”‚ÌŒvŽZ <<
! =======================
!
! --------------
! > á‚È‚µ“yë < 
! --------------
!
      DO I=1,NMAX_01(3)
        HDF_SOIL_H(I,1)
     1     = 1 / (  0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SOIL_INV(I,1)
     1            + 0.5D0 * DEP_SOIL_co(NTYP(I),2) * RAM_SOIL_INV(I,2) )
     1         * FACT_HDF_SOIL_c 
!
        HDF_SOIL_H(I,2)
     1    = 1 / (   0.5D0 * DEP_SOIL_co(NTYP(I),2) * RAM_SOIL_INV(I,2)
     1            + 0.5D0 * DEP_SOIL_co(NTYP(I),3) * RAM_SOIL_INV(I,3))
!    1         * FACT_HDF_SOIL_c 
!
        HDF_SOIL_H(I,3)
     1    = 1 / (   0.5D0 * DEP_SOIL_co(NTYP(I),3) * RAM_SOIL_INV(I,3)
     1            + 0.5D0 * DEP_SOIL4_c            * RAM_SOIL_INV(I,3))
!
      ENDDO
!
! --------------
! > á‚ ‚è“yë <
! --------------
!
      DO I=NMAX_01(2),NMAX
        HDF_SNSL_H(I,1)
     1   = 1. / (  0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SNSL_INV(I,1)
     1           + 0.5D0 * DEP_SOIL_co(NTYP(I),2) * RAM_SNSL_INV(I,2) )
!
        HDF_SNSL_H(I,2)
     1    = 1 / (   0.5D0 * DEP_SOIL_co(NTYP(I),2) * RAM_SNSL_INV(I,2)
     1            + 0.5D0 * DEP_SOIL_co(NTYP(I),3) * RAM_SNSL_INV(I,3))
!
!        HDF_SNSL_H(I,3)
!     1     =  DEP_SOIL_INV_co(NTYP(I),3) * RAM_SNSL(I,3) 
!
        HDF_SNSL_H(I,3)
     1    = 1 / (   0.5D0 * DEP_SOIL_co(NTYP(I),3) * RAM_SNSL_INV(I,3)
     1            + 0.5D0 * DEP_SOIL4_c            * RAM_SNSL_INV(I,3))


      ENDDO
!
! ------------
! > áˆæŒvŽZ < 
! ------------
!
!  S1-L1 ‘w
!
      HDF_SNOW_H(:,:) = 0.
!
      DO I=NMAX_01(2),NMAX_1(3)
        HDF_SNOW_H (I,1)    
     1   = 1/ (   0.5D0 * DEP_SNOW(I,1)          * RAM_SNOW_B_INV (I)
     1          + 0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SNSL_INV   (I,1))      
     1      * FACT_HDF_SNOW_c 
        HDF_SNOW_H (I,ISN) = HDF_SNOW_H (I,1)
      ENDDO 
!
!  S1-S2 ‘w
!
      DO I=NMAX_2(2),NMAX 
        HDF_SNOW_H (I,1)    
     1    = MIN (
     1         1 / (   0.5D0 * DEP_SNOW(I,1) * RAM_SNOW_INV (I,1)
     1               + 0.5D0 * DEP_SNOW(I,2) * RAM_SNOW_INV (I,2) )
     1      , HDF_SKIN_MAX ) 
     1      * FACT_HDF_SNOW_c 
      ENDDO
!
!  S2-L1 ‘w
!
      DO I=NMAX_2(2),NMAX_2(3)
        HDF_SNOW_H (I,2)    
     1   = 1/ (   0.5D0 * DEP_SNOW(I,2)          * RAM_SNOW_B_INV (I)
     1          + 0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SNSL_INV   (I,1))
!    1      * FACT_HDF_SNOW_c 
        HDF_SNOW_H (I,ISN) = HDF_SNOW_H (I,2)     
      ENDDO
!
!  S2-S3 ‘w
!
      DO I=NMAX_3(2),NMAX 
        HDF_SNOW_H (I,2)    
     1    = 1 / (   0.5D0 * DEP_SNOW(I,2) * RAM_SNOW_INV (I,2)
     1            + 0.5D0 * DEP_SNOW(I,3) * RAM_SNOW_INV (I,3) )
!    1      * FACT_HDF_SNOW_c 
      ENDDO
!
!  S3-L1 ‘w
!
      DO I=NMAX_3(2),NMAX_3(3) 
        HDF_SNOW_H (I,3)    
     1   = 1/ (   0.5D0 * DEP_SNOW(I,3)          * RAM_SNOW_B_INV (I)
     1          + 0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SNSL_INV   (I,1))      
!    1      * FACT_HDF_SNOW_c 
        HDF_SNOW_H (I,ISN) = HDF_SNOW_H (I,3)     
      ENDDO
!
      DO I=NMAX_4(2),NMAX 
!  S3-S4 ‘w
        HDF_SNOW_H (I,3)
     1    = 1 / (   0.5D0 * DEP_SNOW(I,3) * RAM_SNOW_INV (I,3)
     1            + 0.5D0 * DEP_SNOW(I,4) * RAM_SNOW_INV (I,4) )
!    1      * FACT_HDF_SNOW_c 
!  S4-L1 ‘w
        HDF_SNOW_H (I,4)  
     1   = 1/ (   0.5D0 * DEP_SNOW(I,4)          * RAM_SNOW_B_INV (I)
     1          + 0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SNSL_INV   (I,1))      
!    1      * FACT_HDF_SNOW_c 
      ENDDO
!
!  L1-L2 
!
      DO I=NMAX_01(2),NMAX 
        HDF_SNOW_H (I,ISN+1)  
     1     = HDF_SNSL_H (I,1)  
     1      * FACT_HDF_SOIL_c 
      ENDDO
!
! ------------------
! > surface ˆæŒvŽZ < 
! ------------------
!
      DO I=1,NMAX_01(3)
! GRSK-L1
        HDF_GRSK (I,1) = 2 * RAM_SOIL(I,1) * DEP_SOIL_INV_co(NTYP(I),1)  
     1                     * FACT_HDF_GRSK_c 
! L1-L2 
        HDF_GRSK (I,2) 
     1     = 1 / (  0.5D0 * DEP_SOIL_co(NTYP(I),1) * RAM_SOIL_INV(I,1)
     1            + 0.5D0 * DEP_SOIL_co(NTYP(I),2) * RAM_SOIL_INV(I,2) )
     1                     * FACT_HDF_SOIL_c 
      ENDDO
!
! SNSK-S1
      DO I=NMAX_01(2),NMAX
        HDF_SNSK  (I,1) = MIN ( 2 * RAM_SNOW(I,1) / DEP_SNOW(I,1) ,
     1                        HDF_SKIN_MAX )  
     1      * FACT_HDF_SNSK_c 
      ENDDO 
!
! S1-S2
      DO I=NMAX_01(2),NMAX
        HDF_SNSK  (I,2)  =   HDF_SNOW_H (I,1)    
     1                      * FACT_HDF_SNOW_c 
      ENDDO  ! áˆê‘w‚È‚ç S1-L1, 2 ‘wˆÈã‚È‚ç S1-S2 
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
        CMARK = 'KDFSNSK' 
        CALL MONIT_SIB ( CMARK, HDF_SNSK(1,1)   ,INS , DELT_SOIL, ONE ) 
        CMARK = 'KDFSN1' 
        CALL MONIT_SIB ( CMARK, HDF_SNOW_H(1,1) ,INS , DELT_SOIL, ONE ) 
        CMARK = 'KDFSN2' 
        CALL MONIT_SIB ( CMARK, HDF_SNOW_H(1,2) ,INS , DELT_SOIL, ONE ) 
        CMARK = 'KDFSN3' 
        CALL MONIT_SIB ( CMARK, HDF_SNOW_H(1,3) ,INS , DELT_SOIL, ONE ) 
!
        CMARK = 'KDFGRSK' 
        CALL MONIT_SIB ( CMARK, HDF_GRSK(1,1) ,ICN , DELT_SOIL, ONE ) 
!
        CMARK = 'KDFSL1' 
        CALL MONIT_SIB ( CMARK, HDF_SOIL_H(1,1) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, HDF_SNSL_H(1,1) ,ICS , DELT_SOIL, ONE ) 
        CMARK = 'KDFSL2' 
        CALL MONIT_SIB ( CMARK, HDF_SOIL_H(1,2) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, HDF_SNSL_H(1,2) ,ICS , DELT_SOIL, ONE ) 
        CMARK = 'KDFSL3' 
        CALL MONIT_SIB ( CMARK, HDF_SOIL_H(1,3) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, HDF_SNSL_H(1,3) ,ICS , DELT_SOIL, ONE ) 
!
        CMARK = 'RAMSNW1' 
        CALL MONIT_SIB ( CMARK, RAM_SNOW(1,1) ,INA , DELT_SOIL, ONE ) 
        CMARK = 'RAMSNW2' 
        CALL MONIT_SIB ( CMARK, RAM_SNOW(1,2) ,INA , DELT_SOIL, ONE ) 
        CMARK = 'RAMSNW3' 
        CALL MONIT_SIB ( CMARK, RAM_SNOW(1,3) ,INA , DELT_SOIL, ONE ) 
        CMARK = 'RAMSNW4' 
        CALL MONIT_SIB ( CMARK, RAM_SNOW(1,4) ,INA , DELT_SOIL, ONE ) 
!
        CMARK = 'RAMSL1' 
        CALL MONIT_SIB ( CMARK, RAM_SOIL(1,1) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, RAM_SNSL(1,1) ,ICS , DELT_SOIL, ONE ) 
        CMARK = 'RAMSL2' 
        CALL MONIT_SIB ( CMARK, RAM_SOIL(1,2) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, RAM_SNSL(1,2) ,ICS , DELT_SOIL, ONE ) 
        CMARK = 'RAMSL3' 
        CALL MONIT_SIB ( CMARK, RAM_SOIL(1,3) ,ICN , DELT_SOIL, ONE ) 
        CALL MONIT_SIB ( CMARK, RAM_SNSL(1,3) ,ICS , DELT_SOIL, ONE ) 
      ENDIF
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_PRESIB_CONDUC_ALL
!
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL  
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
        WRITE(6,*) 'MONIT_REGIST_PRESIB_CONDUC_ALL 02/03 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'KDFSNSK'
      CTITLE = 'RAM SNOW SKIN-SNOW1'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSN1'
      CTITLE = 'RAM SNOW SNOW1-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSN2'
      CTITLE = 'RAM SNOW SNOW2-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSN3'
      CTITLE = 'RAM SNOW SNOW3-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFGRSK'
      CTITLE = 'RAM SOIL SKIN-SOIL1'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSL1'
      CTITLE = 'RAM SNOW SOIL1-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSL2'
      CTITLE = 'RAM SNOW SOIL2-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'KDFSL3'
      CTITLE = 'RAM SNOW SOIL3-'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ----------------------
!
      CMARK  = 'RAMSNW1'
      CTITLE = 'RAM SNOW 1'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSNW2'
      CTITLE = 'RAM SNOW 2'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSNW3'
      CTITLE = 'RAM SNOW 3'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSNW4'
      CTITLE = 'RAM SNOW 4'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSL1'
      CTITLE = 'RAM SOIL 1'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSL2'
      CTITLE = 'RAM SOIL 2'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'RAMSL3'
      CTITLE = 'RAM SOIL 3'
      CUNIT  = 'W/M/K'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      RETURN
      END
!==================================================
      SUBROUTINE PRESIB_CONDUC_SOIL(
     I  NSTR         , NEND        , NTYP     ,  
     I  FR_WTR_SOIL  , FR_GLA_SOIL , 
     O  RAM_SOIL_INV , RAM_SOIL     )
!
! “yë”M“`“±—¦FCox et.al(1999),Climate Dynamics
! ‚½‚¾‚µA–{ƒ‚ƒfƒ‹‚Å‚ÍA“yë’†‚Ì•X‚Æ…‚Í‚¨‚È‚¶–§“x‚ðŽ‚Â‚Ì‚ÅA
! ‘ÌÏ”ä•â³‚Í‚µ‚È‚¢B
!
      USE SIBPRM , ONLY : 
     1      IDX  ,
     1      IDP
!
      USE SIBCON , ONLY :
     1      HDF_SOIL_c  ,  
     1      HDF_WTR_c   ,  
     1      HDF_GLA_c   ,
     1      HDF_AIR_c   
!
      USE SIBVEG , ONLY : 
     1       POR_SOIL_co              ! ‹óŒ„—¦   (ITYP)  
!
      IMPLICIT NONE      

! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: NTYP       (IDX)
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL(IDX,IDP)     ! ‹óŒ„’†‚Ì…‚ÌŠ„‡
!                                                         (‘ÌÏ”ä)
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL(IDX,IDP)     ! ‹óŒ„’†‚Ì•X‚ÌŠ„‡
!                                                        (‘ÌÏ”ä)
!    
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: RAM_SOIL_INV (IDX,IDP)   ! “y‚ÌƒÉ^-
      REAL(8),INTENT(OUT)   :: RAM_SOIL     (IDX,IDP) 
!
! WORK 
!
      REAL(8)               :: RAM_DRY              ! “yë‚ÌŠ£‘‡”M“`“±—¦
      REAL(8)               :: RAM_SAT              ! –O˜a“yë”M“`“±—¦
!
      REAL(8)               :: FR_G      ! “yë•X‚ð—Z‚©‚µ‚½Žž‚Ì–O˜a“x
      REAL(8)               :: FR_W      ! “yë…‚Ì–O˜a“xi‹ŒSiB‚ÌWj
      REAL(8)               :: THETA_S   ! ‹óŒ„—¦(=–O˜aŽž‚Ì‘ÌÏŠÜ…—¦)
      REAL(8)               :: THETA_W
      REAL(8)               :: THETA_G
!
      INTEGER               :: I
      INTEGER               :: K
!
      real(8),save :: fr_min = 1.D-30
!
      LOGICAL, SAVE :: LFIRST 
      DATA LFIRST  / .TRUE. / 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'PRESIB_CONDUC_SOIL 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF 
!
! ”M“`“±—¦‚Ì‹t”ŒvŽZ
!
      DO K=1,IDP
      DO I=NSTR,NEND
!
        FR_W     = FR_WTR_SOIL(I,K)               ! ‹óŒ„’†‚Ì…‚Ìè‚ß‚é‘ÌÏ”ä
        FR_G     = FR_GLA_SOIL(I,K)               ! ‹óŒ„’†‚Ì•X‚Ìè‚ß‚é‘ÌÏ”ä
        THETA_S  = POR_SOIL_co(NTYP(I))           ! –O˜a‚µ‚½ê‡‚Ì…‚Ì‘ÌÏ”ä
        THETA_W  = THETA_S * FR_W / MAX(FR_W + FR_G,FR_MIN)
        THETA_G  = THETA_S - THETA_W 
!
        RAM_DRY  =    HDF_SOIL_c ** (1.0D0-THETA_S)
     1              * HDF_AIR_c  ** THETA_S
        RAM_SAT  =    HDF_SOIL_c ** (1.0D0-THETA_S)
     1              * HDF_WTR_c  ** THETA_W
     1              * HDF_GLA_c  ** THETA_G
!
        RAM_SOIL(I,K) = (   RAM_SAT * ( FR_W + FR_G ) 
     1                    + RAM_DRY * ( 1 - FR_W - FR_G ) )
!
        RAM_SOIL_INV(I,K) = 1. / RAM_SOIL(I,K)
!
      ENDDO
      ENDDO
!
      RETURN
      END

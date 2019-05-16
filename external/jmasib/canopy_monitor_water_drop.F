! CANOPY のフラックス諸量のモニタ。
!  (4) canopy_monitor_water_drop
!      h2o の葉からの落下
!
! =====================================================================
!
      SUBROUTINE CANOPY_MNT_WATER_DROP (
     I    NSTR     , NEND      , 
     I    INSFLAG  , 
     I    WTR_LEAF_DROP        , GLA_LEAF_DROP       , 
     I    TMP_CNP  ,             TMP_GRS  )  
!
      USE SIBPRM , ONLY : IDX , ICG , ICGB 
!
      USE SIBCON , ONLY :
     1   HL_WTR_c , 
     1   HL_W2G_c ,
     1   TMP_FREZ_C , 
     1   CAP_WTR_C , 
     1   CAP_GLA_C  
!
      USE SIB_MONIT , ONLY : MONIT_SIB, ICN, ICS , 
     1                       IMONIT_LEVEL , IMONIT_ALL 
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: INSFLAG 
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
!
      REAL(8),INTENT(IN)   :: WTR_LEAF_DROP(IDX,ICG) 
      REAL(8),INTENT(IN)   :: GLA_LEAF_DROP(IDX,ICG) 
!
      REAL(8),INTENT(IN)   :: TMP_CNP  (IDX) 
      REAL(8),INTENT(IN)   :: TMP_GRS  (IDX) 
!
      CHARACTER(2)      :: CFLAG_NS_G2
!
      REAL(8)           :: WORK  (IDX)
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
      INTEGER           :: I
      INTEGER           :: IFLG_MON
!
      IF ( INSFLAG .EQ.0 ) THEN
        CFLAG_NS_G2 ='NL'
        IFLG_MON    = ICN
      ELSE
        CFLAG_NS_G2 ='SW'
        IFLG_MON    = ICS
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! ==========================
! >> 葉から落下した水の量 <<         
! ==========================
!
        CMARK  = 'WLSBxx'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WTR_LEAF_DROP(1,1), IFLG_MON, ONE , ONE)   
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し下草
          CMARK  = 'WLSBGL'
          CALL MONIT_SIB ( CMARK, WTR_LEAF_DROP(1,2),IFLG_MON, ONE, ONE)   
        ENDIF
!
! ==========================
! >> 葉から落下した氷の量 <<         
! ==========================
!
        CMARK  = 'GLSBxx'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, GLA_LEAF_DROP(1,1), IFLG_MON, ONE, ONE)   
!
        IF   ( INSFLAG.EQ.0 ) THEN   ! 雪無し下草
          CMARK  = 'GLSBGL'
          CALL MONIT_SIB ( CMARK, GLA_LEAF_DROP(1,2),IFLG_MON, ONE,ONE)   
        ENDIF
!
! ==============================
! >> 葉から落下した水に伴う熱 <<
! ==============================
!
        DO I=NSTR,NEND
          WORK(I) = CAP_WTR_c * ( TMP_CNP(I) - TMP_FREZ_c )
     1                * WTR_LEAF_DROP(I,1)
        ENDDO 
        CMARK  = 'ILSBxx'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, ONE , ONE)   
!
        IF ( INSFLAG.EQ.0 ) THEN   ! 雪無し ・ 下草
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_GRS(I) - TMP_FREZ_c )
     1                  * WTR_LEAF_DROP(I,2)
          ENDDO 
          CMARK  = 'ILSBGL'
          CALL MONIT_SIB ( CMARK, WORK ,IFLG_MON, ONE , ONE)   
        ENDIF
!
! ==============================
! >> 葉から落下した氷に伴う熱 <<
! ==============================
!
        DO I=NSTR,NEND
          WORK(I) = ( CAP_GLA_c * ( TMP_CNP(I) - TMP_FREZ_c )
     1                - HL_W2G_C ) 
     1                * GLA_LEAF_DROP(I,1)
        ENDDO 
        CMARK  = 'JLSBxx'
        CMARK(5:6) = CFLAG_NS_G2
        CALL MONIT_SIB ( CMARK, WORK, IFLG_MON, ONE , ONE)   
!
        IF ( INSFLAG.EQ.0 ) THEN   ! 雪無し 
          DO I=NSTR,NEND
            WORK(I) = ( CAP_GLA_c * ( TMP_GRS(I) - TMP_FREZ_c )
     1                  - HL_W2G_C ) 
     1                * GLA_LEAF_DROP(I,2)
          ENDDO 
          CMARK  = 'JLSBGL'
          CALL MONIT_SIB ( CMARK, WORK ,IFLG_MON, ONE , ONE)   
        ENDIF
!
      ENDIF     
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_CANOPY_WATER_DROP 
!
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL 
!
! デバッグ用変数の登録作業
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_WATER_DROP 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ================
! >> 落下した水 <<
! ================
!
! ------
! > 水 <
! ------
!
! キャノピー
!
      CMARK  = 'WLSBNL'
      CTITLE = 'WATER LOSS OF CANOPY NOS'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WLSBSW'
      CTITLE = 'WATER LOSS OF CANOPY SNW'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
      CMARK  = 'WLSBGL'
      CTITLE = 'WATER LOSS OF GRASS'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ------
! > 氷 <
! ------
! 
! キャノピー
!
      CMARK  = 'GLSBNL'
      CTITLE = 'ICE LOSS OF CANOPY NOS ( INTO SNW )'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GLSBSW'
      CTITLE = 'ICE LOSS OF CANOPY SNW'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
      CMARK  = 'GLSBGL'
      CTITLE = 'ICE LOSS OF GRASS ( INTO SNW )'
      CUNIT  = 'KG/M*2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ==================================
! >> 葉から落下した水・氷に伴う熱 <<
! ==================================
!
! ------
! > 水 <
! ------
!
! キャノピー
!
      CMARK  = 'ILSBNL'
      CTITLE = 'HEAT LOSS OF CANOPY ( WATER DROP ) NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'ILSBSW'
      CTITLE = 'HEAT LOSS OF CANOPY ( WATER DROP ) SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
      CMARK  = 'ILSBGL'
      CTITLE = 'HEAT LOSS OF GRASS ( WATER DROP )'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ------
! > 氷 <
! ------
!
! キャノピー 
!
      CMARK  = 'JLSBNL'
      CTITLE = 'HEAT LOSS OF CANOPY ( ICE DROP ) NOS ( INTO SNW )'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JLSBSW'
      CTITLE = 'HEAT LOSS OF CANOPY ( ICE DROP ) SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
      CMARK  = 'JLSBGL'
      CTITLE = 'HEAT LOSS OF GRASS ( ICE DROP INTO SNW )'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      RETURN
      END

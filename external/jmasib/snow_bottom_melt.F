      SUBROUTINE SNOW_BOTTOM_MELT (
     I    NMAX          , NMAX_0        , NMAX_01   , 
     I    CVR_SNW       , ENG_SOIL_TOP  , 
     M    INFO_SNOW_NEXT, 
     M    ENG_SNOW_NEXT , H2O_SNOW_NEXT , 
     M    ENG_BUCKET_INTO_SOIL    , ENG_BUCKET_SNOW_TO_SNSL , 
     M    H2O_BUCKET_INTO_SOIL    , WTR_BUCKET_SNOW_TO_SNSL  )  
!
      USE SIBPRM , ONLY : IDX , ISN  
!
      USE SIBCON , ONLY : 
     1   ITUNE_SNOW_BOTTOM_MELT , 
     1   HL_W2G_c , 
! 20050803
!     1   TMP_FREZ_c , 
     1   CAP_GLA_c  ,
     1   TMP_FREZ_c  , 
     1   ENG_MELT_MAX_c 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, INS , INN , 
     1                       IMONIT_LEVEL, IMONIT_ALL, IMONIT_FLUX_MRI
#endif
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NMAX_0 (3)
      INTEGER,INTENT(IN)    :: NMAX_01(3)
!
      REAL(8),INTENT(IN)    :: CVR_SNW       (IDX)
      REAL(8),INTENT(IN)    :: ENG_SOIL_TOP  (IDX)
!
      INTEGER,INTENT(INOUT) :: INFO_SNOW_NEXT(IDX)
!
      REAL(8),INTENT(INOUT) :: ENG_SNOW_NEXT (IDX,ISN)
      REAL(8),INTENT(INOUT) :: H2O_SNOW_NEXT (IDX,ISN)
!
      REAL(8),INTENT(INOUT) :: ENG_BUCKET_SNOW_TO_SNSL (IDX)
      REAL(8),INTENT(INOUT) :: WTR_BUCKET_SNOW_TO_SNSL (IDX)
!
      REAL(8),INTENT(INOUT) :: ENG_BUCKET_INTO_SOIL    (IDX)
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO_SOIL    (IDX)
!
      INTEGER  ::  I 
      INTEGER  ::  ILEV
! 
      REAL(8)  ::  ENG_MELT  (IDX)   ! 土壌から雪へのエネルギー(下向き正) = 負
      REAL(8)  ::  H2O_MELT  (IDX)   ! 融雪水(下向き正) = 正
!
      REAL(8)  ::  ENG_SOIL_TMP(IDX) 
      REAL(8)  ::  XXX_SNOW_BOT
      LOGICAL,SAVE  ::  LFIRST = .TRUE. 
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK
      REAL(8),PARAMETER :: ONE = 1.D0 
#endif
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_BOTTOM_MELT ITUNE_SNOW_BOTTOM_MELT=' , 
     1             ITUNE_SNOW_BOTTOM_MELT 
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( ITUNE_SNOW_BOTTOM_MELT .EQ. 0 ) THEN
        RETURN
      ENDIF
!
      H2O_MELT    (:) = 0. 
      ENG_MELT    (:) = 0. 
!
      ENG_SOIL_TMP(1:NMAX) = ENG_SOIL_TOP(1:NMAX) 
      DO I=NMAX_01(2),NMAX_01(3)
        ENG_SOIL_TMP(I) = ENG_SOIL_TOP(I) * CVR_SNW(I) 
      ENDDO 
!
! ============
! >> 雪一層 <<
! ============
!   全部融かすこともありうる
!       
      DO I=1,NMAX
        IF ( ENG_SOIL_TMP(I).GT.0 .AND. INFO_SNOW_NEXT(I).EQ.1 ) THEN
          XXX_SNOW_BOT = ENG_SNOW_NEXT(I,1)
!
#ifdef SIB_DEBUG
          IF ( XXX_SNOW_BOT .GT. 1.D-20 ) THEN
            write(6,*) 'snow_bottom_melt error ', i, info_snow_next(i), 
     1                  xxx_snow_bot  
          ENDIF 
#endif    
!
!  全部融け、TF の水が下へ。
!
          IF   ( ENG_SOIL_TMP(I) .GE. - XXX_SNOW_BOT ) THEN
            ENG_MELT(I) = ENG_SNOW_NEXT(I,1)      ! < 0
            H2O_MELT(I) = H2O_SNOW_NEXT(I,1)      ! > 0 
            INFO_SNOW_NEXT(I) = -1                ! かな？ 
            H2O_SNOW_NEXT(I,1) = 0. 
            ENG_SNOW_NEXT(I,1) = 0.
!
!  一部融ける
!
          ELSE 
            ENG_MELT(I) = - ENG_SOIL_TMP(I)                     ! < 0
            H2O_MELT(I) = ENG_MELT(I) / ( ENG_SNOW_NEXT(I,1) )  ! > 0 
     1                    * H2O_SNOW_NEXT(I,1) 
            H2O_SNOW_NEXT(I,1) = H2O_SNOW_NEXT(I,1) - H2O_MELT(I) 
            ENG_SNOW_NEXT(I,1) = ENG_SNOW_NEXT(I,1) - ENG_MELT(I) 
          ENDIF
!
        ENDIF
      ENDDO
!
! ================
! >> 雪二層以上 <<
! ================
!   全部融かしはしない。融雪量に上限値有り。
!       
      DO I=1,NMAX
        IF ( ENG_SOIL_TMP(I).GT.0 .AND. INFO_SNOW_NEXT(I).GE.2 ) THEN
!
          ILEV = INFO_SNOW_NEXT(I)
          XXX_SNOW_BOT = ENG_SNOW_NEXT(I,ILEV)
!
#ifdef SIB_DEBUG
          IF ( XXX_SNOW_BOT .GT. 1.D-20 ) THEN
            write(6,*) 'snow_bottom_melt error ', i, info_snow_next(i), 
     1                  xxx_snow_bot  
          ENDIF 
#endif    
!
          ENG_MELT(I) = MIN ( ENG_SOIL_TMP(I) , 
     1                        - XXX_SNOW_BOT  ,
     1                        - ENG_MELT_MAX_c )              ! < 0 
          H2O_MELT(I) = ENG_MELT(I) / ENG_SNOW_NEXT(I,ILEV)   ! > 0 
     1                  * H2O_SNOW_NEXT(I,ILEV) 
!
          H2O_SNOW_NEXT(I,ILEV) = H2O_SNOW_NEXT(I,ILEV) 
     1                              - H2O_MELT(I) 
          ENG_SNOW_NEXT(I,ILEV) = ENG_SNOW_NEXT(I,ILEV) 
     1                              - ENG_MELT(I) 
!
        ENDIF
      ENDDO
!
! ================
! >> バケツ行き <<
! ================
! 雪なし域なら雪なし土壌へ
! 雪が少しでもあれば雪あり土壌へ
!
      DO I=1,NMAX_0(3) 
        ENG_BUCKET_INTO_SOIL   (I) 
     1      = ENG_BUCKET_INTO_SOIL   (I) + ENG_MELT(I) 
        H2O_BUCKET_INTO_SOIL   (I) 
     1      = H2O_BUCKET_INTO_SOIL   (I) + H2O_MELT(I) 
      ENDDO

      DO I=NMAX_01(2),NMAX 
        ENG_BUCKET_SNOW_TO_SNSL(I) 
     1      = ENG_BUCKET_SNOW_TO_SNSL(I) + ENG_MELT(I) 
        WTR_BUCKET_SNOW_TO_SNSL(I)
     1      = WTR_BUCKET_SNOW_TO_SNSL(I) + H2O_MELT(I) 
      ENDDO
!
! ============
! >> モニタ <<
! ============
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
        CMARK = 'WMLBN'
        CALL MONIT_SIB ( CMARK , H2O_MELT , INN , ONE , ONE ) 
        CMARK = 'WMLBS'
        CALL MONIT_SIB ( CMARK , H2O_MELT , INS , ONE , ONE ) 
        CMARK = 'HMLBN'
        CALL MONIT_SIB ( CMARK , ENG_MELT , INN , ONE , ONE ) 
        CMARK = 'HMLBS'
        CALL MONIT_SIB ( CMARK , ENG_MELT , INS , ONE , ONE ) 
      ENDIF
#endif
!
      RETURN
      END
!============================================================
      SUBROUTINE MR_SNOW_BOTTOM_MELT
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL , IMONIT_ALL , 
     1                      IMONIT_FLUX_MRI   
#endif
!
      IMPLICIT NONE
!
#ifdef MONYOS
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
#endif
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_BOTTOM_MELT 02/03 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL      ) THEN
!
        CMARK  = 'WMLBN'
        CTITLE = 'Snow bottom melt(NOS)'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'WMLBS'
        CTITLE = 'Snow bottom melt(SNW)'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'HMLBN'
        CTITLE = 'Energy of Snow bottom melt(NOS)'
        CUNIT  = 'W/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'HMLBS'
        CTITLE = 'Energy of Snow bottom melt(SNW)'
        CUNIT  = 'W/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END
        

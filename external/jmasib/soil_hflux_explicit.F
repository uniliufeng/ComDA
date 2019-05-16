! 地面内部での温度の式をたてる。
! 99/06/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_HFLUX_EXPLICIT (
     I  NSTR       , NEND          , 
#ifdef MONYOS
     I  INSFLAG    ,
#endif
     I  CVR_SorN_INV , 
     I  TMP_SOIL   , TMP_SOIL4 , 
     I  HDF_SOIL_H , 
     I  ENG_BUCKET_TOP_H , 
     I  DELT_SOIL  , DELT_SOIL_INV , 
     O  HEAT_CONV )
!
!    注意  
!      ・雪なし領域、雪あり領域とも、これを使う。
!      ・上側は熱フラックスが与えられている。 
!      ・水フラックスはまだ入ってこない。
!  
!        D1    1  
!        D2    2  
!        D3    3  
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IDP     
!
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS, INN, INS, 
     1                       IMONIT_LEVEL, IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON
#endif
!
      use sibcon, only : itune_soilfluxbottom_c 
!
      IMPLICIT NONE 
!
! INPUT 
      INTEGER,INTENT(IN)    :: NSTR               ! 開始位置	
      INTEGER,INTENT(IN)    :: NEND               ! 終了位置 
      REAL(8),INTENT(IN)    :: DELT_SOIL          ! SIB用インプリシット時間
      REAL(8),INTENT(IN)    :: DELT_SOIL_INV      ! SIB用インプリシット時間
#ifdef MONYOS
      INTEGER,INTENT(IN)    :: INSFLAG            ! SIB用インプリシット時間
#endif
!
      REAL(8),INTENT(IN)    :: CVR_SorN_INV (IDX)
!
      REAL(8),INTENT(IN)    :: TMP_SOIL   (IDX,IDP)   ! 土壌温度
      REAL(8),INTENT(IN)    :: TMP_SOIL4  (IDX)     
      REAL(8),INTENT(IN)    :: HDF_SOIL_H (IDX,IDP)   ! 土 1-2,2-3 間熱係数
      REAL(8),INTENT(IN)    :: ENG_BUCKET_TOP_H (IDX)   ! 地面上端流入熱量
!                                                     !     * DELT_SOIL 済
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: HEAT_CONV  (IDX,IDP+1)  ! エクスプリシット収束
                                                       ! *DELT はかかっていない
!
! WORK
      INTEGER               :: I
      REAL(8)               :: ENG_BUC (IDX)      ! 地面上端流入熱量
      REAL(8)               :: HFLUX_H (IDX,IDP)  ! エクスプリシット flux
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
#endif
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_HFLUX_EXPLICIT 00/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  ----------------------------------
!  > フラックス（エクスプリシット） <
!  ----------------------------------
!
      DO I=NSTR,NEND
        HFLUX_H(I,1) = HDF_SOIL_H(I,1)*( TMP_SOIL(I,1) - TMP_SOIL(I,2) )      
        HFLUX_H(I,2) = HDF_SOIL_H(I,2)*( TMP_SOIL(I,2) - TMP_SOIL(I,3) )      
        HFLUX_H(I,3) = HDF_SOIL_H(I,3)*( TMP_SOIL(I,3) - TMP_SOIL4(I) )      
     1                 * itune_soilfluxbottom_c 
        ENG_BUC(I)  =  ENG_BUCKET_TOP_H(I)*DELT_SOIL_INV*CVR_SorN_INV(I)
      ENDDO
!
!  --------------------
!  > 熱フラックス収束 <
!  --------------------
!
      DO I=NSTR    ,NEND    
        HEAT_CONV(I,1) = ENG_BUC(I)   - HFLUX_H(I,1)
        HEAT_CONV(I,2) = HFLUX_H(I,1) - HFLUX_H(I,2) 
        HEAT_CONV(I,3) = HFLUX_H(I,2) - HFLUX_H(I,3) 
        HEAT_CONV(I,4) = HFLUX_H(I,3) 
      ENDDO
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
!
! --------------------
! > 雪あり無し域平均 <
! --------------------

      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'FLG0'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H, INN,ONE ,ONE)
          CMARK  = 'FLG1'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1)    , ICN, DELT_SOIL,ONE)  
          CMARK  = 'FLG2'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),     ICN, DELT_SOIL,ONE)  
          CMARK  = 'FLG3'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3),     ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'FLG0'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H, INS,ONE, ONE)
          CMARK  = 'FLG1'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1),     ICS, DELT_SOIL,ONE)  
          CMARK  = 'FLG2'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),     ICS, DELT_SOIL,ONE)  
          CMARK  = 'FLG3'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3),     ICS, DELT_SOIL,ONE)  
        ENDIF
      ENDIF
!
! ----------
! > 領域別 <
! ----------
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'HFLBLN'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H, INN,ONE ,ONE)
          CMARK  = 'HFLN12F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1)    , ICN, DELT_SOIL,ONE)  
          CMARK  = 'HFLN23F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),     ICN, DELT_SOIL,ONE)  
          CMARK  = 'HFLN34F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3),     ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'HFLBLW'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H, INS,ONE, ONE)
          CMARK  = 'HFLS12F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1),     ICS, DELT_SOIL,ONE)  
          CMARK  = 'HFLS23F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),     ICS, DELT_SOIL,ONE)  
          CMARK  = 'HFLS34F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3),     ICS, DELT_SOIL,ONE)  
        ENDIF
!
      ENDIF
#endif
!
      RETURN
      END
!===============================================================
      SUBROUTINE MR_SOIL_HFLUX_EXPLICIT 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL, IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON 
#endif
!
! soil_hflux_eqset で出力するデバッグ用変数の登録作業
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MR_SOIL_HFLUX_EXPLICIT 00/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!  土壌中熱伝導
!
        CMARK  = 'FLG0'
        CTITLE = 'Heat Flux soil surface - soil1'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'FLG1'
        CTITLE = 'Heat Flux soil1 - soil2'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'FLG2'
        CTITLE = 'Heat Flux soil2 - soil3'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'FLG3'
        CTITLE = 'Heat Flux soil3 - bottom'
        CUNIT  = 'W/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'HFLBLN'
      CTITLE = 'HEAT FLUX TO SOIL 1 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN12F'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN23F'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN34F'
      CTITLE = 'HEAT FLUX FROM SOIL 3 TO SOIL 4 NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLBLW'
      CTITLE = 'HEAT FLUX TO SOIL 1 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS12F'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS23F'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS34F'
      CTITLE = 'HEAT FLUX FROM SOIL 3 TO SOIL 4 SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END


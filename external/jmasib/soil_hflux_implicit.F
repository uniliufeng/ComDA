! 地面内部での温度の式をたてる。
! 99/06/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_HFLUX_IMPLICIT_EQSET (
     I  NSTR       , NEND          , 
#ifdef MONYOS
     I  INSFLAG    ,
#endif
     I  CVR_SorN_INV , 
     I  TMP_SOIL   , TMP_SOIL4     , 
     I  PHASE_SOIL , 
     I  CAP_SOIL   , 
     I  HDF_SOIL_H , 
     I  ENG_BUCKET_TOP_H , 
     I  DELT_SOIL  , DELT_SOIL_INV , 
     O  SIB_MTX_SOIL  , HFLUX_H    )

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
      USE SIBCON , ONLY :
     1    DELT_FAC_SOIL_c,  itune_soilfluxbottom_c 

!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS, INN, INS, INA, 
     1                       IMONIT_LEVEL , IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON,
     1                       IMONIT_ALMA 
#endif
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
      REAL(8),INTENT(IN)    :: CAP_SOIL   (IDX,IDP)   ! 土壌熱容量
      REAL(8),INTENT(IN)    :: HDF_SOIL_H (IDX,IDP)   ! 土 1-2,2-3 間熱係数
      REAL(8),INTENT(IN)    :: ENG_BUCKET_TOP_H (IDX)   ! 地面上端流入熱量
!                                                     !     * DELT_SOIL 済
      REAL(8),INTENT(IN)    :: PHASE_SOIL (IDX,IDP)   ! 零度なら 0 , 非 1  
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: SIB_MTX_SOIL   (IDX,IDP,IDP+1) 
!
! WORK
      INTEGER               :: I
      REAL(8)               :: DELT_TMP
      REAL(8),intent(out)   :: HFLUX_H (IDX,0:IDP)    ! エクスプリシット flux
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
        WRITE(6,*) 'SOIL_HFLUX_EQSET 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DELT_TMP = DELT_SOIL * DELT_FAC_SOIL_c
!
!  -----------------------------
!  フラックス（エクスプリシット）
!  -----------------------------
!
      DO I=NSTR,NEND
        HFLUX_H(I,0) = ENG_BUCKET_TOP_H(I)*DELT_SOIL_INV*CVR_SorN_INV(I)
        HFLUX_H(I,1) = HDF_SOIL_H(I,1)*( TMP_SOIL(I,1) - TMP_SOIL(I,2) )      
        HFLUX_H(I,2) = HDF_SOIL_H(I,2)*( TMP_SOIL(I,2) - TMP_SOIL(I,3) )      
        HFLUX_H(I,3) = HDF_SOIL_H(I,3)*( TMP_SOIL(I,3) - TMP_SOIL4(I)  )      
     1                 * itune_soilfluxbottom_c 
      ENDDO
!
      call reset ( SIB_MTX_SOIL , IDX*IDP*(IDP+1) )
!
!  ---------------------------------------------------
!  土壌第一層の式 
!  ---------------------------------------------------
!
      DO 1000 I=NSTR    ,NEND    
        SIB_MTX_SOIL(I,1,1)                         ! 土壌第一層
     1   =   CAP_SOIL(I,1) 
     1     + HDF_SOIL_H(I,1) * DELT_TMP                ! 土壌第二層と
     1         * PHASE_SOIL(I,1)
!
        SIB_MTX_SOIL(I,1,2)                         ! 土壌第二層
     1   = - HDF_SOIL_H(I,1) * DELT_TMP 
     1         * PHASE_SOIL(I,2)
!
        SIB_MTX_SOIL(I,1,4)                                  ! 強制
!    1   =  ENG_BUCKET_TOP_H(I)*DELT_SOIL_INV*CVR_SorN_INV(I)
     1   =  HFLUX_H(I,0) - HFLUX_H(I,1)
 1000 CONTINUE         
!
!  ---------------------------------------------------
!  土壌第二層の式 
!  ---------------------------------------------------
!
      DO 2000 I=NSTR    ,NEND    
!
        SIB_MTX_SOIL(I,2,1)                        ! 土壌第一層   
     1    = - HDF_SOIL_H(I,1) * DELT_TMP 
     1                          * PHASE_SOIL(I,1)
!
        SIB_MTX_SOIL(I,2,2)                        ! 土壌第二層 
     1     =   CAP_SOIL  (I,2) 
     1       + HDF_SOIL_H(I,1) * DELT_TMP             ! 土壌第一層と
     1                         * PHASE_SOIL(I,2)
     1       + HDF_SOIL_H(I,2) * DELT_TMP               ! 土壌第三層と
     1                         * PHASE_SOIL(I,2)
!
        SIB_MTX_SOIL(I,2,3)                       ! 土壌第三層
     1     = - HDF_SOIL_H(I,2) * DELT_TMP
     1                         * PHASE_SOIL(I,3)
!
        SIB_MTX_SOIL(I,2,4) =                     ! 強制
     1      +  HFLUX_H(I,1) -  HFLUX_H(I,2) 
 2000 CONTINUE         
!
!  ---------------------------------------------------
!  土壌第三層の式 
!  ---------------------------------------------------
!
      DO 3000 I=NSTR    ,NEND    
!
        SIB_MTX_SOIL (I,3,2)                       ! 土壌第二層
     1     = - HDF_SOIL_H(I,2) * DELT_TMP
     1                         * PHASE_SOIL(I,2)
!
        SIB_MTX_SOIL (I,3,3)                       ! 土壌第三層
     1     =   CAP_SOIL  (I,3) 
     1       + HDF_SOIL_H(I,2) * DELT_TMP
     1                         * PHASE_SOIL(I,3)
!
        SIB_MTX_SOIL (I,3,4)                       ! 強制
     1     = HFLUX_H(I,2) -  HFLUX_H(I,3) 
 3000 CONTINUE         
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'FLG0'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H,INN, ONE , ONE)
        ELSE
          CMARK  = 'FLG0'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H,INS, ONE , ONE)
        ENDIF
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK  = 'HFLBLN'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H,INN, ONE , ONE)
          CMARK  = 'HFLN12'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1),    ICN, DELT_SOIL,ONE)  
          CMARK  = 'HFLN23'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),    ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'HFLBLW'
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H, INS, ONE , ONE)
          CMARK  = 'HFLS12'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1),     ICS, DELT_SOIL,ONE)  
          CMARK  = 'HFLS23'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2),     ICS, DELT_SOIL,ONE)  
        ENDIF
!
      ENDIF
!
!* ALMA 
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_GH'
        IF ( INSFLAG.EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H,INN, ONE , ONE)
        ELSE
          CALL MONIT_SIB ( CMARK, ENG_BUCKET_TOP_H,INS, ONE , ONE)
        ENDIF
      ENDIF
#endif
!
      RETURN
      END
!===============================================================
! 地中の方程式を解く
      SUBROUTINE SOIL_HFLUX_IMPLICIT_EQSLV (
     I  NSTR  , NEND  , 
     M  SMTX      ) 
!
!     雪なし、雪あり、それぞれで呼ばれる。
!
!     D1  1
!     D2  2
!     D3  3
!
      USE SIBPRM , ONLY : 
     1    IDX      
!
      IMPLICIT NONE 
!
! INPUT
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
!
! INOUT
!
      REAL(8),INTENT(INOUT) :: SMTX(IDX,3,4) 
!
! WORK
      INTEGER               :: I
      REAL(8)               :: A22I 
      REAL(8)               :: A33I 
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_HFLUX_EQSLV 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=NSTR,NEND
!
!a11  a12       a14
!a21  a22  a23  a24
!     a32  a33  a34
!
        A33I = 1. / SMTX(I,3,3) 
        SMTX(I,3,2) = SMTX(I,3,2) * A33I  
        SMTX(I,3,4) = SMTX(I,3,4) * A33I  
!       SMTX(I,3,3) = 1. 
!
        SMTX(I,2,2) = SMTX(I,2,2) - SMTX(I,2,3) * SMTX(I,3,2) 
        SMTX(I,2,4) = SMTX(I,2,4) - SMTX(I,2,3) * SMTX(I,3,4) 
!       SMTX(I,2,3) = 0. 
!
!a11  a12       a14
!a21  a22*  0*  a24*
!     a32*  1*  a34*
!
        A22I = 1. / SMTX(I,2,2) 
        SMTX(I,2,1) = SMTX(I,2,1) * A22I  
        SMTX(I,2,4) = SMTX(I,2,4) * A22I  
!       SMTX(I,2,2) = 1. 
!
        SMTX(I,1,1) = SMTX(I,1,1) - SMTX(I,1,2) * SMTX(I,2,1) 
        SMTX(I,1,4) = SMTX(I,1,4) - SMTX(I,1,2) * SMTX(I,2,4) 
!       SMTX(I,1,2) = 0. 
!
!a11*  0*       a14*
!a21*  1*       a24*
!     a32   1   a34
!
        SMTX(I,1,4) = SMTX(I,1,4) / SMTX(I,1,1)
!       SMTX(I,1,1) = 1. 
!
! 1*            a14*
!a21   1        a24
!     a32   1   a34
!
        SMTX(I,2,4) = SMTX(I,2,4) - SMTX(I,2,1) * SMTX(I,1,4)  
!       SMTX(I,2,1) = 0. 
        SMTX(I,3,4) = SMTX(I,3,4) - SMTX(I,3,2) * SMTX(I,2,4)  
!       SMTX(I,3,2) = 0. 
 1000 CONTINUE
!
      RETURN
      END
!===============================================================
! 土壌中の熱フラックス収束 （水の移動を考慮しない時点）
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SOIL_HFLUX_IMPLICIT_CONV (
     I  NSTR            , NEND            , 
#ifdef MONYOS
     I  INSFLAG         ,
#endif
     I  SIB_MTX_SOIL    , 
     I  HDF_SOIL_H      , 
     I  PHASE_SOIL      , 
     I  DELT_SOIL       , 
     M  HFLUX_H         ,
     O  HEAT_CONV    )  
!
      USE SIBPRM , ONLY :
     1    IDX   ,
     1    IDP  
!
      USE SIBCON , ONLY :
     1   DELT_FAC_SOIL_c
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS, INN, INS, 
     1                       IMONIT_LEVEL , IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON 
#endif
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      REAL(8),INTENT(IN)    :: SIB_MTX_SOIL     (IDX,IDP,IDP+1)
      REAL(8),INTENT(IN)    :: HDF_SOIL_H       (IDX,IDP-1)
      REAL(8),INTENT(IN)    :: PHASE_SOIL       (IDX,IDP)
!
      REAL(8),INTENT(IN)    :: DELT_SOIL                
!
#ifdef MONYOS
      INTEGER,INTENT(IN)    :: INSFLAG            ! SIB用インプリシット時間
#endif
!
! INOUT
!   
      REAL(8),INTENT(INOUT) :: HFLUX_H   (IDX,0:IDP) ! 熱フラックス
      REAL(8),INTENT(OUT)   :: HEAT_CONV (IDX,IDP+1) ! 熱フラックス収束
!
! WORK
!
      INTEGER   :: I
      REAL(8)   :: DELT_TMP
!      REAL(8)   :: HFLX_COR_12 (IDX)
!      REAL(8)   :: HFLX_COR_23 (IDX)
!      REAL(8)   :: HFLX_COR_34 
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
!     REAL(8)           :: WORK(IDX) 
#endif
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'SOIL_HFLUX_CONV 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! エクスプリシット評価に、インプリシット分を加える。
!
      DELT_TMP = DELT_SOIL * DELT_FAC_SOIL_c
!
      DO 1000 I=NSTR,NEND
!
        HFLUX_H(I,1) = HFLUX_H(I,1) +  
     1     HDF_SOIL_H(I,1) * DELT_TMP
     1        * (   SIB_MTX_SOIL(I,1,4) * PHASE_SOIL(I,1) 
     1            - SIB_MTX_SOIL(I,2,4) * PHASE_SOIL(I,2)     )
!
        HFLUX_H(I,2) = HFLUX_H(I,2) +  
     1     HDF_SOIL_H(I,2) * DELT_TMP
     1        * (   SIB_MTX_SOIL(I,2,4) * PHASE_SOIL(I,2) 
     1            - SIB_MTX_SOIL(I,3,4) * PHASE_SOIL(I,3)     )
!
!       HFLUX_H(I,3) = HFLUX_H(I,3)   
!
        HEAT_CONV(I,1) = HFLUX_H(I,0) - HFLUX_H(I,1)
        HEAT_CONV(I,2) = HFLUX_H(I,1) - HFLUX_H(I,2)
        HEAT_CONV(I,3) = HFLUX_H(I,2) - HFLUX_H(I,3)
        HEAT_CONV(I,4) = HFLUX_H(I,3) 
 1000 CONTINUE   
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
!
! 格子平均出力
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
        IF ( INSFLAG.EQ.0) THEN
          CMARK  = 'FLG1'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1), ICN, DELT_SOIL, ONE)  
          CMARK  = 'FLG2'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2), ICN, DELT_SOIL, ONE)  
          CMARK  = 'FLG3'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3), ICN, DELT_SOIL, ONE)  
        ELSE
          CMARK  = 'FLG1'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1), ICS, DELT_SOIL, ONE)  
          CMARK  = 'FLG2'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2), ICS, DELT_SOIL, ONE)  
          CMARK  = 'FLG3'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3), ICS, DELT_SOIL, ONE)  
        ENDIF
      ENDIF
!
! 領域別出力
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
        IF ( INSFLAG.EQ.0) THEN
          CMARK  = 'HFLN12F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1), ICN, DELT_SOIL, ONE)  
          CMARK  = 'HFLN23F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2), ICN, DELT_SOIL, ONE)  
          CMARK  = 'HFLN34F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3), ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'HFLS12F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,1), ICS, DELT_SOIL, ONE)  
          CMARK  = 'HFLS23F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,2), ICS, DELT_SOIL, ONE)  
          CMARK  = 'HFLS34F'
          CALL MONIT_SIB ( CMARK, HFLUX_H(1,3), ICS, DELT_SOIL, ONE)  
        ENDIF
      ENDIF
#endif
!
      RETURN
      END
!===============================================================
      SUBROUTINE MR_SOIL_HFLUX_IMPLICIT 
!
! soil_hflux_eqset で出力するデバッグ用変数の登録作業
!
      USE SIB_MONIT, ONLY :  IMONIT_LEVEL , IMONIT_ALL,
     1                       IMONIT_FLUX_COMMON , IMONIT_ALMA
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
        WRITE(6,*) 'MR_SOIL_HFLUX_IMPLICIT 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  熱フラックス  HFLXG01,HFLXG12,HFLXG23
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
      CTITLE = 'HEAT FLUX TO SOIL 1 NOS explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN12'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 NOS explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN23'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 NOS explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLBLW'
      CTITLE = 'HEAT FLUX TO SOIL 1 SNW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS12'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 SNW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS23'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 SNW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  熱フラックス補正後
!
      CMARK  = 'HFLN12F'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 NOS after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN23F'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 NOS after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLN34F'
      CTITLE = 'HEAT FLUX FROM SOIL 3 TO SOIL 4 NOS after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS12F'
      CTITLE = 'HEAT FLUX FROM SOIL 1 TO SOIL 2 SNW after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS23F'
      CTITLE = 'HEAT FLUX FROM SOIL 2 TO SOIL 3 SNW after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLS34F'
      CTITLE = 'HEAT FLUX FROM SOIL 3 TO SOIL 4 SNW after CORRECT'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_GH'
        CTITLE = 'ALMA Qg '
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF 
!
      RETURN
      END


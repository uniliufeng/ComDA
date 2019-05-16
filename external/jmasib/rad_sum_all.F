! 放射吸収量を計算する。
! 99/05/xx 保坂征宏
! 99/07/26 保坂征宏  キャノピーならびにスキンの吸収量の加算に
!                    雪層への入射放射は加算しないことに注意
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/30 保坂征宏 符号の向きチェック。収束が正を確認。
!
      SUBROUTINE RAD_SUM_ALL (
     I  NSTR       , NEND       , INSFLAG     ,  
#ifdef MONYOS
     I  DELT_CNP   , 
#endif 
     I  NTYP       , 
     I  TMP_CNP    , TMP_GRD    , 
     I  CVR_VEG_ms , CVR_VEG_INV_ms  , 
     I  RAD_LONG_INTO_SIB       , ZLT  ,    
     M  RAD_SHRT   , 
     O  THERMK     , 
     O  RAD_SHRT_INTO_SNOW_LAYER, 
     O  RAD_LONG_SIB_2_ATM , 
     O  RAD_LONG   )
!
!   放射吸収量を計算する。
!     長波放射吸収量の計算。放射計算後。
!          
      USE SIBPRM , ONLY : 
     1   IDX     ,
     1   ITYP    , 
     1   ICG  
!
      USE SIBCON , ONLY :
     1   STEFAN_c      ,
     1   CAP_SNSK_c  , 
     1   CAP_GLA_c_INV ,
     1   TMP_FREZ_c    ,
     1   FACT_TRANS_SR_SNOW_SKIN_c  
!
      USE SIBVEG , ONLY :
     1           CHI_L_co               , ! χL  (ITYP,ICG) 
     1           XLAI_mo                  ! LAI, 格子平均だ!! (ITYP,IMON,ICG)
!
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS,
     1    IMONIT_LEVEL   , 
     1    IMONIT_FLUX_MRI , 
     1    IMONIT_BUDJET , 
     1    IMONIT_ALL    , IMONIT_ALMA 
#endif
!
      IMPLICIT NONE
!
! 機能 : CG それぞれでの短波放射吸収量 RADT に対して、
!        長波放射吸収量を加える。
!        (1-キャノピーでの射出率) である thermk も計算。
!
! 注意 : 雪あり域と雪無し域それぞれについて呼ばれる。
!        雪あり域では TG は TS という意味である。
!
! コメント
!    ・地面 ( 草と裸地 ) は一体として取り扱う。
!    ・1-2-71 にしたがってキャノピー射出率を計算している。
!      EXTINC は G(μbar) ではない。 G(μbar)/μbar を一気に計算している。
!
! 入力
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: NTYP       (IDX) 
      INTEGER,INTENT(IN)    :: INSFLAG               ! 雪無し 0 雪あり 1 
!
!
      REAL(8),INTENT(IN)    :: TMP_CNP    (IDX)      ! キャノピー温度
      REAL(8),INTENT(IN)    :: TMP_GRD    (IDX)      ! 地表面上端温度（スキン）
      REAL(8),INTENT(IN)    :: ZLT        (IDX,ICG)  ! 
      REAL(8),INTENT(IN)    :: CVR_VEG_ms (IDX,ICG)  ! 被覆率((I,1)のみ使用 )
      REAL(8),INTENT(IN)    :: CVR_VEG_INV_ms (IDX,ICG) ! 
      REAL(8),INTENT(INOUT) :: RAD_SHRT   (IDX,ICG)  ! 短波吸収量 BD,WV合計。
      REAL(8),INTENT(IN)    :: RAD_LONG_INTO_SIB (IDX) ! 大気下端長波
!
#ifdef MONYOS
      REAL(8),INTENT(IN)   :: DELT_CNP
#endif
!
! 出力
!
      REAL(8),INTENT(OUT)  ::  THERMK     (IDX)        ! 1 - キャノピー射出率
      REAL(8),INTENT(OUT)   :: RAD_SHRT_INTO_SNOW_LAYER (IDX)  
      REAL(8),INTENT(OUT)   :: RAD_LONG_SIB_2_ATM (IDX)  
!
      REAL(8),INTENT(OUT)  :: RAD_LONG    (IDX,ICG)    ! 
! 
! ワーク
!
      REAL(8)              :: RAD_A2C(IDX)
      REAL(8)              :: RAD_A2G(IDX)
      REAL(8)              :: RAD_C2A(IDX)
      REAL(8)              :: RAD_C2G(IDX)
      REAL(8)              :: RAD_G2A(IDX)
      REAL(8)              :: RAD_G2C(IDX)

      REAL(8)              :: TC4 
      REAL(8)              :: TG4
      REAL(8)              :: AA 
      REAL(8)              :: BB 
      REAL(8)              :: EXTINC
      REAL(8)              :: ZKAT
      REAL(8)              :: FAC1 
      INTEGER              :: I   
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      INTEGER           :: IFLG_MON
#endif
!
      REAL(8),PARAMETER   :: XMU_TF   =  20.D0      ! tmp_frez 
      REAL(8),PARAMETER   :: XMU_GLA  = 100.D0      ! < tmp_frez_c 
      REAL(8),PARAMETER   :: DELT_TEMP_INV = 0.5D0  ! = 1/2K
      REAL(8),PARAMETER   :: ONE  = 1.D0
      REAL(8),PARAMETER   :: ZERO = 0.D0
      REAL(8),SAVE        :: DEP
!
      REAL(8)             :: FACT
      REAL(8)             :: XMU
      REAL(8)             :: TRANS
      LOGICAL,SAVE        :: LFIRST 
!     real(8) , save :: stefan_c = 5.67D-8 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        DEP = CAP_SNSK_c * CAP_GLA_c_INV              ! 質量 kg
     1                    / 200                       ! 厚さ 、単位 m に
!
        WRITE(6,*) 'RAD_SUM_ALL short wave snow skin trans rate ' , 
     1              ' cap_snsk_c = ' , CAP_SNSK_c      ,
     1              ' cap_gla_c  = ' , 1/CAP_GLA_c_INV , 
     1              ' depth(m)   = ' , DEP , 
     1              ' XMU_GLA    = ' , XMU_GLA ,   
     1              ' trns rate  = ' , DEXP ( - XMU_GLA * DEP ) , 
     1              ' fact = ' ,  FACT_TRANS_SR_SNOW_SKIN_c , 
     1              ' so ' , 
     1              FACT_TRANS_SR_SNOW_SKIN_c * DEXP (-XMU_GLA*DEP)  
!
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( INSFLAG.EQ.0 ) THEN
        IFLG_MON = ICN
      ELSE
        IFLG_MON = ICS
      ENDIF
#endif
!
      DO 1000 I=NSTR,NEND
!
        TC4  = TMP_CNP(I)**4 
        TG4  = TMP_GRD(I)**4 
!
!   -----------------------------------
!  >> キャノピーの射出率 (1-THERMK) <<
!   -----------------------------------
!                                               1-2-71 式の 1 - δt
        AA     = 0.500D0-0.633D0*CHI_L_co(NTYP(I),1)
     1                  -0.330D0*CHI_L_co(NTYP(I),1)**2  ! φ1
        BB     = 0.877D0-1.754D0*AA                      ! φ2
        EXTINC = 2 * AA + BB                             ! G(μbar)/μbar 
!       ZKAT   = EXTINC * XLAI_mo(NTYP(I),MON,1) * CVR_VEG_INV_ms(I,1)         
        ZKAT   = EXTINC * ZLT(I,1) * CVR_VEG_INV_ms(I,1)         
        ZKAT   = MIN( 50.0D0  , ZKAT  )
        ZKAT   = MAX( 10.0D-5 , ZKAT  )
        THERMK(I) = DEXP ( -ZKAT )                        ! ？ 
!
!   -----------------------------------------------
!  >> 長波放射で、キャノピーの吸収する割合 FAC1 <<       ( CVR_VEG_ms 込み )
!   -----------------------------------------------
!
        FAC1    = CVR_VEG_ms(I,1) * ( 1.0D0 - THERMK(I) )
!
!   ----------------------------------------
!  >> CG それぞれが吸収する全放射量 RADT <<
!   ----------------------------------------
!
!      CVR_VEG_ms を乗じた単位面積当たりであることに注意。
!
!      短波の吸収量 入力 RADT に加えて、
!      キャノピーのとり分は
!          ( 地面放射の吸収 - キャノピー放射 ) + 長波吸収分
!      地面のとり分は
!          ( キャノピー放射の吸収 - 地面放射 ) + キャノピーで吸収し損ねた長波
!
!      注意。雪では、短波入射 RAD_SHRT(I,2) は SUM には加算しない。
!
        RAD_A2C (I) = RAD_LONG_INTO_SIB(I)  *  FAC1
        RAD_A2G (I) = RAD_LONG_INTO_SIB(I)  *  (1-FAC1)
        RAD_C2A (I) = - TC4 * FAC1      * STEFAN_c
        RAD_C2G (I) =   TC4 * FAC1      * STEFAN_c
        RAD_G2A (I) = - TG4 * (1-FAC1)  * STEFAN_c
        RAD_G2C (I) = - TG4 * FAC1      * STEFAN_c
!
        RAD_LONG_SIB_2_ATM(I) = RAD_C2A(I) + RAD_G2A(I) 
!
        RAD_LONG(I,1) = 
     1             + RAD_A2C(I) + RAD_C2A(I)
     1             - RAD_G2C(I) - RAD_C2G(I) 
        RAD_LONG(I,2) = 
     1             + RAD_A2G(I) + RAD_G2A(I)
     1             + RAD_C2G(I) + RAD_G2C(I) 
!
 1000 CONTINUE
!
! 雪層で吸収する短波放射 ... 雪あり域のみ
!
      IF ( INSFLAG.EQ.1 ) THEN                         ! 雪あり
        DO 2100 I=NSTR,NEND
          FACT = ( TMP_FREZ_c - TMP_GRD(I) ) * DELT_TEMP_INV 
          FACT = MAX ( MIN ( FACT, ONE )  , ZERO )
!    
          XMU =   XMU_TF * (1-FACT) +  XMU_GLA * FACT  
!
          TRANS =  FACT_TRANS_SR_SNOW_SKIN_c * DEXP ( - XMU * DEP ) 
!         TRANS = 0.                                   ! とりあえず
!
          RAD_SHRT_INTO_SNOW_LAYER (I) = RAD_SHRT(I,2) * TRANS  
          RAD_SHRT  (I,2)              = RAD_SHRT(I,2) * (1-TRANS)  
 2100   CONTINUE  
      ENDIF  
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
!
!  下向き長波
!
      IF ( 
     1          IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_BUDJET
     1    ) THEN   
        CMARK = 'DLWBL'                               ! 旧 HLRDTOP 
        CALL MONIT_SIB ( CMARK , RAD_LONG_INTO_SIB , 
     1                   IFLG_MON , DELT_CNP, ONE )
      ENDIF
!
      IF ( 
     1          IMONIT_LEVEL .GE. IMONIT_ALMA 
     1    ) THEN   
        CMARK = 'A_LWDN'                               ! 旧 HLRDTOP 
        CALL MONIT_SIB ( CMARK , RAD_LONG_INTO_SIB , 
     1                   IFLG_MON , DELT_CNP, ONE )
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ------------
! > 雪無し域 <
! ------------
!
!  短波吸収  
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'HSRN' 
          CALL MONIT_SIB ( CMARK , RAD_SHRT     ,ICN, DELT_CNP, ONE )
!
          CMARK = 'HSRG' 
          CALL MONIT_SIB ( CMARK , RAD_SHRT(1,2),ICN, DELT_CNP, ONE )
!
!  長波吸収 ( explicit ) 
!
          CMARK = 'HLRN0' 
          CALL MONIT_SIB ( CMARK , RAD_LONG     ,ICN, DELT_CNP, ONE )
!
          CMARK = 'HLRG0' 
          CALL MONIT_SIB ( CMARK , RAD_LONG(1,2),ICN, DELT_CNP, ONE )
!
!  長波成分別フラックス
!
          CMARK = 'HLRAN'            ! 大気 -> キャノピー          
          CALL MONIT_SIB ( CMARK , RAD_A2C,ICN, DELT_CNP, ONE )
          CMARK = 'HLRNA0'           ! キャノピー -> 大気         
          CALL MONIT_SIB ( CMARK , RAD_C2A,ICN, DELT_CNP, ONE )
!
          CMARK = 'HLRNG0'           ! キャノピー -> 地面          
          CALL MONIT_SIB ( CMARK , RAD_C2G,ICN, DELT_CNP, ONE )
          CMARK = 'HLRGN0'           ! 地面 -> キャノピー          
          CALL MONIT_SIB ( CMARK , RAD_G2C,ICN, DELT_CNP, ONE )
!
          CMARK = 'HLRAG'            ! 大気 -> 地面
          CALL MONIT_SIB ( CMARK , RAD_A2G,ICN, DELT_CNP, ONE )
          CMARK = 'HLRGA0'           ! 地面 -> 大気 
          CALL MONIT_SIB ( CMARK , RAD_G2A,ICN, DELT_CNP, ONE )
!
        ELSE 
!
! --------------
! > 雪ありいき <
! --------------
!
!  短波吸収
!
        CMARK = 'HSRS' 
        CALL MONIT_SIB ( CMARK , RAD_SHRT     ,ICS, DELT_CNP, ONE )
!
        CMARK = 'HSRWTOP' 
        CALL MONIT_SIB ( CMARK , RAD_SHRT(1,2),ICS, DELT_CNP, ONE )
        CALL MONIT_SIB ( CMARK , RAD_SHRT_INTO_SNOW_LAYER, 
     1                                         ICS, DELT_CNP, ONE )
!
        CMARK = 'HSRWSK' 
        CALL MONIT_SIB ( CMARK , RAD_SHRT(1,2),ICS, DELT_CNP, ONE )
!
        CMARK  = 'HSRBSW'                                         
        CALL MONIT_SIB ( CMARK , RAD_SHRT_INTO_SNOW_LAYER, 
     1                   ICS, DELT_CNP, ONE)
!
!  長波吸収
!
          CMARK = 'HLRS0' 
          CALL MONIT_SIB ( CMARK , RAD_LONG     ,ICS, DELT_CNP, ONE )
!
          CMARK = 'HLRW0' 
          CALL MONIT_SIB ( CMARK , RAD_LONG(1,2),ICS, DELT_CNP, ONE )
!
!  長波成分別フラックス
!
          CMARK = 'HLRAS' 
          CALL MONIT_SIB ( CMARK , RAD_A2C,ICS, DELT_CNP, ONE )
          CMARK = 'HLRSA0' 
          CALL MONIT_SIB ( CMARK , RAD_C2A,ICS, DELT_CNP, ONE )
!
          CMARK = 'HLRSW0' 
          CALL MONIT_SIB ( CMARK , RAD_C2G,ICS, DELT_CNP, ONE )
          CMARK = 'HLRWS0' 
          CALL MONIT_SIB ( CMARK , RAD_G2C,ICS, DELT_CNP, ONE )
!
          CMARK = 'HLRAW' 
          CALL MONIT_SIB ( CMARK , RAD_A2G,ICS, DELT_CNP, ONE )
          CMARK = 'HLRWA0' 
          CALL MONIT_SIB ( CMARK , RAD_G2A,ICS, DELT_CNP, ONE )
!
        ENDIF
!
      ENDIF 
#endif
!
      RETURN
      END
!===========================================================
      SUBROUTINE MR_RAD_SUM_ALL
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1    IMONIT_LEVEL  , 
     1    IMONIT_ALL    ,
     1    IMONIT_BUDJET   ,
     1    IMONIT_FLUX_MRI , IMONIT_ALMA 
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
        WRITE(6,*) 'MONIT_REGIST_RAD_SUM_ALL 99/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! 長波トップ ( たまたま格子平均量 )
!
#ifdef MONYOS
      IF ( 
     1          IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_BUDJET
     1    ) THEN   
!       CMARK  = 'HLRDTOP'
        CMARK  = 'DLWBL'
        CTITLE = 'RAD LONG DOWNWARD AT SIB TOP'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      IF ( 
     1          IMONIT_LEVEL .GE. IMONIT_ALMA
     1    ) THEN   
        CMARK  = 'A_LWDN'
        CTITLE = 'ALMA + @LWnet longwave downward '
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF

!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 短波吸収
!
      CMARK  = 'HSRN'
      CTITLE = 'RAD SHORT INTO CANOPY NOS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSRS'
      CTITLE = 'RAD SHORT INTO CANOPY SNW'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSRG'
      CTITLE = 'RAD SHORT INTO GROUND'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSRWTOP'
      CTITLE = 'RAD SHORT INTO SNOW TOP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSRWSK'                                         
      CTITLE = 'SHORT RAD SNOW SKIN'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HSRBSW'                                         
      CTITLE = 'SHORT RAD INTO SNOW LAYER'
      CUNIT  = 'W/M**2'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!
! 長波吸収
!
      CMARK  = 'HLRN0'
      CTITLE = 'RAD LONG INTO CANOPY NOS explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRS0'
      CTITLE = 'RAD LONG INTO CANOPY SNW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRG0'
      CTITLE = 'RAD LONG INTO GROUND explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRW0'
      CTITLE = 'RAD LONG INTO SNOW explicit'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 長波要素別
!
      CMARK  = 'HLRAN' 
      CTITLE = 'LONG RAD FROM ATM TO CNP NOS' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRNA0' 
      CTITLE = 'LONG RAD FROM CNP TO ATM NOS explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRAG' 
      CTITLE = 'LONG RAD FROM ATM TO GROUND' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRGA0' 
      CTITLE = 'LONG RAD FROM GROUND TO ATM explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRNG0' 
      CTITLE = 'LONG RAD FROM CNP NOS TO GROUND explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRGN0' 
      CTITLE = 'LONG RAD FROM GROUND TO CNP NOS explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRAS' 
      CTITLE = 'LONG RAD FROM ATM TO CNP SNW' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRSA0' 
      CTITLE = 'LONG RAD FROM CNP TO ATM SNW explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRAW' 
      CTITLE = 'LONG RAD FROM ATM TO SNOW' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HLRWA0' 
      CTITLE = 'LONG RAD FROM SNOW TO ATM explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRSW0' 
      CTITLE = 'LONG RAD FROM CNP SNW TO SNOW explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK = 'HLRWS0' 
      CTITLE = 'LONG RAD FROM SNOW TO CNP SNW explicit' 
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

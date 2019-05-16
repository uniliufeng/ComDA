! 雪の熱フラックス計算準備  （比熱計算、短波放射吸収量計算）
! 99/07/12 保坂征宏  雪の比熱計算を hcap.F から分離する形で作成
! 99/07/xx 保坂征宏  雪の短波吸収量計算追加
! 99/08/04 保坂征宏  大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SNOW_HFLUX_PRE (
     I    NMAX      , NMAX_01   , NMAX_1  , 
     I    NMAX_2    , NMAX_3    , NMAX_4  , 
     I    DELT_SNOW ,  DELT_SNOW_INV         ,  
     I    WTR_SNOW  , GLA_SNOW  , 
     I    SRAD_BUCKET_SNW       , 
     I    ENG_BUCKET_SNW        , 
     I    DEP_SNOW, TMP_SNOW, 
     I    CVR_SNW        , CVR_SNW_INV    , 
     I    HDF_SNOW_H , 
     O    PHASE_SNOW,  
     O    CAP_SNOW  , 
     O    HFLUX_CONV , SRAD_BUCKET_SNOW_TO_SNSL , 
     O    HDIFF_FLUX_H   )
!
! 注意
!   SRAD_BUCKET_SNW , ENG_BUCKET_SNW は格子平均値である。
!   雪での熱の計算では、雪あり域だけでの平均値を用いる。
!   _01 の格子ではずれが生じるので、このサブルーチンの中で
!       入力時には /CVR_SNW を施す必要がある。
!   
!   なお、_COR がついている変数はそれに /CVR_SNW を施したものである。
!
!   注意
!     ・上からのフラックスは与えられている。      
!
!     ・雪の層数と、TMP_SNSL 等の関係
!       雪 1 層の時    TMP_SNOW(I,2)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,1) = 雪下端と雪下土壌
!       雪 2 層の時    TMP_SNOW(I,3)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,2) = 雪下端と雪下土壌
!       雪 3 層の時    TMP_SNOW(I,4)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,3) = 雪下端と雪下土壌
! 
!     ・T_{frez} に等しい場合はその層は温度変化しないとして評価する。
!       PHASE_SNOW   ! 相変化 on 0 off 1
!
! メモ
!   PHASE_SNOW で、層の数の一つ下の層で零としているのは、フラックス計算の
!   際の便宜のため
!
      USE SIBPRM , ONLY : 
     I    IDX     , 
     I    ISN
!
      USE SIBCON , ONLY : 
     1    CAP_WTR_c  ,                            ! CW 
     1    CAP_GLA_c  ,                            ! CICE 
     1    TMP_FREZ_c
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICS , INS ,
     1                IMONIT_LEVEL , IMONIT_ALL , IMONIT_FLUX_COMMON
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX  
      INTEGER,INTENT(IN)    :: NMAX_01(3)  
      INTEGER,INTENT(IN)    :: NMAX_1(3)    
      INTEGER,INTENT(IN)    :: NMAX_2(3)  
      INTEGER,INTENT(IN)    :: NMAX_3(3)    
      INTEGER,INTENT(IN)    :: NMAX_4(3)  
      REAL(8),INTENT(IN)    :: DELT_SNOW
      REAL(8),INTENT(IN)    :: DELT_SNOW_INV
!
      REAL(8),INTENT(IN)    :: WTR_SNOW  (IDX,ISN)
      REAL(8),INTENT(IN)    :: GLA_SNOW  (IDX,ISN)
      REAL(8),INTENT(IN)    :: CVR_SNW     (IDX)
      REAL(8),INTENT(IN)    :: CVR_SNW_INV   (IDX)
!
      REAL(8),INTENT(IN)    :: SRAD_BUCKET_SNW (IDX)     ! 上端からの短波入射
      REAL(8),INTENT(IN)    :: ENG_BUCKET_SNW  (IDX)      ! 上端からの熱
      REAL(8),INTENT(IN)    :: DEP_SNOW      (IDX,ISN+1) ! 雪の深さ
      REAL(8),INTENT(IN)    :: TMP_SNOW      (IDX,ISN+1) ! 雪の温度
!
      REAL(8),INTENT(IN)    :: HDF_SNOW_H  (IDX,ISN)   ! 雪中熱伝導係数
!     REAL(8),INTENT(IN)    :: HDF_SNOW_BH (IDX)       ! 雪B と地面
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: PHASE_SNOW    (IDX,ISN) ! 零度 -> 0 それ以下->1
!                                                     ! ただし雪がない層では0
!
      REAL(8),INTENT(OUT)  :: CAP_SNOW      (IDX,ISN) ! 雪の比熱
      REAL(8),INTENT(OUT)  :: HFLUX_CONV (IDX,ISN) ! 雪の短波吸収
      REAL(8),INTENT(OUT)  :: SRAD_BUCKET_SNOW_TO_SNSL(IDX)  
!                                                 ! 地面第一層への短波透過
!                                                     !   単位時間当
!
!     REAL(8),INTENT(OUT)  :: ENG_BUC (IDX)       ! 上端からのエネルギー
!                                                 ! 単位時間当、雪有域の平均
!
      REAL(8),INTENT(OUT)  :: HDIFF_FLUX_H   (IDX,0:ISN) 
!                                            ! 層の数と等しい K では、
!                                            ! 地中へのフラックス
!
! WORK
      INTEGER   :: I
      REAL(8)   :: XMU
      REAL(8)   :: ABSORB
      REAL(8)   :: SRAD_FLUX_H   (IDX,0:1) ! 雪の短波フラックス 0->1,1->2
!                                                     !   単位時間当
      REAL(8)   :: SRAD_BUC_DOWN (IDX)     ! 地面第一層への短波透過
!                                                     !   単位時間当
!
      REAL(8),PARAMETER   :: XMU_TF   =  20.D0      ! tmp_frez 
      REAL(8),PARAMETER   :: XMU_GLA  = 100.D0      ! < tmp_frez_c 
      REAL(8),PARAMETER   :: DELT_TEMP_INV = 0.5D0  ! = 1/2K
      REAL(8),PARAMETER   :: ONE = 1.D0 
      REAL(8)             :: FACT 
!
#ifdef MONYOS
!     REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
      
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_HFLUX_PRE POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
! --------------
! > PHASE_SNOW <
! --------------
!
      DO 1100 I=NMAX_01(2),NMAX
        IF ( TMP_SNOW(I,1) .EQ. TMP_FREZ_c ) THEN
          PHASE_SNOW(I,1) = 0
        ELSE
          PHASE_SNOW(I,1) = 1
        ENDIF
 1100 CONTINUE 
      DO 1200 I=NMAX_2(2),NMAX
        IF ( TMP_SNOW(I,2) .EQ. TMP_FREZ_c ) THEN
          PHASE_SNOW(I,2) = 0
        ELSE
          PHASE_SNOW(I,2) = 1
        ENDIF
 1200 CONTINUE 
      DO 1300 I=NMAX_3(2),NMAX
        IF ( TMP_SNOW(I,3) .EQ. TMP_FREZ_c ) THEN
          PHASE_SNOW(I,3) = 0
        ELSE
          PHASE_SNOW(I,3) = 1
        ENDIF
 1300 CONTINUE 
      DO 1400 I=NMAX_4(2),NMAX
        IF ( TMP_SNOW(I,4) .EQ. TMP_FREZ_c ) THEN
          PHASE_SNOW(I,4) = 0
        ELSE
          PHASE_SNOW(I,4) = 1
        ENDIF
 1400 CONTINUE 
!
!  以下はダミー。フラックス計算での便宜のため
!
      DO 1500 I=NMAX_01(2),NMAX_1(3)
        PHASE_SNOW(I,2) = 0
 1500 CONTINUE 
      DO 1600 I=NMAX_2(2),NMAX_2(3)
        PHASE_SNOW(I,3) = 0
 1600 CONTINUE 
      DO 1700 I=NMAX_3(2),NMAX_3(3)
        PHASE_SNOW(I,4) = 0
 1700 CONTINUE 
!
      HFLUX_CONV (:,:) = 0. 
!
! ------------
! > 雪第一層 <
! ------------
!
      DO 2100 I=NMAX_01(2),NMAX           ! 雪があるところ
!
!   熱容量   
!
        CAP_SNOW(I,1) =
     1      (   CAP_GLA_c  * GLA_SNOW(I,1)
     1        + CAP_WTR_c  * WTR_SNOW(I,1) )
     1       * CVR_SNW_INV(I)                 ! fix 01/08/17 格子 B 平均に
!                                         
!
!   上端から入射するエネルギーをフラックスに換算 （雪あり域平均）
!
        HDIFF_FLUX_H(I,0) = ENG_BUCKET_SNW (I) 
     1                       * CVR_SNW_INV(I) * DELT_SNOW_INV 
!
        SRAD_FLUX_H(I,0) = SRAD_BUCKET_SNW(I) 
     1                       * CVR_SNW_INV(I) * DELT_SNOW_INV 
!
!   短波放射フラックス吸収率 
!
        FACT = ( TMP_FREZ_c - TMP_SNOW(I,1) ) * DELT_TEMP_INV 
        FACT = MIN ( FACT, ONE )  
        XMU =   XMU_TF * (1-FACT) +  XMU_GLA * FACT  
!
        ABSORB = 1. - DEXP( -XMU * DEP_SNOW(I,1) )        
        SRAD_FLUX_H(I,1) = SRAD_FLUX_H(I,0) * ( 1 - ABSORB ) ! 一層の下へ
        SRAD_BUC_DOWN (I) =  SRAD_FLUX_H(I,1)
        HFLUX_CONV (I,1) = SRAD_FLUX_H(I,0) - SRAD_FLUX_H(I,1)
!
!   熱伝導フラックス収束  （一層しかない場合は TMP_SNOW(2) = TMP_SNSL(I)）
!
        HDIFF_FLUX_H (I,1) = 
     1        HDF_SNOW_H(I,1) * ( TMP_SNOW (I,1) - TMP_SNOW(I,2) ) 
!
#ifdef MONYOS
        HDIFF_FLUX_H (I,2) = HDIFF_FLUX_H (I,1)  
        HDIFF_FLUX_H (I,3) = HDIFF_FLUX_H (I,1)  
        HDIFF_FLUX_H (I,4) = HDIFF_FLUX_H (I,1)  
#endif
!
 2100 CONTINUE
!
! ------------
! > 雪第二層 <
! ------------
!
      DO 2200 I=NMAX_2(2),NMAX         
!
!   熱容量
!
        CAP_SNOW(I,2) 
     1     =   CAP_GLA_c  * GLA_SNOW(I,2)
     1       + CAP_WTR_c  * WTR_SNOW(I,2) 
!
!   上端から入射するエネルギー = 全て吸収
!
        SRAD_BUC_DOWN (I)   = 0. 
        HFLUX_CONV (I,2) = SRAD_FLUX_H(I,1) 
!
!   熱伝導フラックス収束
!
        HDIFF_FLUX_H (I,2) = 
     1        HDF_SNOW_H(I,2) * ( TMP_SNOW (I,2) - TMP_SNOW(I,3) ) 
!
#ifdef MONYOS
        HDIFF_FLUX_H (I,3) = HDIFF_FLUX_H (I,2)  
        HDIFF_FLUX_H (I,4) = HDIFF_FLUX_H (I,2)  
#endif
 2200 CONTINUE
!
! ------------
! > 雪第三層 <
! ------------
!
      DO 2300 I=NMAX_3(2),NMAX            ! 雪が 3 層以上あるところ
        CAP_SNOW(I,3) 
     1     =   CAP_GLA_c  * GLA_SNOW(I,3)
     1       + CAP_WTR_c  * WTR_SNOW(I,3) 
!
!   熱伝導フラックス収束
!
        HDIFF_FLUX_H (I,3) = 
     1        HDF_SNOW_H(I,3) * ( TMP_SNOW (I,3) - TMP_SNOW(I,4) ) 
!
#ifdef MONYOS
        HDIFF_FLUX_H (I,4) = HDIFF_FLUX_H (I,3)  
#endif
 2300 CONTINUE
!
! ------------
! > 雪第四層 <
! ------------
!
      DO 2400 I=NMAX_4(2),NMAX            ! 雪が 4 層以上あるところ
        CAP_SNOW(I,4) 
     1     =   CAP_GLA_c  * GLA_SNOW(I,4)
     1       + CAP_WTR_c  * WTR_SNOW(I,4) 
!
!   熱伝導フラックス収束
!
        HDIFF_FLUX_H (I,4) = 
     1        HDF_SNOW_H(I,4) * ( TMP_SNOW (I,4) - TMP_SNOW(I,5) ) 
 2400 CONTINUE
!
!   短波の地中への浸透量を、BUCKET に（雪あり・無し域平均、* DELT_SNOW 処理）
!     ( 実際には一層の場合のみ値あり )
!
      DO 3000 I=NMAX_01(2),NMAX
        SRAD_BUCKET_SNOW_TO_SNSL(I) 
     1          = SRAD_BUC_DOWN(I) * DELT_SNOW * CVR_SNW(I)
 3000 CONTINUE
!
#ifdef MONYOS
!
! ==============
! > 格子平均量 <
! ==============
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!  上端でのエネルギー
!
        CMARK  = 'HFLBWS'
        CALL MONIT_SIB ( CMARK , ENG_BUCKET_SNW , INS , ONE , ONE )  
!       CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,0),ICS,DELT_SNOW, ONE )   
!       以上二つは同じはず。
!
!  上端での短波
!     上端
        CMARK  = 'HSRBWS'
        CALL MONIT_SIB ( CMARK , SRAD_BUCKET_SNW , INS , ONE , ONE )  
!
      ENDIF
!
! ========
! > 詳細 <
! ========
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!  下端での短波 = いずれ消す 
!
        CMARK  = 'HSRBWL'
        CALL MONIT_SIB ( CMARK , SRAD_BUCKET_SNOW_TO_SNSL,INS,ONE,ONE)  
      
!
!  熱フラックス （陽的評価）
!
        CMARK  = 'HFLW12' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,1),ICS,DELT_SNOW, ONE )   
        CMARK  = 'HFLW23' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,2),ICS,DELT_SNOW, ONE )  
        CMARK  = 'HFLW34' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,3),ICS,DELT_SNOW, ONE )  
        CMARK  = 'HFLW4L' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,4),ICS,DELT_SNOW, ONE )  
!
      ENDIF
!
!
! 短波フラックス ... いずれ消去
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
        CMARK  = 'HSRW12' 
        CALL MONIT_SIB ( CMARK , SRAD_FLUX_H(1,1),ICS,DELT_SNOW, ONE )  
      ENDIF
#endif
!
      RETURN
      END
!==================================================================
      SUBROUTINE MR_SNOW_HFLUX_PRE
!
! SNOW_PRE で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                IMONIT_LEVEL , IMONIT_ALL , IMONIT_FLUX_COMMON
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
        WRITE(6,*) 'MONIT_REGIST_SNOW_HFLUX_PRE 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! ------------------
! > 雪全体の収支用 <
! ------------------
!
!  エネルギーフラックス
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
      CMARK  = 'HFLBWS'
      CTITLE = 'ENG FLUX FROM CANOPY TO SNOW TOP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  短波放射フラックス
!
      CMARK  = 'HSRBWS'
      CTITLE = 'SRAD FLUX INTO SNOW TOP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'HSRBWL'
      CTITLE = 'SRAD FLUX INTO SNSL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 熱伝導 ( 陽的評価のみ ) 
!
      CMARK  = 'HFLW12'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 1 TO 2 EXP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW23'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 2 TO 3 EXP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW34'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 3 TO 4 EXP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW4L'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 4 TO SNSL EXP'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
! 熱伝導 ( 全体 ) 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
      CMARK  = 'HFLW12F'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 1 TO 2 '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW23F'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 2 TO 3 '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW34F'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 3 TO 4 '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HFLW4LF'
      CTITLE = 'HEAT DIFFUSION FLUX FROM SNOW 4 TO SNSL '
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
!
! 短波放射フラックス
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'HSRW12'
      CTITLE = 'SRAD FLUX FROM SNOW 1 TO 2'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 熱バケツ （雪 -> 土壌） ... snow_hflux_imp/explicit
!
      CMARK  = 'HFLBWL'
      CTITLE = 'HEAT FROM SNOW TO SOIL'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

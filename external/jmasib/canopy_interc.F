! インターセプション
! 99/06/30 保坂征宏
! 99/08/10 保坂征宏  変数名変更
!                    出力要素変更
!                    コメント追加
!                    雨・雪の温度を大気最下層の xx 温度にする
!
      SUBROUTINE CANOPY_INTERC (
     I   NSTR           , NEND           , INSFLAG        , 
     I   NTYP           , 
#ifdef MONYOS
     I   DELT_CNP       ,
#endif
     I   RAIN_LARGE     , RAIN_CONV      , 
     I   CVR_VEG_msx    , CVR_S_or_N     , XLAI_VEG_msx   , 
     I   TMP_REF        , 
     I   TMP_CNP_IN     , TMP_GRSK_IN    , 
     I   FR_WTR_CNP_IN  , FR_GLA_CNP_IN  , 
     I   FR_WTR_GRS_IN  , FR_GLA_GRS_IN  , 
     O   TMP_CNP        , TMP_GRSK       , 
     O   FR_H2O_CNP     , FR_WBH_CNP     , 
     O   FR_H2O_GRS     , FR_WBH_GRS     , 
     O   H2O_SAT_CNP    , H2O_SAT_GRS    , 
     M   H2O_BUCKET_INTO  , ENH_BUCKET_INTO    ) 
!
! -----------------------------------------------------
! ◎canopy_interc の解説
!
!  ○概要 : 降水捕捉の計算
!    大気モデルで計算された降水（H2O とエネルギーをもつ）を
!    キャノピー・下草で捕捉する。捕捉しきれない分は地面に到達する。
!    
!  ○入力される降水量 （大規模凝結 RAIN_LARGE と積雲性降水 RAIN_CONV ） 
!    は、現ステップを 0->Δt として、-2Δt->0 のステップでの降水量。
!    単位mm/m**2 。
!    注意
!      この単位で降水の分布は決定している？（a_pr 等）
!      それ以外の目的ではこの単位は使わない。kg/m**2 を使用。
!      変数は PRECIP である。 
!
! -----------------------------------------------------
!
! インターセプションの計算。
!   インターセプション後のキャノピー・下草の状況計算も行う（比熱計算も含む）。
!
! 注意。
!   ・葉面積指数 ZLT 、保水可能水量、保水量等は
!     植生被覆率で割っていないことに注意。
!   ・保水可能量を越える水があることは無いと仮定。 
!     そして、降水（・降雪）により保水可能量を越えた場合、
!     葉にもともとあった水は落ちず、降水してきた水が下に落ちるとする。
!     キャノピーや下草にもとからあった水が混ざらないため、
!     下草や裸地への降水の温度は、もとの降水の温度と同じである。
!
! 雨・雪の判別とその後
!
!   降ってきた状況では雨または雪のいずれか。
!   大気最下層温度で決める。
!       --- これは xx 温度にする。
!
!   捕捉され損ねた水の温度は、降水の温度と同じとする。
!   つまり
!     ・もともとあったものが落ちるのではなく、降ってきたものが落ちる。）
!     ・キャノピー・下草で一旦捕捉されその後落下した水も、
!       キャノピー・下草の水とは馴染まず、落下時の温度をたもつ。
!     ・したがって捕捉され損ねた水の温度なども、降雪時と同じ
!   とする。
!
! 注意 : B 平均と C平均について
!   B 平均とは、雪あり or なし域全体での平均
!   C 平均とは、植生がある部分 ( CVR_VEG でカウントされる部分 ) のみでの平均
!   降水の遮断では、C 平均で計算し、B 平均にもどすべし
!   しかし、従来 SiB ではそれをやっていない。
!   NEWSIB ではそうするつもりだが、とりあえず現状は保留し、
!   #ifdef CORRECT #endif に書くにとどめる。
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ICG    , 
     1    ITYP
!
      USE SIBCON , ONLY :
     1    CAP_GLA_c ,
     1    CAP_WTR_c ,
     1    CAP_LAI_c , 
     1    TMP_FREZ_c  , 
     1    RHO_WTR_c   ,
     1    HL_W2G_c    , 
     1    CAP_GRSK_c  ,
     1    H2O_SAT_CNP_DEP_c  ,
     1    H2O_SAT_GRS_DEP_c  , 
     1    FR_PREITC_c ,
     1    B_PR_COEF_c   ,
     1    A_PC_COEF_c   ,
     1    A_PL_COEF_c   ,
     1    C_PC_COEF_c   ,
     1    C_PL_COEF_c  

!
      USE SIBVEG , ONLY :
     1           CHI_L_co                  ! χL  (ITYP,ICG) 
!
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS , INA,
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , 
     1     IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR 
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG
!
#ifdef MONYOS
      REAL(8),INTENT(IN)    :: DELT_CNP
#endif
!
      INTEGER,INTENT(IN)    :: NTYP        (IDX) 
      REAL(8),INTENT(IN)    :: CVR_VEG_msx (IDX,ICG) ! 被覆率
      REAL(8),INTENT(IN)    :: CVR_S_or_N  (IDX)     ! 雪被覆率 or 1-被覆率
      REAL(8),INTENT(IN)    :: TMP_REF (IDX) ! 大気最下層温度
      REAL(8),INTENT(IN)    :: XLAI_VEG_msx (IDX,ICG) 
!
      REAL(8),INTENT(IN)    :: RAIN_LARGE    (IDX)  ! 大規模降水降雪量 mm/m**2 
      REAL(8),INTENT(IN)    :: RAIN_CONV     (IDX)  ! 積雲性降水降雪量 mm/m**2 
!
!   予報変数
!
      REAL(8),INTENT(IN)    :: TMP_CNP_IN    (IDX)  ! キャノピー温度
      REAL(8),INTENT(IN)    :: FR_WTR_CNP_IN (IDX)  ! キャノピー水割合
      REAL(8),INTENT(IN)    :: FR_GLA_CNP_IN (IDX)  ! キャノピー氷割合
      REAL(8),INTENT(IN)    :: TMP_GRSK_IN   (IDX)  ! 地面・下草温度
      REAL(8),INTENT(IN)    :: FR_WTR_GRS_IN (IDX)  ! 下草水割合
      REAL(8),INTENT(IN)    :: FR_GLA_GRS_IN (IDX)  ! 下草氷割合
!
!
! OUT
!
!   予報変数 （インターセプションによる変更後）
!
      REAL(8),INTENT(OUT)   :: TMP_CNP    (IDX)   ! 
      REAL(8),INTENT(OUT)   :: TMP_GRSK   (IDX)   ! 
      REAL(8),INTENT(OUT)   :: FR_H2O_CNP (IDX)   ! 
      REAL(8),INTENT(OUT)   :: FR_H2O_GRS (IDX)   ! 
!
!   可能保水量
!
      REAL(8),INTENT(OUT)   :: H2O_SAT_CNP (IDX) 
      REAL(8),INTENT(OUT)   :: H2O_SAT_GRS (IDX) 
!
!
!   水の占める割合
!
      REAL(8),INTENT(OUT)   :: FR_WBH_CNP (IDX)   ! 
      REAL(8),INTENT(OUT)   :: FR_WBH_GRS (IDX)   ! 
!
! INOUT
!
!   バケツ（雪あり/無しを含む格子平均値）
!
      REAL(8),INTENT(INOUT) :: H2O_BUCKET_INTO          (IDX) 
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_INTO          (IDX) 

!
! WORK
!
      REAL(8) :: PRECIP     (IDX)   ! 降水降雪量 ... 上端入射（単位mm/m**2）
      REAL(8) :: PRECIP_W   (IDX)   ! 降水降雪量 ... 上端入射（単位mm/m**2）
      REAL(8) :: PRECIP_G   (IDX)   ! 降水降雪量 ... 上端入射（単位mm/m**2）
!
      REAL(8) :: A_PR_COEF  (IDX)   ! 降水中局在する
      REAL(8) :: C_PR_COEF  (IDX)
!
      REAL(8) :: TMP_PR     (IDX)   ! 降水の温度 
!
      REAL(8) :: PR_DROP_CNP (IDX)  ! キャで一旦補足されるも落下する量 
      REAL(8) :: PR_DROP_GRS (IDX)  ! 下草で一旦補足されるも落下する量 
      REAL(8) :: PR_TRAP_CNP (IDX)  ! 
      REAL(8) :: PR_TRAP_GRS (IDX)  ! 
      REAL(8) :: PR_TRAP_MAX (IDX)  ! 
      REAL(8) :: FR_TRAP_CNP (IDX)  ! 
      REAL(8) :: FR_TRAP_GRS (IDX)  ! 
!
      REAL(8) :: CAP_CNP  
      REAL(8) :: CAP_GRSK  
!
!
      REAL(8) :: CHIV
!
      REAL(8) :: AA
      REAL(8) :: BB
      REAL(8) :: EXTINC
!
      REAL(8) :: ENG_CNP(IDX)
      REAL(8) :: ENX    (IDX)
      REAL(8) :: ENG_GRSK(IDX)  
      REAL(8) :: ENG_DIFF           ! ENG - ENX  
!
      REAL(8) :: FR_WBH_PR  (IDX)   ! 降水・雪中の雪の割合
!
      REAL(8) :: H2O_CNP (IDX)
      REAL(8) :: H2O_GRS (IDX)
!
      REAL(8) :: H2O_SAT_CNP_INV(IDX)
      REAL(8) :: H2O_SAT_GRS_INV(IDX)
!
      REAL(8) :: ENG_TRAP_CNP (IDX)  
      REAL(8) :: ENG_TRAP_GRS (IDX)  
!
      REAL(8) :: WTR_CNP
      REAL(8) :: GLA_CNP
      REAL(8) :: WTR_GRS  
      REAL(8) :: GLA_GRS  
!
      REAL(8),SAVE :: PRECIP_MIN         ! 降水量下限
      REAL(8),SAVE :: CHIL_MIN           ! 降水量下限
      REAL(8),SAVE :: ENG_DIFF_MIN     
      REAL(8),SAVE :: H2O_SAT_CNP_MIN 
      REAL(8),SAVE :: H2O_SAT_GRS_MIN 
      REAL(8),SAVE :: DZERO
! 
      INTEGER :: I
!
      LOGICAL,SAVE  :: LFIRST 
!
#ifdef MONYOS
      CHARACTER(7) :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
      REAL(8)      :: WORK(IDX) 
#endif
!
      DATA PRECIP_MIN  / 1.0D-08 /         ! 単位 mm/step = kg/step 
      DATA CHIL_MIN    / 1.D-2   /
      DATA DZERO       / 0.D0    /
      DATA ENG_DIFF_MIN / 1.D-6 /  
      DATA H2O_SAT_CNP_MIN / 1.D-10 /  
      DATA H2O_SAT_GRS_MIN / 1.D-10 /  
      DATA LFIRST       / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_INTERC POCHAKA 99/08' 
        LFIRST = .FALSE. 
      ENDIF
!
!  ============
!  >> 前処理 <<
!  ============
!    
      DO 100 I= NSTR,NEND
!
!  ----------------------------
!  > 降水・降雪の状況チェック <
!  ----------------------------
!
! 降水の温度                  ... これはいずれ xx 温度に変える
!
        TMP_PR (I) = TMP_REF(I)   
!
! 降水の内の水の割合
!
        IF ( TMP_PR(I) .LT. TMP_FREZ_c ) THEN   ! 零度より低いと氷
          FR_WBH_PR (I) = 0. 
        ELSE                                    ! たまたま零度なら水
          FR_WBH_PR (I) = 1.
        ENDIF
 100  CONTINUE
!
!  ==========================================
!  >> キャノピーインターセプションここから << 
!  ==========================================
!
! ---------------
! >  保水量上限 <
! ---------------
!  これに割合 FR_WTR_CNP or FR_GLA_CNP をかけると水、氷の B 平均になる
!
      DO 900 I=NSTR,NEND 
        H2O_SAT_CNP(I)              ! 単位は kg/m**2 （格子B 平均）
     1      = H2O_SAT_CNP_DEP_c     ! 水換算 0.1mm 
     1        * RHO_WTR_c           ! 1.D3 
     1        * XLAI_VEG_msx(I,1)   ! これで単位面積当たり質量に 
!         
#ifdef SIB_DEBUG
        if ( H2O_SAT_CNP(I) .lt. H2O_SAT_CNP_MIN ) then   
          write(6,*) 'canopy_interc check --- sat_cnp ' , i , 
     1               H2O_SAT_CNP(I) .lt. H2O_SAT_CNP_MIN
        endif
#endif
!
        H2O_SAT_CNP_INV(I) = 1/ MAX( H2O_SAT_CNP(I), H2O_SAT_CNP_MIN )    
!
        H2O_CNP(I) = ( FR_WTR_CNP_IN(I) + FR_GLA_CNP_IN(I) ) 
     1               * H2O_SAT_CNP(I)
!
#ifdef SIB_DEBUG
        IF (  ( FR_WTR_CNP_IN(I) + FR_GLA_CNP_IN(I) ).GT.1.+1.D-8) THEN
          write(6,*) 'canopy_interc error fr_wtr_cnp gt 1 ' , 
     1                 I                , 
     1                 FR_WTR_CNP_IN(I) , FR_GLA_CNP_IN(I) , 
     1                 FR_WTR_CNP_IN(I) + FR_GLA_CNP_IN(I) - 1 
          stop 999
        ENDIF
#endif
 900  CONTINUE  
!
      DO 1000 I=NSTR,NEND
!
! ----------------------------
! > 全降水量の算出、単位変更 <
! ----------------------------
!    格子 B 平均
!
        PRECIP (I) =                                   ! 降水量合計。kg/m**2
     1               (RAIN_CONV(I)+RAIN_LARGE(I)) 
     1               * 1.D-3 * RHO_WTR_c               ! mm/mm*2->kg/m**2  
!    
! --------------------
! > 降水の分布関数形 <
! --------------------
!
        IF ( PRECIP(I) .LT. PRECIP_MIN ) THEN
          A_PR_COEF(I) = A_PL_COEF_c
          C_PR_COEF(I) = C_PL_COEF_c
        ELSE
          A_PR_COEF(I) = (   RAIN_CONV(I) *A_PC_COEF_C 
     1                     + RAIN_LARGE(I)*A_PL_COEF_c ) 
     1                   / ( RAIN_CONV(I) + RAIN_LARGE(I) )
          C_PR_COEF(I) = (   RAIN_CONV(I) *C_PC_COEF_c 
     1                     + RAIN_LARGE(I)*C_PL_COEF_C ) 
     1                   / ( RAIN_CONV(I) + RAIN_LARGE(I) )
        END IF
!
 1000 CONTINUE  
!
! ==================================================
! >> 詳細モニタその 1 ～ キャノピー上端に入る降水 <<    
! ==================================================
!
#ifdef MONYOS
!
!  キャノピーに入る降水フラックス
!
      IF (      IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN 
!
        DO I=NSTR,NEND
          PRECIP_W(I) = PRECIP(I) * FR_WBH_PR(I) 
          PRECIP_G(I) = PRECIP(I) * (1-FR_WBH_PR(I)) 
          WORK(I) 
     1     =  + (    CAP_WTR_c  * PRECIP_W(I)
     1            +  CAP_GLA_c  * PRECIP_G(I) )
     1            * ( TMP_PR (I) - TMP_FREZ_c )  
     1        - HL_W2G_c * PRECIP_G (I)
        ENDDO
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
          IF ( INSFLAG .EQ. 0 ) THEN
            CMARK  = 'WPRAN1'         
            CALL MONIT_SIB ( CMARK , PRECIP_W ,ICN, ONE , ONE ) 
            CMARK  = 'GPRAN1'         
            CALL MONIT_SIB ( CMARK , PRECIP_G ,ICN, ONE , ONE ) 
          ELSE
            CMARK  = 'WPRAS1'         
            CALL MONIT_SIB ( CMARK , PRECIP_W ,ICS, ONE , ONE ) 
            CMARK  = 'GPRAS1'         
            CALL MONIT_SIB ( CMARK , PRECIP_G ,ICS, ONE , ONE ) 
          ENDIF
!
          IF ( INSFLAG .EQ. 0 ) THEN
            CMARK  = 'HPRAN1'
            CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
          ELSE  
            CMARK  = 'HPRAS1'
            CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE ) 
          ENDIF
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          IF ( INSFLAG .EQ. 0 ) THEN
            CMARK  = 'A_RAINF'         
            CALL MONIT_SIB ( CMARK , PRECIP_W ,ICN, ONE , ONE ) 
            CMARK  = 'A_SNOWF'         
            CALL MONIT_SIB ( CMARK , PRECIP_G ,ICN, ONE , ONE ) 
          ELSE
            CMARK  = 'A_RAINF'         
            CALL MONIT_SIB ( CMARK , PRECIP_W ,ICS, ONE , ONE ) 
            CMARK  = 'A_SNOWF'         
            CALL MONIT_SIB ( CMARK , PRECIP_G ,ICS, ONE , ONE ) 
          ENDIF

        ENDIF

      ENDIF
#endif
!
!  --------------------------
!  > キャノピーでの一旦捕捉 <
!  --------------------------
!  ・注意。
!    ここは 格子 C 平均 （植生のある部分 CVR_VEG_msxだけでの平均値）にすべき。
!    現状は B 平均。
!
      DO 1100 I=NSTR,NEND
!       CHIV   = MAX ( CHI_L_co(NTYP(I),1) , CHIL_MIN )
        CHIV   = CHI_L_co(NTYP(I),1) 
        AA     = 0.500D0-0.633D0*CHIV-0.330D0*CHIV*CHIV     ! φ1
        BB     = 0.877D0-1.754D0*AA                         ! φ2
        EXTINC = AA + BB
        FR_TRAP_CNP (I)                       ! 引っかかる割合。
     1       = ( 1. - DEXP( - EXTINC * XLAI_VEG_msx(I,1) 
     1                              / CVR_VEG_msx(I,1) ) )
#ifdef CORRECT 
#else
     1                                * CVR_VEG_msx(I,1)   ! B へいきん
#endif
!
        PR_TRAP_CNP (I)                       ! 一旦捕捉される降水量（kg/m**2）
     1         = PRECIP(I) * FR_TRAP_CNP (I)       ! いまは B へいきん
!
        PR_TRAP_MAX (I)                       ! 捕捉上限量（kg/m**2）
     1         = ( H2O_SAT_CNP(I) - H2O_CNP(I) )   
#ifdef CORRECT 
     1                                / CVR_VEG_msx(I,1)    ! B->C
#else
#endif
!
 1100 CONTINUE
!
!      注意。remark 
!        従来 SiB or CORRECT しない場合と CORRECT した場合の違いについて。
!        従来は FPI * VCOVER(I,1) * PREPIP を LDRAIN への入力としていた。
!        これでは、キャノピーで一旦捕捉される量として、キャノピーの無いと
!        ころを含めて平均した B平均 の捕捉雨量を与えることになる。
!        これでは降水量が不当に少ない(薄められる)ことになり、
!        ldrain への入力としては適切ではない。
!        VCOVER はかけない、C 平均（FPI*PRECIP）を与えて、Cv 平均でトラップ
!        される量を見積もるべきで、最後に B 平均にすべきである。
!
! --------------------------------------------
! > キャノピーで一旦捕捉されたうち落下する量 <
! --------------------------------------------
!     ここも C 平均にすべきだが、とりあえず B 平均
!
      CALL WATER_OUT  (
     I     NSTR        , NEND        , 
     I     PR_TRAP_CNP , PR_TRAP_MAX , 
     I     A_PR_COEF   , B_PR_COEF_c , C_PR_COEF   , 
     O     PR_DROP_CNP )
!
!   キャノピーで捕捉される量  ... B 平均にする ( いまはすでに B へいきん ）
!
      DO 1300 I=NSTR,NEND
        PR_TRAP_CNP (I) = ( PR_TRAP_CNP(I) - PR_DROP_CNP(I) )   
#ifdef CORRECT
     1                    * CVR_VEG_msx(I,1)              !  C->B 
#else
#endif
 1300 CONTINUE  
!
!  ------------------------
!  > キャノピーでの相変化 <
!  ------------------------
!     TMP_FREZ_c での水の状態を基準とする。
!     格子 B 平均
!
      DO 1400 I=NSTR,NEND
!
        WTR_CNP = FR_WTR_CNP_IN(I) * H2O_SAT_CNP(I)  
        GLA_CNP = FR_GLA_CNP_IN(I) * H2O_SAT_CNP(I)
!
!  葉のもともと持つエネルギー
!
        ENG_CNP(I) 
     1     =   (   CAP_LAI_c * XLAI_VEG_msx(I,1)     
     1           + CAP_WTR_c * WTR_CNP               
     1           + CAP_GLA_c * GLA_CNP
     1         )  * ( TMP_CNP_IN(I) - TMP_FREZ_c )    
     1        -  HL_W2G_c * GLA_CNP
!
!  捕捉された水の持つエネルギー
!
        ENG_TRAP_CNP(I) 
     1     =  + (    CAP_WTR_c  * FR_WBH_PR(I) 
     1            +  CAP_GLA_c  * ( 1-FR_WBH_PR (I) )    )
     1            * PR_TRAP_CNP(I) * ( TMP_PR (I) - TMP_FREZ_c )  
     1        - HL_W2G_c * ( 1 - FR_WBH_PR (I) ) * PR_TRAP_CNP(I) 
!
        ENG_CNP(I) = ENG_CNP(I) + ENG_TRAP_CNP(I)
!
        H2O_CNP(I) = WTR_CNP    + GLA_CNP + PR_TRAP_CNP(I)  
        ENX    (I) = ENG_CNP(I) + HL_W2G_c * H2O_CNP(I) 
 1400 CONTINUE  
!
! ------------
! > 相の決定 <
! ------------
!    注意。水が無い場合 ENG_ICE = 0 なので、ELSE にはならない
!
      DO 1410 I=NSTR,NEND
        IF     ( ENG_CNP(I) .GE. 0. ) THEN                      
          FR_WBH_CNP(I) = 1. 
        ELSEIF ( ENX    (I) .LE. 0. ) THEN              
          FR_WBH_CNP(I) = 0. 
        ELSE 
          ENG_DIFF = MAX ( ENX(I) - ENG_CNP(I) , ENG_DIFF_MIN ) 
          FR_WBH_CNP(I) =  ENX(I) / ENG_DIFF 
        ENDIF
 1410 CONTINUE  
!
      DO 1420 I=NSTR,NEND
        WTR_CNP    = H2O_CNP(I) * FR_WBH_CNP(I)
        GLA_CNP    = H2O_CNP(I) * ( 1 - FR_WBH_CNP(I) )
        CAP_CNP    =   CAP_LAI_c * XLAI_VEG_msx(I,1) 
     1               + CAP_WTR_c * WTR_CNP
     1               + CAP_GLA_c * GLA_CNP
!
        TMP_CNP(I) = TMP_FREZ_c 
     1                  + (   MAX ( ENG_CNP(I) , DZERO )  
     1                      + MIN ( ENX(I)     , DZERO ) ) / CAP_CNP
        FR_H2O_CNP(I) = H2O_CNP(I) * H2O_SAT_CNP_INV(I) 
 1420 CONTINUE  
!
!  ------------------------
!  > キャノピーのあと処理 <
!  ------------------------
!     キャノピーを抜ける降水量に変更  ( B 平均 )
!
      DO 1500 I=NSTR,NEND 
        PRECIP (I)      = PRECIP(I) - PR_TRAP_CNP(I)         ! 透過量
 1500 CONTINUE  
!
!  ==========================================
!  >> キャノピーインターセプションここまで << 
!  ==========================================
!
! ==================================================
! >> 詳細モニタその 2 ～ 下草・地面上端に入る降水 <<    
! ==================================================
!
#ifdef MONYOS
!
!  キャノピーを通り抜けた、下草に入る降水フラックス
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
        IF ( INSFLAG .EQ. 0 ) THEN
!
          DO I=NSTR,NEND
            PRECIP_W(I) =  PRECIP(I) * FR_WBH_PR(I) 
            PRECIP_G(I) =  PRECIP(I) * (1-FR_WBH_PR(I) )
            WORK(I) 
     1       =  + (    CAP_WTR_c  * PRECIP_W(I)
     1              +  CAP_GLA_c  * PRECIP_G(I) )
     1              * ( TMP_PR (I) - TMP_FREZ_c )  
     1          - HL_W2G_c * PRECIP_G(I)
          ENDDO
!
          CMARK  = 'WPRAN2'         
          CALL MONIT_SIB ( CMARK , PRECIP_W ,ICN, ONE , ONE ) 
          CMARK  = 'GPRAN2'         
          CALL MONIT_SIB ( CMARK , PRECIP_G ,ICN, ONE , ONE ) 
!
          CMARK  = 'HPRAN2'
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
!
        ENDIF
      ENDIF
#endif
!
!  ====================================
!  >> 下草インターセプションここから <<   （雪無し域のみ）
!  ====================================
!
      IF ( INSFLAG.EQ.0 ) THEN     ! 雪無し域のみ
!
! ---------------
! >  保水量上限 <
! ---------------
!
      DO 2000 I=NSTR,NEND 
        H2O_SAT_GRS(I)
     1       = H2O_SAT_GRS_DEP_c     ! 水換算 0.1mm  
     1         * RHO_WTR_c           ! 1.D3          これで葉 1m^2 当たり質量
     1         * XLAI_VEG_msx(I,2)   !               これで単位面積当たり質量
!
        H2O_SAT_GRS_INV(I) = 1 / MAX( H2O_SAT_GRS(I), H2O_SAT_GRS_MIN )    
        H2O_GRS(I) = ( FR_WTR_GRS_IN(I)+FR_GLA_GRS_IN(I) ) 
     1               * H2O_SAT_GRS(I)
!
#ifdef SIB_DEBUG
        IF (  ( FR_WTR_GRS_IN(I) + FR_GLA_GRS_IN(I) ).GT.1+1.D-8 ) THEN
          write(6,*) 'canopy_interc error fr_wtr_grs gt 1 ' , 
     1                 I , 
     1                 FR_WTR_GRS_IN(I) , FR_GLA_GRS_IN(I) ,
     1                 FR_WTR_GRS_IN(I) + FR_GLA_GRS_IN(I) - 1 
          stop 999
        ENDIF
#endif
!
 2000 CONTINUE  
!
!  下草での一旦捕捉
!    
      DO 2100 I=NSTR,NEND
!       CHIV   = MAX ( CHI_L_co(NTYP(I),2) , CHIL_MIN )
        CHIV   = CHI_L_co(NTYP(I),2) 
        AA     = 0.500D0-0.633D0*CHIV-0.330D0*CHIV*CHIV     ! φ1
        BB     = 0.877D0-1.754D0*AA                         ! φ2
        EXTINC = AA + BB
        FR_TRAP_GRS (I) 
     1        =( 1.0D0 
     1           - DEXP( - EXTINC * XLAI_VEG_msx(I,2) 
     1                           / CVR_VEG_msx(I,2) ) )
#ifdef CORRECT
#else
     1         * CVR_VEG_msx(I,2)                   
#endif
        PR_TRAP_GRS (I) = PRECIP(I) * FR_TRAP_GRS (I)   ! 一旦捕捉される降水量
        PR_TRAP_MAX (I) = (H2O_SAT_GRS(I) - H2O_GRS(I)) ! 捕捉上限量
#ifdef CORRECT
     1                    / CVR_VEG_msx(I,2)                   
#endif
!
 2100 CONTINUE
!
!   下草で一旦捕捉されたうち落下する量
!     C 平均にすべきだが、とりあえず B 平均
!
      CALL WATER_OUT  (
     I     NSTR        , NEND        , 
     I     PR_TRAP_GRS , PR_TRAP_MAX , 
     I     A_PR_COEF   , B_PR_COEF_c , C_PR_COEF   , 
     O     PR_DROP_GRS )
!
!   下草で捕捉される量
!      B 平均にする。
!
      DO 2300 I=NSTR,NEND 
         PR_TRAP_GRS (I) = ( PR_TRAP_GRS(I) - PR_DROP_GRS(I) )
#ifdef CORRECT
     1                     * CVR_VEG_msx(I,2)                    ! ppp いる!
#endif
 2300 CONTINUE  
!
!  ------------------------
!  > 下草・地面での相変化 <
!  ------------------------
!     TMP_FREZ_c での水の状態を基準とする。
!     B 平均
!
      DO 2400 I=NSTR,NEND 
!
!  下草の葉がもともと持つ水                     fix 990917
!
        WTR_GRS = FR_WTR_GRS_IN (I) * H2O_SAT_GRS (I)
        GLA_GRS = FR_GLA_GRS_IN (I) * H2O_SAT_GRS (I)
!
!  下草・地面のもともと持つエネルギー
!
        CAP_GRSK    =   CAP_GRSK_c 
     1                 + CAP_LAI_c * XLAI_VEG_msx(I,2) 
     1                 + CAP_WTR_c * WTR_GRS
     1                 + CAP_GLA_c * GLA_GRS
!
        ENG_GRSK(I) =  CAP_GRSK
     1                   * ( TMP_GRSK_IN(I) - TMP_FREZ_c )    
     1               - HL_W2G_c * GLA_GRS
!
!  捕捉された水の持つエネルギー
!
        ENG_TRAP_GRS (I) 
     1      =       (    CAP_WTR_c  * FR_WBH_PR(I) 
     1                +  CAP_GLA_c  * ( 1 - FR_WBH_PR (I) )      )
     1              * PR_TRAP_GRS(I) * ( TMP_PR (I) - TMP_FREZ_c )  
     1            - HL_W2G_c * ( 1 - FR_WBH_PR (I) ) * PR_TRAP_GRS(I) 
!
        ENG_GRSK(I) = ENG_GRSK(I) + ENG_TRAP_GRS(I)
!
        H2O_GRS(I) = WTR_GRS + GLA_GRS + PR_TRAP_GRS(I)  
        ENX  (I)   = ENG_GRSK(I) + HL_W2G_c * H2O_GRS(I) 
 2400 CONTINUE  
!
! ------------
! > 相の決定 <
! ------------
!
      DO 2410 I=NSTR,NEND
        IF     ( ENG_GRSK(I) .GE. 0. ) THEN                      
          FR_WBH_GRS(I) = 1. 
        ELSEIF ( ENX  (I) .LE. 0. ) THEN              
          FR_WBH_GRS(I) = 0. 
        ELSE 
          ENG_DIFF = MAX ( ENX(I) - ENG_GRSK(I) , ENG_DIFF_MIN ) 
          FR_WBH_GRS(I) = ENX(I) / ENG_DIFF 
        ENDIF
 2410 CONTINUE  
!
      DO 2420 I=NSTR,NEND
        WTR_GRS = H2O_GRS(I) * FR_WBH_GRS(I)
        GLA_GRS = H2O_GRS(I) * ( 1 - FR_WBH_GRS(I) )
        CAP_GRSK    =   CAP_GRSK_c 
     1                 + CAP_LAI_c * XLAI_VEG_msx(I,2) 
     1                 + CAP_WTR_c * WTR_GRS
     1                 + CAP_GLA_c * GLA_GRS
        TMP_GRSK(I) = TMP_FREZ_c 
     1                  + (   MAX ( ENG_GRSK(I) , DZERO )  
     1                      + MIN ( ENX(I)   , DZERO ) ) / CAP_GRSK  
!
        FR_H2O_GRS(I) = H2O_GRS(I) * H2O_SAT_GRS_INV(I) 
 2420 CONTINUE  
!
!  ------------------------
!  > 下草のあと処理 <
!  ------------------------
!     土壌に落ちる降水量 PRECIP （B 平均）
!
      DO 2500 I=NSTR,NEND 
        PRECIP (I) = PRECIP(I) - PR_TRAP_GRS(I)        ! 透過量
 2500 CONTINUE  
!
      ENDIF                             ! 雪無し域のみ
!
!  ====================================
!  >> 下草インターセプションここまで <<   （雪無し域のみ）
!  ====================================
!
! ========================================
! >> 詳細モニタその 3 ～ 裸地に入る降水 <<    
! ========================================
!
#ifdef MONYOS
!
!  土壌・雪面に入る降水フラックス
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
        IF ( INSFLAG .EQ. 0 ) THEN
          DO I=NSTR,NEND
            WORK(I) = PRECIP(I) * FR_WBH_PR(I) 
          ENDDO
          CMARK  = 'WPRBNL'         
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
!
          DO I=NSTR,NEND
            WORK(I) = PRECIP(I) * (1-FR_WBH_PR(I)) 
          ENDDO
          CMARK  = 'GPRBNL'         
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
        ELSE
          DO I=NSTR,NEND
            WORK(I) = PRECIP(I) * FR_WBH_PR(I) 
          ENDDO
          CMARK  = 'WPRBSW'         
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE ) 
!
          DO I=NSTR,NEND
            WORK(I) = PRECIP(I) * (1-FR_WBH_PR(I)) 
          ENDDO
          CMARK  = 'GPRBSW'         
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE ) 
        ENDIF
!
!
        IF ( INSFLAG .EQ. 0 ) THEN
          DO I=NSTR,NEND
            WORK(I) 
     1       =     CAP_WTR_c  * FR_WBH_PR(I) 
     1           * PRECIP(I) * ( TMP_PR (I) - TMP_FREZ_c )  
          ENDDO
          CMARK  = 'IPRBNL'
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
!
          DO I=NSTR,NEND
            WORK(I) 
     1       =   CAP_GLA_c  * ( 1-FR_WBH_PR (I) )    
     1              * PRECIP(I) * ( TMP_PR (I) - TMP_FREZ_c )  
     1           - HL_W2G_c * ( 1 - FR_WBH_PR (I) ) * PRECIP(I)
          ENDDO
          CMARK  = 'JPRBNL'
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE , ONE ) 
        ELSE  
          DO I=NSTR,NEND
            WORK(I) 
     1       =     CAP_WTR_c  * FR_WBH_PR(I) 
     1           * PRECIP(I) * ( TMP_PR (I) - TMP_FREZ_c )  
          ENDDO
          CMARK  = 'IPRBSW'
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE ) 
!
          DO I=NSTR,NEND
            WORK(I) 
     1       =   CAP_GLA_c  * ( 1-FR_WBH_PR (I) )    
     1              * PRECIP(I) * ( TMP_PR (I) - TMP_FREZ_c )  
     1           - HL_W2G_c * ( 1 - FR_WBH_PR (I) ) * PRECIP(I)
          ENDDO
          CMARK  = 'JPRBSW'
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE ) 
        ENDIF
      ENDIF
#endif
!
!  ==============================
!  >> 落下した降水を、バケツへ <<  
!  ==============================
!    注。まだ。
!        ・ 単位のチェック。
!        ・ DELT_CANOPY について 
!
!    雪あり無しを含めた格子平均値にする。
!
      DO 3000 I=NSTR,NEND
        IF ( TMP_PR (I) .GE. TMP_FREZ_c ) THEN 
          H2O_BUCKET_INTO (I) 
     1        = H2O_BUCKET_INTO (I) 
     1          + PRECIP(I) * CVR_S_or_N(I) 
          ENH_BUCKET_INTO(I) 
     1        = ENH_BUCKET_INTO (I) 
     1          + PRECIP(I) * CAP_WTR_c * ( TMP_PR(I) - TMP_FREZ_c )
     1                      * CVR_S_or_N(I)
        ELSE 
          H2O_BUCKET_INTO (I) 
     1        = H2O_BUCKET_INTO (I) 
     1          + PRECIP(I) * CVR_S_or_N(I)
          ENH_BUCKET_INTO(I) 
     1        = ENH_BUCKET_INTO (I) 
     1          + PRECIP(I) 
     1               * ( CAP_GLA_c * ( TMP_PR(I) - TMP_FREZ_c )
     1                   - HL_W2G_c  )     
     1               * CVR_S_or_N(I)
        ENDIF
 3000 CONTINUE
!
! ==========================
! >> 詳細モニタ ～ おまけ <<    
! ==========================
!
#ifdef MONYOS
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'CVRN'         
          CALL MONIT_SIB ( CMARK, CVR_VEG_msx      ,ICN, DELT_CNP, ONE ) 
          CMARK  = 'CVRG'         
          CALL MONIT_SIB ( CMARK, CVR_VEG_msx(1,2) ,ICN, DELT_CNP, ONE ) 
        ELSE
          CMARK  = 'CVRS'         
          CALL MONIT_SIB ( CMARK, CVR_VEG_msx      ,ICS, DELT_CNP, ONE ) 
        ENDIF
!
      ENDIF
#endif
!
      RETURN
      END
!================================================================
      SUBROUTINE MR_CANOPY_INTERC
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , IMONIT_ALMA 
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_INTERC 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! 降水
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
      CMARK  = 'WPRAN1'
      CTITLE = 'PRECIPITATION INTO SIB (RAIN NOS)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GPRAN1'
      CTITLE = 'PRECIPITATION INTO SIB (SNOW NOS)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HPRAN1'
      CTITLE = 'HEAT OF PRECIPITATION INTO SIB (NOS)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WPRAN2'
      CTITLE = 'PRECIPITATION INTO GRASS (RAIN)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GPRAN2'
      CTITLE = 'PRECIPITATION INTO GRASS (SNOW)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HPRAN2'
      CTITLE = 'HEAT OF PRECIPITATION INTO GRASS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WPRBNL'
      CTITLE = 'PRECIPITATION WTR INTO SOIL (NOS)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GPRBNL'
      CTITLE = 'PRECIPITATION ICE INTO SNOW (NOS)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IPRBNL'
      CTITLE = 'HEAT OF PRECIPITATION WTR INTO SOIL (NOS)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JPRBNL'
      CTITLE = 'HEAT OF PRECIPITATION ICE INTO SNOW (NOS)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! ----------------
!
      CMARK  = 'WPRAS1'
      CTITLE = 'PRECIPITATION INTO SIB (SNW RAIN)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GPRAS1'
      CTITLE = 'PRECIPITATION INTO SIB (SNW SNOW)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'HPRAS1'
      CTITLE = 'HEAT OF PRECIPITATION INTO SIB (SNW)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WPRBSW'
      CTITLE = 'PRECIPITATION WTR INTO SNOW (SNW)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'GPRBSW'
      CTITLE = 'PRECIPITATION ICE INTO SNOW (SNW)'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IPRBSW'
      CTITLE = 'HEAT OF PRECIPITATION WTR INTO SNOW (SNW)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'JPRBSW'
      CTITLE = 'HEAT OF PRECIPITATION ICE INTO SNOW (SNW)'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!------------------------------------------------
!
      CMARK  = 'CVRN'
      CTITLE = 'COVERAGE OF VEG .. CANOPY NOS'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'CVRS'
      CTITLE = 'COVERAGE OF VEG .. CANOPY SNW'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'CVRG'
      CTITLE = 'COVERAGE OF VEG .. GRASS NOS'
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF    
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
      CMARK  = 'A_RAINF'
      CTITLE = 'ALMA Rainf'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'A_SNOWF'
      CTITLE = 'ALMA Snowf'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF 
#endif
!
      RETURN
      END

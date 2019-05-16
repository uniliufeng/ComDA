! 雪の短波アルベド・吸収関連定数の設定。 SNWALB 
! 99/05/xx 保坂征宏  第一版
! 99/07/11 保坂征宏      
! 99/07/23 保坂征宏  比率計算削除    
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 
      SUBROUTINE ALBEDO_SNOW  ( 
     I    NMAX        , 
     I    NMAX_01     , 
     I    NTYP        , 
     I    TMP_SNSK    , AGE_SNOW   , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     I    RAIN_1HOUR  ,  
     M    ALB_SNOW_SNW, ALB_SNOW_NOS,
     M    ALB_SNOW_onCNP_SNW  , ALB_SNOW_onCNP_NOS ) 
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IVN     
!
      USE SIBCON , ONLY : 
     1    ITUNE_ALBEDO_SNOW      
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX                   
      INTEGER,INTENT(IN)   :: NMAX_01 (3)                   
      INTEGER,INTENT(IN)   :: NTYP    (IDX)      ! 植生タイプ
!
      REAL(8),INTENT(INOUT):: AGE_SNOW (IDX)     ! 表面の雪の経過時間
      REAL(8),INTENT(IN)   :: TMP_SNSK (IDX)     ! 雪の表面温度

      REAL(8),INTENT(IN)   :: TMP_CNP_NOS (IDX)  ! 雪無しキャノピー温度
      REAL(8),INTENT(IN)   :: TMP_CNP_SNW (IDX)  ! 雪ありキャノピー温度
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_SNW (IDX,IVN)  ! 雪面アルベド
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_NOS (IDX,IVN)  
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_onCNP_SNW (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_onCNP_NOS (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(IN)     :: RAIN_1HOUR (IDX) 
!
      IF ( ITUNE_ALBEDO_SNOW .EQ. 0 ) THEN     ! JMA ルーチン方式
        CALL ALBEDO_SNOW_0   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
      ELSEIF ( ITUNE_ALBEDO_SNOW .EQ. 1 ) THEN   ! 加齢方式その 1 
        CALL ALBEDO_SNOW_1   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , AGE_SNOW   , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
      ELSEIF ( ITUNE_ALBEDO_SNOW .EQ. 2 ) THEN   ! MJ98 方式
        CALL ALBEDO_SNOW_2   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
      ELSEIF (      ITUNE_ALBEDO_SNOW .EQ. 3 
     1         .OR. ITUNE_ALBEDO_SNOW .EQ. 4 ) THEN   ! 加齢方式その 2 (積分)
        CALL ALBEDO_SNOW_3   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    RAIN_1HOUR  , 
     M    AGE_SNOW   , 
!     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     M    ALB_SNOW_SNW, ALB_SNOW_NOS,
     M    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
      ENDIF 
!
      RETURN
      END
!=====================================================================
      SUBROUTINE ALBEDO_SNOW_0   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
! JMA ルーチン方式
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IVN     
!
      USE SIBCON , ONLY : 
     1    TMP_FREZ_c           , 
     1    ALB_SNOW_VIS_c       ,
     1    ALB_SNOW_NIR_c       ,
     1    ALB_SNOW_FACT_TEMP_c    
!
      use sibveg , only : ityp_glacier 
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX                   
      INTEGER,INTENT(IN)   :: NMAX_01 (3)                   
      INTEGER,INTENT(IN)   :: NTYP    (IDX)       ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: TMP_SNSK(IDX)       ! 雪の表面温度

      REAL(8),INTENT(IN)   :: TMP_CNP_NOS (IDX)  ! 雪無しキャノピー温度
      REAL(8),INTENT(IN)   :: TMP_CNP_SNW (IDX)  ! 雪ありキャノピー温度
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: ALB_SNOW_SNW (IDX,IVN)  ! 雪面アルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_NOS (IDX,IVN)  
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_SNW (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_NOS (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
! WORK 
!
      LOGICAL,SAVE :: LFIRST = .TRUE. 
      REAL(8)      ::     DELTG             ! 温度差
      REAL(8)      ::     FWAVE 
      real(8)      ::     ALBSNW

      INTEGER :: I , JVN 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'ALBEDO_SNOW_0 ROUTINE type'
        WRITE(6,*) ' ALBEDO DEPENDS ON WAVE_LENGTH AND TEMP.'
        WRITE(6,*) '  VIS = ' , ALB_SNOW_VIS_c
        WRITE(6,*) '  NIR = ' , ALB_SNOW_NIR_c
        WRITE(6,*) '  FACT_TEMP = ' , ALB_SNOW_FACT_TEMP_c 
        LFIRST = .FALSE. 
      ENDIF
!
      alb_snow_nos(:,:) = 0. 
      alb_snow_snw(:,:) = 0. 
      alb_snow_oncnp_nos(:,:) = 0. 
      alb_snow_oncnp_snw(:,:) = 0. 
!  
      DO JVN = 1, 2                         ! 可視・近赤
!
! ------------
! > 雪あり域 <
! ------------
!
        DO I=NMAX_01(2),NMAX
!
!  ○雪面
!
!     波長依存性
!
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性
!
          DELTG = TMP_FREZ_c - TMP_SNSK(I)           
          IF ( DELTG.LT.0.05D0 ) THEN
            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_C
          ELSE
            ALBSNW = FWAVE 
          END IF
!
!     氷床上の例外処理
!
          IF ( (NTYP(I).EQ.ityp_glacier).AND.(JVN.EQ.1) ) THEN 
            IF ( DELTG.LT.0.05D0 ) THEN
              ALBSNW = 0.95D0 * 0.8D0   
            ELSE
              ALBSNW = 0.95D0 * 1.0D0   
            END IF
          ENDIF
!
!     アルベド値設定
!
          ALB_SNOW_SNW (I,JVN) = ALBSNW 
!
!  ○キャノピー
!
!     波長依存性
!  
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性 ( 0 度に近いときはファクタをかける ) 
!
          DELTG = TMP_FREZ_c - TMP_CNP_SNW(I)         
          IF ( DELTG.LT.0.05D0 ) THEN
            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_c 
          ELSE
            ALBSNW = FWAVE 
          END IF
!
!     アルベド値設定
!
          ALB_SNOW_onCNP_SNW (I,JVN) = ALBSNW  
        ENDDO
!
! ------------
! > 雪なし域 <  ( キャノピーのみ ) 
! ------------
!
        DO I=1,NMAX_01(3)
!
!     波長依存性
! 
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性
!
          DELTG = TMP_FREZ_c - TMP_CNP_NOS(I)           
          IF ( DELTG.LT.0.05D0 ) THEN
            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_c
          ELSE
            ALBSNW = FWAVE 
          END IF
!
!     アルベド値設定
!
          ALB_SNOW_onCNP_NOS (I,JVN) = ALBSNW  
          ALB_SNOW_NOS       (I,JVN) = ALBSNW   !dummy
        ENDDO
!
      ENDDO           
!
      RETURN 
      END SUBROUTINE ALBEDO_SNOW_0
!=====================================================================
      SUBROUTINE ALBEDO_SNOW_1   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , AGE_SNOW   , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
! JMA ルーチン方式
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IVN     
!
      USE SIBCON , ONLY : 
     1    TMP_FREZ_c             , 
     1    AGE_SNOW_MAX1_c        , 
     1    AGE_SNOW_ETIME_INV_c   ,
     1    ALB_SNOW_VIS_INI_c     ,
     1    ALB_SNOW_VIS_FOR_c     ,
     1    ALB_SNOW_NIR_INI_c     ,
     1    ALB_SNOW_NIR_FOR_c     ,
     1    ALB_SNOW_FACT_TEMP_c

!
      use sibveg , only : ityp_glacier 

!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX                   
      INTEGER,INTENT(IN)   :: NMAX_01 (3)                   
      INTEGER,INTENT(IN)   :: NTYP    (IDX)       ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: AGE_SNOW (IDX)      ! 表面の雪の経過時間
      REAL(8),INTENT(IN)   :: TMP_SNSK(IDX)       ! 雪の表面温度

      REAL(8),INTENT(IN)   :: TMP_CNP_NOS (IDX)  ! 雪無しキャノピー温度
      REAL(8),INTENT(IN)   :: TMP_CNP_SNW (IDX)  ! 雪ありキャノピー温度
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: ALB_SNOW_SNW (IDX,IVN)  ! 雪面アルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_NOS (IDX,IVN)  
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_SNW (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_NOS (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
!
! SAVE 
!
      REAL(8),SAVE   ::  ALB_INI(IVN)      ! アルベド初期値
      REAL(8),SAVE   ::  ALB_FOR(IVN)      ! アルベド向かう値
      REAL(8),SAVE   ::  ALB_END(IVN)      ! アルベド最終値
!
! WORK 
!
      REAL(8)      :: DELTG 
      LOGICAL,SAVE :: LFIRST = .TRUE. 
      INTEGER :: I , JVN 
!
      IF ( LFIRST ) THEN
        ALB_INI(1) = ALB_SNOW_VIS_INI_c     ! 可視アルベド初期値    0.95
        ALB_INI(2) = ALB_SNOW_NIR_INI_c     ! 近赤アルベド初期値    0.75 
        ALB_FOR(1) = ALB_SNOW_VIS_FOR_c     ! 可視アルベド向かう値  0.65 
        ALB_FOR(2) = ALB_SNOW_NIR_FOR_c     ! 近赤アルベド向かう値  0.45 
!
        ALB_END(1) = ALB_FOR(1) 
     1              + ( ALB_INI(1) - ALB_FOR(1) ) 
     1             * DEXP ( - AGE_SNOW_MAX1_c * AGE_SNOW_ETIME_INV_c )
        ALB_END(2) = ALB_FOR(2) 
     1              + ( ALB_INI(2) - ALB_FOR(2) ) 
     1             * DEXP ( - AGE_SNOW_MAX1_c * AGE_SNOW_ETIME_INV_c)
!
        write(6,*) 'albedo_snow etime' , 1/ AGE_SNOW_ETIME_INV_c
        write(6,*) 'albedo_snow VIS : ini for end' , 
     1                alb_ini(1), alb_for(1), alb_end(1) 
        write(6,*) 'albedo_snow IR  : ini for end' , 
     1                alb_ini(2), alb_for(2), alb_end(2) 
        LFIRST = .FALSE. 
      ENDIF
!
      alb_snow_nos(:,:) = 0. 
      alb_snow_snw(:,:) = 0. 
      alb_snow_oncnp_nos(:,:) = 0. 
      alb_snow_oncnp_snw(:,:) = 0. 
!
      DO JVN = 1, 2
!
!  ==================
!  >> アルベド計算 <<
!  ==================
!
!  ----------------------------
!  > 雪無し域 ( onCNP のみ ) <
!  ----------------------------
!
        DO I=1,NMAX_01(3)
          ALB_SNOW_onCNP_NOS (I,JVN) = ALB_END(JVN) 
          ALB_SNOW_NOS       (I,JVN) = ALB_END(JVN)    ! ダミー
!
          DELTG = TMP_FREZ_c - TMP_CNP_NOS(I)             
          IF ( DELTG.LT.0.05D0 ) THEN                  ! ぬれ雪
            ALB_SNOW_onCNP_NOS(I,JVN) 
!    1          = ALB_SNOW_onCNP_NOS(I,JVN) * 0.6 
     1          = ALB_SNOW_onCNP_NOS(I,JVN) * ALB_SNOW_FACT_TEMP_c 
          END IF
        ENDDO
!
!  -------------
!  > 雪あり域  <
!  -------------
!     新雪なら時間について exp で減少する。旧雪なら定数。
!
        DO I=NMAX_01(2),NMAX
!
!  ベース
!
          IF ( AGE_SNOW(I) .GT. AGE_SNOW_MAX1_c ) THEN
            ALB_SNOW_SNW(I,JVN) = ALB_END(JVN)
          ELSE 
            ALB_SNOW_SNW(I,JVN) = ALB_FOR(JVN) 
     1                + ( ALB_INI(JVN) - ALB_FOR(JVN) ) 
     1                  * DEXP ( - AGE_SNOW(I) * AGE_SNOW_ETIME_INV_c )
          ENDIF
!
!  キャノピー
!
          ALB_SNOW_onCNP_SNW(I,JVN) = ALB_SNOW_SNW(I,JVN)
          DELTG = TMP_FREZ_c - TMP_CNP_SNW(I)             
          IF ( DELTG.LT.0.05D0 ) THEN                   ! ぬれ雪
            ALB_SNOW_onCNP_SNW(I,JVN) 
!    1             = ALB_SNOW_onCNP_SNW(I,JVN) * 0.6D0           
     1             = ALB_SNOW_onCNP_SNW(I,JVN) * ALB_SNOW_FACT_TEMP_c           
          END IF
!
!  雪面
!
          DELTG   =  TMP_FREZ_c - TMP_SNSK(I)           
          IF ( DELTG.LT.0.05D0 ) THEN                 
!           ALB_SNOW_SNW(I,JVN) = ALB_SNOW_SNW(I,JVN) * 0.6            
            ALB_SNOW_SNW(I,JVN) = ALB_SNOW_SNW(I,JVN) * 
     1                                         ALB_SNOW_FACT_TEMP_c            
          END IF
!
! ----------------------
! > 氷床アルベドの補正 <
! ----------------------
!    雪面可視のみ。キャノピーはそもそもあまりないだろうから、補正なし。
!
          IF ( ( NTYP(I).EQ.ityp_glacier).AND.(JVN.EQ.1) ) THEN 
            IF ( DELTG.LT.0.05D0 ) THEN
              ALB_SNOW_SNW(I,JVN) = 0.95D0 * 0.8D0   
            ELSE
              ALB_SNOW_SNW(I,JVN) = 0.95D0 * 1.0D0   
            END IF
          ENDIF
!
        ENDDO
!
      ENDDO
!
      RETURN 
      END SUBROUTINE ALBEDO_SNOW_1
!=====================================================================
      SUBROUTINE ALBEDO_SNOW_2   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    TMP_CNP_SNW , TMP_CNP_NOS, 
     O    ALB_SNOW_SNW, ALB_SNOW_NOS,
     O    ALB_SNOW_onCNP_SNW       , ALB_SNOW_onCNP_NOS        ) 
!
!  For MRI/JMA98 ( 行本 L3SiB チューニング ) 
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IVN     
!
      USE SIBCON , ONLY : 
     1    TMP_FREZ_c           , 
     1    ALB_SNOW_VIS_c       ,
     1    ALB_SNOW_NIR_c       ,
     1    ALB_SNOW_FACT_TEMP_c   , 
     1    ALB_SNOW_VIS_FOR_c     ,
     1    ALB_SNOW_NIR_FOR_c     
!
      use sibveg , only : ityp_glacier 

!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX                   
      INTEGER,INTENT(IN)   :: NMAX_01 (3)                   
      INTEGER,INTENT(IN)   :: NTYP    (IDX)       ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: TMP_SNSK(IDX)       ! 雪の表面温度

      REAL(8),INTENT(IN)   :: TMP_CNP_NOS (IDX)  ! 雪無しキャノピー温度
      REAL(8),INTENT(IN)   :: TMP_CNP_SNW (IDX)  ! 雪ありキャノピー温度
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: ALB_SNOW_SNW (IDX,IVN)  ! 雪面アルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_NOS (IDX,IVN)  
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_SNW (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(OUT)  :: ALB_SNOW_onCNP_NOS (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
!
! WORK 
!
      LOGICAL,SAVE :: LFIRST = .TRUE. 
      REAL(8)      :: ALBXX 
      REAL(8)      :: ALBSNW
      REAL(8)      :: FWAVE
      INTEGER :: I , JVN 
    
!
      IF ( LFIRST ) THEN  
        WRITE(6,*) 'ALBEDO_SNOW_2 : ' , 
     1             ' ALBEDO DEPENDS ON WAVE_LENGTH AND TEMP.'
        WRITE(6,*) '  VIS = ' , ALB_SNOW_VIS_c
        WRITE(6,*) '  NIR = ' , ALB_SNOW_NIR_c
        WRITE(6,*) '  FACT_TEMP = ' , ALB_SNOW_FACT_TEMP_c 
        LFIRST = .FALSE. 
      ENDIF
!
      alb_snow_nos(:,:) = 0. 
      alb_snow_snw(:,:) = 0. 
      alb_snow_oncnp_nos(:,:) = 0. 
      alb_snow_oncnp_snw(:,:) = 0. 
!
! ------------
! > 雪あり域 <
! ------------
!
      DO JVN=1,IVN
        DO I=NMAX_01(2),NMAX
!
!  ○雪面
!
!     波長依存性
!
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性
!
!!          DELTG = TMP_FREZ_c - TMP_SNSK(I)           
!!          IF ( DELTG.LT.0.05D0 ) THEN
!!            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_C
!!          ELSE
!!            ALBSNW = FWAVE 
!!          END IF
          ALBXX = 1.D0 - ( TMP_SNSK(I) -269.1D0 ) * 0.05
          ALBXX = MAX(0.7D0,MIN(1.D0,ALBXX))
          ALBSNW = FWAVE * ALBXX 
!
!     氷床上の例外処理
!
!!          IF ( (NTYP(I).EQ.25).AND.(JVN.EQ.1) ) THEN 
!!            IF ( DELTG.LT.0.05D0 ) THEN
!!              ALBSNW = 0.95D0 * 0.8D0   
!!            ELSE
!!              ALBSNW = 0.95D0 * 1.0D0   
!!            END IF
!!          ENDIF
!
          IF ( NTYP(I).EQ.ityp_glacier ) THEN
            IF (JVN.EQ.1) THEN 
              ALBSNW = 0.95D0 
            ELSE
              ALBXX = 0.4D0 + (270.D0-TMP_SNSK(I))*0.2D0
              ALBSNW = MAX ( 0.4D0, MIN(0.7D0,ALBXX))
            END IF
          ENDIF
!
!     アルベド値設定
!
          ALB_SNOW_SNW (I,JVN) = ALBSNW 
!
!  ○キャノピー
!
!     波長依存性
!  
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性 ( 0 度に近いときはファクタをかける ) 
!
!!          DELTG = TMP_FREZ_c - TMP_CNP_SNW(I)         
!!          IF ( DELTG.LT.0.05D0 ) THEN
!!            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_c 
!!          ELSE
!!            ALBSNW = FWAVE 
!!          END IF
!
          ALBXX = 1.D0 - ( TMP_SNSK(I) -269.1D0 ) * 0.05
          ALBXX = MAX(0.7D0,MIN(1.D0,ALBXX))
          ALBSNW = FWAVE * ALBXX 
!
          IF ( NTYP(I).EQ.ityp_glacier ) THEN
            IF (JVN.EQ.1) THEN 
              ALBSNW = 0.95D0 
            ELSE
              ALBXX = 0.4D0 + (270.D0-TMP_SNSK(I))*0.2D0
              ALBSNW = MAX ( 0.4D0, MIN(0.7D0,ALBXX))
            END IF
          ENDIF
!
!     アルベド値設定
!
          ALB_SNOW_onCNP_SNW (I,JVN) = ALBSNW  
        ENDDO
!
! ------------
! > 雪なし域 <  ( キャノピーのみ ) 
! ------------
!
        DO I=1,NMAX_01(3)
!
!     波長依存性
! 
          IF      ( JVN .EQ. 1 ) THEN  ! 可視
            FWAVE = ALB_SNOW_VIS_c
          ELSEIF  ( JVN .EQ. 2 ) THEN  ! 近赤外
            FWAVE = ALB_SNOW_NIR_c
          ENDIF
!
!     温度依存性
!
!!          DELTG = TMP_FREZ_c - TMP_CNP_NOS(I)           
!!          IF ( DELTG.LT.0.05D0 ) THEN
!!            ALBSNW = FWAVE * ALB_SNOW_FACT_TEMP_c
!!          ELSE
!!            ALBSNW = FWAVE 
!!          END IF
!
          ALBXX = 1.D0 - ( TMP_SNSK(I) -269.1D0 ) * 0.05
          ALBXX = MAX(0.7D0,MIN(1.D0,ALBXX))
          ALBSNW = FWAVE * ALBXX 
!
          IF ( NTYP(I).EQ.ityp_glacier ) THEN
            IF (JVN.EQ.1) THEN 
              ALBSNW = 0.95D0 
            ELSE
              ALBXX = 0.4D0 + (270.D0-TMP_SNSK(I))*0.2D0
              ALBSNW = MAX ( 0.4D0, MIN(0.7D0,ALBXX))
            END IF
          ENDIF
!
!     アルベド値設定
!
          ALB_SNOW_onCNP_NOS (I,JVN) = ALBSNW  
          ALB_SNOW_NOS       (I,JVN) = ALBSNW   !dummy
        ENDDO
!
      ENDDO           
!
      RETURN 
      END SUBROUTINE ALBEDO_SNOW_2
!=====================================================================
      SUBROUTINE ALBEDO_SNOW_3   ( 
     I    NMAX        , NMAX_01    , 
     I    NTYP        , 
     I    TMP_SNSK    , 
     I    RAIN_1HOUR  , 
     M    AGE_SNOW    , 
!     I    TMP_CNP_SNW , TMP_CNP_NOS , 
     M    ALB_SNOW_SNW, ALB_SNOW_NOS,
     M    ALB_SNOW_onCNP_SNW        , ALB_SNOW_onCNP_NOS        ) 
!
!  T.Aoki et al. (2003) 方式
!    年齢は、若い ( 1 時間以下 ) かどうかの情報のみを使用する。
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    IVN     
!
      USE SIBCON , ONLY : 
     1    TMP_FREZ_c             , 
     1    ALB_SNOW_VIS_FOR_c     ,
     1    ALB_SNOW_NIR_FOR_c     ,
     1    ITUNE_ALBEDO_SNOW      ,
     1    RAIN_1HOUR_FRESH_c     ,
     1    AGE_SNOW_MAX2_c     
!
!
      use sibveg , only : ityp_glacier 

! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX                   
      INTEGER,INTENT(IN)   :: NMAX_01 (3)                   
      INTEGER,INTENT(IN)   :: NTYP    (IDX)      ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: TMP_SNSK (IDX)     ! 雪の表面温度
!
      REAL(8),INTENT(IN)   :: RAIN_1HOUR (IDX)  

!      REAL(8),INTENT(IN)   :: TMP_CNP_NOS (IDX)  ! 雪無しキャノピー温度
!      REAL(8),INTENT(IN)   :: TMP_CNP_SNW (IDX)  ! 雪ありキャノピー温度
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_SNW (IDX,IVN)  ! 雪面アルベド
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_NOS (IDX,IVN)  
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_onCNP_SNW (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
      REAL(8),INTENT(INOUT)  :: ALB_SNOW_onCNP_NOS (IDX,IVN) 
!                                                 ! キャの葉の雪のアルベド
!
      REAL(8),INTENT(INOUT)  :: AGE_SNOW (IDX) 

!
      REAL(8),SAVE :: DELT_CALC_ALBEDO = 3600.D0   ! unit sec 
                                       ! 1 時間に一度の計算を仮定
!
! SAVE 
!
!      REAL(8),SAVE   ::  ALB_STR(IVN)      ! 初期値
      REAL(8),SAVE   ::  ALB_END(IVN)      ! 最小値
!
      LOGICAL,SAVE   :: LFIRST = .TRUE. 
!
! WORK 
!
      REAL(8)      :: DEL 
      REAL(8)      :: ALBSNW 
      INTEGER      :: I , JVN 
      real(8),parameter :: zero = 0.D0 
      real(8),parameter :: one  = 1.D0 
      REAL(8)      :: TMP
      REAL(8)      :: ALBXX
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'ALBEDO_SNOW_3 Aoki type ' , itune_albedo_snow 
        LFIRST = .FALSE. 
        ALB_END(1) = ALB_SNOW_VIS_FOR_c 
        ALB_END(2) = ALB_SNOW_NIR_FOR_c 
!
        write(6,*) '   visible goto ' , alb_end(1) 
        write(6,*) '   ir      goto ' , alb_end(2) 
!        ALB_STR(1) = ALB_SNOW_VIS_INI_c 
!        ALB_STR(2) = ALB_SNOW_NIR_INI_c 
      ENDIF
!
      DO I=1,NMAX_01(2)-1
        AGE_SNOW     (I) = AGE_SNOW_MAX2_c
      ENDDO
      DO I=NMAX_01(2),NMAX 
        AGE_SNOW     (I) = AGE_SNOW     (I)  + DELT_CALC_ALBEDO
        AGE_SNOW     (I) = MIN ( AGE_SNOW(I) , AGE_SNOW_MAX2_c)
      ENDDO
!
! 計算方式
!   基本的には予報方程式
!   年齢が若ければ ( 1 時間以内 ) 新値にリセットする。
!     
      DO JVN = 1, 2                         ! 可視・近赤
!
! ------------
! > 雪なし域 <  ( キャノピーのみ。下限値を与える。) 
! ------------
!
        DO I=1,NMAX_01(3)
!
!     波長依存性
! 
          ALBSNW = ALB_END(JVN) 
!
!     温度依存性 ( 濡れ雪 ) 
!         DELTG = TMP_FREZ_c - TMP_CNP_NOS(I)           
!         IF ( DELTG.LT.0.05D0 ) THEN
!           ALBSNW = ALBSNW * ALB_SNOW_FACT_TEMP_c
!         END IF
!
!     アルベド値設定
!
          ALB_SNOW_onCNP_NOS (I,JVN) = ALBSNW  
          ALB_SNOW_NOS       (I,JVN) = ALBSNW   !dummy
!
        ENDDO
!
! ------------
! > 雪あり域 <
! ------------
!
        DO I=NMAX_01(2),NMAX
!
!  ○雪面
!     
!      アルベド減率 ( 1 時間当たり )    
!
          IF     ( TMP_SNSK(I) - TMP_FREZ_c .LT. -20 ) THEN
            DEL = - 0.321D-4 * ( DELT_CALC_ALBEDO / 3600.D0 )
          ELSEIF ( TMP_SNSK(I) - TMP_FREZ_c .LT. -10 ) THEN
            DEL = - 1.60 D-4 * ( DELT_CALC_ALBEDO / 3600.D0 )
          ELSEIF ( TMP_SNSK(I) - TMP_FREZ_c .LT. - 5 ) THEN
            DEL = - 4.33 D-4 * ( DELT_CALC_ALBEDO / 3600.D0 )
          ELSE
            DEL = - 7.18 D-4 * ( DELT_CALC_ALBEDO / 3600.D0 )
          ENDIF
          DEL = DEL * ( 0.4 + JVN * 0.4 )   ! 可視 1 は 0.8 , 近赤 2 は 1.2
!
!      時間積分
!
          ALB_SNOW_SNW (I,JVN) = ALB_SNOW_SNW (I,JVN) + DEL 
!
!      下限処理
!
          ALB_SNOW_SNW (I,JVN) 
     1         = MAX ( ALB_SNOW_SNW (I,JVN) , ALB_END(JVN) )
!
!      新雪時はリセット
!        可視 : -2.5 度で 0.9 , -7.5 度 0.95
!        近赤 : -2.5 度で 0.7 , -7.5 度 0.75
!
!         IF ( AGE_SNOW (I) .LT. DELT_CALC_ALBEDO - 0.1 ) THEN
          IF ( RAIN_1HOUR (I) .GT. RAIN_1HOUR_FRESH_c ) THEN
!
!               TMP : -2.5 より高い → 0 , -7.5 より低い → 1 
!               可視(1) 0.95or0.9 から、近赤(2) : 0.75or0.7 から
!
            TMP =  - ( TMP_SNSK(I) - TMP_FREZ_c + 2.5 )/5  
            TMP = MIN ( MAX ( TMP , ZERO ) , ONE ) 

            ALB_SNOW_SNW(I,JVN) = TMP * 0.05 + 1.1 - JVN * 0.2  
            AGE_SNOW    (I)     = 0. 
!
          ENDIF 
!
          ALB_SNOW_onCNP_SNW (I,JVN) = ALB_SNOW_SNW(I,JVN) 
!
        ENDDO
!
!     氷床上の例外処理 ( ITUNE_ALBEDO_SNOW .EQ. 4 のときのみ ) 
!
        IF ( ITUNE_ALBEDO_SNOW .EQ. 4 ) THEN
        DO I=NMAX_01(2),NMAX
          IF ( NTYP(I).EQ.ityp_glacier ) THEN
            IF (JVN.EQ.1) THEN 
              ALBSNW = 0.95D0 
            ELSE
              ALBXX = 0.4D0 + (270.D0-TMP_SNSK(I))*0.2D0
              ALBSNW = MAX ( 0.4D0, MIN(0.7D0,ALBXX))
            END IF
          ENDIF
        ENDDO
        ENDIF
!
      ENDDO           
!
      RETURN
      END SUBROUTINE ALBEDO_SNOW_3   


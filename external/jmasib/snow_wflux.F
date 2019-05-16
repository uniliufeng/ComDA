! 旧雪層を解く 第一段。ただし旧雪のとり方はこれまでとあまり変えない場合。
! 99/07/10 保坂征宏   新規作成
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SNOW_WFLUX (
     I   NMAX         , 
     I   NMAX_01      , NMAX_2    , NMAX_3    , 
     I   H2O_PRECIPI  , ENG_PRECIPI   , RHO_PRECIPI_INV       , 
     I   H2O_SNOW     , ENG_SNOW      , RHO_SNOW_INV        , 
     I   I_SNOW_LEV   , 
     M   IFLG_SLEV_CHANGE      , 
     O   ENG_SNOW_NEXT , H2O_SNOW_NEXT       , 
     O   RHO_SNOW_INV_NEXT   , 
     O   ENH_BUCKET_SNOW_TO_SNSL  , WTR_BUCKET_SNOW_TO_SNSL  ,    
     O   INFO_SNOW_NEXT         )  
!
!   ・バケツの単位はこれでよい。水周り計算では格子平均値であるため。
!     
      USE SIBPRM , ONLY :
     1   IDX  ,
     1   ISN 
!
      USE SIBCON , ONLY :
     1   HL_W2G_c        , 
     1   H2O_SNOW1_MIN_c , 
     1   H2O_SNOW1_MAX_c , 
     1   H2O_SNOW2_MIN_c , 
     1   H2O_SNOW2_MAX_c , 
     1   H2O_SNOW3_MIN_c , 
     1   H2O_SNOW3_MAX_c  
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_01(3)
      INTEGER,INTENT(IN)   :: NMAX_2 (3)
      INTEGER,INTENT(IN)   :: NMAX_3 (3)
!     INTEGER,INTENT(IN)   :: NMAX_4 (3)
!
      REAL(8),INTENT(IN)   :: ENG_SNOW         (IDX,ISN)
      REAL(8),INTENT(IN)   :: H2O_SNOW         (IDX,ISN)
      REAL(8),INTENT(IN)   :: RHO_SNOW_INV     (IDX,ISN)  
      REAL(8),INTENT(IN)   :: ENG_PRECIPI      (IDX)
      REAL(8),INTENT(IN)   :: H2O_PRECIPI      (IDX)
      REAL(8),INTENT(IN)   :: RHO_PRECIPI_INV  (IDX)       ! 降雪の密度
!
      INTEGER,INTENT(IN)   :: I_SNOW_LEV       (IDX)

!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: H2O_SNOW_NEXT     (IDX,ISN)
      REAL(8),INTENT(OUT)  :: ENG_SNOW_NEXT     (IDX,ISN)
      REAL(8),INTENT(OUT)  :: RHO_SNOW_INV_NEXT (IDX,ISN)
!
      REAL(8),INTENT(OUT)  :: WTR_BUCKET_SNOW_TO_SNSL (IDX)
      REAL(8),INTENT(OUT)  :: ENH_BUCKET_SNOW_TO_SNSL (IDX)
!
      INTEGER,INTENT(OUT)    :: INFO_SNOW_NEXT   (IDX)
!
! INOUT
!
      INTEGER,INTENT(INOUT):: IFLG_SLEV_CHANGE  (IDX) 
!
! WORK
      INTEGER  :: I 
      REAL(8)  :: H2O_1_TO_2 (IDX)   
      REAL(8)  :: H2O_2_TO_3 (IDX)
      REAL(8)  :: H2O_3_TO_4 (IDX)
!
      REAL(8)  :: H2O_TO_SOIL(IDX)  
!
      REAL(8),PARAMETER :: ZERO_IDX   (IDX) = 0.
!
      REAL(8)  :: FR
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_OLD_LAYER POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
      H2O_TO_SOIL(:) = 0. 
!
! ==================
! >> 第一層を解く <<
! ==================
!  降雪を混ぜる
!
      CALL MK_SNOW_LAYER (
     I   NMAX_01(2)        , NMAX             , 1           ,
     I   ENG_PRECIPI       , H2O_PRECIPI      , 
     I   H2O_SNOW1_MIN_c   , H2O_SNOW1_MAX_c  , 
     I   ENG_SNOW          , H2O_SNOW         ,    
     O   H2O_1_TO_2        , 
     M   IFLG_SLEV_CHANGE  , H2O_TO_SOIL      )
!
! ==================
! >> 第二層を解く <<
! ==================
!
      CALL MK_SNOW_LAYER (
     I   NMAX_2(2)         , NMAX             , 2          ,
     I   ZERO_IDX          , H2O_1_TO_2       , 
     I   H2O_SNOW2_MIN_c   , H2O_SNOW2_MAX_c , 
     I   ENG_SNOW          , H2O_SNOW      ,    
     O   H2O_2_TO_3        , 
     M   IFLG_SLEV_CHANGE  , H2O_TO_SOIL      )
!
! ==================
! >> 第三層を解く <<
! ==================
!
      CALL MK_SNOW_LAYER (
     I   NMAX_3(2)         , NMAX             , 3         ,
     I   ZERO_IDX          , H2O_2_TO_3       , 
     I   H2O_SNOW3_MIN_c   , H2O_SNOW3_MAX_c , 
     I   ENG_SNOW          , H2O_SNOW      ,    
     O   H2O_3_TO_4        , 
     M   IFLG_SLEV_CHANGE  , H2O_TO_SOIL      )
!
! ==============
! >> 更新作業 <<   
! ==============
!    ただし、層切り直しではない場合
!
      INFO_SNOW_NEXT(:) = -1                    ! 初期値は雪無し
!
      ENH_BUCKET_SNOW_TO_SNSL(:) = 0.
      WTR_BUCKET_SNOW_TO_SNSL(:) = 0.
      H2O_SNOW_NEXT(:,:) = 0.
      ENG_SNOW_NEXT(:,:) = 0.
      RHO_SNOW_INV_NEXT (:,:) = 0. 
!
      DO I=NMAX_01(2),NMAX
        IF ( IFLG_SLEV_CHANGE(I).EQ.0 ) THEN
!
! ----------
! > 層の数 < 
! ----------
!
          INFO_SNOW_NEXT(I) = I_SNOW_LEV(I)     ! 層数を代入 
!
! ----------
! > 第一層 < 
! ----------
!
!    密度の計算。水の落下は無視して比を見積もる。（いい加減?）
!
          FR  = H2O_PRECIPI(I) / ( H2O_PRECIPI(I) + H2O_SNOW(I,1) ) 
!                
          RHO_SNOW_INV_NEXT (I,1) 
     1            =   FR      * RHO_PRECIPI_INV (I)
     1              + (1-FR)  * RHO_SNOW_INV    (I,1)
!
!     水とエネルギーの計算
!
          H2O_SNOW_NEXT (I,1) = H2O_SNOW (I,1) 
     1                         + H2O_PRECIPI(I) - H2O_1_TO_2(I)       
          ENG_SNOW_NEXT (I,1) = ENG_SNOW (I,1) 
     1                         + ENG_PRECIPI(I) 
!
       ENDIF 
      ENDDO 
!
!
! --------------
! > 第二層 < 
! --------------
!
!    密度の計算。水の落下は無視して比を見積もる。（いい加減?）
!      ----> 第二層の密度のままということだ。
!
      DO I=NMAX_2(2),NMAX
        IF ( IFLG_SLEV_CHANGE(I).EQ.0 ) THEN

          RHO_SNOW_INV_NEXT (I,2) = RHO_SNOW_INV (I,2)
!
!    水とエネルギーの計算
!
          H2O_SNOW_NEXT (I,2) = H2O_SNOW(I,2) 
     1                         + H2O_1_TO_2(I) - H2O_2_TO_3(I)       
          ENG_SNOW_NEXT (I,2) = ENG_SNOW(I,2) 
        ENDIF   
      ENDDO
!

!
! ----------
! > 第三層 < 
! ----------
!
!    密度の計算。水の落下は無視して比を見積もる。（いい加減?）
!      ----> 第三層の密度のままということだ。
!
      DO I=NMAX_3(2),NMAX
        IF ( IFLG_SLEV_CHANGE(I).EQ.0 ) THEN
          RHO_SNOW_INV_NEXT (I,3) = RHO_SNOW_INV(I,3)
!
!    水とエネルギーの計算
!
          H2O_SNOW_NEXT (I,3) = H2O_SNOW(I,3) 
     1                         + H2O_2_TO_3(I) - H2O_3_TO_4(I)       
          ENG_SNOW_NEXT (I,3) = ENG_SNOW(I,3) 
        ENDIF
      ENDDO
!
! ================================
! >> 土壌への水フラックスの計算 <<
! ================================
!        ( 下に落さないことを前提としていれば、不要 )
!
      DO I=NMAX_01(2),NMAX_3(3)
        IF ( IFLG_SLEV_CHANGE(I).EQ.0 ) THEN
          WTR_BUCKET_SNOW_TO_SNSL(I) = H2O_TO_SOIL(I)             
          ENH_BUCKET_SNOW_TO_SNSL(I) = 0. 
        ENDIF
      ENDDO
!
      RETURN
      END
!===============================================
      SUBROUTINE MK_SNOW_LAYER (
     I   NSTR   , NEND    , KLEV      ,
     I   ENG_FROM_UPPER   , H2O_FROM_UPPER   , 
     I   H2O_SNOW_MIN     , H2O_SNOW_MAX     , 
     I   ENG_SNOW         , H2O_SNOW         ,    
     O   H2O_TO_LOWER     , 
     M   IFLG_SLEV_CHANGE , H2O_TO_SOIL     )
!
! 上から来た水（最上層では +氷）と既存の氷+水を混ぜる。
! 飽和含水量以上の水がある場合はそれを下に落す。
! 
! ・注意
!     氷は下には落さない。
!     氷を落さねばならないような状況では、層の切り直しが必要。
!
      USE SIBPRM , ONLY :
     1   IDX      , 
     1   ISN  
! 
      USE SIBCON , ONLY : 
     1   HL_W2G_c        , 
     1   HL_W2G_INV_c    , 
     1   FR_SNOWWTR_MAX_c         ! 雪層において全 H2O 中で持てる液水の割合
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: KLEV                 ! 扱うレベル
!
      REAL(8),INTENT(IN)   :: H2O_SNOW_MIN          ! 旧雪 KLEV 層での下限
      REAL(8),INTENT(IN)   :: H2O_SNOW_MAX          ! 旧雪 KLEV 層での上限
!
      REAL(8),INTENT(IN)   :: ENG_FROM_UPPER (IDX)    
      REAL(8),INTENT(IN)   :: H2O_FROM_UPPER (IDX)    
      REAL(8),INTENT(IN)   :: ENG_SNOW       (IDX,ISN-1)
      REAL(8),INTENT(IN)   :: H2O_SNOW       (IDX,ISN-1)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)    :: H2O_TO_LOWER (IDX)    
      REAL(8),INTENT(INOUT)  :: H2O_TO_SOIL  (IDX)    
!
! INOUT
!
      INTEGER,INTENT(INOUT)  :: IFLG_SLEV_CHANGE (IDX) 
!
! WORK
!
      INTEGER  :: I
      REAL(8)  :: H2O_0
      REAL(8)  :: H2O_1
      REAL(8)  :: ENG_0
      REAL(8)  :: ENX_0
      REAL(8)  :: ENX_2
      REAL(8)  :: H2O_DEL
      REAL(8),PARAMETER  :: ZERO = 0. 
!
      DO 1000 I=NSTR,NEND 
      IF ( IFLG_SLEV_CHANGE(I).EQ.0 ) THEN  
!
!  エネルギー・水の和をとる
!
        ENG_0  = ENG_SNOW(I,KLEV) + ENG_FROM_UPPER(I) 
        H2O_0  = H2O_SNOW(I,KLEV) + H2O_FROM_UPPER(I)  
!
        ENX_0  = ENG_0 + HL_W2G_c * H2O_0
        ENX_2  = - ENG_0 * FR_SNOWWTR_MAX_c 
!
! その層が融ける場合。
! あるいは、その層からかなりの量の水が下に落ちる場合   更新フラグ A 
!
        IF     ( ENG_0.GT.0. .OR. ENX_0.GT.ENX_2*4 ) THEN
          IFLG_SLEV_CHANGE (I) = 3
!
! その層からある程度の量の水が下に落ちる場合   更新フラグ B
!
        ELSEIF ( ENX_0 .GT. ENX_2*2 ) THEN
          IFLG_SLEV_CHANGE (I) = 2
!
! その層からは水がちょっと下に落ちる、または落ちない。
!   --> 更新フラグは設定せず、そのまま進む。
!
!      落ちる   ... ENX_0 - ENX_2 > 0 の時。このとき ENX_1 = ENX_2 になる。
!      落ちない ... ENX_0 - ENX_2 < 0 の時。このとき ENX_1 = ENX_0 。
!
        ELSE
!
          H2O_DEL = MAX ( ENX_0 - ENX_2 , ZERO ) * HL_W2G_INV_c
          H2O_1 = H2O_0 - H2O_DEL
!
!   この H2O_1 が上下限に収まれば問題なし? 
!         H2O_SNOW , ENG_SNOW のゼロクリア
!         H2O_SNOW_NEXT      等への代入はしない。
!
          IF (       H2O_1 .GE. H2O_SNOW_MIN 
     1         .AND. H2O_1 .LE. H2O_SNOW_MAX ) THEN 
            H2O_TO_LOWER (I)   = H2O_DEL
            H2O_TO_SOIL  (I)   = H2O_DEL
!
!   この H2O_1 が上下限に収まらなければ、層の切り直しが必要。更新フラグ C 
!
          ELSE
            IFLG_SLEV_CHANGE (I) = 1  
          ENDIF
!
        ENDIF
      ENDIF
 1000 CONTINUE
!
      RETURN
      END

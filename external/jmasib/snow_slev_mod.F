! 旧雪層の切替え ～ 不良度が高いものについて
! 99/07/10 保坂征宏  第一版
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/31 保坂征宏 ENG_DROP 等の設定追加
!
      SUBROUTINE SNOW_SLEV_MOD (
     I    NMAX      , 
     I    ENG_PRECIPI  , H2O_PRECIPI  , RHO_PRECIPI_INV       , 
     I    ENG_SNOW     , H2O_SNOW     , RHO_SNOW_INV          , 
     M    ENG_SNOW_NEXT, H2O_SNOW_NEXT, RHO_SNOW_INV_NEXT     ,  
     M    ENH_BUCKET_SNOW_TO_SNSL     , WTR_BUCKET_SNOW_TO_SNSL   , 
     M    IFLG_SLEV_CHANGE            , INFO_SNOW_NEXT    )
!
! ・メモ
!   最初で水はしぼる。落とすものは全てここで落としてしまう。
!   あとは、質量に応じて層をとるだけ。ここでは水としては落ちないことに注意。
!
! ・注意
!   作業領域のゼロクリアがきちんとできていないかも知れない。要チェック。
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ISN   
!
      USE SIBCON , ONLY : 
     1    HL_W2G_c        , 
     1    H2O_SNOW_1_2a_c  , 
     1    H2O_SNOW_2a_2b_c , 
     1    H2O_SNOW_2b_3_c  , 
     1    H2O_SNOW1_STD_c  , 
     1    H2O_SNOW3_STD_c 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , INA ,
     1                       IMONIT_LEVEL , IMONIT_ALL ,
     1                       IMONIT_FLUX_MRI , IMONIT_ALMA
#endif
!
      IMPLICIT NONE    
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NMAX 
!
      REAL(8),INTENT(IN)   :: ENG_PRECIPI   (IDX) 
      REAL(8),INTENT(IN)   :: H2O_PRECIPI   (IDX) 
      REAL(8),INTENT(IN)   :: ENG_SNOW      (IDX,ISN) 
      REAL(8),INTENT(IN)   :: H2O_SNOW      (IDX,ISN) 
!
      REAL(8),INTENT(IN)   :: RHO_PRECIPI_INV (IDX)    
      REAL(8),INTENT(IN)   :: RHO_SNOW_INV    (IDX,ISN)
!
! INOUT 
!
!
      REAL(8),INTENT(INOUT) :: ENG_SNOW_NEXT  (IDX,ISN) 
      REAL(8),INTENT(INOUT) :: H2O_SNOW_NEXT  (IDX,ISN) 
      REAL(8),INTENT(INOUT) :: RHO_SNOW_INV_NEXT  (IDX,ISN)
!
      INTEGER,INTENT(INOUT) :: INFO_SNOW_NEXT (IDX)  ! 層数
      INTEGER,INTENT(INOUT) :: IFLG_SLEV_CHANGE (IDX)  ! フラグ
!
      REAL(8),INTENT(INOUT) :: ENH_BUCKET_SNOW_TO_SNSL (IDX) 
      REAL(8),INTENT(INOUT) :: WTR_BUCKET_SNOW_TO_SNSL (IDX) 
!
!
! WORK 
!
      INTEGER    :: I
!
      REAL(8)    :: WGT(0:3,1:3)
!
      REAL(8),PARAMETER :: DUMMY (IDX) = 0.
!
      REAL(8)    :: H2O_TO_1_NEXT  (IDX) 
      REAL(8)    :: H2O_TO_2_NEXT  (IDX) 
      REAL(8)    :: H2O_TO_3_NEXT  (IDX) 
!
      REAL(8)    :: H2O_SUM    (IDX)
!
      REAL(8)    :: H2O_TO_1   (IDX)
      REAL(8)    :: H2O_TO_2   (IDX)
      REAL(8)    :: H2O_TO_3   (IDX)
!
      REAL(8)    :: ENG_P_TO_1 (IDX) 
      REAL(8)    :: H2O_P_TO_1 (IDX) 
      REAL(8)    :: ENG_1_TO_2 (IDX) 
      REAL(8)    :: H2O_1_TO_2 (IDX) 
      REAL(8)    :: ENG_2_TO_3 (IDX) 
      REAL(8)    :: H2O_2_TO_3 (IDX) 
!
      REAL(8)    :: ENG_PRE_M (IDX)          ! しぼったあと。
      REAL(8)    :: H2O_PRE_M (IDX) 
      REAL(8)    :: ENG_SNOW_M (IDX,ISN) 
      REAL(8)    :: H2O_SNOW_M (IDX,ISN) 
!
      REAL(8)    :: BUNBO_P
      REAL(8)    :: BUNBO_1
      REAL(8)    :: BUNBO_2
      REAL(8)    :: BUNBO_3
!
      REAL(8),PARAMETER  :: SMALL   = 1.D-40
      REAL(8),PARAMETER  :: ZERO  = 0.D0
!
#ifdef MONYOS
      REAL(8),PARAMETER  :: ONE = 1.D0 
!     REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
#ifdef SIB_DEBUG
      REAL(8)   :: ENG_SNOW_SUM_SAVE  (IDX)
      REAL(8)   :: H2O_SNOW_SUM_SAVE  (IDX)
      REAL(8)   :: ENG_SNOW_SUM_TMP   (IDX)
      REAL(8)   :: H2O_SNOW_SUM_TMP   (IDX)
      REAL(8)   :: ENG_SNOW_SUM_NEW   (IDX)
      REAL(8)   :: H2O_SNOW_SUM_NEW   (IDX)
      REAL(8)   :: H2O_ERR_MAX
      REAL(8)   :: ENG_ERR_MAX
#endif
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_SLEV_MOD POCHAKA 99/11/02' 
        LFIRST = .FALSE.
      ENDIF
!
!  ==============
!  >> 暫定処理 <<
!  ==============
!   現状（991019）で、良質な層の切り分け直しでよい場合(1,2)についても
!   この方法、つまり相当たちの悪い場合でも何とかできる方法、をとっている。
!
      DO 100 I=1,NMAX
        IF (      IFLG_SLEV_CHANGE(I) .EQ. 1 
     1       .OR. IFLG_SLEV_CHANGE(I) .EQ. 2 ) THEN
          IFLG_SLEV_CHANGE(I) = 3
        ENDIF
 100  CONTINUE  
!
#ifdef SIB_DEBUG
      DO I=1,NMAX
        IF ( IFLG_SLEV_CHANGE(I) .EQ. 3 ) THEN
          H2O_SNOW_SUM_SAVE(I) =  H2O_PRECIPI (I) 
     1                            + H2O_SNOW(I,1)
     1                            + H2O_SNOW(I,2)
     1                            + H2O_SNOW(I,3)
        ENDIF
      ENDDO
#endif      
!
!  ======================
!  >> 上から順にしぼる <<
!  ======================
!    しぼる操作で RHO_SNOW は変わらないと仮定する。
!    保水量上限までは持てることにしたい。( 今は全てしぼって落している ) 
!
! 降水
!
      CALL SIBORU (
     I    1             , NMAX          , 
     I    IFLG_SLEV_CHANGE              , 3         , 
     I    ENG_PRECIPI   , H2O_PRECIPI   , 
     I    DUMMY         , DUMMY         ,  
     O    ENG_PRE_M     , H2O_PRE_M     ,        ! しぼった後
     O    ENG_P_TO_1    , H2O_P_TO_1      )      ! 落ちる量
!
! 第一層
!
      CALL SIBORU (
     I    1                  , NMAX               , 
     I    IFLG_SLEV_CHANGE                        , 3         , 
     I    ENG_SNOW(1,1)      , H2O_SNOW(1,1)      , 
     I    ENG_P_TO_1         , H2O_P_TO_1         ,  
     O    ENG_SNOW_M(1,1)    , H2O_SNOW_M(1,1) , 
     O    ENG_1_TO_2       , H2O_1_TO_2    ) 
!
! 第二層
!
      CALL SIBORU (
     I    1                  , NMAX               , 
     I    IFLG_SLEV_CHANGE                        , 3         , 
     I    ENG_SNOW  (1,2)    , H2O_SNOW(1,2)   , 
     I    ENG_1_TO_2         , H2O_1_TO_2       ,  
     O    ENG_SNOW_M(1,2)    , H2O_SNOW_M(1,2) , 
     O    ENG_2_TO_3         , H2O_2_TO_3    ) 
!
! 第三層
!
      CALL SIBORU (
     I    1                  , NMAX               , 
     I    IFLG_SLEV_CHANGE                        , 3         , 
     I    ENG_SNOW(1,3)      , H2O_SNOW(1,3)      , 
     I    ENG_2_TO_3         , H2O_2_TO_3         ,  
     O    ENG_SNOW_M(1,3)    , H2O_SNOW_M(1,3)    , 
     O    ENH_BUCKET_SNOW_TO_SNSL  , WTR_BUCKET_SNOW_TO_SNSL     ) 
!
#ifdef SIB_DEBUG
      DO I=1,NMAX
        IF ( IFLG_SLEV_CHANGE(I) .EQ. 3 ) THEN
          H2O_SNOW_SUM_TMP(I) =     H2O_PRE_M (I) 
     1                            + H2O_SNOW_M(I,1)
     1                            + H2O_SNOW_M(I,2)
     1                            + H2O_SNOW_M(I,3)
     1                            + WTR_BUCKET_SNOW_TO_SNSL(I) 
        ENDIF
      ENDDO
#endif      

!
!  ======================================
!  >> H2O の鉛直和を計算、層数を決める <<
!  ======================================
!
! ----------------
! > H2O の鉛直和 <
! ----------------
!
      DO 1000 I=1,NMAX
        IF ( IFLG_SLEV_CHANGE(I).EQ.3 ) THEN 
          H2O_TO_1 (I)  = H2O_PRE_M(I) + H2O_SNOW_M(I,1)    ! N までの和
          H2O_TO_2 (I)  = H2O_TO_1 (I) + H2O_SNOW_M(I,2)  
          H2O_TO_3 (I)  = H2O_TO_2 (I) + H2O_SNOW_M(I,3)  
          H2O_SUM  (I)  = H2O_TO_3 (I)
        ENDIF
 1000 CONTINUE    
!
! ----------------
! > 層の数を決定 <
! ----------------
!
      DO 1100 I=1,NMAX
!
        IF ( IFLG_SLEV_CHANGE(I).EQ.3 )  THEN
!
!    零層
          IF      ( H2O_SUM(I) .LE. 0. ) THEN
            INFO_SNOW_NEXT (I) = - 1
            H2O_TO_1_NEXT (I) = 0.
            H2O_TO_2_NEXT (I) = 0.
            H2O_TO_3_NEXT (I) = 0.
!
!    一層
          ELSEIF  ( H2O_SUM(I) .LT. H2O_SNOW_1_2a_c ) THEN
            INFO_SNOW_NEXT  (I) = 1
            H2O_TO_1_NEXT (I) = H2O_SUM(I) 
            H2O_TO_2_NEXT (I) = H2O_SUM(I) 
            H2O_TO_3_NEXT (I) = H2O_SUM(I) 
!  
!    二層、  ただし 1 と 2 とで層厚おなじ 
          ELSEIF ( H2O_SUM(I) .LT. H2O_SNOW_2a_2b_c ) THEN
            INFO_SNOW_NEXT  (I) = 2
            H2O_TO_1_NEXT (I) = H2O_SUM(I) * 0.5
            H2O_TO_2_NEXT (I) = H2O_SUM(I) 
            H2O_TO_3_NEXT (I) = H2O_SUM(I) 
!
!    二層、 
          ELSEIF ( H2O_SUM(I) .LT. H2O_SNOW_2b_3_c ) THEN
            INFO_SNOW_NEXT  (I) = 2
            H2O_TO_1_NEXT (I) = H2O_SNOW1_STD_c 
            H2O_TO_2_NEXT (I) = H2O_SUM(I) 
            H2O_TO_3_NEXT (I) = H2O_SUM(I) 
!
!    三層、 
          ELSE 
            INFO_SNOW_NEXT  (I) = 3
            H2O_TO_1_NEXT (I) = H2O_SNOW1_STD_c 
            H2O_TO_2_NEXT (I) = H2O_SUM(I) - H2O_SNOW3_STD_c
            H2O_TO_3_NEXT (I) = H2O_SUM(I) 
          ENDIF
        ENDIF
!
 1100 CONTINUE  
!
! ========================
! >> 割合テーブルの作成 <<
! ========================
!
! WGT(K_OLD,K_NEW) として、
! ENG(K_NEW) = SUM_K_OLD ENG(K_OLD) * WGT(K_OLD,KNEW)
!
      DO 2000 I=1,NMAX 
        IF ( IFLG_SLEV_CHANGE(I).EQ.3 )  THEN
          BUNBO_P = 1/MAX ( H2O_PRE_M(I)                , SMALL ) 
          BUNBO_1 = 1/MAX ( H2O_TO_1 (I) - H2O_PRE_M(I) , SMALL ) 
          BUNBO_2 = 1/MAX ( H2O_TO_2 (I) - H2O_TO_1 (I) , SMALL ) 
          BUNBO_3 = 1/MAX ( H2O_TO_3 (I) - H2O_TO_2 (I) , SMALL ) 
!
!  降水の分配
!
          WGT( 0,1)= MAX (    MIN(H2O_PRE_M(I),H2O_TO_1_NEXT(I))
     1                     ,  ZERO ) * BUNBO_P 
!
          WGT( 0,2)= MAX (    MIN(H2O_PRE_M(I),H2O_TO_2_NEXT(I))
     1                     -                   H2O_TO_1_NEXT(I)
     1                     ,  ZERO ) * BUNBO_P 
!
          WGT( 0,3)= MAX (    MIN(H2O_PRE_M(I),H2O_TO_3_NEXT(I))
     1                     -                   H2O_TO_2_NEXT(I)
     1                     ,  ZERO ) * BUNBO_P 
!
! 雪第一層の分配
!
          WGT(1,1) = MAX (    MIN(H2O_TO_1 (I),H2O_TO_1_NEXT(I))
     1                     -      H2O_PRE_M(I)
     1                     ,  ZERO ) * BUNBO_1 
          WGT(1,2) = MAX (    MIN(H2O_TO_1 (I),H2O_TO_2_NEXT(I))
     1                     -  MAX(H2O_PRE_M(I),H2O_TO_1_NEXT(I))
     1                     ,  ZERO ) * BUNBO_1 
          WGT(1,3) = MAX (    MIN(H2O_TO_1 (I),H2O_TO_3_NEXT(I))
     1                     -  MAX(H2O_PRE_M(I),H2O_TO_2_NEXT(I))
     1                     ,  ZERO ) * BUNBO_1 
!
! 雪第二層の分配
!
          WGT(2,1) = MAX (    MIN(H2O_TO_2(I),H2O_TO_1_NEXT(I))
     1                     -      H2O_TO_1(I)
     1                     ,  ZERO ) * BUNBO_2 
          WGT(2,2) = MAX (    MIN(H2O_TO_2(I),H2O_TO_2_NEXT(I))
     1                     -  MAX(H2O_TO_1(I),H2O_TO_1_NEXT(I))
     1                     ,  ZERO ) * BUNBO_2 
          WGT(2,3) = MAX (    MIN(H2O_TO_2(I),H2O_TO_3_NEXT(I))
     1                     -  MAX(H2O_TO_1(I),H2O_TO_2_NEXT(I))
     1                     ,  ZERO ) * BUNBO_2 
!
! 雪第三層の分配
!
          WGT(3,1) = MAX (    MIN(H2O_TO_3(I),H2O_TO_1_NEXT(I))
     1                     -      H2O_TO_2(I)
     1                     ,  ZERO ) * BUNBO_3 
          WGT(3,2) = MAX (    MIN(H2O_TO_3(I),H2O_TO_2_NEXT(I))
     1                     -  MAX(H2O_TO_2(I),H2O_TO_1_NEXT(I))
     1                     ,  ZERO ) * BUNBO_3 
          WGT(3,3) = MAX (    MIN(H2O_TO_3(I),H2O_TO_3_NEXT(I))
     1                     -  MAX(H2O_TO_2(I),H2O_TO_2_NEXT(I))
     1                     ,  ZERO ) * BUNBO_3 
!
! 各層の水・エネルギー・密度の計算
!
          H2O_SNOW_NEXT(I,1) = H2O_TO_1_NEXT (I)
          H2O_SNOW_NEXT(I,2) = H2O_TO_2_NEXT (I) - H2O_TO_1_NEXT (I)
          H2O_SNOW_NEXT(I,3) = H2O_TO_3_NEXT (I) - H2O_TO_2_NEXT (I)
!
          ENG_SNOW_NEXT(I,1) =
     1             ENG_PRE_M (I)   * WGT( 0,1)  
     1           + ENG_SNOW_M(I,1) * WGT( 1,1)   
     1           + ENG_SNOW_M(I,2) * WGT( 2,1)   
     1           + ENG_SNOW_M(I,3) * WGT( 3,1)   
!
          ENG_SNOW_NEXT(I,2) =
     1             ENG_PRE_M (I)   * WGT( 0,2)  
     1           + ENG_SNOW_M(I,1) * WGT( 1,2)   
     1           + ENG_SNOW_M(I,2) * WGT( 2,2)   
     1           + ENG_SNOW_M(I,3) * WGT( 3,2)   
!
          ENG_SNOW_NEXT(I,3) =
     1             ENG_PRE_M (I)   * WGT( 0,3)  
     1           + ENG_SNOW_M(I,1) * WGT( 1,3)   
     1           + ENG_SNOW_M(I,2) * WGT( 2,3)   
     1           + ENG_SNOW_M(I,3) * WGT( 3,3)   
!
          BUNBO_1 = 1/MAX (    WGT( 0,1) * H2O_PRE_M(I)       
     1                      +  WGT( 1,1) * H2O_SNOW_M(I,1)   
     1                      +  WGT( 2,1) * H2O_SNOW_M(I,2)   
     1                      +  WGT( 3,1) * H2O_SNOW_M(I,3) 
     1                      , SMALL )   
          BUNBO_2 = 1/MAX (    WGT( 0,2) * H2O_PRE_M(I)       
     1                      +  WGT( 1,2) * H2O_SNOW_M(I,1)   
     1                      +  WGT( 2,2) * H2O_SNOW_M(I,2)   
     1                      +  WGT( 3,2) * H2O_SNOW_M(I,3) 
     1                      , SMALL )   
          BUNBO_3 = 1/MAX (    WGT( 0,3) * H2O_PRE_M(I)       
     1                      +  WGT( 1,3) * H2O_SNOW_M(I,1)   
     1                      +  WGT( 2,3) * H2O_SNOW_M(I,2)   
     1                      +  WGT( 3,3) * H2O_SNOW_M(I,3) 
     1                      , SMALL )   
!
          RHO_SNOW_INV_NEXT(I,1) =
     1      (   RHO_PRECIPI_INV(I)   * WGT( 0,1) * H2O_PRE_M (I)       
     1        + RHO_SNOW_INV   (I,1) * WGT( 1,1) * H2O_SNOW_M(I,1)   
     1        + RHO_SNOW_INV   (I,2) * WGT( 2,1) * H2O_SNOW_M(I,2)     
     1        + RHO_SNOW_INV   (I,3) * WGT( 3,1) * H2O_SNOW_M(I,3) )      
     1       * BUNBO_1
!
          RHO_SNOW_INV_NEXT(I,2) =
     1      (   RHO_PRECIPI_INV(I)   * WGT( 0,2) * H2O_PRE_M(I)       
     1        + RHO_SNOW_INV   (I,1) * WGT( 1,2) * H2O_SNOW_M(I,1)   
     1        + RHO_SNOW_INV   (I,2) * WGT( 2,2) * H2O_SNOW_M(I,2)     
     1        + RHO_SNOW_INV   (I,3) * WGT( 3,2) * H2O_SNOW_M(I,3) )      
     1       * BUNBO_2
!
          RHO_SNOW_INV_NEXT(I,3) =
     1      (   RHO_PRECIPI_INV(I)   * WGT( 0,3) * H2O_PRE_M(I)       
     1        + RHO_SNOW_INV   (I,1) * WGT( 1,3) * H2O_SNOW_M(I,1)   
     1        + RHO_SNOW_INV   (I,2) * WGT( 2,3) * H2O_SNOW_M(I,2)     
     1        + RHO_SNOW_INV   (I,3) * WGT( 3,3) * H2O_SNOW_M(I,3) )      
     1       * BUNBO_3
        ENDIF
 2000 CONTINUE   
!
!!    rho_snow_inv_next(:,:) = 1.D0/200.D0 
!
#ifdef SIB_DEBUG
      DO I=1,NMAX
        IF ( IFLG_SLEV_CHANGE(I) .EQ. 3 ) THEN
          H2O_SNOW_SUM_NEW(I) =   + H2O_SNOW_NEXT(I,1)
     1                            + H2O_SNOW_NEXT(I,2)
     1                            + H2O_SNOW_NEXT(I,3)
     1                            + WTR_BUCKET_SNOW_TO_SNSL(I) 
!
          IF ( 
     1         ABS(H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I)).GT.1.D-9 
     1       ) THEN
            write(6,*) 'SNOW_SLEV_MOD ERROR SIBORU ' , I , 
     1         H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_SAVE(I) , 
     1         H2O_SNOW_SUM_SAVE(I),H2O_SNOW_SUM_TMP(I)  
          ENDIF 
!
          IF ( 
     1         ABS(H2O_SNOW_SUM_TMP(I)-H2O_SNOW_SUM_NEW(I)).GT.1.D-9 
     1       ) THEN
            write(6,*) 'SNOW_SLEV_MOD ERROR BUNPAI ' , I , 
     1         H2O_SNOW_SUM_NEW(I)-H2O_SNOW_SUM_TMP(I) , 
     1         H2O_SNOW_SUM_TMP(I),H2O_SNOW_SUM_NEW(I)  
!
            write(6,*) 'TMP ' , 
     1                              H2O_PRE_M (I) ,
     1                            + H2O_SNOW_M(I,1),
     1                            + H2O_SNOW_M(I,2),
     1                            + H2O_SNOW_M(I,3),
     1                            + WTR_BUCKET_SNOW_TO_SNSL(I) 
            write(6,*) 'NEW ' , 
     1           H2O_SNOW_NEXT(I,1) ,
     1           H2O_SNOW_NEXT(I,2) , 
     1           H2O_SNOW_NEXT(I,3) , 
     1           WTR_BUCKET_SNOW_TO_SNSL(I)

          ENDIF 

        ENDIF
      ENDDO
#endif      

!
#ifdef MONYOS
!
! 格子平均フラックス量
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_MRI ) THEN
! 融解雪量 =  水バケツ （雪 -> 土壌） 
        CMARK = 'SNMLT'
        CALL MONIT_SIB ( CMARK , WTR_BUCKET_SNOW_TO_SNSL,INA, 
     1                   ONE , ONE )                            
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 水バケツ （雪 -> 土壌）
!
        CMARK = 'WFLBWL'
        CALL MONIT_SIB ( CMARK , WTR_BUCKET_SNOW_TO_SNSL,INA, 
     1                   ONE , ONE )                            
!
! 水移動に伴うエネルギー移動 （雪 -> 土壌）
!
        CMARK  = 'IFLBWL' 
        CALL MONIT_SIB ( CMARK , ENH_BUCKET_SNOW_TO_SNSL,INA, 
     1                   ONE , ONE )                               ! ver2

      ENDIF
!
! * ALMA
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          CMARK  = 'A_QST'         ! WFLBWL に同じ
          CALL MONIT_SIB ( CMARK , WTR_BUCKET_SNOW_TO_SNSL,INA, 
     1                   ONE , ONE )              
      ENDIF              

#endif
!
      RETURN
      END
!==========================================================
      SUBROUTINE SIBORU (
     I    NSTR      , NEND       , 
     I    IFLG_SLEV_CHANGE       , IFLG_VAL   , 
     I    ENG_ORG   , H2O_ORG    , 
     I    ENG_ADD   , H2O_ADD    ,             
     O    ENG_MOD   , H2O_MOD    ,                ! 変更後
     O    ENG_DROP  , H2O_DROP   )                ! 下に落ちる分
!
      USE SIBPRM , ONLY :
     1    IDX    
!
      USE SIBCON , ONLY : 
     1    HL_W2G_c    ,
     1    HL_W2G_INV_c    ,
     1    FR_SNOWWTR_MAX_c
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: IFLG_SLEV_CHANGE(IDX)
      INTEGER,INTENT(IN)    :: IFLG_VAL
      REAL(8),INTENT(IN)    :: ENG_ORG(IDX)
      REAL(8),INTENT(IN)    :: ENG_ADD(IDX)
      REAL(8),INTENT(IN)    :: H2O_ORG(IDX)
      REAL(8),INTENT(IN)    :: H2O_ADD(IDX)
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT) :: ENG_MOD (IDX)
      REAL(8),INTENT(INOUT) :: ENG_DROP(IDX)
      REAL(8),INTENT(INOUT) :: H2O_MOD (IDX)
      REAL(8),INTENT(INOUT) :: H2O_DROP(IDX)
!
! WORK
!
      INTEGER  :: I
      REAL(8)  :: ENG_1
      REAL(8)  :: H2O_1
      REAL(8)  :: ENX_1
      REAL(8)  :: ENX_2
!
! ここは条件を満たすところだけの計算にしたい
!
      DO 1000 I = NSTR,NEND

        IF ( IFLG_SLEV_CHANGE(I) .EQ. IFLG_VAL ) THEN 
!
          ENG_1 = ENG_ORG(I) + ENG_ADD(I)
          H2O_1 = H2O_ORG(I) + H2O_ADD(I)
!
          ENX_1 = ENG_1 + HL_W2G_c * H2O_1  
          ENX_2 = - ENG_1 * FR_SNOWWTR_MAX_c
!
          IF     ( ENG_1 .GE.0. ) THEN     ! 全部融ける
            ENG_DROP(I) = ENG_1
            H2O_DROP(I) = H2O_1
!
          ELSEIF ( ENX_1 .LE. ENX_2 ) THEN 
            ENG_DROP(I) = 0. 
            H2O_DROP(I) = 0. 
!
          ELSE
            ENG_DROP(I) = 0. 
            H2O_DROP(I) = ( ENX_1 - ENX_2 ) * HL_W2G_INV_c
!
          ENDIF
!
          ENG_MOD(I) = ENG_1 - ENG_DROP(I)
          H2O_MOD(I) = H2O_1 - H2O_DROP(I)
!
        ENDIF 
 1000 CONTINUE  
!
      RETURN
      END
!==============================================================
      SUBROUTINE MR_SNOW_SLEV_MOD
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL ,
     1                       IMONIT_FLUX_MRI  , IMONIT_ALMA 
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
        WRITE(6,*) 'MONIT_REGIST_SNOW_SLEV_MOD 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_MRI ) THEN
!   融雪量 = WFLBWL  
        CMARK  = 'SNMLT'
        CTITLE = 'Snow Melting Rate'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'WFLBWL'
      CTITLE = 'WATER FROM SNOW TO BS'
      CUNIT  = 'KG/M**2/S'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IFLBWL'
      CTITLE = 'ENG WITH WATER FROM SNOW TO BS'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_QST'         ! WFLBLW と同じ
        CTITLE = 'ALMA Qst : Water from Snow to Soil'
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF

#endif
!
      RETURN
      END



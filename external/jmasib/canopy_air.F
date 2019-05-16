! サブルーチン CUT  キャノピー空気状況の計算
! 99/06/0x 大泉
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_AIR ( 
     I  NSTR           , NEND           , INSFLAG      , 
     I  ICHECK         , 
     I  FLAG_DRY_CNP   , FLAG_DRY_SKIN  , 
     I  CVR_VEG_msx    , CVR_VEG_INV_msx, 
     I  R_AC_BAR       , R_AC_BAR_INV   , 
     I  R_AG           , R_AG_INV       , 
     I  R_STM_BAR      , R_BB           , R_AIR_INV    , 
     I  XFR_H2O_CNP    , XFR_H2O_GRS    , 
     I  E_CNP_SAT      , E_REF          , E_SKIN_SAT    ,
     I  FAC_BSH        , 
     O  E_AIR          , 
     O  R_ACX_BAR      , R_AGX_BAR      , 
     O  R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O  R_AB_INV   )
!
!=================================================================
! 蒸発関連の変数  E_AIR, R_ACX_BAR, R_AGX_BAR の計算 + RCP の計算
!================================================================
!
!  していること
!    具体的には
!      キャノピー空間での水蒸気圧             E_AIR 
!      キャノピーの葉の気孔抵抗 
!          + 葉の表面とキャノピー空間間の抵抗 R_ACX_BAR 
!      下草      の葉の気孔抵抗 
!          + 葉の表面とキャノピー空間間の抵抗 R_AGX_BAR
!    を求める。
!
!    これらは、キャノピー空間が過飽和かどうかで状況が変わり、値も大きく変わる。
!    過飽和かどうかのフラグ 
!       FLAG_DRY_CNP,FLAG_DRY_SKIN もここで確認されることに注意。
!
!  プログラマーへの注意
!
!  ・雪あり、雪無しとも、これで対応。
!
!  ・ R_AGX_BAR は以前の RG とは違うので注意。
!     R_AGX_BAR = RG + RD にあたる。 R_STM_BAR + R_AG 
!
!  ・ COC, COG1, COG2 について
!       COC  キャノピーからの蒸発・蒸散の並列抵抗の逆数。
!       COG1 地面・下草からの蒸発・蒸散の並列抵抗の逆数。
!       COG2 地面・下草からの蒸発・蒸散の並列抵抗の逆数、
!            ただし fac_bsh はかけない。
!
!  ・バグ
!      雪での下草?? の処理。どうするかを決めてから変形すること。
!
      USE SIBPRM , ONLY :
     1   IDX  , 
     1   ICG 
!
      IMPLICIT NONE
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT) :: E_AIR (IDX)   ! キャノピー空間の水蒸気圧
      REAL(8),INTENT(INOUT) :: R_ACX_BAR   (IDX)  
!                                             ! キャノピー気孔抵抗+境界層抵抗
      REAL(8),INTENT(INOUT) :: R_ACX_BAR_INV (IDX) 
      REAL(8),INTENT(INOUT) :: R_AGX_BAR   (IDX)   ! 下草気孔抵抗+境界層抵抗
      REAL(8),INTENT(INOUT) :: R_AGX_BAR_INV (IDX)
      REAL(8),INTENT(INOUT) :: R_AB_INV (IDX)
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: INSFLAG
      INTEGER,INTENT(IN)   :: ICHECK    (IDX)
      REAL(8),INTENT(IN)   :: FLAG_DRY_CNP   (IDX)  
      REAL(8),INTENT(IN)   :: FLAG_DRY_SKIN  (IDX)  
!
      REAL(8),INTENT(IN)   :: XFR_H2O_CNP(IDX)      ! キャの濡れた部分の割合
      REAL(8),INTENT(IN)   :: XFR_H2O_GRS(IDX)      ! 下草の濡れた部分の割合
!
      REAL(8),INTENT(IN)   :: R_STM_BAR (IDX,ICG)  ! 気孔抵抗
!
      REAL(8),INTENT(IN)   :: R_AC_BAR     (IDX)   ! キャ・キャ空間間
      REAL(8),INTENT(IN)   :: R_AC_BAR_INV (IDX)   !   境界層抵抗 
      REAL(8),INTENT(IN)   :: R_AG         (IDX)   ! 下草・地面とキャ空間間の
      REAL(8),INTENT(IN)   :: R_AG_INV     (IDX)   !   境界層抵抗 RD 
! 
      REAL(8),INTENT(IN)   :: R_BB         (IDX)   ! 裸地からの蒸発抵抗RSOIL
!
      REAL(8),INTENT(IN)   :: CVR_VEG_msx     (IDX,ICG)  ! 植生被覆率
      REAL(8),INTENT(IN)   :: CVR_VEG_INV_msx (IDX,ICG)  ! 植生被覆率
!
      REAL(8),INTENT(IN)   :: FAC_BSH   (IDX)      ! f_h ( 式 .. )土壌相対湿度
!
      REAL(8),INTENT(IN)   :: R_AIR_INV (IDX)      ! キャ空間と大気最下層間抵抗
!
      REAL(8),INTENT(IN)   :: E_CNP_SAT   (IDX)      ! TC での飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_SKIN_SAT  (IDX)      ! TG での飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_REF       (IDX)      ! 参照レベルでの水蒸気圧
!
! WORK 
!
      REAL(8)   :: COC  (IDX) 
      REAL(8)   :: COG1 (IDX) 
      REAL(8)   :: COG2 (IDX) 
!
      INTEGER   :: I 
      REAL(8)   :: D2
      REAL(8)   :: TOP
      REAL(8)   :: XNUM
!
!   ---------------
!   >> RC の計算 <<
!   ---------------
!
      DO 100 I = NSTR,NEND
      IF ( ICHECK(I).EQ.1 ) THEN
!
!  キャノピー結露せず
!
        IF ( FLAG_DRY_CNP(I).GT.0.5D0 ) THEN  
          R_ACX_BAR    (I) 
!Z   1       = R_STM_BAR(I,1) + 2.0D0 * R_AC_BAR(I)
     1       = R_STM_BAR(I,1) + R_AC_BAR(I)
!
          R_ACX_BAR_INV(I) 
     1       = 1. / R_ACX_BAR (I) 
!
          COC (I) 
     1       =   (1.0D0 -XFR_H2O_CNP(I)) * R_ACX_BAR_INV(I) 
!Z   1          +   XFR_H2O_CNP(I)       * R_AC_BAR_INV (I) * 0.5
     1          +   XFR_H2O_CNP(I)       * R_AC_BAR_INV (I) 
!
!  キャノピー結露
!
        ELSE                        
!
          R_ACX_BAR     (I) = 1.D30
          R_ACX_BAR_INV (I) = 0. 
!Z        COC  (I)          = R_AC_BAR_INV (I) * 0.5 
          COC  (I)          = R_AC_BAR_INV (I) 
!
        END IF  
      ENDIF
 100  CONTINUE   
!
! =======================
! >> R_AGX_BAR の計算 <<
! =======================
!          結露の意味は、全面が水または雪で覆われるという意味。
!          雪ならば結露とおなじことに注意。
!  ------------
!  > 雪無し域 <
!  ------------
!
      IF ( INSFLAG.EQ.0 ) THEN       
!
        DO 200 I = NSTR,NEND
        IF ( ICHECK(I).EQ.1 ) THEN
!
!   地面結露無し
!
!!        IF ( FLAG_DRY_SKIN(I).GT.0.5D0 .AND. INSFLAG.EQ.1 ) THEN  
          IF ( FLAG_DRY_SKIN(I).GT.0.5D0 ) THEN                    !010725
            XNUM = XFR_H2O_GRS(I) * R_AG_INV(I) 
     1                            * CVR_VEG_msx(I,2)            ! 下草蒸発
!
            R_AGX_BAR(I) 
     1           = R_STM_BAR (I,2) 
     1             + R_AG(I) * CVR_VEG_INV_msx(I,2)
!
            R_AGX_BAR_INV (I) 
     1           = 1 / R_AGX_BAR(I) 
!
            R_AB_INV (I) 
     1           = 1 / ( R_BB(I)+R_AG(I) ) 
!
            COG2 (I)     
     1       =   XNUM                                          ! 下草蒸発
     1         + (1.0D0 - XFR_H2O_GRS(I)) * R_AGX_BAR_INV (I)  ! 下草蒸散
     1         + ( 1.0D0 - CVR_VEG_msx(I,2))  * R_AB_INV(I)    ! 裸地蒸発
!
            COG1(I)      
     1         =
     1          +  XNUM                                        ! 下草蒸発
     1          +  (1.0D0 - XFR_H2O_GRS(I)) * R_AGX_BAR_INV(I) ! 下草蒸散
     1          + ( 1.0D0 - CVR_VEG_msx(I,2)) * R_AB_INV(I)    ! 裸地蒸発
     1                                        * FAC_BSH (I)
!
!   地面結露 
!
          ELSE      
!
!Z          R_AB_INV      (I) = 1 / (   R_BB(I)+R_AG(I) ) 
            R_AB_INV      (I) = 1 / ( 0*R_BB(I)+R_AG(I) ) 
!
            COG1(I) =  CVR_VEG_msx(I,2)        * R_AB_INV(I)    ! 下草蒸発
     1              + (1.0D0-CVR_VEG_msx(I,2)) * R_AB_INV(I)    ! 裸地蒸発
!!   1                                         * FAC_BSH(I)     ! ketsuro
!                                                          
            COG2(I) =  CVR_VEG_msx(I,2)        * R_AB_INV(I)    ! 下草蒸発
     1              + (1.0D0-CVR_VEG_msx(I,2)) * R_AB_INV(I)    ! 裸地蒸発
!
            R_AGX_BAR     (I) = 1.D30
            R_AGX_BAR_INV (I) = 0.
          END IF                                         
        ENDIF
 200    CONTINUE
!
!  ------------
!  > 雪あり域 <
!  ------------
!      雪無しの場合の結露条件とおなじ。かつ、R_BB = 0
!
      ELSE                                  ! 雪あり
        DO 300 I = NSTR,NEND
        IF ( ICHECK(I).EQ.1 ) THEN
          COG1(I) = R_AG_INV(I)                                 ! 雪面蒸発
          COG2(I) = COG1(I)
!
          R_AGX_BAR     (I) = 1.D30
          R_AGX_BAR_INV (I) = 0.
          R_AB_INV      (I) = R_AG_INV(I) 
        ENDIF
 300    CONTINUE  
!
      ENDIF
!
! ====================
! >> E_AIR の計算 <<
! ====================
!
      DO 400 I = NSTR,NEND
      IF ( ICHECK(I).EQ.1 ) THEN
        D2     = R_AIR_INV(I)  + COC(I) + COG2(I)
        TOP    =   COC(I)   * E_CNP_SAT(I) 
     1           + E_REF(I) * R_AIR_INV(I) 
     1           + COG1(I)  * E_SKIN_SAT(I)
        E_AIR (I) = TOP / D2
      END IF
  400 CONTINUE
!
      RETURN
      END

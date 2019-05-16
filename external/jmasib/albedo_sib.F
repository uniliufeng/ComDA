! 放射伝達方程式を解く
! 99/05/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE ALBEDO_SIB  ( 
     I   MON         , 
     I   NSTR        , NEND       , INSFLAG    , 
     I   ALB_SNOW_onCNP           , FR_GLA_CNP , 
     I   ALB_SKIN    , NTYP       ,  
     I   CVR_VEG_ms  ,  
     I   DEP_SNOW    , ZMEAN       , ZLT , 
     O   EXTK        , ABSORB_SHRT   , 
     O   REACH_GRD   , ALB_SIB        )
!
! 陸面の
!    短波アルベドと、
!    SiB 上端に 1 入射する短波の、キャノピー・地面での吸収率
! とを評価する。
!
! ただし雪あり域と雪無し域とは別々に呼ばれる。
!
      USE SIBPRM , ONLY : 
     1    IDX    ,  
     1    ITYP   , 
     1    IMON   , 
     1    ICG    , 
     1    IVN    , 
     1    ILD    , 
     1    IBD    , 
     1    IBS  
!
      USE SIBVEG , ONLY :
     1            TRAN_co               ,   ! 透過特性 (ITYP,ICG,IWV,ILD)
     1            REF_co                ,   ! 反射特性 (ITYP,ICG,IWV,ILD)
     1            CHI_L_co              ,   ! χL  (ITYP,ICG) 
     1            Z2_mo                 ,   ! キャノピー上端高さ (ITYP,IMON)
     1            Z1_mo                 ,   ! キャノピー下端高さ (ITYP,IMON)
     1            XLAI_mo               ,   ! LAI, 格子平均 (ITYP,IMON,ICG)
     1            GREEN_mo                  ! 生きている割合 (ITYP,IMON,ICG) 
!
      IMPLICIT NONE 
!
! INPUT START 
      INTEGER    NSTR                      ! 陸面格子数
      INTEGER    NEND                      ! 陸面格子数
      INTEGER    MON                       ! 月
      INTEGER,INTENT(IN)   :: INSFLAG      ! 地面上雪無し 0 雪あり 1 
!
      REAL(8),INTENT(IN)   :: DEP_SNOW(IDX)            ! 葉のもつ水
      REAL(8),INTENT(IN)   :: ZMEAN   (IDX)            ! cos ( 天頂角 ) = μ
      INTEGER,INTENT(IN)   :: NTYP    (IDX)            ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: CVR_VEG_ms(IDX,ICG)      ! 植生被覆率 
      REAL(8),INTENT(IN)   :: ALB_SKIN       (IDX,IVN) ! 裸地 or 雪面アルベド
      REAL(8),INTENT(IN)   :: ALB_SNOW_onCNP (IDX,IVN) !  木の葉の雪アルベド
      REAL(8),INTENT(IN)   :: FR_GLA_CNP ( IDX )       ! キャノピーの雪の割合
!
! OUTPUT  START 
!
      REAL(8),INTENT(INOUT)  :: 
     1        EXTK        (IDX,ICG,IVN,IBD)       ,  ! 消散係数 
     1        ABSORB_SHRT (IDX,ICG,IVN,IBD)       ,  ! 地上下向短波率
     1        ALB_SIB     (IDX,    IVN,IBD)       ,  ! アルベド
     1        REACH_GRD   (IDX,        IBD)          ! 地面到達率
!
      REAL(8),intent(in) :: ZLT   (IDX,ICG)           ! 葉面積指数
!
! WORK START 
!
      REAL(8)  ::   ALBEDO (IDX,ICG,IVN,IBD)         !CG・L・BD別アルベド  
!
      REAL(8)  ::   RMEW    (IDX)     ! μ                   F
      REAL(8)  ::   RMEWI   (IDX)  
      REAL(8)  ::   RMEWBAR (IDX)     ! μbar d              ZMEW
      REAL(8)  ::   RMEWBARI(IDX)            
      REAL(8)  ::   GG      (IDX)     ! G(μ)                PROJ
      REAL(8)  ::   RKK     (IDX)     ! K = G / μ           EXTKB
      REAL(8)  ::   HH    (IDX)       ! h                    PSI
      REAL(8)  ::   OMEGA (IDX)       ! ω                   SCAT
      REAL(8)  ::   BETOMG(IDX)       ! βω                 UPSCAT
      REAL(8)  ::   BETAO (IDX)       ! β_0                 UPSCAT
!
      REAL(8)  ::   SS    (IDX)       ! s   ☆
      REAL(8)  ::   TT    (IDX)       ! t   ☆
      REAL(8)  ::   TTI   (IDX)       ! 1/t   ☆
      REAL(8)  ::   UU    (IDX)       ! u   ☆
      REAL(8)  ::   VV    (IDX)       ! v   ☆
      REAL(8)  ::   STK   (IDX)       ! s^2-t^2-K^2 ☆
      REAL(8)  ::   STKI  (IDX)       ! 1/ (s^2-t^2-K^2) ☆
      REAL(8)  ::   ALP1  (IDX)       ! α1 ☆
      REAL(8)  ::   ALP2  (IDX)
      REAL(8)  ::   ALP3  (IDX)
      REAL(8)  ::   ALP4  (IDX)
      REAL(8)  ::   ALP5  (IDX)
      REAL(8)  ::   ALP6  (IDX)
      REAL(8)  ::   A23X  (IDX)
      REAL(8)  ::   A23Y  (IDX)
      REAL(8)  ::   TRA1 (IDX)  ! 可視直達、大気下端下向きで葉に達したうち
!                             ! 葉を直達のまま通り抜ける割合
!                             ! 葉での反射は考慮。
      REAL(8)  ::   TRA2 (IDX)  ! 可視散乱、大気下端下向きで葉に達したうち
!                             ! 葉を通り抜ける割合
      REAL(8)  ::   TRA3 (IDX)  ! 可視直達、大気下端下向きで葉に達したうち
!                             ! 葉で散乱されて通り抜ける割合
!     REAL(8)  ::   EKL   ( IDX )            ! exp(KL)           
      REAL(8)  ::   EHL   ( IDX )            ! exp(hL)            
      REAL(8)  ::   EKLI  ( IDX )            ! exp(-KL)     EK  
      REAL(8)  ::   EHLI  ( IDX )            ! exp(-hL)     EPSI
      REAL(8)  ::   CHIV  ( IDX ) 
      REAL(8)  ::   REFF1 ( IDX )
      REAL(8)  ::   REFF2 ( IDX )
      REAL(8)  ::   TRAN1 ( IDX )
      REAL(8)  ::   TRAN2 ( IDX )
      REAL(8)  ::   AA    ( IDX )
      REAL(8)  ::   BB    ( IDX )
      REAL(8)  ::   ACSS  ( IDX )
      REAL(8)  ::   FAC   ( IDX )   ! 雪でキャノピーが減る割合
      REAL(8)  ::   RHL   ( IDX )
      REAL(8)  ::   RKL   ( IDX )
      REAL(8)  ::   ROSB  ( IDX )
      REAL(8)  ::   ROSD  ( IDX )
      REAL(8)  ::   ZAT   ( IDX )
      REAL(8)  ::   SCOV1 ( IDX )
!
      INTEGER ::   I , JVN , IRD  , IVEG
!
      REAL(8),PARAMETER ::   DZERO  = 0.D0
      REAL(8),PARAMETER ::   D095   = 0.95D0
      REAL(8),PARAMETER ::   D05    = 0.5D0
      LOGICAL,SAVE ::   LFIRST
!
      DATA LFIRST  / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'ALBEDO_SIB 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF
!
!  ========================
!  >>> 入力データ微修正 <<<
!  ========================
!
!   --------------------
!  >> 天頂角下限処理 <<
!  --------------------
!                                                     ! 収支合うの？
      DO 1000 I=NSTR,NEND
        RMEW(I)  = MAX ( ZMEAN(I), 0.01746D0 )         ! 天頂角 μ
        RMEWI(I) = 1. / RMEW(I)
 1000 CONTINUE
!
!  -----------------------------
!  > 雪がある場合の LAI の補正 <
!  -----------------------------
!       注。雪無し域計算では DEP_SNOW = 0 
!           上限値は 0.95 にしている。
!
      DO 1300 I=NSTR,NEND
        FAC(I)= ( DEP_SNOW(I)         - Z1_mo(NTYP(I),MON) )
     1          / ( Z2_mo(NTYP(I),MON)- Z1_mo(NTYP(I),MON) ) 
        FAC(I)= MAX ( DZERO , FAC(I) )
        FAC(I)= MIN ( D095  , FAC(I) )
 1300 CONTINUE
!
!  =====================================
!! >> 波長別処理ここから ( 最後まで ) <<
!  =====================================
!
#include "noparallel_fixedform"
      DO 2000 JVN = 1, 2
!
!  =====================================================
!! >> キャノピー(IVEG=2)・地面(IVEG=1) 別処理ここから <<
!  =====================================================
!!        注意 : 地面から先に計算する。
!
#include "noparallel_fixedform"
      DO 2100 IVEG =2, 1, -1
!
!  --------------------------
!! >>> 散乱係数 ω の計算 <<<
!  --------------------------
!
      IF ( IVEG.EQ.1 ) THEN    ! 木の上
!
      DO 3000 I=NSTR,NEND
!
        SCOV1 (I) = MIN ( D05 , FR_GLA_CNP(I) )  ! 木の葉の半分は見えている 
!
!   反射率 α の計算 
!
        REFF1(I)=                                ! 生きている葉   α生
     1           ( 1.0D0-SCOV1(I) ) 
     1                   * REF_co (NTYP(I),IVEG,JVN,1)  
     1              + SCOV1(I)      
     1                   * ALB_SNOW_onCNP(I,JVN)                  
!
        REFF2(I)=                                  ! 死んだ葉       α死
     1           ( 1.0D0-SCOV1(I) ) 
     1                   * REF_co (NTYP(I),IVEG,JVN,2) 
     1              + SCOV1(I)      
     1                   * ALB_SNOW_onCNP(I,JVN)                  
!
! 
!  透過率 δ の計算 
!
        TRAN1(I) =                                    ! 生きた葉 δ生
     1       (  ( 1.0D0 - SCOV1(I) )
     1           +  SCOV1(I)
     1                * ( 1.0D0-ALB_SNOW_onCNP(I,JVN) ) ) 
     1       * TRAN_co(NTYP(I),IVEG,JVN,1)           
!
        TRAN2(I) =                                    ! 死んだ葉 δ死
     1       (  ( 1.0D0 - SCOV1(I) )
     1        + SCOV1(I)
     1                * ( 1.0D0-ALB_SNOW_onCNP(I,JVN))*0.9D0)      
     1       * TRAN_co(NTYP(I),IVEG,JVN,2)              
!                                                 ! *0.9 は不明
! 
!  散乱係数 ω
!
        OMEGA(I) =                                      ! ω
     1        GREEN_mo(NTYP(I),MON,IVEG)      * ( TRAN1(I) + REFF1(I) )       
     1   + (1.0D0-GREEN_mo(NTYP(I),MON,IVEG)) * ( TRAN2(I) + REFF2(I) )
!
 3000 CONTINUE
!
      ELSEIF ( IVEG.EQ.2 .AND. INSFLAG.EQ.0 ) THEN   ! 地面で雪無し
!
!                   旧radalb DO LOOP 440 で SCOV1=0 を与えた式にする
!
        DO 3001 I=NSTR,NEND
!
!   反射率 α の計算 
!
          REFF1(I)= REF_co (NTYP(I),IVEG,JVN,1)     ! 生きている葉   α生
          REFF2(I)= REF_co (NTYP(I),IVEG,JVN,2)     ! 死んだ葉       α死
! 
!  透過率 δ の計算  
!
          TRAN1(I) =  TRAN_co(NTYP(I),IVEG,JVN,1)   ! 生きた葉 δ生
!         TRAN2(I) =  TRAN_co(NTYP(I),IVEG,JVN,1)   ! 死んだ葉 δ死 
          TRAN2(I) =  TRAN_co(NTYP(I),IVEG,JVN,2)   ! 死んだ葉 δ死
!
!  散乱係数 ω
!
          OMEGA(I) =                                      ! ω
     1       GREEN_mo(NTYP(I),MON,IVEG)       * ( TRAN1(I) + REFF1(I) )
     1   + (1.0D0-GREEN_mo(NTYP(I),MON,IVEG)) * ( TRAN2(I) + REFF2(I) )       
!
 3001   CONTINUE
      ENDIF          
!
!  ここから先、雪の地面計算はしない   ------------------ AAAA 
!
      IF ( .NOT. ( IVEG.EQ.2 .AND. INSFLAG.EQ.1 ) ) THEN   ! 雪面計算★以外★
!
!   -------------------------------
!!  >> 下面のアルベド ROSB, ROSD <<
!   -------------------------------
!
      IF     ( IVEG.EQ.2 .AND. INSFLAG.EQ.0 ) THEN   ! 土の地面
        DO 3010 I=NSTR,NEND
!         ZAT (I) = XLAI_mo(NTYP(I),MON,IVEG) 
          ZAT (I) = ZLT(I,IVEG) 
     1                      / CVR_VEG_ms(I,IVEG)     ! L = Lg/Vg
          ROSB(I) = ALB_SKIN (I,JVN)                ! 裸地・直達光アルベドA^b
          ROSD(I) = ALB_SKIN (I,JVN)                ! 裸地・直達光アルベドA^b
 3010   CONTINUE
      ELSEIF ( IVEG.EQ.1 ) THEN                      ! キャノピー
        DO 3020 I=NSTR,NEND
          ZAT (I) = ZLT(I,IVEG) 
!         ZAT (I) = XLAI_mo(NTYP(I),MON,IVEG) 
     1                      / CVR_VEG_ms(I,IVEG)     ! L = Lc/Vc
!         ZAT (I) = ZAT(I) * (1.0D0 -FAC(I))         ! L ( 1-FAC )
          ROSB(I) = ALBEDO(I,2,JVN,1)              ! 直達光  A^μ
          ROSD(I) = ALBEDO(I,2,JVN,2)              ! 散乱光  A^μ
 3020  CONTINUE
      END IF
!
!    ---------------------------
!!   >> 葉の向きのパラメタχL <<
!    ---------------------------
!
      DO 3100 I=NSTR,NEND
        
        CHIV(I) = CHI_L_co(NTYP(I),IVEG)
        AA(I) = 0.500D0 - 0.633D0 * CHIV(I)- 0.33D0 * CHIV(I) * CHIV(I) 
        BB(I) = 0.877D0 - 1.754D0 * AA(I)                               
!
        RMEWBAR (I) = 1.0D0
     1             / BB(I)
     1             * (1.0D0 -AA(I)/BB(I)*LOG((AA(I)+BB(I))/AA(I)))   !μbar
!!         ただし  式とは異なる？
        RMEWBARI(I) = 1.D0 / RMEWBAR(I)
 3100 CONTINUE
!
      DO 3200 I=NSTR,NEND
        GG   (I) = AA(I) + BB(I) * RMEW(I)    ! 方位分布 G(μ)   ( F は μ )
        RKK  (I) = GG(I) * RMEWI(I)           ! 光学的深さ K = G / μ
 3200 CONTINUE
!
!  --------------------------
!  > 消散係数の計算（出力）<
!  --------------------------
!
      DO 3210 I=NSTR,NEND
        EXTK(I, IVEG, JVN, 1 ) 
     1        = GG(I) * SQRT(1.0D0-OMEGA(I)) * RMEWI(I)    ! G*sqrt(1-ω)/μ
        EXTK(I, IVEG, JVN, 2 ) 
     1        =         SQRT(1.0D0-OMEGA(I)) * RMEWBARI(I) !sqrt(1-ω)/μbar
 3210 CONTINUE
!
      DO 3220 I=NSTR,NEND
        ACSS (I) = OMEGA(I)* 0.5D0 *GG(I)/( GG(I)+RMEW(I)*BB(I) )
        ACSS (I) = ACSS(I)*(1.0D0-RMEW(I)*AA(I)/(GG(I)+RMEW(I)*BB(I))     
     1        * LOG(( GG(I)+RMEW(I)*(BB(I)+AA(I)) )/( RMEW(I)*AA(I) )))   
!                  !   "ω/2 \int_0^1 μ'G'/(μG'+μ'G) d μ'"  ではない？
!
        BETOMG(I) =  
     1            GREEN_mo(NTYP(I),MON,IVEG)       * TRAN1(I)
     1      + ( 1.0D0 -GREEN_mo(NTYP(I),MON,IVEG)) * TRAN2(I)
        BETOMG(I) =   0.50D0 *( OMEGA(I)+( OMEGA(I)- 2.0D0 *BETOMG(I) )
     1              * 0.25D0 *( 1.0D0 -CHIV(I) )**2 )    
!               ! βω ？
!               ! ただし式  とは違う。χL と θbar の間は仮定あり。
!
        BETAO(I)   =   ( 1.0D0  + RMEWBAR(I)*RKK(I) )         ! 分子
     1              / ( OMEGA(I) *RMEWBAR(I)*RKK(I) ) 
     1              * ACSS(I)                               
!               ! β0 1-2-61 ？ ただし "ACCS = ω/2 \int_0^1 d μ'" なら？
 3220 CONTINUE
!
      DO 3230 I=NSTR,NEND
        SS   (I) = ( 1.0D0 - OMEGA(I) + BETOMG(I) ) * RMEWBARI(I)
        TT   (I) = BETOMG(I) * RMEWBARI(I)
        UU   (I) = OMEGA(I)  * RKK(I) * BETAO(I)  
        VV   (I) = OMEGA(I)  * RKK(I) * ( 1.0D0 -BETAO(I) )
        STK  (I) = SS(I)**2 - TT(I)**2 - RKK(I)**2       ! s^2-t^2-K^2
        IF ( ABS(STK(I)).LT.1.D-10 ) THEN
          OMEGA(I) = OMEGA(I) * 0.98D0       ! インチキ臭いエラー処理
          SS   (I) = ( 1.0D0 - OMEGA(I) + BETOMG(I) ) * RMEWBARI(I)
          TT   (I) = BETOMG(I) * RMEWBARI(I)
          UU   (I) = OMEGA(I)  * RKK(I) * BETAO(I)  
          VV   (I) = OMEGA(I)  * RKK(I) * ( 1.0D0 -BETAO(I) )
          STK  (I) = SS(I)**2 - TT(I)**2 - RKK(I)**2       ! s^2-t^2-K^2
        ENDIF 
        STKI (I) = 1 / STK(I)
        TTI  (I) = 1 / TT(I)
        HH   (I) = SQRT ( SS(I)**2 - TT(I)**2 ) 
 3230 CONTINUE
!
      DO 3240 I=NSTR,NEND
        RHL (I) = MIN ( HH  (I)*ZAT(I), 50.0D0 )   ! hL 
        RKL (I) = MIN ( RKK (I)*ZAT(I), 50.0D0 )   ! KL
        EHLI(I) = DEXP ( - RHL(I) )                 ! exp(-hL)
        EKLI(I) = DEXP ( - RKL(I) )                 ! exp(-KL)
        EHL (I) = 1. / EHLI(I)
!       EKL (I) = 1. / EKLI(I)
!
        A23X (I) = ( SS(I) - TT(I) * ROSD(I) + HH(I) ) * EHL (I)
        A23Y (I) = ( SS(I) - TT(I) * ROSD(I) - HH(I) ) * EHLI(I) 
 3240 CONTINUE
!
      ELSE                            ! 雪面 
        DO I=NSTR,NEND
          EXTK(I, IVEG, JVN, 1 ) =10.  ! 0 以外の値をダミーでいれておく
          EXTK(I, IVEG, JVN, 2 ) =100.  ! 同上
        ENDDO 
      ENDIF
!  ここまで、雪の地面計算は、計算しない   ------------------ AAAA 
!
!  -------------
!  >>  散乱光 <<
!  -------------
!
!  ここから、雪の地面計算は、計算しない   ------------------ BBBB
!
      IF ( .NOT. ( IVEG.EQ.2 .AND. INSFLAG.EQ.1 ) ) THEN  ! 雪面計算★以外★
!
!  ---------------------
!  >>  散乱光アルベド <<
!  ---------------------
!
      DO 4000 I=NSTR,NEND
!
        ALP1(I) = 0.
        ALP2(I) =   A23X(I) / ( A23X(I) - A23Y(I) )
        ALP3(I) = 1 - ALP2(I) 
        ALP4(I) = 0.
        ALP5(I) = ( SS(I) - HH(I) ) * TTI(I) * ALP2(I)
        ALP6(I) = ( SS(I) + HH(I) ) * TTI(I) * ALP3(I)
!
        ALBEDO (I,IVEG,JVN,2) = ALP5(I) + ALP6(I) 
        TRA2(I) = ALP2(I) * EHLI(I) + ALP3(I) * EHL(I) 
!
 4000 CONTINUE
!
!  ここまで、雪の地面計算は、計算しない   ------------------ BBBB
!
      ENDIF
!
!  -----------------------------------------------------
!! >>  散乱光アルベド ～ 土の地面は裸地と草の効果考慮 <<
!  -----------------------------------------------------
!
      IF ( IVEG.EQ.2 .AND. INSFLAG.EQ.0 ) THEN                ! 土の地面    
        DO 4100 I=NSTR,NEND
          ALBEDO   (I,2,JVN,2)                              ! 裸地+草 A^μ
     1      =    ROSD(I)            * ( 1.0D0 - CVR_VEG_ms(I,2))
     1         + ALBEDO (I,2,JVN,2) *  CVR_VEG_ms(I,2)
 4100   CONTINUE
      END IF
!
      IF ( IVEG.EQ.2 .AND. INSFLAG.EQ.1 ) THEN                ! 雪面    
        DO 4200 I=NSTR,NEND
          ALBEDO   (I,2,JVN,2) = ALB_SKIN(I,JVN)
 4200   CONTINUE
      END IF
!
!  ここから、雪の地面計算は、計算しない   ------------------ CCCC
!
      IF ( .NOT. ( IVEG.EQ.2 .AND. INSFLAG.EQ.1 ) ) THEN   ! 雪面計算★以外★

!  ---------------------
!! >>  直達光アルベド <<
!  ---------------------
!
      DO 5000 I=NSTR,NEND
!
        ALP1(I) = ( (SS(I) + RKK(I))*VV(I) + UU(I)*TT(I) ) * STKI(I)     
        ALP4(I) = ( (SS(I) - RKK(I))*UU(I) + VV(I)*TT(I) ) * STKI(I)     
!
        ALP2(I) = ( 
     1           TT(I) * ( -ALP4(I)+ROSD(I)*ALP1(I)+ROSB(I) )*EKLI(I)      
     1    + ((SS(I)+HH(I))*ALP1(I) - TT(I)*ROSD(I)*ALP1(I) ) *EHL (I) )  
     1             / ( A23Y(I) - A23X(I) )
!
        ALP3(I) =  - ALP1(I) - ALP2(I)

        ALP5(I) = ( SS(I) - HH(I) ) * TTI(I) * ALP2(I)
        ALP6(I) = ( SS(I) + HH(I) ) * TTI(I) * ALP3(I)
!
        ALBEDO (I,IVEG,JVN,1) = ALP4(I) + ALP5(I) + ALP6(I) 
!
        TRA1  (I) = EKLI(I)           ! exp(-KL) 
        TRA3  (I) =                   ! 散乱光のうち葉を透過される割合
     1                ALP1 (I) * EKLI(I)    ! α1 exp(-KL)
     1              + ALP2 (I) * EHLI(I)    ! α2 exp(-hL)
     1              + ALP3 (I) * EHL (I)    ! α3 exp( hL)
!
 5000 CONTINUE
!
      ENDIF
!    
!  ここまで、雪の地面計算は、計算しない   ------------------ CCCC
!
!  -----------------------------------------------
!! >>  直達光アルベド ～ 地面は裸地・草効果考慮 <<
!  -----------------------------------------------
!
      IF ( IVEG.EQ.2 .AND. INSFLAG.EQ.0 ) THEN      ! 土の地面
        DO 5100 I=NSTR,NEND
          ALBEDO   (I,2,JVN,1)                             ! 裸地 + 草
     1      =  ROSB(I)            * ( 1.0D0 -CVR_VEG_ms(I,2) )
     1       + ALBEDO (I,2,JVN,1) * CVR_VEG_ms(I,2)
 5100   CONTINUE
      END IF
!
      IF ( IVEG.EQ.2 .AND. INSFLAG.EQ.1 ) THEN       ! 雪面    
         DO 5200 I=NSTR,NEND
           ALBEDO   (I,2,JVN,1) = ALB_SKIN(I,JVN)
 5200   CONTINUE
      END IF
!
 2100 CONTINUE
!
!  =====================================================
!! >> キャノピー(IVEG=2)・裸地(IVEG=1) 別処理ここまで <<
!  =====================================================
!
!  ======================================================
!! >> 以下、出力量 ( ABSORB_SHRT, ALB_SIB , REACH_GRD ) 計算 <<
!  ======================================================
!
!  -----------------------------------------------------------
!! >> 地上での下向き短波中の、吸収 absorb する割合 ABSORB_SHRT <<
!  -----------------------------------------------------------
      DO 7000 I=NSTR,NEND
!
        ABSORB_SHRT(I,2,JVN,1)=
     1        ( 1.0D0 -CVR_VEG_ms(I,1)) * ( 1.0D0 - ALBEDO(I,2,JVN,1))
     2        +CVR_VEG_ms(I,1) *( TRA1(I)*( 1.0D0 - ALBEDO(I,2,JVN,1))
     3                       +TRA3(I)*( 1.0D0 - ALBEDO(I,2,JVN,2)))
!
!         直達光のうち地面で吸収する割合
!         =   地面に直接到達したうち、地面で反射されない分
!           + キャノピーを直達光のまま透過して、地面で反射されない分
!           + キャノピーで散乱されながらも透過し、地面で反射されない分
!
        ABSORB_SHRT (I,2,JVN,2)=
     1        ( 1.0D0 -CVR_VEG_ms(I,1))   * (1.0D0 -ALBEDO(I,2,JVN,2))
     2        +CVR_VEG_ms(I,1) * TRA2(I)* (1.0D0 -ALBEDO(I,2,JVN,2))
!
!         散乱光のうち地面で吸収する割合
!         =   地面に直接到達したうち、反射されない分
!           + キャノピーに到達した後透過し、地面で反射されない分
!
!         キャノピーに到達した散乱光がさらに散乱された分は？
!
        ABSORB_SHRT(I,1,JVN,1) = CVR_VEG_ms(I,1) * (  
     1                     ( 1.0D0 -ALBEDO(I,1,JVN,1) )
     &         - TRA1(I) * ( 1.0D0 -ALBEDO(I,2,JVN,1) )
     &         - TRA3(I) * ( 1.0D0 -ALBEDO(I,2,JVN,2) ) )
!

!
!!         直達光のうちキャノピーで吸収される割合
!!         =   キャノピーに到達したうち、
!!           ( キャノピーで反射されない分
!!             - キャノピーを直達光のまま透過され地面で吸収される分
!!             - キャノピーで散乱光になり透過され地面で吸収される分 )
!
!!         注 : 地面で反射されれば、キャノピーで吸収される。
!
        ABSORB_SHRT(I,1,JVN,2) = CVR_VEG_ms(I,1) * (
     1                     ( 1.0D0 -ALBEDO(I,1,JVN,2) )
     &       - TRA2(I) * ( 1.0D0 -ALBEDO(I,2,JVN,2) ) )
!
!!         散乱光のうちキャノピーで吸収される割合
!!         =   キャノピーに到達したうち、
!!           ( キャノピーで反射されない分
!!             - キャノピーで透過され地面で吸収される分 )
!
!!         注 : 地面で反射されれば、キャノピーで吸収される。
!
 7000 CONTINUE
!
!  --------------------------
!! >> 全アルベド ALB_SIB 計算 <<
!  --------------------------
!
      DO 7100 IRD=1,2
      DO 7100 I=NSTR,NEND
        ALB_SIB(I,JVN,IRD) = 
     1         ( 1.0D0 - CVR_VEG_ms(I,1)) * ALBEDO(I,2,JVN,IRD)  ! 地面
     1      +    CVR_VEG_ms(I,1)          * ALBEDO(I,1,JVN,IRD)  ! キャノピー
 7100 CONTINUE
!
!  --------------------------------------------------------------
!! >> 可視で大気下端を抜けたうちで地面上端に届く割合 REACH_GRD <<
!  --------------------------------------------------------------
!!    キャノピー上端で反射される分は tranc? で考慮。
!     SAVING OF EXTINCTION COEFFICIENTS ( PAR ) FOR STOMAT CALCULATION
!
      IF ( JVN.EQ.1 ) THEN                               !  可視のみ
        DO 7200 I=NSTR,NEND
          REACH_GRD(I,1) =                                  ! 直達
     1             1.0D0 -CVR_VEG_ms(I,1)                          ! 地面
     1             + CVR_VEG_ms(I,1) * ( TRA1(I)+TRA3(I) )     ! キャノピー
          REACH_GRD(I,2) =                                         ! 散乱
     1             1.0D0 -CVR_VEG_ms(I,1)                          ! 地面
     1             + CVR_VEG_ms(I,1) * TRA2(I)                   ! キャノピー
 7200   CONTINUE
!   
      END IF
!
 2000 CONTINUE      !  <-----------------------  波長別処理ここまで
!
      RETURN
      END

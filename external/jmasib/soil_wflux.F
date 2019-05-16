! 土壌中の水フラックスの評価
! 99/07/06 保坂征宏 
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/09/22 保坂征宏 XFR_WTR の下限処理。
!                   凍った時 XFR_WTR(I,1)**(-BEE_co(J))でまずいので。
! 99/09/23 保坂征宏 dep_h_inv(i,3) の式を変える
! 99/09/24 保坂征宏 方程式の対角成分に問題あり。FACT_WGT_h2o 導入
! 99/09/24 保坂征宏 拡散係数の与え方を変更。（JCN_COEF_WFLUX_SOIL=0 の時。）
!
      SUBROUTINE SOIL_WFLUX (
     I   NSTR             , NEND             , 
     I   INSFLAG          , 
     I   FR_WTR_SOIL      , 
     I   FR_GLA_SOIL      , 
     I   FACT_WGT_H2O      , 
     I   DELT_SOIL        , 
     I   NTYP             , 
     I   TMP_SOIL         , 
     O   WFLUX_H          , WTR_INC ) 
! /*
! イメージ
!   TMP_FREZ よりも低温ならば 土 + 空気 + 氷
!                             ただし水が浸透してくることはある
!   TMP_FREZ ならば           土 + 空気 + 水
!   TMP_FREZ よりも高温ならば 土 + 空気 + 氷 + 水
!
! 問題点。
!   ・ポテンシャルの下限処理、透水係数の下限処理、...  XK, POT_SOIL
!     また下限時には重力ポテンシャル差でのフラックスをなくしたい。
!     1 を変数化して、0 or 1 にする?
!
! 要チェック
!  ・上端から入って来る水は、降水 - 裸地蒸発であるので、必ずしも正ではない。
!    これに対する対策はたてられているのか？
!    とくに、所有する水よりも多くの液水の蒸発を要求された場合の措置は？
!
! 注意。拡散係数の与え方について。99/09/24
!
! メモ 1 : 方程式と解 ～ 陽的に解く場合
!
!   空隙中を水が占める割合 FR_WTR = w とする。
!   水の質量は FACT_WGT_H2O * w である。ただし FACT_WGT_H2O=RHO_WTR_C*POR 
!
!   方程式は、
!      水の質量の時間変化率 = 質量フラックス収束
!   の形である。つまり
!      d (FACT_WGT_H2O*w)/dt = 上からの水フラックス - 下への水フラックス
!   である。
!   陽に評価される水フラックス WFLUX_H_exp は、
!      WFLUX_H_exp = XK_H * ( POT_SOIL - POT_SOIL ) / DEP_H * RHO_WTR_c -- *E
!   という形である。
!
!   解くべき方程式は
!      d (FACT_WGT_H2O*w)/dt =   WFLUX_H_exp_上 - WFLUX_H_exp_下
!   である。このサブルーチンでは *E で表される各フラックスを計算して、
!   上位のサブルーチンに返す。
!
! メモ 2 : 方程式と解 ～ 陰的に解く場合
!
!   陰的に解く場合は、フラックスとして補正項 WFLUX_H_cor が加わる。
!
!   水フラックスの補正量は、POT_SOIL と w の関係
!      POT_SOIL = POT_SAT_SOIL * w'**(-B)  
!      w'       = w/(1-空隙中の氷の割合)
!   すなわち
!      POT_SOIL = POT_SAT_SOIL  * (1-氷の割合)**B * w**(-B) 
!   を反映して、以下のようになる。ただし Δw = dw/dt * DELT 
!   WFLUX_H_cor_上
!   = XK_H / DEP_H * POT_SAT_SOIL * (1-氷の割合)**B * (-B) * w**(-B-1)  
!          * (Δw_上の層 - Δw_当該層)  * RHO_WTR_C 
!   = XK_H / DEP_H * POT_SAT_SOIL                   * (-B) * w'**(-B)/w 
!          * (Δw_上の層 - Δw_当該層)  * RHO_WTR_C 
!   = XK_H / DEP_H * POT_SOIL                       * (-B)           /w * d
!          * (Δw_上の層 - Δw_当該層)  * RHO_WTR_C 
!   = XK_H / DEP_H * POT_SOIL * (-B)/w * (Δw_上の層 - Δw_当該層) * RHO_WTR_C 
!   になる。
!
!   方程式は
!      d (FACT_WGT_H2O*w)/dt =   WFLUX_H_exp_上 - WFLUX_H_exp_下
!                              + WFLUX_H_cor_上 - WFLUX_H_cor_下 
!   すなわち
!      FACT_WGT_H2O * dw/dt  
!          = (   WFLUX_H_exp_上 - WFLUX_H_exp_下
!              + WFLUX_H_cor_上 - WFLUX_H_cor_下 ) 
!   である。
!
!   ここでは簡単のために当該層の補正分のみを記すことにしよう。
!       WFLUX_H_cor_上 
!         = - XK_H / DEP_H * POT_SOIL * (-B) /w * dw/dt * DELT * RHO_WTR_C  
!   ここで
!       XC_H  = XK_H / DEP_H 
!       XC    = POT_SOIL / (-B) / w * DELT * RHO_WTR_C 
!       DF_xx = XC_H * XC        
!   と書くと、
!       WFLUX_H_cor_上 = - XC_H * XC * dw/dt 
!                      = - DF_xx     * dw/dt 
!   である。
!
!   方程式は
!      ( FACT_WGT_H2O + DF_xx ) dw/dt = 陽的フラックス収束
!   であり、陰的補正項はこの解 AA を用いて
!       WFLUX_H_cor_上 = - DF_xx     * AA 
!   のように書くことができる。
!
! メモ 2 ～ 第 3 層とその下のやりとり
!   3 層目とその下とのフラックスについての基本的な考え方は、
!     ・下には抜ける。ランオフ
!       3 層目の拡散係数（外挿したということ）、
!       2 層目・3 層目間のポテンシャル差（外挿したということ）、
!       地面の傾き効果の積。
!     ・上に来るかどうかは考え方次第。デフォルトでは無し。
!       ある場合はランオフと同じ評価方法。
!
! */       
!
!
      USE SIBPRM , ONLY : 
     1   IDX   ,
     1   ITYP  , 
     1   IDP  
!
      USE SIBCON , ONLY :
     1   DELT_FAC_SOIL_c  , 
!    1   JCN_COEF_WFLUX_SOIL_c  ,
     1   TMP_FREZ_c        ,
     1   RHO_WTR_C        
     1   , JCN_RUNOFF_GRAV_c            ! default=1 2004/08/10 Hosaka 
     1   , COEF_SOIL_WFLUX_MIN_c        ! default=0 
     1   , POT_SOIL_MIN_c              
!
      USE SIBVEG , ONLY : 
     1           SINX_co                 , ! 土壌層の傾き  (ITYP)
     1           POR_SOIL_co            ,  ! 空隙率   (ITYP)  
     1           BEE_co                 ,  ! 土壌水分依存性係数 (ITYP)
     1           DEP_SOIL_co             , ! 土壌層の深さ    (ITYP,IDP)
     1           POT_SAT_SOIL_co        ,  ! 飽和時の土壌ポテンシャル (ITYP)
     1           WCOEF_SAT_SOIL_co         ! 飽和透水係数    (ITYP)
!
      USE COM_RUNCONF_SIB0109 , ONLY :
     1   JCN_RRUNOFF_B        ,
     1   JCN_SOIL_WATER_EXPLICIT       ! implicit -> 0 , explicit -> 1 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS,
     1                      IMONIT_LEVEL, IMONIT_ALL , 
     1                       IMONIT_FLUX_COMMON
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG         ! 雪無し 0, 雪下 1 
!
!  定数
!
      INTEGER,INTENT(IN)    :: NTYP      (IDX)
!
      REAL(8),INTENT(IN)    :: DELT_SOIL       ! 時間刻
!
      REAL(8),INTENT(IN)    :: FR_WTR_SOIL (IDX,IDP)  ! 
      REAL(8),INTENT(IN)    :: FR_GLA_SOIL (IDX,IDP)
!
      REAL(8),INTENT(IN)    :: TMP_SOIL    (IDX,IDP)   ! 
!
!
      REAL(8),INTENT(IN)    :: FACT_WGT_H2O  (IDX,IDP)    ! 空隙量x水の比重
!
! OUTPUT
!
      REAL(8),INTENT(OUT)   :: WFLUX_H(IDX,IDP) ! 水フラックス
      REAL(8),INTENT(OUT)   :: WTR_INC(IDX,IDP) ! 水増分
!
! SAVE 
!
      REAL(8),SAVE    :: DEP_H_INV (ITYP,IDP)     ! 1/層平均厚さ
!
! WORK
!
      REAL(8)   :: XPOT_SOIL(IDX,IDP)       ! 土壌ポテンシャル
      REAL(8)   :: XFR_WTR  (IDX,IDP) ! W 空隙中で液水の占める体積比
!     REAL(8)   :: XPOROS   (IDX,IDP) ! 空隙率 （1- 土 - 氷）
      REAL(8)   :: XK    (IDX,IDP)    ! 透水係数 （フル）
      REAL(8)   :: XK_H  (IDX,IDP)    ! 透水係数 （ハーフ）
      REAL(8)   :: EPSILON(IDX)       ! q_{3+1/2} 有無フラグ     
      REAL(8)   :: WFLUX_H_EXP (IDX,IDP)
      REAL(8)   :: WFLUX_H_COR (IDX,IDP)
!
      REAL(8)   :: AA    (IDX,3,4)    ! 方程式の係数
      REAL(8)   :: A11I
      REAL(8)   :: A22I
      REAL(8)   :: A33I
      REAL(8)   :: XC_H    (IDX,IDP)  ! 方程式の係数用作業変数
      REAL(8)   :: XC      (IDX,IDP)  ! 方程式の係数用作業変数
      REAL(8)   :: DF12_DW1(IDX)    ! 方程式の係数用作業変数
      REAL(8)   :: DF12_DW2(IDX)    ! 方程式の係数用作業変数
      REAL(8)   :: DF23_DW2(IDX)    ! 方程式の係数用作業変数
      REAL(8)   :: DF23_DW3(IDX)    ! 方程式の係数用作業変数
      REAL(8)   :: DF34_DW2(IDX)    ! 方程式の係数用作業変数
      REAL(8)   :: DF34_DW3(IDX)    ! 方程式の係数用作業変数
!
      REAL(8)   :: DELS    (IDX)     
      REAL(8)   :: WMAX 
      REAL(8)   :: WMIN
      REAL(8)   :: PSMAX 
      REAL(8)   :: PSMIN
      REAL(8)   :: BOT 
!
      REAL(8)   :: PHASE_SOIL(IDX,IDP) ! 零度以上 1 、以下 0 
!                                      ! 零度以下ではポテンシャルは変わらない
!
      REAL(8)   :: WP3_1 
      REAL(8)   :: WP3_2 
      REAL(8)   :: WP3_3 
      REAL(8)   :: WM1_1  
      REAL(8)   :: WM1_2  
      REAL(8)   :: WM1_3  
!
      REAL(8)   :: AVK          (IDX,IDP)  
      REAL(8)   :: DIV   
      REAL(8)   :: AVKMIN 
      REAL(8)   :: AVKMAX 
      REAL(8)   :: SAT 
!
      REAL(8)   :: DELT_SL
!
      REAL(8)   :: XFR_GLA
      REAL(8)   :: XFACT        (IDX,IDP)
      REAL(8)   :: XFACT_INV    (IDX,IDP)
      REAL(8)   :: XFACT_WGT_H2O(IDX,IDP)
!     
      INTEGER   :: I
      INTEGER   :: J
      INTEGER   :: K
!
      REAL(8),PARAMETER :: DIV_MIN = 1.D-7    ! 昔と同じ値
      REAL(8),PARAMETER :: ZERO   = 0.D0
      REAL(8),PARAMETER :: ONE = 1.D0 
      REAL(8),PARAMETER :: XFR_WTR_MIN = 1.D-2 ! 各種係数評価に使う飽和度下限
      REAL(8),PARAMETER :: XFR_GLA_MAX = 0.9999D0
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
#endif
!
      LOGICAL,SAVE :: LFIRST
!
      DATA LFIRST / .TRUE. /
!
!
! ==================================
! >> 一度計算すればいい定数の設定 << 
! ==================================
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SOIL_WFLUX 99/08 POCHAKA'
        LFIRST = .FALSE. 
!
!       注意。dep_h_inv(I,3) の与え方について。  99/09/24 
!             3,4 層目間の水フラックスを評価する際に、
!             4 層目は 3 層目と同じ厚さだと考え、
!             3,4, 層目の各々の中央での値の差から評価する、という立場。
!
        DO 100 J = 1,ITYP
          DEP_H_INV(J,1) = 2. / ( DEP_SOIL_co(J,1)+ DEP_SOIL_co(J,2) ) 
          DEP_H_INV(J,2) = 2. / ( DEP_SOIL_co(J,2)+ DEP_SOIL_co(J,3) ) 
          DEP_H_INV(J,3) = 1. / DEP_SOIL_co(J,3)   
 100    CONTINUE   
!
      ENDIF
!
      DELT_SL      = DELT_SOIL * DELT_FAC_SOIL_c    ! 意味は?  000731
!
!     write(6,*) 'soil_wflux KKK delt_fac ' , 
!    1           DELT_FAC_SOIL_1 , DELT_FAC_SOIL_2 , delt_soil
!
! =======================================
! >> ポテンシャル変化フラグ PHASE_SOIL << 
! =======================================
!  零度以下では水分ポテンシャルは変化しない
!
      PHASE_SOIL(:,:) = ONE 
! 
      DO K=1,IDP
      DO I=NSTR,NEND
        IF ( TMP_SOIL(I,K) .LT. TMP_FREZ_c ) THEN
          PHASE_SOIL(I,K) = 0.
        ENDIF
      ENDDO
      ENDDO
!
! ======================
! >> 透水係数等の計算 << 
! ======================
!
      DO 1000 K=1,IDP
      DO 1000 I=NSTR,NEND
        J = NTYP(I)
!
! 氷の存在による飽和度等にかける係数の設定
!       XFACT     = (1-XFR_GLA1) 
!       XFACT_INV = 1/(1-XFR_GLA1) 
!       XFACT     * XFR_W = FR_W 
!
        XFR_GLA  = MIN ( FR_GLA_SOIL (I,K) , XFR_GLA_MAX )  ! max で 0.9999
        XFACT    (I,K) = 1 - XFR_GLA                        ! min で 0.0001
        XFACT_INV(I,K) = 1./XFACT(I,K)                      ! max で 10000 
!
! 飽和度 XFR_WTR  （ 空隙中で水の占める体積比 ）
!
        XFR_WTR  (I,K) = FR_WTR_SOIL (I,K) * XFACT_INV(I,K)
!
#ifdef SIB_DEBUG 
        if ( xfr_wtr(I,K)-1 .GT. 1.D-7 ) THEN
          write(6,*) 'soil_wflux error xfr_wtr' , i , K ,xfr_wtr(I,K),
     1               XFR_GLA , XFACT(I,K) ,  FR_WTR_SOIL (I,K) , 
     1               XFACT_INV(I,K) , 
     1               FR_WTR_SOIL (I,K) + FR_GLA_SOIL (I,K) 
          stop 999
        endif
#endif
!
        XFR_WTR  (I,K) = MAX ( XFR_WTR  (I,K) , XFR_WTR_MIN )  ! min で 0.01 
!
! 空隙率 XPOROS （ 空隙の割合。1 - 土 - 氷 )
!       XPOROS (I,K) = POR_SOIL_co(J) * (1.-FR_GLA_SOIL(I,K)) 
!
! 土壌水ポテンシャル XPOT_SOIL
!   問題点... 下限処理。水が無い場合は有限値に。
!             また水が少ない場合、重力ポテンシャルが効かないようにしたい。
!
        XPOT_SOIL(I,K) = POT_SAT_SOIL_co(J)*XFR_WTR(I,K)**(-BEE_co(J))    
        XPOT_SOIL(I,K) = MAX ( XPOT_SOIL(I,K) , POT_SOIL_MIN_c ) 
!
        XFACT_WGT_H2O(I,K) = FACT_WGT_H2O(I,K) * XFACT(I,K) 
!
! 透水係数 XK 
!    問題点... 下限処理。水が無い場合は 0 にする。
!
        XK(I,K) = WCOEF_SAT_SOIL_co(J) 
     1                * XFR_WTR(I,K)**(2.*BEE_co(J)+3.)    
 1000 CONTINUE  
!
! ========================
! >> 方程式の係数の設定 <<
! ========================
!
!     IF ( JCN_COEF_WFLUX_SOIL_c .EQ. 0 ) THEN   ! MILLY AND EAGLESON (1982)?
      DO 4000 I=NSTR,NEND 
        J = NTYP(I)
        WP3_1     = XFR_WTR (I,1)**(2*BEE_co(J)+3.) 
        WP3_2     = XFR_WTR (I,2)**(2*BEE_co(J)+3.) 
        WP3_3     = XFR_WTR (I,3)**(2*BEE_co(J)+3.) 
!
        WM1_1     = XFR_WTR (I,1)**(-BEE_co(J)) 
        WM1_2     = XFR_WTR (I,2)**(-BEE_co(J)) 
        WM1_3     = XFR_WTR (I,3)**(-BEE_co(J)) 
!
!   上限の評価 （ 陰的補正項係数 DELS = dψ/dw > 0 ) 
!
        WMAX      = MAX ( XFR_WTR (I,1), XFR_WTR (I,2), XFR_WTR (I,3) )
        WMIN      = MIN ( XFR_WTR (I,1), XFR_WTR (I,2), XFR_WTR (I,3) )
        PSMAX     = MAX ( WM1_1        , WM1_2        , WM1_3         )
        PSMIN     = MIN ( WM1_1        , WM1_2        , WM1_3         )
        BOT       = MAX ( 1.D-7 , WMAX - WMIN ) 
        DELS(I)   = - ( PSMAX - PSMIN ) / BOT * POT_SAT_SOIL_co(J)   
!                   ^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^
!                      正                          負
!
!   第一層と第二層
!
        AVK(I,1)  = WP3_1 * WM1_1 - WP3_2 * WM1_2
        DIV       = WM1_1 - WM1_2
!
        AVKMIN    = MIN ( WP3_1, WP3_2 )                                
        AVKMAX    = MAX ( WP3_1, WP3_2 )  
! 
        IF ( ABS(DIV).GT.DIV_MIN ) THEN   ! 普通通る
          AVK(I,1) = - AVK(I,1) * BEE_co(J) / ( (BEE_co(J)+3) * DIV )  
        ENDIF
!
        AVK(I,1)  = MIN ( AVK(I,1) , AVKMAX )
        AVK(I,1)  = MAX ( AVK(I,1) , AVKMIN )
        SAT       = WCOEF_SAT_SOIL_co(J) 
     1               * ( 1 - MAX(FR_GLA_SOIL(I,1),FR_GLA_SOIL(I,2)) )
        XK_H(I,1) = MAX ( SAT * AVK(I,1) , COEF_SOIL_WFLUX_MIN_c ) 
!
!   第二層と第三層
!
        AVK(I,2)  = WP3_2 * WM1_2 - WP3_3 * WM1_3 
        DIV       = WM1_2 - WM1_3
        AVKMIN    = MIN( WP3_2, WP3_3 )                                
        AVKMAX    = MAX( WP3_2, WP3_3 )  
!
        IF ( ABS(DIV).GT.DIV_MIN ) THEN
          AVK(I,2)= - AVK(I,2) * BEE_co(J) / ( (BEE_co(J)+3) * DIV )  
        ENDIF
!
        AVK(I,2)  = MIN ( AVK(I,2) , AVKMAX )
        AVK(I,2)  = MAX ( AVK(I,2) , AVKMIN )
        SAT       = WCOEF_SAT_SOIL_co(J) 
     1             * ( 1 - MAX(FR_GLA_SOIL(I,2),FR_GLA_SOIL(I,3)) )
        XK_H(I,2) = MAX ( SAT * AVK(I,2) , COEF_SOIL_WFLUX_MIN_c ) 
!
!   第三層とその下
!
        AVK (I,3) = XK  (I,3)       ! ただの K
        XK_H(I,3) = XK  (I,3)       ! *SAT しなくていいのかな?? いいはず?
 4000 CONTINUE
!
!     ELSE
!       DO 4100 I=NSTR,NEND
!         J = NTYP(I) 
!         XK_H(I,1) 
!    1       = ( DEP_SOIL_co(J,1) + DEP_SOIL_co(J,2) ) 
!    1         * XK(I,1) * XK(I,2)
!    1        / ( DEP_SOIL_co(J,1)* XK(I,2)+DEP_SOIL_co(J,2)* XK(I,1) )
!         XK_H(I,2) 
!    1       = ( DEP_SOIL_co(J,2) + DEP_SOIL_co(J,3) ) 
!    1         * XK(I,2) * XK(I,3)
!    1        / ( DEP_SOIL_co(J,2)* XK(I,3)+DEP_SOIL_co(J,3)* XK(I,2) )     
!         XK_H(I,3) = XK  (I,3) 
!4100   CONTINUE  
!     ENDIF
!
!     層間の量   XC_H = XK_H / DEP_H 
!     層の量     XC   = XPOT_SOIL * (- BEE_co(J) ) / w * m dt * RHO_WTR_C 
!
      DO 4200 I=NSTR,NEND
        J = NTYP(I) 
!
!  XC_H の計算
!     密度 x ( AVK x K_s ) / 厚さ          > 0 
!
        XC_H  (I,1) = RHO_WTR_C * XK_H(I,1) * DEP_H_INV(J,1) 
        XC_H  (I,2) = RHO_WTR_C * XK_H(I,2) * DEP_H_INV(J,2)
        XC_H  (I,3) = RHO_WTR_C * XK_H(I,3) * DEP_H_INV(J,3) 
     1                                      * SINX_co(J) 
!
        WFLUX_H_EXP(I,1) 
     1     = - XC_H(I,1) 
     1        * ( XPOT_SOIL(I,2) - XPOT_SOIL(I,1) ) 
!            + RHO_WTR_C * AVK(I,1)                 ! 重力ポテンシャル
!
        WFLUX_H_EXP(I,2) 
     1     = - XC_H(I,2) 
     1        * ( XPOT_SOIL(I,3) - XPOT_SOIL(I,2) )  
!            + RHO_WTR_C * AVK(I,2)                 ! 重力ポテンシャル
!
        WFLUX_H_EXP(I,3) 
     1     = - XC_H(I,3) 
     1        * ( XPOT_SOIL(I,3) - XPOT_SOIL(I,2) )  
!    1        + RHO_WTR_C * AVK(I,3) * sinx_co(j) 
!
!  陰的補正項関係
!     dψ/dw 
!     とりあえずは各層ともおなじ値 DELS にしてみる
!
        XC    (I,1) = DELS (I) 
        XC    (I,2) = DELS (I) 
        XC    (I,3) = DELS (I) 
!
!
 4200 CONTINUE  
!
      IF ( JCN_RUNOFF_GRAV_c .GE. 1 ) THEN 
        DO I=NSTR,NEND
          WFLUX_H_EXP(I,3) = WFLUX_H_EXP(I,3) 
     1        + RHO_WTR_C * AVK(I,3) * sinx_co(j) 
                                                     ! 重力ポテンシャル
        ENDDO
      ENDIF
!
      IF ( JCN_RUNOFF_GRAV_c .GE. 2 ) THEN 
        DO I=NSTR,NEND
          WFLUX_H_EXP(I,1) = WFLUX_H_EXP(I,1) 
!    1         + RHO_WTR_C * AVK(I,1)               
     1         + RHO_WTR_C * XK_H(I,1)               
          WFLUX_H_EXP(I,2) = WFLUX_H_EXP(I,2) 
!    1         + RHO_WTR_C * AVK(I,2)               
     1         + RHO_WTR_C * XK_H(I,2)               
        ENDDO
      ENDIF
!
! 第 3 層目への下からの流入
!     JCN_RRUNOFF_B = 0 （上への逆流禁止）で、WFLUX_H が負（上向き）の場合
!
      DO 4300 I=NSTR,NEND 
        EPSILON(I) = 1.
 4300 CONTINUE  
!
      IF ( JCN_RRUNOFF_B .EQ. 0 ) THEN
        DO 4310 I=NSTR,NEND 
          IF ( WFLUX_H_EXP(I,3) .LT. 0. ) THEN
            EPSILON(I) = 0. 
            WFLUX_H_EXP(I,3) = 0. 
          ENDIF
 4310   CONTINUE  
      ENDIF
!
!  ================================================
!  >> フラックスの仮決定 ～ インプリシットの場合 <<
!  ================================================
!
      IF ( JCN_SOIL_WATER_EXPLICIT.EQ.0 ) THEN
!
! ------------------------
! > インプリシット寄与分 < 
! ------------------------
!
!   例
!      DF12_DW1 ... W1 が増えたときに 1->2 へのフラックスの増分
!        XC_H(I,1)       = 1 層と 2 層との間の定数
!        XC  (I,1)       = W1 が増えたときの 1 層目のポテンシャル増分
!        PHASE_SOIL(I,1) = W1 が増えたとしてポテンシャルが増えるか?
!      
      DO I=NSTR,NEND
        DF12_DW1(I) =   XC_H (I,1) * XC (I,1) 
!    1                             * PHASE_SOIL(I,1) 
        DF12_DW2(I) = - XC_H (I,1) * XC (I,2) 
!    1                             * PHASE_SOIL(I,2)         
        DF23_DW2(I) =   XC_H (I,2) * XC (I,2) 
!    1                             * PHASE_SOIL(I,2)         
        DF23_DW3(I) = - XC_H (I,2) * XC (I,3) 
!    1                             * PHASE_SOIL(I,3)         
        DF34_DW2(I) =   XC_H (I,3) * XC (I,2) 
!    1                                        * PHASE_SOIL(I,2) 
!    1                                        * EPSILON(I)
        DF34_DW3(I) = - XC_H (I,3) * XC (I,3) 
!    1                                        * PHASE_SOIL(I,3)         
!    1                                        * EPSILON(I) 
      ENDDO
!
!  ------------------
!  > 方程式を建てる <
!  ------------------
!
      DO 5000 I=NSTR,NEND
        AA(I,1,1) 
     1     =   XFACT_WGT_H2O (I,1)  / DELT_SL 
     1       + DF12_DW1 (I)
        AA(I,1,2) 
     1     =   DF12_DW2 (I)
        AA(I,1,3) =   0. 
        AA(I,1,4) 
     1     =                - WFLUX_H_EXP(I,1) 
!
        AA(I,2,1) 
     1     = - DF12_DW1 (I)
        AA(I,2,2) 
     1     =   XFACT_WGT_H2O (I,2) / DELT_SL  
     1       - DF12_DW2 (I) + DF23_DW2(I)  
        AA(I,2,3) 
     1     =                + DF23_DW3(I)  
        AA(I,2,4) 
     1     =   WFLUX_H_EXP(I,1) - WFLUX_H_EXP(I,2) 
!
        AA(I,3,1) =   0. 
        AA(I,3,2) 
     1     = - DF23_DW2(I) 
!    1       + DF34_DW2(I)               ! 下端排水陰的補正
        AA(I,3,3) 
     1     =   XFACT_WGT_H2O (I,3) / DELT_SL 
     1       - DF23_DW3(I) 
!    1       + DF34_DW3(I)               ! 下端排水陰的補正
        AA(I,3,4) 
     1     =   WFLUX_H_EXP(I,2) - WFLUX_H_EXP(I,3) 
 5000 CONTINUE  
!
! ----------------
! > 方程式を解く <
! ----------------
!
      DO 6000 I=NSTR,NEND
        A33I = 1./AA(I,3,3) 
        AA(I,3,2) = AA(I,3,2) * A33I
        AA(I,3,4) = AA(I,3,4) * A33I
!       AA(I,3,3) = 1. 
!
        AA(I,2,2) = AA(I,2,2) - AA(I,3,2) * AA(I,2,3)     
        AA(I,2,4) = AA(I,2,4) - AA(I,3,4) * AA(I,2,3)     
!       AA(I,2,3) = 0. 
!
        A22I = 1./AA(I,2,2)
        AA(I,2,1) = AA(I,2,1) * A22I 
        AA(I,2,4) = AA(I,2,4) * A22I 
!       AA(I,2,2) = 1.      
!
        AA(I,1,1) = AA(I,1,1) - AA(I,2,1) * AA(I,1,2)     
        AA(I,1,4) = AA(I,1,4) - AA(I,2,4) * AA(I,1,2)     
!       AA(I,1,2) = 0. 
!
        A11I = 1./AA(I,1,1)
        AA(I,1,4) = AA(I,1,4) * A11I                      ! dW1/dt * Δt
!       AA(I,1,1) = 1. 
!
        AA(I,2,4) = AA(I,2,4) - AA(I,1,4) * AA(I,2,1)     ! dW2/dt * Δt
!       AA(I,2,1) = 0. 
!
        AA(I,3,4) = AA(I,3,4) - AA(I,2,4) * AA(I,3,2)     ! dW3/dt * Δt
!       AA(I,3,2) = 0. 
 6000 CONTINUE
!
! ---------------------- 
! > 水フラックスの決定 <
! ---------------------- 
!
      DO 7000 I=NSTR,NEND
!
! フラックス WFLUX_H の計算 ... 下向き正
!
        WFLUX_H_COR(I,1) 
     1     = 
     1       + (   DF12_DW1(I) * AA(I,1,4) 
     1           + DF12_DW2(I) * AA(I,2,4) ) 
!!   1         * DELT_SL
!
        WFLUX_H_COR(I,2) 
     1     = 
     1       + (   DF23_DW2(I) * AA(I,2,4) 
     1           + DF23_DW3(I) * AA(I,3,4) ) 
!!   1         * DELT_SL
!
        WFLUX_H_COR(I,3) 
     1     = 0
!    1       + (   DF34_DW2(I) * AA(I,2,4) 
!    1           + DF34_DW3(I) * AA(I,3,4) ) 
!              * DELT_SL
!
        WFLUX_H(I,1) = WFLUX_H_EXP(I,1) + WFLUX_H_COR(I,1) 
!
        WFLUX_H(I,2) = WFLUX_H_EXP(I,2) + WFLUX_H_COR(I,2) 
!
        WFLUX_H(I,3) = WFLUX_H_EXP(I,3) + WFLUX_H_COR(I,3) 
!
 7000 CONTINUE  
!
!  ==================================================
!  >> フラックスの仮決定 ～ エクスプリシットの場合 <<
!  ==================================================
!
      ELSE 
        DO I=NSTR,NEND
          WFLUX_H(I,1) = WFLUX_H_EXP(I,1) 
          WFLUX_H(I,2) = WFLUX_H_EXP(I,2) 
          WFLUX_H(I,3) = WFLUX_H_EXP(I,3) 
        ENDDO
      ENDIF
!
!  ==========================
!  >> 第三層の下からの流入 <<
!  ==========================
!     評価結果が WFLUX_H(I,3) < 0 、つまり下からの流入を要請する場合に、
!     それにこたえるかどうか。
!
      IF ( JCN_RRUNOFF_B .EQ. 0 ) THEN
        DO 8000 I=NSTR,NEND
          WFLUX_H(I,3) = MAX ( WFLUX_H(I,3) , ZERO )
 8000   CONTINUE 
      ENDIF 
!
!  増加する水量計算
!
      DO I=NSTR,NEND 
        WTR_INC(I,1) = (              - WFLUX_H(I,1) ) 
     1                 * DELT_SL 
        WTR_INC(I,2) = ( WFLUX_H(I,1) - WFLUX_H(I,2) ) 
     1                 * DELT_SL
        WTR_INC(I,3) = ( WFLUX_H(I,2) - WFLUX_H(I,3) )
     1                 * DELT_SL
      ENDDO
!
!  ============================
!  >> フラックス量のチェック <<
!  ============================
!    水フラックス計算で、実際にある以上のものを流している場合は補正
!
!
! ================
! >> 詳細モニタ <<
! ================
!
! 地下ランオフ、および、01 以外の水フラックスに伴う熱は soil_flux_conv で。
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!
! 水フラックス 1 -> 2 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN12'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,1), 
     1                     ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS12'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,1), 
     1                     ICS, DELT_SOIL,ONE)  
        ENDIF
!
! 水フラックス 2 -> 3 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN23'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,2), 
     1                     ICN   , DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS23'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,2), 
     1                     ICS   , DELT_SOIL,ONE)  
        ENDIF
!
! 水フラックス 2 -> 3 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN34'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,3), 
     1                     ICN   , DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS34'
          CALL MONIT_SIB ( CMARK , WFLUX_H_EXP(1,3), 
     1                     ICS   , DELT_SOIL,ONE)  
        ENDIF
!
! 水フラックス 1 -> 2 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN12C'
          CALL MONIT_SIB ( CMARK , WFLUX_H_COR(1,1), 
     1                     ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS12C'
          CALL MONIT_SIB ( CMARK , WFLUX_H_COR(1,1), 
     1                     ICS, DELT_SOIL,ONE)  
        ENDIF
!
! 水フラックス 2 -> 3 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN23C'
          CALL MONIT_SIB ( CMARK , WFLUX_H_COR(1,2), 
     1                     ICN   , DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS23C'
          CALL MONIT_SIB ( CMARK , WFLUX_H_COR(1,2), 
     1                     ICS   , DELT_SOIL,ONE)  
        ENDIF
!
!
! COEF XK
!
        CMARK  = 'XK1'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK(1,1),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK(1,1),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'XK2'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK(1,2),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK(1,2),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'XK3'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK(1,3),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK(1,3),ICS, DELT_SOIL,ONE)  
        ENDIF
!
! COEF XK_H 
!
        CMARK  = 'XKH1'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK_H(1,1),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK_H(1,1),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'XKH2'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK_H(1,2),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK_H(1,2),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'XKH3'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XK_H(1,3),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XK_H(1,3),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'AVK1'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , AVK(1,1),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , AVK(1,1),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'AVK2'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , AVK(1,2),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , AVK(1,2),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'AVK3'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , AVK(1,3),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , AVK(1,3),ICS, DELT_SOIL,ONE)  
        ENDIF
!
        CMARK  = 'POTL1'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,1),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,1),ICS, DELT_SOIL,ONE)  
        ENDIF
        CMARK  = 'POTL2'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,2),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,2),ICS, DELT_SOIL,ONE)  
        ENDIF
        CMARK  = 'POTL3'
        IF ( INSFLAG .EQ.0 ) THEN
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,3),ICN, DELT_SOIL,ONE)  
        ELSE 
          CALL MONIT_SIB ( CMARK , XPOT_SOIL(1,3),ICS, DELT_SOIL,ONE)  
        ENDIF
!
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!
! 水フラックス 1 -> 2 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN12F'
          CALL MONIT_SIB ( CMARK , WFLUX_H    (1,1), 
     1                     ICN, DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS12F'
          CALL MONIT_SIB ( CMARK , WFLUX_H    (1,1), 
     1                     ICS  , DELT_SOIL,ONE)  
        ENDIF
!
! 水フラックス 2 -> 3 
!
        IF ( INSFLAG .EQ.0 ) THEN
          CMARK  = 'WFLN23F'
          CALL MONIT_SIB ( CMARK , WFLUX_H    (1,2), 
     1                     ICN   , DELT_SOIL,ONE)  
        ELSE
          CMARK  = 'WFLS23F'
          CALL MONIT_SIB ( CMARK , WFLUX_H    (1,2), 
     1                     ICS   , DELT_SOIL,ONE)  
        ENDIF


!
      ENDIF
#endif 
!
      RETURN
      END
!===============================================================
      SUBROUTINE MR_SOIL_WFLUX
!
! soil_wflux で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1         IMONIT_LEVEL, IMONIT_ALL, IMONIT_FLUX_COMMON
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
        WRITE(6,*) 'MONIT_REGIST_SOIL_WFLUX 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!-------------------------------------------------------------
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!
!  1->2 の水フラックス   （熱は別途）
!
      CMARK  = 'WFLN12'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 NOS explicit'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS12'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 SNW explicit'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  2->3 の水フラックス   （熱は別途）
!
      CMARK  = 'WFLN23'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 NOS explicit'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS23'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 SNW explicit'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLN34'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 NOS explicit'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS34'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 SNW explicit'  
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLN12C'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 NOS explicit'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS12C'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 SNW explicit'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!  2->3 の水フラックス   （熱は別途）
!
      CMARK  = 'WFLN23C'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 NOS explicit'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS23C' 
      CTITLE = 'H2O FLUX SOIL 2 TO 3 SNW explicit'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
!
      CMARK  = 'XK1'
      CTITLE = 'H2O FLUX COEF 1-2'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'XK2'
      CTITLE = 'H2O FLUX COEF 2-3'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'XK3'
      CTITLE = 'H2O FLUX COEF 3-4'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'XKH1'
      CTITLE = 'H2O FLUX COEF 1-2'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'XKH2'
      CTITLE = 'H2O FLUX COEF 2-3'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'XKH3'
      CTITLE = 'H2O FLUX COEF 3-4'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'AVK1'
      CTITLE = 'H2O FLUX COEF 1-2'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'AVK2'
      CTITLE = 'H2O FLUX COEF 2-3'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      CMARK  = 'AVK3'
      CTITLE = 'H2O FLUX COEF 3-4'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'POTL1'
      CTITLE = 'potential L1'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'POTL2'
      CTITLE = 'potential L2'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'POTL3'
      CTITLE = 'potential L3'
      CUNIT  = '-'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )

!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
      CMARK  = 'WFLN12F'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 NOS after CORRECT'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS12F'
      CTITLE = 'H2O FLUX SOIL 1 TO 2 SNW after CORRECT'
      CUNIT  = 'kg/m**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLN23F'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 NOS after CORRECT'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WFLS23F'
      CTITLE = 'H2O FLUX SOIL 2 TO 3 SNW after CORRECT'
      CUNIT  = 'kg/m**2/s'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END














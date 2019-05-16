! 気孔抵抗の計算
! 99/07/19 保坂征宏 大幅な変更
! 99/07/30 保坂征宏 抵抗で *_bar にする
! 99/08/03          根の抵抗も *_bar にする
!                   蒸散の計算（ポテンシャルの計算）は *_bar ですることに。
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/10/29 保坂征宏 土壌ポテンシャルが低いときにも蒸散するバグ発見。
! 2001 6/18  かなりの手直し  M.Oh'izumi   

      SUBROUTINE CANOPY_STRES (
     I   NSTR        , NEND        , INSFLAG      ,
     I   ITR         , ITR_MAX     , 
     I   IFLG_AFT_MOD_FLG  , IFLG_LAST , IFLG_MOD_FLG  , 
     I   NTYP        , 
     I   R_AC_BAR    , R_AG        ,
     I   E_CNP_SAT   , E_SKIN_SAT  , E_AIR      , E_AIR_SAT      , 
     I   FUNC_TMP    , 
     I   POT_SOIL    ,
     I   D_ZERO_msx  , 
     I   R_ROOT_BAR_INV  ,               CVR_VEG_INV_msx   , 
     I   FAC_EVP     , R_STM_BAR_PAR    , 
     I   XFR_H2O_CNP  , XFR_H2O_GRS  ,
#ifdef MONYOS
     I   DELT_CNP , 
#endif
     O   POT_LEAF    , R_STM_BAR   , 
!    M   FAC_TRNS_A  , FAC_TRNS_B  , KLEV_TRNS_MAX , 
     M   FAC_TRNS_A  , FAC_TRNS_B  , 
     M   FLAG_TR     , FLAG_TR_ALL   )


C    注意：蒸散がある時に、負の符号がつく

!
! 出力・変更
!   POT_LEAF   (IDX,ICG)        葉のポテンシャル
!   R_STM_BAR  (IDX,ICG)        気孔抵抗
!   FAC_TRNS_A (IDX,IDP  ,ICG)  蒸散分配係数
!   FAC_TRNS_B (IDX,IDP  ,ICG)  蒸散分配係数
!   FLAG_TR    (IDX,IDP  ,ICG)  蒸散発生フラグ 1 あり 0 なし
!
! 蒸散が起こり得ない場合としては、
!   ・結露条件を満たす              ... どの土壌層からも蒸散なし
!   ・土壌のポテンシャルが十分低い
!
! このサブルーチンですること
!   蒸散が起こるときの葉のポテンシャルと気孔抵抗を計算する。
!   後に蒸散フラックスが求まったとき、層ごとの分配をするための係数を計算
!
! ややこしいこと
!   蒸散が起こらない層がある場合に、その層がないものとして計算する必要あり。
!
      USE SIBPRM , ONLY :
     1    IDX   ,
     1    ITYP  ,
     1    ICG   ,
     1    IDP
!
      USE SIBCON , ONLY :
     1   RHO_WTR_c ,
     1   HL_WTR_c
!
      USE SIBVEG , ONLY :
     1           DEP_RMEAN_co            , ! 層内根平均深さ  (ITYP,IDP,ICG) 
     1           F_RST_DEFI_co           , ! 気孔抵抗の飽差依存定数 (ITYP,ICG) 
     1           F_RST_POTLMAX_co        , ! 葉のポテンシャル上限 (ITYP,ICG) 
     1           F_RST_POTLMIN_co        , ! 葉のポテンシャル下限 (ITYP,ICG) 
     I           R_PLANT_co                ! 茎の抵抗        (ITYP,ICG)
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS,
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: ITR
      INTEGER,INTENT(IN)   :: ITR_MAX
!
      INTEGER,INTENT(IN)   :: INSFLAG           ! 雪無し 0 雪あり 1 
      INTEGER,INTENT(IN)   :: IFLG_AFT_MOD_FLG  ! FLAG_TR を変更した直後フラグ 
      INTEGER,INTENT(IN)   :: IFLG_LAST         ! 最終フラグ
      INTEGER,INTENT(IN)   :: IFLG_MOD_FLG      ! FLAG_TR を変更フラグ
!
      REAL(8),INTENT(IN)   :: R_AC_BAR (IDX)    ! R_AC(キャ-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_AG     (IDX)    ! R_AG(地面-キャ空間抵抗)
      REAL(8),INTENT(IN)   :: R_ROOT_BAR_INV(IDX,IDP,ICG) 
!                                               ! 各層の(透水＋吸水)抵抗
!
#ifdef MONYOS
      REAL(8),INTENT(IN)   :: DELT_CNP
#endif
!
      REAL(8),INTENT(IN)   :: E_CNP_SAT (IDX)    ! キャノピー温度飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_AIR     (IDX)    ! キャノピー水蒸気圧
      REAL(8),INTENT(IN)   :: E_AIR_SAT (IDX)    ! キャノピー水蒸気圧
      REAL(8),INTENT(IN)   :: E_SKIN_SAT(IDX)    ! 地面温度飽和水蒸気圧
!
!     REAL(8),INTENT(IN)   :: ICHECK  (IDX)     ! 計算する 1 , しない 0
!
      REAL(8),INTENT(IN)   :: FUNC_TMP(IDX,ICG) ! 気孔抵抗の温度依存性
!
      REAL(8),INTENT(IN)   :: POT_SOIL (IDX,IDP)  ! 土壌のmatric pot
      REAL(8),INTENT(IN)   :: D_ZERO_msx (IDX)    ! 零面変位   (旧；Ｄ)
!
!     REAL(8),INTENT(IN)   :: CVR_VEG_ms (IDX,ICG)  ! 植生被覆率
      REAL(8),INTENT(IN)   :: CVR_VEG_INV_msx (IDX,ICG)  ! 植生被覆率
      REAL(8),INTENT(IN)   :: FAC_EVP  (IDX)      !  ρε/P_s 
!                                                 ! 文書の XX 
!
      REAL(8),INTENT(IN)   :: R_STM_BAR_PAR  (IDX,ICG) ! 
!
      INTEGER,INTENT(IN)   :: NTYP   (IDX)       ! 植生インデックス

      REAL(8),INTENT(IN)   :: XFR_H2O_CNP (IDX)   ! キャノピーの濡れ領域の割合
      REAL(8),INTENT(IN)   :: XFR_H2O_GRS (IDX)   ! 下草の濡れ領域の割合
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: R_STM_BAR  (IDX,ICG)    ! 気孔抵抗
      REAL(8),INTENT(OUT)  :: POT_LEAF   (IDX,ICG)    ! 葉の水分ポテンシャル
!
!
! INOUT
!
      REAL(8),INTENT(INOUT)  :: FAC_TRNS_A (IDX,IDP,ICG)  
      REAL(8),INTENT(INOUT)  :: FAC_TRNS_B (IDX,IDP,ICG)  
!     INTEGER,INTENT(INOUT)  :: KLEV_TRNS_MAX(IDX,ICG)  
      REAL(8),INTENT(INOUT)  :: FLAG_TR     (IDX,IDP,ICG)  ! 
      REAL(8),INTENT(INOUT)  :: FLAG_TR_ALL (IDX,    ICG)  ! 
!
! SAVE
!
      REAL(8),SAVE  :: AA (IDX,ICG)             ! 方程式の係数
      REAL(8),SAVE  :: BB (IDX,ICG)
      REAL(8),SAVE  :: SUM_ROOT_INV (IDX,ICG)
      REAL(8),SAVE  :: SUM_POT      (IDX,ICG)
!
! WORK
!
      REAL(8)     :: POT_SOIL_WRK(IDX,IDP,ICG)
      REAL(8)     :: TRNS_SOIL(IDX,IDP,ICG)
      REAL(8)     :: TRNS_ALL (IDX,    ICG)
!
      REAL(8)     :: XFR_DRY  (IDX,ICG)     ! 葉の湿っていない割合
!              
      REAL(8)     :: DIFF_EAX(IDX,ICG)     ! e*(Tc) or e*(Tg) - e(Ta) 旧 EPOT 
      REAL(8)     :: DEFI                  ! 飽差 ( deficit ) 
!
      REAL(8)     :: R_AX_BAR    (IDX,ICG)     ! 抵抗 キャ空気と、葉または地面
      REAL(8)     :: R_STM_BAR_1 (IDX,ICG)     ! 気孔抵抗作業
      REAL(8)     :: R_STM_BAR_ORG1(IDX,ICG)   ! 気孔抵抗作業
      REAL(8)     :: R_STM_BAR_ORG2(IDX,ICG)   ! 気孔抵抗作業
      REAL(8)     :: R_STM_BAR_ORG3(IDX,ICG)   ! 気孔抵抗作業
      REAL(8)     :: R_STM_BAR_PAR1(IDX,ICG)   ! 気孔抵抗作業
      REAL(8)     :: R_STM_work  (IDX)         ! 気孔抵抗作業
!
      REAL(8)     :: CCI                   ! 方程式の係数
      REAL(8)     :: DD (IDX,ICG)
      REAL(8)     :: EE
!
      REAL(8)     :: ACOEF                 ! 二次方程式の係数
      REAL(8)     :: BCOEF
      REAL(8)     :: CCOEF
!
      REAL(8)     :: TERM                  ! 二次方程式の解の公式のルート
! 
      REAL(8)     :: FUNC_DPOT(IDX,ICG)   ! 土壌水分依存性
      REAL(8)     :: FUNC_DPOT_ORG(IDX,ICG)   ! 土壌水分依存性

      REAL(8)     :: FUNC_TMP_X
      REAL(8)     :: FUNC_TMP_MIN
!
      REAL(8)     :: FUNC_DEFI    (IDX,ICG)
      REAL(8)     :: FUNC_DEFI_ORG(IDX,ICG)
      REAL(8)     :: FUNC_DEFI_MIN
!     REAL(8)     :: FUNC_DEFI_X
!
      REAL(8)     :: DIFF               
!     REAL(8)     :: PM
!
      REAL(8),SAVE :: R_STM_BAR_MAX 
!
      REAL(8),SAVE :: POT_MIN_ADD 
      REAL(8),SAVE :: DZERO

      INTEGER      :: IVEG
      REAL(8)      :: FVEG              ! 下草 IVEG=2 では考慮しない
      INTEGER      :: I
!
      REAL(8)      :: FDP_MAX
      REAL(8)      :: FDP_MIN
      REAL(8),PARAMETER   :: SMALL = 1.D-30
!
#ifdef MONYOS 
      CHARACTER(7) :: CMARK 
      REAL(8),PARAMETER :: ONE = 1.D0 
      REAL(8)      :: ONE_WORK(IDX) 
      REAL(8)      :: two_WORK(IDX) 
      INTEGER      :: IFLG_MON
#endif
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      DATA DZERO /0.D0/
      DATA POT_MIN_ADD   / 1.D-5 /
      DATA FUNC_TMP_MIN  / 1.D-3 /              ! ??
      DATA FUNC_DEFI_MIN / 1.D-1 /              ! 01/07/28 =L3SiB 
!
      DATA FDP_MIN  / 1.D-4 /
      DATA FDP_MAX  / 1.D0 /
!      
C------ after Routine SiB   2001 6/25 ----------------

      DATA R_STM_BAR_MAX  / 1.D30 /              ! ??
!      DATA R_STM_BAR_MAX / 1.D5 /              ! ??
!
      r_stm_bar     (:,:) = 1.D50 
      r_stm_bar_par1(:,:) = 1.D50 
      r_stm_bar_org1(:,:) = 1.D50 
      r_stm_bar_org2(:,:) = 1.D50 
      r_stm_bar_org3(:,:) = 1.D50 
!
!!      if ( insflag.eq.0 ) 
!!     1   write(6,*) 'ISTEP canopy_stres' , 4,R_STM_BAR_PAR(4,1) 
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_STRES 99/07 POCHAKA'  
        LFIRST = .FALSE.
      ENDIF
!
#ifdef MONYOS 
      IF ( INSFLAG.EQ.0 ) THEN
        IFLG_MON=ICN
      ELSE
        IFLG_MON=ICS
      ENDIF
#endif

!
!  =====================================================================
!  >> 蒸散可能層のとりかたを変えたとき（ FLAG_TR が変わった後）計算 <<
!  =====================================================================
!
!  FLAG_TR には、蒸散が可能な層では 1 、可能でない層では 0 が格納済。
!  可能な層だけでの、
!      根の抵抗の逆数          の和  R_ROOT_BAR_INV 
!      ポテンシャル差/根の抵抗 の和  SUM_POT 
!  と、
!      方程式の係数 AA と BB   
!  を計算する。
!
      IF ( IFLG_AFT_MOD_FLG .EQ.1 ) THEN
#include "noparallel_fixedform"
        DO 1000 IVEG=1,2-INSFLAG 
!
        DO 1000 I=NSTR,NEND
          IF ( FLAG_TR_ALL(I,IVEG).GT.0.5d0  ) THEN
!
            SUM_ROOT_INV(I,IVEG)  
     1      =   R_ROOT_BAR_INV(I,1,IVEG) * FLAG_TR(I,1,IVEG)
     1        + R_ROOT_BAR_INV(I,2,IVEG) * FLAG_TR(I,2,IVEG)
     1        + R_ROOT_BAR_INV(I,3,IVEG) * FLAG_TR(I,3,IVEG)
!
            SUM_POT(I,IVEG)
     1      =   R_ROOT_BAR_INV(I,1,IVEG) 
     1            * ( POT_SOIL(I,1) - DEP_RMEAN_co(NTYP(I),1,IVEG) )
     1            * FLAG_TR(I,1,IVEG)
     2        + R_ROOT_BAR_INV(I,2,IVEG) 
     2            * ( POT_SOIL(I,2) - DEP_RMEAN_co(NTYP(I),2,IVEG) )
     2            * FLAG_TR(I,2,IVEG)
!    3        + R_ROOT_BAR_INV(I,3,IVEG) 
     3        + R_ROOT_BAR_INV(I,3,IVEG)                              ! bug
     3            * ( POT_SOIL(I,3) - DEP_RMEAN_co(NTYP(I),3,IVEG) )
     3            * FLAG_TR(I,3,IVEG)
!
            FVEG = FLOAT(2-IVEG)  
!
            AA(I,IVEG) 
     1         = ( SUM_POT(I,IVEG) 
     1              - ( FVEG * D_ZERO_msx(I) 
     1                       * SUM_ROOT_INV(I,IVEG) ) ) 
     1            / (   1.0D0 
     1               + R_PLANT_co(NTYP(I),IVEG)
     1                       * SUM_ROOT_INV(I,IVEG) )      

C-------- bug 2001 6/18 MO ------------------------------------------
C       根の抵抗と蒸散ドキュメント（大泉原案、保坂編集、平成１１年７月２６日）
C 　　　に誤り有り。水の密度をかける必要有り


            aa(i,iveg)=rho_wtr_c*aa(i,iveg)
!----------------------------------------------------------------
            BB(I,IVEG) 
     1          = - SUM_ROOT_INV(I,IVEG) 
     1                / ( 1.0D0 
     1                     + R_PLANT_co(NTYP(I),IVEG) 
     1                       * SUM_ROOT_INV(I,IVEG)   )
!----------------------------------------------------------------

            bb(i,iveg)=rho_wtr_c*bb(i,iveg)

!----------------------------------------------------------------

          ENDIF
 1000   CONTINUE
      ENDIF
!
#include "noparallel_fixedform"
      DO 1100 IVEG=1,2-INSFLAG
!
      DO 1100 I=NSTR,NEND
        POT_SOIL_WRK(I,1,IVEG) = POT_SOIL(I,1) 
     1                           - DEP_RMEAN_co(NTYP(I),1,IVEG) 
        POT_SOIL_WRK(I,2,IVEG) = POT_SOIL(I,2) 
     1                           - DEP_RMEAN_co(NTYP(I),2,IVEG) 
        POT_SOIL_WRK(I,3,IVEG) = POT_SOIL(I,3) 
     1                           - DEP_RMEAN_co(NTYP(I),3,IVEG) 
 1100 CONTINUE  
!
!  ======================================
!  >> 結露状態（蒸散無し）のときの対処 <<
!  ======================================
!    結露（e_sat(Tc) < e_air 等のとき ）の場合蒸散は発生せず、
!    全面への結露 = 負の蒸発 になる。
!
!    注意。
!      不用意にこれをすると、本当は結露条件でない場合に、全面からの
!      蒸発が起こってしまう。
!
!    葉の乾き度 XFR_DRY 、抵抗 R_AX_BAR 、ポテンシャルの差
!    
#include "noparallel_fixedform"
      DO 2000 IVEG=1,2-INSFLAG
!
      DO 2000 I=NSTR,NEND
!
!  葉の乾き度、キャノピー空気の飽和度計算
!
        IF ( IVEG.EQ.1 ) THEN
          XFR_DRY  (I,IVEG) = 1.0D0 - XFR_H2O_CNP(I)  
          DIFF_EAX(I,IVEG) = E_CNP_SAT(I)  - E_AIR(I)  
        ELSE
          XFR_DRY  (I,IVEG) = 1.0D0 - XFR_H2O_GRS(I)
          DIFF_EAX(I,IVEG) = E_SKIN_SAT(I) - E_AIR(I) 
        END IF
!
!  全面濡れていたり、結露条件なら蒸散はない
!
        IF(XFR_DRY(I,IVEG).EQ.0.d0.OR.DIFF_EAX(I,IVEG).LE.0.d0) THEN
!
          FLAG_TR    (I,1,IVEG) = 0.d0 
          FLAG_TR    (I,2,IVEG) = 0.d0 
          FLAG_TR    (I,3,IVEG) = 0.d0 
          FLAG_TR_ALL(I  ,IVEG) = 0.d0 
!
          FAC_TRNS_A (I,1,IVEG) = 0.d0
          FAC_TRNS_B (I,1,IVEG) = 0.d0 
          FAC_TRNS_A (I,2,IVEG) = 0.d0
          FAC_TRNS_B (I,2,IVEG) = 0.d0 
          FAC_TRNS_A (I,3,IVEG) = 0.d0
          FAC_TRNS_B (I,3,IVEG) = 0.d0 
!
        ENDIF
 2000 CONTINUE  
!
!  ==========================
!  >> 蒸散ありのときの対処 <<
!  ==========================
!
      FUNC_DEFI     (:,:) = 0.
      FUNC_DEFI_ORG (:,:) = 0.
      FUNC_DPOT     (:,:) = 0.
      FUNC_DPOT_ORG (:,:) = 0.
!
#include "noparallel_fixedform"
      DO 2100 IVEG=1,2-INSFLAG
!
      DO 2100 I=NSTR,NEND
!
        IF ( FLAG_TR_ALL(I,IVEG).GT.0.5d0 ) THEN
!
!  ------------------------------
!  > 気孔抵抗の各種依存性の評価 <
!  ------------------------------
!
!  飽差依存性
!
!Z bug    DEFI        = MAX ( DZERO, E_AIR_SAT(I) - E_AIR_SAT(I) )     ! 飽差
          DEFI        = MAX ( DZERO, E_AIR_SAT(I) - E_AIR(I)     )     ! 飽差
!
C-------コメント  2001 6/18 MO------------------------------------------------
Cーーーーー　　現ルーチン版の以下の式に変更
         FUNC_DEFI(I,IVEG)
     1        = (1.0D0 +  F_RST_DEFI_co(NTYP(I),IVEG) * DEFI )
     &        / (1.0D0 +  3.d0 * F_RST_DEFI_co(NTYP(I),IVEG) * DEFI )
C-----　　以下はオリジナル（ルーチンではコメント化）の定式
C          FUNC_DEFI   = 1.0D0 -  F_RST_DEFI_co(NTYP(I),IVEG) * DEFI 
C------------------------------------------------------------------------------

          FUNC_DEFI_ORG(I,IVEG) = FUNC_DEFI(I,IVEG) 
          FUNC_DEFI(I,IVEG) = MAX( FUNC_DEFI_MIN , FUNC_DEFI(I,IVEG) )
!             L3SiB に従い、下限値処理 ( 0.1 ) 
!
!  温度依存性
!
!Z        FUNC_TMP_X  = MAX( FUNC_TMP_MIN  , FUNC_TMP(I,IVEG) )
          FUNC_TMP_X  = 1.d0 
!
!  飽差・温度依存性・ CVR_VEG を考慮した気孔抵抗
!
          r_stm_bar_par1(I,IVEG) = R_STM_BAR_PAR (I,IVEG)   ! モニタ
!
!!      if ( insflag.eq.0 .and. i.eq.4 .and. iveg.eq. 1) 
!!     1   write(6,*) 'ISTEP canopy_stres_2' , R_STM_BAR_PAR1(i,iveg) 
          
!
          R_STM_BAR_1 (I,IVEG) 
     1        = R_STM_BAR_PAR (I,IVEG) 
     1          / ( FUNC_DEFI(I,IVEG) * FUNC_TMP_X  )
!
          R_STM_BAR_1 (I,IVEG) 
     1          = MIN ( R_STM_BAR_1(I,IVEG) , R_STM_BAR_MAX )
!
          IF ( IVEG.EQ.1 ) THEN
            R_AX_BAR(I,IVEG) = R_AC_BAR(I) 
          ELSE
            R_AX_BAR(I,IVEG) = R_AG(I) * CVR_VEG_INV_msx(I,2)
          END IF
!
!  ----------------------
!  > 二次方程式を立てる <
!  ----------------------
!
          CCI = FAC_EVP(I) * XFR_DRY (I,IVEG)
!
          DD(I,IVEG) 
     1       =   R_STM_BAR_1 (I,IVEG) 
     1              * (   F_RST_POTLMAX_co(NTYP(I),IVEG) 
     1                  - F_RST_POTLMIN_co(NTYP(I),IVEG) )
     1          - R_AX_BAR(I,IVEG) * F_RST_POTLMIN_co(NTYP(I),IVEG)
          EE = - DIFF_EAX(I,IVEG) * F_RST_POTLMIN_co(NTYP(I),IVEG)
!
!  二次方程式の係数設定
!
          ACOEF = BB(I,IVEG) * R_AX_BAR(I,IVEG)           ! BB<0 , R_AX_BAR>0 
          BCOEF = ( - DIFF_EAX(I,IVEG) * CCI 
     1              + ( AA(I,IVEG) * R_AX_BAR(I,IVEG) 
     1                  + BB(I,IVEG) * DD(I,IVEG) ) )      
          CCOEF = AA(I,IVEG) * DD(I,IVEG) - EE * CCI
!
!  解のルートの中身
!
          TERM    = MAX ( DZERO , BCOEF*BCOEF - 4.0D0*ACOEF*CCOEF )
!
!  --------------------
!  > 二次方程式を解く <
!  --------------------
!  解。ただし POT_LEAF が大きい方をとる。　（保坂）

C       2001 6/19 MO memo
C                 　以下の式は代数的に大きい方を採用する方式である。
C                   正値になった場合は、強制的にF_RST_POTLMAX_coまで下げる。
C                   意味のある負値が２個得られた場合は、蒸散させ易い方を
C                   採用することになる。


!       ACOEF < 0 に注意。
!       POT_LEAF < F_RST_POTLMIN_co(NTYP(I),IVEG) なら蒸散は負
!       POT_LEAF > F_RST_POTLMAX_co(NTYP(I),IVEG) なら解き直し
!
          POT_LEAF(I,IVEG) 
     1       = ( - BCOEF - SQRT ( TERM ) ) / ( 2.0D0 * ACOEF )
!

C--------  check quardatic eq. root   -- 2001 6/19 -------------
C          後でコメント化せよ

C       if( POT_LEAF(I,IVEG) .gt.
C     1               F_RST_POTLMAX_co(NTYP(I),IVEG)) then
C           write(*,*) ' - Pot_leaf =', POT_LEAF(I,IVEG) 
C           write(*,*) ' + Pot_leaf =',  
C     1       ( - BCOEF + SQRT ( TERM ) ) / ( 2.0D0 * ACOEF )
C           write(*,*) ' Pot_MAX  =', F_RST_POTLMAX_co(NTYP(I),IVEG) 
C       endif
C---------------------------------------------------------------- 
!
!  F_RST_POTLMIN_co より小さいときは、気孔抵抗が大きくて、実質的に蒸散無し
!
          POT_LEAF(I,IVEG) 
     1       = MAX ( POT_LEAF(I,IVEG) , 
     1               F_RST_POTLMIN_co(NTYP(I),IVEG) + POT_MIN_ADD  )
!
!  F_RST_POTLMAX_co より大きいときは、F_RST_POTLMAX_co にする 
!
          POT_LEAF(I,IVEG) 
     1       = MIN ( POT_LEAF   (I,IVEG) , 
     1               F_RST_POTLMAX_co(NTYP(I),IVEG) )
!
!  気孔抵抗に土壌水ポテンシャル依存性をかける ( f(ψl) の逆数 )
!
!        FUNC_DPOT(I,IVEG)
!    1        =   ( F_RST_POTLMAX_co(NTYP(I),IVEG) 
!    1                  - F_RST_POTLMIN_co(NTYP(I),IVEG) )
!    1           / ( POT_LEAF   (I,IVEG)       
!    1                  - F_RST_POTLMIN_co(NTYP(I),IVEG) ) 
!
!  気孔抵抗に土壌水ポテンシャル依存性をかける ( f(ψl) )
!
         FUNC_DPOT(I,IVEG)
     1        =   
     1             ( POT_LEAF   (I,IVEG)       
     1                  - F_RST_POTLMIN_co(NTYP(I),IVEG) ) 
     1           / ( F_RST_POTLMAX_co(NTYP(I),IVEG) 
     1                  - F_RST_POTLMIN_co(NTYP(I),IVEG) )
!
         FUNC_DPOT_ORG(I,IVEG) = FUNC_DPOT(I,IVEG) 
!
         FUNC_DPOT(I,IVEG) = MIN(FDP_MAX,MAX(FDP_MIN,FUNC_DPOT(I,IVEG)))
!        
         R_STM_BAR(I,IVEG) = R_STM_BAR_1 (I,IVEG) / FUNC_DPOT(I,IVEG) 
         R_STM_BAR_ORG1(I,IVEG) = R_STM_BAR (I,IVEG) 
!
        ENDIF 

 2100 CONTINUE
!
! ====================
! >> 蒸散量の仮評価 <<
! ====================
!   蒸散量を仮に評価する。
!   その蒸散量が負である場合、
!   土壌層からの蒸散はないと判定して FLAG_TR = 0 にする
!
      IF ( IFLG_MOD_FLG.EQ.1 .OR. IFLG_LAST.EQ.1 ) THEN
        CALL RESET ( TRNS_SOIL,IDX*IDP*ICG )
!
#include "noparallel_fixedform"
        DO 3000 IVEG=1,2-INSFLAG 
!
        DO 3000 I=NSTR,NEND

          IF ( FLAG_TR_ALL(I,IVEG).GT.0.5d0  ) THEN
!
! --------------------------------
! > 土壌層ごとの蒸散量を評価する <
! --------------------------------
!      FLAG_FR を調べるために計算
!
            TRNS_ALL(I,IVEG)
     1        = - FAC_EVP(I) * ( E_CNP_SAT(I)-E_AIR(I) ) 
     1                       * XFR_DRY(I,IVEG)  
     1            / ( R_STM_BAR(I,IVEG) + R_AX_BAR(I,IVEG) ) 
!
!  ---------------------------------
!  > 層毎の蒸散量 TRNS_SOIL の評価 <
!  ---------------------------------
!      FLAG_TR を調べるために計算
!      注意。
!        FAC_TRNS_A, FAC_TRNS_B は後に蒸散合計から各層毎の
!        蒸散量分配を評価するのに使う。
!        A の三層での和は 0 , B の和は 1 になるはずで、
!        各層での蒸散の和は全層の蒸散に等しい。
!        ... いずれは、正確に等しくなるようにするために、もっとも蒸散が
!            大きい層の A,B を他から求めるように調節しよう。
!            ずれがある場所については警告。
! 
!!        現時点のプログラムでは、FLAG_TR がこの後で変更されるし、
!!        canopy_air の FC,FG との整合性もきちんととっていないので、
!!        整合的ではない。とりあえずの処置である。
!
            FAC_TRNS_A(I,1,IVEG) 
     1        =  ( POT_SOIL_WRK(I,1,IVEG) 
     1             - SUM_POT(I,IVEG) /( SUM_ROOT_INV(I,IVEG)+SMALL) ) 
     1           * R_ROOT_BAR_INV(I,1,IVEG)
     1           * FLAG_TR(I,1,IVEG)    
!
C-------- bug 2001 6/18 MO --------------------------------------------------
C       根の抵抗と蒸散ドキュメント（大泉原案、保坂編集、平成１１年７月２６日）
C 　　　に誤り有り。水の密度をかける必要有り


            FAC_TRNS_A(I,1,IVEG) = FAC_TRNS_A(I,1,IVEG)*rho_wtr_c

C-------------------------------------------------------------------
            FAC_TRNS_A(I,2,IVEG) 
     1        =  ( POT_SOIL_WRK(I,2,IVEG) 
     1             - SUM_POT(I,IVEG) / (SUM_ROOT_INV(I,IVEG)+SMALL) ) 
     1           * R_ROOT_BAR_INV(I,2,IVEG)
     1           * FLAG_TR(I,2,IVEG)    
!
C----------------------------------------------------------------


            FAC_TRNS_A(I,2,IVEG) = FAC_TRNS_A(I,2,IVEG)*rho_wtr_c

C------------------------------------------------------------------
             FAC_TRNS_A(I,3,IVEG) 
     1        =  ( POT_SOIL_WRK(I,3,IVEG) 
     1             - SUM_POT(I,IVEG) / (SUM_ROOT_INV(I,IVEG)+SMALL) ) 
     1           * R_ROOT_BAR_INV(I,3,IVEG)
     1           * FLAG_TR(I,3,IVEG)    
!
C-------------------------------------------------------------------

            FAC_TRNS_A(I,3,IVEG) = FAC_TRNS_A(I,3,IVEG)*rho_wtr_c

C-------------------------------------------------------------------
#ifdef SIB_DEBUG
            DIFF =     FAC_TRNS_A(I,1,IVEG) 
     1               + FAC_TRNS_A(I,2,IVEG) 
     1               + FAC_TRNS_A(I,3,IVEG) 
            IF ( ABS(DIFF).GT.1.D-15 ) THEN
              write(*,*) 'canopy_stress warning diff is too large A' , 
     1                 i , IVEG ,  DIFF ,
     1                 FAC_TRNS_A(I,1,IVEG) ,
     1                 FAC_TRNS_A(I,2,IVEG) ,
     1                 FAC_TRNS_A(I,3,IVEG) ,
     1                 POT_SOIL_WRK(I,1,IVEG) ,
     1                 POT_SOIL_WRK(I,2,IVEG) ,
     1                 POT_SOIL_WRK(I,3,IVEG) ,
     1                 - SUM_POT(I,IVEG) /(SUM_ROOT_INV(I,IVEG)+SMALL),
     1                   SUM_POT(I,IVEG) , SUM_ROOT_INV(I,IVEG)  ,
     1                 R_ROOT_BAR_INV(I,1,IVEG) ,
     1                 R_ROOT_BAR_INV(I,2,IVEG) ,
     1                 R_ROOT_BAR_INV(I,3,IVEG)
            ENDIF
            IF ( ABS(DIFF).GT.1.D-1 ) THEN
              write(*,*) 'canopy_stress error diff is too large A'  
              stop 999
            ENDIF
#endif
!
            FAC_TRNS_B(I,1,IVEG) 

C----- bug fixでBBに水の密度をかけているので、もう１度元に戻す
C     1        =   - ( 1/BB(I,IVEG) + R_PLANT_co(NTYP(I),IVEG) ) 
     1        = ( - rho_wtr_c/(BB(I,IVEG)+small)
     1            - R_PLANT_co(NTYP(I),IVEG) ) 
     1           * R_ROOT_BAR_INV(I,1,IVEG)   
     1           * FLAG_TR(I,1,IVEG)    
!
            FAC_TRNS_B(I,2,IVEG) 
C     1        =   - ( 1/BB(I,IVEG) + R_PLANT_co(NTYP(I),IVEG) ) 
     1        = ( - rho_wtr_c/(BB(I,IVEG)+small)
     1            - R_PLANT_co(NTYP(I),IVEG) ) 
     1           * R_ROOT_BAR_INV(I,2,IVEG)   
     1           * FLAG_TR(I,2,IVEG)    
!
            FAC_TRNS_B(I,3,IVEG) 
C     1        =   - ( 1/BB(I,IVEG) + R_PLANT_co(NTYP(I),IVEG) ) 
     1        = ( - rho_wtr_c/(BB(I,IVEG)+small)
     1            - R_PLANT_co(NTYP(I),IVEG) ) 
     1           * R_ROOT_BAR_INV(I,3,IVEG)   
     1           * FLAG_TR(I,3,IVEG)    
!
#ifdef SIB_DEBUG
            DIFF =     FAC_TRNS_B(I,1,IVEG) 
     1               + FAC_TRNS_B(I,2,IVEG) 
     1               + FAC_TRNS_B(I,3,IVEG) - 1.D0 
!
            IF (       ABS(DIFF)     .GT.1.D-10 
     1           .AND. ABS(DIFF+1.D0).GT.1.D-10 ) THEN
               write(*,*) 'canopy_stress warning diff is too large B' , 
     1                 i , IVEG ,  DIFF ,
     1                 FAC_TRNS_B(I,1,IVEG) ,
     1                 FAC_TRNS_B(I,2,IVEG) ,
     1                 FAC_TRNS_B(I,3,IVEG) ,
     1                 rho_wtr_c / (BB(I,IVEG)+small) ,  
     1                 rho_wtr_c , BB(I,IVEG) ,  
     1                 R_PLANT_co(NTYP(I),IVEG)  ,
     1                 R_ROOT_BAR_INV(I,1,IVEG)  , 
     1                 R_ROOT_BAR_INV(I,2,IVEG)  , 
     1                 FLAG_TR(I,1,IVEG)         ,
     1                 FLAG_TR(I,2,IVEG)         
            ENDIF
            IF (       ABS(DIFF)     .GT.1.D-1
     1           .AND. ABS(DIFF+1.D0).GT.1.D-1 ) THEN
              write(*,*) 'canopy_stress error diff is too large B'  
              stop 999
            ENDIF

#endif
!
            TRNS_SOIL(I,1,IVEG) 
     1        = - FAC_TRNS_A(I,1,IVEG) 
     1          + FAC_TRNS_B(I,1,IVEG) * TRNS_ALL(I,IVEG) 
!
            TRNS_SOIL(I,2,IVEG) 
     1        = - FAC_TRNS_A(I,2,IVEG) 
     1          + FAC_TRNS_B(I,2,IVEG) * TRNS_ALL(I,IVEG) 
!
            TRNS_SOIL(I,3,IVEG) 
     1        = - FAC_TRNS_A(I,3,IVEG) 
     1          + FAC_TRNS_B(I,3,IVEG) * TRNS_ALL(I,IVEG) 

          ENDIF
 3000   CONTINUE     
!
      ENDIF
!
!  ====================
!  >> FLAG_TR の変更 <<
!  ====================
!
!   FLAG_TR の変更
!
      IF ( IFLG_MOD_FLG.EQ.1 .AND. IFLG_LAST.NE.1 ) THEN
#include "noparallel_fixedform"
        DO 4000 IVEG=1,2-INSFLAG 
!
        DO 4000 I=NSTR,NEND
          IF ( FLAG_TR_ALL(I,IVEG).GT.0.5d0 ) THEN

           IF (      TRNS_SOIL(I,1,IVEG).GT.0.D0 
!   以下の部分の意味が不明。
!  ３層の根の間で水のやりとりをしないようにということか？
     1       .OR. POT_LEAF(I,IVEG)-POT_SOIL_WRK(I,1,IVEG).GT.100.d0
     1          ) THEN
              FLAG_TR   (I,1,IVEG) = 0.d0 
              FAC_TRNS_A(I,1,IVEG) = 0.d0
              FAC_TRNS_B(I,1,IVEG) = 0.d0 
            ENDIF

            IF (      TRNS_SOIL(I,2,IVEG).GT.0.D0 
     1        .OR. POT_LEAF(I,IVEG)-POT_SOIL_WRK(I,2,IVEG).GT.100.d0
     1          ) THEN
              FLAG_TR   (I,2,IVEG) = 0.d0 
              FAC_TRNS_A(I,2,IVEG) = 0.d0
              FAC_TRNS_B(I,2,IVEG) = 0.d0 
            ENDIF

            IF (      TRNS_SOIL(I,3,IVEG).GT.0.D0
     1        .OR. POT_LEAF(I,IVEG)-POT_SOIL_WRK(I,3,IVEG).GT.100.d0
     1          ) THEN
              FLAG_TR   (I,3,IVEG) = 0.d0 
              FAC_TRNS_A(I,3,IVEG) = 0.d0
              FAC_TRNS_B(I,3,IVEG) = 0.d0 
            ENDIF

            FLAG_TR_ALL(I,IVEG)                  
     1         = 1.d0 -  (1.d0-FLAG_TR(I,1,IVEG)) 
     1                 * (1.d0-FLAG_TR(I,2,IVEG)) 
     1                 * (1.d0-FLAG_TR(I,3,IVEG)) 
          ENDIF
 4000   CONTINUE     
!
!
      ENDIF

!
!  =========================================================
!  >> どの層からも蒸散がない場合は FAC_TRNS_* を零にする。<<
!  =========================================================
!
!
#include "noparallel_fixedform"
      DO IVEG=1,2-INSFLAG 
!
      DO I=NSTR,NEND
!
        if (       abs(flag_tr_all(I,IVEG)) .GT. 0.001 
     1       .AND. abs(flag_tr_all(I,IVEG)-1) .GT. 0.001 )then
          write(6,*) 'canopy_stres flag_tr_all error' , i , iveg, 
     1                flag_tr_all(I,IVEG),
     1                 FLAG_TR(I,1,IVEG), 
     1                 FLAG_TR(I,2,IVEG), 
     1                 FLAG_TR(I,3,IVEG) 
          stop 1    
        endif
#ifdef SX6_OPTION
      end do
      end do        
      DO IVEG=1,2-INSFLAG 
      DO I=NSTR,NEND
#endif
!
        R_STM_BAR_ORG2(I,IVEG) = R_STM_BAR (I,IVEG) 
        IF ( FLAG_TR_ALL(I,IVEG).LT.0.5d0  ) THEN
          FAC_TRNS_A(I,1,IVEG) = 0.d0
          FAC_TRNS_B(I,1,IVEG) = 0.d0 
          FAC_TRNS_A(I,2,IVEG) = 0.d0
          FAC_TRNS_B(I,2,IVEG) = 0.d0 
          FAC_TRNS_A(I,3,IVEG) = 0.d0
          FAC_TRNS_B(I,3,IVEG) = 0.d0 
          R_STM_BAR (I,IVEG)   = R_STM_BAR_MAX          
          R_STM_BAR_ORG3(I,IVEG) = R_STM_BAR (I,IVEG) 
#ifdef MONYOS
          if ( iveg.eq.1 ) two_work(i) = 1. 
#endif
        ELSE 
          R_STM_BAR_ORG3(I,IVEG) = R_STM_BAR (I,IVEG) 
#ifdef MONYOS
          if ( iveg.eq.1 ) two_work(i) = 0. 
#endif
        ENDIF
      ENDDO
      ENDDO 
!
! -----------------------------------------------
! > FAC_TRNS_* があまりにもひどくないことの確認 <
! -----------------------------------------------
!
#ifdef SIB_DEBUG
      IF ( IFLG_LAST.EQ.1 ) THEN
#include "noparallel_fixedform"
      DO IVEG=1,2-INSFLAG
!
      DO I=NSTR,NEND
        DIFF = (   FAC_TRNS_A(I,1,IVEG) + FAC_TRNS_A(I,2,IVEG) 
     1           + FAC_TRNS_A(I,3,IVEG) ) * FLAG_TR_ALL(I,IVEG)
        IF ( ABS(DIFF) .GT. 1.D-8 ) THEN 
          write(6,*) 'canopy_stres sum fac_trns_a' , i, IVEG, diff 
        ENDIF
!
        DIFF = (   FAC_TRNS_B(I,1,IVEG) + FAC_TRNS_B(I,2,IVEG) 
     1           + FAC_TRNS_B(I,3,IVEG) -1.D0 ) * FLAG_TR_ALL(I,IVEG)
        IF ( ABS(DIFF) .GT. 1.D-12 ) THEN 
          write(6,*) 'canopy_stres sum fac_trns_b' , i, IVEG, diff,
     1                FAC_TRNS_B(I,1,IVEG) , FAC_TRNS_B(I,2,IVEG) ,
     1                FAC_TRNS_B(I,3,IVEG) , FLAG_TR_ALL(I,IVEG)
        ENDIF
      ENDDO
      ENDDO
      ENDIF
#endif
!
!  ============================
!  >> 各層からの蒸散の仮評価 <<
!  ============================
!    ただの目安の評価。
!    目安の評価が実際のそれらとどれだけあっているのかの、確認用。
!
#ifdef MONYOS
      IF ( IFLG_LAST.EQ.1  .AND. ITR.EQ.ITR_MAX ) THEN 
      IF ( IMONIT_LEVEL.GE.IMONIT_ALL ) THEN
!
! 蒸散仮評価
!
        CMARK = 'WTRCL1Z' 
        CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,1,1),IFLG_MON,DELT_CNP,ONE ) 
        CMARK = 'WTRCL2Z' 
        CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,2,1),IFLG_MON,DELT_CNP,ONE ) 
        CMARK = 'WTRCL3Z' 
        CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,3,1),IFLG_MON,DELT_CNP,ONE ) 
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'WTRGL1Z' 
          CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,1,2),IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'WTRGL2Z' 
          CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,2,2),IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'WTRGL3Z' 
          CALL MONIT_SIB ( CMARK,TRNS_SOIL(1,3,2),IFLG_MON,DELT_CNP,ONE) 
        ENDIF
!
! 気孔抵抗
!
        CMARK = 'RSTSRC'  ! 短波依存性
        r_stm_work(:) = 1/R_STM_BAR_PAR1(:,1)
        CALL MONIT_SIB ( CMARK,R_STM_work,IFLG_MON,DELT_CNP,ONE) 
!
!!        if ( insflag.eq.0 ) 
!!     1   write(6,*) 'ISTEP canopy_stres_3' , 4 , 
!!     1              R_STM_BAR_PAR1(4,1),R_STM_work(4) 
!
        CMARK = 'RSTECC'  ! 飽差依存性
        CALL MONIT_SIB ( CMARK,FUNC_DEFI(1,1),IFLG_MON,DELT_CNP,ONE ) 
        CMARK = 'RSTEOC'  ! 飽差依存性
        CALL MONIT_SIB ( CMARK,FUNC_DEFI_ORG(1,1),IFLG_MON,DELT_CNP,ONE) 
        CMARK = 'RSTPCC'  ! 土壌ポテンシャル依存性
        CALL MONIT_SIB ( CMARK,FUNC_DPOT(1,1),IFLG_MON,DELT_CNP,ONE ) 
        CMARK = 'RSTPOC'  ! 土壌ポテンシャル依存性
        CALL MONIT_SIB ( CMARK,FUNC_DPOT_ORG(1,1),IFLG_MON,DELT_CNP,ONE) 
        CMARK = 'RSTOC'   ! 全体
        r_stm_work(:) = 1/R_STM_BAR_ORG1(:,1)
        CALL MONIT_SIB ( CMARK,R_STM_work,        IFLG_MON,DELT_CNP,ONE) 
        CMARK = 'RSTNC'   ! 全体
        r_stm_work(:) = 1/R_STM_BAR_ORG2(:,1)
        CALL MONIT_SIB ( CMARK,R_STM_work,        IFLG_MON,DELT_CNP,ONE) 
        CMARK = 'RSTMC'   ! 全体
        r_stm_work(:) = 1/R_STM_BAR_ORG3(:,1)
        CALL MONIT_SIB ( CMARK,R_STM_work,        IFLG_MON,DELT_CNP,ONE) 
        CMARK = 'RSTCC'   ! 全体
        r_stm_work(:) = 1/R_STM_BAR(:,1)
        CALL MONIT_SIB ( CMARK,R_STM_work,IFLG_MON,DELT_CNP,ONE) 
!
        IF ( INSFLAG.EQ.0 ) THEN
          CMARK = 'RSTSRG'
          r_stm_work(:) = 1/R_STM_BAR_PAR1(:,2)
          CALL MONIT_SIB ( CMARK,R_STM_work,IFLG_MON,DELT_CNP,ONE)      
          CMARK = 'RSTECG'
          CALL MONIT_SIB ( CMARK,FUNC_DEFI(1,2),IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTEOG'
          CALL MONIT_SIB ( CMARK,FUNC_DEFI_ORG(1,2),
     1                                          IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTPCG'
          CALL MONIT_SIB ( CMARK,FUNC_DPOT(1,2),IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTPOG'
          CALL MONIT_SIB ( CMARK,FUNC_DPOT_ORG(1,2),
     1                                          IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTOG'
          r_stm_work(:) = 1/R_STM_BAR_ORG1(:,2)
          CALL MONIT_SIB ( CMARK,R_STM_work,
     1                                          IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTNG'
          r_stm_work(:) = 1/R_STM_BAR_ORG2(:,2)
          CALL MONIT_SIB ( CMARK,R_STM_work,
     1                                          IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTMG'
          r_stm_work(:) = 1/R_STM_BAR_ORG3(:,2)
          CALL MONIT_SIB ( CMARK,R_STM_work,
     1                                          IFLG_MON,DELT_CNP,ONE) 
          CMARK = 'RSTCG'
          r_stm_work(:) = 1/R_STM_BAR(:,2)
          CALL MONIT_SIB ( CMARK,R_STM_work,IFLG_MON,DELT_CNP,ONE) 
        ENDIF
!
!
      ENDIF
        IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          CMARK = 'A_CCOND' 
          r_stm_work(:) = 1/R_STM_BAR(:,1)
          CALL MONIT_SIB ( CMARK,R_STM_work,IFLG_MON,DELT_CNP,ONE) 
        ENDIF 
      ENDIF
#endif
!
      RETURN
      END
!=====================================================
      SUBROUTINE FLAG_TR_RESET (
     I   NSTR     , NEND       , INSFLAG , 
     I   TMP_SOIL , 
!    M   FLAG_TR  , FLAG_TR_ALL, KLEV_TRNS_MAX )
     M   FLAG_TR  , FLAG_TR_ALL )
!
      USE SIBPRM , ONLY :
     1   IDX   ,
     1   ICG   , 
     1   IDP
!
      USE SIBCON , ONLY :
     1   TMP_FREZ_c    
!
      IMPLICIT NONE 
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG 
      REAL(8),INTENT(IN)    :: TMP_SOIL(IDX,IDP)
!
! INOUT
!
      REAL(8),INTENT(INOUT) :: FLAG_TR       (IDX,IDP,ICG)   
      REAL(8),INTENT(INOUT) :: FLAG_TR_ALL   (IDX    ,ICG)   
!     INTEGER,INTENT(INOUT) :: KLEV_TRNS_MAX (IDX    ,ICG)   
!
! WORK
!
      INTEGER  :: I
      INTEGER  :: IVEG
      REAL(8),PARAMETER :: ONE =1.D0
!
#include "noparallel_fixedform"
      DO 1000 IVEG=1,2-INSFLAG 
!
      DO 1000 I=NSTR,NEND
        FLAG_TR      (I,1,IVEG) = 1.d0 
        FLAG_TR      (I,2,IVEG) = 1.d0 
        FLAG_TR      (I,3,IVEG) = 1.d0 
        FLAG_TR_ALL  (I,IVEG)   = 1.d0 
!       KLEV_TRNS_MAX(I,IVEG)   = 1          ! 単なる初期値
 1000 CONTINUE  
!
! ------------------------------------
! > おまけ。FLAG_TR に初期値を与える <
! ------------------------------------
!  土中温度が零度以下は蒸散なし
!
      IF ( INSFLAG .EQ. 0 ) THEN
        DO 2000 I=NSTR,NEND 
          FLAG_TR(I,1,1) = 1.d0    
          FLAG_TR(I,2,1) = 1.d0 
          FLAG_TR(I,3,1) = 1.d0 
          FLAG_TR(I,1,2) = 1.d0 
          FLAG_TR(I,2,2) = 1.d0 
          FLAG_TR(I,3,2) = 1.d0 
          IF ( TMP_SOIL(I,1)   .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,1,1) = 0.d0 
            FLAG_TR(I,1,2) = 0.d0 
          ENDIF
          IF ( TMP_SOIL(I,2) .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,2,1) = 0.d0 
            FLAG_TR(I,2,2) = 0.d0 
          ENDIF
          IF ( TMP_SOIL(I,3) .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,3,1) = 0.d0 
            FLAG_TR(I,3,2) = 0.d0 
          ENDIF
          FLAG_TR_ALL(I,1) =  
     1          MIN ( FLAG_TR(I,1,1) + FLAG_TR(I,2,1) + FLAG_TR(I,3,1) 
     1               , ONE )
          FLAG_TR_ALL(I,2) =  
     1          MIN ( FLAG_TR(I,1,2) + FLAG_TR(I,2,2) + FLAG_TR(I,3,2) 
     1               , ONE )
 2000   CONTINUE  
      ELSE                            ! 雪では下草からの蒸散は無し
        DO 2100 I=NSTR,NEND 
          FLAG_TR  (I,1,1) = 1.d0 
          FLAG_TR  (I,2,1) = 1.d0 
          FLAG_TR  (I,3,1) = 1.d0 
          IF ( TMP_SOIL(I,1) .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,1,1) = 0.d0 
          ENDIF
          IF ( TMP_SOIL(I,2) .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,2,1) = 0.d0 
          ENDIF
          IF ( TMP_SOIL(I,3) .LE. TMP_FREZ_c ) THEN 
            FLAG_TR(I,3,1) = 0.d0 
          ENDIF
          FLAG_TR_ALL(I,1) = 
     1          MIN ( FLAG_TR(I,1,1) + FLAG_TR(I,2,1) + FLAG_TR(I,3,1) 
     1               , ONE )
!
          FLAG_TR  (I,1,2) = 0.d0 
          FLAG_TR  (I,2,2) = 0.d0 
          FLAG_TR  (I,3,2) = 0.d0 
          FLAG_TR_ALL(I,2) = 0.d0 
 2100   CONTINUE  
      ENDIF
!
      RETURN
      END
!================================================================
      SUBROUTINE MR_CANOPY_STRES
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL     , IMONIT_ALMA
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
        WRITE(6,*) 'MONIT_REGIST_CANOPY_STRES 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL.GE.IMONIT_ALL ) THEN
!
!  蒸散仮評価
!
      CMARK  = 'WTRCL1Z'
      CTITLE = 'WATER TRNS FROM SOIL1 VIA CANOPY (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK  = 'WTRCL2Z'
      CTITLE = 'WATER TRNS FROM SOIL2 VIA CANOPY (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK  = 'WTRCL3Z'
      CTITLE = 'WATER TRNS FROM SOIL3 VIA CANOPY (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK  = 'WTRGL1Z'
      CTITLE = 'WATER TRNS FROM SOIL1 VIA GRASS (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK  = 'WTRGL2Z'
      CTITLE = 'WATER TRNS FROM SOIL2 VIA GRASS (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK  = 'WTRGL3Z'
      CTITLE = 'WATER TRNS FROM SOIL3 VIA GRASS (KARI)'
      CUNIT  = 'KG/M**2/S'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
! 気孔抵抗 ( 雪あり無し平均 )
!
!   短波依存性
!
      CMARK = 'RSTSRG'
      CTITLE = 'registance srad grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTSRC'
      CTITLE = 'registance srad canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!   飽差依存性
!
      CMARK = 'RSTECG'
      CTITLE = 'registance edef grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTECC'
      CTITLE = 'registance edef canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTEOG'
      CTITLE = 'registance edef grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTEOC'
      CTITLE = 'registance edef canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!   土壌水分ポテンシャル依存性
!
      CMARK = 'RSTPCG'
      CTITLE = 'registance potdef grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTPCC'
      CTITLE = 'registance potdef canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTPOG'
      CTITLE = 'registance potdef grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTPOC'
      CTITLE = 'registance potdef canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
!   気孔抵抗
!
      CMARK = 'RSTCG'
      CTITLE = 'registance grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTCC'
      CTITLE = 'registance canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTOG'
      CTITLE = 'registance grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTOC'
      CTITLE = 'registance canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTNG'
      CTITLE = 'registance grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTNC'
      CTITLE = 'registance canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTMG'
      CTITLE = 'registance grass'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      CMARK = 'RSTMC'
      CTITLE = 'registance canopy'
      CUNIT  = '1'
      CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_CCOND' 
        CTITLE = 'ALMA CCond (Maybe) : registance stomatol (canopy)'
        CUNIT  = '1'
        CALL MONIT_REGIST_SIB ( CMARK , CTITLE , CUNIT )
      ENDIF 

#endif
!
      RETURN
      END

! 根関連の抵抗計算
! 99/07/19 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_ROOT (
     I   NSTR            , NEND          , INSFLAG     ,
     I   NTYP            , MON           , 
     I   FR_WTR_SOIL     , 
     I   FR_GLA_SOIL     , 
     I   TMP_SOIL        , 
     O   POT_SOIL        , R_ROOT_BAR_INV, 
     O   HUM_GRSK        , R_BB        )   ! 雪 INSFLAG=1 では 1/0
!
! LOOP 中一度行う。
!
      USE SIBPRM , ONLY :
     1   IDX   ,
     1   ITYP  ,
     1   ICG   ,
     1   IDP
!
      USE SIBCON , ONLY :
     1   GRAV_c  , 
     1   R_VAP_c 
!
      USE SIBVEG , ONLY :
     1           POT_SAT_SOIL_co        ,  ! 飽和時の土壌ポテンシャル (ITYP)
     1           BEE_co                 ,  ! 土壌水分依存性係数 (ITYP)
     1           WCOEF_SAT_SOIL_co      ,  ! 飽和透水係数    (ITYP)
     1           DEP_ROOT_co            ,  ! 層内根の長さ    (ITYP,IDP,ICG) 
     1           ALPHAF_BAR_mo          ,  ! 根抵抗関連係数α(ITYP,IMON,ICG) 
     1           BETAF_BAR_mo              ! 根抵抗関連係数α(ITYP,IMON,ICG) 
!
      IMPLICIT NONE
!
! 機能 :根の抵抗計算、土壌・CG の水分ポテンシャルの計算。
!
! 出力
!
      REAL(8),INTENT(OUT) :: POT_SOIL(IDX,IDP)
!                               ! I層目の土壌マトリックポテンシャル
      REAL(8),INTENT(OUT) :: R_ROOT_BAR_INV (IDX,IDP,ICG)
!                               ! 各層のC,G 別の(根・透水）抵抗の和の逆数
!
      REAL(8)  :: XFR_WTR
!                               ! 土壌含水率：飽和度（全層） : 水/(1-土-氷) 
!
      REAL(8),INTENT(OUT) :: HUM_GRSK    (IDX) 
!                               ! 雪無し土壌第一層相対湿度 旧 HRR 
      REAL(8),INTENT(OUT) :: R_BB        (IDX) 
!                               ! 雪無し土壌、土壌表層の抵抗
!                                                   
! 入力
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: MON 
!
      INTEGER,INTENT(IN)   :: INSFLAG          ! 雪無し域 0 雪あり域 1
!
      REAL(8),INTENT(IN)   :: FR_WTR_SOIL(IDX,IDP) ! 土壌含水率：飽和度（全層）
      REAL(8),INTENT(IN)   :: FR_GLA_SOIL(IDX,IDP) ! 土壌含水率：飽和度（全層）
!
      REAL(8),INTENT(IN)   :: TMP_SOIL   (IDX)     ! 土壌表層温度（雪無し）
!
      INTEGER,INTENT(IN)   :: NTYP  (IDX)          ! 植生タイプ
!
!
! WORK
!
      REAL(8)   :: WWTR               ! 作業時の土壌含水率
      REAL(8)   :: XPOT_SOIL          ! 作業時の土壌ポテンシャル
      REAL(8)   :: WCOEF_SOIL(IDX,IDP)
      INTEGER   :: I
      INTEGER   :: J
      INTEGER   :: K
      INTEGER   :: IVEG
!
      REAL(8)   :: DMIN_FRW 
      REAL(8)   :: DMIN_SLC
!
      REAL(8),save   :: ONE = 1.D0 
!
#ifdef SIB_DEBUG
      integer :: imsg_num
      real(8) :: xmsg_sq
#endif

!
      LOGICAL, SAVE :: LFIRST 
      DATA LFIRST  / .TRUE. / 
!

      DATA DMIN_FRW / 1.D-3 /
      DATA DMIN_SLC / 1.D-10/
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY__ROOT 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF 
!
!  ========================================
! >> 土壌のマトリックポテンシャル PHSOIL <<
!  ========================================
!                             1-2-81 に従う
!                             ただし評価に使う W は下限 0.1m に設定。
!
      DO 1000 K=1,IDP
      DO 1000 I=NSTR,NEND
        IF (FR_GLA_SOIL(I, K) == 1.0) THEN
          XFR_WTR = 0.0d0
        ELSE
          XFR_WTR         = FR_WTR_SOIL(I,K) / (1.-FR_GLA_SOIL(I,K)) 
        ENDIF
        WWTR            = MAX ( DMIN_FRW , XFR_WTR )
        POT_SOIL(I,K)   = POT_SAT_SOIL_co(NTYP(I)) 
     1                        * WWTR**(-BEE_co(NTYP(I)))
        WCOEF_SOIL(I,K) = WCOEF_SAT_SOIL_co(NTYP(I))
     1                        * WWTR**(2.0D0*BEE_co(NTYP(I))+3.0D0)    
        WCOEF_SOIL(I,K) = MAX ( DMIN_SLC , WCOEF_SOIL(I,K) )
 1000 CONTINUE  
!
! -----------------------------
! > 根の抵抗の逆数 R_ROOT_INV <                           ! 土 + 根の表面
! -----------------------------
!
#include "indep_fixedform"
      DO 2000 IVEG=1,2-INSFLAG                 ! 植生ごと
      DO 2000 J=1,IDP                          ! 各層ごと
      DO 2000 I=NSTR,NEND                      ! 各格子ごと
         R_ROOT_BAR_INV (I,J,IVEG) 
     1        = DEP_ROOT_co(NTYP(I),J,IVEG)
     1          / (   BETAF_BAR_mo (NTYP(I),MON,IVEG) 
     1              + ALPHAF_BAR_mo(NTYP(I),MON,IVEG) 
     1                  / WCOEF_SOIL(I,J) ) 
 2000 CONTINUE
!
#ifdef SIB_DEBUG
      imsg_num = 0 
#endif
!
      IF ( INSFLAG.EQ.0 ) THEN
      DO 3000 I=NSTR,NEND
        WWTR  = MAX ( DMIN_FRW , 
     1                FR_WTR_SOIL(I,1) + FR_GLA_SOIL(I,1) )
!
! 雪無し土壌表層の空気の相対湿度。
!     1-2-130 式の上の式.
!     ここでの土壌水分ポテンシャルは、氷を土とみなす POT_SOIL とは異なる
!     ことに注意。
!     水・氷共存の場合のみならず、氷の場合も昇華の計算に使うことになる。

        XPOT_SOIL    = POT_SAT_SOIL_co(NTYP(I)) 
     1                    * WWTR**(-BEE_co(NTYP(I)))
!20051209
        if ( ISNAN(TMP_SOIL(I)) ) then
	     HUM_GRSK (I) = 1

	  else   
	      HUM_GRSK (I) = DEXP ( XPOT_SOIL * GRAV_c 
     1                       / ( TMP_SOIL(I) * R_VAP_c ) ) 
	  endif
!
#ifdef SIB_DEBUG
        if ( hum_grsk(i).gt.1 ) then
          imsg_num = imsg_num + 1
          xmsg_sq  = xmsg_sq + (hum_grsk(i)-1.)**2   
        endif
#endif
!
        hum_grsk(i) = min ( hum_grsk(i) , one ) 
!
! 雪無し土壌表面での抵抗
!      1-2-129 式の上の式
!
        R_BB (I) = 30. + 3.5 * WWTR**(-2.3) 
 3000 CONTINUE   
      ELSE 
        HUM_GRSK(NSTR:NEND) = 1. 
        R_BB    (NSTR:NEND) = 0. 
      ENDIF
!
#ifdef SIB_DEBUG
      if ( imsg_num.gt.0 ) then
         write(6,*) 'presib_root hum gt 1 ' , imsg_num , xmsg_sq
      endif
#endif
!
      RETURN
      END

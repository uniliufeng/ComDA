! 雪がある場合の大気中の輸送特性の変更（キャノピー空間）
! 99/06/0x 大泉
! 99/07/13 保坂征宏  引数変更。 
! 99/07/27 保坂征宏  XLAI 追加
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_MOD_AIR_SNW (
     I  NSTR             , NEND        , 
     I  DEP_SNOW         , NTYP         , MON         , 
     I  FAC_AG_INV_ms    , FAC_AC_BAR_INV_ms         , 
     I  XLAI_VEG_ms      , 
     I  D_ZERO_ms        , Z0_ms            , 
     I  CVR_VEG_ms       , CVR_VEG_INV_ms , 
     I  WTR_SNOW         , GLA_SNOW       , 
     O  Z0_msx           , 
     O  FAC_AG_INV_msx   , FAC_AC_BAR_INV_msx ,
     O  CVR_VEG_msx      , CVR_VEG_INV_msx    , 
     O  XLAI_VEG_msx     , D_ZERO_msx         , 
     O  HUM_SNSK         , FR_WBH_SNW         
     1    ) 
!
! やっていること
!   雪の深さの、キャノピーの高さに対する比が大きいほど
!     キャノピー空間の、雪が無いときに対する割合 FAC は小さくなり、
!     ZERO_CS が雪の深さ DEP_SNOW に近く、
!     粗度長が小さく、
!     積雪表面とキャノピー空間間の抵抗が小さく（長さが短くなるので）、
!     キャノピーとキャノピー空間間の抵抗定数が大きく
!           （やりとりできる葉が少なくなるので。
!             ここで葉の分布は一様と仮定した。）
!   なる。
!
!   ただし... 
!     ・雪が多くキャノピー高さ Z2 の 0.95 倍を越えるときは、便宜上、雪
!       の高さが Z2*0.95 しかないとする。 
!       これは Sellers ら(1986) のとっている方法である。
!       キャノピーが埋もれると、埋もれた箇所の扱いを例外処理する必要がある。
!       これをしたくないからであろう。
!
!     ・ファクタの評価方法は基本的にいい加減である。
!       雪があることによる効果を適当に、かつ無い場合（FAC=1）と連続的にな
!       るようにいれている。
!
!  要チェック
!    DEP_SNOW から作られる作業変数 SDEP をとっておかなくてよいか?  
!
      USE SIBPRM , ONLY : 
     1  IDX   , 
     1  ISN   , 
     1  ICG
!
      USE SIBVEG , ONLY : 
     1           Z2_mo              ,     ! キャノピー上端高さ (ITYP,IMON)
     1           Z1_mo                    ! キャノピー下端高さ (ITYP,IMON)
!
      USE SIBCON , ONLY :
     1            CVR_VEG_INV_MAX_c    
!
!     use com_debug , only : ipoint_new 
     
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR 
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: MON
      INTEGER,INTENT(IN)   :: NTYP   (IDX)         ! 植生タイプ
!
      REAL(8),INTENT(IN)   :: DEP_SNOW(IDX,ISN+1)! 積雪深( m ) 
!
      REAL(8),INTENT(IN)   :: FAC_AG_INV_ms (IDX)! 雪面とキャ間空力抵抗特性
!                                                ! に関する定数 旧 RDC 
      REAL(8),INTENT(IN)   :: FAC_AC_BAR_INV_ms (IDX) 
!                                                ! 積雪上のキャとキャ空間間の
!                                                ! 抵抗に関係する定数 旧 RBC
!
      REAL(8),INTENT(IN)   :: XLAI_VEG_ms    (IDX,ICG)
!
      REAL(8),INTENT(IN)   :: CVR_VEG_ms     (IDX,ICG)
      REAL(8),INTENT(IN)   :: CVR_VEG_INV_ms (IDX,ICG)
!
      REAL(8),INTENT(IN)   :: WTR_SNOW       (IDX)
      REAL(8),INTENT(IN)   :: GLA_SNOW       (IDX)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: FAC_AG_INV_msx (IDX)! 雪面とキャ間空力抵抗特性
!                                                ! に関する定数 旧 RDC 
      REAL(8),INTENT(OUT)  :: FAC_AC_BAR_INV_msx (IDX) 
!                                                ! 積雪上のキャとキャ空間間の
!                                                ! 抵抗に関係する定数 旧 RBC
!
      REAL(8),INTENT(OUT)  :: XLAI_VEG_msx (IDX,ICG)
      REAL(8),INTENT(OUT)  :: HUM_SNSK     (IDX)
      REAL(8),INTENT(OUT)  :: FR_WBH_SNW   (IDX)
!
      REAL(8),INTENT(OUT)  :: CVR_VEG_msx    (IDX,ICG)
      REAL(8),INTENT(OUT)  :: CVR_VEG_INV_msx(IDX,ICG)
!

!
! INOUT
!
      REAL(8),INTENT(IN)   :: D_ZERO_ms  (IDX)    ! 雪面上の零面変位 (m)
      REAL(8),INTENT(INOUT):: D_ZERO_msx (IDX)    ! 雪面上の零面変位 (m)
      REAL(8),INTENT(IN)   :: Z0_ms      (IDX)    ! 雪面上の粗度長   (m)
      REAL(8),INTENT(OUT)  :: Z0_msx     (IDX)    ! 雪面上の粗度長   (m)
!
! WORK 
!
      INTEGER   :: I
      REAL(8)   :: XZ
      REAL(8)   :: SDEP
      REAL(8)   :: FAC
      REAL(8),PARAMETER :: ZERO = 0.D0  
      REAL(8)   :: H2O_SNOW
! 
      REAL(8),PARAMETER :: H2O_MIN = 1.D-10 
!
      LOGICAl,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_MOD_AIR_SNOW 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 100 I = NSTR, NEND
#ifdef CORRECT
        XZ     = Z2_mo ( NTYP(I) , MON ) - Z1_mo ( NTYP(I) , MON ) 
        SDEP   = MIN ( DEP_SNOW(I,5)- Z1_mo ( NTYP(I) , MON )
     1                 , XZ * 0.95D0 )  
#else
        XZ     = Z2_mo  ( NTYP(I) , MON )
        SDEP   = MIN ( DEP_SNOW(I,5)  , XZ * 0.95D0 )  
#endif
        SDEP   = MAX ( SDEP , ZERO )  
        FAC    = 1.0D0 - SDEP / XZ                ! 雪あり 0.05 - 雪無し 1
!
        D_ZERO_msx (I)        = SDEP   + FAC * D_ZERO_ms(I)
        Z0_msx     (I)        = Z0_ms (I)  * FAC
!
!#ifdef CORRECT
        FAC_AG_INV_msx (I)    = FAC_AG_INV_ms (I)    / FAC
        FAC_AC_BAR_INV_msx(I) = FAC_AC_BAR_INV_ms(I) * FAC 
!#else
!        FAC_AG_INV_msx (I)    = FAC_AG_INV_ms (I)    * FAC
!        FAC_AC_BAR_INV_msx(I) = FAC_AC_BAR_INV_ms(I) / FAC 
!#endif
!
!  -----------------------------------------------
!
        XZ     = Z2_mo ( NTYP(I) , MON ) - Z1_mo ( NTYP(I) , MON ) 
        SDEP   = MIN ( DEP_SNOW(I,5)- Z1_mo ( NTYP(I) , MON )
     1                 , XZ * 0.95D0 )  
        SDEP   = MAX ( SDEP , ZERO )  
        FAC    = 1.D0 - SDEP /  XZ 
!
!
        XLAI_VEG_msx (I,1)    = XLAI_VEG_ms (I,1) * FAC
        XLAI_VEG_msx (I,2)    = 0. 
!
!  -----------------------------------------------
!
        HUM_SNSK      (I)     = 1.  
!
        CVR_VEG_msx (I,1)      = CVR_VEG_ms (I,1)
        CVR_VEG_msx (I,2)      = 0. 
        CVR_VEG_INV_msx (I,1)  = CVR_VEG_INV_ms (I,1)
        CVR_VEG_INV_msx (I,2)  = CVR_VEG_INV_MAX_c   
!
        H2O_SNOW      = MAX ( WTR_SNOW(I)+GLA_SNOW(I) , H2O_MIN )
        FR_WBH_SNW(I) = WTR_SNOW(I) / H2O_SNOW  
  100 CONTINUE
!
      RETURN
!
      END


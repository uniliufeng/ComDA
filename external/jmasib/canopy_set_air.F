! 蒸発・昇華関連計算   氷と水の割合
! 99/09/28 保坂征宏 canopy_get_ea から分割作成。
! 99/10/01 保坂征宏 裸地上での蒸発・昇華率の変更。水・氷そのままの割合にする
!
      SUBROUTINE CANOPY_SET_AIR  (
     I    NSTR         , NEND         , INSFLAG     , 
     I    FR_H2O_CNP   , FR_WBH_CNP   , FR_H2O_GRS  , FR_WBH_GRS  , 
     I    FR_WBH_BS    ,    
     I    FLAG_DRY_CNP , FLAG_DRY_SKIN, 
     O    XFR_WTR_CNP  , XFR_GLA_CNP  , XFR_WTR_GRS , XFR_GLA_GRS , 
     O    XFR_WTR_BS   , XFR_GLA_BS   )
!
! 注意
!     以前はインチキながらも、DEG_SAT_DTG(I) 
!     これは本当はインチキ（本当は零なのは温度変化量）だが、
!     やむなくこういう対処をする。
!
      USE SIBPRM , ONLY :
     1   IDX , 
     1   IDP , 
     1   ICG , 
     1   ITYP 
!
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR 
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG 
!
      REAL(8),INTENT(IN)   :: FR_H2O_CNP (IDX)   ! キャノピーの濡れ領域の割合
      REAL(8),INTENT(IN)   :: FR_H2O_GRS (IDX)   ! 下草の濡れ領域の割合
      REAL(8),INTENT(IN)   :: FR_WBH_CNP (IDX)   ! 
      REAL(8),INTENT(IN)   :: FR_WBH_GRS (IDX)   ! 
      REAL(8),INTENT(IN)   :: FR_WBH_BS  (IDX)   ! 
      REAL(8),INTENT(IN)   :: FLAG_DRY_SKIN  (IDX) 
      REAL(8),INTENT(IN)   :: FLAG_DRY_CNP   (IDX) 
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)  :: XFR_WTR_CNP (IDX)
      REAL(8),INTENT(OUT)  :: XFR_GLA_CNP (IDX)
      REAL(8),INTENT(OUT)  :: XFR_WTR_GRS (IDX)
      REAL(8),INTENT(OUT)  :: XFR_GLA_GRS (IDX)
      REAL(8),INTENT(OUT)  :: XFR_WTR_BS  (IDX)
      REAL(8),INTENT(OUT)  :: XFR_GLA_BS  (IDX)
!
!  WORK 
!
      INTEGER  ::  I 
      REAL(8)  ::  XFR_TMP 
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_SET_AIR 00/02 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  --------------------------------------
!  > 蒸発散計算に使う水・氷被覆率の設定 <
!  --------------------------------------
!               XFR_WTR_CNP 等
!               蒸散あり FLAG_DRY_CNP = 1 
!
! キャノピー
!
      DO 8000 I=NSTR,NEND
        XFR_TMP =  (1-FLAG_DRY_CNP(I)) 
     1              + FLAG_DRY_CNP(I)     * FR_H2O_CNP(I) 
        XFR_WTR_CNP(I) = XFR_TMP * FR_WBH_CNP(I)
        XFR_GLA_CNP(I) = XFR_TMP * (1-FR_WBH_CNP(I)) 
!
 8000 CONTINUE  
!
!  下草と裸地            
!        要 check : CVR_GRS や 1-CVR_GRS をかける必要の有無
!        かける必要なしと思う。
! 
      IF ( INSFLAG.EQ.0) THEN
      DO 8100 I=NSTR,NEND
        XFR_TMP =  (1-FLAG_DRY_SKIN(I)) 
     1              + FLAG_DRY_SKIN(I)     * FR_H2O_GRS(I) 
        XFR_WTR_GRS(I) = XFR_TMP * FR_WBH_GRS(I)
        XFR_GLA_GRS(I) = XFR_TMP * (1-FR_WBH_GRS(I)) 
        XFR_WTR_BS (I) = FR_WBH_BS(I)        
        XFR_GLA_BS (I) = 1-FR_WBH_BS(I)         
 8100 CONTINUE  
      ENDIF
!
!  雪上。
!
      IF ( INSFLAG .EQ. 1 ) THEN
        DO 8200 I=NSTR,NEND 
          XFR_WTR_BS (I) = FR_WBH_BS(I)        
          XFR_GLA_BS (I) = 1 - FR_WBH_BS(I)         
 8200   CONTINUE
      ENDIF
!
      RETURN
      END
!====================================================================
      SUBROUTINE CANOPY_AIR_ALL (
     I    NSTR       , NEND       , INSFLAG    , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , FR_H2O_GRS , 
     I    E_CNP_SAT  , E_REF      , E_SKIN_SAT   ,
     I    HUM_REL_SKIN   , 
     M    FLAG_DRY_CNP , FLAG_DRY_SKIN     , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
!
      USE SIBPRM , ONLY :
     1   IDX , 
     1   IDP , 
     1   ICG , 
     1   ITYP 
!
!
      IMPLICIT NONE
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT):: E_AIR (IDX)   ! キャノピー空間の水蒸気圧
      REAL(8),INTENT(OUT)  :: R_ACX_BAR   (IDX)  
!                                             ! キャノピー気孔抵抗+境界層抵抗
      REAL(8),INTENT(OUT)  :: R_ACX_BAR_INV (IDX) 
      REAL(8),INTENT(OUT)  :: R_AGX_BAR   (IDX)   ! 下草気孔抵抗+境界層抵抗
      REAL(8),INTENT(OUT)  :: R_AGX_BAR_INV (IDX)
      REAL(8),INTENT(OUT)  :: R_AB_INV (IDX)
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NSTR
      INTEGER,INTENT(IN)   :: NEND
      INTEGER,INTENT(IN)   :: INSFLAG
      REAL(8),INTENT(INOUT):: FLAG_DRY_CNP  (IDX) 
      REAL(8),INTENT(INOUT):: FLAG_DRY_SKIN (IDX)   
!
      REAL(8),INTENT(IN)   :: HUM_REL_SKIN (IDX)
!
      REAL(8),INTENT(IN)   :: FR_H2O_CNP(IDX)      ! キャの濡れた部分の割合
      REAL(8),INTENT(IN)   :: FR_H2O_GRS(IDX)      ! 下草の濡れた部分の割合
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
      REAL(8),INTENT(IN)   :: R_AIR_INV (IDX)      ! キャ空間と大気最下層間抵抗
!
      REAL(8),INTENT(IN)   :: E_CNP_SAT   (IDX)      ! TC での飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_SKIN_SAT  (IDX)      ! TG での飽和水蒸気圧
      REAL(8),INTENT(IN)   :: E_REF       (IDX)      ! 参照レベルでの水蒸気圧
!
! INOUT
!
      INTEGER :: ICHECK    (IDX)
!
!  WORK 
!
      INTEGER  ::  I 
      INTEGER  ::  IDEWCO  
!
      REAL(8)  ::  EGF     (IDX)   
      REAL(8)  ::  ECF     (IDX)   
!
      REAL(8)  ::  DONE
      REAL(8)  ::  DEWC 
      REAL(8)  ::  DEWG 
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  DONE / 1.D0 /
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_AIR_ALL 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO I=NSTR,NEND 
        ICHECK(I) = 1 
      ENDDO
!
! ------------------
! > E_AIR の計算 <
! ------------------
!
!     FLAG_DRY_CNP=1,FLAG_DRY_SKIN=1 
!       -> FLAG_DRY_CNP=1,FLAG_DRY_SKIN=0 
!          -> FLAG_DRY_CNP=0,FLAG_DRY_SKIN=1 
!             -> FLAG_DRY_CNP=0,FLAG_DRY_SKIN=0  
!     未飽和（蒸散あり） = 1 
!
      IDEWCO = 0 
!
 3000 CONTINUE                ! <<<----------------------------------  C
!
#ifdef SIB_DEBUG 
      IF ( IDEWCO .EQ. 4 ) THEN
        write(6,*) 'canopy_set_air error idewco=4 nstr nend icheck ',
     1             nstr,nend,
     1             icheck(nstr:nend)
        stop 999
      ENDIF
#endif
!
        IDEWCO = IDEWCO + 1
!
        CALL CANOPY_AIR ( 
     I    NSTR       , NEND       , INSFLAG    , ICHECK     , 
     I    FLAG_DRY_CNP  , FLAG_DRY_SKIN        , 
     I    CVR_VEG_msx  , CVR_VEG_INV_msx, 
     I    R_AC_BAR   , R_AC_BAR_INV , R_AG     , R_AG_INV    , 
     I    R_STM_BAR  , R_BB       , R_AIR_INV  , 
     I    FR_H2O_CNP , FR_H2O_GRS , 
     I    E_CNP_SAT  , E_REF      , E_SKIN_SAT   ,
     I    HUM_REL_SKIN     , 
     O    E_AIR      , 
     O    R_ACX_BAR  , R_AGX_BAR  , R_ACX_BAR_INV  , R_AGX_BAR_INV  , 
     O    R_AB_INV )
!
      DO 3200 I=NSTR,NEND
        IF ( ICHECK(I).EQ.1 ) THEN    
          ECF (I) = SIGN ( DONE , E_CNP_SAT (I)-E_AIR(I) ) ! 過飽和-1 未飽和1 
          EGF (I) = SIGN ( DONE , E_SKIN_SAT(I)-E_AIR(I) )
          DEWC    = FLAG_DRY_CNP (I)  * 2.0D0 - 1.0D0    ! 未飽和 1 
          DEWG    = FLAG_DRY_SKIN(I)  * 2.0D0 - 1.0D0    
          ECF (I) = ECF(I) * DEWC
          EGF (I) = EGF(I) * DEWG
        END IF
 3200 CONTINUE
!
#ifdef CORRECT
      DO 3300 I=NSTR,NEND 
        IF ( ( ECF(I).GT.0.0D0 .AND. EGF(I).GT.0.0D0 )
     1        .OR. ( IDEWCO.EQ.4 )                         ) THEN
          ICHECK(I)=0
        END IF
 3300 CONTINUE
#else
      DO 3300 I=NSTR,NEND 
        IF ( ECF(I).GT.0.0D0 .AND. EGF(I).GT.0.0D0 )  THEN 
          ICHECK(I)=0
        END IF
 3300 CONTINUE
#endif
!
      DO 3400 I=NSTR,NEND            ! 雪では、FLAG_DRY_SKIN = 0. のみでよい
        IF ( ICHECK(I).EQ.1 ) THEN 
          IF     (     FLAG_DRY_CNP (I).GT.0.5D0 
     1           .AND. FLAG_DRY_SKIN(I).GT.0.5D0 ) THEN     
            FLAG_DRY_CNP(I)  = 1.0D0
            FLAG_DRY_SKIN(I) = 0.0D0
          ELSEIF (     FLAG_DRY_CNP(I).GT.0.5D0
     1           .AND. FLAG_DRY_SKIN(I).LT.0.5D0 ) THEN     
            FLAG_DRY_CNP(I)  = 0.0D0
            FLAG_DRY_SKIN(I) = 1.0D0
          ELSEIF (     FLAG_DRY_CNP(I) .LT.0.5D0
     1           .AND. FLAG_DRY_SKIN(I).GT.0.5D0 ) THEN     
            FLAG_DRY_CNP(I)  = 0.0D0
            FLAG_DRY_SKIN(I) = 0.0D0
          ELSEIF (     FLAG_DRY_CNP(I) .LT.0.5D0
     1           .AND. FLAG_DRY_SKIN(I).LT.0.5D0 ) THEN     
            FLAG_DRY_CNP (I) = 1.0D0
            FLAG_DRY_SKIN(I) = 1.0D0
          ENDIF
        ENDIF
 3400 CONTINUE
!
#include "noparallel_fixedform"
      DO  I=NSTR,NEND            ! 雪では、FLAG_DRY_SKIN = 0. のみでよい
        IF ( ICHECK(I).EQ.1 ) THEN 
          GOTO 3000                   ! --------------- >>> C
        ENDIF
      ENDDO

!
      RETURN
      END







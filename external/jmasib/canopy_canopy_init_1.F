! キャノピー計算で必要な、キャノピー空間量に関する情報、初期値を作成する
! 99/07/25 保坂征宏 新規作成
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_CANOPY_INIT_1 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , MON         , 
     I   TMP_REF     , U_REF       , V_REF       , Z_REF      ,   
     I   Z0_msx      , D_ZERO_msx  , CVR_VEG_msx ,
     M   R_STM_BAR_PAR , 
     O   TMP_AIR     , R_AIR_INV   , R_AIR_CONST_INV  , 
     O   U2          , FAC_CORTMP     , 
     O   A_COEF      , 
     O   R_STM_BAR   , FLAG_DRY_CNP   , FLAG_DRY_SKIN           )
!
! 注意。
!    Z0 については雪の存在による変更処理を施した値を用いるので、
!    canopy_mod_air_snw の後に計算する必要がある。
!    そのため、雪あり部分と雪無し部分それぞれで呼ぶ必要がある。
!
      USE SIBPRM , ONLY :
     1   IDX  , 
     1   ICG
!
      USE SIBVEG , ONLY : 
     1           Z2_mo                   , ! キャノピー上端高さ (ITYP,IMON)
     1           Z1_mo                     ! キャノピー上端高さ (ITYP,IMON)
!
      USE SIBCON , ONLY : 
     1   GRAV_c                     , 
     1   VKRMN_c                     , 
     1   VKRMN_INV_c
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: NTYP(IDX)
      INTEGER,INTENT(IN)    :: MON
!
      INTEGER,INTENT(IN)    :: INSFLAG
!
      REAL(8),INTENT(IN)    :: TMP_REF     (IDX)
      REAL(8),INTENT(IN)    :: U_REF       (IDX)
      REAL(8),INTENT(IN)    :: V_REF       (IDX)
      REAL(8),INTENT(IN)    :: Z_REF       (IDX)
      REAL(8),INTENT(IN)    :: FAC_CORTMP  (IDX)
      REAL(8),INTENT(IN)    :: Z0_msx      (IDX)
      REAL(8),INTENT(IN)    :: D_ZERO_msx  (IDX)
      REAL(8),INTENT(IN)    :: CVR_VEG_msx      (IDX,ICG)  !
      REAL(8),INTENT(INOUT) :: R_STM_BAR_PAR    (IDX,ICG)  !
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: U2        (IDX)   ! 
      REAL(8),INTENT(OUT)   :: TMP_AIR   (IDX)   ! 
      REAL(8),INTENT(OUT)   :: R_AIR_INV (IDX)   ! 
      REAL(8),INTENT(OUT)   :: R_AIR_CONST_INV(IDX)   ! 
      REAL(8),INTENT(OUT)   :: A_COEF         (IDX)   ! 
      REAL(8),INTENT(OUT)   :: R_STM_BAR      (IDX,ICG)  !
      REAL(8),INTENT(OUT)   :: FLAG_DRY_SKIN  (IDX)  !
      REAL(8),INTENT(OUT)   :: FLAG_DRY_CNP   (IDX)  !
!
! WORK
!
      INTEGER    :: I
      INTEGER    :: IVEG
!  
      REAL(8)    :: UMIN_TUNE 
      REAL(8)    :: Z0MIN_TUNE 
      REAL(8)    :: ABS_U_REF     ! 
      REAL(8)    :: Z0            ! 上下限処理した、モデルで使う Z0 、かな??
!                                 ! 実際にはここで使うだけ。
!
      REAL(8)    :: Z_REFER       ! Z_REF - D_ZERO 
      REAL(8)    :: C_M_CONST     ! Cm の一定部分、つまり a^2
      REAL(8)    :: C_M_FAC       ! Cm の RiB 依存性部分。わからないので 1
      REAL(8)    :: A_TMP         ! 1-2-95 の分母
!
      LOGICAL,SAVE ::  LFIRST
!
      DATA UMIN_TUNE   / 2.D0   /
      DATA Z0MIN_TUNE  / 0.01D0 /
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'CANOPY_CANOPY_INIT 99/07 POCHKA' 
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=NSTR,NEND 
!
        ABS_U_REF  = MAX ( SQRT(U_REF(I)**2+V_REF(I)**2) , UMIN_TUNE )
!
        Z0         = MAX ( Z0_msx(I) , Z0MIN_TUNE )
        Z_REFER    = MAX ( 15.0D0*Z0 , Z_REF (I) - D_ZERO_msx(I) )
!
        A_TMP      = 1 / LOG ( Z_REFER / Z0 )   
        A_COEF (I) =  VKRMN_c * A_TMP 
!
        U2  (I)  
     1     =  ABS_U_REF 
     1          * LOG ( (Z2_mo(NTYP(I),MON)-D_ZERO_msx(I))/Z0_msx(I) )
     1          * A_TMP 
!
        C_M_CONST = A_COEF(I)**2 
        C_M_FAC   = 1. 
!
        R_AIR_CONST_INV (I)  =  C_M_CONST * ABS_U_REF  
        R_AIR_INV       (I)  =  R_AIR_CONST_INV (I) * C_M_FAC
!
        TMP_AIR(I) = TMP_REF (I)  * FAC_CORTMP(I) 
 1000 CONTINUE
!
      DO IVEG=1,2-INSFLAG
        DO I=NSTR,NEND
          R_STM_BAR_PAR(I,IVEG) = R_STM_BAR_PAR(I,IVEG)
     1                                   /CVR_VEG_msx(I,IVEG)     
          R_STM_BAR    (I,IVEG) = R_STM_BAR_PAR(I,IVEG)     
!!          if (i.eq.4 .and. iveg.eq.1 ) 
!!     1           write(6,*) 'ISTEP canopy_1 ' , i , 
!!     1                      R_STM_BAR_PAR(I,IVEG)
!!     1                      , CVR_VEG_msx(I,IVEG)     
        ENDDO
      ENDDO
!
      DO I=NSTR,NEND 
        FLAG_DRY_SKIN(I) = 1.  
        FLAG_DRY_CNP(I)  = 1. 
      ENDDO 
!
      RETURN
      END

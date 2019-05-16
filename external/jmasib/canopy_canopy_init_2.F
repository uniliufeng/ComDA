! キャノピー計算で必要な、キャノピー空間量に関する情報、初期値を作成する
! 99/07/25 保坂征宏 新規作成
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_CANOPY_INIT_2 (
     I   NSTR        , NEND        , INSFLAG     , 
     I   NTYP        , 
     I   TMP_CNP     , TMP_GRSK    , 
     I   U_REF       , V_REF       , Z_REF      ,   
     I   Z0_msx      , D_ZERO_msx  ,
     I   TMP_REF     , FAC_CORTMP , 
     O   RIB_CONST   ,
     O   FUNC_TMP    )   

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
     1           F_RST_T_OPT_co          , ! 最適温度  (ITYP,ICG) 
     1           F_RST_T_MIN_co          , ! 最低限界温度 (ITYP,ICG) 
     1           F_RST_T_MAX_co            ! 最高限界温度 (ITYP,ICG)
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
!
      INTEGER,INTENT(IN)    :: INSFLAG
!
      REAL(8),INTENT(IN)    :: TMP_CNP     (IDX)
      REAL(8),INTENT(IN)    :: TMP_GRSK    (IDX)
      REAL(8),INTENT(IN)    :: U_REF       (IDX)
      REAL(8),INTENT(IN)    :: V_REF       (IDX)
      REAL(8),INTENT(IN)    :: Z_REF       (IDX)
      REAL(8),INTENT(IN)    :: Z0_msx      (IDX)
      REAL(8),INTENT(IN)    :: D_ZERO_msx  (IDX)
      REAL(8),INTENT(IN)    :: FAC_CORTMP  (IDX)   ! 
      REAL(8),INTENT(IN)    :: TMP_REF     (IDX)   ! 
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: RIB_CONST   (IDX)   ! 
      REAL(8),INTENT(OUT)   :: FUNC_TMP    (IDX,ICG) ! 気孔抵抗の温度依存性
!
! WORK
!
      INTEGER    :: I
!  
      REAL(8)    :: UMIN_TUNE 
      REAL(8)    :: Z0MIN_TUNE 
      REAL(8)    :: ABS_U_REF     ! 
      REAL(8)    :: Z0            ! 上下限処理した、モデルで使う Z0 、かな??
!                                 ! 実際にはここで使うだけ。
!
      REAL(8)    :: Z_REFER       ! Z_REF - D_ZERO 
      REAL(8)    :: T1
      REAL(8)    :: T2
      REAL(8)    :: T3
!
      REAL(8),PARAMETER  :: ZERO = 0.D0
!
      LOGICAL,SAVE ::  LFIRST
!
      DATA UMIN_TUNE   / 2.D0   /
      DATA Z0MIN_TUNE  / 0.01D0 /
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'CANOPY_CANOPY_INIT_2 00/01 POCHKA' 
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
!  RIB_CONST 
!       佐藤文書 1-2-96 の (Ta-κT1) 以外。
!       SiB9603 では TMP_AIR ではなく TMP_REF で割っていた。これを修正した。
!
!     010725 , TMP_AIR ではなく κT_1 であるところの
!
        RIB_CONST(I) =  - GRAV_c * Z_REFER 
     1                  / ( TMP_REF(I) * FAC_CORTMP(I) * ABS_U_REF**2 )        
!
!  FUNC_TMP 
!
        T1  =   MAX( TMP_CNP(I) - F_RST_T_MIN_co(NTYP(I),1) , ZERO )
     1           / ( F_RST_T_OPT_co(NTYP(I),1) 
     1                   - F_RST_T_MIN_co(NTYP(I),1) )
!
        T2  =   MAX( F_RST_T_MAX_co(NTYP(I),1) - TMP_CNP(I) , ZERO )
     1           / ( F_RST_T_MAX_co(NTYP(I),1) 
     1                  - F_RST_T_OPT_co(NTYP(I),1) )
!
        T3  =   ( F_RST_T_MAX_co(NTYP(I),1)  
     1               - F_RST_T_OPT_co(NTYP(I),1) )
     1           / ( F_RST_T_OPT_co(NTYP(I),1) 
     1                  - F_RST_T_MIN_co(NTYP(I),1) )
        FUNC_TMP(I,1) = T1 * T2**T3 
 1000 CONTINUE
!
      IF ( INSFLAG .EQ. 0 ) THEN
      DO 2000 I=NSTR,NEND 
!
        T1  =   MAX ( TMP_GRSK(I) - F_RST_T_MIN_co(NTYP(I),2) , ZERO )
     1           / ( F_RST_T_OPT_co(NTYP(I),2) 
     1                   - F_RST_T_MIN_co(NTYP(I),2) )
        T2  =   MAX( F_RST_T_MAX_co(NTYP(I),2) - TMP_CNP(I) , ZERO ) 
     1           / ( F_RST_T_MAX_co(NTYP(I),2) 
     1                  - F_RST_T_OPT_co(NTYP(I),2) )
        T3  =   ( F_RST_T_MAX_co(NTYP(I),2)  
     1               - F_RST_T_OPT_co(NTYP(I),2) )
     1           / ( F_RST_T_OPT_co(NTYP(I),2) 
     1                  - F_RST_T_MIN_co(NTYP(I),2) )
        FUNC_TMP(I,2) = T1 * T2**T3
 2000 CONTINUE  
      ENDIF
!
      RETURN
      END

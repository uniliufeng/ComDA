! 放射フラックス （各種温度依存性の微係数）計算
!
! 99/06/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 99/08/30 保坂征宏 DRAD_DX の符号変更。収束が正になるようにする。
!
!
      SUBROUTINE CANOPY_RAD (
     I  NSTR       , NEND        , 
     I  TMP_CNP    , TMP_SKIN    , 
     I  THERMK     , 
     I  CVR_VEG_msx, 
!    I  FLAG_TMPMOD_CNP  , FLAG_TMPMOD_SKIN  , 
     O  DRAD_DX    )
!
! 注意。
!     RC ... キャ・地面・雪面射出長波放射ネットで、キャが出す放射
!     RG ... キャ・地面・雪面射出長波放射ネットで、地面が出す放射
!     RS ... キャ・地面・雪面射出長波放射ネットで、雪面が出す放射
! 
!     RC    には CVR_VEG(I,1) がかかっている！
!     RG,TS には CVR_SNW はかかっていない。
!
!     雪のあるところ、無いところ共通。
!
!     DRAD_DX は、受け取る向き（収束）が正!
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ICG    , 
     1    ICGR
!
      USE SIBCON , ONLY : 
     1    STB4_C    
!
      IMPLICIT NONE 
!
! INPUT
!
      REAL(8),INTENT(IN)    :: TMP_CNP   (IDX)      ! キャノピー温度
      REAL(8),INTENT(IN)    :: TMP_SKIN  (IDX)      ! 地面・雪面上端温度
      REAL(8),INTENT(IN)    :: THERMK    (IDX)      ! キャ長波透過率
      REAL(8),INTENT(IN)    :: CVR_VEG_msx(IDX,ICG)
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
!
!     REAL(8),INTENT(IN)    :: FLAG_TMPMOD_CNP  (IDX)
!     REAL(8),INTENT(IN)    :: FLAG_TMPMOD_SKIN (IDX)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)   :: DRAD_DX     (IDX,ICGR,ICG) 
!
! WORK
!    
      REAL(8)               :: S4_TMP_C3 (IDX)
      REAL(8)               :: S4_TMP_G3 (IDX)
      REAL(8)               :: FAC       (IDX)     ! キャ吸収率 * Cv 
!
      INTEGER               :: I
!
      LOGICAL,SAVE :: LFIRST
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_IR 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=NSTR,NEND
!
! キャノピーの長波放射吸収率
!
        FAC (I)   = CVR_VEG_msx(I,1) * (1.0D0-THERMK (I))
!
! 作業変数
!
        S4_TMP_C3 (I) = STB4_c * TMP_CNP (I)**3
        S4_TMP_G3 (I) = STB4_c * TMP_SKIN (I)**3 
!
!     TC 変化によるフラックス変化 
!
        DRAD_DX(I,1,1) = 0                       ! TC 変化で C-C の変化
        DRAD_DX(I,2,1) =   FAC(I) * S4_TMP_C3(I) ! TC 変化で C-G の変化
        DRAD_DX(I,3,1) = - FAC(I) * S4_TMP_C3(I) ! TC 変化で A-C の変化
!
!     TG 変化によるフラックス変化 
!
        DRAD_DX(I,1,2) = - FAC(I) * S4_TMP_G3(I) ! TG 変化で C-G の変化
        DRAD_DX(I,2,2) = 0. 
        DRAD_DX(I,3,2) = - (1-FAC(I))* S4_TMP_G3(I)! TG 変化で A-G の変化
!
 1000 CONTINUE 
!
#ifdef KILLRAD
      call reset ( drad_dx_tmp , idx * icg * icg ) ! zzz
#endif
!
      RETURN
      END




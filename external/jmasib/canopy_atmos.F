! キャノピー計算に必要な、大気参照レベルの情報を作る。
! 99/07/25 保坂征宏 新規作成
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_ATMOS (
     I   NMAX          , 
     I   P_SFC         , TMP_REF       ,  Q_REF     , 
     O   RHO_AIR       , RHO_CP   , 
     O   E_REF         , 
     O   FAC_EVP       , FAC_EVP_2 )
!
! 注意。
!   雪有り無しに関係ない量の計算
!
!   fac_evp について。
!     後で使うときに 水蒸気圧 E_*  とのかけ算の形で使う。
!     E_* の単位も、P_SFC の単位も HPa = である。
!
! 水蒸気圧 e ( 単位 HPa ) の評価には、
!    e ～ pq/ε という近似式を用いている
!
      USE SIBPRM , ONLY :
     1   IDX
!
      USE SIBCON , ONLY : 
     1   R_DRY_INV_c                ,  ! 乾燥空気気体定数
     1   HL_WTR_c                   ,  ! 
     1   CP_AIR_c                   ,  ! 
     1   EPSFAC_c                      ! 0.622(無次元)=分子量比 H2O/乾燥
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX
!
      REAL(8),INTENT(IN)    :: P_SFC    (IDX)
      REAL(8),INTENT(IN)    :: TMP_REF  (IDX)
      REAL(8),INTENT(IN)    :: Q_REF    (IDX)
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: RHO_AIR   (IDX)   ! 
      REAL(8),INTENT(OUT)   :: RHO_CP    (IDX)   ! 
      REAL(8),INTENT(OUT)   :: E_REF   (IDX)     ! 
      REAL(8),INTENT(OUT)   :: FAC_EVP   (IDX)   ! ρε/P_s
      REAL(8),INTENT(OUT)   :: FAC_EVP_2 (IDX)   ! P_s/ε   e ～ FAC_EVP_2 q
!
      INTEGER    :: I
!  
      LOGICAL,SAVE ::  LFIRST
!
      DATA LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'CANOPY_ATMOS 99/07 POCHKA' 
        LFIRST = .FALSE. 
      ENDIF
!
      DO 1000 I=1,NMAX
        RHO_AIR  (I) = P_SFC(I) * 100.0D0 * R_DRY_INV_c / TMP_REF(I)  
        RHO_CP   (I) = RHO_AIR(I) * CP_AIR_c
        FAC_EVP  (I) = RHO_AIR(I) * EPSFAC_c / P_SFC(I)
        FAC_EVP_2(I) = P_SFC(I) / EPSFAC_c 
        E_REF    (I) = FAC_EVP_2(I) * Q_REF(I)  
 1000 CONTINUE
!
      RETURN
      END

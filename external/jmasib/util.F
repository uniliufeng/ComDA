! SiB 用ユーティリティーサブルーチンパッケージ
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE RESET ( WORK , IMAX )
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: IMAX      
!
      REAL(8),INTENT(INOUT) :: WORK(IMAX)
!
      INTEGER               :: I
!
      DO 100 I=1,IMAX
        WORK(I) = 0.D0
 100  CONTINUE
!
      RETURN
      END
!===========================================================
      SUBROUTINE COPY ( WORK_IN , WORK_OUT , IMAX )
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: IMAX      
      REAL(8),INTENT(IN)    :: WORK_IN  (IMAX)
      REAL(8),INTENT(OUT)   :: WORK_OUT (IMAX)
!
      INTEGER               :: I
!
      DO 100 I=1,IMAX
        WORK_OUT(I) = WORK_IN(I) 
 100  CONTINUE
!
      RETURN
      END
!=============================================================
      SUBROUTINE IRESET ( IWORK , IMAX )
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: IMAX      
!
      INTEGER,INTENT(INOUT) :: IWORK(IMAX)
!
      INTEGER               :: I
!
      DO 100 I=1,IMAX
        IWORK(I) = 0
 100  CONTINUE
!
      RETURN
      END
!==========================================================
      SUBROUTINE CON_SET ( WORK , IMAX , CON )
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: IMAX      
      REAL(8),INTENT(IN)    :: CON 
!
      REAL(8),INTENT(INOUT) :: WORK(IMAX)
!
      INTEGER               :: I
!
      DO 100 I=1,IMAX
        WORK(I) = CON
 100  CONTINUE
!
      RETURN
      END
!==========================================================
      SUBROUTINE SET_UNDEF ( WORK , IMAX )
!
      IMPLICIT NONE 
!
      INTEGER,INTENT(IN)    :: IMAX      
!
      REAL(8),INTENT(INOUT) :: WORK(IMAX)
!
      INTEGER               :: I
!
      DO 100 I=1,IMAX
        WORK(I) = -9.99D33
 100  CONTINUE
!
      RETURN
      END

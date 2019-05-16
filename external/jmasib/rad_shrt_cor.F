! 地表面（SiB）が受け取る短波放射量の、時間補正
! 99/05/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE RAD_SHRT_COR (
     I  NMAX    , NMAX_01 , 
     I  ZMEAN   , ZTEMP   , 
     M  RAD_SHRT_NOS      ,  RAD_SHRT_SNW       )  
!
! 地表面（SiB）が受け取る短波放射量の、時間補正
!   短波放射は一時間に一度しか計算しないが、
!   SiB は毎ステップ計算するために必要になる。
!
! 注意点
!   ・RAD_CG はこの時点では、短波寄与分のみ。あとで長波分が足しこまれる。
!     また、雪面については、この時点ではゼロでない値が入っている。
!     （いずれ RAD_SNW や RAD_DS1 に振り分けられ、ゼロになる）
!
!   ・従来の SiB では同じバンド内の陸格子を取り出して順につめればよかったが、
!     大泉・保坂 SiB では雪無し域と雪あり域を別々にする必要あり。
!     （雪あり域の RAD_CG,RADN と雪無し域の RAD_CG,RADN とを別々に得るため。）
!
! 課題 （バグ??? ）
!   補正の仕方がエネルギー保存的でない?! 
!   つまり、
!     大気短波放射計算で一時間に地面が吸収したであろう放射エネルギー量と
!     地表面が実際に受け取った放射エネルギー量にずれがある。
!   要チェック
!   
!
      USE SIBPRM , ONLY : 
     1   IDX    , 
     1   ICG    
!
      IMPLICIT NONE 
!
! 入力
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_01(3)
!
      REAL(8),INTENT(IN)   :: ZMEAN  (IDX)
      REAL(8),INTENT(IN)   :: ZTEMP  (IDX)
!
! OUTPUT 
!
      REAL(8),INTENT(INOUT):: RAD_SHRT_NOS (IDX,ICG) ! 短波吸収量 BD,WV合計。
      REAL(8),INTENT(INOUT):: RAD_SHRT_SNW (IDX,ICG) ! 短波吸収量 BD,WV合計。
!
! WORK 
!
      INTEGER   :: N 
      REAL(8),PARAMETER  :: D01 = 0.01D0
!

      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD_SHRT_COR 99/07 POCHAKA '
        LFIRST = .FALSE.
      ENDIF 
!
!   --------------------
!   >>  陸面入射短波  <<
!   --------------------
!     これはもともと XRAD_CG が SiB 格子である。
!     雪あり時は下草はない。
!
      DO 1000 N=1,NMAX_01(3)
        RAD_SHRT_NOS(N,1  ) = RAD_SHRT_NOS (N,1) 
     1                * ZTEMP(N) / MAX(D01,ZMEAN(N))    ! 補正  ppp
        RAD_SHRT_NOS(N,2  ) = RAD_SHRT_NOS (N,2) 
     1                * ZTEMP(N) / MAX(D01,ZMEAN(N))    ! 補正  ppp
 1000 CONTINUE
!
      DO 2000 N=NMAX_01(2),NMAX
        RAD_SHRT_SNW(N,1  ) = RAD_SHRT_SNW (N,1) 
     1                * ZTEMP(N) / MAX(D01,ZMEAN(N))    ! 補正  ppp
        RAD_SHRT_SNW(N,2  ) = RAD_SHRT_SNW (N,2) 
     1                * ZTEMP(N) / MAX(D01,ZMEAN(N))    ! 補正
 2000 CONTINUE
!
      RETURN
      END

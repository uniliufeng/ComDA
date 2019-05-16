!  機能 : 大気下端での長波放射フラックスを放射格子から陸面格子に詰め替え。
!         長波放射計算ステップのみ
!  99/05/xx 保坂征宏
!  99/07/21 保坂征宏 変数名、サブルーチン名変更
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE RAD_SHRT_RGRD_2_SIBGRD (
     I IMSKLT_RX , IMSKLT_RY , NMAX      ,
     I RVISB     , RVISD     , RNIRB     , RNIRD     , 
     O RAD_SHRT_INTO_SIB   )
!
!  機能 : 大気下端での長波放射フラックスを放射格子から陸面格子に詰め替え。
!         長波放射計算ステップのみ
!
!  注意 : 従来との対応としては、DLWB を消去。これはもともと全面データだっ
!         たため。
!
      USE PRM , ONLY : 
     1    IRAD  , 
     1    JRAD  , 
     1    JLPHY  
!
      USE SIBPRM , ONLY :
     1    IDX    ,
     1    IVN    ,
     1    IBD   
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   ::  IMSKLT_RX ( IDX )      ! SiB -> 物理
      INTEGER,INTENT(IN)   ::  IMSKLT_RY ( IDX )
      INTEGER,INTENT(IN)   ::  NMAX 
!
      REAL(8),INTENT(IN)   ::  RVISB   ( IRAD , JRAD /JLPHY )   
      REAL(8),INTENT(IN)   ::  RVISD   ( IRAD , JRAD /JLPHY)
      REAL(8),INTENT(IN)   ::  RNIRB   ( IRAD , JRAD /JLPHY)
      REAL(8),INTENT(IN)   ::  RNIRD   ( IRAD , JRAD /JLPHY)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  ::  
     1   RAD_SHRT_INTO_SIB ( IDX, IVN, IBD ) ! 大気下端下向放射
!
! WORK
      INTEGER              :: N
      INTEGER              :: IR
      INTEGER              :: JR
!
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD_SHRT_RGRD_2_SIBGRD 99/07 POCHAKA '
        LFIRST = .FALSE.
      ENDIF 

!
      DO 1000 N=1,NMAX
        IR  = IMSKLT_RX(N)
        JR  = IMSKLT_RY(N)
        RAD_SHRT_INTO_SIB(N,1,1) = RVISB (IR,JR)               ! 可視・直達
        RAD_SHRT_INTO_SIB(N,1,2) = RVISD (IR,JR)               ! 可視・散乱
        RAD_SHRT_INTO_SIB(N,2,1) = RNIRB (IR,JR)               ! 近赤・直達
        RAD_SHRT_INTO_SIB(N,2,2) = RNIRD (IR,JR)               ! 近赤・散乱
 1000 CONTINUE
!
      RETURN
      END
! ===============================================================
      SUBROUTINE RAD_LONG_RGRD_2_SIBGRD (
     I IMSKLT_RX , IMSKLT_RY , NMAX      ,
     I DLWB      , 
     O RAD_LONG_INTO_SIB   )
!
!  機能 : 大気下端での長波放射フラックスを放射格子から陸面格子に詰め替え。
!         長波放射計算ステップのみ
!
!  注意 : 従来との対応としては、DLWB を消去。これはもともと全面データだっ
!         たため。
!
      USE PRM , ONLY : 
     1    IRAD  , 
     1    JRAD  ,
     1    JLPHY 
!
      USE SIBPRM , ONLY :
     1    IDX    
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   ::  IMSKLT_RX ( IDX )      ! SiB -> 物理
      INTEGER,INTENT(IN)   ::  IMSKLT_RY ( IDX )
      INTEGER,INTENT(IN)   ::  NMAX 
!
      REAL(8),INTENT(IN)   ::  DLWB   ( IRAD , JRAD/JLPHY )   
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  ::  
     1   RAD_LONG_INTO_SIB ( IDX ) ! 大気下端下向放射
!
! WORK
      INTEGER              :: N
      INTEGER              :: IR
      INTEGER              :: JR
!
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD_LONG_RGRD_2_SIBGRD 99/07 POCHAKA '
        LFIRST = .FALSE.
      ENDIF 
!
      DO 1000 N=1,NMAX
        IR  = IMSKLT_RX(N)
        JR  = IMSKLT_RY(N)
        RAD_LONG_INTO_SIB(N) = DLWB (IR,JR) 
 1000 CONTINUE
!
      RETURN
      END

      SUBROUTINE WATER_OUT ( 
     I     NSTR       , NEND       , 
     I     PR_TRAP    , PR_TRAP_MAX, 
     I     A_PR_COEF  , B_PR_COEF  , C_PR_COEF  , 
     O     PR_DROP     )
!
! 植生・土壌第一層で一旦トラップした PR_TRAP のうち、
! 可能貯水容量 PR_TRAP_MAX を越える分 PR_DROP を落下させる。
!
      USE SIBPRM , ONLY :
     1    IDX 
!
      IMPLICIT NONE
!
! INPUT 
      INTEGER,INTENT(IN)    :: NSTR 
      INTEGER,INTENT(IN)    :: NEND
      REAL(8),INTENT(IN)    :: PR_TRAP     (IDX)  ! 一旦捕捉される雨量
      REAL(8),INTENT(IN)    :: PR_TRAP_MAX (IDX)  ! 捕捉できる上限値
                                                  !    = 飽和 - 現有
      REAL(8),INTENT(IN)    :: A_PR_COEF   (IDX)   
      REAL(8),INTENT(IN)    :: B_PR_COEF 
      REAL(8),INTENT(IN)    :: C_PR_COEF   (IDX)   
!
! OUTPUT 
      REAL(8),INTENT(OUT)   :: PR_DROP     (IDX)  ! 落下する量
!
! WORK 
      INTEGER    :: I
!
      REAL(8)    :: ARG 
      REAL(8)    :: X_SAT
      REAL(8)    :: ARG_MIN     
      REAL(8)    :: DZERO
      REAL(8)    :: DONE
      REAL(8)    :: TRAP_MIN     ! PR_TRAP_MAX がこれ以下の場合は全て落下
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA ARG_MIN  / 1.D-12 /   ! B_PR_COEF より小さいこと。
      DATA TRAP_MIN / 1.D-12 /   
      DATA DZERO    / 0.D0 /
      DATA DONE     / 1.D0 /
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'LDRAIN 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      DO 100 I=NSTR,NEND
!
! 雨量が少なければ全て捕捉 
!      （ゼロ割りのエラー処理のため。既に飽和している場合は最後）
!
        IF ( PR_TRAP(I) .LT. 1.0D-20 ) THEN
          PR_DROP (I) = 0.0D0
!
! 雨量が少しはあれば
        ELSE
          ARG = ( PR_TRAP_MAX(I) / PR_TRAP(I) - C_PR_COEF(I) ) 
     1           / A_PR_COEF(I)
!
          IF ( ARG .GE. ARG_MIN ) THEN                  ! 雨量が多くない
            X_SAT = - 1.0D0 / B_PR_COEF * LOG ( ARG )
            X_SAT = MIN ( X_SAT , DONE  )
            X_SAT = MAX ( X_SAT , DZERO )
          ELSE                                          ! 雨量が多い
            X_SAT  = 1.0D0
          END IF
!
          PR_DROP(I) 
     1      =  PR_TRAP(I) * A_PR_COEF(I) / B_PR_COEF 
     1                    * ( 1.0D0 - DEXP ( - B_PR_COEF * X_SAT ) )  
     1         + ( C_PR_COEF(I) * PR_TRAP(I) - PR_TRAP_MAX(I) ) * X_SAT     
!
        END IF
!
! もし飽和していたら

        IF ( PR_TRAP_MAX (I) .LT. TRAP_MIN ) THEN
          PR_DROP (I) = PR_TRAP (I)
        ENDIF
  100 CONTINUE
!
      RETURN
      END

      SUBROUTINE RADSHRT_MOD (
     I    NMAX              , ILIST             , JL         , 
     I    CVR_SNW           , CVR_SNW_NEXT      , 
     I    INFO_SNOW         , INFO_SNOW_NEXT    , 
     M    RAD_SHRT_NOS_SAVE   , R_STM_BAR_PAR_NOS_SAVE  , 
     M    RAD_SHRT_SNW_SAVE   , R_STM_BAR_PAR_SNW_SAVE  ) 
!
!  雪の被覆率が変わったときで短波放射計算の直前以外のステップで、
!  短波放射関連の諸量を修正する。
!
!   ・短波フラックスについては、
!     (1) 未定義域に値ができたときはコピー
!     (2) 面積が変わった場合は平均値が変わらないように修正
!   ・気孔抵抗については
!     (1) 未定義域に値ができたときはコピー
!
!  という対策をとる。
!
      USE PRM  , ONLY :
     1    IJPHY , 
     1    JLPHY , 
     1    ISPT  
!
      USE SIBPRM , ONLY : 
     1     IDX   ,
     1     ICG 
!
! INPUT
!
      INTEGER,INTENT(IN)   ::  NMAX
      INTEGER,INTENT(IN)   ::  JL
      INTEGER,INTENT(IN)   ::  ILIST  (IDX)    ! SiB -> 全SiB 
!
      REAL(8),INTENT(IN)   ::  CVR_SNW       (IDX)
      REAL(8),INTENT(IN)   ::  CVR_SNW_NEXT  (IDX)
      INTEGER,INTENT(IN)   ::  INFO_SNOW     (IDX)
      INTEGER,INTENT(IN)   ::  INFO_SNOW_NEXT(IDX)
!
! INOUT
!
      REAL(8),INTENT(INOUT):: RAD_SHRT_NOS_SAVE  (ISPT*IJPHY,JLPHY,ICG)  
      REAL(8),INTENT(INOUT):: RAD_SHRT_SNW_SAVE  (ISPT*IJPHY,JLPHY,ICG)  
!
      REAL(8),INTENT(INOUT):: 
     1            R_STM_BAR_PAR_NOS_SAVE(ISPT*IJPHY,JLPHY,ICG) ,
     1            R_STM_BAR_PAR_SNW_SAVE(ISPT*IJPHY,JLPHY,ICG) 
!
! WORK
!
      INTEGER  :: I
      REAL(8)  :: FAC_NOS_NOS       (IDX)
      REAL(8)  :: FAC_SNW_SNW       (IDX)
      REAL(8)  :: RAD_SHRT_NOS      (IDX,ICG)
      REAL(8)  :: RAD_SHRT_SNW      (IDX,ICG)
!
      REAL(8),PARAMETER  :: R_STM_MAX = 1.D10 
!
! ================================
! >> 雪無し格子に雪が生じたとき <<
! ================================
!
#include "indep_fixedform"
      DO I = 1,NMAX
        J = ILIST(I) 
        IF ( INFO_SNOW_NEXT(I).GE.0 .AND. INFO_SNOW(I).LT.0 ) THEN
          RAD_SHRT_SNW_SAVE (J,JL,1) = RAD_SHRT_NOS_SAVE (J,JL,1)  
          RAD_SHRT_SNW_SAVE (J,JL,2) = RAD_SHRT_NOS_SAVE (J,JL,2)  
!
          R_STM_BAR_PAR_SNW_SAVE(J,JL,1)
     1                               = R_STM_BAR_PAR_NOS_SAVE(J,JL,1)
        ENDIF
      ENDDO
!
! ======================================
! >> 全面雪格子に雪無し域が生じたとき <<
! ======================================
!    下草は困るので、大きい値を与えてしまおう。
!
#include "indep_fixedform"
      DO I = 1,NMAX
        J = ILIST(I) 
        IF ( INFO_SNOW_NEXT(I).LE.0 .AND. INFO_SNOW(I).GE.1 ) THEN
          RAD_SHRT_NOS_SAVE (J,JL,1) = RAD_SHRT_SNW_SAVE (J,JL,1)  
          RAD_SHRT_NOS_SAVE (J,JL,2) = RAD_SHRT_SNW_SAVE (J,JL,2)  
          R_STM_BAR_PAR_NOS_SAVE(J,JL,1)
     1                               = R_STM_BAR_PAR_SNW_SAVE(J,JL,1)
          R_STM_BAR_PAR_NOS_SAVE(J,JL,2)
     1                               = R_STM_MAX 
        ENDIF
      ENDDO
!
! ============================
! >> 短波入射セーブ値の変更 <<
! ============================
!    全格子平均が変わらないように補正する。
!
      DO I=1,NMAX
        IF     ( CVR_SNW_NEXT(I) .EQ. CVR_SNW(I) ) THEN
          FAC_NOS_NOS(I) = 1.                                ! 無し -> 無し
          FAC_SNW_SNW(I) = 1.                                ! あり -> あり
        ELSEIF ( CVR_SNW_NEXT(I) .GT. CVR_SNW(I) ) THEN  ! 雪が増える
          FAC_NOS_NOS(I) = 1.                                ! 無し -> 無し
          FAC_SNW_SNW(I) = CVR_SNW(I) / CVR_SNW_NEXT(I)      ! あり -> あり
        ELSEIF ( CVR_SNW_NEXT(I) .LT. CVR_SNW(I) ) THEN  ! 雪が減る
          FAC_NOS_NOS(I) = ( 1 - CVR_SNW(I) ) 
     1                       / ( 1 - CVR_SNW_NEXT(I) )       ! 無し -> 無し
          FAC_SNW_SNW(I) = 1.                                ! あり -> あり
        ENDIF
      ENDDO
!
#include "indep_fixedform"
      DO I=1,NMAX
        J = ILIST(I) 
        IF ( CVR_SNW_NEXT(I) .NE. CVR_SNW(I) ) THEN  
          RAD_SHRT_NOS(I,1)                                   
     1       =   RAD_SHRT_NOS_SAVE (J,JL,1) * FAC_NOS_NOS(I)          
     1         + RAD_SHRT_SNW_SAVE (J,JL,1) * ( 1 - FAC_NOS_NOS(I) )  
          RAD_SHRT_NOS(I,2)                                   
     1       =   RAD_SHRT_NOS_SAVE (J,JL,2) * FAC_NOS_NOS(I)          
     1         + RAD_SHRT_SNW_SAVE (J,JL,2) * ( 1 - FAC_NOS_NOS(I) )  
!
          RAD_SHRT_SNW(I,1)                                   
     1       =   RAD_SHRT_NOS_SAVE (J,JL,1) * ( 1 - FAC_SNW_SNW(I) )  
     1         + RAD_SHRT_SNW_SAVE (J,JL,1) * FAC_SNW_SNW(I)          
          RAD_SHRT_SNW(I,2)                                   
     1       =   RAD_SHRT_NOS_SAVE (J,JL,2) * ( 1 - FAC_SNW_SNW(I) )  
     1         + RAD_SHRT_SNW_SAVE (J,JL,2) * FAC_SNW_SNW(I)          
        ENDIF
      ENDDO
!
      DO I=1,NMAX
        J = ILIST(I) 
        IF ( CVR_SNW_NEXT(I) .NE. CVR_SNW(I) ) THEN  
          RAD_SHRT_NOS_SAVE (J,JL,1) = RAD_SHRT_NOS      (I,1)    
          RAD_SHRT_NOS_SAVE (J,JL,2) = RAD_SHRT_NOS      (I,2)    
          RAD_SHRT_SNW_SAVE (J,JL,1) = RAD_SHRT_SNW      (I,1)    
          RAD_SHRT_SNW_SAVE (J,JL,2) = RAD_SHRT_SNW      (I,2)    
        ENDIF
      ENDDO
!
      RETURN
      END

      MODULE NMAX_INFO
!
      USE PRM , ONLY :
     1   IDIM  ,
     1   JDIM  ,
     1   JPHY2 ,
     1   JLPHY , 
     1   IJPHY , 
     1   INTRI , 
     1   INTRJ , 
     1   IRAD  , 
     1   ISPT
!
      use sibprm , only : idx 
!
      implicit none
!
      INTEGER, SAVE    ::  ILIST_NOSNOW    (IDX,JLPHY) 
      INTEGER, SAVE    ::  IMSKLT_NOSNOW   (IDX,JLPHY) 
      INTEGER, SAVE    ::  IMSKLT_RX_NOSNOW(IDX,JLPHY) 
      INTEGER, SAVE    ::  IMSKLT_RY_NOSNOW(IDX,JLPHY)
      INTEGER, SAVE    ::  NMAX_NOSNOW     (    JLPHY)     
      INTEGER, SAVE    ::  NTYP_NOSNOW     (IDX,JLPHY)
      INTEGER, SAVE    ::  NMAX_0_NOSNOW   (3  ,JLPHY)     
      INTEGER, SAVE    ::  NMAX_X_NOSNOW   (3  ,JLPHY)     
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      CONTAINS 
!----------------------------------------------------------
      SUBROUTINE NMAX_INFO_INI (
     I    IMASK             ) 
!
      USE SIBPRM , ONLY :
     1   IDX     
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: IMASK            (ISPT*IJPHY,JLPHY) 
!
      INTEGER :: IJ_SIB 
      INTEGER :: IJ_PHY
      INTEGER :: I_PHY 
      INTEGER :: J_PHY 
      INTEGER :: I_RAD
      INTEGER :: J_RAD
!
      INTEGER :: NM 
!
      INTEGER :: JL
      INTEGER :: I
!
#ifdef SIB_DEBUG
      INTEGER :: IMASK_TMP
      CHARACTER(1) :: CMASK ( IDIM , JDIM )
      INTEGER :: IS
      INTEGER :: ISTR 
      INTEGER :: IEND
#endif
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'NMAX_INFO 01/09 hosaka'
        LFIRST = .FALSE.
      ENDIF
!
! ---------
! > 表示 <
! ---------
!
#ifdef DDDDD
!
      WRITE(6,*) 'config_set : map information start'

#include "noparallel_fixedform"
      DO JL     =1,JLPHY
#include "noparallel_fixedform"
      DO IJ_PHY =1,IJPHY

        I_PHY = MOD ( IJ_PHY-1 , IDIM ) + 1 
        J_PHY = (IJ_PHY-1)/IDIM + 1 + JPHY2 * (JL-1)  
        IS    = 1 + (IJ_PHY-1)*ISPT

        IF     (       IMASK(IS,JL).GE.0 
     1           .AND. IMASK(IS,JL).LE.12 ) THEN
          IMASK_TMP = IMASK(IS,JL) 
        ELSEIF ( IMASK(IS,JL).LE.25 ) THEN
          IMASK_TMP = IMASK(IS,JL) - 12 
        ELSE
          write(6,*) 'config_set error IMASK' , 
     1               I_PHY , J_PHY , IMASK(IS,JL)  
          stop 999
        ENDIF
!        
        IF  (( IMASK_TMP .GT. 0 ).and.( IMASK_TMP .LT. 10 )) THEN 
          WRITE( CMASK(I_PHY,J_PHY) , '(I1)' ) IMASK_TMP
        ELSEIF ( IMASK_TMP .EQ. 10 ) THEN  
          CMASK(I_PHY,J_PHY) = 'A' 
        ELSEIF ( IMASK_TMP .EQ. 11 ) THEN  
          CMASK(I_PHY,J_PHY) = 'B' 
        ELSEIF ( IMASK_TMP .EQ. 12 ) THEN  
          CMASK(I_PHY,J_PHY) = 'C' 
        ELSEIF ( IMASK_TMP .EQ. 13 ) THEN  
          CMASK(I_PHY,J_PHY) = 'D' 
        ELSEIF ( IMASK_TMP .EQ. 0 ) THEN  
          CMASK(I_PHY,J_PHY) = ' ' 
        ELSE
          WRITE(6,*) 'CONFIG_SET ERROR I J IMASK' , 
     1               I_PHY , J_PHY , IMASK(IS,JL) , IMASK_TMP
          STOP 999
        ENDIF
!
      ENDDO
      ENDDO
!
      ISTR = 0
 1000 CONTINUE
      IEND = MIN ( ISTR+80 , IDIM )
#include "noparallel_fixedform"
      DO J_PHY=1,JDIM 
        WRITE(6,*) (CMASK(I,J_PHY),I=ISTR+1,IEND)
      ENDDO
      IF ( IEND .LT. IDIM ) GOTO 1000
!
      WRITE(6,*) 'config_set : map information end'
!
#endif 
!
      write(6,*) 'nmax_info imask start'
      do jl=1,jlphy
       write(6,*) jl , imask(:,jl)
      enddo
      write(6,*) 'nmax_info imask end'
!
! -------------------
! > *_NOSNOW の設定 <
! -------------------
!
      ILIST_NOSNOW    (:,:) = 0
      IMSKLT_NOSNOW   (:,:) = 0
      IMSKLT_RX_NOSNOW(:,:) = 0
      IMSKLT_RY_NOSNOW(:,:) = 0
      NTYP_NOSNOW     (:,:) = 0
      NMAX_NOSNOW     (  :) = 0
!
#include "noparallel_fixedform"
      DO JL=1,JLPHY
!
        NM = 0                        ! SiB 格子カウンタ 0<NM<IDX
#include "noparallel_fixedform"
        DO IJ_SIB = 1,IDX
          IJ_PHY = (IJ_SIB-1) / ISPT + 1 
!
          I_PHY   = MOD ( IJ_PHY - 1 , IDIM ) + 1
          J_PHY   = ( IJ_PHY - 1 ) / IDIM + 1 
!
          I_RAD   = (I_PHY-1)/INTRI + 1 
          J_RAD   = (J_PHY-1)/INTRJ + 1 
!
          IF ( IMASK(IJ_SIB,JL).GT.0 ) THEN
            NM = NM + 1 
            ILIST_NOSNOW     (NM,JL) = IJ_SIB       ! SiB -> 全SiB
            IMSKLT_NOSNOW    (NM,JL) = IJ_PHY       ! SiB -> 物理
            IMSKLT_RX_NOSNOW (NM,JL) = I_RAD
            IMSKLT_RY_NOSNOW (NM,JL) = J_RAD 
            NTYP_NOSNOW      (NM,JL) = IMASK(IJ_SIB,JL)  
          ENDIF
        ENDDO
!
        NMAX_NOSNOW          (JL)    = NM
        NMAX_0_NOSNOW        (1,JL)  = NM
        NMAX_0_NOSNOW        (2,JL)  =  1
        NMAX_0_NOSNOW        (3,JL)  = NM
        NMAX_X_NOSNOW        (1,JL)  =  0
        NMAX_X_NOSNOW        (2,JL)  = NM+1
        NMAX_X_NOSNOW        (3,JL)  = NM
!
#ifdef SIB_DEBUG
        write(6,*) 'config_set jl= ' , jl
        write(6,*) 'imask=         ' , (imask(i,jl),i=1,idx)
        write(6,*) 'imsklt_nosnow= ' , (imsklt_nosnow(i,jl),i=1,nm)
#endif
!
      ENDDO
!
      RETURN
      END SUBROUTINE NMAX_INFO_INI 
! --------------------------------------------------------     
      SUBROUTINE NMAX_INFO_RUN_NOSNOW (
     I   JL         , 
     O   NMAX       , NTYP       ,     
     O   IMSKLT     , IMSKLT_RX  , IMSKLT_RY  , ILIST )
!
      IMPLICIT NONE
!
! INPUT START 
!
      INTEGER,INTENT(IN)    ::  JL 
!
! OUTPUT START
!
      INTEGER,INTENT(OUT)   ::  ILIST    ( IDX )   ! SiB -> 全 SiB 
      INTEGER,INTENT(OUT)   ::  IMSKLT   ( IDX )   ! SiB -> 物理
      INTEGER,INTENT(OUT)   ::  IMSKLT_RX( IDX )   ! SiB -> 放射X
      INTEGER,INTENT(OUT)   ::  IMSKLT_RY( IDX )   ! SiB -> 放射Y
      INTEGER,INTENT(OUT)   ::  NTYP     ( IDX )   ! 植生タイプ（陸面格子）
      INTEGER,INTENT(OUT)   ::  NMAX
!
#ifdef SIB_DEBUG
      IF ( LFIRST ) THEN
        WRITE(6,*) 'NMAX_INFO_RUN_NOSNOW : ERROR '
        WRITE(6,*) 'SHOULD BE CALLED AFTER CALL NMAX_INFO_INI'
        STOP 999
      ENDIF
#endif
!
      NMAX         = NMAX_NOSNOW     (  JL) 
      NTYP     (:) = NTYP_NOSNOW     (:,JL)
      ILIST    (:) = ILIST_NOSNOW    (:,JL)
      IMSKLT   (:) = IMSKLT_NOSNOW   (:,JL)
      IMSKLT_RX(:) = IMSKLT_RX_NOSNOW(:,JL)
      IMSKLT_RY(:) = IMSKLT_RY_NOSNOW(:,JL)
!
      RETURN    
      END SUBROUTINE NMAX_INFO_RUN_NOSNOW
!-----------------------------------------------------------------
      SUBROUTINE NMAX_INFO_RUN_SNOW (
     I   IMASK      , INFO_SNOW  , JL         , 
     O   NMAX       , NTYP       ,     
     O   NMAX_0     , NMAX_01    , NMAX_1     , 
     O   NMAX_2     , NMAX_3     , NMAX_4     ,     
     O   IMSKLT     , IMSKLT_RX  , IMSKLT_RY  , ILIST )
!
      IMPLICIT NONE
!
! INPUT START 
!
      INTEGER,INTENT(IN)    ::  JL 
      INTEGER,INTENT(IN)    ::  IMASK    (ISPT*IJPHY)! 植生タイプ
      INTEGER,INTENT(IN)    ::  INFO_SNOW(ISPT*IJPHY)! 雪の層数
!
!
!
! OUTPUT START
!
      INTEGER,INTENT(OUT)   ::  ILIST    ( IDX )   ! SiB -> 全 SiB 
      INTEGER,INTENT(OUT)   ::  IMSKLT   ( IDX )   ! SiB -> 物理
      INTEGER,INTENT(OUT)   ::  IMSKLT_RX( IDX )   ! SiB -> 放射X
      INTEGER,INTENT(OUT)   ::  IMSKLT_RY( IDX )   ! SiB -> 放射Y
      INTEGER,INTENT(OUT)   ::  NTYP     ( IDX )   ! 植生タイプ（陸面格子）
      INTEGER,INTENT(OUT)   ::  NMAX
      INTEGER,INTENT(OUT)   ::  NMAX_0 ( 3 )     ! 雪無し 個数・開始・終了
      INTEGER,INTENT(OUT)   ::  NMAX_01( 3 )     ! 雪01   個数・開始・終了
      INTEGER,INTENT(OUT)   ::  NMAX_1 ( 3 )     ! 雪1    個数・開始・終了
      INTEGER,INTENT(OUT)   ::  NMAX_2 ( 3 )     ! 雪2    個数・開始・終了
      INTEGER,INTENT(OUT)   ::  NMAX_3 ( 3 )     ! 雪3    個数・開始・終了
      INTEGER,INTENT(OUT)   ::  NMAX_4 ( 3 )     ! 雪4    個数・開始・終了
!
      INTEGER               :: I 
!
#ifdef SIB_DEBUG
      IF ( LFIRST ) THEN
        WRITE(6,*) 'NMAX_INFO_RUN_SNOW : ERROR '
        WRITE(6,*) 'SHOULD BE CALLED AFTER CALL NMAX_INFO_INI'
        STOP 999
      ENDIF
#endif
!
! 初期化
!
      NMAX       = 0
      NMAX_0 (:) = 0
      NMAX_01(:) = 0
      NMAX_1 (:) = 0
      NMAX_2 (:) = 0
      NMAX_3 (:) = 0
      NMAX_4 (:) = 0
!
! 全ての場所で雪がない場合
!
#include "noparallel_fixedform"
      DO 10 I=1,ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).GE.0 ) GOTO 50   
!                                                ! 雪が一箇所でもあれば飛べ
 10   CONTINUE
      NMAX         = NMAX_NOSNOW     (  JL) 
      NMAX_0   (:) = NMAX_0_NOSNOW   (:,JL) 
      NMAX_01  (:) = NMAX_X_NOSNOW   (:,JL) 
      NMAX_1   (:) = NMAX_X_NOSNOW   (:,JL) 
      NMAX_2   (:) = NMAX_X_NOSNOW   (:,JL) 
      NMAX_3   (:) = NMAX_X_NOSNOW   (:,JL) 
      NMAX_4   (:) = NMAX_X_NOSNOW   (:,JL) 
      NTYP     (:) = NTYP_NOSNOW     (:,JL)
      ILIST    (:) = ILIST_NOSNOW    (:,JL)
      IMSKLT   (:) = IMSKLT_NOSNOW   (:,JL)
      IMSKLT_RX(:) = IMSKLT_RX_NOSNOW(:,JL)
      IMSKLT_RY(:) = IMSKLT_RY_NOSNOW(:,JL)
!
#ifdef SIB_DEBUG
      write(6,*) 'set_nmax nosnow ' , nmax , ' imsklt ' , 
     1            (imsklt(i),i=1,nmax)
      write(6,*) 'set_nmax nosnow ' , nmax , ' ilist  ' , 
     1            (ilist (i),i=1,nmax)
      write(6,*) 'set_nmax nosnow ' , nmax , ' imsklt_rx ' , 
     1            (imsklt_rx(i),i=1,nmax)
      write(6,*) 'set_nmax nosnow ' , nmax , ' imsklt_ry ' , 
     1            (imsklt_ry(i),i=1,nmax)
#endif
!
      RETURN    
!
! 雪がある場合の設定ここから
!
 50   CONTINUE
!
!       ILIST  ... SIB 格子 -> 全 SiB 格子
!       IMSKLT ... SIB 格子 -> 物理格子
!
! 雪無し 
!
      NMAX_0(2) = 1
#include "noparallel_fixedform"
      DO 100 I = 1, ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).LT.0 ) THEN
          NMAX      = NMAX+1
          NMAX_0(1) = NMAX_0(1)+1
          NTYP     (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  100 CONTINUE
      NMAX_0(3) = NMAX 
!
! 雪10 
!
      NMAX_01(2) = NMAX+1
#include "noparallel_fixedform"
      DO 200 I = 1, ISPT*IJPHY 
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).EQ.0 ) THEN
          NMAX       = NMAX       + 1
          NMAX_01(1) = NMAX_01(1) + 1
          NTYP  (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  200 CONTINUE
      NMAX_01(3) = NMAX
!
! 雪1

      NMAX_1(2) = NMAX + 1 
#include "noparallel_fixedform"
      DO 400 I = 1, ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).EQ.1 ) THEN
          NMAX     = NMAX+1
          NMAX_1(1)= NMAX_1(1)+1
          NTYP  (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  400 CONTINUE
      NMAX_1(3) = NMAX
!
! 2 層
!
      NMAX_2 (2) = NMAX + 1 
#include "noparallel_fixedform"
      DO 510 I = 1, ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).EQ.2 ) THEN 
          NMAX       = NMAX       + 1
          NMAX_2(1) = NMAX_2(1) + 1
          NTYP  (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  510 CONTINUE
      NMAX_2 (3) = NMAX
!
! 3 層
!

      NMAX_3 (2) = NMAX + 1 
#include "noparallel_fixedform"
      DO 610 I = 1, ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).EQ.3 ) THEN
          NMAX       = NMAX       + 1
          NMAX_3(1) = NMAX_3(1) + 1
          NTYP  (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  610 CONTINUE
      NMAX_3 (3) = NMAX
!
! 雪4
!
      NMAX_4(2) = NMAX+1
      DO 700 I = 1, ISPT*IJPHY
        IF ( IMASK(I).GE.1 .AND. INFO_SNOW(I).EQ.4 ) THEN
          NMAX      = NMAX      + 1
          NMAX_4(1) = NMAX_4(1) + 1
          NTYP  (NMAX) = MOD(IMASK(I),100)
!
          ILIST    (NMAX) = I
          IMSKLT_RX(NMAX) = MOD((I-1)/INTRI,IRAD)+1
          IMSKLT_RY(NMAX) = ((I-1)/IDIM)/INTRJ+1
!
        END IF
  700 CONTINUE
      NMAX_4(3) = NMAX
!
      DO 800 I=1,NMAX
        IMSKLT(I) = (ILIST(I)-1)/ISPT + 1 
 800  CONTINUE   
!
! デバッグ時のみ
!
#ifdef SIB_DEBUG
      IF ( NMAX .LT. NMAX_NOSNOW(JL) ) THEN
        write(6,*) 'setnmax error : NMAX_NOSNOW = ' , NMAX_NOSNOW(JL) , 
     1                            '  nmax = ' , nmax 
        stop 999
      ENDIF
#endif
!
#ifdef SIB_DEBUG
      write(6,*) 'set_nmax snow  ' , nmax , ' imsklt ' , 
     1            (imsklt(i),i=1,nmax)
      write(6,*) 'set_nmax snow  ' , nmax , ' ilist  ' , 
     1            (ilist (i),i=1,nmax)
      write(6,*) 'set_nmax snow  ' , nmax , ' imsklt_rx ' , 
     1            (imsklt_rx(i),i=1,nmax)
      write(6,*) 'set_nmax snow  ' , nmax , ' imsklt_ry ' , 
     1            (imsklt_ry(i),i=1,nmax)
      write(6,*) 'set_nmax snow  ' , 
     1            'nmax_0  ' , nmax_0  ,
     1            'nmax_01 ' , nmax_01 ,
     1            'nmax_1  ' , nmax_1  ,
     1            'nmax_2  ' , nmax_2  ,
     1            'nmax_3  ' , nmax_3  ,
     1            'nmax_4  ' , nmax_4
#endif
!
      RETURN
      END SUBROUTINE NMAX_INFO_RUN_SNOW
!----------------------------------------------------------------
      END MODULE NMAX_INFO


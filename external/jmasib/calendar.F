! cvdate4.F - GSM に由来する日時計算
!
! subroutine CDSET()
!       com_runconf を参照して日数テーブルを初期化する
! entry GETID(IDIN, IDOUT, KTLAG)
!       IDIN に KTLAG を加えた IDOUT を出力
! entry GETKT(IDIN1, IDIN2, IFLAG, KTLAG)
!       IDIN1 から IDIN2 までの時間差を KTLAG を出力
!       IFLAG =4  時間単位 , 3 日単位
! entry GETSBD(IDIN, MONS, MONE, WGT) 
!       月ごとに与えられるデータの補間のための重みを計算
!
! ID で始まる名前の変数は大きさ5の整数型の配列で、以下のように日時を格納する。
!
!       id(1): 西暦年
!       id(2): 月
!	id(3): 日
!	id(4): 時
!	id(5): 未使用 (0で初期化)

!==================================================================
! COM_RUNCONF モジュールの日付情報を参照して日数テーブルを初期化する。
!
      MODULE CALENDAR 
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: IYEARM = 200

      INTEGER,PARAMETER :: 
     +     MON(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
!    +     MON(12) = (/ 30,30,30,30,30,30,30,35,30,30,30,30 /)

      INTEGER,PARAMETER :: IYSTR = 1900
!
      INTEGER,SAVE :: MDAY ( IYEARM, 0:25 ) 
!!           0    --- 西暦
!!           1-12 --- 月の日数
!!          13-24 --- その年の 1/1 0 時から前月の終りまでの時間数
!!          25    --- IYSTR 1/1 0h からその年の 1/1 0h までの時間数
!
      LOGICAL,SAVE :: FIRST = .TRUE.
!
      CONTAINS 
!==========================================================
      SUBROUTINE  CALENDAR_INI
!
      USE COM_RUNCONF_SIB0109 , ONLY : JCN_SIB_ISLSCP, JCN_SIB_ATM_CYCLE 
!
      IMPLICIT NONE
!
      INTEGER :: II
      INTEGER :: J
      INTEGER :: IURU

!
      IF ( FIRST ) THEN
!
        WRITE(6,*) 'CALENDAR_INI : CALENDAR IS SET FROM ' , IYSTR , 
     1              ' TO ' , IYSTR+IYEARM-1 
!
!!  年 MDAY(II,0) と 各月の日数 MDAY(II,1-12) の設定
!
        DO 100 II=1,IYEARM
          MDAY(II,0) = IYSTR + II - 1 
          DO 110 J=1,12
            MDAY(II,J) = MON(J)
 110      CONTINUE
!
!     閏年処理
!       ISLSCP の場合は 1988 のみ
          IF ( JCN_SIB_ISLSCP .EQ. 1 ) THEN 
            IF ( MDAY(II,0) .EQ. 1988 ) THEN        ! 1988 だけ閏年あり
              MDAY(II,2) = MDAY(II,2) + 1 
            ENDIF
!       それ以外
          ELSE 
!           ある特定の年でサイクルする場合は強制的に年をリセット
            IF( JCN_SIB_ATM_CYCLE.NE.0 ) THEN
                  IURU=JCN_SIB_ATM_CYCLE
            ELSE
                  IURU=MDAY(II,0)
            ENDIF
            IF ( (MOD(IURU,400).EQ.0) .OR. 
     -           (MOD(IURU,4).EQ.0.AND.MOD(IURU,100).NE.0) )THEN
              MDAY(II,2) = MDAY(II,2) + 1 
            ENDIF
          ENDIF
 100    CONTINUE
!
!! 前月までの時間の設定
!
        DO 200 II=1,IYEARM
          IF ( II .EQ. 1 ) THEN
            MDAY(II,1+24) = 0
          ELSE
            MDAY(II,1+24) = MDAY(II-1,1+24)       ! 前々年まで
     1                     + MDAY(II-1,24)        ! 前年 12/1 まで      
     1                     + MDAY(II-1,12)*24     ! 前年 12 月分
          ENDIF
!
          MDAY(II,1+12) = 0 
!
          DO 210 J=2,12 
            MDAY(II,J+12) = MDAY(II,J-1+12) + MDAY(II,J-1)*24
  210     CONTINUE
  200   CONTINUE
!
#ifdef DEBUG
        WRITE(6,*) '  CDSET : EXAMPLES ARE FOLLOWING. ' 
        DO II =1,IYEARM/10 
          WRITE(6,*) '   YEAR IS ' , MDAY(II,0), 
     1               ' DAYS ARE ' , (MDAY(II,J),J=1,25)
        ENDDO
        DO II =IYEARM*9/10 , IYEARM
          WRITE(6,*) '   YEAR IS ' , MDAY(II,0), 
     1               ' DAYS ARE ' , (MDAY(II,J),J=1,25)
        ENDDO
#endif
! 20050516
! /* DEBUG */
        FIRST = .FALSE. 
!
      ELSE 
        WRITE(6,*) 'CALENDAR_INI : ALREADY SET , SO RETURN '
      ENDIF
!
      RETURN
      END SUBROUTINE CALENDAR_INI
!-----------------------------------------------------------------
      SUBROUTINE CALENDAR_RUN_GETID  ( ID , IDNEW , KTLAG )
!
      implicit none
!
      integer,intent(in)  :: id(5) 
      integer,intent(out) :: idnew(5) 
      integer,intent(in)  :: ktlag
!
      integer  :: ii 
      integer  :: kt1
      integer  :: kt2
      integer  :: m

!
!! 入力
!!    ID   (5) :  もとになる日時  ( 年, 月, 日, 時, ？ )
!!    KTLAG    :  時間差  
!
!! 出力      
!!    IDNEW(5) :  計算される日時  ( 年, 月, 日, 時, ？ )
!                 IDNEW = ID  +  KTLAG (HOUR) 
!
!! < 初期設定 >
!
#ifdef PCHECK
      IF ( FIRST ) THEN
        WRITE(6,*) 'GETID : CALENDAR WAS NOT SET.' 
        STOP 999
      ENDIF  
#endif
!
      IDNEW(1) = ID(1)
      IDNEW(5) = 0                          ! 意味無し。
!
!! < ID の 1/1 からの通算時間 >             ! 例 IYSTR = 1970 で 
!                                           ! 92/12/31/18 に + 12 する
!
      II  = ID(1) - IYSTR + 1               ! ID(1)=1992 -> II= 23
      KT1 = ID(4) + ( ID(3)-1 ) * 24               
     1            + MDAY(II,ID(2)+12)       ! KT1 =18+30*24+334*24
!                                           !               ^^^11/30 まで 
!! < IDNEW の 1/1 からの通算時間 > 
!
      KT2 = KT1 + KTLAG                     ! KT2 = 6+365*24
!
!! < IDNEW(1-4) を決める > 
!
 1100 CONTINUE 
!
      IF ( KT2 .LT. 0 ) THEN   ! 前の年 ( またはそれより前 )             
        IDNEW(1) = IDNEW(1) - 1                    
        II       = II - 1                          
        KT2      = KT2 + ( MDAY(II+1,25) - MDAY(II,25) ) 
        GOTO 1100
!
      ELSEIF ( KT2 .GE. ( MDAY(II+1,25) - MDAY(II,25) )  ) THEN !翌年 ( 以降 )
        IDNEW(1) = IDNEW(1) + 1 
        II       = II + 1 
        KT2      = KT2  - ( MDAY(II,25) - MDAY(II-1,25) )
        GOTO 1100 
!
      ELSE                      ! その年
        DO 1200 M = 2, 12 
          IF ( KT2 .LT. MDAY(II,M+12) ) THEN   ! その前の月
            IDNEW(2) = M-1 
            GOTO 1300 
          ENDIF
 1200   CONTINUE
        IDNEW(2) = 12         
!
 1300   CONTINUE
!
        KT2 = KT2 - MDAY(II,IDNEW(2)+12)       ! その月の 1 日からの時間
        IDNEW(3) = KT2 / 24 + 1 
        IDNEW(4) = MOD ( KT2 , 24 ) 
      ENDIF
!
      RETURN      
      END SUBROUTINE CALENDAR_RUN_GETID
!-----------------------------------------------------------------
      SUBROUTINE CALENDAR_RUN_GETKT ( 
     I     ID    , IDNEW , IFLAG   , 
     O     KTLAG )
!
!! 入力
!!    ID   (5) :  もとになる日時  ( 年, 月, 日, 時, ？ )
!!    IDNEW(5) :  欲しい日時      ( 年, 月, 日, 時, ？ )
!!    IFLAG     :  フラグ. 時間単位 4   日単位 3
!!         ただし 1/1 12h と 1/2 1h の日差、1/2 14h の日差はともに 1 日
!
!! 出力      
!!    KTLAG    :  時間差  
!                 KTLAG = IDNEW - ID
!
      implicit none
!
      integer,intent(in) :: id(5)
      integer,intent(in) :: idnew(5)
      integer,intent(in) :: iflag
      integer,intent(out):: ktlag
!
      integer :: i1
      integer :: kt1
      integer :: i2
      integer :: kt2
      integer :: kt3
!
#ifdef PCHECK
      IF ( FIRST ) THEN
        WRITE(6,*) 'GETKT : CALENDAR WAS NOT SET.' 
        STOP 999
      ENDIF  
#endif
!
!! < ID の 1/1 からの通算時間 >     
!                                   
      I1  = ID(1) - IYSTR + 1       
      KT1 = ID(4) + ( ID(3)-1 ) * 24               
     1      + MDAY(I1,ID(2)+12) 
!
!! < IDNEW の 1/1 からの通算時間 >     
!                                   
      I2  = IDNEW(1) - IYSTR + 1       
      KT2 = IDNEW(4) + ( IDNEW(3)-1 ) * 24               
     1         + MDAY(I2,IDNEW(2)+12)
!
!! < ID の 1/1 と IDNEW の 1/1 との時間差 ( 年が違う場合 0 でない ) >
!
      KT3 = MDAY(I2,25) - MDAY(I1,25)
!
!! < KTLAG の計算 > 
!
      IF ( IFLAG .EQ. 4 ) THEN          ! 時間単位 
        KTLAG = KT3 + KT2 - KT1 
      ELSEIF ( IFLAG .EQ. 3 ) THEN
        KTLAG = KT3 / 24 + KT2 / 24 - KT1 / 24 
#ifdef PCHECK
      ELSE 
        WRITE(6,*) 'GETKT : NOT SUPPORT IFLAG = ' , IFLAG 
        STOP 999
#endif
      ENDIF
!
      RETURN        
      END SUBROUTINE CALENDAR_RUN_GETKT       
!-------------------------------------------------------------------
      SUBROUTINE CALENDAR_RUN_GETSBD  (
     I    IDNEW , 
     O    MONS  , MONE , IYEARS , IYEARE , WGT ) 
!
      implicit none
!
      integer,intent(in)  :: idnew(5)
      integer,intent(out) :: mons
      integer,intent(out) :: mone
      integer,intent(out) :: iyears
      integer,intent(out) :: iyeare
      real(8),intent(out) :: wgt
!
      integer :: iy
!
!       月ごとに与えられるデータの補間のための重みを計算する。
!       IDNEW を与えると、直前の月央 (IYEARS 年 MONS 月) と
!       直後の月央 (IYEARE 年 MONE 月) を求め、この2者を含む線形補間で
!       (1.0 - WGT) * (MONS 月のデータ) + WGT * (MONE 月のデータ)
!	のように IDNEW のデータを得るための係数 WGT を求める。
!       月央は月末の15日前であり、1 月は 30 日であると仮定されている。
!
!! WGT は, MONE の月の割合. 
!! 重み自体はちょっと変な計算で求める. いずれ変更. 
!! 
#ifdef PCHECK
      IF ( FIRST ) THEN
        WRITE(6,*) 'GETSBD : CALENDAR WAS NOT SET.' 
        STOP 999
      ENDIF  
#endif
!
      IY = IDNEW(1) - IYSTR + 1 
      IYEARS = IDNEW(1)
      IYEARE = IDNEW(1)
!
      IF ( MDAY(IY,IDNEW(2)) - IDNEW(3) .LE. 15 ) THEN  ! 月の後半
        MONS = IDNEW(2) 
        MONE = IDNEW(2)+1 
        WGT  = ( 15 - ( MDAY(IY,IDNEW(2)) - IDNEW(3) ) ) / 30.0D0 
!       WGT  = ( 15 - ( MDAY(IY,IDNEW(2)) - IDNEW(3) ) ) / 30.0 
      ELSE                                              ! 月の前半
        MONE = IDNEW(2) 
        MONS = IDNEW(2)-1
        WGT  = ( 45 - ( MDAY(IY,IDNEW(2)) - IDNEW(3) ) ) / 30.0D0 
!       WGT  = ( 45 - ( MDAY(IY,IDNEW(2)) - IDNEW(3) ) ) / 30.0 
      ENDIF
!
      IF ( MONE .EQ. 13 ) THEN
        MONE   = 1 
        IYEARE = IYEARE + 1
      ELSEIF ( MONS .EQ. 0 ) THEN
        MONS   = 12 
        IYEARS = IYEARS - 1 
      ENDIF
!
      RETURN
!-------------------------------------------------------------------
      END SUBROUTINE CALENDAR_RUN_GETSBD
      END MODULE CALENDAR
!===================================================================

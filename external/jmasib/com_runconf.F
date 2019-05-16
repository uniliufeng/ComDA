      MODULE COM_RUNCONF          ! オフライン用
      

!! 実行ジョブ設定情報 ( ジョブ内で一定 )

!! 下の設定値はdefaultの値。
!! com_runconf_iniで、いくつかの値が置き換わる。

! run configuration information module file
!
!--- FLG0 
!
!     INTEGER,SAVE :: 
!    1   JCN_KT0_START = 0   ! KT=0 からのスタートなら 1 
!     NAMELIST / NAMFLG0 / JCN_KT0_START 
!
!--- FLG1 ありなしフラグ ( 物理関係 )

      INTEGER,SAVE :: 
     1   JCNLSM   = 1    ,    ! 1 L3SiB ( 大泉版 ) , 2 SiB0007 
     1   JCNTTICE = 1  
!
      INTEGER,SAVE ::          ! モニタ
     1   JCNMMM   = 0,         ! 柴田モニタ               980730 保坂
     1   JCNGRP   = 1          ! EP-FLUX and 
!                              ! GRIPS MONITOR, i.e,   zonal mean daily values
!                                     T(90N),U(60N), T(90S),U(60S)
!
      INTEGER, SAVE :: 
     1  JCNIMNT    = 999, ! 月毎モニタ出力間隔。 980804 保坂 旧 INTMNT
     1  JCNIMNTDAY =  24, ! 日毎モニタ出力間隔。
     1  JCNIMNT6HR =   6

      NAMELIST / SIB0109_COM_RUNCONF /       
     1   JCNIMNT, JCNIMNTDAY, JCNIMNT6HR
!
      CONTAINS
!
! ===============================================================
!
      SUBROUTINE COM_RUNCONF_INI
!
        use com_stdinout_unit , only : iunit_in , iunit_out  
!
       !--- ネームリスト読み込みによる変更 
!
        REWIND(iunit_in)
        READ  (iunit_in ,SIB0109_COM_RUNCONF)
        WRITE (iunit_out,SIB0109_COM_RUNCONF)

      END SUBROUTINE COM_RUNCONF_INI
!====================================================================
      END MODULE COM_RUNCONF



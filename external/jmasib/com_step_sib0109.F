! module com_step 
! 99/07/20 保坂征宏  モデルのうち必要なものを切り出す。
!                    ICN_SIB_CALC_SOIL_SNOW と ICNSW_NEXT は新規作成
!
!  各モデル対応の COM_STEP_SIB0007_RUN_* が必要
!  とにかく名前は 31 文字までに留めましょう
!
      MODULE COM_STEP_SIB0109 
!
      IMPLICIT NONE
!
!                      1234567890123456789012345678901
      INTEGER,SAVE  :: ICN_SIB0109_CALC_SOIL_SNOW    ! SiB で土壌・雪の計算
      INTEGER,SAVE  :: ICN_SIB0109_SNOW_RESET     =0 ! 雪層をリセット
      INTEGER,SAVE  :: ICN_SIB0109_AFTER_CALC_SOIL_SNW  
                                                     ! SiB で土壌・雪の計算した
                                                     ! 次のステップ
      INTEGER,SAVE  :: ICN_SIB0109_SW                ! 短波計算するステップ
      INTEGER,SAVE  :: ICN_SIB0109_LW                ! 長波計算するステップ
!
      INTEGER,SAVE  :: ICN_SIB0109_MONITOR           ! モニタ出力するステップ
!
      REAL(8),SAVE  :: SCN_SIB0109_DELT            ! 時間刻
      REAL(8),SAVE  :: SCN_SIB0109_DELT_ATM        ! 大気の時間刻 
                                                   !  or リープフロッグのとび
      REAL(8),SAVE  :: SCN_SIB0109_DELT_SR         ! 短波の時間刻 ( 1 hour ? )
!
      INTEGER,SAVE  :: ICN_SIB0109_FIRST_STEP        ! ジョブの第一ステップ 
      INTEGER,SAVE  :: ICN_SIB0109_LAST_STEP         ! ジョブの最終ステップ 
!
      contains 
!===================================================================
      SUBROUTINE COM_STEP_SIB0109_RUN  (
     I    FSECP ) 
!
      USE COM_RUNCONF_SIB0109 , ONLY : 
     I    JCN_IWL_SKIP      ,
     I    RCN_DELT_CNP      ,  
     I    RCN_DELT_SOIL     ,
     1    RCN_DELT_SR     
!
      USE COM_STEP , ONLY :
     I    ICNSW   , 
     I    ICNIR   , 
     I    ICNMNTMON , 
     I    SCN_DELT  , 
     I    SCN_2DELTX ,
!    I    SCNDTSW   ,
     I    ICNFST  ,
     I    ICNLAST

!
!
      IMPLICIT NONE
!
      real(8),intent(in) :: fsecp
      integer  :: is 
!
      logical , save :: lfirst = .true. 
!
!   土壌・雪を計算した直後フラグ
!
      IF ( LFIRST ) THEN
        ICN_SIB0109_AFTER_CALC_SOIL_SNW = 1 
        lfirst = .false. 
      ELSE 
        ICN_SIB0109_AFTER_CALC_SOIL_SNW = ICN_SIB0109_CALC_SOIL_SNOW 
      ENDIF 
!
      IF     ( 
     1              JCN_IWL_SKIP.EQ.1 .OR. JCN_IWL_SKIP.EQ.-1 
     1         .OR. JCN_IWL_SKIP.EQ.2 .OR. JCN_IWL_SKIP.EQ.-2 ) THEN 
        ICN_SIB0109_CALC_SOIL_SNOW = 1
!
      ELSEIF ( JCN_IWL_SKIP.EQ.3 .OR. JCN_IWL_SKIP.EQ.-3 ) THEN  
!
        IS = ( FSECP + 0.01 ) / RCN_DELT_SOIL
        IF ( ABS(FSECP-IS*RCN_DELT_SOIL) .LT. 0.1 ) THEN
          ICN_SIB0109_CALC_SOIL_SNOW = 1
        ELSE
          ICN_SIB0109_CALC_SOIL_SNOW = 0
        ENDIF
!
      ENDIF
!
      ICN_SIB0109_SW = ICNSW 
      ICN_SIB0109_LW = ICNIR
!
#ifdef SIB_DEBUG
      write(6,*) 'com_step_sib0109_run_stepini SW LR SS' , 
     1      ICN_SIB0109_SW , 
     1      ICN_SIB0109_LW , 
     1      ICN_SIB0109_CALC_SOIL_SNOW
#endif
!
      ICN_SIB0109_MONITOR = ICNMNTMON 
!
      SCN_SIB0109_DELT     = SCN_DELT 
!
      SCN_SIB0109_DELT_ATM = SCN_2DELTX 
!
!     SCN_SIB0109_DELT_SR  = SCNDTSW 
      SCN_SIB0109_DELT_SR  = RCN_DELT_SR 
!
      ICN_SIB0109_FIRST_STEP = ICNFST
      ICN_SIB0109_LAST_STEP  = ICNLAST
!
      end subroutine com_step_sib0109_run
!===================================================================
      END MODULE COM_STEP_SIB0109


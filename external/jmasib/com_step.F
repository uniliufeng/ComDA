! module com_step 
! 99/07/20 保坂征宏  モデルのうち必要なものを切り出す。
! 00/08/24 保坂征宏  mj98 にあるもののみ残す。残りは com_step_sib0007 へ
!
!  ステップ毎に変わる変数/フラグ
!  サブルーチン time_step で設定されている。
!
      MODULE COM_STEP 
!
      IMPLICIT NONE
!
      REAL(8),SAVE  :: SCN_DELT                ! 時間刻そのもの 
      REAL(8),SAVE  :: SCN_2DELTX              ! リープフロッグのとび
!     REAL(8),SAVE  :: SCNDTSW                 ! 短波放射用 DELT 時間重み
      INTEGER,SAVE  :: ICNSW                   ! 短波放射計算
      INTEGER,SAVE  :: ICNIR                   ! 長波放射計算
      INTEGER,SAVE  :: ICNMNTMON               ! 月モニタを出力する
      INTEGER,SAVE  :: ICNFST                  ! ラン始めのステップ 
      INTEGER,SAVE  :: ICNLAST                 ! ラン最終ステップ
!
      END MODULE 






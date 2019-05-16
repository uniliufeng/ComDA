! 99/07/07 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
! 2000-09-11 豊田英司  水の上限処理の追加に伴う抜本書き換え
!
      SUBROUTINE SOIL_MOD (
     I   NMAX, FACT_WGT_H2O, INFO_SNOW_NEXT, 
     M   TMP_SOIL, TMP_SNSL, 
     M   FR_WTR_SOIL, FR_WTR_SNSL, FR_GLA_SOIL, FR_GLA_SNSL) 
!
! 土壌水比・氷比が物理的にありえない値になっている場合の処置をする。
! 温度についても将来何らかの処置がなされるかもしれない。
!
! 新しいステップの積雪被覆率で第1層の水・氷を再分配したあとで行う。
!
! 現状では以下の処置をしている。ここで含有率とは空隙に対する比である。
! 積雪域に対する処理をまとめてもよいが、abc の順番を変えてはいけない。
!
!   a 土壌液水含有率 が負である場合強制的に 0 にする
!   b 土壌氷含有率 が負である場合強制的に 0 にする
!   c 土壌水分(液水+氷)総含有率 が 1 を超える場合 1 にする
!
!   d 土壌総水分(液水+氷)含有率 が
!     が 0.05 正しくは SIBCON::FR_SOIL_H2O_MIN) より小さいときに、
!     同じ温度の H2O を与えて 0.05 になるようにする
!
! 土壌総水分に関する強制補正時はなるべく現在の存在度で分配している。
! かつては
!
!   e 土壌氷含有率が 0.95 (正しくは SIBCON::FR_SOIL_GLA_MAX) を超えた
!     ときに 0.95 になるようにする
!
! もやっていた。これは今のところやっていない。
!
! 高速化に関して速度測定は実行していないが、旧 SOIL_MOD にくらべて
! 修正対象のかき集めをしないようにしているので、有利だと考えられる。
! SR8000 にて配列代入文が DO ループと同様に疑似ベクトル化・並列化の
! 対象となることは確認している。
!
      USE SIBPRM , ONLY :
     1   IDX   , 
     1   IDP   
!
      USE SIBCON , ONLY :
     1    FR_SOIL_H2O_MIN_c  , 
     1    FR_SOIL_GLA_MAX_c  , 
     1    RHO_WTR_c        , 
     1    CAP_GLA_c        , 
     1    CAP_WTR_c        , 
     1    TMP_FREZ_c       , 
     1    HL_W2G_c         
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB_NEXT , ICN , ICS ,
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER, INTENT(IN):: NMAX
      REAL(8), INTENT(IN):: FACT_WGT_H2O(IDX, IDP)  ! 空隙量x水の比重
      INTEGER, INTENT(IN):: INFO_SNOW_NEXT(IDX)
!
! INOUT
!
      ! 土壌の温度
      REAL(8), INTENT(INOUT) :: TMP_SOIL   (IDX,IDP) 
      REAL(8), INTENT(INOUT) :: TMP_SNSL   (IDX,IDP) 
      ! 土壌水比
      REAL(8), INTENT(INOUT) :: FR_WTR_SOIL(IDX,IDP) 
      REAL(8), INTENT(INOUT) :: FR_WTR_SNSL(IDX,IDP) 
      ! 土壌氷比
      REAL(8), INTENT(INOUT) :: FR_GLA_SOIL(IDX,IDP) 
      REAL(8), INTENT(INOUT) :: FR_GLA_SNSL(IDX,IDP) 
!
! WORK
!
      REAL(8):: H2O_MOD_SOIL(IDX,IDP)
      REAL(8):: H2O_MOD_SNSL(IDX,IDP)
      REAL(8):: ENG_MOD_SOIL(IDX,IDP)
      REAL(8):: ENG_MOD_SNSL(IDX,IDP)
      ! 動作フラグ
      LOGICAL:: BARESOIL_FLAG (IDX)    ! 無雪域をもつ格子だけ真
      LOGICAL:: UNDERSNOW_FLAG(IDX)    ! 積雪域をもつ格子だけ真
      ! 配列探査用 DO 変数
      INTEGER:: I
      INTEGER:: K
      ! 配列探査中に演算途中結果を置く変数。POPTION TLOCAL してよい
      REAL(8):: WTR_MOD  ! 水分補正量 W/m2
      REAL(8):: GLA_MOD  ! 氷補正量 W/m2
      REAL(8):: H2O_MOD  ! 氷+水補正量 W/m2
      REAL(8):: FR_H2O  ! 氷+水 の対空隙比 0--1
      REAL(8):: FR_H2O_MOD  ! 氷+水 の対空隙比の補正量
      REAL(8):: FR_WTR_MOD
      REAL(8):: FR_GLA_MOD
      REAL(8):: LIQUID_RATIO  ! 氷+水補正を配分するための 水/(氷+水)
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
#endif
      REAL(8), PARAMETER :: UNITY = 1.0_8
!
#ifdef SIB_DEBUG
      INTEGER            :: IERR_COUNT
      REAL(8)            :: XERR_SUM
#endif
!
      LOGICAL, SAVE:: LFIRST = .TRUE.
      IF (LFIRST) THEN 
        WRITE(6, *) 'SOIL_MOD 2000-09-13 TOYODA'
        LFIRST = .FALSE. 
      ENDIF

      H2O_MOD_SOIL(:,:) = 0.0_8
      H2O_MOD_SNSL(:,:) = 0.0_8
      ENG_MOD_SOIL(:,:) = 0.0_8
      ENG_MOD_SNSL(:,:) = 0.0_8

!
! 動作フラグ
!
      BARESOIL_FLAG (1:NMAX) = (INFO_SNOW_NEXT(1: NMAX) <= 0)
      UNDERSNOW_FLAG(1:NMAX) = (INFO_SNOW_NEXT(1: NMAX) >= 0)

! ====================
! 負の水または氷を排除する処理
! ====================
! いちおう DO ループで書いておく。修正が多発したら考えよう。
!
! ----------
! > 無雪域 <
! ----------
!
#ifdef SIB_DEBUG
      IERR_COUNT = 0 
      XERR_SUM   = 0.
#endif      
!
      DO, K = 1, IDP
      DO, I = 1, NMAX
        IF (.NOT. BARESOIL_FLAG(I)) CYCLE   
!
        IF ( FR_WTR_SOIL(I,K) < 0.0_8 ) THEN
          WTR_MOD = -FR_WTR_SOIL(i, K) * FACT_WGT_H2O(i, K)
          FR_WTR_SOIL(i, K) = 0.0_8
          H2O_MOD_SOIL(i, K) = WTR_MOD
          ENG_MOD_SOIL(i, K) = WTR_MOD * CAP_WTR_C *
     1      (TMP_SOIL(i, K) - TMP_FREZ_c)
#ifdef SIB_DEBUG
          IERR_COUNT = IERR_COUNT + 1
          XERR_SUM   = XERR_SUM   + WTR_MOD**2
#endif
        ENDIF
!
        IF (FR_GLA_SOIL(i, K) < 0.0_8) THEN
          GLA_MOD = -FR_GLA_SOIL(i, K) * FACT_WGT_H2O(i, K)
          FR_GLA_SOIL(i, K) = 0.0_8
          H2O_MOD_SOIL(i, K) = H2O_MOD_SOIL(i, K) + GLA_MOD
          ENG_MOD_SOIL(i, K) = ENG_MOD_SOIL(i, K) + GLA_MOD *
     *      (-HL_W2G_C + CAP_GLA_C * (TMP_SOIL(i, K) - TMP_FREZ_c))
!    *      (-HL_W2G_C * CAP_GLA_C * (TMP_SOIL(i, K) - TMP_FREZ_c))
#ifdef SIB_DEBUG
          IERR_COUNT = IERR_COUNT + 1
          XERR_SUM   = XERR_SUM   + GLA_MOD**2
#endif
        ENDIF
      ENDDO
      ENDDO

! ----------
! > 有雪域 <
! ----------
! 
      DO, K = 1, IDP
      DO, i = 1, NMAX
        IF (.NOT. UNDERSNOW_FLAG(i)) CYCLE
!
        IF (FR_WTR_SNSL(i,k) < 0.0_8) THEN
          WTR_MOD = -FR_WTR_SNSL(i,k) * FACT_WGT_H2O(i,k)
          FR_WTR_SNSL (i,k) = 0.0_8
          H2O_MOD_SNSL(i,k) = WTR_MOD
          ENG_MOD_SNSL(i,k) = WTR_MOD * CAP_WTR_C *
     *      (TMP_SNSL(i,k) - TMP_FREZ_c)
#ifdef SIB_DEBUG
          IERR_COUNT = IERR_COUNT + 1
          XERR_SUM   = XERR_SUM   + WTR_MOD**2
#endif
        ENDIF
!
        IF (FR_GLA_SNSL(i,k) < 0.0_8) THEN
          GLA_MOD = -FR_GLA_SNSL(i,k) * FACT_WGT_H2O(i,k)
          FR_GLA_SNSL(i,k) = 0.0_8
          H2O_MOD_SNSL(i,k) = H2O_MOD_SNSL(i,k) + GLA_MOD
          ENG_MOD_SNSL(i,k) = ENG_MOD_SNSL(i,k) + GLA_MOD *
     1      (-HL_W2G_C + CAP_GLA_C * (TMP_SNSL(i,k) - TMP_FREZ_c))
!    1      (-HL_W2G_C * CAP_GLA_C * (TMP_SNSL(i,k) - TMP_FREZ_c))
#ifdef SIB_DEBUG
          IERR_COUNT = IERR_COUNT + 1
          XERR_SUM   = XERR_SUM   + GLA_MOD**2
#endif
        ENDIF
      ENDDO
      ENDDO
!
#ifdef SIB_DEBUG
      IF ( IERR_COUNT .GT. 0 ) THEN
        WRITE(6,*) 'SOIL_MOD H2O IS LACKED' , IERR_COUNT , XERR_SUM   
      ENDIF
#endif

! ===================
! H2O 総量規制 (上限)
! ===================

! 無雪域
!
      DO, K = 1, IDP
      DO, i = 1, NMAX

        IF (.NOT. BARESOIL_FLAG(i)) CYCLE
!
        FR_H2O = FR_WTR_SOIL(i, K) + FR_GLA_SOIL(i, K)
        IF (FR_H2O <= UNITY) CYCLE
! 補正開始
        FR_H2O_MOD = UNITY - FR_H2O
        H2O_MOD = FACT_WGT_H2O(i, K) * FR_H2O_MOD
        H2O_MOD_SOIL(i, K) = H2O_MOD_SOIL(i, K) + H2O_MOD
! 注。FR_H2O は大きすぎるくらいだからゼロ割は起こさない
        LIQUID_RATIO = FR_WTR_SOIL(i, K) / FR_H2O
        FR_WTR_MOD = FR_H2O_MOD * LIQUID_RATIO
        FR_GLA_MOD = FR_H2O_MOD * (UNITY - LIQUID_RATIO)

        GLA_MOD = FR_GLA_MOD * FACT_WGT_H2O(i, K)                ! 000915
!
        ENG_MOD_SOIL(i, K) = ENG_MOD_SOIL(i, K) +
     +    (TMP_SOIL(i, K) - TMP_FREZ_c) * FACT_WGT_H2O(i, K) *
     *    (FR_WTR_MOD * CAP_WTR_C + FR_GLA_MOD * CAP_GLA_C) -
     -    GLA_MOD * HL_W2G_C
        FR_WTR_SOIL(i, K) = FR_WTR_SOIL(i, K) + FR_WTR_MOD
        FR_GLA_SOIL(i, K) = FR_GLA_SOIL(i, K) + FR_GLA_MOD
#ifdef SIB_DEBUG
        write(6, *) 'soil_mod: excessive soil H2O reduced i=', i,
     ,    ' k=', k, ' fr_wtr_mod=', fr_wtr_mod,
     ,    ' fr_gla_mod=', fr_gla_mod , 
     1    'ENGMOD = ' ,   ENG_MOD_SOIL(i,K) 
#endif
      ENDDO
      ENDDO

      ! 有雪域
      DO, K = 1, IDP
      DO, i = 1, NMAX
        ! 補正すべきでない場合は次格子へ
        IF (.NOT. UNDERSNOW_FLAG(i)) CYCLE
        FR_H2O = FR_WTR_SNSL(i,k) + FR_GLA_SNSL(i,k)
        IF (FR_H2O <= UNITY) CYCLE
        ! 補正開始
        FR_H2O_MOD = UNITY - FR_H2O
        H2O_MOD = FACT_WGT_H2O(i,k) * FR_H2O_MOD
        H2O_MOD_SNSL(i,k) = H2O_MOD_SNSL(i,k) + H2O_MOD
        ! FR_H2O は大きすぎるくらいだからゼロ割は起こさない
        LIQUID_RATIO = FR_WTR_SNSL(i,k) / FR_H2O
        FR_WTR_MOD = FR_H2O_MOD * LIQUID_RATIO
        FR_GLA_MOD = FR_H2O_MOD * (UNITY - LIQUID_RATIO)
!
        GLA_MOD    = FR_GLA_MOD * FACT_WGT_H2O(i, k)              ! 000915
!
        ENG_MOD_SNSL(i,k) = ENG_MOD_SNSL(i,k) +
     +    (TMP_SNSL(i,k) - TMP_FREZ_c) * FACT_WGT_H2O(i, k) *
     *    (FR_WTR_MOD * CAP_WTR_C + FR_GLA_MOD * CAP_GLA_C) -
     -    GLA_MOD * HL_W2G_C
        FR_WTR_SNSL(i,k) = FR_WTR_SNSL(i,k) + FR_WTR_MOD
        FR_GLA_SNSL(i,k) = FR_GLA_SNSL(i,k) + FR_GLA_MOD
#ifdef SIB_DEBUG
        write(6, *) 'soil_mod: excessive snsl H2O reduced i=', i,
     ,    ' fr_wtr_mod=', fr_wtr_mod, ' fr_gla_mod=', fr_gla_mod,
     ,    ' fr_wtr=', fr_wtr_snsl(i,k), ' fr_gla=', fr_gla_snsl(i,k) ,
     1    'H2O_MOD = ' ,   H2O_MOD_SNSL(i,k) ,
     1    'ENG_MOD = ' ,   ENG_MOD_SNSL(i,k) ,
     1    (TMP_SNSL(i,k) - TMP_FREZ_c) , FACT_WGT_H2O(i,k) , 
     1     FR_WTR_MOD * CAP_WTR_C , FR_GLA_MOD * CAP_GLA_C , 
     1     GLA_MOD * HL_W2G_C

#endif
      ENDDO
      ENDDO

! ===================
! H2O 総量規制 (下限 ... 廃止されるかも)
! ===================
!     上限とのコーディング上の主な相違は水氷の分配比の
!     決定におけるゼロ割抑止にある。

      ! 無雪域
      DO, K = 1, IDP
      DO, i = 1, NMAX
        ! 補正すべきでない場合は次格子へ
        IF (.NOT. BARESOIL_FLAG(i)) CYCLE
        FR_H2O = FR_WTR_SOIL(i, K) + FR_GLA_SOIL(i, K)
        IF (FR_H2O >= FR_SOIL_H2O_MIN_c) CYCLE
        ! 補正開始
        FR_H2O_MOD = FR_SOIL_H2O_MIN_c - FR_H2O
        H2O_MOD = FACT_WGT_H2O(i, K) * FR_H2O_MOD
        H2O_MOD_SOIL(i, K) = H2O_MOD_SOIL(i, K) + H2O_MOD
        ! 分配比の決定
        IF (FR_H2O /= 0.0_8) THEN
          LIQUID_RATIO = FR_WTR_SOIL(i, K) / FR_H2O
        ELSE IF (TMP_SOIL(i, K) > TMP_FREZ_c) THEN
          LIQUID_RATIO = UNITY
        ELSE IF (TMP_SOIL(i, K) == TMP_FREZ_c) THEN
          LIQUID_RATIO = 0.5_8
        ELSE
          LIQUID_RATIO = 0.0_8
        ENDIF
!
        FR_WTR_MOD = FR_H2O_MOD * LIQUID_RATIO
        FR_GLA_MOD = FR_H2O_MOD * (UNITY - LIQUID_RATIO)
        GLA_MOD    = FR_GLA_MOD * FACT_WGT_H2O(i, k)             
        ENG_MOD_SOIL(i, K) = ENG_MOD_SOIL(i, K) +
     +    (TMP_SOIL(i, K) - TMP_FREZ_c) * FACT_WGT_H2O(i, K) *
     *    (FR_WTR_MOD * CAP_WTR_C + FR_GLA_MOD * CAP_GLA_C) -
     -    GLA_MOD * HL_W2G_C
        FR_WTR_SOIL(i, K) = FR_WTR_SOIL(i, K) + FR_WTR_MOD
        FR_GLA_SOIL(i, K) = FR_GLA_SOIL(i, K) + FR_GLA_MOD
#ifdef SIB_DEBUG
        write(6, *) 'soil_mod: soil H2O supplement i=', i,
     ,    ' k=', k,
     ,    ' fr_wtr_mod=', fr_wtr_mod, ' fr_gla_mod=', fr_gla_mod
#endif
      ENDDO
      ENDDO

      ! 有雪域
      DO, K = 1, IDP
      DO, i = 1, NMAX
        ! 補正すべきでない場合は次格子へ
        IF (.NOT. UNDERSNOW_FLAG(i)) CYCLE
        FR_H2O = FR_WTR_SNSL(i,k) + FR_GLA_SNSL(i,k)
        IF (FR_H2O >= FR_SOIL_H2O_MIN_c) CYCLE
        ! 補正開始
        FR_H2O_MOD = FR_SOIL_H2O_MIN_c - FR_H2O
        H2O_MOD = FACT_WGT_H2O(i,k) * FR_H2O_MOD
        H2O_MOD_SNSL(i,k) = H2O_MOD_SNSL(i,k) + H2O_MOD
        ! 分配比の決定
        IF (FR_H2O /= 0.0_8) THEN
          LIQUID_RATIO = FR_WTR_SNSL(i,k) / FR_H2O
        ELSE IF (TMP_SOIL(i,k) > TMP_FREZ_c) THEN
          LIQUID_RATIO = UNITY
        ELSE IF (TMP_SOIL(i,k) == TMP_FREZ_c) THEN
          LIQUID_RATIO = 0.5_8
        ELSE
          LIQUID_RATIO = 0.0_8
        ENDIF
        !
        FR_WTR_MOD = FR_H2O_MOD * LIQUID_RATIO
        FR_GLA_MOD = FR_H2O_MOD * (UNITY - LIQUID_RATIO)
        GLA_MOD    = FR_GLA_MOD * FACT_WGT_H2O(i, k)             
        ENG_MOD_SNSL(i,k) = ENG_MOD_SNSL(i,k) +
     +    (TMP_SNSL(i,k) - TMP_FREZ_c) * FACT_WGT_H2O(i,k) *
     *    (FR_WTR_MOD * CAP_WTR_C + FR_GLA_MOD * CAP_GLA_C) -
     -    GLA_MOD * HL_W2G_C
        FR_WTR_SNSL(i,k) = FR_WTR_SNSL(i,k) + FR_WTR_MOD
        FR_GLA_SNSL(i,k) = FR_GLA_SNSL(i,k) + FR_GLA_MOD
#ifdef SIB_DEBUG
        write(6, *) 'soil_mod: snsl H2O supplement ik=', i,k,
     ,    ' fr_wtr_mod=', fr_wtr_mod, ' fr_gla_mod=', fr_gla_mod
#endif
      ENDDO
      ENDDO

! ================
! >> 詳細モニタ <<
! ================

#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! 土壌水の補正量 ver2 
!
        CMARK = 'WADLN1'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SOIL(1,1), ICN, UNITY, UNITY) 
        CMARK = 'WADLS1'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SNSL(1,1), ICS, UNITY, UNITY) 
!
        CMARK = 'WADLN2'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SOIL(1,2), ICN, UNITY, UNITY)     
        CMARK = 'WADLS2'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SNSL(1,2), ICS, UNITY, UNITY)     
!
        CMARK = 'WADLN3'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SOIL(1,3), ICN, UNITY, UNITY)
        CMARK = 'WADLS3'
        CALL MONIT_SIB_NEXT(CMARK, H2O_MOD_SNSL(1,3), ICS, UNITY, UNITY)
!
! エネルギーの補正量 
!
        CMARK = 'HADLN1C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SOIL(1,1), ICN, UNITY, UNITY)
        CMARK = 'HADLS1C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SNSL(1,1), ICS, UNITY, UNITY)
!
        CMARK = 'HADLN2C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SOIL(1,2), ICN, UNITY, UNITY)
        CMARK = 'HADLS2C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SNSL(1,2), ICS, UNITY, UNITY)
!
        CMARK = 'HADLN3C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SOIL(1,3), ICN, UNITY, UNITY)
        CMARK = 'HADLS3C' 
        CALL MONIT_SIB_NEXT(CMARK, ENG_MOD_SNSL(1,3), ICS, UNITY, UNITY)
      ENDIF
!
#endif

      END SUBROUTINE SOIL_MOD
!===============================================================
      SUBROUTINE MR_SOIL_MOD
!
! soil_mod で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
      CHARACTER(7):: CMARK 
      CHARACTER(32):: CTITLE 
      CHARACTER(13):: CUNIT
!
      LOGICAL, SAVE::  LFIRST 
!
!
      DATA LFIRST /.TRUE./
!
      IF (LFIRST) THEN
        WRITE(6, *) 'MONIT_REGIST_SOIL_MOD 99/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!  水の補正量   
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN   
!
      CMARK = 'WADLN1'
      CTITLE = 'H2O MODIFY IN SOIL 1 NOS (SOILMOD)'
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'WADLS1'
      CTITLE = 'H2O MODIFY IN SOIL 1 SNW (SOILMOD)'
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'WADLN2'
      CTITLE = 'H2O MODIFY IN SOIL 2 NOS (SOILMOD) '
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'WADLS2'
      CTITLE = 'H2O MODIFY IN SOIL 2 SNW (SOILMOD) '
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'WADLN3'
      CTITLE = 'H2O MODIFY IN SOIL 3 NOS (SOILMOD) '
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'WADLS3'
      CTITLE = 'H2O MODIFY IN SOIL 3 SNW (SOILMOD) '
      CUNIT = 'KG/M**2/S'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
! エネルギー補正
!
      CMARK = 'HADLN1C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 1 NOS (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'HADLS1C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 1 SNW (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'HADLN2C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 2 NOS (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'HADLS2C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 2 SNW (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'HADLN3C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 3 NOS (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      CMARK = 'HADLS3C'
      CTITLE = 'ENG WITH H2O MODIFY IN SOIL 3 SNW (SOILMOD) '
      CUNIT = 'W/M**2'
      CALL monit_regist_sib(CMARK, CTITLE, CUNIT)
!
      ENDIF
#endif
!
      RETURN
      END SUBROUTINE MR_SOIL_MOD

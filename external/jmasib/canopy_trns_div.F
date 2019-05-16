! 蒸散の、土壌層毎の分配を計算する。
! 99/07/25 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE CANOPY_TRNS_DIV (
     I   NSTR         , NEND        , INSFLAG    , 
     I   TRNS_TOT     , CVR_S_or_N  , DELT_CNP   ,  
     I   FAC_TRNS_A   , FAC_TRNS_B  , 
!    I   KLEV_TRNS_MAX ,
#ifdef MONYOS
     I   ITR          , ITR_MAX     , 
     I   TMP_SOIL     , 
#endif
     M   TRNS_BUCKET )
!
! 注意
!   蒸散フラックスの向きは上から下を正とするので、蒸散の値は負。
!
! 課題。
!   割り振り方の問題で値がわずかにだと思うが、正になる場合あり。
!   その場合、最も蒸散が多い層に押しつけている。
!
      USE SIBPRM , ONLY :
     1   IDX     , 
     1   ICG     , 
     1   IDP
!
      USE SIBCON , ONLY : 
     1   TMP_FREZ_c  ,
     1   CAP_WTR_c   
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN,ICS ,
     1     IMONIT_LEVEL        , 
     1     IMONIT_ALL          , IMONIT_ALMA 
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG
!
      REAL(8),INTENT(IN)    :: DELT_CNP 
      REAL(8),INTENT(IN)    :: CVR_S_or_N   (IDX)
      REAL(8),INTENT(IN)    :: TRNS_TOT     (IDX,ICG)       ! 蒸散フラックス
      REAL(8),INTENT(INOUT) :: FAC_TRNS_A   (IDX,IDP,ICG)
      REAL(8),INTENT(INOUT) :: FAC_TRNS_B   (IDX,IDP,ICG)
!
#ifdef MONYOS
      INTEGER,INTENT(IN)    :: ITR
      INTEGER,INTENT(IN)    :: ITR_MAX
      REAL(8),INTENT(IN)   :: TMP_SOIL (IDX,IDP) ! 地面第一層温度
#endif

!
! INOUT
!
      REAL(8),INTENT(INOUT) :: TRNS_BUCKET (IDX,IDP) ! 蒸散吸い上げ量
!                               土壌水温は積算中一定なので、量だけモニタする
!
! WORK
!
      REAL(8) :: TRNS (IDX,IDP,ICG) 
      INTEGER :: KMAX_TRNS(IDX,ICG) 
      REAL(8) :: TRNS_MAX (IDX,ICG) 
      REAL(8) :: TRNS_SUM (IDX,ICG) 
!
      INTEGER :: I
      INTEGER :: JDP 
      INTEGER :: JCG
#ifdef SIB_DEBUG
      REAL(8) :: DIFF
#endif
!
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      REAL(8)           :: WORK(IDX) 
      CHARACTER(7) :: CMARK 
#endif
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_TRNS_DIV 99/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
! ==========================
! >> 蒸散が零の場合の処理 <<
! ==========================
!
      DO JCG = 1, 2-INSFLAG
      DO I   = NSTR,NEND
!
#ifdef SIB_DEBUG
        IF  ( TRNS_TOT(I,JCG) .GT. 0.D0 ) THEN
          write(6,*) 'canopy_trns_div error trns_tot>0 ' , i , jcg , 
     1               trns_tot(i,jcg)
          stop 999
        ENDIF
!
        DIFF =   FAC_TRNS_B(I,1,JCG) 
     1         + FAC_TRNS_B(I,2,JCG) 
     1         + FAC_TRNS_B(I,3,JCG) 
        IF  ( ABS( DIFF - 1.D0 ) .GE. 1.D-1 ) THEN
          IF ( .NOT. 
     1               (        ABS(DIFF).LE.1.D-1 
     1                  .AND. ABS(TRNS_TOT(I,JCG)).LE.1.D-10 )
     1       ) THEN
            write(6,*) 'canopy_trns_div error : ',
     1                 ' sum_fact_b.NE.1 and trns>0',
     1                i , jcg , 
     1                fac_trns_b(i,1,jcg)   ,
     1                fac_trns_b(i,2,jcg)   ,
     1                fac_trns_b(i,3,jcg)   ,
     1                TRNS_TOT(I,JCG)    
            stop 999
          ENDIF
        ENDIF 

#endif
!
        IF  ( TRNS_TOT(I,JCG) .EQ. 0.D0 ) THEN
!         TRNS_TOT  (I,  JCG) = 0.D0
          FAC_TRNS_A(I,1,JCG) = 0
          FAC_TRNS_A(I,2,JCG) = 0
          FAC_TRNS_A(I,3,JCG) = 0
        ENDIF
      ENDDO 
      ENDDO 
!
! ====================
! >> 各層ごとの分配 <<
! ====================
!
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        TRNS(I,JDP,JCG) 
     1     =  - FAC_TRNS_A(I,JDP,JCG) 
     1        + FAC_TRNS_B(I,JDP,JCG) * TRNS_TOT(I,JCG) 
      ENDDO
      ENDDO
      ENDDO
!
! ================
! >> エラー処理 <<
! ================
!
#ifdef SIB_DEBUG
!
!  とりあえず蒸散が正の層があったら警告
!
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        if ( trns(i,jdp,jcg) .gt.0 ) then
          write(6,*) 'canopy_trns_div warning : tr > 0 : ' , 
     1               ' It will be modified.' ,  
     1            'I' , i , 'IDP' , jdp , 'ICG' , jcg  , 
     1            'TOTAL' , trns_tot(i,jcg),
     1            '1st'   , trns(i,1,jcg) ,
     1            '2nd'   , trns(i,2,jcg) ,
     1            '3rd'   , trns(i,3,jcg) 
        endif
!
      ENDDO
      ENDDO
      ENDDO
#endif
!
! -------------------------------------------------------
! > sum TRNS と TRNS_TOT の差を一番大きい層に押しつける <
! -------------------------------------------------------
!
      TRNS_MAX  (:,:) = 1.D40
      KMAX_TRNS (:,:) = 0
      TRNS_SUM  (:,:) = 0.D0
!
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        IF ( TRNS(I,JDP,JCG) .LT. TRNS_MAX(I,JCG) ) THEN
          KMAX_TRNS (I,JCG) = JDP 
          TRNS_MAX  (I,JCG) = TRNS(I,JDP,JCG)
	  ELSE
	    TRNS(I,JDP,JCG) = 0.
	    KMAX_TRNS (I,JCG) = JDP 
          TRNS_MAX  (I,JCG) = TRNS(I,JDP,JCG)
	  ENDIF
        TRNS_SUM(I,JCG) = TRNS_SUM(I,JCG) + TRNS(I,JDP,JCG)   
      ENDDO
      ENDDO
      ENDDO
!
      DO JCG=1, 2-INSFLAG
      DO I=NSTR,NEND 
        TRNS(I,KMAX_TRNS(I,JCG),JCG)
     1    = TRNS(I,KMAX_TRNS(I,JCG),JCG) 
     1      + ( TRNS_TOT(I,JCG) - TRNS_SUM(I,JCG) ) 
      ENDDO
      ENDDO
!
! --------------------------------------------
! > 正の層の分を負の層に、値に応じて比例配分 <
! --------------------------------------------
!
      TRNS_SUM     (:,:) = 0.D0
!
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        IF ( TRNS(I,JDP,JCG) .LT. 0. ) THEN
          TRNS_SUM(I,JCG) = TRNS_SUM(I,JCG) + TRNS(I,JDP,JCG)   
        ELSE
          TRNS(I,JDP,JCG) = 0. 
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        IF   ( TRNS(I,JDP,JCG) .LT. 0. ) THEN
          TRNS(I,JDP,JCG) = TRNS(I,JDP,JCG) 
     1                      * ( TRNS_TOT(I,JCG) / TRNS_SUM(I,JCG) )
        ENDIF
      ENDDO
      ENDDO
      ENDDO
!
#ifdef SIB_DEBUG
      DO JCG = 1, 2-INSFLAG
      DO JDP = 1, IDP
      DO I   = NSTR,NEND
        if ( trns(i,jdp,jcg) .gt.0 ) then
          write(6,*) 'canopy_trns_div error i dep cg ' , 
     1            i , jdp , jcg  , trns_tot(i,jcg),
     1            KMAX_TRNS(I,JCG) , 
     1            trns(i,1,jcg) ,
     1            trns(i,2,jcg) ,
     1            trns(i,3,jcg) 
!         stop 999
        endif
      ENDDO
      ENDDO
      ENDDO
#endif
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
      IF ( ITR.EQ.ITR_MAX) THEN

      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
        IF ( INSFLAG .EQ. 0 ) THEN 
          CMARK = 'WTRBNL1' 
          CALL MONIT_SIB ( CMARK , TRNS(1,1,1),ICN, DELT_CNP , ONE )
          CMARK = 'WTRBNL2' 
          CALL MONIT_SIB ( CMARK , TRNS(1,2,1),ICN, DELT_CNP , ONE )
          CMARK = 'WTRBNL3' 
          CALL MONIT_SIB ( CMARK , TRNS(1,3,1),ICN, DELT_CNP , ONE )
          CMARK = 'WTRBGL1' 
          CALL MONIT_SIB ( CMARK , TRNS(1,1,2),ICN, DELT_CNP , ONE )
          CMARK = 'WTRBGL2' 
          CALL MONIT_SIB ( CMARK , TRNS(1,2,2),ICN, DELT_CNP , ONE )
          CMARK = 'WTRBGL3' 
          CALL MONIT_SIB ( CMARK , TRNS(1,3,2),ICN, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,1)-TMP_FREZ_c ) 
     1                          * TRNS(I,1,1)
          ENDDO
          CMARK = 'ITRBNL1' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,2)-TMP_FREZ_c ) 
     1                          * TRNS(I,2,1)
          ENDDO
          CMARK = 'ITRBNL2' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,3)-TMP_FREZ_c ) 
     1                          * TRNS(I,3,1)
          ENDDO
          CMARK = 'ITRBNL3' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,1)-TMP_FREZ_c ) 
     1                          * TRNS(I,1,2)
          ENDDO
!
          CMARK = 'ITRBGL1' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,2)-TMP_FREZ_c ) 
     1                          * TRNS(I,2,2)
          ENDDO
          CMARK = 'ITRBGL2' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,3)-TMP_FREZ_c ) 
     1                          * TRNS(I,3,2)
          ENDDO
          CMARK = 'ITRBGL3' 
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
! --------------
! >>  雪あり <<   :  TMP_SOIL = TMP_SNSL である。
! --------------
!
        ELSE
!
          CMARK = 'WTRBSL1' 
          CALL MONIT_SIB ( CMARK , TRNS(1,1,1),ICS, DELT_CNP , ONE )
          CMARK = 'WTRBSL2' 
          CALL MONIT_SIB ( CMARK , TRNS(1,2,1),ICS, DELT_CNP , ONE )
          CMARK = 'WTRBSL3' 
          CALL MONIT_SIB ( CMARK , TRNS(1,3,1),ICS, DELT_CNP , ONE )
!
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,1)-TMP_FREZ_c ) 
     1                          * TRNS(I,1,1)
          ENDDO
          CMARK = 'ITRBSL1' 
          CALL MONIT_SIB ( CMARK , WORK ,ICS, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,2)-TMP_FREZ_c ) 
     1                          * TRNS(I,2,1)
          ENDDO
          CMARK = 'ITRBSL2' 
          CALL MONIT_SIB ( CMARK , WORK ,ICS, DELT_CNP , ONE )
!
          DO I=NSTR,NEND
            WORK(I) = CAP_WTR_c * ( TMP_SOIL(I,3)-TMP_FREZ_c ) 
     1                          * TRNS(I,3,1)
          ENDDO
          CMARK = 'ITRBSL3' 
          CALL MONIT_SIB ( CMARK , WORK ,ICS, DELT_CNP , ONE )
        ENDIF
!
      ENDIF
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK = 'A_TRVEG'
        IF ( INSFLAG .EQ.0 ) THEN     ! 雪なし
          DO I=NSTR,NEND
            WORK(I) = ( TRNS(I,1,1) + TRNS(I,2,1) + TRNS(I,3,1) 
     1                + TRNS(I,1,2) + TRNS(I,2,2) + TRNS(I,3,2) )
          ENDDO
          CALL MONIT_SIB ( CMARK , WORK ,ICN, DELT_CNP , ONE )
!
        ELSEIF ( INSFLAG .EQ.1 ) THEN     ! 雪なし
          DO I=NSTR,NEND
            WORK(I) = ( TRNS(I,1,1) + TRNS(I,2,1) + TRNS(I,3,1) )
          ENDDO
          CALL MONIT_SIB ( CMARK , WORK ,ICS, DELT_CNP , ONE )
        ENDIF 
      ENDIF
#endif
!   
! ====================
! >> バケツへの格納 <<
! ====================
!
      DO 3000 JCG = 1, 2-INSFLAG
      DO 3000 JDP = 1, IDP
* ここで並列化すべし
      DO 3000 I = NSTR,NEND
        TRNS_BUCKET(I,JDP) 
     1      = TRNS_BUCKET(I,JDP) 
     1        + TRNS(I,JDP,JCG) * CVR_S_or_N(I) * DELT_CNP
 3000 CONTINUE
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_CANOPY_TRNS_DIV
!
! canopy_trns_div で出力するデバッグ用変数の登録作業
!
      USE SIBPRM ,  ONLY : 
     1   IDP 
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1     IMONIT_LEVEL    , 
     1     IMONIT_ALL      , IMONIT_ALMA
#endif
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      INTEGER      :: K 
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_CANOPY_TRNS_DIV 99/11 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
#include "noparallel_fixedform"
      DO K=1,IDP      
!
! 雪なしキャノピー
!
        CMARK  = 'WTRBNLx'
        CTITLE = 'TRNS WATER SOIL x VIA CANOPY NOS'
        CUNIT  = 'KG/M**2/S'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(17:17),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 雪ありキャノピー
!
        CMARK  = 'WTRBSLx'
        CTITLE = 'TRNS WATER SOIL x VIA CANOPY SNW'
        CUNIT  = 'KG/M**2/S'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(17:17),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
        CMARK  = 'WTRBGLx'
        CTITLE = 'TRNS WATER SOIL x VIA GRASS'
        CUNIT  = 'KG/M**2/S'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(17:17),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 雪なしキャノピー
!
        CMARK  = 'ITRBNLx'
        CTITLE = 'ENG WITH TRNS WATER SOIL x VIA CANOPY NOS'
        CUNIT  = 'W/M**2'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(26:26),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 雪ありキャノピー
!
        CMARK  = 'ITRBSLx'
        CTITLE = 'ENG WITH TRNS WATER SOIL x VIA CANOPY SNW'
        CUNIT  = 'W/M**2'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(26:26),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 下草
!
        CMARK  = 'ITRBGLx'
        CTITLE = 'ENG WITH TRNS WATER SOIL x VIA GRASS'
        CUNIT  = 'W/M**2'
        WRITE(CMARK ( 7: 7),'(I1)') K        
        WRITE(CTITLE(26:26),'(I1)') K        
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDDO
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_TRVEG'
        CTITLE = 'ALMA - TVeg Transipiration from Cnp + Grass' 
        CUNIT  = 'kg/m**2/s'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
#endif
!
      RETURN
      END

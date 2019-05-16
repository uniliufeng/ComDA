! 雪の水まわり計算後のあと処理
! 99/07/12 保坂征宏 第一版
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
    
      SUBROUTINE SNOW_WFLUX_POST ( 
     I    NMAX             , 
     I    H2O_SNOW_NEXT    , ENG_SNOW_NEXT    ,
#ifdef MONYOS 
     I    DELT_CNP         , DELT_SNOW        ,
#endif
     M    INFO_SNOW_NEXT   , 
     O    CVR_SNW_NEXT     , CVR_SNW_INV_NEXT , 
     O    CVR_NOS_NEXT     , CVR_NOS_INV_NEXT , 
     O    TMP_SNOW_NEXT    , WTR_SNOW_NEXT    , GLA_SNOW_NEXT   , 
     O    LIST_SNOW        , ICNT_SNOW )
!
! INFO_SNOW_NEXT 
!      INFO_SNOW_NEXT(I) 
!                          -2  雪なし少雪あり
!                          -1  雪無し
!                           0  01
!                           1  1 層 
!                           2  2 層 
!                           3  3 層 
!                           4  4 層 
!
      USE SIBPRM , ONLY :
     1    IDX   , 
     1    ISN 
!
      USE SIBCON , ONLY :
     1    HL_W2G_c        , 
     1    HL_W2G_INV_c    , 
     1    CAP_WTR_c     ,
     1    CAP_GLA_c     , 
     1    CAP_GLA_c_INV , 
     1    TMP_FREZ_c    , 
     1    H2O_SNOW_01_1_c   ,       ! 雪がこれ以下なら 01 
     1    H2O_SNOW_01_1_INV_c ,     ! 
     1    CVR_SNW_INV_MAX_c            ! CVR_SNW 逆数上限。ダミー。
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB, INA ,
     1                     IMONIT_LEVEL, IMONIT_ALL, 
     1                     IMONIT_MEAN_MRI ,
     1                     IMONIT_SNAP_NPD
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)   :: NMAX 
#ifdef MONYOS
      REAL(8),INTENT(IN)   :: DELT_CNP
      REAL(8),INTENT(IN)   :: DELT_SNOW
#endif
!
      REAL(8),INTENT(IN)   :: H2O_SNOW_NEXT (IDX,ISN)  
      REAL(8),INTENT(IN)   :: ENG_SNOW_NEXT (IDX,ISN)  
!

!     REAL(8),INTENT(IN)   :: CVR_SNW     (IDX)  
!
! OUTPUT
!
      INTEGER,INTENT(INOUT) :: INFO_SNOW_NEXT       (IDX)
!
      REAL(8),INTENT(OUT)  :: CVR_SNW_NEXT     (IDX)  
      REAL(8),INTENT(OUT)  :: CVR_NOS_NEXT     (IDX)  
      REAL(8),INTENT(OUT)  :: CVR_SNW_INV_NEXT (IDX)  
      REAL(8),INTENT(OUT)  :: CVR_NOS_INV_NEXT (IDX)  
!
      REAL(8),INTENT(OUT)  :: TMP_SNOW_NEXT (IDX,ISN)
      REAL(8),INTENT(OUT)  :: WTR_SNOW_NEXT (IDX,ISN)
      REAL(8),INTENT(OUT)  :: GLA_SNOW_NEXT (IDX,ISN)
!
      INTEGER,INTENT(OUT)  :: LIST_SNOW (IDX,0:4)
      INTEGER,INTENT(OUT)  :: ICNT_SNOW (    0:4)
!
! WORK
!
      INTEGER   :: I 
      INTEGER   :: K 
      REAL(8)   :: ENX
      REAL(8),PARAMETER :: ONE = 1.D0 
      LOGICAL,SAVE :: LFIRST 
!
#ifdef MONYOS
!     REAL(8),PARAMETER  :: ONE = 1.D0 
      REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_WFLUX_POST POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
!  ==============
!  >> 雪の更新 <<
!  ==============
!
      TMP_SNOW_NEXT(:,:) = 0. 
      GLA_SNOW_NEXT(:,:) = 0. 
      WTR_SNOW_NEXT(:,:) = 0. 
!
!  注意。
!    この時点で、INFO_SNOW_NEXT は
!    雪無し・少雪バケツ なら 負、01 は 1 である必要あり。
!
      DO K=1,ISN
      DO I=1,NMAX
        IF ( INFO_SNOW_NEXT(I).GE.K ) THEN   
          ENX = ENG_SNOW_NEXT(I,K) + HL_W2G_c * H2O_SNOW_NEXT(I,K) 
!
          IF ( ENX .GT. 0. ) THEN                             !  水あり
            TMP_SNOW_NEXT (I,K) = TMP_FREZ_c 
            GLA_SNOW_NEXT (I,K) = - ENG_SNOW_NEXT(I,K) * HL_W2G_INV_c 
            WTR_SNOW_NEXT (I,K) = H2O_SNOW_NEXT(I,K) 
     1                              - GLA_SNOW_NEXT(I,K)
          ELSE
            TMP_SNOW_NEXT(I,K) 
     1          =   TMP_FREZ_c 
     1            + ENX / H2O_SNOW_NEXT(I,K) * CAP_GLA_c_INV       
            GLA_SNOW_NEXT(I,K) = H2O_SNOW_NEXT(I,K)
          ENDIF
        ENDIF 
      ENDDO
      ENDDO
!
!  ======================
!  >> 雪の被覆率の更新 <<
!  ======================
!   01 のときはここで INFO_SNOW_NEXT=0 にする
!
      DO 3000 I=1,NMAX 
        IF ( INFO_SNOW_NEXT(I) .GE. 0       ) THEN  ! 雪無し
          IF ( H2O_SNOW_NEXT(I,1) .GE. H2O_SNOW_01_1_c ) THEN  ! 1 以上
            CVR_SNW_NEXT(I)     = 1. 
            CVR_NOS_NEXT(I)     = 0. 
            CVR_SNW_INV_NEXT(I) = 1. 
            CVR_NOS_INV_NEXT(I) = CVR_SNW_INV_MAX_c   
!
          ELSE                                            ! 01  
            CVR_SNW_NEXT(I) 
     1       = MIN ( H2O_SNOW_NEXT(I,1) * H2O_SNOW_01_1_INV_c 
     1               , ONE )
            CVR_NOS_NEXT(I) = 1. - CVR_SNW_NEXT(I) 
            CVR_SNW_INV_NEXT(I)  
     1         = MIN ( 1. / CVR_SNW_NEXT(I) , CVR_SNW_INV_MAX_c ) 
            CVR_NOS_INV_NEXT(I)  
     1         = MIN ( 1. / CVR_NOS_NEXT(I) , CVR_SNW_INV_MAX_c ) 
            INFO_SNOW_NEXT(I)   = 0
          ENDIF
!
        ELSEIF ( INFO_SNOW_NEXT(I) .LT. 0       ) THEN  ! 雪無し
          CVR_SNW_NEXT(I)     = 0. 
          CVR_NOS_NEXT(I)     = 1. 
          CVR_SNW_INV_NEXT(I) = CVR_SNW_INV_MAX_c    
          CVR_NOS_INV_NEXT(I) = 1. 
!
        ENDIF
 3000 CONTINUE  
!
!  --------------------------
!  > 詰め替え用テーブル作成 <    
!  --------------------------
!   きっと必要だよねえ... 一応
!
      ICNT_SNOW(0) = 0
      ICNT_SNOW(1) = 0
      ICNT_SNOW(2) = 0
      ICNT_SNOW(3) = 0
      ICNT_SNOW(4) = 0
      DO 5000 I=1,NMAX
        IF     ( INFO_SNOW_NEXT(I) .EQ. 0 ) THEN  ! 01
          ICNT_SNOW(0)  = ICNT_SNOW(0) + 1 
          LIST_SNOW(ICNT_SNOW(0),0) = I  
        ELSEIF ( INFO_SNOW_NEXT(I) .EQ. 1 ) THEN  ! 1 
          ICNT_SNOW(1)  = ICNT_SNOW(1) + 1 
          LIST_SNOW(ICNT_SNOW(1),1) = I  
        ELSEIF ( INFO_SNOW_NEXT(I) .EQ. 2 ) THEN  ! 2 
          ICNT_SNOW(2) = ICNT_SNOW(2) + 1 
          LIST_SNOW(ICNT_SNOW(2),2) = I  
        ELSEIF ( INFO_SNOW_NEXT(I) .EQ. 3 ) THEN  ! 3 
          ICNT_SNOW(3) = ICNT_SNOW(3) + 1 
          LIST_SNOW(ICNT_SNOW(3),3) = I  
        ELSEIF ( INFO_SNOW_NEXT(I) .EQ. 4 ) THEN  ! 4 
          ICNT_SNOW(4) = ICNT_SNOW(4) + 1 
          LIST_SNOW(ICNT_SNOW(4),4) = I  
#ifdef SIB_DEBUG
        ELSEIF ( 
     1                  INFO_SNOW_NEXT(I) .NE. -1    ! 0  
     1           .AND.  INFO_SNOW_NEXT(I) .NE. -2    ! 少雪バケツ
     1         ) THEN
          WRITE(6,*) 'SNOW_WFLUX_POST ERROR INFO_SNOW_NEXT ' , 
     1               I , INFO_SNOW_NEXT(I) 
          STOP 999
#endif
        ENDIF
 5000 CONTINUE  
!
! ================
! >> 詳細モニタ <<
! ================
!
#ifdef MONYOS
      IF (     IMONIT_LEVEL .GE. IMONIT_MEAN_MRI 
     1    .OR. IMONIT_LEVEL .GE. IMONIT_ALL      ) THEN
        CMARK = 'CVRSNWx'
        DO K = 2,4
          WORK(:) = 0.
          WRITE(CMARK(7:7),'(I1)') K 
          DO I=1,NMAX
            IF ( INFO_SNOW_NEXT(I) .EQ. K ) THEN
              WORK(I) = 1.
            ENDIF
          ENDDO
          CALL MONIT_SIB ( CMARK , WORK, INA , DELT_SNOW , ONE )         
        ENDDO
      ENDIF
!
      IF (     IMONIT_LEVEL .GE. IMONIT_SNAP_NPD 
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1 
     1   ) THEN
        DO K = 1,ISN
          CMARK = 'SnWex'
          WRITE(CMARK(5:5),'(I1)') K
          CALL MONIT_SIB ( CMARK , H2O_SNOW_NEXT(1,K),INA,DELT_CNP,ONE)        
        ENDDO
      ENDIF
#endif
!
      RETURN
      END
!============================================================
      SUBROUTINE MR_SNOW_WFLUX_POST
!
! SNOW_WFLUX_POST で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL, IMONIT_ALL, 
     1                       IMONIT_MEAN_MRI ,
     1                       IMONIT_SNAP_NPD 
#endif
!
      USE SIBPRM , ONLY : ISN 
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      INTEGER        :: K
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_SNOW_WFLUX_POST 00/08 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI 
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
      CMARK  = 'CVRSNW2'
      CTITLE = 'FRACTION OF EXISTING SNOW LEV2' 
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'CVRSNW3'
      CTITLE = 'FRACTION OF EXISTING SNOW LEV3' 
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'CVRSNW4'
      CTITLE = 'FRACTION OF EXISTING SNOW LEV4' 
      CUNIT  = '0-1'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF (   IMONIT_LEVEL .GE. IMONIT_SNAP_NPD ) THEN
        DO K = 1,ISN
          CMARK = 'SnWex'
          CTITLE = 'Snow Equivalent Water LEV x' 
          WRITE(CMARK ( 5: 5),'(I1)') K
          WRITE(CTITLE(27:27),'(I1)') K
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
        ENDDO
      ENDIF
#endif
!
      RETURN
      END

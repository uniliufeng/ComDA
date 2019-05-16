!  雪での温度の式をたてる。
!  99/06/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SNOW_HFLUX_IMPLICIT_EQSET (
     I  NMAX      , NMAX_01   , 
     I  NMAX_2    , NMAX_3    , NMAX_4    , 
     I  DELT_SNOW , 
     I  PHASE_SNOW       , 
     I  CAP_SNOW         , HDIFF_FLUX_H      , 
     I  HDF_SNOW_H       ,
     O  SIB_MTX_SNOW     , 
     M  HFLUX_CONV )

!  雪での温度の式をたてる。
!
!        S1    1
!        S2    2
!        S3    3
!        S4    4
! 
!   注意
!     ・上からのフラックスは与えられている。      
!
!     ・雪の層数と、TMP_SNSL 等の関係
!       雪 1 層の時    TMP_SNOW(I,2)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,1) = 雪下端と雪下土壌
!       雪 2 層の時    TMP_SNOW(I,3)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,2) = 雪下端と雪下土壌
!       雪 3 層の時    TMP_SNOW(I,4)   = TMP_SNSL(I)  
!                      HDF_SNOW_H(I,3) = 雪下端と雪下土壌
! 
!     ・T_{frez} に等しい場合はその層は温度変化しないとして評価する。
!       PHASE_SNOW   ! 相変化 on 0 off 1
!
!
      USE SIBPRM , ONLY : 
     1    IDX    , 
     1    ISN      
!
      USE SIBCON , ONLY :
     1    HL_W2G_c , 
     1    DELT_FAC_SNOW_c
!
      IMPLICIT NONE 
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NMAX
      INTEGER,INTENT(IN)    :: NMAX_01(3)
!     INTEGER,INTENT(IN)    :: NMAX_02(3)
!     INTEGER,INTENT(IN)    :: NMAX_1 (3)
      INTEGER,INTENT(IN)    :: NMAX_2 (3)
      INTEGER,INTENT(IN)    :: NMAX_3 (3)
      INTEGER,INTENT(IN)    :: NMAX_4 (3)
!
!     REAL(8),INTENT(IN)    :: TMP_SNOW    (IDX,ISN)
!     REAL(8),INTENT(IN)    :: TMP_SNSL    (IDX)
!
      REAL(8),INTENT(IN)    :: PHASE_SNOW  (IDX,ISN)
!
      REAL(8),INTENT(IN)    :: DELT_SNOW
!
      REAL(8),INTENT(IN)    :: CAP_SNOW    (IDX,ISN)
!
      REAL(8),INTENT(IN)    :: HDIFF_FLUX_H(IDX,0:ISN) 
!
      REAL(8),INTENT(IN)    :: HDF_SNOW_H  (IDX,ISN)   ! 雪中 1-2,2-3,3-4  
!
! OUTPUT 
!
      REAL(8),INTENT(OUT)   :: SIB_MTX_SNOW(IDX,ISN,ISN+1) 
      REAL(8),INTENT(INOUT) :: HFLUX_CONV(IDX,ISN) 
!
! WORK
      INTEGER               :: I
      REAL(8)               :: DELT_TMP 
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_HFLUX_IMPLICIT_EQSET POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
      DELT_TMP = DELT_SNOW * DELT_FAC_SNOW_c
!
      CALL RESET ( SIB_MTX_SNOW , IDX*ISN*(ISN+1) )
!
!  ---------------------------------------------------
!  雪第一層の式 
!  ---------------------------------------------------
!
      DO 1100 I=NMAX_01(2),NMAX
        SIB_MTX_SNOW(I,1,1)                       ! 雪第一層
     1     =   CAP_SNOW(I,1) 
     1       + HDF_SNOW_H(I,1) * DELT_TMP            ! 雪第二層 or 土壌と
     1         * PHASE_SNOW(I,1)
!
        SIB_MTX_SNOW(I,1,2)                       ! 雪第二層
     1     = - HDF_SNOW_H(I,1) * DELT_TMP 
     1         * PHASE_SNOW(I,2)                     ! 一層しかない場合は 零
!
        SIB_MTX_SNOW(I,1,5) =                     ! 強制 （単位時間当たり）
     1      + HFLUX_CONV  (I,1)                      ! 短波分
     1      + HDIFF_FLUX_H(I,0)
     1      - HDIFF_FLUX_H(I,1) 
 1100 CONTINUE         
!
!  ---------------------------------------------------
!  雪第二層の式 
!  ---------------------------------------------------
!
      DO 1200 I=NMAX_2(2),NMAX
!
        SIB_MTX_SNOW(I,2,1)                        ! 第一層
     1      =  - HDF_SNOW_H(I,1) * DELT_TMP             
     1              * PHASE_SNOW(I,1)
!
        SIB_MTX_SNOW(I,2,2)                        ! 第二層 
     1    =   CAP_SNOW(I,2) 
     1      + HDF_SNOW_H(I,1) * DELT_TMP               ! 第一層と
     1              * PHASE_SNOW(I,2)
     1      + HDF_SNOW_H(I,2) * DELT_TMP               ! 第三層 or 土壌と
     1              * PHASE_SNOW(I,2)
!
        SIB_MTX_SNOW(I,2,3)                         ! 第三層
     1      =  - HDF_SNOW_H(I,2) * DELT_TMP             
     1              * PHASE_SNOW(I,3)                  ! 二層しかない場合零
!
        SIB_MTX_SNOW(I,2,5) =                       ! 強制
     1      + HFLUX_CONV    (I,2)                      ! 短波分
     1      + HDIFF_FLUX_H  (I,1)
     1      - HDIFF_FLUX_H  (I,2)
 1200 CONTINUE         
!
!  ---------------------------------------------------
!  雪第三層の式 
!  ---------------------------------------------------
!
      DO 1300 I=NMAX_3(2),NMAX
!
        SIB_MTX_SNOW(I,3,2)                          ! 第二層
     1      =  - HDF_SNOW_H(I,2) * DELT_TMP             
     1              * PHASE_SNOW(I,2)
!
        SIB_MTX_SNOW(I,3,3)                          ! 第三層 
     1    =   CAP_SNOW(I,3) 
     1      + HDF_SNOW_H(I,2) * DELT_TMP                ! 第二層と
     1              * PHASE_SNOW(I,3)
     1      + HDF_SNOW_H(I,3) * DELT_TMP                ! 第四層 or 土壌と
     1              * PHASE_SNOW(I,3)
!
        SIB_MTX_SNOW(I,3,4)                          ! 第四層
     1      =  - HDF_SNOW_H(I,3) * DELT_TMP             
     1              * PHASE_SNOW(I,4)                   ! 三層しかない場合零
!
        SIB_MTX_SNOW(I,3,5) =                        ! 強制
     1      + HFLUX_CONV    (I,3)                       ! 短波分
     1      + HDIFF_FLUX_H  (I,2)
     1      - HDIFF_FLUX_H  (I,3)
 1300 CONTINUE         
!
!  ---------------------------------------------------
!  雪第四層の式 
!  ---------------------------------------------------
!
      DO 1400 I=NMAX_4(2),NMAX
!
        SIB_MTX_SNOW(I,4,3)                           ! 第三層
     1      =  - HDF_SNOW_H(I,3) * DELT_TMP             
     1              * PHASE_SNOW(I,3)
!
        SIB_MTX_SNOW(I,4,4)                           ! 第四層 
     1    =   CAP_SNOW(I,4) 
     1      + HDF_SNOW_H(I,3) * DELT_TMP                   ! 第三層と
     1              * PHASE_SNOW(I,4)
     1      + HDF_SNOW_H(I,4) * DELT_TMP                   ! 地面と
     1              * PHASE_SNOW(I,4)
!
        SIB_MTX_SNOW(I,4,5) =                         ! 強制
     1      + HFLUX_CONV    (I,4)                          ! 短波分
     1      + HDIFF_FLUX_H  (I,3)
     1      - HDIFF_FLUX_H  (I,4)
 1400 CONTINUE         
!
      RETURN
      END
!===============================================================
! 雪内の熱の方程式を解く
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス

      SUBROUTINE SNOW_HFLUX_IMPLICIT_EQSLV (
     I  NMAX_01   , 
     I  NMAX_1    , NMAX_2    , NMAX_3    , NMAX_4    , 
     M  SMTX      ) 
!
!     S1  1
!     S2  2
!     S3  3
!     S4  4
!
      USE SIBPRM , ONLY : 
     1    IDX      
!
      IMPLICIT NONE 
!
! INPUT
      INTEGER,INTENT(IN)    :: NMAX_01(3)
      INTEGER,INTENT(IN)    :: NMAX_1 (3)
      INTEGER,INTENT(IN)    :: NMAX_2 (3)
      INTEGER,INTENT(IN)    :: NMAX_3 (3)
      INTEGER,INTENT(IN)    :: NMAX_4 (3)
!
! INOUT
!
      REAL(8),INTENT(INOUT) :: SMTX(IDX,4,5) 
!
! WORK
      INTEGER               :: I
      REAL(8)               :: A22I 
      REAL(8)               :: A33I 
      REAL(8)               :: A44I 
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_HFLUX_IMPLICIT_EQSLV POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF

!
! -----------
! >> 一層 << 
! -----------
!
      DO 1000 I=NMAX_01(2),NMAX_1(3)
        SMTX(I,1,5) = SMTX(I,1,5) / SMTX(I,1,1) 
!       SMTX(I,1,1) = 1.
 1000 CONTINUE
!
! -----------
! >> 二層 << 
! -----------
!
      DO 2000 I=NMAX_2(2),NMAX_2(3)
!
!a11  a12           a15
!a21  a22           a25
!
        A22I = 1. / SMTX(I,2,2) 
        SMTX(I,2,1) = SMTX(I,2,1) * A22I  
        SMTX(I,2,5) = SMTX(I,2,5) * A22I  
!       SMTX(I,2,2) = 1. 
!
        SMTX(I,1,1) = SMTX(I,1,1) - SMTX(I,1,2) * SMTX(I,2,1) 
        SMTX(I,1,5) = SMTX(I,1,5) - SMTX(I,1,2) * SMTX(I,2,5) 
!       SMTX(I,1,2) = 0. 
!
!a11*  0*           a15*
!a21*  1*           a25*
!
        SMTX(I,1,5) = SMTX(I,1,5) / SMTX(I,1,1)
!       SMTX(I,1,1) = 1. 
!
! 1*                a15*
!a21   1            a25
!
        SMTX(I,2,5) = SMTX(I,2,5) - SMTX(I,2,1) * SMTX(I,1,5)  
!       SMTX(I,2,1) = 0. 
 2000 CONTINUE
!
! -----------
! >> 三層 << 
! -----------
!
      DO 3000 I=NMAX_3(2),NMAX_3(3)
!
!a11  a12           a15
!a21  a22  a23      a25
!     a32  a33      a35
!
        A33I = 1. / SMTX(I,3,3) 
        SMTX(I,3,2) = SMTX(I,3,2) * A33I  
        SMTX(I,3,5) = SMTX(I,3,5) * A33I  
!       SMTX(I,3,3) = 1. 
!
        SMTX(I,2,2) = SMTX(I,2,2) - SMTX(I,2,3) * SMTX(I,3,2) 
        SMTX(I,2,5) = SMTX(I,2,5) - SMTX(I,2,3) * SMTX(I,3,5) 
!       SMTX(I,2,3) = 0. 
!
!a11  a12           a15
!a21  a22*  0*      a25*
!     a32*  1*      a35*
!
        A22I = 1. / SMTX(I,2,2) 
        SMTX(I,2,1) = SMTX(I,2,1) * A22I  
        SMTX(I,2,5) = SMTX(I,2,5) * A22I  
!       SMTX(I,2,2) = 1. 
!
        SMTX(I,1,1) = SMTX(I,1,1) - SMTX(I,1,2) * SMTX(I,2,1) 
        SMTX(I,1,5) = SMTX(I,1,5) - SMTX(I,1,2) * SMTX(I,2,5) 
!       SMTX(I,1,2) = 0. 
!
!a11*  0*           a15*
!a21*  1*           a25*
!     a32   1       a35
!
        SMTX(I,1,5) = SMTX(I,1,5) / SMTX(I,1,1)
!       SMTX(I,1,1) = 1. 
!
! 1*                a15*
!a21   1            a25
!     a32   1       a35
!
        SMTX(I,2,5) = SMTX(I,2,5) - SMTX(I,2,1) * SMTX(I,1,5)  
!       SMTX(I,2,1) = 0. 
        SMTX(I,3,5) = SMTX(I,3,5) - SMTX(I,3,2) * SMTX(I,2,5)  
!       SMTX(I,3,2) = 0. 
 3000 CONTINUE
!
! -----------
! >> 四層 << 
! -----------
!
      DO 4000 I=NMAX_4(2),NMAX_4(3)
!
!a11  a12            a15
!a21  a22  a23       a25
!     a32  a33  a34  a35
!          a43  a44  a45
!
        A44I = 1. / SMTX(I,4,4) 
        SMTX(I,4,3) = SMTX(I,4,3) * A44I  
        SMTX(I,4,5) = SMTX(I,4,5) * A44I  
!       SMTX(I,4,4) = 1. 
!
        SMTX(I,3,3) = SMTX(I,3,3) - SMTX(I,3,4) * SMTX(I,4,3) 
        SMTX(I,3,5) = SMTX(I,3,5) - SMTX(I,3,4) * SMTX(I,4,5) 
!       SMTX(I,3,4) = 0. 
!
!a11  a12            a15
!a21  a22  a23       a25
!     a32  a33*  0   a35*
!          a43*  1   a45*
!
        A33I = 1. / SMTX(I,3,3) 
        SMTX(I,3,2) = SMTX(I,3,2) * A33I  
        SMTX(I,3,5) = SMTX(I,3,5) * A33I  
!       SMTX(I,3,3) = 1. 
!
        SMTX(I,2,2) = SMTX(I,2,2) - SMTX(I,2,3) * SMTX(I,3,2) 
        SMTX(I,2,5) = SMTX(I,2,5) - SMTX(I,2,3) * SMTX(I,3,5) 
!       SMTX(I,2,3) = 0. 

!
!a11  a12            a15
!a21  a22*  0*       a25*
!     a32*  1*       a35*
!          a43   1   a45
!
        A22I = 1. / SMTX(I,2,2) 
        SMTX(I,2,1) = SMTX(I,2,1) * A22I  
        SMTX(I,2,5) = SMTX(I,2,5) * A22I  
!       SMTX(I,2,2) = 1. 
!
        SMTX(I,1,1) = SMTX(I,1,1) - SMTX(I,1,2) * SMTX(I,2,1) 
        SMTX(I,1,5) = SMTX(I,1,5) - SMTX(I,1,2) * SMTX(I,2,5) 
!       SMTX(I,1,2) = 0. 
!
!a11*  0*            a15*
!a21*  1*            a25*
!     a32   1        a35
!          a43   1   a45
!
        SMTX(I,1,5) = SMTX(I,1,5) / SMTX(I,1,1)
!       SMTX(I,1,1) = 1. 
!
! 1*                 a15*
!a21   1             a25
!     a32   1        a35
!          a43   1   a45
!
        SMTX(I,2,5) = SMTX(I,2,5) - SMTX(I,2,1) * SMTX(I,1,5)  
!       SMTX(I,2,1) = 0. 
        SMTX(I,3,5) = SMTX(I,3,5) - SMTX(I,3,2) * SMTX(I,2,5)  
!       SMTX(I,3,2) = 0. 
        SMTX(I,4,5) = SMTX(I,4,5) - SMTX(I,4,3) * SMTX(I,3,5)  
!       SMTX(I,4,3) = 0. 
 4000 CONTINUE
!
      RETURN
      END
!================================================================
! 雪の各層の得る熱フラックス収束の計算 （陽的評価に、陰的評価補正分を加える）
! 土壌への熱フラックスの計算
! 99/07/11 保坂征宏 第一版
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SNOW_HFLUX_IMPLICIT_CORR (
     I  NMAX      , NMAX_01   , NMAX_1    , 
     I  NMAX_2    , NMAX_3    , NMAX_4    , 
     I  SIB_MTX_SNOW          , 
     I  HDF_SNOW_H            , 
     I  DELT_SNOW             , PHASE_SNOW           , 
     I  CVR_SNW              ,               
     I  SRAD_BUCKET_SNOW_TO_SNSL     ,  
     M  HDIFF_FLUX_H          , 
     O  ENG_BUCKET_SNOW_TO_SNSL      , 
     M  HFLUX_CONV             )  
!
      USE SIBPRM , ONLY :
     1   IDX  , 
     1   ISN    
!
      USE SIBCON , ONLY :
     1   DELT_FAC_SNOW_c
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICS ,
     1                IMONIT_LEVEL, IMONIT_ALL ,
     1                IMONIT_FLUX_COMMON 
#endif
!
      IMPLICIT NONE
!
! INPUT 
!
      INTEGER,INTENT(IN)   :: NMAX
      INTEGER,INTENT(IN)   :: NMAX_01(3)    
      INTEGER,INTENT(IN)   :: NMAX_1 (3)    
      INTEGER,INTENT(IN)   :: NMAX_2 (3)    
      INTEGER,INTENT(IN)   :: NMAX_3 (3)    
      INTEGER,INTENT(IN)   :: NMAX_4 (3)    
!
      REAL(8),INTENT(IN)   :: SIB_MTX_SNOW (IDX,ISN,ISN+1)
      REAL(8),INTENT(IN)   :: HDF_SNOW_H   (IDX,ISN)
      REAL(8),INTENT(IN)   :: PHASE_SNOW   (IDX,ISN)
      REAL(8),INTENT(IN)   :: DELT_SNOW  
!
      REAL(8),INTENT(IN)   :: CVR_SNW      (IDX)
!
      REAL(8),INTENT(IN)   :: SRAD_BUCKET_SNOW_TO_SNSL(IDX) 
!
      REAL(8),INTENT(INOUT):: HDIFF_FLUX_H(IDX,0:ISN) 
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: ENG_BUCKET_SNOW_TO_SNSL(IDX) 
!
! INOUT 
! 
      REAL(8),INTENT(INOUT):: HFLUX_CONV(IDX,ISN)  ! 雪の各層が受け取る熱
!
! WORK
!
      INTEGER :: I
      REAL(8) :: DELT_TMP
      REAL(8) :: HDIFF_FLUX_H_CORR 
      REAL(8) :: ENG_BUC_SNOW_TO_SNSL(IDX) 
!
!
#ifdef MONYOS
      REAL(8),PARAMETER  :: ONE = 1.D0 
!     REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
#endif
!
      LOGICAL,SAVE :: LFIRST 
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_HFLUX_IMPLICIT_CORR POCHAKA 99/08/09' 
        LFIRST = .FALSE.
      ENDIF
!
! エクスプリシット評価に、インプリシット分を加える。
!
      DELT_TMP = DELT_SNOW * DELT_FAC_SNOW_c
!
! 第一層
!     一層しかない場合は PHASE_SNOW(I,2) = 0
!    
      DO 1100 I=NMAX_01(2),NMAX
        HDIFF_FLUX_H_CORR = 
     1         HDF_SNOW_H(I,1) * DELT_TMP 
     1           * ( SIB_MTX_SNOW(I,1,5) * PHASE_SNOW(I,1) 
     1                - SIB_MTX_SNOW(I,2,5) * PHASE_SNOW(I,2) )
        HDIFF_FLUX_H(I,1) = HDIFF_FLUX_H(I,1) + HDIFF_FLUX_H_CORR
 1100 CONTINUE  
!
! 第二層
!     二層しかない場合は PHASE_SNOW(I,3) = 0
!    
      DO 1200 I=NMAX_2(2),NMAX
        HDIFF_FLUX_H_CORR = 
     1         HDF_SNOW_H(I,2) * DELT_TMP 
     1           * ( SIB_MTX_SNOW(I,2,5) * PHASE_SNOW(I,2) 
     1                - SIB_MTX_SNOW(I,3,5) * PHASE_SNOW(I,3) )
        HDIFF_FLUX_H(I,2) = HDIFF_FLUX_H(I,2) + HDIFF_FLUX_H_CORR
 1200 CONTINUE  
!
! 第三層
!     三層しかない場合は PHASE_SNOW(I,4) = 0
!    
      DO 1300 I=NMAX_3(2),NMAX
        HDIFF_FLUX_H_CORR = 
     1         HDF_SNOW_H(I,3) * DELT_TMP 
     1           * ( SIB_MTX_SNOW(I,3,5) * PHASE_SNOW(I,3) 
     1                - SIB_MTX_SNOW(I,4,5) * PHASE_SNOW(I,4) )
        HDIFF_FLUX_H(I,3) = HDIFF_FLUX_H(I,3) + HDIFF_FLUX_H_CORR
 1300 CONTINUE  
!
! 第四層
!    
      DO 1400 I=NMAX_4(2),NMAX
        HDIFF_FLUX_H_CORR = 
     1         HDF_SNOW_H(I,4) * DELT_TMP 
     1           *  SIB_MTX_SNOW(I,4,5) * PHASE_SNOW(I,4) 
        HDIFF_FLUX_H(I,4) = HDIFF_FLUX_H(I,4) + HDIFF_FLUX_H_CORR
 1400 CONTINUE  
!
      DO I=NMAX_01(2),NMAX_1(3)
        HDIFF_FLUX_H(I,2) = HDIFF_FLUX_H(I,1) 
        HDIFF_FLUX_H(I,3) = HDIFF_FLUX_H(I,1) 
        HDIFF_FLUX_H(I,4) = HDIFF_FLUX_H(I,1) 
      ENDDO 
      DO I=NMAX_2(2),NMAX_2(3)
        HDIFF_FLUX_H(I,3) = HDIFF_FLUX_H(I,2) 
        HDIFF_FLUX_H(I,4) = HDIFF_FLUX_H(I,2) 
      ENDDO 
      DO I=NMAX_3(2),NMAX_3(3)
        HDIFF_FLUX_H(I,4) = HDIFF_FLUX_H(I,3) 
      ENDDO 
!
      DO I=NMAX_01(2),NMAX 
        HFLUX_CONV(I,1) = HFLUX_CONV(I,1) 
     1                    + HDIFF_FLUX_H(I,0) - HDIFF_FLUX_H(I,1)  
        HFLUX_CONV(I,2) = HFLUX_CONV(I,2) 
     1                    + HDIFF_FLUX_H(I,1) - HDIFF_FLUX_H(I,2)  
        HFLUX_CONV(I,3) = HFLUX_CONV(I,3) 
     1                    + HDIFF_FLUX_H(I,2) - HDIFF_FLUX_H(I,3)  
        HFLUX_CONV(I,4) = HFLUX_CONV(I,4) 
     1                    + HDIFF_FLUX_H(I,3) - HDIFF_FLUX_H(I,4)  
      ENDDO
!
!  ========================
!  > 地面への熱フラックス <  （下向き正）
!  ========================
!
!  --------------------
!  > 熱拡散フラックス <
!  --------------------
!
      DO 2100 I=NMAX_01(2),NMAX_1(3)
        ENG_BUC_SNOW_TO_SNSL(I) 
     1      =  HDIFF_FLUX_H(I,1)   
 2100 CONTINUE  
!
      DO 2200 I=NMAX_2(2),NMAX_2(3)
        ENG_BUC_SNOW_TO_SNSL(I)  
     1      =    HDIFF_FLUX_H(I,2) 
 2200 CONTINUE  
!
      DO 2300 I=NMAX_3(2),NMAX_3(3)
        ENG_BUC_SNOW_TO_SNSL(I)  
     1      =    HDIFF_FLUX_H(I,3) 
 2300 CONTINUE  
!
      DO 2400 I=NMAX_4(2),NMAX_4(3)
        ENG_BUC_SNOW_TO_SNSL(I)  
     1      =    HDIFF_FLUX_H(I,4) 
 2400 CONTINUE  
!
!  --------------------------------
!  > 雪あり・雪なし共存域での補正 <
!  --------------------------------
!     ここまでの計算では ENG_BUC_SNOW_TO_SNSL は雪あり域だけ
!     での平均量であった。
!     ここで格子平均にする。
!     これらがそのまま土壌への入力になる。
!
      ENG_BUCKET_SNOW_TO_SNSL(:) = 0.
!
      DO 3000 I=NMAX_01(2),NMAX
        ENG_BUCKET_SNOW_TO_SNSL(I) = 
     1         + ENG_BUC_SNOW_TO_SNSL(I) 
     1               * CVR_SNW(I) * DELT_SNOW
     1         + SRAD_BUCKET_SNOW_TO_SNSL(I)    
 3000 CONTINUE  
!
!  --------------
!  > 詳細モニタ <
!  --------------
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
!  雪あり土壌への熱
!
        CMARK  = 'HFLBWL' 
        CALL MONIT_SIB ( CMARK , ENG_BUC_SNOW_TO_SNSL,ICS,DELT_SNOW,ONE)  
      ENDIF
!
! 熱伝導フラックス  （インプリシット補正項を加算）
!
      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
        CMARK  = 'HFLW12F' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,1), ICS, 
     1                   DELT_SNOW, ONE)  
!
        CMARK  = 'HFLW23F' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,2), ICS, 
     1                   DELT_SNOW, ONE)  
!
        CMARK  = 'HFLW34F' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,3), ICS, 
     1                   DELT_SNOW, ONE)  
!
        CMARK  = 'HFLW4LF' 
        CALL MONIT_SIB ( CMARK , HDIFF_FLUX_H(1,4), ICS, 
     1                   DELT_SNOW, ONE)  
!
      ENDIF
#endif
!
      RETURN
      END

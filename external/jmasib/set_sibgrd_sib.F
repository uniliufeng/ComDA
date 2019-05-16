! SiB 予報変数を SiB 格子に詰める
! 99/07/xx 保坂征宏
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
      SUBROUTINE SET_SIBGRD_SIB (
     I   JL      , ILIST    , 
     I   NMAX    , NMAX_01  , NMAX_1  , NMAX_2   , NMAX_3   , NMAX_4   ,     
     I   TMP_CNP_NOS_PHY    , FR_WTR_CNP_NOS_PHY , FR_GLA_CNP_NOS_PHY , 
     I   TMP_CNP_SNW_PHY    , FR_WTR_CNP_SNW_PHY , FR_GLA_CNP_SNW_PHY , 
     I   TMP_GRSK_PHY       , FR_WTR_GRS_PHY     , FR_GLA_GRS_PHY     , 
     I   TMP_SNSK_PHY       , 
     I   INFO_SNOW_PHY      , 
     I   TMP_SOIL_PHY       , FR_WTR_SOIL_PHY    , FR_GLA_SOIL_PHY    , 
     I   TMP_SNSL_PHY       , FR_WTR_SNSL_PHY    , FR_GLA_SNSL_PHY    , 
     I   TMP_SNOW_PHY       , WTR_SNOW_PHY       , GLA_SNOW_PHY       , 
     I   RHO_SNOW_INV_PHY   , AGE_SNOW_PHY       , 
     I   ENG_SNOW_BUCKET_PHY, H2O_SNOW_BUCKET_PHY, 
     I   TMP_SOIL4_PHY      ,
     I   ALB_SNOW_SNW_PHY   , RAIN_1HOUR_PHY     , 
     O   TMP_CNP_NOS_SIB    , FR_WTR_CNP_NOS_SIB , FR_GLA_CNP_NOS_SIB , 
     O   TMP_CNP_SNW_SIB    , FR_WTR_CNP_SNW_SIB , FR_GLA_CNP_SNW_SIB , 
     O   TMP_GRSK_SIB       , FR_WTR_GRS_SIB     , FR_GLA_GRS_SIB     , 
     O   TMP_SNSK_SIB       , 
     O   INFO_SNOW_SIB      , 
     O   TMP_SOIL_SIB       , FR_WTR_SOIL_SIB    , FR_GLA_SOIL_SIB    , 
     O   TMP_SNSL_SIB       , FR_WTR_SNSL_SIB    , FR_GLA_SNSL_SIB    ,  
     O   TMP_SNOW_SIB       , WTR_SNOW_SIB       , GLA_SNOW_SIB       , 
     O   RHO_SNOW_INV_SIB   , AGE_SNOW_SIB       , 
     O   ENG_SNOW_BUCKET_SIB, H2O_SNOW_BUCKET_SIB ,
     O   TMP_SOIL4_SIB      , 
     O   ALB_SNOW_SNW_SIB   , RAIN_1HOUR_SIB)
!
      USE PRM , ONLY :
     I   IJPHY , 
     1   JLPHY , 
     1   ISPT 
!
      USE SIBPRM , ONLY : 
     1   IDX  , 
     1   ICG  , 
     1   IDP  , 
     1   ISN  , IVN
!
      USE SIBCON, ONLY :
     1   AGE_SNOW_MAX2_c   , 
     1   AGE_SNOW_MAX2_M_c 
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: JL
      INTEGER,INTENT(IN)  :: NMAX
      INTEGER,INTENT(IN)  :: NMAX_01(3)
      INTEGER,INTENT(IN)  :: NMAX_1 (3)
      INTEGER,INTENT(IN)  :: NMAX_2 (3)
      INTEGER,INTENT(IN)  :: NMAX_3 (3)
      INTEGER,INTENT(IN)  :: NMAX_4 (3)
      INTEGER,INTENT(IN)  :: ILIST              (IDX)  ! SiB 格子 -> 全 SiB 
!
      REAL(8),INTENT(IN)  :: TMP_CNP_NOS_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(IN)  :: FR_WTR_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(IN)  :: FR_GLA_CNP_NOS_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(IN)  :: TMP_CNP_SNW_PHY    (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(IN)  :: FR_WTR_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(IN)  :: FR_GLA_CNP_SNW_PHY (ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(IN)  :: TMP_GRSK_PHY       (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(IN)  :: FR_WTR_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(IN)  :: FR_GLA_GRS_PHY     (ISPT*IJPHY,JLPHY)   
      REAL(8),INTENT(IN)  :: TMP_SNSK_PHY       (ISPT*IJPHY,JLPHY)
!
      REAL(8),INTENT(IN)  :: TMP_SOIL_PHY       (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(IN)  :: TMP_SOIL4_PHY      (ISPT*IJPHY,JLPHY)
      REAL(8),INTENT(IN)  :: FR_WTR_SOIL_PHY    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(IN)  :: FR_GLA_SOIL_PHY    (ISPT*IJPHY,JLPHY,IDP)
      REAL(8),INTENT(IN)  :: TMP_SNSL_PHY       (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),INTENT(IN)  :: FR_WTR_SNSL_PHY    (ISPT*IJPHY,JLPHY,IDP)    
      REAL(8),INTENT(IN)  :: FR_GLA_SNSL_PHY    (ISPT*IJPHY,JLPHY,IDP)    
!
      REAL(8),INTENT(IN)  :: TMP_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(IN)  :: WTR_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(IN)  :: GLA_SNOW_PHY       (ISPT*IJPHY,JLPHY,ISN)
      REAL(8),INTENT(IN)  :: RHO_SNOW_INV_PHY   (ISPT*IJPHY,JLPHY,ISN)    
!
      REAL(8),INTENT(IN)  :: AGE_SNOW_PHY       (ISPT*IJPHY,JLPHY)
!
      REAL(8),INTENT(IN)  :: ENG_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY)    
      REAL(8),INTENT(IN)  :: H2O_SNOW_BUCKET_PHY(ISPT*IJPHY,JLPHY)    
!
      REAL(8),INTENT(IN)  :: ALB_SNOW_SNW_PHY   (ISPT*IJPHY,JLPHY,IVN)    
      REAL(8),INTENT(IN)  :: RAIN_1HOUR_PHY     (ISPT*IJPHY,JLPHY)    
!
      INTEGER,INTENT(IN)  :: INFO_SNOW_PHY      (ISPT*IJPHY,JLPHY) 
!
! OUTPUT
!
      REAL(8),INTENT(OUT) :: TMP_CNP_NOS_SIB    (IDX)
      REAL(8),INTENT(OUT) :: FR_WTR_CNP_NOS_SIB (IDX)    
      REAL(8),INTENT(OUT) :: FR_GLA_CNP_NOS_SIB (IDX)    
      REAL(8),INTENT(OUT) :: TMP_CNP_SNW_SIB    (IDX)
      REAL(8),INTENT(OUT) :: FR_WTR_CNP_SNW_SIB (IDX)    
      REAL(8),INTENT(OUT) :: FR_GLA_CNP_SNW_SIB (IDX)    
      REAL(8),INTENT(OUT) :: TMP_GRSK_SIB       (IDX)
      REAL(8),INTENT(OUT) :: FR_WTR_GRS_SIB     (IDX)   
      REAL(8),INTENT(OUT) :: FR_GLA_GRS_SIB     (IDX)   
      REAL(8),INTENT(OUT) :: TMP_SNSK_SIB       (IDX)
!
      REAL(8),INTENT(OUT) :: TMP_SOIL_SIB       (IDX,IDP)
      REAL(8),INTENT(OUT) :: FR_WTR_SOIL_SIB    (IDX,IDP)
      REAL(8),INTENT(OUT) :: FR_GLA_SOIL_SIB    (IDX,IDP)
      REAL(8),INTENT(OUT) :: TMP_SNSL_SIB       (IDX,IDP)    
      REAL(8),INTENT(OUT) :: FR_WTR_SNSL_SIB    (IDX,IDP)    
      REAL(8),INTENT(OUT) :: FR_GLA_SNSL_SIB    (IDX,IDP)    
!
      REAL(8),INTENT(OUT) :: TMP_SNOW_SIB       (IDX,ISN+1) ! 雪の温度
!                                                           ! 最後は TMP_SNSL 
      REAL(8),INTENT(OUT) :: WTR_SNOW_SIB       (IDX,ISN)
      REAL(8),INTENT(OUT) :: GLA_SNOW_SIB       (IDX,ISN)
      REAL(8),INTENT(OUT) :: RHO_SNOW_INV_SIB   (IDX,ISN)    
!
      REAL(8),INTENT(OUT) :: AGE_SNOW_SIB       (IDX)
!
      REAL(8),INTENT(OUT) :: ENG_SNOW_BUCKET_SIB(IDX)
      REAL(8),INTENT(OUT) :: H2O_SNOW_BUCKET_SIB(IDX)
!
      REAL(8),INTENT(OUT) :: TMP_SOIL4_SIB      (IDX)
!
      REAL(8),INTENT(OUT)  :: ALB_SNOW_SNW_SIB  (IDX,IVN)
      REAL(8),INTENT(OUT)  :: RAIN_1HOUR_SIB   (IDX)
!
      INTEGER,INTENT(OUT) :: INFO_SNOW_SIB      (IDX) 
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: J
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SET_SIBGRD_SIB 99/07 POCHAKA'
        LFIRST = .FALSE.
      ENDIF

!
      IF ( NMAX.LT.1 ) RETURN
!
! ------------
! > チェック <
! ------------
!
! 雪の年齢
!
#ifdef SIB_DEBUG
      DO I=1,NMAX_01(3)
        J = ILIST(I)
        IF ( AGE_SNOW_PHY(J,JL) .LT. AGE_SNOW_MAX2_M_c ) THEN
          WRITE(6,*) 'SET_SIBGRD_SIB ERROR AGE_SNOW_PHY ' , 
     1      J , JL , AGE_SNOW_PHY(J,JL) 
        ENDIF
      ENDDO  
#endif
!
        TMP_CNP_NOS_SIB   (:)    = 0.
        FR_WTR_CNP_NOS_SIB(:)    = 0.
        FR_GLA_CNP_NOS_SIB(:)    = 0.
        TMP_GRSK_SIB      (:)    = 0. 
        FR_WTR_GRS_SIB    (:)    = 0.
        FR_GLA_GRS_SIB    (:)    = 0. 
!
        TMP_CNP_SNW_SIB   (:)    = 0.
        FR_WTR_CNP_SNW_SIB(:)    = 0.
        FR_GLA_CNP_SNW_SIB(:)    = 0.
        TMP_SNSK_SIB      (:)    = 0. 
!
        AGE_SNOW_SIB      (:)    = AGE_SNOW_MAX2_c  
        TMP_SNOW_SIB      (:,:)  = 0.
        WTR_SNOW_SIB      (:,:)  = 0.
        GLA_SNOW_SIB      (:,:)  = 0.
        RHO_SNOW_INV_SIB  (:,:)  = 0.
!
        TMP_SOIL_SIB      (:,:)  = 0.
        FR_WTR_SOIL_SIB   (:,:)  = 0.
        FR_GLA_SOIL_SIB   (:,:)  = 0.
!
        TMP_SNSL_SIB      (:,:)  = 0.
        FR_WTR_SNSL_SIB   (:,:)  = 0.
        FR_GLA_SNSL_SIB   (:,:)  = 0.
!
! --------
! > 全域 <
! --------
!
      DO 1000 I=1,NMAX        
        J = ILIST(I)
        INFO_SNOW_SIB      (I) = INFO_SNOW_PHY      (J,JL) 
        H2O_SNOW_BUCKET_SIB(I) = H2O_SNOW_BUCKET_PHY(J,JL)
        ENG_SNOW_BUCKET_SIB(I) = ENG_SNOW_BUCKET_PHY(J,JL)
        TMP_SOIL4_SIB      (I) = TMP_SOIL4_PHY      (J,JL)
        ALB_SNOW_SNW_SIB   (I,:) = ALB_SNOW_SNW_PHY   (J,JL,:) 
        RAIN_1HOUR_SIB     (I) = RAIN_1HOUR_PHY     (J,JL)  
 1000 CONTINUE  
!
! ------------
! > 雪なし域 <
! ------------
!
      DO 2000 I=1,NMAX_01(3)
        J = ILIST(I)
        TMP_CNP_NOS_SIB   (I)    = TMP_CNP_NOS_PHY   (J,JL) 
        FR_WTR_CNP_NOS_SIB(I)    = FR_WTR_CNP_NOS_PHY(J,JL)  
        FR_GLA_CNP_NOS_SIB(I)    = FR_GLA_CNP_NOS_PHY(J,JL)  
!
        TMP_GRSK_SIB      (I)    = TMP_GRSK_PHY      (J,JL) 
        FR_WTR_GRS_SIB    (I)    = FR_WTR_GRS_PHY    (J,JL)  
        FR_GLA_GRS_SIB    (I)    = FR_GLA_GRS_PHY    (J,JL)  
!
        TMP_SOIL_SIB      (I,1)  = TMP_SOIL_PHY      (J,JL,1)
        TMP_SOIL_SIB      (I,2)  = TMP_SOIL_PHY      (J,JL,2)
        TMP_SOIL_SIB      (I,3)  = TMP_SOIL_PHY      (J,JL,3)
        FR_WTR_SOIL_SIB   (I,1)  = FR_WTR_SOIL_PHY   (J,JL,1)
        FR_WTR_SOIL_SIB   (I,2)  = FR_WTR_SOIL_PHY   (J,JL,2)
        FR_WTR_SOIL_SIB   (I,3)  = FR_WTR_SOIL_PHY   (J,JL,3)
        FR_GLA_SOIL_SIB   (I,1)  = FR_GLA_SOIL_PHY   (J,JL,1)
        FR_GLA_SOIL_SIB   (I,2)  = FR_GLA_SOIL_PHY   (J,JL,2)
        FR_GLA_SOIL_SIB   (I,3)  = FR_GLA_SOIL_PHY   (J,JL,3)
!
 2000 CONTINUE
!
! ------------
! > 雪あり域 <
! ------------
!
      DO 3000 I=NMAX_01(2),NMAX   
        J = ILIST(I)
        TMP_CNP_SNW_SIB   (I)    = TMP_CNP_SNW_PHY   (J,JL) 
        FR_WTR_CNP_SNW_SIB(I)    = FR_WTR_CNP_SNW_PHY(J,JL)  
        FR_GLA_CNP_SNW_SIB(I)    = FR_GLA_CNP_SNW_PHY(J,JL)  
        TMP_SNSK_SIB      (I)    = TMP_SNSK_PHY      (J,JL) 
!
        AGE_SNOW_SIB      (I)    = AGE_SNOW_PHY      (J,JL)
        TMP_SNOW_SIB      (I,1)  = TMP_SNOW_PHY      (J,JL,1)
        WTR_SNOW_SIB      (I,1)  = WTR_SNOW_PHY      (J,JL,1)
        GLA_SNOW_SIB      (I,1)  = GLA_SNOW_PHY      (J,JL,1)
        RHO_SNOW_INV_SIB  (I,1)  = RHO_SNOW_INV_PHY  (J,JL,1)
!
        TMP_SNSL_SIB      (I,1)  = TMP_SNSL_PHY      (J,JL,1)
        TMP_SNSL_SIB      (I,2)  = TMP_SNSL_PHY      (J,JL,2)
        TMP_SNSL_SIB      (I,3)  = TMP_SNSL_PHY      (J,JL,3)
        FR_WTR_SNSL_SIB   (I,1)  = FR_WTR_SNSL_PHY   (J,JL,1)
        FR_WTR_SNSL_SIB   (I,2)  = FR_WTR_SNSL_PHY   (J,JL,2)
        FR_WTR_SNSL_SIB   (I,3)  = FR_WTR_SNSL_PHY   (J,JL,3)
        FR_GLA_SNSL_SIB   (I,1)  = FR_GLA_SNSL_PHY   (J,JL,1)
        FR_GLA_SNSL_SIB   (I,2)  = FR_GLA_SNSL_PHY   (J,JL,2)
        FR_GLA_SNSL_SIB   (I,3)  = FR_GLA_SNSL_PHY   (J,JL,3)

 3000 CONTINUE
!
!    メモ  ここから下は、土壌・雪計算ステップのみ必要なのかも知れない
!
! 雪第二層 
!
      DO 3100 I=NMAX_2(2),NMAX  
        J = ILIST(I)
        TMP_SNOW_SIB      (I,2)  = TMP_SNOW_PHY      (J,JL,2)
        WTR_SNOW_SIB      (I,2)  = WTR_SNOW_PHY      (J,JL,2)
        GLA_SNOW_SIB      (I,2)  = GLA_SNOW_PHY      (J,JL,2)
        RHO_SNOW_INV_SIB  (I,2)  = RHO_SNOW_INV_PHY  (J,JL,2)
 3100 CONTINUE
!
! 雪第三層 
!
      DO 3200 I=NMAX_3(2),NMAX  
        J = ILIST(I)
        TMP_SNOW_SIB      (I,3)  = TMP_SNOW_PHY      (J,JL,3)
        WTR_SNOW_SIB      (I,3)  = WTR_SNOW_PHY      (J,JL,3)
        GLA_SNOW_SIB      (I,3)  = GLA_SNOW_PHY      (J,JL,3)
        RHO_SNOW_INV_SIB  (I,3)  = RHO_SNOW_INV_PHY  (J,JL,3)
 3200 CONTINUE  
!
! 雪第四層
!
      DO 3300 I=NMAX_4(2),NMAX 
        J = ILIST(I)
        TMP_SNOW_SIB      (I,4)  = TMP_SNOW_PHY      (J,JL,4)
        WTR_SNOW_SIB      (I,4)  = WTR_SNOW_PHY      (J,JL,4)
        GLA_SNOW_SIB      (I,4)  = GLA_SNOW_PHY      (J,JL,4)
        RHO_SNOW_INV_SIB  (I,4)  = RHO_SNOW_INV_PHY  (J,JL,4)
 3300 CONTINUE
!
!  雪の温度には、最下層の下に SNSL の温度の層を与える
!      （雪の）
!
      DO 4000 I=NMAX_01(2),NMAX_1(3)
        TMP_SNOW_SIB(I,2) = TMP_SNSL_SIB(I,1)
 4000 CONTINUE  
      DO 4100 I=NMAX_2(2) ,NMAX_2(3)
        TMP_SNOW_SIB(I,3) = TMP_SNSL_SIB(I,1)
 4100 CONTINUE  
      DO 4200 I=NMAX_3(2) ,NMAX_3(3)
        TMP_SNOW_SIB(I,4) = TMP_SNSL_SIB(I,1)
 4200 CONTINUE  
      DO 4300 I=NMAX_4(2) ,NMAX_4(3)
        TMP_SNOW_SIB(I,5) = TMP_SNSL_SIB(I,1)
 4300 CONTINUE  
!
!
      RETURN
      END

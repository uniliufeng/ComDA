      SUBROUTINE MONIT_SNAP_SNOW (
     I    NMAX    , 
     I    IMEAN           , INFO_MEAN , 
     I    DELT_SNOW       , DELT_CNP , 
     I    INFO_SNOW_NEXT   , 
     I    CVR_SNW_NEXT     , 
     I    WTR_SNOW_NEXT       , GLA_SNOW_NEXT       , 
     I    TMP_SNOW_NEXT       , DEP_SNOW_NEXT       , 
     I    H2O_SNOW_BUCKET_NEXT , 
     I    ENG_SNOW_BUCKET_NEXT , 
     I    RHO_SNOW_INV_NEXT    
     1   )  
!
      USE SIBPRM , ONLY : 
     1    IDX , 
     1    ISN , 
     1    IDP , 
     1    ICG 
!
      USE SIBCON , ONLY :
     1    CAP_LAI_c  ,
     1    CAP_WTR_c  ,
     1    CAP_GLA_c  , 
     1    TMP_FREZ_c , 
     1    HL_W2G_c   , 
     1    RHO_WTR_c  , 
     1    RHO_GLA_c  , 
     1    H2O_SAT_CNP_DEP_c ,
     1    H2O_SAT_GRS_DEP_c
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB_NEXT, INS,INA,INN ,ICS,
     1                       IMONIT_LEVEL , IMONIT_ALL , 
     1                       IMONIT_MEAN_MRI,
     1                       IMONIT_SNAP_NPD,
     1                       IMONIT_ALMA 
!
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif
    
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)    :: IMEAN        ! 平均値も出す
      INTEGER,INTENT(IN)    :: INFO_MEAN    ! 負なら delt_mean=DELT_CNP  
      INTEGER,INTENT(IN)    :: NMAX       
!
!     INTEGER,INTENT(IN)    :: NMAX_01(3)
      INTEGER,INTENT(IN)    :: INFO_SNOW_NEXT  (IDX)
      REAL(8),INTENT(IN)    :: CVR_SNW_NEXT    (IDX) 
      REAL(8),INTENT(IN)    :: DELT_SNOW
      REAL(8),INTENT(IN)    :: DELT_CNP
!
      REAL(8),INTENT(IN)    :: DEP_SNOW_NEXT        (IDX,ISN+1) 
      REAL(8),INTENT(IN)    :: WTR_SNOW_NEXT        (IDX,ISN) 
      REAL(8),INTENT(IN)    :: GLA_SNOW_NEXT        (IDX,ISN) 
      REAL(8),INTENT(IN)    :: RHO_SNOW_INV_NEXT    (IDX,ISN) 
      REAL(8),INTENT(IN)    :: TMP_SNOW_NEXT        (IDX,ISN) 
      REAL(8),INTENT(IN)    :: H2O_SNOW_BUCKET_NEXT      (IDX) 
      REAL(8),INTENT(IN)    :: ENG_SNOW_BUCKET_NEXT      (IDX) 
!
! WORK
!     
#ifdef MONYOS
      REAL(8),PARAMETER :: ONE = 1. 
      CHARACTER(7) :: CMARK 
      REAL(8)      :: WORK   (IDX)
      REAL(8)      :: WORK2  (IDX)
      REAL(8)      :: DELT_MEAN 
      real(8)      :: delt_snap
      REAL(8)      :: BUNSI
      REAL(8)      :: BUNBO
      REAL(8),parameter :: SMALL = 1.D-30
#endif
!
      INTEGER   ::  I
      INTEGER   ::  K
      INTEGER   ::  KK
!
      LOGICAL,SAVE :: LFIRST
!
      DATA  LFIRST /.TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_SNAP_SNOW 99/10 POCHAKA'

        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( INFO_MEAN .LT. 0 ) THEN
        DELT_MEAN = DELT_CNP 
      ELSE
        DELT_MEAN = DELT_SNOW
      ENDIF
      DELT_SNAP   = DELT_CNP
!
!
! ========================
! >> 雪、収支チェック用 <<  
! ========================
! 注意... 雪の質量等は、もともと CVR_SNOW がかかっている!!! 
!
! メモ... 3 ( 雪がある場合だけ足す ) ではなく -1 ( 常に足す ) でもいい
!         気がする。
!
! 注意。
!    ・瞬間値については DELT_CNP をかける。
!       -- 出力直前に DELT_CNP で割るため。      
!    ・時間平均値は DELT_MEAN をかける。
!       -- 出力直前に積算時間で割るため。
!
! 積雪水
!
      IF ( (      IMONIT_LEVEL .GE. IMONIT_SNAP_NPD 
     1       .OR. IMONIT_LEVEL .GE. IMONIT_ALL        )
     1      .AND. ICN_SIB0109_MONITOR.EQ. 1 
     1   ) THEN
        DO K=1,4
         CMARK = 'WTRWxX' 
         WRITE(CMARK(5:5),'(I1)') K 
         CALL MONIT_SIB_NEXT (CMARK,WTR_SNOW_NEXT(1,K),INS,DELT_CNP,ONE)   
        ENDDO
      ENDIF
!
! 積雪氷
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .AND. ICN_SIB0109_MONITOR.EQ.1
     1   ) THEN
        DO K=1,4
         CMARK = 'GLAWxX' 
         WRITE(CMARK(5:5),'(I1)') K 
         CALL MONIT_SIB_NEXT (CMARK,GLA_SNOW_NEXT(1,K),INS,DELT_CNP,ONE)      
        ENDDO
!
! 積雪エネルギー
!
        DO K=1,4
          CMARK = 'ENGWxX' 
          WRITE(CMARK(5:5),'(I1)') K 
          DO I=1,NMAX
            WORK(I) = 
     1              (   CAP_GLA_c * GLA_SNOW_NEXT(I,K)  
     1                + CAP_WTR_c * WTR_SNOW_NEXT(I,K) ) 
     1              *  ( TMP_SNOW_NEXT(I,K) - TMP_FREZ_c )
     1              - HL_W2G_c * GLA_SNOW_NEXT(I,K)  
          ENDDO
          CALL MONIT_SIB_NEXT (CMARK,WORK,INS,DELT_CNP,ONE)      
        ENDDO
!
! 少雪バケツ
!
        CMARK = 'ENGLSNW' 
        CALL MONIT_SIB_NEXT (CMARK,ENG_SNOW_BUCKET_NEXT,
     1                             INN,DELT_CNP,ONE)      
        CMARK = 'H2OLSNW' 
        CALL MONIT_SIB_NEXT (CMARK,H2O_SNOW_BUCKET_NEXT,
     1                               INN,DELT_CNP,ONE)      
      ENDIF
!
! ==========================
! >> 雪、クイックルック用 <<
! ==========================
! 注意... 被覆率を乗じていないので、まさに瞬間値。
!         雪がなくても足している。
!         描画で /cvrsnwm などはしないこと。
!
      IF (  (     IMONIT_LEVEL .GE. IMONIT_ALL 
     1       .OR. IMONIT_LEVEL .GE. IMONIT_SNAP_NPD )
     1     .AND. ICN_SIB0109_MONITOR.EQ.1
     1    ) THEN
!
! 積雪温度
!
        DO K=1,ISN
          CMARK = 'TMPWxX' 
          WRITE(CMARK(5:5),'(I1)') K 
          CALL MONIT_SIB_NEXT(CMARK,TMP_SNOW_NEXT(1,K),
     1                        INA,DELT_CNP,ONE)      
        ENDDO
!
! 積雪密度
!
        DO K=1,4
          KK = K
          CMARK = 'RHOWxX' 
          IF ( K .EQ. 1 ) KK = 0
          CALL RESET ( WORK , IDX ) 
          DO I=1,NMAX
            IF ( INFO_SNOW_NEXT(I) .GE. KK ) THEN
              WORK(I) = 1./RHO_SNOW_INV_NEXT(I,K)
            ENDIF
          ENDDO
          WRITE(CMARK(5:5),'(I1)') K 
          CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_CNP,ONE)      
        ENDDO
!
! 積雪被覆率
!
        CMARK = 'CVRSNWX' 
        CALL MONIT_SIB_NEXT(CMARK,CVR_SNW_NEXT,INA,DELT_CNP,ONE)      
!
      ENDIF
!
! 積雪深
!
      IF (       IMONIT_LEVEL .GE. IMONIT_ALL 
     1     .AND. ICN_SIB0109_MONITOR.EQ.1
     1   ) THEN
!
        DO K=1,ISN
          CMARK = 'DEPWxX' 
          WRITE(CMARK(5:5),'(I1)') K 
          CALL MONIT_SIB_NEXT(CMARK,DEP_SNOW_NEXT(1,K),INA,DELT_CNP,ONE)      
        ENDDO
!
! 積雪状況
!
        WORK(:) = -9.99E33
        DO I=1,NMAX
          IF     ( INFO_SNOW_NEXT(I).EQ.-2) THEN  ! sukunai
            WORK(I) = 0.25
          ELSEIF ( INFO_SNOW_NEXT(I).EQ.-1) THEN
            WORK(I) = 0
          ELSEIF ( INFO_SNOW_NEXT(I).EQ. 0) THEN
            WORK(I) = 0.5
          ELSEIF ( INFO_SNOW_NEXT(I).GE. 1) THEN
            WORK(I) = INFO_SNOW_NEXT(I)
          ENDIF 
        ENDDO
        CMARK  = 'INFSNWX'
        CALL MONIT_SIB_NEXT(CMARK,WORK,INA,DELT_CNP,ONE)      
!
      ENDIF
!
! ==================
! >> 雪、その他用 <<
! ==================
!   雪がないときも合わせた平均値。/cvrsnwm 等はしないこと。
!   ( お絵描き時に必要 ) 
!
      IF ( 
     1          IMONIT_LEVEL .GE. IMONIT_MEAN_MRI
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL 
     1     ) THEN
!
        IF ( IMEAN .EQ. 1 ) THEN
!
! 積雪被覆率
!
          CMARK = 'CVRSNWM' 
          CALL MONIT_SIB_NEXT (CMARK,CVR_SNW_NEXT,INA,DELT_MEAN,ONE)      
!
! 積雪深  
!
          CMARK  = 'DEPSNWM'
          CALL MONIT_SIB_NEXT (CMARK,DEP_SNOW_NEXT(1,5),
!    1                         INA,DELT_MEAN,ONE)      
     1                         ICS,DELT_MEAN,ONE)      
!
! 積雪温度
!
        DO K=1,ISN
          CMARK = 'TMPWxM' 
          WRITE(CMARK(5:5),'(I1)') K 
          CALL MONIT_SIB_NEXT(CMARK,TMP_SNOW_NEXT(1,K),
!    1                        INA,DELT_CNP,ONE)      
     1                        ICS,DELT_CNP,ONE)      
        ENDDO
!
! 積雪密度
!
        DO K=1,4
          KK = K
          CMARK = 'RHOWxM' 
          IF ( K .EQ. 1 ) KK = 0
          CALL RESET ( WORK , IDX ) 
          DO I=1,NMAX
            IF ( INFO_SNOW_NEXT(I) .GE. KK ) THEN
              WORK(I) = 1./RHO_SNOW_INV_NEXT(I,K)
            ENDIF
          ENDDO
          WRITE(CMARK(5:5),'(I1)') K 
!         CALL MONIT_SIB_NEXT (CMARK,WORK,INA,DELT_CNP,ONE)      
          CALL MONIT_SIB_NEXT (CMARK,WORK,ICS,DELT_CNP,ONE)      
        ENDDO
!
! 積雪質量
!
        DO K=1,ISN
!
          WORK(:) = 0 
          DO I=1,NMAX
            WORK(I) = GLA_SNOW_NEXT(I,K) + WTR_SNOW_NEXT(I,K)
          ENDDO

          CMARK = 'WGTWxM' 
          WRITE(CMARK(5:5),'(I1)') K 
          CALL MONIT_SIB_NEXT(CMARK,WORK, 
!    1                        INA,DELT_CNP,ONE)      
     1                        ICS,DELT_CNP,ONE)      
        ENDDO
!
! 積雪平均温度 K 
!
          WORK (:) = 0.
          WORK2(:) = 0.
          DO I=1,NMAX
            BUNBO =   GLA_SNOW_NEXT(I,1) + WTR_SNOW_NEXT(I,1)
     1              + GLA_SNOW_NEXT(I,2) + WTR_SNOW_NEXT(I,2)
     1              + GLA_SNOW_NEXT(I,3) + WTR_SNOW_NEXT(I,3)
     1              + GLA_SNOW_NEXT(I,4) + WTR_SNOW_NEXT(I,4)
            BUNSI =   TMP_SNOW_NEXT(I,1) 
     1                  * ( GLA_SNOW_NEXT(I,1) + WTR_SNOW_NEXT(I,1) )
     1              + TMP_SNOW_NEXT(I,2) 
     1                  * ( GLA_SNOW_NEXT(I,2) + WTR_SNOW_NEXT(I,2) )
     1              + TMP_SNOW_NEXT(I,3) 
     1                  * ( GLA_SNOW_NEXT(I,3) + WTR_SNOW_NEXT(I,3) )
     1              + TMP_SNOW_NEXT(I,4) 
     1                  * ( GLA_SNOW_NEXT(I,4) + WTR_SNOW_NEXT(I,4) )
            WORK (I) = BUNSI/MAX(BUNBO,SMALL) 
            WORK2(I) = BUNBO  
          ENDDO
!
          CMARK  = 'TMPSNWM'
          CALL MONIT_SIB_NEXT (CMARK,WORK ,
!    1                         INA,DELT_MEAN,ONE)      
     1                         ICS,DELT_MEAN,ONE)      
!
! 積雪質量 kg/m**2 
!  
          CMARK  = 'WGTSNWM'
          CALL MONIT_SIB_NEXT (CMARK,WORK2,
!    1                         INA,DELT_MEAN,ONE)      
     1                         ICS,DELT_MEAN,ONE)      
!
        ENDIF
!
      ENDIF
!
! * for ALMA 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
! A_SWEx  ( 積雪水当量 ) 
!
        DO K=1,ISN
          CMARK = 'A_SWEx' 
          WRITE(CMARK(6:6),'(I1)') K 
          WORK(:) = 0 
          DO I=1,NMAX
            WORK(I) = GLA_SNOW_NEXT(I,K) + WTR_SNOW_NEXT(I,K)
          ENDDO
!         CALL MONIT_SIB_NEXT(CMARK,WORK,INA,DELT_CNP,ONE)      
          CALL MONIT_SIB_NEXT(CMARK,WORK,ICS,DELT_CNP,ONE)      
        ENDDO
!
! A_DPSNx (積雪厚さ) 
!         
        DO K=1,ISN
          CMARK = 'A_DPSNx' 
          WRITE(CMARK(7:7),'(I1)') K 
!         CALL MONIT_SIB_NEXT(CMARK,DEP_SNOW_NEXT(1,K),INA,DELT_CNP,ONE)      
          CALL MONIT_SIB_NEXT(CMARK,DEP_SNOW_NEXT(1,K),ICS,DELT_CNP,ONE)      
        ENDDO
!
! A_FACSN ( 積雪被覆率 ) 
!
        CMARK = 'A_FACSN' 
        CALL MONIT_SIB_NEXT(CMARK,CVR_SNW_NEXT,INA,DELT_CNP,ONE)      
!
! A_SWEAX  ( 積雪水当量 鉛直積分 snap) 
!
        WORK(:) = H2O_SNOW_BUCKET_NEXT(:)  
        DO K=1,ISN
          DO I=1,NMAX
            WORK(I) = WORK(I) + GLA_SNOW_NEXT(I,K) + WTR_SNOW_NEXT(I,K)
          ENDDO
        ENDDO
        CMARK = 'A_SWEAX' 
        CALL MONIT_SIB_NEXT(CMARK,WORK,INA,DELT_SNAP,ONE)      
!
      ENDIF
!
#endif
!
      RETURN
      END
!=========================================================
      SUBROUTINE MR_MONIT_SNAP_SNOW
!
      USE SIBPRM , ONLY :
     1    ISN , 
     1    IDP 
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY :
     1                       IMONIT_LEVEL , IMONIT_ALL ,
     1                       IMONIT_MEAN_MRI,
     1                       IMONIT_SNAP_NPD,
     1                       IMONIT_ALMA 
#endif

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
        WRITE(6,*) 'MONIT_REGIST_MONIT_SNAP_SNOW 99/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF (      IMONIT_LEVEL .GE. IMONIT_MEAN_MRI
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN 
        CMARK  = 'CVRSNWM'
        CTITLE = 'SNOW COVER (TMEAN)'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'DEPSNWM'
        CTITLE = 'SNOW DEPTH (TMEAN)'
        CUNIT  = 'M'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'TMPSNWM'
        CTITLE = 'SNOW TEMP (TMEAN)'
        CUNIT  = 'K'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'WGTSNWM'
        CTITLE = 'SNOW MASS (TMEAN)'
        CUNIT  = 'kg/m**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'TMPWxM'
        CTITLE = 'SNOW TEMP LEVEL x (MEAN)'
        CUNIT  = 'K'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(17:17),'(I1)') K 
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
        ENDDO
!
        CMARK  = 'RHOWxM'
        CTITLE = 'SNOW DENS LEVEL x (MEAN)'
        CUNIT  = 'K'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(17:17),'(I1)') K 
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
        ENDDO

!
! 雪のもつ水 + 氷
!
        CMARK  = 'WGTWxM'
        CTITLE = 'SNOW WEIGHT LEVEL x (MEAN)'
        CUNIT  = 'KG/M**2'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(19:19),'(I1)') K 
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
        ENDDO
!
      ENDIF
!
      IF (      IMONIT_LEVEL .GE. IMONIT_SNAP_NPD
     1     .OR. IMONIT_LEVEL .GE. IMONIT_ALL       ) THEN 
!
! 雪の被覆率
!
        CMARK  = 'CVRSNWX'
        CTITLE = 'SNOW COVER (TSNAP)'
        CUNIT  = '0-1'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! 雪の温度
!
        CMARK  = 'TMPWxX'
        CTITLE = 'SNOW TEMP LEVEL x (SNAP)'
        CUNIT  = 'K'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(17:17),'(I1)') K 
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
        ENDDO
!
! 雪のもつ水 
!
        CMARK  = 'WTRWxX'
        CTITLE = 'SNOW WATER LEVEL x (SNAP)'
        CUNIT  = 'KG/M**2'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(18:18),'(I1)') K 
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
        ENDDO
!
! 雪密度
!
        CMARK  = 'RHOWxX'
        CTITLE = 'SNOW DENSITY LEVEL x (TSNAP)'
        CUNIT  = 'KG/M**3'
        DO K = 1,ISN 
          WRITE(CMARK (5:5)  ,'(I1)') K 
          WRITE(CTITLE(20:20),'(I1)') K 
          CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
        ENDDO
!
      ENDIF    
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
! 積雪状況
!
      CMARK  = 'INFSNWX'
      CTITLE = 'Info snow (snap)'
      CUNIT  = '1'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! ==========================
! >> 雪、収支関連スナップ <<
! ==========================
!
! 雪のもつエネルギー 
!
      CMARK  = 'ENGWxX'
      CTITLE = 'SNOW ENERGY LEVEL x (SNAP)'
      CUNIT  = 'J/M**2'
      DO K = 1,ISN 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        WRITE(CTITLE(19:19),'(I1)') K 
        CALL monit_regist_sib_SNAP ( CMARK , CTITLE , CUNIT )
      ENDDO
!
!
! 雪のもつ氷 
!
      CMARK  = 'GLAWxX'
      CTITLE = 'SNOW GLA LEVEL x (SNAP)'
      CUNIT  = 'KG/M**2'
      DO K = 1,ISN 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        WRITE(CTITLE(16:16),'(I1)') K 
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDDO
!
! 少雪バケツ
!
      CMARK  = 'ENGLSNW'
      CTITLE = 'ENH / LITTLE SNOW BUCKET (SNAP)'
      CUNIT  = 'J/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'H2OLSNW'
      CTITLE = 'H2O / LITTLE SNOW BUCKET (SNAP)'
      CUNIT  = 'KG/M**2'
      CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
!
! ================
! >> 雪、その他 <<
! ================
!
! 積雪深
!
      CMARK  = 'DEPWxX'
      CTITLE = 'SNOW DEPTH LEVEL x (TSNAP)'
      CUNIT  = 'M'
      DO K = 1,ISN 
        WRITE(CMARK (5:5)  ,'(I1)') K 
        WRITE(CTITLE(18:18),'(I1)') K 
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )
      ENDDO
!
      ENDIF
!
! * For ALMA 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
!
        CMARK  = 'A_SWEx'
        CTITLE = 'ALMA SWE'
        CUNIT  = 'kg/m**2'
        DO K = 1,ISN 
          WRITE(CMARK (6:6)  ,'(I1)') K 
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
        ENDDO
!
        CMARK  = 'A_DPSNx'
        CTITLE = 'ALMA SnowDepth'
        CUNIT  = 'M'
        DO K = 1,ISN 
          WRITE(CMARK (7:7)  ,'(I1)') K 
          CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
        ENDDO
!
        CMARK = 'A_FACSN'
        CTITLE = 'ALMA SnowFac'
        CUNIT  = '0-1'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'A_SWEAX'
        CTITLE = 'ALMA Dt DelSWE ( Tsnap )'
        CUNIT  = 'kg/m**2'
        CALL monit_regist_sib_snap ( CMARK , CTITLE , CUNIT )

      ENDIF
#endif
!
      RETURN
      END

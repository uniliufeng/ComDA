!   各種短波放射吸収量入射量を加え RAD_CG を算出する。
! 99/05/xx 保坂征宏
! 99/07/23 保坂征宏 引数変更、snw,nos をまとめて計算するようにする。
! 99/08/04 保坂征宏 大泉さんに渡すバージョンフィックス
!
!  一時間に一度だけ呼ばれる
!  注意。モニタ量は、daytime だけの値を 1 時間の重みをつけて積算している。
!        よって日の出・日の入付近は、モニタ値は過大評価されている。
!        そこで出力マークには、Z をつけた
!
      SUBROUTINE RAD_SHRT_ADD (
     I    NMAX            , NMAX_01         , 
     I    ABSORB_SHRT_NOS , ABSORB_SHRT_SNW , 
     I    RAD_INTO_SIB    , 
     I    DAYTIME_1HR       , 
     O    RAD_SHRT_NOS    , RAD_SHRT_SNW    ) 
!
! キャノピー・地面・雪層それぞれについて、
!   各種短波放射吸収量入射量を加え RAD_CG を算出する。
!
! 直達光・散乱光、可視・近赤の和をとる。
! ( 透過する分は含む？吸収分のみ？ どちら？)
!
      USE SIBPRM , ONLY : 
     1   IDX  , 
     1   ICG  , 
     1   IVN  , 
     1   IBD   
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN,ICS,INA, 
     1    IMONIT_LEVEL  , 
     1    IMONIT_FLUX_MRI   , 
     1    IMONIT_BUDJET ,
     1    IMONIT_ALL , IMONIT_ALMA
#endif
!
!      USE COM_STEP_SIB0109, only :
!     1    SCN_SIB0109_DELT_SR   


!
      IMPLICIT NONE
!
! 入力
!
      INTEGER,INTENT(IN)    :: NMAX_01(3)
      INTEGER,INTENT(IN)    :: NMAX
!
      REAL(8),INTENT(IN)    :: DAYTIME_1HR(IDX) 
      REAL(8),INTENT(IN)    :: 
     1      ABSORB_SHRT_NOS  (IDX,ICG,IVN,IBD)   ,  ! 短波の吸収率
     1      ABSORB_SHRT_SNW  (IDX,ICG,IVN,IBD)   ,  ! 短波の吸収率
     1      RAD_INTO_SIB     (IDX,    IVN,IBD)      ! 大気下端下向flux
!
! 出力
!
      REAL(8),INTENT(OUT)  ::  
     1      RAD_SHRT_NOS     (IDX,ICG)           ,  ! 短波吸収量
     1      RAD_SHRT_SNW     (IDX,ICG)              ! 短波吸収量
!
! ワーク
!
      INTEGER               :: I
      INTEGER               :: IVEG    
      REAL(8)               :: WORK(IDX)
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE  =  1.D0 
!     REAL(8),PARAMETER :: ONEM = -1.D0 
#endif
!
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'RAD_SHRT_ADD 99/07 POCHAKA '
        LFIRST = .FALSE.
      ENDIF 
!
      DO 1000 IVEG = 1, ICG
      DO 1000 I    = 1,NMAX_01(3)
        RAD_SHRT_NOS(I,IVEG) 
     1     =  ABSORB_SHRT_NOS(I,IVEG,1,1) * RAD_INTO_SIB(I,1,1)   ! 可視直達
     1      + ABSORB_SHRT_NOS(I,IVEG,2,1) * RAD_INTO_SIB(I,2,1)   ! 近赤直達
     2      + ABSORB_SHRT_NOS(I,IVEG,1,2) * RAD_INTO_SIB(I,1,2)   ! 可視散乱
     3      + ABSORB_SHRT_NOS(I,IVEG,2,2) * RAD_INTO_SIB(I,2,2)   ! 近赤散乱
 1000 CONTINUE
!
! 雪 の IVEG=2 は、雪上端に入る短波。
!
      DO 2000 IVEG = 1, ICG
      DO 2000 I    = NMAX_01(2),NMAX 
        RAD_SHRT_SNW(I,IVEG) 
     1     =  ABSORB_SHRT_SNW(I,IVEG,1,1) * RAD_INTO_SIB(I,1,1)   ! 可視直達
     1      + ABSORB_SHRT_SNW(I,IVEG,2,1) * RAD_INTO_SIB(I,2,1)   ! 近赤直達
     2      + ABSORB_SHRT_SNW(I,IVEG,1,2) * RAD_INTO_SIB(I,1,2)   ! 可視散乱
     3      + ABSORB_SHRT_SNW(I,IVEG,2,2) * RAD_INTO_SIB(I,2,2)   ! 近赤散乱
 2000 CONTINUE
!
! --------------
! > 詳細モニタ <
! --------------
!
#ifdef MONYOS
!
      IF ( 
     1             IMONIT_LEVEL .GE. IMONIT_ALL   
     1        .OR. IMONIT_LEVEL .GE. IMONIT_FLUX_MRI 
     1        .OR. IMONIT_LEVEL .GE. IMONIT_ALMA
     1   ) THEN
!
! 入射
!    下向き短波入射量 DSWBL は別に積算
!
        DO I    = 1,NMAX
          WORK(I)
     1    = (  RAD_INTO_SIB(I,1,1)                               ! 可視直達
     2       + RAD_INTO_SIB(I,2,1)                               ! 近赤直達
     3       + RAD_INTO_SIB(I,1,2)                               ! 可視散乱
     4       + RAD_INTO_SIB(I,2,2) )                             ! 近赤散乱
     5       * DAYTIME_1HR(I)   
        ENDDO
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
          CMARK = 'HSRDNTP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICN, SCN_SIB0109_DELT_SR, ONE )
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE, ONE )
          CMARK = 'HSRDSTP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICS, SCN_SIB0109_DELT_SR, ONE )
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE, ONE )
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_BUDJET ) THEN
          CMARK = 'DSWBL' 
!         CALL MONIT_SIB ( CMARK , WORK ,INA, SCN_SIB0109_DELT_SR, ONE )
          CALL MONIT_SIB ( CMARK , WORK ,INA, ONE, ONE )
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
          CMARK = 'A_SWDN' 
!         CALL MONIT_SIB ( CMARK , WORK ,INA, SCN_SIB0109_DELT_SR, ONE )
          CALL MONIT_SIB ( CMARK , WORK ,INA, ONE, ONE )
        ENDIF
!
      ENDIF
!
! 反射 ( 非吸収分 )  
!    上向き短波量 USWBL は別に積算
!
      IF ( 
     1           IMONIT_LEVEL .GE. IMONIT_ALL   
     1      .OR. IMONIT_LEVEL .GE. IMONIT_FLUX_MRI 
     1      .OR. IMONIT_LEVEL .GE. IMONIT_ALMA
     1   ) THEN
        DO I    = 1,NMAX_01(3)
          WORK(I)
     1     =   ( 1 - ABSORB_SHRT_NOS(I,1,1,1) 
     1             - ABSORB_SHRT_NOS(I,2,1,1) ) 
     1         * RAD_INTO_SIB(I,1,1)                               ! 可視直達
     2       + ( 1 - ABSORB_SHRT_NOS(I,1,2,1) 
     2             - ABSORB_SHRT_NOS(I,2,2,1) ) 
     2         * RAD_INTO_SIB(I,2,1)                               ! 近赤直達
     3       + ( 1 - ABSORB_SHRT_NOS(I,1,1,2) 
     3             - ABSORB_SHRT_NOS(I,2,1,2) ) 
     3         * RAD_INTO_SIB(I,1,2)                               ! 可視散乱
     4       + ( 1 - ABSORB_SHRT_NOS(I,1,2,2) 
     4             - ABSORB_SHRT_NOS(I,2,2,2) ) 
     4         * RAD_INTO_SIB(I,2,2)                               ! 近赤散乱
          WORK(I)=WORK(I) * DAYTIME_1HR (I) 
        ENDDO
!
        IF ( IMONIT_LEVEL .GE. IMONIT_BUDJET ) THEN
          CMARK = 'USWBL' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICN, SCN_SIB0109_DELT_SR,-ONE)
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE ,-ONE)
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
          CMARK = 'HSRUNTP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICN, SCN_SIB0109_DELT_SR, ONE)
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE, ONE)
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN 
          CMARK = 'A_SWUP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICN, SCN_SIB0109_DELT_SR, ONE)
          CALL MONIT_SIB ( CMARK , WORK ,ICN, ONE, ONE)
        ENDIF
!
        DO I    = NMAX_01(2),NMAX
          WORK(I) 
     1     =   ( 1 - ABSORB_SHRT_SNW(I,1,1,1) 
     1             - ABSORB_SHRT_SNW(I,2,1,1) ) 
     1         * RAD_INTO_SIB(I,1,1)                               ! 可視直達
     2       + ( 1 - ABSORB_SHRT_SNW(I,1,2,1) 
     2             - ABSORB_SHRT_SNW(I,2,2,1) ) 
     2         * RAD_INTO_SIB(I,2,1)                               ! 近赤直達
     3       + ( 1 - ABSORB_SHRT_SNW(I,1,1,2) 
     3             - ABSORB_SHRT_SNW(I,2,1,2) ) 
     3         * RAD_INTO_SIB(I,1,2)                               ! 可視散乱
     4       + ( 1 - ABSORB_SHRT_SNW(I,1,2,2) 
     4             - ABSORB_SHRT_SNW(I,2,2,2) ) 
     4         * RAD_INTO_SIB(I,2,2)                               ! 近赤散乱
          WORK(I) = WORK(I) * DAYTIME_1HR(I) 
        ENDDO
!
        IF ( IMONIT_LEVEL .GE. IMONIT_BUDJET ) THEN
          CMARK = 'USWBL' 
!          CALL MONIT_SIB ( CMARK , WORK ,ICS, SCN_SIB0109_DELT_SR, -ONE)
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE, -ONE)
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
          CMARK = 'HSRUSTP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICS, SCN_SIB0109_DELT_SR, ONE )
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE, ONE )
        ENDIF
!
        IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN 
          CMARK = 'A_SWUP' 
!         CALL MONIT_SIB ( CMARK , WORK ,ICS, SCN_SIB0109_DELT_SR, ONE)
          CALL MONIT_SIB ( CMARK , WORK ,ICS, ONE , ONE)
        ENDIF
!
      ENDIF
!
#endif
      RETURN
      END
!=============================================================
      SUBROUTINE MR_RAD_SHRT_ADD 
!
      USE SIB_MONIT , ONLY : 
     1    IMONIT_LEVEL  , IMONIT_ALL ,
     1    IMONIT_BUDJET ,
     1    IMONIT_FLUX_MRI , IMONIT_ALMA
!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_RAD_SHRT_ADD 00/09 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_BUDJET ) THEN
!
! 下向き短波
!
        CMARK  = 'DSWBL'           
        CTITLE = 'RAD SHRT DOWNWARD AT SIB TOP'    
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
! 上向き短波
!
        CMARK = 'USWBL'           
        CTITLE = 'RAD SHRT UPWARD   AT SIB TOP'   ! 
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
        CMARK = 'HSRDNTP' 
        CTITLE = 'RAD SHRT DOWNWARD AT SIB TOP (NOSNOW)'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK = 'HSRUNTP' 
        CTITLE = 'RAD SHRT UPWARD   AT SIB TOP (NOSNOW)'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK = 'HSRDSTP' 
        CTITLE = 'RAD SHRT DOWNWARD AT SIB TOP (SNOW)'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK = 'HSRUSTP' 
        CTITLE = 'RAD SHRT UPWARD   AT SIB TOP (SNOW)'
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF 
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALMA ) THEN
        CMARK  = 'A_SWDN'           
        CTITLE = 'ALMA + Swnet and / Albedo (Shortwave down)'    
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK = 'A_SWUP'           
        CTITLE = 'ALMA - Swnet and * Albedo (Shortwave up)'    
        CUNIT  = 'W/M**2'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
      ENDIF
!
      RETURN
      END

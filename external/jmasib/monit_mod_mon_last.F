! 
! 注意
!   ・月末のステップで、全ての計算が終わった後に行なう。
!     雪や土壌の計算は今の月のパラメタを使って計算しているが、
!     このステップの計算後の値は次の月の頭の値であるので、
!     月変わり時に LAI 等が変わったことに対する対処をしないと
!     いけないのである。
!
      SUBROUTINE MONIT_MOD_MON_LAST (
     I   NMAX             , NTYP             , 
     I   MON              , MON_NEXT         ,
     I   CVR_SNW_NEXT     , 
     I   DEP_SNOW_NEXT    , 
     I   TMP_CNP_NOS_NEXT , FR_WTR_CNP_NOS_NEXT , FR_GLA_CNP_NOS_NEXT ,        
     I   TMP_CNP_SNW_NEXT , FR_WTR_CNP_SNW_NEXT , FR_GLA_CNP_SNW_NEXT ,
     I   TMP_GRSK_NEXT    , FR_WTR_GRS_NEXT     , FR_GLA_GRS_NEXT    ) 
! 
      USE SIBPRM , ONLY :
     1   IDX  , 
     1   ISN  ,  
     1   ITYP  
!
      USE SIBVEG , ONLY :
     1           XLAI_mo                    , !              (ITYP,IMON,ICG)
     1           D_ZERO_mo                  , ! 零面変位 d   (ITYP,IMON)  
     1           Z2_mo                      , 
     1           Z1_mo                 
!
      USE SIBCON , ONLY :
     1           TMP_FREZ_c      ,
     1           HL_W2G_c        ,
     1           CAP_GLA_c       ,
     1           CAP_WTR_c       ,
     1           CAP_LAI_c       ,
     1           RHO_WTR_c       ,
     1           H2O_SAT_CNP_DEP_c ,
     1           H2O_SAT_GRS_DEP_c 

!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB_NEXT, ICN, ICS ,
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif
!
      IMPLICIT NONE
!
      REAL(8),INTENT(IN) :: CVR_SNW_NEXT        (IDX)  ! 雪面被覆率
!
      REAL(8),INTENT(IN) :: TMP_CNP_NOS_NEXT    (IDX)  ! キャノピー温度
      REAL(8),INTENT(IN) :: TMP_CNP_SNW_NEXT    (IDX)  ! キャノピー温度
      REAL(8),INTENT(IN) :: FR_WTR_CNP_SNW_NEXT (IDX)  ! キャノピー水割合
      REAL(8),INTENT(IN) :: FR_GLA_CNP_SNW_NEXT (IDX)  ! キャノピー氷割合
      REAL(8),INTENT(IN) :: FR_WTR_CNP_NOS_NEXT (IDX)  ! キャノピー水割合
      REAL(8),INTENT(IN) :: FR_GLA_CNP_NOS_NEXT (IDX)  ! キャノピー氷割合
      REAL(8),INTENT(IN) :: TMP_GRSK_NEXT       (IDX)  ! 地面温度
      REAL(8),INTENT(IN) :: FR_WTR_GRS_NEXT     (IDX)  ! 地面水割合
      REAL(8),INTENT(IN) :: FR_GLA_GRS_NEXT     (IDX)  ! 地面氷割合
!
      REAL(8),INTENT(IN) :: DEP_SNOW_NEXT (IDX,ISN+1)   ! 積雪層の厚さ(m)
!
      INTEGER,INTENT(IN) :: NMAX 
      INTEGER,INTENT(IN) :: MON
      INTEGER,INTENT(IN) :: MON_NEXT
      INTEGER,INTENT(IN) :: NTYP(IDX) 
!
! WORK
!     
      INTEGER     :: I
!
      REAL(8)     :: SDEP
      REAL(8)     :: XZ 
      REAL(8)     :: FAC 
      REAL(8)     :: SDEP_NEXT
      REAL(8)     :: XZ_NEXT
      REAL(8)     :: FAC_NEXT
      REAL(8)     :: XLAI_CNP     
      REAL(8)     :: XLAI_CNP_NEXT
      REAL(8)     :: XLAI_GRS     
      REAL(8)     :: XLAI_GRS_NEXT
      REAL(8)     :: ADD_LAI_CNP_SNW (IDX) 
      REAL(8)     :: ADD_LAI_CNP_NOS (IDX) 
      REAL(8)     :: ADD_LAI_GRS     (IDX) 
      REAL(8),PARAMETER  :: ZERO = 0.D0 

!
      REAL(8),PARAMETER  :: ONE = 1.D0 
      REAL(8)            :: WORK(IDX) 
      CHARACTER(7)       :: CMARK
!
      LOGICAL,SAVE :: LFIRST
      DATA LFIRST / .TRUE. /
!
      IF ( LFIRST ) THEN 
        WRITE(6,*) 'MONIT_MOD_MON_LAST 00/07 POCHAKA'
        write(6,*) cvr_snw_next(1)
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ============================
! >> 補正量・変化量のモニタ <<
! ============================
!
!  面積などの変化率 （全格子に対する割合、全格子・全LAI に対する割合 ）
!
      DO I = 1,NMAX
        XZ        = Z2_mo  ( NTYP(I) , MON      ) 
     1              - Z1_mo  ( NTYP(I) , MON      )
        XZ_NEXT   = Z2_mo  ( NTYP(I) , MON_NEXT )
     1              - Z1_mo  ( NTYP(I) , MON_NEXT )
!
        SDEP      = MIN ( DEP_SNOW_NEXT(I,5)
     1                       - Z1_mo  ( NTYP(I) , MON     )
     1                   , XZ     * 0.95D0 )
        SDEP      = MAX ( SDEP , ZERO )
        SDEP_NEXT = MIN ( DEP_SNOW_NEXT(I,5) 
     1                       - Z1_mo  ( NTYP(I) , MON_NEXT )
     1                   , XZ_NEXT * 0.95D0 )
        SDEP_NEXT = MAX ( SDEP_NEXT , ZERO )
!
        FAC       = 1.0D0 - SDEP      / XZ    
        FAC_NEXT  = 1.0D0 - SDEP_NEXT / XZ_NEXT    
!
        XLAI_CNP       = XLAI_mo ( NTYP(I) , MON      , 1 ) 
        XLAI_CNP_NEXT  = XLAI_mo ( NTYP(I) , MON_NEXT , 1 ) 
!
        XLAI_GRS       = XLAI_mo ( NTYP(I) , MON      , 2 ) 
        XLAI_GRS_NEXT  = XLAI_mo ( NTYP(I) , MON_NEXT , 2 ) 
!
        ADD_LAI_CNP_SNW (I) 
     1       = XLAI_CNP_NEXT  * FAC_NEXT - XLAI_CNP * FAC
        ADD_LAI_CNP_NOS (I) 
     1       = XLAI_CNP_NEXT - XLAI_CNP
        ADD_LAI_GRS (I) 
     1       = XLAI_GRS_NEXT  - XLAI_GRS
      ENDDO
!
!  雪ありキャノピーでの水の変化量 
!
      DO I=1,NMAX 
        WORK(I) 
     1    = ADD_LAI_CNP_SNW(I) * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1      * ( FR_WTR_CNP_SNW_NEXT (I) + FR_GLA_CNP_SNW_NEXT (I) )
      ENDDO
      CMARK = 'WADSMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICS, ONE , ONE )  
!
!  雪無しキャノピーでの水の変化量
!
      DO I=1,NMAX
        WORK(I) 
     1    = ADD_LAI_CNP_NOS(I) * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1      * ( FR_WTR_CNP_NOS_NEXT (I) + FR_GLA_CNP_NOS_NEXT (I) )
      ENDDO
      CMARK = 'WADNMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICN, ONE , ONE )  
!
!  雪無し下草での水の変化量
!
      DO I=1,NMAX
        WORK(I) 
     1    = ADD_LAI_GRS(I) * RHO_WTR_c * H2O_SAT_GRS_DEP_c 
     1      * ( FR_WTR_GRS_NEXT (I) + FR_GLA_GRS_NEXT (I) )
      ENDDO
      CMARK = 'WADGMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICN, ONE , ONE )  
!
!  雪ありキャノピーでのエネルギーの変化量
!
      DO I=1,NMAX
        WORK(I) 
     1    = (    CAP_LAI_c 
     1         + CAP_WTR_c * FR_WTR_CNP_SNW_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1         + CAP_GLA_c * FR_GLA_CNP_SNW_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1       )
     1       * ADD_LAI_CNP_SNW(I) 
     1       * ( TMP_CNP_SNW_NEXT(I) - TMP_FREZ_c )
     1      - HL_W2G_c * FR_GLA_CNP_SNW_NEXT (I) 
     1                 * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1                 * ADD_LAI_CNP_SNW(I) 
      ENDDO
      CMARK = 'IADSMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICS, ONE , ONE )  
!
!  雪無しキャノピーでのエネルギーの変化量
!
      DO I=1,NMAX
        WORK(I) 
     1    = (    CAP_LAI_c 
     1         + CAP_WTR_c * FR_WTR_CNP_NOS_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1         + CAP_GLA_c * FR_GLA_CNP_NOS_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1      )
     1       * ADD_LAI_CNP_NOS(I) 
     1       * ( TMP_CNP_NOS_NEXT(I) - TMP_FREZ_c )
     1      - HL_W2G_c * FR_GLA_CNP_NOS_NEXT (I) 
     1                 * RHO_WTR_c * H2O_SAT_CNP_DEP_c 
     1                 * ADD_LAI_CNP_NOS(I) 
      ENDDO
      CMARK = 'IADNMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICN, ONE , ONE )  
!
!  雪無し地面・下草でのエネルギーの変化量
!      ... マットに貯めてるエネルギーは収支には関係無いことに注意。
!          あれは、GRSK を大きく変化させないための重しにすぎない。
!          （その熱は全て地中に抜ける）
! 
      DO I=1,NMAX
        WORK(I) 
     1    = (    CAP_LAI_c 
     1         + CAP_WTR_c * FR_WTR_GRS_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_GRS_DEP_c 
     1         + CAP_GLA_c * FR_GLA_GRS_NEXT (I) 
     1                     * RHO_WTR_c * H2O_SAT_GRS_DEP_c 
     1      )
     1       * ADD_LAI_GRS(I) 
     1       * ( TMP_GRSK_NEXT(I) - TMP_FREZ_c )
     1      - HL_W2G_c * FR_GLA_GRS_NEXT (I) 
     1                 * RHO_WTR_c * H2O_SAT_GRS_DEP_c 
     1                 * ADD_LAI_GRS(I) 
      ENDDO
      CMARK = 'IADGMI' 
      CALL MONIT_SIB_NEXT ( CMARK , WORK ,ICN, ONE , ONE )  
!
      ENDIF
#endif
!       
      RETURN
      END
!==================================================================
      SUBROUTINE MR_MONIT_MOD_MON_LAST 
!
! SNOWMOD_MAIN で出力するデバッグ用変数の登録作業
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : 
     1                       IMONIT_LEVEL , IMONIT_ALL
#endif

!
      IMPLICIT NONE
!
      CHARACTER( 7)  :: CMARK 
      CHARACTER(32)  :: CTITLE 
      CHARACTER(13)  :: CUNIT
!
      LOGICAL,SAVE ::  LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_MOD_MON_LAST 00/07 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN
!
! ---------------------------------------------------------
!
      CMARK  = 'WADSMI' 
      CTITLE = 'H2O ADD ON SNW-CANOPY BY MONTH INI'
      CUNIT  = 'KG/S/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'WADNMI' 
      CTITLE = 'H2O ADD ON NOS-CANOPY BY MONTH INI'
      CUNIT  = 'KG/S/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )

      CMARK  = 'WADGMI' 
      CTITLE = 'H2O ADD ON GRASS BY MONTH INI'
      CUNIT  = 'KG/S/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IADSMI' 
      CTITLE = 'ENERGY WITH H2O ADD ON SNW-CANOPY BY MONTH INI'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IADNMI' 
      CTITLE = 'ENERGY WITH H2O ADD ON NO-CANOPY BY MONTH INI'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      CMARK  = 'IADGMI' 
      CTITLE = 'ENERGY WITH H2O ADD ON GRASS BY MONTH INI'
      CUNIT  = 'W/M**2'
      CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      RETURN
      END

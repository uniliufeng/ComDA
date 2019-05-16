      SUBROUTINE SNOW_RHO_INTG (
     I    NMAX  , NMAX_01 , NMAX_2 , NMAX_3 , NMAX_4 , 
     I    RHO_SNOW_INV    , 
     I    DELT_SNOW       , 
!    I    WTR_SNOW        , GLA_SNOW      , TMP_REF   ,  
     I    WTR_SNOW        , GLA_SNOW      , TMP_SNOW  ,  
     O    RHO_SNOW_INV_NEXT )
!
      USE SIBPRM , ONLY :
     1    IDX , 
     1    ISN 
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: NMAX 
      INTEGER,INTENT(IN)  :: NMAX_01(3) 
      INTEGER,INTENT(IN)  :: NMAX_2 (3) 
      INTEGER,INTENT(IN)  :: NMAX_3 (3) 
      INTEGER,INTENT(IN)  :: NMAX_4 (3) 
      REAL(8),INTENT(IN)  :: RHO_SNOW_INV     (IDX,ISN)
      REAL(8),INTENT(IN)  :: WTR_SNOW         (IDX,ISN)
      REAL(8),INTENT(IN)  :: GLA_SNOW         (IDX,ISN)
!     REAL(8),INTENT(IN)  :: TMP_REF          (IDX)
      REAL(8),INTENT(IN)  :: TMP_SNOW         (IDX,ISN)
      REAL(8),INTENT(IN)  :: DELT_SNOW
!
! OUTPUT
!
      REAL(8),INTENT(OUT) :: RHO_SNOW_INV_NEXT(IDX,ISN)
!
! WORK
!
      REAL(8) :: WGT       (IDX)       ! ÇªÇÃëwÇÃíÜâõÇ≈ÇÃê√ïXà≥ÅiíPà  PaÅj
      
!
      LOGICAL,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_RHO_INTG 99/10/15 POCHAKA'
        LFIRST = .FALSE.
      ENDIF
!
      CALL RESET ( WGT , IDX )
      CALL RESET ( RHO_SNOW_INV_NEXT , IDX*ISN )
!
      CALL SNOW_RHO_INTG_SUB (
     I    NMAX_01(2)         , NMAX               , 1         , 
     I    RHO_SNOW_INV(1,1)  , 
     I    DELT_SNOW          , 
     I    WTR_SNOW(1,1)      , GLA_SNOW(1,1)      , TMP_SNOW(1,1)  ,  
     M    WGT                , 
     O    RHO_SNOW_INV_NEXT(1,1) )
!
      CALL SNOW_RHO_INTG_SUB (
     I    NMAX_2 (2)         , NMAX               , 2         , 
     I    RHO_SNOW_INV(1,2)  , 
     I    DELT_SNOW          , 
     I    WTR_SNOW(1,2)      , GLA_SNOW(1,2)      , TMP_SNOW(1,2) ,  
     M    WGT                , 
     O    RHO_SNOW_INV_NEXT(1,2) )
!
      CALL SNOW_RHO_INTG_SUB (
     I    NMAX_3 (2)         , NMAX               , 3         , 
     I    RHO_SNOW_INV(1,3)  , 
     I    DELT_SNOW          , 
     I    WTR_SNOW(1,3)      , GLA_SNOW(1,3)      , TMP_SNOW(1,3)  ,  
     M    WGT                , 
     O    RHO_SNOW_INV_NEXT(1,3) )
!
      CALL SNOW_RHO_INTG_SUB (
     I    NMAX_4(2)         , NMAX               , 4          , 
     I    RHO_SNOW_INV(1,4)  , 
     I    DELT_SNOW          , 
     I    WTR_SNOW(1,4)      , GLA_SNOW(1,4)      , TMP_SNOW(1,4)  ,  
     M    WGT                , 
     O    RHO_SNOW_INV_NEXT(1,4) )
!
      RETURN
      END
!===========================================================
      SUBROUTINE SNOW_RHO_INTG_SUB (
     I    NSTR            , NEND          , KLEV      , 
     I    RHO_SNOW_INV    , 
     I    DELT_SNOW       , 
     I    WTR_SNOW        , GLA_SNOW      , TMP_SNOW   ,  
     M    WGT             , 
     O    RHO_SNOW_INV_NEXT )
!
      USE SIBPRM , ONLY :
     1    IDX  
!
      USE SIBCON , ONLY :
     1    TMP_FREZ_c  , 
     1    GRAV_c      ,
     1    RHO_SNOW_INTG_THR_c ,
     1    RHO_SNOW_INTG2_THR_c,
     1    RHO_GLA_c,
     1    RHO_WTR_c
!
!      USE COM_RUNCONF_SIB0109 , ONLY :
!     1   JCN_CALC_SNOW_ETA
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN)  :: NSTR
      INTEGER,INTENT(IN)  :: NEND
      INTEGER,INTENT(IN)  :: KLEV
!
      REAL(8),INTENT(IN)  :: RHO_SNOW_INV     (IDX)
      REAL(8),INTENT(IN)  :: WTR_SNOW         (IDX)
      REAL(8),INTENT(IN)  :: GLA_SNOW         (IDX)
      REAL(8),INTENT(IN)  :: TMP_SNOW          (IDX)
      REAL(8),INTENT(IN)  :: DELT_SNOW
!
! INOUT     
!
      REAL(8),INTENT(INOUT) :: WGT     (IDX)  ! ÇªÇÃëwíÜâõïîÇ≈ÇÃê√ïXà≥ÅiPaÅj
!
! OUTPUT
!
      REAL(8),INTENT(INOUT) :: RHO_SNOW_INV_NEXT(IDX)
!
! WORK
!
      INTEGER :: I
      REAL(8) :: ETA       (IDX)
      REAL(8) :: RHO_SNOW
      REAL(8) :: WGT_LEV 
!     REAL(8),PARAMETER :: RHO_THR1 =  78.
!     REAL(8),PARAMETER :: RHO_THR2 =  0.0025     ! = 1/400
      REAL(8),PARAMETER :: FAC_1A   =    1.42D5     
      REAL(8),PARAMETER :: FAC_1R   =    0.0663D0
      REAL(8),PARAMETER :: FAC_1T   =  - 0.0958D0
      REAL(8),PARAMETER :: FAC_2A   =    3.44D6 
      REAL(8),PARAMETER :: FAC_2R   =    0.0253D0
      REAL(8),PARAMETER :: FAC_2T   =  - 0.0958D0  
!
!
!-------  ä‹êÖÇÃîSê´ó¶Ç÷ÇÃâeãøÇ™î≤ÇØÇƒÇ¢ÇΩÅiHIRAI & Oh'izumi)
!     integer,parameter  :: JCN_CALC_SNOW_ETA = 1 
      REAL(8) :: D_a, W_vol,A_coef
      REAL(8),parameter :: fac_A_W_vol = -0.092
!---------------------------------------------------------------------
!
      LOGICAL,SAVE :: LFIRST 
!
      DATA LFIRST /.TRUE./
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'SNOW_RHO_INTG_SUB 99/10/15 POCHAKA' , NSTR,NEND
        LFIRST = .FALSE.
      ENDIF
!
! ê·ëwÇÃà≥èkîSê´åWêîÇ…ëÃêœä‹êÖó¶ÇÃåWêî(A(Wvol))ÇÃåWêîÇ™î≤ÇØÇƒÇ¢ÇΩ
! ï ç˚46çÜP59(7)Ç…è]Ç¢ïsãÔçáÇíºÇ∑ÅB   Hirai 2004/05/18 
!
        DO I=NSTR,NEND
            RHO_SNOW = 1.0D0 / RHO_SNOW_INV(I)
          IF ( RHO_SNOW .LT. RHO_SNOW_INTG_THR_c ) THEN   ! êVê·
            ETA(I) = FAC_1A * DEXP(   FAC_1R * RHO_SNOW
     1                        + FAC_1T * ( TMP_SNOW(I) - TMP_FREZ_c ) )
          ELSE IF ( RHO_SNOW .LT. RHO_SNOW_INTG2_THR_c) then       ! ãåê·1
!
            D_a=(1.0/rho_snow - 1.0/rho_gla_c)*(gla_snow(i)+wtr_snow(i))
            W_vol=(wtr_snow(i)/rho_wtr_c)/((wtr_snow(i)/rho_wtr_c)+D_a)
            A_coef=DEXP( fac_A_W_vol * W_vol * 100.D0 )
!                                       W_volÇÃíPà ÇÕ[%]Ç≈çsÇ§

            ETA(I) =A_coef* FAC_2A * DEXP(   FAC_2R * RHO_SNOW
     1                        + FAC_2T * ( TMP_SNOW(I) - TMP_FREZ_c ) )
          ELSE                                                     ! ãåê·2
            ETA(I) = FAC_2A * DEXP(   FAC_2R * RHO_SNOW
     1                        + FAC_2T * ( TMP_SNOW(I) - TMP_FREZ_c ) )
          ENDIF
        ENDDO

!
      DO I=NSTR,NEND
        WGT_LEV = ( WTR_SNOW(I) + GLA_SNOW(I) ) * GRAV_c 
!
        RHO_SNOW_INV_NEXT(I) 
     1          = RHO_SNOW_INV(I)  
     1            - ( WGT(I) + WGT_LEV / 2 ) / ETA(I)     
     1              * RHO_SNOW_INV(I)          
     1              * DELT_SNOW
!
        WGT(I) = WGT(I) + WGT_LEV 
!
#ifdef SIB_DEBUG
        IF (       1/RHO_SNOW_INV_NEXT(I) .LT. 30. 
     1        .OR. 1/RHO_SNOW_INV_NEXT(I) .GT. 999 ) THEN
          write(6,*) 'SNOW_RHO_INTG : WARNING ' , 
     1               ' RHO_SNOW_INV_NEXT IS NOT GOOD. ' , 
     1               ' I K NEXT_INV INV RHO ' ,
     1               I , KLEV , 
     1               RHO_SNOW_INV_NEXT(I) , RHO_SNOW_INV(I),
     1               1/RHO_SNOW_INV(I) , 1/RHO_SNOW_INV_NEXT(I) 
         STOP 999
        ENDIF
#endif
!
      ENDDO
!
      RETURN
      END

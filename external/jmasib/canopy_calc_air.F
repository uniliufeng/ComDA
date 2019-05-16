      SUBROUTINE CANOPY_CALC_AIR (
     I    NSTR      , NEND       , MON       , NTYP       , CVR_NorS  , 
#ifdef MONYOS
     I    DELT_CNP  , INSFLAG    ,
#endif
     I    TA        , EA         , TG         , 
     I    TMP_REF   , Q_REF      , Z_REF      , P_SFC     , 
     I    U_REF     , V_REF      , 
     I    Z0_msx    , D_ZERO_msx , FAC_CORTMP , U2        , 
     O    TS        , QS         , US         , VS  ) 
!
      USE SIBPRM , ONLY : IDX 
      USE COM_RUNCONF_SIB0109 , ONLY : 
     1    JCN_SIB_UPPER_MODEL 
!
      USE SIBVEG , ONLY :
     1           Z2_mo
!
#ifdef MONYOS
      USE SIB_MONIT , ONLY : MONIT_SIB , ICN, ICS ,
     1                       imonit_level , 
     1                       imonit_mean_mri ,
     1                       imonit_snap_npd
      USE COM_STEP_SIB0109 , ONLY : ICN_SIB0109_MONITOR
#endif 
!
      IMPLICIT NONE
!
      REAL(8), INTENT(INOUT)  :: TS(IDX)
      REAL(8), INTENT(INOUT)  :: QS(IDX)
      REAL(8), INTENT(INOUT)  :: US(IDX)
      REAL(8), INTENT(INOUT)  :: VS(IDX)
!
#ifdef MONYOS
      REAL(8), INTENT(IN)     :: DELT_CNP 
      INTEGER, INTENT(IN)     :: INSFLAG 
#endif
      INTEGER, INTENT(IN)     :: MON 
      INTEGER, INTENT(IN)     :: NSTR
      INTEGER, INTENT(IN)     :: NEND
      INTEGER, INTENT(IN)     :: NTYP      (IDX)
      REAL(8), INTENT(IN)     :: CVR_NorS  (IDX)
      REAL(8), INTENT(IN)     :: EA        (IDX)
      REAL(8), INTENT(IN)     :: TG        (IDX)
      REAL(8), INTENT(IN)     :: TA        (IDX)
      REAL(8), INTENT(IN)     :: P_SFC     (IDX)
      REAL(8), INTENT(IN)     :: Z0_msx    (IDX)
      REAL(8), INTENT(IN)     :: D_ZERO_msx(IDX)
      REAL(8), INTENT(IN)     :: TMP_REF   (IDX)
      REAL(8), INTENT(IN)     :: FAC_CORTMP(IDX)
      REAL(8), INTENT(IN)     :: Z_REF     (IDX)
      REAL(8), INTENT(IN)     :: U_REF     (IDX)
      REAL(8), INTENT(IN)     :: V_REF     (IDX)
      REAL(8), INTENT(IN)     :: Q_REF     (IDX)
      REAL(8), INTENT(IN)     :: U2        (IDX)
!
      REAL(8)    :: QA 
      REAL(8)    :: FMLTH 
      REAL(8)    :: FMLTM
      REAL(8)    :: Z0S
      REAL(8)    :: ZS
!
      REAL(8) ,parameter :: vs2_small = 1.D-10
!
      INTEGER    :: I  
!
#ifdef MONYOS
      CHARACTER(7) :: CMARK
      REAL(8),PARAMETER  :: ONE = 1.D0
#endif
!
      IF ( JCN_SIB_UPPER_MODEL .EQ. 1 ) THEN         ! mj98 
#include "indep_fixedform"
        DO I=NSTR,NEND
          QA         =0.622*EA(I)/P_SFC(I)
          IF(Z2_mo(NTYP(I),MON).GE.2.0) THEN           ! for AMIP 9807 Ohizumi
            TS(I) = TS(I) + TA(I) * CVR_NorS(I)
            QS(I) = QS(I) + QA    * CVR_NorS(I)
          ELSE
            FMLTH   = LOG( MAX(Z0_msx(I),2.0-D_ZERO_msx(I))/Z0_msx(I) )
     1              / LOG(    (Z_REF(I)     -D_ZERO_msx(I))/Z0_msx(I) )
            TS(I)=  TS(I)  
     1              + (TA(I)+(TMP_REF(I)*FAC_CORTMP(I)-TA(I))*FMLTH)
     1                /(1.0+(FAC_CORTMP(I)-1.0)*FMLTH) 
     1                * CVR_NorS(I) 
            QS(I)= QS(I) + ( QA +(Q_REF(I)-QA)*FMLTH ) * CVR_NorS(I) 
          END IF
!
          IF(Z2_mo(NTYP(I),MON).GE.10.0) THEN
            US(I)= US(I)  
     1             + U_REF(I)*U2(I)
     1                       /(SQRT(U_REF(I)**2+V_REF(I)**2)+vs2_small)
     1               * CVR_NorS(I) 
            VS(I)= VS(I) 
     1             + V_REF(I)*U2(I)
     1                       /(SQRT(U_REF(I)**2+V_REF(I)**2)+vs2_small)
     1               * CVR_NorS(I) 
          ELSE
            FMLTM  = LOG( MAX(Z0_msx(I),10.0-D_ZERO_msx(I))/Z0_msx(I))
     1              /LOG(    (Z_REF(I)      -D_ZERO_msx(I))/Z0_msx(I))
            US(I) = US(I) + U_REF(I)*FMLTM * CVR_NorS(I) 
            VS(I) = VS(I) + V_REF(I)*FMLTM * CVR_NorS(I) 
          END IF
        ENDDO
!
      ELSE                     ! –{’¡ GSM 
!
        Z0S = 0.01
        ZS = 2.D0
#include "indep_fixedform"
        DO I=NSTR,NEND
          QA         =0.622*EA(I)/(P_SFC(I)-0.378*EA(I))
          FMLTH      = LOG(ZS/Z0S)
     1               / LOG(Z_REF(I)/Z0S)
          TS(I) = TS(I) 
     1           + (TG(I)+(TMP_REF(I)*FAC_CORTMP(I)-TG(I))*FMLTH)
     1             /(1.0+(FAC_CORTMP(I)-1.0)*FMLTH)
     1             * CVR_NorS(I) 
          QS(I) = QS(I) 
     1            + ( QA +(Q_REF(I)-QA)*FMLTH )  
     1             * CVR_NorS(I) 
!
          IF(Z2_mo(NTYP(I),MON).GE.10.0) THEN
            US(I) = US(I)  
     1              + U_REF(I)*U2(I)
     1                        /(SQRT(U_REF(I)**2+V_REF(I)**2)+vs2_small)
     1                * CVR_NorS(I) 
            VS(I) = VS(I) 
     1              + V_REF(I)*U2(I)
     1                        /(SQRT(U_REF(I)**2+V_REF(I)**2)+vs2_small)
     1                * CVR_NorS(I) 
          ELSE
            FMLTM =LOG( MAX(Z0_msx(I),10.0 -D_ZERO_msx(I))/Z0_msx(I) )
     1            /LOG(           (Z_REF(I)-D_ZERO_msx(I))/Z0_msx(I) )
            US(I) = US(I) 
     1              + U_REF(I)*FMLTM * CVR_NorS(I) 
            VS(I) = VS(I) 
     1              + V_REF(I)*FMLTM * CVR_NorS(I) 
          ENDIF
        ENDDO
      ENDIF
!
#ifdef MONYOS
      IF ( IMONIT_LEVEL .GE. IMONIT_MEAN_MRI ) THEN  ! time-mean 
        CMARK = 'Z0' 
        IF ( INSFLAG .EQ. 0 ) THEN
          CALL MONIT_SIB ( CMARK , Z0_msx ,ICN, DELT_CNP , ONE )      
        ELSE
          CALL MONIT_SIB ( CMARK , Z0_msx ,ICS, DELT_CNP , ONE )      
        ENDIF
      ENDIF
!
      IF (       IMONIT_LEVEL .GE. IMONIT_SNAP_NPD 
     1     .AND. ICN_SIB0109_MONITOR .EQ. 1         ) THEN  ! snap 
        CMARK = 'Z0' 
        IF ( INSFLAG .EQ. 0 ) THEN
          CALL MONIT_SIB ( CMARK , Z0_msx ,ICN, DELT_CNP , ONE )      
        ELSE
          CALL MONIT_SIB ( CMARK , Z0_msx ,ICS, DELT_CNP , ONE )      
        ENDIF
      ENDIF

#endif
!
      RETURN
      END


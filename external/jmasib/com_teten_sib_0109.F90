MODULE COM_TETEN_sib_0109

  IMPLICIT NONE

!--- COMMON/CTETEN/
  REAL(8),SAVE :: TABLE(25000)

!--- COMMON/DTETEN/
  REAL(8),SAVE :: DTABLE(25000)

!--- COMMON/CLATENT/
  REAL(8),SAVE :: CTABLE(25000)
  REAL(8),SAVE :: DCL
  REAL(8),SAVE :: TEMP0
  REAL(8),SAVE :: TEMPI

!--- COMMON/RLQIC/
  REAL(8),SAVE :: RTABLE(25000)
  REAL(8),SAVE :: TLI1
  REAL(8),SAVE :: TLI2

!--- COMMON/COMEVP/
  REAL(8),SAVE :: CEV(4001)
  REAL(8),SAVE :: DFW
  REAL(8),SAVE :: RDFW

!-------------------------------------------------

  CONTAINS

    SUBROUTINE COM_TETEN_sib_0109_INI

    USE PRM_PHCONST
    USE COM_RUNCONF

    REAL(8),PARAMETER :: XB = 21.18123D0
    REAL(8),PARAMETER :: XC = 5418.0D0
    REAL(8),PARAMETER :: B  = 19.480254D0
    REAL(8),PARAMETER :: C  = 4304.412D0
    REAL(8),PARAMETER :: BI = 23.684812D0
    REAL(8),PARAMETER :: CI = 5803.3203D0
    REAL(8),PARAMETER :: TSAT  = 29.55D0
    REAL(8),PARAMETER :: TSATI = 7.85D0

    LOGICAL,SAVE :: FIRST = .TRUE.

    REAL(8) :: DTEMP
    REAL(8) :: DICE
    REAL(8) :: HICE
    REAL(8) :: X
    REAL(8) :: RR
    REAL(8) :: TBL1
    REAL(8) :: DTBL1
    REAL(8) :: TBL2
    REAL(8) :: DTBL2
    REAL(8) :: FWMX
    REAL(8) :: FWMN
    REAL(8) :: FW

    INTEGER :: I
    INTEGER :: IFWM

    IF ( .NOT. FIRST ) THEN
      WRITE(6,*) 'TETEN : NOT FIRST TIME CALL , SO RETURN ONLY '
      RETURN
    ENDIF

    FIRST = .FALSE.    
!
    TEMP0 = 273.15
    TEMPI = 233.15
    DTEMP = TEMP0-TEMPI
    TLI1  = TEMPI
    TLI2  = TEMP0
!P  IF(ICE.EQ.1) THEN
    IF(JCNTTICE.EQ.1) THEN
      DICE = 3.33E5
    ELSE
      DICE = 0.0
    ENDIF
    HICE = HL + DICE
!  DL/DT
    DCL = -DICE/DTEMP
!X  CLBYCP = HL/CP
!X  CLBYCPI = HICE/CP
!P  IF(ICE.EQ.0) THEN
    IF(JCNTTICE.EQ.0) THEN
      DO I = 1,25000
        X = 123.2D0 + 0.01D0*I
        TABLE(I) = 0.622*DEXP(B-C/(X-TSAT))
        DTABLE(I) = TABLE(I)*C/(X-TSAT)**2
        CTABLE(I) = HL
        RTABLE(I) = 1.
      END DO
    ELSE
      DO I = 1,25000
        X = 123.2D0 + 0.01D0*I
        IF(X.GE.TEMP0) THEN
          TABLE(I) = 0.622*DEXP(B-C/(X-TSAT))
          DTABLE(I) = TABLE(I)*C/(X-TSAT)**2
          CTABLE(I) = HL
          RTABLE(I) = 1.
        ELSEIF(X.LE.TEMPI) THEN
          TABLE(I) = 0.622*DEXP(BI-CI/(X-TSATI))
          DTABLE(I) = TABLE(I)*CI/(X-TSATI)**2
          CTABLE(I) = HICE
          RTABLE(I) = 0.
        ELSE
          RR = (TEMP0-X)/DTEMP
          CTABLE(I) = HL*(1.0-RR) + HICE*RR
          TBL1 = 0.622*DEXP(B-C/(X-TSAT))
          DTBL1 = TBL1*C/(X-TSAT)**2
          TBL2 = 0.622*DEXP(BI-CI/(X-TSATI))
          DTBL2 = TBL2*CI/(X-TSATI)**2
          TABLE(I)  = TBL1*(1.D0-RR)+TBL2*RR
          DTABLE(I) = DTBL1*(1.D0-RR)+DTBL2*RR+(TBL1-TBL2)/DTEMP
          RTABLE(I) = 1.D0-RR
! ##                         3-JI KANSUU : RTABLE(I) = TT*TT*(3.-2.*TT)
! ##                                       WHERE   TT = 1.D0-RR
        ENDIF
      END DO
    ENDIF
    FWMX = 5.0
    FWMN = 0.0
    IFWM = 4001
    DFW = (FWMX-FWMN)/(IFWM-1.)
    RDFW=1./DFW
!..#include "nofval_freeform"
    DO I = 1,IFWM
      FW = FWMN + DFW*(I-1.)
      CEV(I) = 8.*RGAS*(1.6+23.2*(FW)**0.167)*(FW)**0.467
    END DO
    CEV(1) = 0.0

  END SUBROUTINE COM_TETEN_sib_0109_INI

END MODULE COM_TETEN_sib_0109

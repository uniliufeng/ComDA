      SUBROUTINE SETBACK_TO_ATM  (
     I   ILIST       , NMAX           , 
     I   GMT_OUT     , GMQ_OUT        , GMUV_OUT      ,  
     I   RAD_LONG_SIB_2_ATM_SIB       , 
     O   TMTX2L_PHY  , QMTX2L_PHY     , UMTX2L_PHY    , VMTX2L_PHY    , 
     O   RAD_LONG_SIB_2_ATM_PHY       )
!
      USE PRM , ONLY :
     1   IJPHY  
!
      USE SIBPRM , ONLY : 
     1   IDX   
!
      IMPLICIT NONE
!
! INPUT
!
      INTEGER,INTENT(IN)  :: NMAX
!     INTEGER,INTENT(IN)  :: NMAX_01(3)
      INTEGER,INTENT(IN)  :: ILIST    ( IDX )
!
!     REAL(8),INTENT(IN)  :: CVR_SNW  ( IDX )
!     REAL(8),INTENT(IN)  :: CVR_NOS  ( IDX )
!
      REAL(8),INTENT(IN)  :: GMT_OUT  ( IDX,2 )
      REAL(8),INTENT(IN)  :: GMQ_OUT  ( IDX,2 )
      REAL(8),INTENT(IN)  :: GMUV_OUT ( IDX,3 )
      REAL(8),INTENT(IN)  :: RAD_LONG_SIB_2_ATM_SIB (IDX)
!
! OUTPUT
!
      REAL(8),INTENT(OUT)  :: TMTX2L_PHY ( IJPHY )
      REAL(8),INTENT(OUT)  :: QMTX2L_PHY ( IJPHY )
      REAL(8),INTENT(OUT)  :: UMTX2L_PHY ( IJPHY )
      REAL(8),INTENT(OUT)  :: VMTX2L_PHY ( IJPHY )
      REAL(8),INTENT(OUT)  :: RAD_LONG_SIB_2_ATM_PHY (IJPHY)
!
! WORK
!
      INTEGER   :: I
      INTEGER   :: J
!
      IF ( NMAX.LT.1 ) RETURN
!
! これでよいかどうか、online にするときは要チェック!!! 
!
      TMTX2L_PHY (:) = 0.
      QMTX2L_PHY (:) = 0.
      UMTX2L_PHY (:) = 0.
      VMTX2L_PHY (:) = 0.
      RAD_LONG_SIB_2_ATM_PHY (:) = 0.
!
      DO I=1,NMAX 
        J = ILIST(I)
        TMTX2L_PHY (J)   = GMT_OUT (I,2) 
        QMTX2L_PHY (J)   = GMQ_OUT (I,2) 
        UMTX2L_PHY (J)   = GMUV_OUT(I,2) 
        VMTX2L_PHY (J)   = GMUV_OUT(I,3) 
        RAD_LONG_SIB_2_ATM_PHY (J) = - RAD_LONG_SIB_2_ATM_SIB (I)  
      ENDDO
!
      RETURN
      END SUBROUTINE SETBACK_TO_ATM

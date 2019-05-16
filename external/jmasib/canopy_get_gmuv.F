! 大気運動量拡散行列 GMUV の計算
! 99/08/12 保坂征宏
!
      SUBROUTINE CANOPY_GET_GMUV (
     I    NSTR       , NEND      , INSFLAG   ,
     I    RIB        , DELT_ATM  , DELT_CNP  , 
     I    U_REF      , V_REF     , Z_REF     , P_DEL      , 
     I    R_AIR_CONST_INV        , A_COEF    , RHO_AIR    , 
     I    Z0_msx     , 
     I    GMUV_IN    , 
     O    U_FLUX     , V_FLUX     , 
     O    GMUV   )
!
!   要チェック
!      ・U_FLUX と V_FLUX の単位は？
!      ・DELT_A = DELT_ATM * 2 は修正が必要。
!                               第一ステップは異なる。
!                               時間刻変更時もチェック要。
!
!-----------------------------------------------------------
!
      USE SIBPRM , ONLY : 
     1    IDX 
!
      USE SIBCON , ONLY :
     1    GRAV_c , 
     1    DELT_FAC_ATM_c   
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : MONIT_SIB, ICN, ICS,
     1                      IMONIT_LEVEL, IMONIT_ALL ,
     1                      IMONIT_FLUX_COMMON 
#endif
!
      IMPLICIT NONE 
!
! INPUT 
!
      INTEGER,INTENT(IN)    :: NSTR   
      INTEGER,INTENT(IN)    :: NEND
      INTEGER,INTENT(IN)    :: INSFLAG
      REAL(8),INTENT(IN)    :: RIB       (IDX)
      REAL(8),INTENT(IN)    :: U_REF     (IDX)
      REAL(8),INTENT(IN)    :: V_REF     (IDX)
      REAL(8),INTENT(IN)    :: Z_REF     (IDX)
      REAL(8),INTENT(IN)    :: Z0_msx    (IDX)
      REAL(8),INTENT(IN)    :: GMUV_IN   (IDX,3)
      REAL(8),INTENT(IN)    :: R_AIR_CONST_INV  (IDX)
      REAL(8),INTENT(IN)    :: RHO_AIR   (IDX)
      REAL(8),INTENT(IN)    :: A_COEF    (IDX)
      REAL(8),INTENT(IN)    :: P_DEL     (IDX)
      REAL(8),INTENT(IN)    :: DELT_ATM 
      REAL(8),INTENT(IN)    :: DELT_CNP

	REAL(8)               :: tmp
!
! OUTPUT
!
!
! INOUT
!
      REAL(8),INTENT(INOUT) :: GMUV          (IDX,2)
!
! WORK
!
      REAL(8),intent(out)   :: U_FLUX (IDX)
      REAL(8),intent(out)   :: V_FLUX (IDX)
!
      REAL(8) :: CUI           (IDX) 
      REAL(8) :: AA            (IDX) 
      REAL(8) :: BB            (IDX) 
      REAL(8) :: CC            (IDX) 
      REAL(8) :: DELT_A 
      REAL(8) :: R_MOM_AIR_INV 
      INTEGER :: I  
!
      REAL(8),SAVE :: CUI_MIN
      LOGICAL,SAVE :: LFIRST
!
#ifdef MONYOS
      CHARACTER(7)      :: CMARK 
      REAL(8),PARAMETER :: ONE  = 1.D0 
      INTEGER           :: IFLG_MON
#endif    
!
      DATA CUI_MIN / 0.1D0  /
      DATA LFIRST  / .TRUE. /
!
      IF ( LFIRST ) THEN
        WRITE(6,*) 'CANOPY_GET_GMUV_1 POCAHKA 99/08/12'  
        LFIRST = .FALSE. 
      ENDIF
!
!     DELTT = DELT_ATM * 3                   ! 3 は仮
      DELT_A = DELT_ATM * DELT_FAC_ATM_c
!
! 運動量フラックス計算
!
      DO 1000 I=NSTR,NEND
	  IF ( RIB(I).GT.0.0D0 ) THEN
          CUI(I)  = 1.0D0/(1.0D0+10.0D0*RIB(I)/SQRT(1.0D0+5.0D0*RIB(I)))
          CUI(I)  = MAX(CUI(I),CUI_MIN)
        ELSE
	    tmp = -Z_REF(I)/Z0_msx(I)*RIB(I)
	    if (tmp.le.0.0) then
	     tmp = 0.0001
	    endif
	    if (ISNAN(tmp)) then
	     tmp = 0.0001
	    endif 
	     CUI(I)  =  1.0D0-10.0D0*RIB(I)/( 1.0D0+75.0D0*A_COEF(I)**2
     1               *SQRT(tmp) )
	     CUI(I)  = MAX(CUI(I),CUI_MIN)
        END IF
 1000 CONTINUE
!
      DO 2000 I=NSTR,NEND 
        R_MOM_AIR_INV = CUI(I) * R_AIR_CONST_INV(I) 
        U_FLUX(I)     = - RHO_AIR(I) * R_MOM_AIR_INV * U_REF(I) 
        V_FLUX(I)     = - RHO_AIR(I) * R_MOM_AIR_INV * V_REF(I) 
!
! 方程式を建てる
!
!       AA(I) = P_DEL    (I) * 100 / GRAV_c * GMUV_IN(I,1)       
        AA(I) = GMUV_IN(I,1)       
     1          + RHO_AIR(I) * R_MOM_AIR_INV * DELT_A 

!       BB(I) = P_DEL    (I) * 100 / GRAV_c * GMUV_IN(I,2)       
        BB(I) = GMUV_IN(I,2)       
     1          + U_FLUX(I) 
!
!       CC(I) = P_DEL    (I) * 100 / GRAV_c * GMUV_IN(I,3)       
        CC(I) = GMUV_IN(I,3)       
     1          + V_FLUX(I) 
!
! 方程式を解く
!
        GMUV(I,1) = BB(I) / AA(I) 
        GMUV(I,2) = CC(I) / AA(I) 
!
! フラックス補正する
!
        U_FLUX(I) 
     1     = U_FLUX(I) 
     1       - RHO_AIR(I) * R_MOM_AIR_INV * DELT_A * GMUV(I,1) ! 010726 
!    1       + RHO_AIR(I) * R_MOM_AIR_INV * DELT_A * GMUV(I,1)
        V_FLUX(I) 
     1     = V_FLUX(I) 
     1       - RHO_AIR(I) * R_MOM_AIR_INV * DELT_A * GMUV(I,2) ! 010726 
!    1       + RHO_AIR(I) * R_MOM_AIR_INV * DELT_A * GMUV(I,2)
!
 2000 CONTINUE
!
!
!  ================
!  >> モニタする <<
!  ================
! 
#ifdef MONYOS
      IF ( INSFLAG .EQ.0 ) THEN
        IFLG_MON    = ICN
      ELSE
        IFLG_MON    = ICS
      ENDIF
!
!      IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN 
!        CMARK  = 'UMOM'
!        CALL MONIT_SIB ( CMARK , U_FLUX ,IFLG_MON, DELT_CNP , ONE )      
!        CMARK  = 'VMOM'
!        CALL MONIT_SIB ( CMARK , V_FLUX ,IFLG_MON, DELT_CNP , ONE )      
!      ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
        IF ( INSFLAG .EQ. 0 ) THEN
          CMARK  = 'MUFLG'
          CALL MONIT_SIB ( CMARK , U_FLUX ,IFLG_MON, DELT_CNP , ONE )      
          CMARK  = 'MVFLG'
          CALL MONIT_SIB ( CMARK , V_FLUX ,IFLG_MON, DELT_CNP , ONE )      
        ELSE
          CMARK  = 'MUFLW'
          CALL MONIT_SIB ( CMARK , U_FLUX ,IFLG_MON, DELT_CNP , ONE )      
          CMARK  = 'MVFLW'
          CALL MONIT_SIB ( CMARK , V_FLUX ,IFLG_MON, DELT_CNP , ONE )      
        ENDIF
      ENDIF
#endif
!
      RETURN
      END 
!================================================================
      SUBROUTINE     MR_CANOPY_GET_GMUV 
!
!
#ifdef MONYOS
      USE SIB_MONIT, ONLY : 
     1                      IMONIT_LEVEL, IMONIT_ALL ,
     1                      IMONIT_FLUX_COMMON 
#endif
!
! デバッグ用変数の登録作業
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
#ifdef MONYOS
      IF ( LFIRST ) THEN
        WRITE(6,*) 'MONIT_REGIST_CANOPY_GET_GMUV' , 
     1             ' 01/10 POCHAKA'
        LFIRST = .FALSE. 
      ENDIF
!
!     IF ( IMONIT_LEVEL .GE. IMONIT_FLUX_COMMON ) THEN
!
!  運動量フラックス ( canopy_get_gmuv )
!     CMARK='UMOM' = MFLUG + MFLUW    登録は海と共通
!     CMARK='VMOM' = MFLVG + MFLVW    登録は海と共通
!
!     ENDIF
!
      IF ( IMONIT_LEVEL .GE. IMONIT_ALL ) THEN 
!
        CMARK  = 'MUFLG'
        CTITLE = 'Momentum flux (U) (NOS)'
        CUNIT  = 'N/M**2/S'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'MVFLG'
        CTITLE = 'Momentum flux (V) (NOS)'
        CUNIT  = 'N/M**2/S'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'MUFLW'
        CTITLE = 'Momentum flux (U) (SNW)'
        CUNIT  = 'N/M**2/S'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
        CMARK  = 'MVFLW'
        CTITLE = 'Momentum flux (V) (SNW)'
        CUNIT  = 'N/M**2/S'
        CALL monit_regist_sib ( CMARK , CTITLE , CUNIT )
!
      ENDIF
#endif
!
      END SUBROUTINE MR_CANOPY_GET_GMUV 


